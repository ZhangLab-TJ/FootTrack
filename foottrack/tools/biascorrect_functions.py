import os
import sys
import gc
import numpy as np
import multiprocessing as mp
import time
from datetime import datetime
import matplotlib
matplotlib.use('Agg')	#non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit
import pickle
import pyBigWig
import math
import random
#Bio-specific packages
import pysam

#Internal functions and classes
from foottrack.utils.sequences import SequenceMatrix, GenomicSequence
from foottrack.utils.signals import fast_rolling_math
from foottrack.utils.utilities import * 
from foottrack.utils.regions import OneRegion, RegionList
from foottrack.utils.ngs import OneRead, ReadList
from foottrack.utils.logger import foottrackLogger

#Catch warnings from curve_fit
import warnings
from scipy.optimize import OptimizeWarning
warnings.simplefilter("error", OptimizeWarning)

#--------------------------------------------------------------------------------------------------#
class enzBias:
	""" Class for storing information about estimated bias """
	def __init__(self, L=10, stype="DWM"):
		self.stype = stype	#type of score matrix
		self.bias = {"forward": SequenceMatrix.create(L, self.stype),
					 "reverse": SequenceMatrix.create(L, self.stype),
					 "both": SequenceMatrix.create(L, self.stype)}
		self.no_crs = 0
	def join(self, obj):
		""" Join counts from enzBias obj with other enzBias obj """
		self.bias["forward"].add_counts(obj.bias["forward"])
		self.bias["reverse"].add_counts(obj.bias["reverse"])
		self.bias["both"].add_counts(obj.bias["both"])
		self.no_crs += obj.no_crs
	def to_pickle(self, f):
		""" Pickle an enzBias object to a .pickle file """
		handle = open(f, "wb")
		pickle.dump(self, handle)
		handle.close()
		return(self)
	def from_pickle(self, f):
		""" Read anq enzBias object from a .pickle file """
		handle = open(f, "rb")
		self = pickle.load(handle)
		return(self)
	
class CrList(list):
	def __init__(self, data_array=[]):
		super().__init__(data_array)
		self.bw_obj = None
	@classmethod
	def from_bw(cls, bw_obj, region):
		""" Create a CrList from a bigWig file for a specific region """
		new_instance = cls()
		chrom, start, end = region
		values = bw_obj.values(chrom, start, end)
		for i, val in enumerate(values):
			if val is not None and not math.isnan(val):
				new_instance.append(OneCr(chrom, start + i, val))
		return new_instance
	def split_strands(self):
		""" Split CrList into forward/reverse crs """
		forward_cr = [cr for cr in self if cr.is_reverse == False]
		reverse_cr = [cr for cr in self if cr.is_reverse == True]
		return [CrList(forward_cr), CrList(reverse_cr)]
	def signal(self, region):
		chrom, reg_start, reg_end = region
		reg_len = reg_end - reg_start
		values = np.full(reg_len, np.nan)
		for cr_obj in self:
			cut = cr_obj.pos
			if reg_start < cut <= reg_end:
				values[cut - reg_start] = cr_obj.value
		return values
	def bias(self, region):
		chrom, reg_start, reg_end = region
		reg_len = reg_end - reg_start
		bias = np.full(reg_len, np.nan)
		for cr_obj in self:
			cut = cr_obj.pos
			if reg_start < cut <= reg_end:
				bias[cut - reg_start] = cr_obj.bias
		return bias
	def close_bw(self):
		""" Close the bigWig file """
		self.bw_obj.close()

class OneCr():
	""" Class representing a single record from bigWig file """
	def __init__(self, chrom, pos, value):
		self.chrom = chrom
		self.pos = pos
		self.value = value
		self.is_reverse = None	# Initialize strand as unknown
		self.kmer = None  # Initialize kmer as None
		self.bias = None
	def get_strand(self, genomic_sequence):
		""" Infer strand information based on the genomic sequence at a position """
		seq_start = genomic_sequence.region.start
		seq_end = genomic_sequence.region.end
		i = self.pos - seq_start
		sequence = genomic_sequence.sequence[i:i+1]  # Assuming genomic_sequence.sequence uses 0-based indexing
		if sequence == 2:
			self.is_reverse = False	# Not reverse strand
		elif sequence == 3:
			self.is_reverse = True	# Is reverse strand			
	def get_kmer(self, genomic_sequence, k_flank):
		""" Get the kmer sequence around a given position """
		seq_start = genomic_sequence.region.start
		seq_end = genomic_sequence.region.end
		if self.pos > seq_start + k_flank + 1 and self.pos <= seq_end - k_flank:
			if self.is_reverse == False:
				i = self.pos - seq_start
				self.kmer = genomic_sequence.sequence[i-k_flank:i+k_flank+1]
			else:
				i = seq_end - self.pos - 1
				self.kmer = genomic_sequence.revcomp[i-k_flank:i+k_flank+1]
		else:
			self.kmer = np.array([4] * (k_flank*2 + 1))
	def get_bias(self, bias_obj):
		if self.is_reverse == False:
			strand = "forward"
		else:
			strand = "reverse"
		bias_array = bias_obj.bias[strand].score_sequence(self.kmer)
		if len(bias_array) == 0:
			self.bias = np.nan
		else:
			length = len(bias_array)
			middle_index = length // 2
			self.bias = bias_array[middle_index]

#--------------------------------------------------------------------------------------------------#

def count_crs(regions_list, params):
	""" Count crs from bw within regions (counts position of pos to prevent double-counting) """

	bw_f = params.bw
	bw_obj = pyBigWig.open(bw_f)

	log_q = params.log_q
	logger = foottrackLogger("", params.verbosity, log_q) #sending all logger calls to log_q

	#Count per region
	cr_count = 0
	logger.spam("Started counting region_chunk ({0} -> {1})".format("_".join([str(element) for element in regions_list[0]]), "_".join([str(element) for element in regions_list[-1]])))
	
	for region in regions_list:
		cr_lst = CrList.from_bw(bw_obj, region)
		logger.spam("- {0} ({1} crs)".format(region, len(cr_lst)))
		for cr in cr_lst:  
			if cr.pos > region.start and cr.pos <= region.end:  #only crs within borders
				cr_count += 1
				
	logger.spam("Finished counting region_chunk ({0} -> {1})".format("_".join([str(element) for element in regions_list[0]]), "_".join([str(element) for element in regions_list[-1]])))
	bw_obj.close()

	return(cr_count)

#--------------------------------------------------------------------------------------------------#

def bias_estimation(regions_lst, params):
	""" Estimates bias of insertions within regions """

	#Info on run
	bw_f = params.bw
	fasta_f = params.genome
	k_flank = params.k_flank
	bg_shift = params.bg_shift
	L = 2 * k_flank + 1

	logger = foottrackLogger("", params.verbosity, params.log_q) 	#sending all logger calls to log_q

	#Open objects for reading
	bw_obj = pyBigWig.open(bw_f)
	fasta_obj = pysam.FastaFile(fasta_f)
	chrom_lengths = bw_obj.chroms()	#Chromosome boundaries from bw_obj

	bias_obj = enzBias(L, params.score_mat)

	strands = ["forward", "reverse"]

	#Estimate bias at each region
	for region in regions_lst:

		cr_lst = CrList().from_bw(bw_obj, region)

		## Kmer cutting bias ##
		if len(cr_lst) > 0:

			#Extract sequence
			extended_region = region.extend_reg(k_flank + bg_shift)	#Extend to allow full kmers
			extended_region.check_boundary(chrom_lengths, "cut")
			sequence_obj = GenomicSequence(extended_region).from_fasta(fasta_obj)
			for cr in cr_lst:
				cr.get_strand(sequence_obj)
				cr.get_kmer(sequence_obj, k_flank)
			
			#Split crs forward/reverse
			for_lst, rev_lst = cr_lst.split_strands()
			cr_lst_strand = {"forward": for_lst, "reverse": rev_lst} 
			logger.spam("Region: {0}. Forward crs: {1}. Reverse crs: {2}".format(region, len(for_lst), len(rev_lst)))

			for strand in strands:
				for cr in cr_lst_strand[strand]:
					if cr.pos > region.start and cr.pos < region.end:  	#only crs within borders
						no_cut = cr.value
						cr.get_kmer(sequence_obj, k_flank)
						bias_obj.bias[strand].add_sequence(cr.kmer, no_cut)
						bias_obj.bias[strand].add_background(cr.kmer, no_cut)
						bias_obj.no_crs += no_cut

	bw_obj.close()
	fasta_obj.close()

	return(bias_obj) 	#object containing information collected on bias 	

#--------------------------------------------------------------------------------------------------#

def relu(x, a, b):
	""" a and b are components of a linear curve (y=a*x+b) """
	y = np.maximum(0.0, a*x + b)
	return(y)

#--------------------------------------------------------------------------------------------------#

def bias_correction(regions_lst, params, bias_obj, golobal_cr):
	""" Corrects bias in poss (from bwfile) using estimated bias """

	logger = foottrackLogger("", params.verbosity, params.log_q) 	#sending all logger calls to log_q

	bw_f = params.bw
	fasta_f = params.genome
	k_flank = params.k_flank
	L = 2 * k_flank + 1
	w = params.window
	f = int(w/2.0)
	qs = params.qs

	f_extend = k_flank + f

	strands = ["forward", "reverse"]
	pre_bias = {strand: SequenceMatrix.create(L, "PWM") for strand in strands}
	post_bias = {strand: SequenceMatrix.create(L, "PWM") for strand in strands}

	#Open bwfile and fasta
	bw_obj = pyBigWig.open(bw_f)
	fasta_obj = pysam.FastaFile(fasta_f)

	out_signals = {}

	#Go through each region
	for region_obj in regions_lst:

		region_obj.extend_reg(f_extend)
		reg_len = region_obj.get_length()	#length including flanking
		reg_key = (region_obj.chrom, region_obj.start+f_extend, region_obj.end-f_extend)	#output region
		out_signals[reg_key] = {"bias":{}, "corrected_lc":{}, "corrected_gl":{}, "expected_lc":{}, "expected_gl":{}}

		#Get pos positions for each cr
		cr_lst = CrList().from_bw(bw_obj, region_obj)
		#Get sequence in this region
		sequence_obj = GenomicSequence(region_obj).from_fasta(fasta_obj)
		for cr in cr_lst:
			cr.get_strand(sequence_obj)
			cr.get_kmer(sequence_obj, k_flank)
			cr.get_bias(bias_obj)
		logger.spam("Cr {0} crs from region {1}".format(len(cr_lst), region_obj))

		for_lst, rev_lst = cr_lst.split_strands()
		cr_lst_strand = {"forward": for_lst, "reverse": rev_lst}

		for strand in strands:

			########################################
			####### Uncorrected crs  and bias ######
			########################################

			uncorrected_signal = cr_lst_strand[strand].signal(region_obj)
			bias_log = cr_lst_strand[strand].bias(region_obj)
			bias = np.power(2, bias_log)
			out_signals[reg_key]["bias"][strand] = bias

			#################################
			###### Correction of crs ######
			#################################

			# local background
			signal_mean = fast_rolling_math(uncorrected_signal, w, "mean")
			expected_lc = signal_mean * bias
			corrected_lc = uncorrected_signal - expected_lc
			
			out_signals[reg_key]["expected_lc"][strand] = expected_lc
			out_signals[reg_key]["corrected_lc"][strand] = corrected_lc

			# global background
			expected_gl = golobal_cr * bias
			corrected_gl = uncorrected_signal - expected_gl
			out_signals[reg_key]["expected_gl"][strand] = expected_gl
			out_signals[reg_key]["corrected_gl"][strand] = corrected_gl

			###############################################################

			# ######## Correct signal ########
			
			# out_signals[reg_key]["expected"][strand] = np.zeros_like(uncorrected_signal)
			# out_signals[reg_key]["corrected"][strand] = uncorrected_signal / bias_probas

		#######################################
		########   Verify correction   ########
		#######################################
		
		#Verify correction across all crs

			for idx in range(k_flank,reg_len - k_flank): 

				orig = uncorrected_signal[idx]
				correct = corrected_lc[idx]

				if not math.isnan(orig) and not math.isnan(correct):	#if one is nan, don't add to pre/post bias
					if strand == "forward":
						i = idx
						kmer = sequence_obj.sequence[i-k_flank:i+k_flank+1]
					else:
						i = reg_len - idx - 1
						kmer = sequence_obj.revcomp[i-k_flank:i+k_flank+1]

					#Save kmer for bias correction verification
					pre_bias[strand].add_sequence(kmer, orig)
					post_bias[strand].add_sequence(kmer, correct)

		#######################################
		########    Write to queue    #########
		#######################################
		#Set size back to original
		for track in out_signals[reg_key]:
			for strand in out_signals[reg_key][track]:
				out_signals[reg_key][track][strand] = out_signals[reg_key][track][strand][f_extend:-f_extend]

		#Calculate "both" if split_strands == False
		if params.split_strands == False:
			for track in out_signals[reg_key]:
				forward = out_signals[reg_key][track]["forward"]
				reverse = out_signals[reg_key][track]["reverse"]				
				# sum_without_nans = np.nan_to_num(forward) + np.nan_to_num(reverse)
				sum_without_nans = np.nansum([forward, reverse], axis=0)
				combined_sum = np.where(np.isnan(forward) & np.isnan(reverse), np.nan, sum_without_nans)
				out_signals[reg_key][track]["both"] = combined_sum

		#Send to queue
		strands_to_write = ["forward", "reverse"] if params.split_strands == True else ["both"]
		for track in out_signals[reg_key]:

			#Send to writer per strand
			for strand in strands_to_write:
				key = "{0}:{1}".format(track, strand)

				if key in qs: #only write the signals where the files were initialized
					logger.spam("Sending {0} signal from region {1} to writer queue".format(key, reg_key))
					qs[key].put((key, reg_key, out_signals[reg_key][track][strand]))

		#Sent to qs - delete from this process
		out_signals[reg_key] = None
	
	bw_obj.close()
	fasta_obj.close()

	gc.collect()

	return([pre_bias, post_bias])

def bias_generate(regions_lst, params, bias_obj):
	""" Corrects bias in poss (from bwfile) using estimated bias """

	logger = foottrackLogger("", params.verbosity, params.log_q) 	#sending all logger calls to log_q

	bw_f = params.bw
	fasta_f = params.genome
	k_flank = params.k_flank
	L = 2 * k_flank + 1
	w = params.window
	f = int(w/2.0)
	qs = params.qs

	f_extend = k_flank + f

	strands = ["forward", "reverse"]

	#Open bwfile and fasta
	bw_obj = pyBigWig.open(bw_f)
	fasta_obj = pysam.FastaFile(fasta_f)

	out_signals = {}

	#Go through each region
	for region_obj in regions_lst:

		region_obj.extend_reg(f_extend)
		reg_key = (region_obj.chrom, region_obj.start+f_extend, region_obj.end-f_extend)	#output region
		out_signals[reg_key] = {"bias":{}}

		#Get pos positions for each cr
		cr_lst = CrList().from_bw(bw_obj, region_obj)
		#Get sequence in this region
		sequence_obj = GenomicSequence(region_obj).from_fasta(fasta_obj)
		for cr in cr_lst:
			cr.get_strand(sequence_obj)
			cr.get_kmer(sequence_obj, k_flank)
			cr.get_bias(bias_obj)
		logger.spam("Cr {0} crs from region {1}".format(len(cr_lst), region_obj))

		for_lst, rev_lst = cr_lst.split_strands()
		cr_lst_strand = {"forward": for_lst, "reverse": rev_lst}

		for strand in strands:

			########################################
			####### Uncorrected crs  and bias ######
			########################################

			bias_log = cr_lst_strand[strand].bias(region_obj)
			bias = np.power(2, bias_log)
			out_signals[reg_key]["bias"][strand] = bias 

		#######################################
		########    Write to queue    #########
		#######################################
		#Set size back to original
		for track in out_signals[reg_key]:
			for strand in out_signals[reg_key][track]:
				out_signals[reg_key][track][strand] = out_signals[reg_key][track][strand][f_extend:-f_extend]

		#Calculate "both" if split_strands == False
		if params.split_strands == False:
			for track in out_signals[reg_key]:
				forward = out_signals[reg_key][track]["forward"]
				reverse = out_signals[reg_key][track]["reverse"]				
				# sum_without_nans = np.nan_to_num(forward) + np.nan_to_num(reverse)
				sum_without_nans = np.nansum([forward, reverse], axis=0)
				combined_sum = np.where(np.isnan(forward) & np.isnan(reverse), np.nan, sum_without_nans)
				out_signals[reg_key][track]["both"] = combined_sum

		#Send to queue
		strands_to_write = ["forward", "reverse"] if params.split_strands == True else ["both"]
		for track in out_signals[reg_key]:

			#Send to writer per strand
			for strand in strands_to_write:
				key = "{0}:{1}".format(track, strand)

				if key in qs: #only write the signals where the files were initialized
					logger.spam("Sending {0} signal from region {1} to writer queue".format(key, reg_key))
					qs[key].put((key, reg_key, out_signals[reg_key][track][strand]))

		#Sent to qs - delete from this process
		out_signals[reg_key] = None

	bw_obj.close()
	fasta_obj.close()
	gc.collect()

####################################################################################################
######################################## Plot functions ############################################
####################################################################################################

colors = {0:"green", 1:"red", 2:"blue", 3:"darkkhaki"}
names = {0:"A", 1:"T", 2:"C", 3:"G"}

def plot_pssm(matrix, title):
	""" Plot pssm in matrix """
	#Make figure
	fig, ax = plt.subplots()
	fig.suptitle(title, fontsize=16, weight="bold")
	#Formatting of x axis
	length = matrix.shape[1]
	flank = int(length/2.0)		
	xvals = np.arange(length)  # each position corresponds to i in mat	
	#Customize minor tick labels
	xtick_pos = xvals[:-1] + 0.5
	xtick_labels = list(range(-flank, flank+1))
	ax.xaxis.set_major_locator(ticker.FixedLocator(xvals))
	ax.xaxis.set_major_formatter(ticker.FixedFormatter(xtick_labels))	
	ax.xaxis.set_minor_locator(ticker.FixedLocator(xtick_pos))			#locate minor ticks between major ones (cutsites)
	ax.xaxis.set_minor_formatter(ticker.NullFormatter())
	#Make background grid on major ticks
	plt.grid(color='0.8', which="minor", ls="--", axis="x")
	plt.xlim([0, length-1])
	# plt.ylim([-0.26, 0.26])
	plt.xlabel('Position from cytosine')
	plt.ylabel('PSSM score')
	######## Plot data #######
	#Plot PSSM / bias motif
	for nuc in range(4):
		plt.plot(xvals, matrix[nuc,:], color=colors[nuc], label=names[nuc])
	#Cutsite-line
	#Finish up
	plt.legend(loc="lower right")
	plt.tight_layout()
	fig.subplots_adjust(top=0.88, hspace=0.5)
	return(fig)

#----------------------------------------------------------------------------------------------------#
def plot_correction(pre_mat, post_mat, title):
	""" Plot comparison of pre-correction and post-correction matrices """

	#Make figure
	fig, ax = plt.subplots()
	fig.suptitle(title, fontsize=16, weight="bold")

	L = pre_mat.shape[1]
	flank = int(L/2.0)		
	xvals = np.arange(L)  # each position corresponds to i in mat
	
	#Customize minor tick labels
	xtick_pos = xvals[:-1] + 0.5
	xtick_labels = list(range(-flank, flank))		#-flank - flank without 0
	ax.xaxis.set_major_locator(ticker.FixedLocator(xvals))
	ax.xaxis.set_major_formatter(ticker.FixedFormatter(xtick_labels))	
	ax.xaxis.set_minor_locator(ticker.FixedLocator(xtick_pos))			#locate minor ticks between major ones (conversion rate)
	ax.xaxis.set_minor_formatter(ticker.NullFormatter())
	
	#PWMs for all mats
	pre_pwm = pre_mat
	post_pwm = post_mat

	#Pre correction
	for nuc in range(4):
		yvals = [pre_pwm[nuc, m] for m in range(L)]
		plt.plot(xvals, yvals, linestyle="--", color=colors[nuc], linewidth=1, alpha=0.5)
	
	#Post correction
	for nuc in range(4):
		yvals = [post_pwm[nuc, m] for m in range(L)]
		plt.plot(xvals, yvals, color=colors[nuc], linewidth=2, label=names[nuc])
	
	plt.xlim([0, L-1])
	plt.xlabel('Position from cutsite')
	plt.ylabel('Nucleotide frequency')

	#Set legend
	plt.plot([0],[0], linestyle="--", linewidth=1, color="black", label="pre-correction")
	plt.plot([0],[0], color="black", label="post-correction")
	plt.legend(loc="lower right", prop={'size':6})

	plt.tight_layout()
	fig.subplots_adjust(top=0.88, hspace=0.4)

	return(fig)
