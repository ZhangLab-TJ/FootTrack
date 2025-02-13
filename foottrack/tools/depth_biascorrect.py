import argparse
import collections
from collections import OrderedDict
import copy
from copy import deepcopy
import csv
import gc
import math
import multiprocessing as mp
import os
import pickle
import random
import re
import subprocess
import sys
import time
from datetime import datetime
from itertools import *
import warnings

import matplotlib
matplotlib.use('Agg')  # non-interactive backend
import numpy as np
import pandas as pd
import pysam
import scipy.stats
from scipy.optimize import curve_fit, OptimizeWarning
from scipy.stats import binom
from statsmodels.stats.multitest import fdrcorrection

warnings.simplefilter("error", OptimizeWarning)

from foottrack.parsers import add_depth_biascorrect_arguments
from foottrack.tools.biascorrect_functions import *
from foottrack.utils.logger import *
from foottrack.utils.ngs import OneRead, ReadList
from foottrack.utils.regions import OneRegion, RegionList
from foottrack.utils.sequences import *
from foottrack.utils.signals import fast_rolling_math
from foottrack.utils.utilities import *

# Class and Function

class OneCrd():
	def __init__(self, chrom, pos, strand, conversion, convert_N, total_N):
		self.chrom = chrom
		self.pos = pos
		self.strand = strand
		self.conversion = conversion
		self.convert_N = convert_N
		self.total_N = total_N
	def get_strand(self):
		if self.strand == "+":
			self.is_reverse = False # Not reverse strand
		else:
			self.is_reverse = True # Is reverse strand
	def get_kmer(self, genomic_sequence, k_flank):
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
	def get_field(self, field_name):
		return getattr(self, field_name, np.nan)

class CrdList(list):
	def __init__(self, data_array=[]):
		super().__init__(data_array)
		self.csv_obj = None
	@classmethod
	def from_csv(cls, meth_file, region):
		new_instance = cls()
		chrom, start, end = region
		region_data = meth_file[(meth_file['chrom'] == chrom) & (meth_file['start'] >= start) & (meth_file['start'] < end)]
		for index, row in region_data.iterrows():
			one_crd = OneCrd(
				chrom=row['chrom'],
				pos=row['start'],
				strand=row['strand'],
				conversion=row['conversion'],
				convert_N=row['convert_N'],
				total_N=row['total_N']
			)
			new_instance.append(one_crd)
		return new_instance
	def split_strands(self):
		forward_crd = [crd for crd in self if crd.is_reverse == False]
		reverse_crd = [crd for crd in self if crd.is_reverse == True]
		return [CrdList(forward_crd), CrdList(reverse_crd)]
	def signal(self, region, type):
		chrom, reg_start, reg_end = region
		reg_len = reg_end - reg_start
		values = np.full(reg_len, np.nan)
		for cdr_obj in self:
			cut = cdr_obj.pos
			if reg_start < cut <= reg_end:
				values[cut - reg_start] = cdr_obj.get_field(type)
		return values
	def bias(self, region):
		chrom, reg_start, reg_end = region
		reg_len = reg_end - reg_start
		bias = np.full(reg_len, np.nan)
		for crd_obj in self:
			cut = crd_obj.pos
			if reg_start < cut <= reg_end:
				bias[cut - reg_start] = crd_obj.bias
		return bias

def correct_and_pval(meth_f, regions_lst, w, bias_f, fasta_f, k_flank, split_strands, standard, mode):

	strands = ["forward", "reverse"]
	L = 2 * k_flank + 1
	f = int(w/2.0)
	f_extend = k_flank + f
	pre_bias = {strand: SequenceMatrix.create(L, "PWM") for strand in strands}
	post_bias = {strand: SequenceMatrix.create(L, "PWM") for strand in strands}
	
	# Open fasta, bias
	fasta_obj = pysam.FastaFile(fasta_f)
	bias_obj = enzBias().from_pickle(bias_f)

	out_signals = {}
	for region_obj in regions_lst: #Go through each region
		region_obj.extend_reg(f_extend)
		reg_len = region_obj.get_length()
		reg_key = (region_obj.chrom, region_obj.start+f_extend, region_obj.end-f_extend)
		
		# Initialize out_signals for the region
		out_signals[reg_key] = {
			"p_value": {}, 
			"fdr": {}, 
			"bias": {}, 
			"expected": {}, 
			"corrected": {}
		}

		crd_lst = CrdList.from_csv(meth_f, region_obj) #Get pos positions for each cr
		sequence_obj = GenomicSequence(region_obj).from_fasta(fasta_obj) #Get sequence in this region
		for crd in crd_lst:
			crd.get_strand()
			crd.get_kmer(sequence_obj, k_flank)
			crd.get_bias(bias_obj)
		
		for_lst, rev_lst = crd_lst.split_strands()
		crd_lst_strand = {"forward": for_lst, "reverse": rev_lst}
		
		for strand in strands:

			########################################
			####### Uncorrected crs and bias #######
			########################################
			uncorrected_signal = crd_lst_strand[strand].signal(region_obj, 'conversion')
			bias_log = crd_lst_strand[strand].bias(region_obj)
			bias = np.power(2, bias_log)
			out_signals[reg_key]["bias"][strand] = bias

			#################################
			###### Correction of crds #######
			#################################
			if mode == "local":
				signal_mean = fast_rolling_math(uncorrected_signal, w, "mean")
				expected = signal_mean * bias
				corrected = uncorrected_signal - expected
				out_signals[reg_key]["expected"][strand] = expected
				out_signals[reg_key]["corrected"][strand] = corrected
			elif mode == "global":
				expected = standard * bias
				corrected = uncorrected_signal - expected
				out_signals[reg_key]["expected"][strand] = expected
				out_signals[reg_key]["corrected"][strand] = corrected

			###################################
			######## Verify correction ########
			###################################
			for idx in range(k_flank,reg_len - k_flank): 
				orig = uncorrected_signal[idx]
				correct = corrected[idx]
				if not math.isnan(orig) and not math.isnan(correct): #if one is nan, don't add to pre/post bias
					if strand == "forward":
						i = idx
						kmer = sequence_obj.sequence[i-k_flank:i+k_flank+1]
					else:
						i = reg_len - idx - 1
						kmer = sequence_obj.revcomp[i-k_flank:i+k_flank+1]
					#Save kmer for bias correction verification
					pre_bias[strand].add_sequence(kmer, orig)
					post_bias[strand].add_sequence(kmer, correct)

			###################################
			######## Calculate p-value ########
			###################################
			N_total_array = crd_lst_strand[strand].signal(region_obj, 'total_N')
			N_convert_array = crd_lst_strand[strand].signal(region_obj, 'convert_N')
			expected_array = out_signals[reg_key]["expected"][strand]

			region_pvalue = np.full(len(N_convert_array), np.nan)
			valid_mask = (~np.isnan(N_total_array)) & (N_total_array != 0) & (~np.isnan(expected_array))

			for i in range(len(valid_mask)):
				if valid_mask[i]:  # Only process valid data points
					n_covert = int(N_convert_array[i])
					n_total = int(N_total_array[i])
					p_exp = np.clip(expected_array[i], 0, 1)
					pvalue = scipy.stats.binomtest(n_covert, n_total, p=p_exp, alternative='less').pvalue # Perform the binomial test
					region_pvalue[i] = pvalue # Assign the computed p-value to the region_pvalue array

			out_signals[reg_key]["p_value"][strand] = region_pvalue
			out_signals[reg_key]["fdr"][strand] = np.full(len(N_convert_array), np.nan)

		####################################
		######## not split strands #########
		####################################
		#Set size back to original
		for track in out_signals[reg_key]:
			for strand in out_signals[reg_key][track]:
				out_signals[reg_key][track][strand] = out_signals[reg_key][track][strand][f_extend:-f_extend]
		#Calculate "both" if split_strands == False
		if split_strands == False:
			for track in out_signals[reg_key]:
				forward = out_signals[reg_key][track]["forward"]
				reverse = out_signals[reg_key][track]["reverse"]
				sum_without_nans = np.nansum([forward, reverse], axis=0)
				sum = np.where(np.isnan(forward) & np.isnan(reverse), np.nan, sum_without_nans)
				out_signals[reg_key][track]["both"] = sum

	return(out_signals)

def calculate_fdr(df, pval_column='p_value', alpha=0.05, method='indep'):
	pvals = df[pval_column].values
	_, fdr_corrected = fdrcorrection(pvals, alpha=alpha, method=method)
	df = df.copy()
	df['fdr'] = fdr_corrected
	return df

# run_depth_biascorrect

def run_depth_biascorrect(args):
	"""
	Main function to run correction and p-value calculation.
	
	Parameters:
	- args: parsed command-line arguments
	"""
	# Extract arguments
	input_f = args.input
	fasta_f = args.genome
	bed_f = args.peaks
	bias_f = args.bias
	w = args.window
	k_flank = args.k_flank
	prefix = args.prefix
	split_strands = args.split_strands
	verbosity = args.verbosity
	cores = args.cores
	mode = args.mode

	# Print info on run
	logger = foottrackLogger("DepthBiascorrect", verbosity)
	logger.begin()
	parser = add_depth_biascorrect_arguments(argparse.ArgumentParser())
	logger.arguments_overview(parser, args)

	# Open bed, meth_file
	peak_regions = RegionList().from_bed(bed_f)
	regions_lst = peak_regions.apply_method(OneRegion.split_region, 50000)
	output_regions_chunks = regions_lst.chunks(args.split)
	names = ['chrom', 'start', 'end', 'strand', 'conversion', 'convert_N', 'total_N']
	logger.info("Reading info from .txt file")
	meth_f = pd.read_csv(input_f, sep='\t', header=0, names=names, compression='infer',low_memory=False)
	if args.standard:
		standard = args.standard
	else:
		standard = meth_f['conversion'].mean()
	logger.info("Read done")

	logger.info("Start calculate correct signal and p-value")
	worker_pool = mp.Pool(processes=cores)
	task_list = [worker_pool.apply_async(correct_and_pval, args=[meth_f, chunk, w, bias_f, fasta_f, k_flank, split_strands, standard, mode]) for chunk in output_regions_chunks]
	worker_pool.close()
	monitor_progress(task_list, logger, "Correction progress:")	#does not exit until tasks in task_list finished
	results = [task.get() for task in task_list]
	# Combine all results into a single dictionary
	global_out_signals = {}
	for out_signals in results:
		global_out_signals.update(out_signals)
	logger.info("Calculate done")

	# logger.info("Start calculate correct signal and p-value")
	# with mp.Pool(processes=cores) as pool:
	# 	task_list = pool.starmap(correct_and_pval, [(meth_f, chunk, w, bias_f, fasta_f, k_flank, split_strands) for chunk in output_regions_chunks])
	# # Combine all results into a single dictionary
	# global_out_signals = {}
	# for out_signals in task_list:
	# 	global_out_signals.update(out_signals)
	# logger.info("Calculate done")

	# # one core
	# logger.info("Start calculate correct singal and p-value")
	# global_out_signals = correct_and_pval(meth_f, peak_regions, w, bias_f, fasta_f, k_flank, split_strands)
	# logger.info("Calculate done")

	# output
	output_data = {
		"forward": [],
		"reverse": [],
		"both": []
	}

	logger.info("Start calculate FDR")
	for reg_key, tracks in global_out_signals.items():
		for strand in ["forward", "reverse", "both"]:
			corrected_array = tracks["corrected"].get(strand, np.array([]))
			p_value_array = tracks["p_value"].get(strand, np.array([]))
			fdr_value_array = tracks["fdr"].get(strand, np.array([]))

			for idx in range(len(corrected_array)):
				corrected = corrected_array[idx]
				p_value = p_value_array[idx]
				fdr = fdr_value_array[idx]
				
				if not np.isnan(corrected) and not np.isnan(p_value):
					position = reg_key[1] + idx
					data_entry = {
						"region": reg_key[0],
						"start_position": position,
						"end_position": position + 1,
						"corrected": corrected,
						"p_value": p_value,
						"fdr": fdr
					}
					output_data[strand].append(data_entry)

	forward_df = pd.DataFrame(output_data["forward"])
	reverse_df = pd.DataFrame(output_data["reverse"])
	both_df = pd.DataFrame(output_data["both"])
	# Calculate and fill FDR for each DataFrame
	if not forward_df.empty:
		forward_df = calculate_fdr(forward_df)
	if not reverse_df.empty:
		reverse_df = calculate_fdr(reverse_df)
	if not both_df.empty:
		both_df = calculate_fdr(both_df)
	logger.info("Calculate done")

	logger.info("output")
	if split_strands == False:
		output_file = prefix + ".txt"
		both_df.to_csv(output_file, sep='\t', index=False)
	else:
		output_file_forward = prefix + "_forward.txt"
		forward_df.to_csv(output_file_forward, sep='\t', index=False)
		output_file_reverse = prefix + "_reverse.txt"
		reverse_df.to_csv(output_file_reverse, sep='\t', index=False)
	logger.info("DepthBiascorrect Done!")

#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser = add_depth_biascorrect_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	run_depth_biascorrect(args)