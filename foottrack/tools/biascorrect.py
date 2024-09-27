#!/usr/bin/env python

"""
BiasCorrect.py: Estimates cFOOT-seq bias and corrects cr counts from .bw and .fasta input

"""

#--------------------------------------------------------------------------------------------------------#
#----------------------------------------- Import libraries ---------------------------------------------# 
#--------------------------------------------------------------------------------------------------------#

import os
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 

import sys
import argparse
import numpy as np
import multiprocessing as mp
from copy import deepcopy

from collections import OrderedDict
import itertools
import matplotlib
matplotlib.use("Agg")  #non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#Bio-specific packages
import pyBigWig
import pysam

#Internal functions and classes
from foottrack.parsers import add_biascorrect_arguments
from foottrack.tools.biascorrect_functions import *
from foottrack.utils.utilities import *
from foottrack.utils.regions import OneRegion, RegionList
from foottrack.utils.sequences import *
from foottrack.utils.logger import foottrackLogger

#np.seterr(divide='raise', invalid='raise')

#--------------------------------------------------------------------------------------------------------#
#-------------------------------------- Main pipeline function ------------------------------------------# 
#--------------------------------------------------------------------------------------------------------#

def run_biascorrect(args):

	"""
	Function for bias correction of input .bw files
	Calls functions in BiasCorrect_functions and several internal classes
	"""

	#Test if required arguments were given:
	if args.bw == None:
		sys.exit("Error: No .bw-file given")
	if args.genome == None:
		sys.exit("Error: No .fasta-file given")

	#Adjust some parameters depending on input
	args.prefix = os.path.splitext(os.path.basename(args.bw))[0] if args.prefix == None else args.prefix
	args.outdir = os.path.abspath(args.outdir) if args.outdir != None else os.path.abspath(os.getcwd())

	#Set output bigwigs based on input
	tracks = ["bias", "corrected_lc", "corrected_gl", "expected_lc", "expected_gl"]
	tracks = [track for track in tracks if track not in args.track_off] 	# switch off printing

	if args.split_strands == True:
		strands = ["forward", "reverse"]
	else:
		strands = ["both"]

	output_bws = {}
	for track in tracks:
		output_bws[track] = {}
		for strand in strands:
			elements = [args.prefix, track] if strand == "both" else [args.prefix, track, strand]
			output_bws[track][strand] = {"fn": os.path.join(args.outdir, "{0}.bw".format("_".join(elements)))}

	#Set all output files
	bigwigs = [output_bws[track][strand]["fn"] for (track, strand) in itertools.product(tracks, strands)]
	figures_f = os.path.join(args.outdir, "{0}_biascorrect.pdf".format(args.prefix))
	
	output_files = bigwigs + [figures_f]
	output_files = list(OrderedDict.fromkeys(output_files)) 	#remove duplicates due to "both" option

	strands = ["forward", "reverse"]

	#----------------------------------------------------------------------------------------------------#
	# Print info on run
	#----------------------------------------------------------------------------------------------------#

	logger = foottrackLogger("BiasCorrect", args.verbosity)
	logger.begin()

	parser = add_biascorrect_arguments(argparse.ArgumentParser())
	logger.arguments_overview(parser, args)
	logger.output_files(output_files)

	args.cores = check_cores(args.cores, logger)

	#----------------------------------------------------------------------------------------------------#
	# Test input file availability for reading 
	#----------------------------------------------------------------------------------------------------#

	logger.info("----- Processing input data -----")

	logger.debug("Testing input file availability")
	check_files([args.bw, args.genome, args.peaks], "r")

	logger.debug("Testing output directory/file writeability")
	make_directory(args.outdir)
	check_files(output_files, "w")

	#Open pdf for figures
	figure_pdf = PdfPages(figures_f, keep_empty=False)

	#----------------------------------------------------------------------------------------------------#
	# Cr information in bw/fasta
	#----------------------------------------------------------------------------------------------------#

	logger.info("Reading info from .bw file")
	bwfile = pyBigWig.open(args.bw)
	bw_chrom_info = {chrom: bwfile.chroms(chrom) for chrom in bwfile.chroms()}
	bw_references = list(bwfile.chroms().keys())
	logger.debug("bw_chrom_info: {0}".format(bw_chrom_info))

	logger.info("Reading info from .fasta file")
	fastafile = pysam.FastaFile(args.genome)
	fasta_chrom_info = dict(zip(fastafile.references, fastafile.lengths))
	logger.debug("fasta_chrom_info: {0}".format(fasta_chrom_info))
	fastafile.close()

	#Compare chrom lengths
	chrom_in_common = set(bw_chrom_info.keys()).intersection(fasta_chrom_info.keys())
	for chrom in chrom_in_common:
		bwlen = bw_chrom_info[chrom]
		fastalen = fasta_chrom_info[chrom]
		if bwlen != fastalen:
			logger.warning("(Fastafile)\t{0} has length {1}".format(chrom, fasta_chrom_info[chrom]))
			logger.warning("(Bamfile)\t{0} has length {1}".format(chrom, bw_chrom_info[chrom]))
			sys.exit("Error: .bw and .fasta have different chromosome lengths. Please make sure the genome file is similar to the one used in mapping.")

	#Subset bw_references to those for which there are sequences in fasta
	chrom_not_in_fasta = set(bw_references) - set(fasta_chrom_info.keys())
	if len(chrom_not_in_fasta) > 1:
		logger.warning("The following contigs in --bw did not have sequences in --fasta: {0}. NOTE: These contigs will be skipped in calculation and output.".format(chrom_not_in_fasta))

	bw_references = [ref for ref in bw_references if ref in fasta_chrom_info]
	chrom_in_common = [ref for ref in chrom_in_common if ref in bw_references]

	#Drop mitochrondrial (or other chroms) from list
	chrom_in_common_orig = chrom_in_common
	chrom_in_common = [chrom for chrom in chrom_in_common if chrom not in args.drop_chroms]
	bw_references = [chrom for chrom in bw_references if chrom in chrom_in_common]

	dropped = set(chrom_in_common_orig) - set(chrom_in_common)
	if len(dropped) > 0:
		logger.info("The following contigs were dropped from analysis because they were found in '--drop-chroms': {0}".format(list(dropped)))
	else:
		logger.warning("No additional chromosomes were removed. Consider using '--drop-chroms' to remove mitochondrial and/or other unwanted contigs.")

	#Check if any contigs were left; else exit
	if len(chrom_in_common) == 0:
		logger.error("No common contigs left to run BiasCorrect on. Please check that '--bw' and '--fasta' are matching.")
		sys.exit()

	#----------------------------------------------------------------------------------------------------#
	# Cr regions from bedfiles
	#----------------------------------------------------------------------------------------------------#

	logger.info("Processing input/output regions")

	#Chromosomes included in analysis
	genome_regions = RegionList().from_list([OneRegion([chrom, 0, bw_chrom_info[chrom]]) for chrom in chrom_in_common]) #full genome length
	logger.debug("CHROMS\t{0}".format("; ".join(["{0} ({1})".format(reg.chrom, reg.end) for reg in genome_regions])))
	genome_bp = sum([region.get_length() for region in genome_regions])

	# calculate general conversion rate
	chrom_lengths = bwfile.chroms()
	signal_sum = 0
	total_length = 0
	for chrom in chrom_lengths:
		intervals = bwfile.intervals(chrom)
		for interval in intervals:
			start, end, value = interval
			length = end - start
			signal_sum += value * length
			total_length += length
	golobal_cr = signal_sum / total_length

	# Process peaks
	if args.peaks != None:
		peak_regions = RegionList().from_bed(args.peaks)
	else:
		chroms = bwfile.chroms().keys()
		chrom_effective_ranges = []
		for chrom in chroms:
			intervals = bwfile.intervals(chrom)
			if intervals:
				start_positions = [interval[0] for interval in intervals]
				end_positions = [interval[1] for interval in intervals]
				effective_start = min(start_positions)
				effective_end = max(end_positions)
				chrom_effective_ranges.append([chrom, effective_start, effective_end])
			else:
				chrom_effective_ranges.append([chrom, None, None])
		output_bed_file = os.path.join(args.outdir, args.prefix + "_effective_ranges.bed")
		with open(output_bed_file, 'w') as bed_file:
			for chrom, start, end in chrom_effective_ranges:
				if start is not None and end is not None:
					bed_file.write(f"{chrom}\t{start}\t{end}\n")
		peak_regions = RegionList().from_bed(output_bed_file)
		bwfile.close()

	#### Statistics about regions ####
	blacklist_regions = RegionList([])

	regions_dict = {"genome": genome_regions, "peak_regions": peak_regions}
	sub="peak_regions"
	regions_sub = regions_dict[sub]
	regions_sub.subtract(blacklist_regions)
	regions_sub = regions_sub.apply_method(OneRegion.split_region, 50000)
	regions_sub.keep_chroms(chrom_in_common)
	regions_dict[sub] = regions_sub

	genome_bp = sum([region.get_length() for region in regions_dict["genome"]])
	for key in regions_dict:
		total_bp = sum([region.get_length() for region in regions_dict[key]])
		logger.stats("{0}: {1} regions | {2} bp | {3:.2f}% coverage".format(key, len(regions_dict[key]), total_bp, total_bp/genome_bp*100))

	#Estallish variables for regions to be used
	input_regions = regions_dict["peak_regions"]
	output_regions = regions_dict["peak_regions"]

	#----------------------------------------------------------------------------------------------------#
	# Estimate normalization factors
	#----------------------------------------------------------------------------------------------------#

	#Setup logger queue
	logger.debug("Setting up listener for log")
	logger.start_logger_queue()
	args.log_q = logger.queue

	#----------------------------------------------------------------------------------------------------#
	logger.comment("")
	logger.info("----- Estimating normalization factors -----")
	correct_factor = 1.0
	logger.stats("CORRECTION_FACTOR:\t{0:.5f}".format(correct_factor))

	#----------------------------------------------------------------------------------------------------#
	# Estimate sequence bias
	#----------------------------------------------------------------------------------------------------#

	logger.comment("")

	if args.bias_pkl is None:

		logger.info("Started estimation of sequence bias...")

		input_region_chunks = input_regions.chunks(args.split)	#split to 100 chunks (also decides the step of output)
		out_lst = run_parallel(bias_estimation, input_region_chunks, [args], args.cores, logger)	#Output is list of enzBias objects

		#Join objects
		estimated_bias = out_lst[0]		#initialize object with first output
		for output in out_lst[1:]:
			estimated_bias.join(output)		#bias object contains bias/background SequenceMatrix objects

		logger.debug("Bias estimated\tno_crs: {0}".format(estimated_bias.no_crs))

	else:
		logger.info("Loading sequence bias from '--bias-pkl' file...")
		estimated_bias = enzBias().from_pickle(args.bias_pkl)

	#----------------------------------------------------------------------------------------------------#
	# Join estimations from all chunks of regions
	#----------------------------------------------------------------------------------------------------#

	bias_obj = estimated_bias
	bias_obj.correction_factor = correct_factor

	### Bias motif ###
	logger.info("Finalizing bias motif for scoring")
	for strand in strands:
		bias_obj.bias[strand].prepare_mat()
		logger.debug("Saving pssm to figure pdf")
		fig = plot_pssm(bias_obj.bias[strand].pssm, "enzyme bias ({0} strand)".format(strand))
		figure_pdf.savefig(fig)

	#Write bias motif to pickle
	out_f = os.path.join(args.outdir, args.prefix + "_enzBias.pickle")
	logger.debug("Saving bias object to pickle ({0})".format(out_f))
	bias_obj.to_pickle(out_f)
	
	#----------------------------------------------------------------------------------------------------#
	# Correct cr bias and write to bigwig
	#----------------------------------------------------------------------------------------------------#

	logger.comment("")
	logger.info("----- Correcting crs from .bw within output regions -----")

	output_regions.loc_sort(bw_references)		#sort in order of references
	output_regions_chunks = output_regions.chunks(args.split)
	# no_tasks = float(len(output_regions_chunks))
	chunk_sizes = [len(chunk) for chunk in output_regions_chunks]
	logger.debug("All regions chunked: {0} ({1})".format(len(output_regions), chunk_sizes))

	### Create key-file linking for bigwigs 
	key2file = {}
	for track in output_bws:
		for strand in output_bws[track]:
			filename = output_bws[track][strand]["fn"]
			key = "{}:{}".format(track, strand)
			key2file[key] = filename

	#Start correction/write cores
	n_bigwig = len(key2file.values())
	writer_cores = min(n_bigwig, max(1,int(args.cores*0.1)))	#at most one core per bigwig or 10% of cores (or 1)
	worker_cores = max(1, args.cores - writer_cores) 				
	logger.debug("Worker cores: {0}".format(worker_cores))
	logger.debug("Writer cores: {0}".format(writer_cores))

	worker_pool = mp.Pool(processes=worker_cores)
	writer_pool = mp.Pool(processes=writer_cores)
	manager = mp.Manager()

	#Start bigwig file writers
	writer_tasks = []
	header = [(chrom, bw_chrom_info[chrom]) for chrom in bw_references]
	key_chunks = [list(key2file.keys())[i::writer_cores] for i in range(writer_cores)]
	qs_list = []
	qs = {}
	for chunk in key_chunks:
		logger.debug("Creating writer queue for {0}".format(chunk))

		q = manager.Queue()
		qs_list.append(q)

		files = [key2file[key] for key in chunk]
		writer_tasks.append(writer_pool.apply_async(bigwig_writer, args=(q, dict(zip(chunk, files)), header, output_regions, args)))	 #, callback = lambda x: finished.append(x) print("Writing time: {0}".format(x)))
		for key in chunk:
			qs[key] = q

	args.qs = qs
	writer_pool.close() #no more jobs applied to writer_pool

	#Start correction
	logger.debug("Starting correction")
	task_list = [worker_pool.apply_async(bias_correction, args=[chunk, args, bias_obj, golobal_cr]) for chunk in output_regions_chunks]

	worker_pool.close()
	monitor_progress(task_list, logger, "Correction progress:")	#does not exit until tasks in task_list finished
	results = [task.get() for task in task_list]

	#Get all results 
	pre_bias = results[0][0]	#initialize with first result
	post_bias = results[0][1]	#initialize with first result

	for result in results[1:]:
		pre_bias_chunk = result[0]
		post_bias_chunk = result[1]

		for direction in strands:
			pre_bias[direction].add_counts(pre_bias_chunk[direction])
			post_bias[direction].add_counts(post_bias_chunk[direction])

	#Stop all queues for writing
	logger.debug("Stop all queues by inserting None")
	for q in qs_list:
		q.put((None, None, None))

	#Fetch error codes from bigwig writers
	logger.debug("Fetching possible errors from bigwig_writer tasks")
	results = [task.get() for task in writer_tasks]	#blocks until writers are finished

	logger.debug("Joining bigwig_writer queues")
	
	qsum = sum([q.qsize() for q in qs_list])
	while qsum != 0:
		qsum = sum([q.qsize() for q in qs_list])
		logger.spam("- Queue sizes {0}".format([(key, qs[key].qsize()) for key in qs]))
		time.sleep(0.5)

	#Waits until all queues are closed
	writer_pool.join() 
	worker_pool.terminate()
	worker_pool.join()

	#Stop multiprocessing logger	
	logger.stop_logger_queue()

	#----------------------------------------------------------------------------------------------------#
	# Information and verification of corrected cr frequencies
	#----------------------------------------------------------------------------------------------------#		

	logger.comment("")
	logger.info("Verifying bias correction")

	#Calculating variance per base
	for strand in strands:

		pre_bias[strand].prepare_mat()
		post_bias[strand].prepare_mat()

		pre_var = np.mean(np.var(pre_bias[strand].bias_pwm, axis=1)[:4])   #mean of variance per nucleotide
		post_var = np.mean(np.var(post_bias[strand].bias_pwm, axis=1)[:4])
		logger.stats("BIAS\tpre-bias variance {0}:\t{1:.7f}".format(strand, pre_var))
		logger.stats("BIAS\tpost-bias variance {0}:\t{1:.7f}".format(strand, post_var))

		#Plot figure
		# # pssm
		# fig_title = "corrected enzyme bias \n({0} strand)".format(strand)
		# figure_pdf.savefig(plot_pssm(post_bias[strand].pssm, fig_title))
		# Nucleotide frequencies
		fig_title = "Nucleotide frequencies \n({0} strand)".format(strand)
		figure_pdf.savefig(plot_correction(pre_bias[strand].bias_pwm, post_bias[strand].bias_pwm, fig_title))
	
	#----------------------------------------------------------------------------------------------------#
	# Finish up
	#----------------------------------------------------------------------------------------------------#
	plt.close('all')
	figure_pdf.close()
	logger.stats("general conversion rate :\t{}".format(golobal_cr))
	logger.end()

#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_biascorrect_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	run_biascorrect(args)