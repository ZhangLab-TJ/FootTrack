#!/usr/bin/env python

import argparse
from foottrack.utils.utilities import format_help_description, restricted_float, add_underscore_options
from foottrack.utils.logger import add_logger_args

#--------------------------------------------------------------------------------------------------------#
def add_biascorrect_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=35, width=90)

	description = "BiasCorrect corrects the cutsite-signal from cFOOT-seq with regard to the underlying sequence preference of deaminase.\n\n"
	description += "Usage:\nFootTrack BiasCorrect --bw <reads.bw> --genome <genome.fa> --peaks <peaks.bed>\n\n"
	description += "Output files:\n"
	description += "\n".join(["- <outdir>/<prefix>_{0}.bw".format(track) for track in ["uncorrected", "bias", "expected", "corrected"]]) + "\n"
	description += "- <outdir>/<prefix>_biascorrect.pdf"
	parser.description = format_help_description("BiasCorrect", description)

	parser._action_groups.pop()	#pop -h

	#Required arguments
	reqargs = parser.add_argument_group('Required arguments')
	reqargs.add_argument('-b', '--bw', metavar="<bw>", help="A .bw-file containing reads to be corrected")
	reqargs.add_argument('-g', '--genome', metavar="<fasta>", help="A .fasta-file containing whole genomic sequence")
	reqargs.add_argument('-p', '--peaks', metavar="<bed>", help="A .bed-file containing analysis peak regions")

	#Optional arguments
	optargs = parser.add_argument_group('Optional arguments')
	optargs.add_argument('--extend', metavar="<int>", type=int, help="Extend output regions with basepairs upstream/downstream (default: 100)", default=100)
	optargs.add_argument('--split-strands', help="Write out tracks per strand", action="store_true")
	optargs.add_argument('--norm-off', help="Switches off normalization based on number of reads", action='store_true')
	optargs.add_argument('--track-off', metavar="<track>", help="Switch off writing of individual .bigwig-tracks (bias/corrected_lc/corrected_gl/expected_lc/expected_gl) lc=local background, gl=global background", nargs="*", default=[])
	optargs.add_argument('--drop-chroms', metavar="<chrom>", help="Drop any chromosomes in the list from the correction. The default is to drop the mitochrondrial chromosome. Default: ['chrM', 'chrMT', 'M', 'MT', 'Mito']", nargs="*", default=['chrM', 'chrMT', 'M', "MT", "Mito"])

	optargs = parser.add_argument_group('Advanced BiasCorrect arguments (no need to touch)')
	optargs.add_argument('--k_flank', metavar="<int>", help="Flank +/- of cutsite to estimate bias from (default: 12)", type=int, default=12)
	optargs.add_argument('--bg_shift', metavar="<int>", type=int, help="Read shift for estimation of background frequencies (default: 100)", default=100)
	optargs.add_argument('--window', metavar="<int>", help="Window size for calculating expected signal (default: 100)", type=int, default=100)
	# optargs.add_argument('--score_mat', metavar="<mat>", help="Type of matrix to use for bias estimation (PWM/DWM) (default: DWM)", choices=["PWM", "DWM"], default="DWM")
	optargs.add_argument('--score_mat', metavar="<mat>", help="Type of matrix to use for bias estimation (PWM) (default: PWM)", choices=["PWM"], default="PWM")
	optargs.add_argument('--bias-pkl', metavar="<obj>", help="Path to a pre-calculated enzBias.pkl-object, as output from a previous BiasCorrect run (default: None). Can be used to bypass the internal bias estimation.", default=None)

	runargs = parser.add_argument_group('Run arguments')
	runargs.add_argument('--prefix', metavar="<prefix>", help="Prefix for output files (default: same as .bw file)")
	runargs.add_argument('--outdir', metavar="<directory>", help="Output directory for files (default: current working directory)", default="")
	runargs.add_argument('--cores', metavar="<int>", type=int, help="Number of cores to use for computation (default: 1)", default=1)
	runargs.add_argument('--split', metavar="<int>", type=int, help="Split of multiprocessing jobs (default: 100)", default=100)
	
	runargs = add_logger_args(runargs)

	return(parser)

#--------------------------------------------------------------------------------------------------------#
def add_scorebigwig_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "ScoreBigwig calculates scores (such as footprint-scores) from bigwig files (such as cFOOT-seq conversion rate calculated using the BiasCorrect tool). " 
	description += "NOTE: ScoreBigwig replaces the previous FootprintScores command.\n\n"
	description += "Usage: ScoreBigwig --signal <conversion rate.bw> --regions <regions.bed> --output <output.bw>\n\n"
	description += "Output:\n- <output.bw>"
	parser.description = format_help_description("ScoreBigwig", description)
	
	parser._action_groups.pop()	#pop -h

	#Required arguments
	required = parser.add_argument_group('Required arguments')
	required.add_argument('-s', '--signal', metavar="<bigwig>", help="A .bw file of cFOOT-seq cutsite signal")
	required.add_argument('-o', '--output', metavar="<bigwig>", help="Full path to output bigwig")			
	required.add_argument('-r', '--regions', metavar="<bed>", help="Genomic regions to run footprinting within")

	optargs = parser.add_argument_group('Optional arguments')
	optargs.add_argument('--score', metavar="<score>", choices=["footprint", "sum", "mean", "none"], help="Type of scoring to perform on conversion rate (footprint/sum/mean/none) (default: footprint)", default="footprint")
	optargs.add_argument('--absolute', action='store_true', help="Convert bigwig signal to absolute values before calculating score")
	optargs.add_argument('--extend', metavar="<int>", type=int, help="Extend input regions with bp (default: 100)", default=100)
	optargs.add_argument('--smooth', metavar="<int>", type=int, help="Smooth output signal by mean in <bp> windows (default: no smoothing)", default=1)
	optargs.add_argument('--min-limit', metavar="<float>", type=float, help="Limit input bigwig value range (default: no lower limit)") 		#default none
	optargs.add_argument('--max-limit', metavar="<float>", type=float, help="Limit input bigwig value range (default: no upper limit)") 		#default none

	footprintargs = parser.add_argument_group('Parameters for score == footprint')
	footprintargs.add_argument('--fp-min', metavar="<int>", type=int, help="Minimum footprint width (default: 10)", default=10)
	footprintargs.add_argument('--fp-max', metavar="<int>", type=int, help="Maximum footprint width (default: 10)", default=10)
	footprintargs.add_argument('--flank-min', metavar="<int>", type=int, help="Minimum range of flanking regions (default: 20)", default=20)
	footprintargs.add_argument('--flank-max', metavar="<int>", type=int, help="Maximum range of flanking regions (default: 20)", default=20)
	
	sumargs = parser.add_argument_group('Parameters for score == sum')
	sumargs.add_argument('--window', metavar="<int>", type=int, help="The window for calculation of sum (default: 100)", default=100)

	runargs = parser.add_argument_group('Run arguments')
	runargs.add_argument('--cores', metavar="<int>", type=int, help="Number of cores to use for computation (default: 1)", default=1)
	runargs.add_argument('--split', metavar="<int>", type=int, help="Split of multiprocessing jobs (default: 100)", default=100)
	runargs = add_logger_args(runargs)

	return(parser)

#--------------------------------------------------------------------------------------------------------#
def add_bindetect_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=35, width=90)
	description = "BINDetect takes motifs, signals (footprints) and genome as input to estimate bound transcription factor binding sites and differential binding between conditions. "
	description += "The underlying method is a modified motif enrichment test to see which motifs have the largest differences in signal across input conditions. "
	description += "The output is an in-depth overview of global changes as well as the individual binding site signal-differences.\n\n"
	description += "Usage:\nFootTrack BINDetect --signals <bigwig1> (<bigwig2> (...)) --motifs <motifs.txt> --genome <genome.fasta> --peaks <peaks.bed>\n\n"
	description += "Output files:\n- <outdir>/<prefix>_figures.pdf\n- <outdir>/<prefix>_results.{txt,xlsx}\n- <outdir>/<prefix>_distances.txt\n"
	description += "- <outdir>/<TF>/<TF>_overview.{txt,xlsx} (per motif)\n- <outdir>/<TF>/beds/<TF>_all.bed (per motif)\n"
	description += "- <outdir>/<TF>/beds/<TF>_<condition>_bound.bed (per motif-condition pair)\n- <outdir>/<TF>/beds/<TF>_<condition>_unbound.bed (per motif-condition pair)\n\n"
	parser.description = format_help_description("BINDetect", description)

	parser._action_groups.pop()	#pop -h
	
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--signals', metavar="<bigwig>", help="Signal per condition (.bigwig format)", nargs="*")
	required.add_argument('--peaks', metavar="<bed>", help="Peaks.bed containing open chromatin regions across all conditions")
	required.add_argument('--motifs', metavar="<motifs>", help="Motif file(s) in pfm/jaspar/meme format", nargs="*")
	required.add_argument('--genome', metavar="<fasta>", help="Genome .fasta file")

	optargs = parser.add_argument_group('Optional arguments')
	optargs.add_argument('--cond-names', metavar="<name>", nargs="*", help="Names of conditions fitting to --signals (default: prefix of --signals)")
	optargs.add_argument('--peak-header', metavar="<file>", help="File containing the header of --peaks separated by whitespace or newlines (default: peak columns are named \"_additional_<count>\")")
	optargs.add_argument('--naming', metavar="<string>", help="Naming convention for TF output files ('id', 'name', 'name_id', 'id_name') (default: 'name_id')", choices=["id", "name", "name_id", "id_name"], default="name_id")
	optargs.add_argument('--motif-pvalue', metavar="<float>", type=lambda x: restricted_float(x, 0, 1), help="Set p-value threshold for motif scanning (default: 1e-4)", default=0.0001)
	optargs.add_argument('--bound-pvalue', metavar="<float>", type=lambda x: restricted_float(x, 0, 1), help="Set p-value threshold for bound/unbound split (default: 0.001)", default=0.001)
	#optargs.add_argument('--volcano-diff-thresh', metavar="<float>", help="", default=0.2)	#not yet implemented
	#optargs.add_argument('--volcano-p-thresh', metavar="<float>", help="", default=0.05)	#not yet implemented

	optargs.add_argument('--pseudo', type=float, metavar="<float>", help="Pseudocount for calculating log2fcs (default: estimated from data)", default=None)
	optargs.add_argument('--time-series', action='store_true', help="Will only compare signals1<->signals2<->signals3 (...) in order of input, and skip all-against-all comparison.")
	optargs.add_argument('--skip-excel', action='store_true', help="Skip creation of excel files - for large datasets, this will speed up BINDetect considerably")
	optargs.add_argument('--output-peaks', metavar="<bed>", help="""Gives the possibility to set the output peak set differently than the input --peaks.
													 				This will limit all analysis to the regions in --output-peaks. 
																	NOTE: --peaks must still be set to the full peak set!""")
	optargs.add_argument('--norm-off', action='store_true', help="Turn off normalization of footprint scores across conditions")

	runargs = parser.add_argument_group("Run arguments")
	runargs.add_argument('--outdir', metavar="<directory>", help="Output directory to place TFBS/plots in (default: bindetect_output)", default="bindetect_output")
	optargs.add_argument('--prefix', metavar="<prefix>", help="Prefix for overview files in --outdir folder (default: bindetect)", default="bindetect")
	runargs.add_argument('--cores', metavar="<int>", type=int, help="Number of cores to use for computation (default: 1)", default=1)
	runargs.add_argument('--split', metavar="<int>", type=int, help="Split of multiprocessing jobs (default: 100)", default=100)
	runargs.add_argument('--debug', action='store_true', help="Creates an additional '_debug.pdf'-file with debug plots")	#creates extra output for debugging
	
	runargs = add_logger_args(runargs)

	return(parser)