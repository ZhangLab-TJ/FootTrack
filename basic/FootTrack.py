#!/usr/bin/env python

import sys
import argparse
from argparse import SUPPRESS
import textwrap
import importlib

#Import extra dependencies
try:
	import MOODS
except:
	sys.exit("ERROR: Package MOODS is not installed and is needed by FootTrack. You can install it via conda using:\n"
			"$ conda install moods -c bioconda\n\n"
			"Or directly from source:\n"
			"$ wget https://github.com/jhkorhonen/MOODS/releases/download/v1.9.3/MOODS-python-1.9.3.tar.gz\n"
			"$ tar xzvf MOODS-python-1.9.3.tar.gz\n"
			"$ cd  MOODS-python-1.9.3\n"
			"$ python setup.py install"
			)

#Import parsers from foottrack
from foottrack.parsers import *
from foottrack import __version__ as FootTrack_VERSION

#Ignore gimmemotifs plot warning
import warnings
import matplotlib
matplotlib_version = tuple([int(i) for i in matplotlib.__version__.split(".")])

try:
	if matplotlib_version < (3, 8):
		import matplotlib.cbook
		warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)
	else:
		warnings.filterwarnings("ignore", category=matplotlib.MatplotlibDeprecationWarning)  # after 3.8, cbook.mplDeprecation is deprecated
except Exception:
	pass  # error in filtering warnings; not critical


def main():

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	all_parser_info = {"Tools for footprinting analysis":
							{
							"BiasCorrect":{"help":"Correct reads with regards to enzyme sequence bias", "add_arguments": add_biascorrect_arguments, "function": "foottrack.tools.biascorrect.run_biascorrect"},
							"ScoreBigwig":{"help":"Calculate scores such as footprints from conversion rate", "add_arguments": add_scorebigwig_arguments, "function": "foottrack.tools.score_bigwig.run_scorebigwig", "replaces":"FootprintScores"},
							"BINDetect":{"help":"Detect TF binding from footprints and motifs", "add_arguments": add_bindetect_arguments, "function": "foottrack.tools.bindetect.run_bindetect"},
							}
						}

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#					

	parser = argparse.ArgumentParser("FootTrack", usage=SUPPRESS)
	parser._action_groups.pop()
	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=25, width=90)
	parser.description = textwrap.dedent('''
										 ______________________________________________________________________________
										|                                                                              |
										|                                ~ FootTrack ~                                 |
										|                   Transcription factor Occupancy prediction                  |
										|                       By Investigation of cFOOT-seq Signal                   |
										|______________________________________________________________________________|
									
										Usage: FootTrack <TOOLNAME> [arguments]

										''')

	subparsers = parser.add_subparsers(title=None, metavar="")
	
	#Add all tools to parser
	all_tool_parsers = {}				
	for group in all_parser_info:
		parser.description += group + ":\n"

		info = all_parser_info[group]
		for tool in info:
			parser.description += "   {0}{1}{2}\n".format(tool, info[tool].get("space", "\t\t"), info[tool]["help"])
			subparser = subparsers.add_parser(tool, usage=SUPPRESS)
			subparser = info[tool]["add_arguments"](subparser)
			subparser.set_defaults(module=info[tool]["function"])
			all_tool_parsers[tool.lower()] = subparser

			#Add version to subparser
			subparser.add_argument("--version", action='version', version=FootTrack_VERSION)
			subparser = add_underscore_options(subparser)

			#Add parser for old tool names
			if "replaces" in info[tool]:
				replace_tool = info[tool]["replaces"]
				subparser = subparsers.add_parser(replace_tool, usage=SUPPRESS)
				subparser = info[tool]["add_arguments"](subparser)
				subparser.set_defaults(module=info[tool]["function"])
				all_tool_parsers[replace_tool.lower()] = subparser
			
		parser.description += "\n"

	parser.description += "For help on each tool, please run: FootTrack <TOOLNAME> --help\n"
	
	#Add version number to upper FootTrack parser 
	parser.description += "For version number: FootTrack --version"
	parser.add_argument("--version", action='version', version=FootTrack_VERSION)

	#If no args, print help for top-level FootTrack
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	#if args are pointing to specific tool, and no other arguments are given, print help for this
	if sys.argv[1].lower() in all_tool_parsers and len(sys.argv) == 2:
		chosen_tool = sys.argv[1]
		if chosen_tool != "DownloadData":	#Downloaddata can be run without options
			all_tool_parsers[chosen_tool.lower()].print_help()
			sys.exit()
	
	args = parser.parse_args()

	#Depending on subparser chosen, load main script entry and run
	function_str = args.module
	mod_name, func_name = function_str.rsplit(".", 1)
	module = importlib.import_module(mod_name)	#load specific module
	func = getattr(module, func_name)
	args.func = func

	#Run specified function with arguments
	args.func(args)		

	
if __name__ == "__main__":
    main()