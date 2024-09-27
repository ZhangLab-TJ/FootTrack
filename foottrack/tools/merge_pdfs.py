#!/usr/bin/env python

"""
Small utility for mering pdfs
"""

from PyPDF2 import PdfMerger, PdfReader
import argparse
import sys
import os

#Internal functions
from foottrack.parsers import add_mergepdf_arguments
from foottrack.utils.utilities import *

#--------------------------------------------------------------------------------------------------------#
def run_mergepdf(args):

	check_required(args, ["input", "output"])
	print("Number of input files: {0}".format(len(args.input)))

	#Preliminary checks
	print("Checking read/write status")
	check_files(args.input, action="r")
	check_files([args.output], action="w")

	#Join pdfs
	print("Starting to merge PDFs")
	merger = PdfMerger(strict=False)
	for pdf in args.input:
		if os.stat(pdf).st_size != 0:	#only join files containing plots
			merger.append(PdfReader(pdf))
	
	print("Writing merged file: {0}".format(args.output))
	merger.write(args.output)

	print("PDFs merged successfully!")


#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_mergepdf_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	run_mergepdf(args)
