#!/usr/bin/env python

"""
Class for dealing with logger-function used across all FootTrack tools

"""

import sys
import os
from datetime import datetime
import logging
import logging.handlers
import multiprocessing as mp
import time

#-------------------------------------------------------------------------------------------#

def add_logger_args(args):
	""" Function for adding FootTrack-wide verbosity to command-line parsers """
	args.add_argument('--verbosity', metavar="<int>", help="Level of output logging (0: silent, 1: errors/warnings, 2: info, 3: stats, 4: debug, 5: spam) (default: 3)", choices=[0,1,2,3,4,5], default=3, type=int)
	return(args)

class foottrackLogger(logging.Logger):
	""" foottrackLogger is an instance of a logging.Logger with special functions for formatting and creating automatic logging """

	logger_levels = {
						0: 0,
						1: logging.WARNING,							#also includes errors
						2: logging.INFO, 							#info 
						3: int((logging.INFO+logging.DEBUG)/2),		#statistics
						4: logging.DEBUG,							#debugging info
						5: logging.DEBUG - 5						#spam-level debug
					}

	def __init__(self, tool_name="FootTrack", level=3, queue=None):

		self.tool_name = tool_name		#name of tool within FootTrack
		logging.Logger.__init__(self, self.tool_name)

		if level == 0:
			self.disabled = True

		####### Setup custom levels #######
		#Create custom level for comments (Same level as errors/warnings)
		comment_level = foottrackLogger.logger_levels[1] + 1
		logging.addLevelName(comment_level, "comment") #log_levels[lvl])
		setattr(self, 'comment', lambda *args: self.log(comment_level, *args))

		#Create custom level for stats (between info and debug)
		stats_level = foottrackLogger.logger_levels[3]
		logging.addLevelName(stats_level, "STATS") #log_levels[lvl])
		setattr(self, 'stats', lambda *args: self.log(stats_level, *args))
		
		#Create custom level for spamming debug messages
		spam_level = foottrackLogger.logger_levels[5]
		logging.addLevelName(spam_level, "SPAM") #log_levels[lvl])
		setattr(self, 'spam', lambda *args: self.log(spam_level, *args))

		#Set level
		self.level = foottrackLogger.logger_levels[level]
		
		######### Setup formatter ########

		#Setup formatting
		self.formatter = FootTrackFormatter()
		self.setLevel(self.level)
		
		########## Setup streaming #########
		##Log file stream
		#if log_f != None:
		#	log = logging.FileHandler(log_f, "w")
		#	log.setLevel(self.level)
		#	log.setFormatter(self.formatter)
		#	self.addHandler(log)

		if queue == None:
			#Stdout stream
			con = logging.StreamHandler(sys.stdout)		#console output
			con.setLevel(self.level)
			con.setFormatter(self.formatter)
			self.addHandler(con)
		else:
			h = logging.handlers.QueueHandler(queue)  	# Just the one handler needed
			self.handlers = []
			self.addHandler(h)

		#Lastly, initialize time
		self.begin_time = datetime.now()
		self.end_time = None
		self.total_time = None


	def begin(self):
		""" Begin logging by writing comments about the current run """
		from foottrack import __version__ as FootTrack_VERSION

		self.cmd = "FootTrack " + " ".join(sys.argv[1:])

		#Print info on run
		self.comment("# FootTrack {0} {1} (run started {2})".format(FootTrack_VERSION, self.tool_name, self.begin_time))
		self.comment("# Working directory: {0}".format(os.getcwd()))
		self.comment("# Command line call: {0}\n".format(self.cmd))

	def stop(self):
		""" Stop without printing status """
		
		self.end_time = datetime.now()
		self.total_time = self.end_time - self.begin_time
		
	def end(self):
		""" End logging - write out how long it took """

		self.end_time = datetime.now()
		self.total_time = self.end_time - self.begin_time

		self.comment("")	#creates empty line; only for pretty output
		self.info("Finished {0} run (total time elapsed: {1})".format(self.tool_name, self.total_time))


	def start_logger_queue(self):
		""" start process for listening and handling through the main logger queue """

		self.debug("Starting logger queue for multiprocessing")
		self.queue = mp.Manager().Queue()
		self.listener = mp.Process(target=self.main_logger_process)
		self.listener.start()
		

	def stop_logger_queue(self):
		""" Stop process for listening """

		self.debug("Waiting for listener to finish")
		self.queue.put(None)
		while self.listener.exitcode != 0:
			self.debug("Listener exitcode is: {0}. Waiting for exitcode = 0.".format(self.listener.exitcode))
			time.sleep(0.1)

		self.debug("Joining listener")
		self.listener.join()


	def main_logger_process(self):

		self.debug("Started main logger process")
		while True:
			try:
				record = self.queue.get()
				if record is None:
					break
				self.handle(record) 	#this logger is coming from the main process

			except EOFError:
				self.error("Multiprocessing logger lost connection to queue - probably due to an error raised from a child process.")
				break

		return(1)

	def arguments_overview(self, parser, args):
		""" Creates string of arguments and options to print using logger """

		content = ""
		content += "# ----- Input parameters -----\n"
		for group in parser._action_groups:

			group_actions = group._group_actions
			#print(args)
			if len(group_actions) > 0:
				#content += "# ----- {0} -----\n".format(group.title)
				for option in group_actions:
					if option.help != "==SUPPRESS==": #only show if not suppressed
						name = option.dest
						attr = getattr(args, name, None)
						content += "# {0}:\t{1}\n".format(name, attr)
				#content += "\n"
		self.comment(content + "\n")

	def output_files(self, outfiles):
		""" Print out list of output files"""

		self.comment("# ----- Output files -----")
		for outf in outfiles:
			if outf != None:
				self.comment("# {0}".format(outf))
		self.comment("\n")
		
class FootTrackFormatter(logging.Formatter):
	""" Formatter class used in foottrackLogger """
	default_fmt = logging.Formatter("%(asctime)s (%(process)d) [%(levelname)s]\t%(message)s", "%Y-%m-%d %H:%M:%S")
	comment_fmt = logging.Formatter("%(message)s")

	def format(self, record):

		#Comments
		if record.levelname == "comment":
			return self.comment_fmt.format(record)
		elif record.levelno != 0:
			return self.default_fmt.format(record)	
		else:
			return
