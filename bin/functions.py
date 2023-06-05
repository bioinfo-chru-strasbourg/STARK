"""
Generic functions to use with STARK

Author: Samuel Nicaise (2019)
"""

from __future__ import division
from __future__ import print_function

import argparse
import glob
import math
import os
import random
import re
import subprocess
import time

from os.path import join as osj

def assert_file_exists_and_is_readable(filePath):
	assert os.path.isfile(filePath) and os.access(filePath, os.R_OK), \
			"[ERROR] File "+filePath+" doesn't exist or isn't readable"

def find_any_samplesheet(runDir, fromResDir):
	"""
	Input:
		runDir: String ; path/to/run where you want to find a samplesheet
		fromResDir: boolean ; if True the run dir is organized as a STARK result dir,
							else as a STARK repository dir (i.e. patient dirs contains a STARK dir)
	Output:
		String with samplesheet path or "NO_SAMPLESHEET_FOUND"

	Method:
	1) look up recursively all files named SampleSheet.csv in the runDir
	2) check if file path follows an expected samplesheet name and location
		(the latter depends on if we're in a STARK result or repository dir,
		defined by the bool fromResDir)
	3) first correct file path is returned
	"""
	p = subprocess.Popen("find -L "+runDir+" -maxdepth 3 -name *SampleSheet.csv", stdout=subprocess.PIPE, shell=True)
	out = p.stdout.readlines()
	for ss in out:
		ss = ss.decode("utf-8").strip()
		if fromResDir:
			r = re.match(runDir.rstrip("/")+"/(.*)/(.*).SampleSheet.csv", ss)
		else:
			r = re.match(runDir.rstrip("/")+"/(.*)/STARK/(.*).SampleSheet.csv", ss)
		if r is None:
			continue
		elif r.group(1) == r.group(2): #checks if (.*) == (.*)
			return ss
	return "NO_SAMPLESHEET_FOUND"

def get_sample_list_from_samplesheet(samplesheetPath):
	"""
	Returns a python list containing all sample names in a samplesheet.
	"""
	assert samplesheetPath != "NO_SAMPLESHEET_FOUND", \
			"[ERROR] find_any_samplesheet() couldn't find any samplesheet. Check if the --fromResultDir argument is set correctly."
	assert_file_exists_and_is_readable(samplesheetPath)
	inDataTable = False
	sampleList = []
	with open(samplesheetPath, "r") as f:
			for l in f:
				if not inDataTable:
					if l.startswith("Sample_ID,"):
						inDataTable = True
				else:
					if "," in l:
						sampleList.append(l.strip().split(",")[0])
	#if there are spaces in samplesheet names, change them to "_" because that's what demultiplexing.sh will do
	#otherwise the fastq won't be found when looking in the DEM dir
	for i in range(len(sampleList)):
		if " " in sampleList[i]:
			sampleList[i] = sampleList[i].replace(" ", "_")
	return sampleList

def get_descriptions_from_samplesheet(samplesheetPath):
	"""
	Returns a python list containing all tags from samplesheet.
	"""
	assert samplesheetPath != "NO_SAMPLESHEET_FOUND", \
			"[ERROR] find_any_samplesheet() couldn't find any samplesheet. Check if the --fromResultDir argument is set correctly."
	assert_file_exists_and_is_readable(samplesheetPath)
	inDataTable = False
	tagList = []
	with open(samplesheetPath, "r") as f:
		for l in f:
			if not inDataTable:
				if l.startswith("Sample_ID,"):
					inDataTable = True
					descriptionIndex = l.strip().split(",").index("Description")
			else:
				if "," in l:
					tagList.append(l.strip().split(",")[descriptionIndex])
	return tagList


def tags_and_types_to_lists(tags, tagTypes="", associator="#", separator="!"):
	"""
	Input:
		tags: String containing the tags, usually a Description field from a Samplesheet
		tagTypes: String containing the space-separated tag types you want to get in the output.
				If left empty, return all the tags.
	Output:
		2 same length lists of Strings in a tuple:
		(tagValuesList, tagTypesList)
	"""
	tagValuesList = []
	tagTypesList = []
	if associator in tags:
		if separator in tags:
			if tagTypes == "":
				for tag in tags.split(separator):
					if tag == "":
						continue # in case someone put a '!' at the end of the description field

					#python split()'s second argument defines the maximum number of splits
					#so split('#', 1) splits only based on the first #
					#examples:
					# '#a#b'.split('#')[1] --> 'a'
					# '#a#b'.split('#', 1)[1] --> 'a#b'

					#so here, if tag = "test#v1#v2":
					# tagValuesList.append("v1 v2")
					# tagTypesList.append("test")
					if associator in tag:
						tagValuesList.append(tag.split(associator, 1)[1].replace(associator, " "))
						tagTypesList.append(tag.split(associator)[0])
		else:
			tagValuesList.append(tags.split(associator, 1)[1].replace(associator, " "))
			tagTypesList.append(tags.split(associator)[0])

	return (tagValuesList, tagTypesList)

def get_lockfile_path(lockfile_prefix, target):
	return lockfile_prefix + os.path.basename(target)

def is_rule_startable(lockfile_prefix, target, max_jobs, shell_mode=False):
	"""
	Input:
		run_dir: run directory, used because it is a common path accessible for all rules
		target: targeted file, used to create the lockfile suffix
		max_jobs: maximum number of instances of the rule that can be executed at once
	
	Output:
		print stating if the rule instance invoking this function can be allowed to start
			1 if yes
			0 if no

	Usage:
		Invoke from a make function with a format such as:
		while <target> does not exist:
			if limit_rule_concurrency == True:
				create target
			else:
				sleep 10

	Python functions can be invoked from make with a syntax such as:
	python -c "import sys; sys.path.append('/home1/bin/STARK/0.9.18.4.rnaseq/bin'); import functions; functions.limit_rule_concurrency(args)"

	"""
	if not os.path.exists(os.path.dirname(lockfile_prefix)):
		raise FileNotFoundError(os.path.dirname(lockfile_prefix))

	if len(glob.glob(osj(lockfile_prefix + "*"))) < int(max_jobs): #int cast in case of user error
		lock_file = get_lockfile_path(lockfile_prefix, target)
		open(lock_file, 'a').close()
		if shell_mode:
			print("1")
		return True
	else:
		if shell_mode:
			print("0")
		return False

def rule_finished(lockfile_prefix, target):
	"""
	Inputs: 
		See is_job_startable
	Output: 
		None, cleans run_dir of target lockfile so other jobs can start
	"""
	lock_file = get_lockfile_path(lockfile_prefix, target)
	if os.path.exists(lock_file):
		os.remove(lock_file)

def randomized_sleep(string):
	"""
	Multiple simultaneous executions of launch_when_possible() can result in a number of lockfiles/started processes over the defined maximum.
	To limit this issue, sleep for a random time based on target filename to spread out the executions.
	"""
	# ord gives the unicode point for each char in filename --> converts filename in a "random" number
	unicode_points_sum = sum([ord(x) for x in string])
	# we want to avoid that similar filenames lead to similar wait times
	# but wait times should stay small overall
	# the random function allows to spread wait times while staying between 1 and 15 seconds max
	random.seed(unicode_points_sum)
	wait_time = random.random() * 15 #the arbitrary *15 increases spread
	time.sleep(wait_time)


def launch_when_possible(cmd, lockfile_prefix, target, max_jobs, current_dir, timeout):
	"""
	current_dir allows to emulate launching commands from the user's current directory, even if this python script is elsewhere.
	timeout is in hours
	"""
	randomized_sleep(target)
	cmd = "cd " + current_dir + " && " + cmd
	starting_time = time.time()

	while True:
		try:
			if is_rule_startable(lockfile_prefix, target, max_jobs):
				print("Launching: ", cmd)
				subprocess.call(cmd, shell=True)
				rule_finished(lockfile_prefix, target)
				return
		finally:
			#remove lockfile if script is killed
			rule_finished(lockfile_prefix, target)

		if time.time() - starting_time > timeout * 3600:
			raise RuntimeError("Timeout: couldn't create target: " + target)
		time.sleep(10)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers(help="sub-command help")
	parser_launch = subparsers.add_parser(
		"launch", help="launch command, ensuring rule is not started more than X times at once"
	)
	parser_startable = subparsers.add_parser(
		"is_rule_startable", help="returns 1 if rule is startable, 0 otherwise"
	)
	parser_finished = subparsers.add_parser(
		"rule_finished", help="signals that rule is finished so others can start"
	)

	parser_launch.add_argument("-c", "--cmd", type=str, required=True)

	for p in (parser_launch, parser_startable, parser_finished):
		p.add_argument("-l", "--lockfile_prefix", type=str, required=True)
		p.add_argument("-t", "--target", type=str, required=True)
	
	for p in (parser_launch, parser_startable):
		p.add_argument("-m", "--max_jobs", type=int, required=True)
		p.add_argument("--shell_mode", action='store_true', default=False)

	parser_launch.add_argument("-d", "--current_dir", type=str, default=os.path.abspath("."))
	parser_launch.add_argument("-z", "--timeout", type=int, default=48)

	args = parser.parse_args()

	if args.max_jobs < 1:
		raise ValueError("max_jobs argument should be >= 1")

	if "max_jobs" in args:
		if "cmd" in args:
			launch_when_possible(args.cmd,
								args.lockfile_prefix,
								args.target,
								args.max_jobs,
								args.current_dir,
								args.timeout)
		else:
			is_rule_startable(args.lockfile_prefix,
							args.target,
							args.max_jobs,
							args.shell_mode)
	else:
		rule_finished(args.lockfile_prefix, args.lock_name)