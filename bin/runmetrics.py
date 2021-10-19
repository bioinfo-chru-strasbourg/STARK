"""
@GOAL
Goal of this script is to create a run metrics tsv for one run. It will be launched by STARK every run.
It will also create secondary files with gene/design/amplicon specific metrics.

######## STARK 0.9.18 documentation
Generates up to 4 tsv files, for each design and panel.
Mainly used to regroup content from <sample>.bwamem.metrics files
- global_coverage.tsv : summarized stats with one line for each sample in the run
- coverage.tsv : compiled content of all HsMetrics.per_target_coverage.flags files
- depth.tsv : compiled content of all genes.genes.txt files
- amplicon_coverage.tsv : compiled content from all .per_amplicon_coverage.flags files


########## STARK 0.9.17 documentation (no longer relevant for STARK 0.9.18, useful to understand the structure of this script)
@INPUTS AND OUTPUTS
Metrics to fetch according to julien/nadege/jean/samuel meeting (27/05/2019):
for each sample, both total and %:
- total reads <-- in bwamem.bam.flagstat
- mapped reads <-- in bwamem.bam.flagstat
- duplicate reads <-- in bwamem.bam.flagstat
- on target reads (ratio based on non-duplicate reads) <-- (*)
- and % depth at 5, 10, 20, 30X <-- computed from depthbed

+ specify what bed is used

(*): 2 values exist in stark report:
	- %bases on target in hsmetrics file ; not what was discussed
	- %reads on target made with the filtered bam: view -F 1024 -F 4 -q 10
		The reads on target are computed with samtools view cleaned.bam -L sample.bed
		and the results are found in bwamem.bam.metrics/<sample>.bwamem.<on/off>.nbreads
The latter is exactly what was discussed with the biologists; except for the -F4 -q 10 flags.
I'm taking this value for now, but this should be discussed.

All stark files listed above exist in Stark 0.9.17 and below outputs. The script is fitted to use the new 0.9.18 files, but the old data file names are commented out for an easy rollback.

The secondary cov metrics files created here will be copies of their equivalent in the report, but grouped by sample in a single file.

@LAUNCH EXAMPLE
[New main_routine] Made to be launched by Stark using "metrics" file to identify metrics directories and data files. Example featuring a full test run and tags (data files on bioinfo@INT):
python /home1/TOOLS/tools/toolkit/current/lib/sam/forStark/runmetrics.py --metricsFileList `echo /home1/datas/WORK_DIR_SAM/TOOLSstark18RES/results/VALIDATION_RUN_TEST_HUSTUMSOL_tagsArgs/P1335_R1/P1335_R1.bwamem.bam.metrics/metrics /home1/datas/WORK_DIR_SAM/TOOLSstark18RES/results/VALIDATION_RUN_TEST_HUSTUMSOL_tagsArgs/P1338_R1/P1338_R1.bwamem.bam.metrics/metrics /home1/datas/WORK_DIR_SAM/TOOLSstark18RES/results/VALIDATION_RUN_TEST_HUSTUMSOL_tagsArgs/P1376_R1/P1376_R1.bwamem.bam.metrics/metrics /home1/datas/WORK_DIR_SAM/TOOLSstark18RES/results/VALIDATION_RUN_TEST_HUSTUMSOL_tagsArgs/P1406_R1/P1406_R1.bwamem.bam.metrics/metrics /home1/datas/WORK_DIR_SAM/TOOLSstark18RES/results/VALIDATION_RUN_TEST_HUSTUMSOL_tagsArgs/P1408_R1/P1408_R1.bwamem.bam.metrics/metrics /home1/datas/WORK_DIR_SAM/TOOLSstark18RES/results/VALIDATION_RUN_TEST_HUSTUMSOL_tagsArgs/P1431B_R1/P1431B_R1.bwamem.bam.metrics/metrics /home1/datas/WORK_DIR_SAM/TOOLSstark18RES/results/VALIDATION_RUN_TEST_HUSTUMSOL_tagsArgs/P1439_R1/P1439_R1.bwamem.bam.metrics/metrics | tr " " ","` --outputPrefix /home1/datas/WORK_DIR_SAM/TOOLSstark18RES/results/VALIDATION_RUN_TEST_HUSTUMSOL_tagsArgs/analysis.V20190724-150053.report.

@TODO
DONE: the legends of coverage files could be improved.
DONE: colonne nom du run
TODO: stats globales pour une maladie/appli [beyond the scope of this script]
DONE: variables environnement pour legende (MAPQ...) si applicable
DONE: gestion tags (CTUM...)
DONE: nom des fichiers (avec date, et changer run_metrics en reads.metrics)
DONE: how do I deal with so many .genes???
- reads.metrics: add columns for every .gene (or use only the design values <-- I did that)
- one more design/amplicon type file for every .gene file
--> adds only one more file per .gene, and should be doable without rewriting everything. Find the .genes list common to all samples by looking up each <sample>/DATA/<sample>.list.gene ; and do one metrics file with each common .genes file
############################

@Author: Samuel Nicaise (29/05/2019)
"""

from __future__ import division
from __future__ import print_function

import argparse
import os
import re
import subprocess
from os.path import join as osj

from functions import assert_file_exists_and_is_readable, tags_and_types_to_lists, \
					get_descriptions_from_samplesheet, find_any_samplesheet

def get_run_path_from_metrics_file(metricsFile, fromResDir):
	"""
	Example: "/home1/datas/WORK_DIR_SAM/TOOLSstark18RES/results/test_bam_bed_metrics/TEST/TEST.bwamem.bam.metrics/metrics"
	should return "/home1/datas/WORK_DIR_SAM/TOOLSstark18RES/results/test_bam_bed_metrics/"
	"""
	if fromResDir:
		r = re.match("(.*)/(.*)/(.*)/metrics", metricsFile)
	else:
		r = re.match("(.*)/(.*)/STARK/(.*)/metrics", metricsFile)
	assert r.group(3).startswith(r.group(2)) #way to check if the pattern fits a patient name
	return r.group(1)

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

def get_total_mapped_duplicate_reads(sampleDir, sample, aligner):
	"""
	Fetches total reads, mapped reads (total and %), duplicate reads (total and %) from the flagstat file in sampleDir.
	Returns them as a list in that order.
	"""
	flagstatFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+".flagstat")
	assert_file_exists_and_is_readable(flagstatFile)
	with open(flagstatFile, "r") as f:
		col1 = [int(line.split()[0]) for line in f]
	if col1[0] != 0:
		return [col1[0], col1[4], format(col1[4]/col1[0]*100, ".2f"), col1[3], format(col1[3]/col1[0]*100, ".2f")]
	else:
		return [col1[0], col1[4], "NA", col1[3], "NA"]

def get_on_target_reads(sampleDir, sample, aligner, bed):
	"""
	Fetches on target reads (total and %). Returns them as a list.
	More explanations in the docstring at the top of this script.
	"""
	if bed == "design":
		#for stark before 0.9.18d
		#onReadsFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+".on.nbreads")
		# onReadsFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+"."+sample+".genes_from_manifest.on.nbreads")
		# onReadsFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+"."+sample+"."+aligner+".design.bed.on.nbreads")
		onReadsFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+"."+sample+"."+aligner+".design.bed.on.target")
		#for stark before 0.9.18d
		#offReadsFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+".off.nbreads")
		# offReadsFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+"."+sample+".genes_from_manifest.off.nbreads")
		offReadsFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+"."+sample+"."+aligner+".design.bed.off.target")
	else: #bed is a .genes file
		onReadsFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+"."+sample+"."+bed+".on.target")
		offReadsFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+"."+sample+"."+bed+".off.target")

	with open(onReadsFile, "r") as fon:
		on = int(fon.readline().rstrip())
	with open(offReadsFile, "r") as foff:
		off = int(foff.readline().rstrip())
	if (on+off) != 0:
		return [on, format(on*100/(on+off), ".2f")]
	else:
		return [on, "NA"]

def get_cov_criteria():
	if "COVERAGE_CRITERIA" in os.environ.keys():
		covCriteria = os.environ["COVERAGE_CRITERIA"]
	else:
		covCriteria = "1,5,10,20,30,100"
	maxCriteria = str(max([int(v) for v in covCriteria.split(",")]))
	return covCriteria, maxCriteria

def deprecated_get_depth_metrics(sampleDir, sample, aligner):
	"""
	Computes % depth at 1X,5X... (defined by an envt variable or a default value) based on depthbed. Returns them as a list in that order.
	Uses the awk magic formula from stark_report.sh to keep the same values as the pdf report.

	DEPRECATED as the depthbed file no longer exists in newer stark 0.9.18 versions.
	"""
	#for stark before 0.9.18d
	#depthFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+".depthbed")
	# depthFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+"."+sample+".genes_from_manifest.depthbed")
	depthFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+"."+sample+"."+aligner+".design.bed.depthbed")
	assert_file_exists_and_is_readable(depthFile)
	covCriteria, maxCriteria = get_cov_criteria()
	cmd = """cat %s | awk -v MDP=%s  '{SUM++} { if ($3>MDP) {DP[MDP]++} else {DP[$3]++} } END { for (i=MDP; i>=0; i-=1) {print i" "DP[i]" SUM"SUM}}' | sort -g -r | awk -v COVS=%s '{SUM+=$2} {CUM[$1]=SUM} {split(COVS,C,",")}  END  { for (j in C) {print C[j]" "CUM[C[j]]" SUM "(CUM[C[j]]/SUM)} }' | sort -g""" % (depthFile, maxCriteria, covCriteria)
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
	out = p.stdout.readlines()
	covMetrics = []
	for l in out:
		l = l.decode("utf-8").strip()
		covMetrics.append(format(float(l.split()[3])*100, ".2f"))
	return covMetrics

def get_depth_metrics(sampleDir, sample, aligner, bed):
	"""
	coverage can now be found in files such as TEST.bwamem.TEST.from_design.genes.coverage
	whose content is like:
	#Depth  CoveredBases    TotalBases      Percent
	1X      274     11413   0.0240077
	5X      274     11413   0.0240077
	10X     274     11413   0.0240077
	20X     274     11413   0.0240077
	30X     274     11413   0.0240077
	50X     254     11413   0.0222553
	100X    0       11413   0
	200X    0       11413   0
	300X    0       11413   0
	This function extracts the % values and returns them as a list.
	"""
	if bed == "design":
		depthFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+"."+sample+"."+aligner+".design.bed.coverage")
	else: #.genes file
		depthFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+"."+sample+"."+bed+".coverage")
	covCriteria, maxCriteria = get_cov_criteria()
	if "," in covCriteria:
		covCriteriaList = covCriteria.split(",")
	else:
		covCriteriaList = [covCriteria]
	covMetrics = []
	with open(depthFile, "r") as f:
		for l in f:
			l = l.strip().split()
			if l[0][:-1] in covCriteriaList:
				covMetrics.append(format(float(l[3])*100, ".5f"))
	assert len(covMetrics) == len(covCriteriaList), "[ERROR] Content of the environment variable COVERAGE_CRITERIA did not fit the content of the coverage file "+depthFile
	return covMetrics

def get_sample_metrics(runPath, sample, fromResDir, aligner, bed):
	"""
	Fetches all the metrics defined in the docstring at the top of this script for <run>.reads.metrics.
	"""
	if fromResDir:
		sampleDir = osj(runPath, sample)
	else:
		sampleDir = osj(runPath, sample, "STARK")
	metricsList = get_total_mapped_duplicate_reads(sampleDir, sample, aligner)
	metricsList += get_on_target_reads(sampleDir, sample, aligner, bed)
	metricsList += get_depth_metrics(sampleDir, sample, aligner, bed)
	return metricsList

def get_tags_list(sampleList, sampleDirList, samplesheet=""):
	"""
	Returns a list of tags string in the order corresponding to sampleList

	Tries to fetch tags through STARK analysis and sample tags files.
	If it fails and a samplesheet is provided, fetches them from there instead.
	"""
	try:
		tagsList = []
		for sample, sampleDir in zip(sampleList, sampleDirList):
			tags = ""
			with open(osj(sampleDir, sample+".analysis.tag"), "r") as f:
				for l in f:
					tags += l.strip()
			with open(osj(sampleDir, sample+".tag"), "r") as f:
				for l in f:
					if tags != "":
						tags += "!"+l.strip()
					else:
						tags = l.strip()
			tagsList.append(tags)
		# print(tagsList)
		return tagsList
		# return ["APP#SWAG!#TUMCELL" for v in sampleList]
	except IOError as e:
		print(e)
		print("[WARNING] Couldn't use <sample>.tag and <sample>.analysis.tag to determine tags, trying to use samplesheet instead.")
		assert samplesheet != "", "[ERROR] Tried to use samplesheet to determine tags, but no samplesheet was provided."
		if os.path.getsize(samplesheet) == 0:
			#then no tags to be found
			return ["" for v in sampleList]
			# return ["APP#SWAG!#TUMCELL" for v in sampleList]
		else:
			#get both samples and tags to make sure the tags returned
			#are in the same order as the sampleList provided in input
			samples = get_sample_list_from_samplesheet(samplesheet)
			tags = get_descriptions_from_samplesheet(samplesheet)
			return [tags[samples.index(s)] for s in sampleList]

def build_cov_metrics_header(dataFileList, tagsList):
	"""Builds the header for write_cov_metrics_file()"""
	#Base
	header = "#Run\tSample"
	
	# #DEPRECATED - used to allow to have a column per tag
	# #appending tag types
	# allTagTypes = []
	# for tags in tagsList:
		# values, types = tags_and_types_to_lists(tags)
		# allTagTypes += types
	# seen = set()
	# #remove duplicate types for the header, while still keeping them ordered
	# allTagTypes = [t for t in allTagTypes if not(t in seen or seen.add(t))]

	# untypedValuesExist = False
	# for t in allTagTypes:
		# if t == '':
			# untypedValuesExist = True
		# else:
			# header += "\t" + t
	# if untypedValuesExist:
		# header += "\t" + "Other tags"
	if any([len(t) > 0 for t in tagsList]):
		header += "\t" + "Tags"

	#appending data file header
	with open(dataFileList[0], "r") as fin:
		headerTmp = fin.readline().strip()
	if headerTmp.startswith("#"):
		headerTmp = headerTmp[1:]
	header += "\t" + headerTmp
	header += "\n"

	# return header, allTagTypes
	return header

def write_cov_metrics_file(name, run, sampleList, dataFileList, legend, tagList):
	"""
	Input:
	- name: metrics file name (including path)
	- sampleList: list of sample base names
	- dataFileList : python list of the files to concatenate
	- legend: text at the start of the file commented with ##. Usually a multiline docstring.
	- tags: list of tag strings (should be the same length as sample)
	"""
	with open(name, "w") as fout:
		fout.write(legend)

		#Create header with the first file
		# header, headerTags = build_cov_metrics_header(dataFileList, tagList)
		header = build_cov_metrics_header(dataFileList, tagList)
		fout.write(header)

		for tags, sample, dataFile in zip(tagList, sampleList, dataFileList):
			tagString = ""
			# values, types = tags_and_types_to_lists(tags)
			# for h in headerTags:
				# if h in types:
					# tagString += "\t" + values[types.index(h)]
				# else:
					# tagString += "\t" + ""
			if len(tags) > 0:
				tagString = "\t" + tags
			with open(dataFile, "r") as fin:
				next(fin) #skips the header
				for l in fin:
					fout.write(os.path.basename(run)+"\t"+sample+tagString+"\t"+l)

def str_to_bool(v):
	"""
	Returns a boolean based on the String v in input. Used for argument parsing
	"""
	if v.lower() in ('yes', 'true', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('One arg is expected to be a String convertible to a boolean value.')

def reverse_dict_of_lists(dic):
	"""
	>>> reverse_dict_of_lists({'s1': ['g1'] ; 's2' : ['g1', 'g2'})
	{'s1': ['g1', 'g2'] ; 's2': ['g1']}
	"""
	r = {}
	for k in dic.keys():
		for v in dic[k]:
			if v in r.keys():
				r[v].append(k)
			else:
				r[v] = [k]
	return r

def make_global_coverage_file(run, regionsType, regionsFileList, metricsFileList, sampleList ,alignerList, genesSampleDic, outputPrefix, fromResDir):
	"""
	Generates the files containing samtools flagstat data and mean Cov 1X/5X/...X over all samples in a run

	design.<nomdesign>.global_coverage.metrics
	or
	panel.<nompanel>.global_coverage.metrics
	depending on input files

	where regionsFileList is either [design_name] or [name1.genes, name2.genes, ...]
	and regionsType indicates if it's a design or a panel

	sampleGenesDic indicates if a given regionFile is associated to a given sample.
	"""
	assert regionsType in ("Design", "Panel")
	assert len(regionsFileList) > 0

	#create matrix
	print("Getting " + regionsType+ " global_coverage data...")
	runMetrics = {}
	for regionsFile in regionsFileList:
		for metricsFile, sample, aligner in zip(metricsFileList, sampleList ,alignerList):
			# runMetrics[sample]=[] #faster execution to debug everything else
			#distinction needed here because design metrics files are named "design" in STARK, while .genes metrics used the full .genes name
			if regionsType == "Design":
				runMetrics[sample] = get_sample_metrics(run, sample, fromResDir, aligner, "design")
				runMetrics[sample].insert(0, regionsFile)
			elif regionsType == "Panel":
				if sample in genesSampleDic[regionsFile]:
					runMetrics[sample] = get_sample_metrics(run, sample, fromResDir, aligner, regionsFile)
					runMetrics[sample].insert(0, regionsFile)

		#write matrix
		if regionsType == "Design":
			finalTsv = osj(outputPrefix + regionsType + ".global_coverage.tsv")
		elif regionsType == "Panel":
			finalTsv = osj(outputPrefix + regionsType + "." + regionsFile + ".global_coverage.tsv")
		print("Writing " + finalTsv)
		with open(finalTsv, "w") as f:
			covHeader = "\t".join(["Cov "+v+"X" for v in get_cov_criteria()[0].split(",")])
			f.write("## Run Metrics\n"
				"##\n"
				"## On-target reads are reads that are aligned within the regions defined in " + regionsType + " and that aren't duplicates, unmapped or with low mapping quality (MAPQ < 10)\n"
				"## Cov 30X is the % of coverage with at least 30X read depth ; in the regions defined in " + regionsType +".\n"
				"##\n"
				"#Run\tSample\t" + regionsType + "\tTotal reads\tMapped reads\t% Mapped reads\tDuplicate reads\t% Duplicate reads\tOn-target reads\t% On-target reads\t"+covHeader+"\n")
			#Could be converted to "," for French Excel visualization with str(v).replace(".", ",") ; triggered by a variable from a STARK app

			if regionsType == "Design":
				for sample in sampleList:
					f.write(os.path.basename(run)+"\t"+sample+"\t"+"\t".join([str(v) for v in runMetrics[sample]])+"\n")
			else:
				for sample in genesSampleDic[regionsFile]:
					f.write(os.path.basename(run)+"\t"+sample+"\t"+"\t".join([str(v) for v in runMetrics[sample]])+"\n")

def make_depth_file(run, regionsType, regionsFileList, metricsFileList, sampleList, sampleDirList, alignerList, tagsList, outputPrefix):
	"""
	Previously named "metrics.design.metrics"
	Made to generate:
	- design.<nomdesign>.coverage.metrics
	- panel.<nompanel>.coverage.metrics
	"""
	assert regionsType in ("Design", "Panel")
	assert len(regionsFileList) > 0

	for regionsFile in regionsFileList:
		finalSampleList = []
		finalTagsList = []
		dataFileList = []

		#create matrix
		for metricsFile, sample, sampleDir, aligner, tags in zip(metricsFileList, sampleList, sampleDirList, alignerList, tagsList):
			#before latest 0.9.18 version
			# dataFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+".HsMetrics.per_target_coverage")
			if regionsType == "Design":
				dataFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+"."+sample+"."+aligner+".design.bed.HsMetrics.per_target_coverage.flags")
				dataFileList.append(dataFile)
				finalSampleList.append(sample)
				finalTagsList.append(tags)
			elif regionsType == "Panel":
				dataFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+"."+sample+"."+regionsFile+".HsMetrics.per_target_coverage.flags")
				if os.path.exists(dataFile):
					print("iamhere")
					dataFileList.append(dataFile)
					finalSampleList.append(sample)
					finalTagsList.append(tags)

		#write matrix
		if regionsType == "Design":
			finalTsv = osj(outputPrefix + regionsType + ".depth.tsv")
		elif regionsType == "Panel":
			finalTsv = osj(outputPrefix + regionsType + "." + regionsFile + ".depth.tsv")
		print("Writing " + finalTsv)
		legend=("## Depth metrics - per region in " + regionsType + "\n"
				"##\n"
				"## Compilation of sample statistics computed by STARK.\n"
				"##\n")
		write_cov_metrics_file(finalTsv, run, finalSampleList, dataFileList, legend, finalTagsList)

def make_coverage_file(run, regionsType, regionsFileList, metricsFileList, sampleList, sampleDirList, alignerList, tagsList, outputPrefix):
	"""
	fichier metrics.sample.<nompanel>.metrics
	- ajouter design.<nomdesign>.depth.metrics
	- renommer existant en panel.<nompanel>.depth.metrics
	"""
	assert regionsType in ("Design", "Panel")
	assert len(regionsFileList) > 0

	for regionsFile in regionsFileList:
		finalSampleList = []
		finalTagsList = []
		dataFileList = []
		legend=("## Coverage metrics - - per region in " + regionsType + " " + regionsFile +"\n"
		"##\n"
		"## Compilation of sample statistics computed by STARK.\n"
		"##\n")

		for metricsFile, sample, sampleDir, aligner, tags in zip(metricsFileList, sampleList, sampleDirList, alignerList, tagsList):
			if regionsType == "Design":
				dataFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+"."+sample+"."+aligner+".design.bed.genes.txt")
				dataFileList.append(dataFile)
				finalSampleList.append(sample)
				finalTagsList.append(tags)
			else:
				dataFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+"."+sample+"."+regionsFile+".genes.txt")
				if os.path.exists(dataFile):
					dataFileList.append(dataFile)
					finalSampleList.append(sample)
					finalTagsList.append(tags)

		if regionsType == "Design":
			finalTsv = osj(outputPrefix + regionsType + ".coverage.tsv")
		elif regionsType == "Panel":
			finalTsv = osj(outputPrefix + regionsType + "." + regionsFile + ".coverage.tsv")
		print("Writing " + finalTsv)
		write_cov_metrics_file(finalTsv, run, sampleList, dataFileList, legend, tagsList)

def make_amplicon_file(run, metricsFileList, sampleList, sampleDirList, alignerList, tagsList, outputPrefix):
	"""
	relevant files only exist with design
	"""
	dataFileList = []
	for metricsFile, sample, sampleDir, aligner in zip(metricsFileList, sampleList, sampleDirList, alignerList):
		dataFile = osj(sampleDir, sample+"."+aligner+".bam.metrics", sample+"."+aligner+".HsMetrics.per_amplicon_coverage.flags")
		assert_file_exists_and_is_readable(dataFile)
		dataFileList.append(dataFile)
	finalTsv = osj(outputPrefix+"Design.amplicon_coverage.tsv")
	legend=("## Coverage metrics - per amplicon\n"
		"##\n"
		"## The regions listed here correspond to the amplicons listed in the manifest.\n"
		"##\n")
	print("Writing " + finalTsv)
	write_cov_metrics_file(finalTsv, run, sampleList, dataFileList, legend, tagsList)

def main_routine(metricsFileList, outputPrefix, fromResDir):
	"""
	metricsFileList is a list of <run>/.../<sample>.<aligner>.bam.metrics/metrics files ; it defines the sample list and input directories
	outputPrefix is the output prefix of all files that will be written
	fromResDir is a boolean defining if input files are from a STARK result dir (True) or from a STARK repository (False)
	"""
	#extract parameters
	#####
	print("fromResDir:", fromResDir)
	if "," in metricsFileList:
		#applying os.path.normpath now so I don't have to bother with // in paths later
		metricsFileList = [os.path.normpath(v) for v in metricsFileList.split(',')]
	else:
		metricsFileList = [os.path.normpath(metricsFileList)]
	run = get_run_path_from_metrics_file(metricsFileList[0], fromResDir)
	sampleList = [v.replace(run+"/", "").split("/")[0] for v in metricsFileList]
	#to get the aligner: need to get e.g. "bwamem" from a path like
	#/home1/datas/WORK_DIR_SAM/TOOLSstark18RES/results/test_bam_bed_metrics/TEST/TEST.bwamem.bam.metrics/metrics
	alignerList = [os.path.dirname(mF).split("/")[-1].split(".")[1] for mF in metricsFileList]
	if fromResDir:
		sampleDirList = [osj(run, s) for s in sampleList]
	else:
		sampleDirList = [osj(run, s, "STARK") for s in sampleList]
	tagsList = get_tags_list(sampleList, sampleDirList, samplesheet=find_any_samplesheet(run, fromResDir))

	#find .genes (Panel) files
	#####
	sampleGenesDic = {}
	genesSetList = []
	for metricsFile in metricsFileList:
		genesSet = set()
		sample = metricsFile.replace(run+"/", "").split("/")[0]
		genesListFile = osj(os.path.dirname(os.path.dirname(metricsFile)), sample+".list.genes")
		with open(genesListFile, "r") as f:
			for l in f:
				#removing sample name if present in .genes name
				genesSet.add(l.strip().replace(sample+".", ""))
		genesSetList.append(genesSet)
		sampleGenesDic[sample] = list(genesSet)
	###Deprecated: intersect all sets to find the .genes files that are common to all samples
	#genesList = list(set.intersection(*genesSetList))
	#print(".genes files common to all samples:", genesList)
	genesList = list(set.union(*genesSetList))
	#reverse dict
	genesSampleDic = reverse_dict_of_lists(sampleGenesDic)


	#design files
	#TODO: deal with design names in the case there are multiples
	make_global_coverage_file(run, "Design", ["Design"], metricsFileList, sampleList, alignerList, genesSampleDic, outputPrefix, fromResDir)
	make_depth_file(run, "Design", ["Design"], metricsFileList, sampleList, sampleDirList, alignerList, tagsList, outputPrefix)
	make_coverage_file(run, "Design", ["Design"], metricsFileList, sampleList, sampleDirList, alignerList, tagsList, outputPrefix)
	#panels files
	make_global_coverage_file(run, "Panel", genesList, metricsFileList, sampleList, alignerList, genesSampleDic, outputPrefix, fromResDir)
	make_depth_file(run, "Panel", genesList, metricsFileList, sampleList, sampleDirList, alignerList, tagsList, outputPrefix)
	make_coverage_file(run, "Panel", genesList, metricsFileList, sampleList, sampleDirList, alignerList, tagsList, outputPrefix)


	#amplicon coverage metrics
	#####
	# Check if the first amplicon metrics file is not empty...
	testAligner = alignerList[0]
	skipAmpliconMetrics = False
	if fromResDir:
		testFile = osj(run, sampleList[0], sampleList[0]+"."+testAligner+".bam.metrics", sampleList[0]+"."+testAligner+".HsMetrics.per_amplicon_coverage.flags")
	else:
		testFile = osj(run, sampleList[0], "STARK", sampleList[0]+"."+testAligner+".bam.metrics", sampleList[0]+"."+testAligner+".HsMetrics.per_amplicon_coverage.flags")
	with open(testFile, "r") as f:
		if f.readline().rstrip() == "":
			skipAmpliconMetrics = True
	# ... If it isn't, go on to create the metrics file
	if not skipAmpliconMetrics:
		make_amplicon_file(run, metricsFileList, sampleList, sampleDirList, alignerList, tagsList, outputPrefix)

	print("runmetrics.py done")

if __name__=="__main__":
	parser = argparse.ArgumentParser(prog='python runmetrics.py')
	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument("-m", "--metricsFileList", type=str, help="path/to/metrics,path/to/metrics")
	parser.add_argument("-f", "--fromResultDir", type=str_to_bool, default=True, help="specifies if the directories in the folder are organised like in STARK results (default) or like in STARK repository (then enter False/No/N/0)")
	parser.add_argument("-o", "--outputPrefix", type=str, default="", help="output prefix (including dir). Files go into run dir by default if not specified.")
	args = parser.parse_args()

	main_routine(args.metricsFileList, args.outputPrefix, args.fromResultDir)
