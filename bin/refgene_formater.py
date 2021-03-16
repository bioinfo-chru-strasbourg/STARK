# -*- coding: utf-8 -*-

"""
Goal: automatically download and reformat RefGene all_fields file into the format that we commonly use with STARK. 

Link for more info about the columns of the all_fields file: 
https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=847920841_aCwaNipMVAu126E4Wu7dY5naakdP&clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=refSeqComposite&hgta_table=refGene&hgta_regionType=genome&position=chrX%3A15%2C578%2C261-15%2C621%2C068&hgta_outputType=bed&hgta_outFileName=RefGene_full.bed

Refgene all_fields:
chr1    48999844        48999965        NM_032785_cds_0_0_chr1_48999845_r       0       -

Output file:
For Stark17:
chr1    14408   14409   NR_046018,exon3 DDX11L1 +
For Stark18:
chr1    67000041    67000051    SGIP1    NM_001308203    +    exon1
"""

from __future__ import division
from __future__ import print_function

import argparse
import os
import subprocess
import time

def load_config_dic():
	"""
	Could do anything, like loading a json file or catching environment variables from a STARK application. 
	Right now, it's hardcoded url for the UCSC website. 
	"""
	config = {}
	config["UCSC_Refgene_url"] = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"
	config["NCBI_RefSeqGene_url"] = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ncbiRefSeq.txt.gz"
	return config

def get_stdout(cmd, debug=False):
	"""
	input is a string of a shell cmd
	returns a list of strings, where every string is a line from the stdout.
	"""
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	out = p.stdout.readlines()
	if debug:
		print(out)
	outList = [l.rstrip() for l in out]
	err = p.stderr.readlines()
	if debug:
		print(err)
	errList = [l.rstrip() for l in err]
	# assert len(errList) == 0, print("ERROR: len(errList)>0. cmd=", cmd, "errList=", errList)
	# wget uses stderr for its normal operations, so I can't do this check that I usually do. 
	return outList

def reformat_stark17(refgene_file, output, include_utr_5, include_utr_3, include_chrM, include_non_canonical_chr, include_non_coding_transcripts):
	"""
	Core of this module. Takes refgene_file as input, and writes the output file. 
	Other options are boolean determining parameters for file creation. 
	"""
	with open (refgene_file, "r") as fi: 
		with open(output, "w") as fo:
			for l in fi:
				if l.startswith("#"):
					continue
				l = l.strip().split()
				
				#exclude non-coding transcripts
				if not include_non_coding_transcripts:
					if not l[1].startswith("NM_"):
						continue
				
				#exclude unwanted chr
				if include_chrM and include_non_canonical_chr:
					pass
				elif include_chrM and not include_non_canonical_chr:
					if l[2] not in {"chr1": None, "chr10": None, "chr11": None, "chr12": None, "chr13": None, "chr14": None, "chr15": None, "chr16": None, "chr17": None, "chr18": None, "chr19": None, "chr2": None, "chr20": None, "chr21": None, "chr22": None, "chr3": None, "chr4": None, "chr5": None, "chr6": None, "chr7": None, "chr8": None, "chr9": None, "chrX": None, "chrY": None, "chrM": None, "chrMT": None}:
						continue
				elif not include_chrM and not include_non_canonical_chr:
					if l[2] not in {"chr1": None, "chr10": None, "chr11": None, "chr12": None, "chr13": None, "chr14": None, "chr15": None, "chr16": None, "chr17": None, "chr18": None, "chr19": None, "chr2": None, "chr20": None, "chr21": None, "chr22": None, "chr3": None, "chr4": None, "chr5": None, "chr6": None, "chr7": None, "chr8": None, "chr9": None, "chrX": None, "chrY": None}:
						continue
				
				# packing exons together to be able to number them (need to start from the end if negative strand, if using IGV as standard)
				exonStarts = l[9][:-1].split(",")
				exonEnds = l[10][:-1].split(",")
				exonsToWrite = []
				
				#deal with utr
				# l[6] = cds start ; l[7] = cds end. To remove UTR on both ends, consider only exons between those columns and truncate them if needed.
				if include_utr_5 and include_utr_3:
					for i in range(len(exonStarts)):
						exonsToWrite.append([l[2], exonStarts[i], exonEnds[i], l[1], l[12], l[3]])
				
				elif not include_utr_5 and include_utr_3:
					for i in range(len(exonStarts)):
						if int(exonStarts[i]) < int(l[6]):
							exonsToWrite.append([l[2], l[6], exonEnds[i], l[1], l[12], l[3]])
						else:
							exonsToWrite.append([l[2], exonStarts[i], exonEnds[i], l[1], l[12], l[3]])
				
				elif include_utr_5 and not include_utr_3:
					for i in range(len(exonStarts)):
						if int(exonEnds[i]) > int(l[7]):
							exonsToWrite.append([l[2], exonStarts[i], l[7], l[1], l[12], l[3]])
						else:
							exonsToWrite.append([l[2], exonStarts[i], exonEnds[i], l[1], l[12], l[3]])
				
				elif not include_utr_5 and not include_utr_3:
					for i in range(len(exonStarts)):
						if int(exonStarts[i]) >= int(l[6]) and int(exonEnds[i]) <= int(l[7]):
							exonsToWrite.append([l[2], exonStarts[i], exonEnds[i], l[1], l[12], l[3]])
						elif int(exonStarts[i]) < int(l[6]) and int(l[6]) <= int(exonEnds[i]) <= int(l[7]):
							exonsToWrite.append([l[2], l[6], exonEnds[i], l[1], l[12], l[3]])
						elif int(l[6])<= int(exonStarts[i]) <= int(l[7]) and int(exonEnds[i]) > int(l[7]):
							exonsToWrite.append([l[2], exonStarts[i], l[7], l[1], l[12], l[3]])
						elif int(exonStarts[i]) <= int(l[6]) <= int(exonEnds[i]) and int(exonStarts[i]) <= int(l[7]) <= int(exonEnds[i]):
							exonsToWrite.append([l[2], l[6], l[7], l[1], l[12], l[3]])
				
				if l[3] == "-":
					k = len(exonsToWrite)
					for i in range(len(exonsToWrite)):
						exonsToWrite[i][3] = exonsToWrite[i][3] + ",exon"+str(k)
						k-=1
						fo.write("\t".join(exonsToWrite[i])+"\n")
				else:
					k = 1
					for i in range(len(exonsToWrite)):
						exonsToWrite[i][3] = exonsToWrite[i][3] + ",exon"+str(k)
						k+=1
						fo.write("\t".join(exonsToWrite[i])+"\n")

def reformat_stark18(refgene_file, output, include_utr_5, include_utr_3, include_chrM, include_non_canonical_chr, include_non_coding_transcripts):
	"""
	Core of this module. Takes refgene_file as input, and writes the output file. 
	Other options are boolean determining parameters for file creation. 
	
	Difference with stark17: slightly different columns.
	"""
	with open (refgene_file, "r") as fi: 
		with open(output, "w") as fo:
			for l in fi:
				if l.startswith("#"):
					continue
				l = l.strip().split()
				
				#exclude non-coding transcripts
				if not include_non_coding_transcripts:
					if not l[1].startswith("NM_"):
						continue
				
				#exclude unwanted chr
				if include_chrM and include_non_canonical_chr:
					pass
				elif include_chrM and not include_non_canonical_chr:
					if l[2] not in {"chr1": None, "chr10": None, "chr11": None, "chr12": None, "chr13": None, "chr14": None, "chr15": None, "chr16": None, "chr17": None, "chr18": None, "chr19": None, "chr2": None, "chr20": None, "chr21": None, "chr22": None, "chr3": None, "chr4": None, "chr5": None, "chr6": None, "chr7": None, "chr8": None, "chr9": None, "chrX": None, "chrY": None, "chrM": None, "chrMT": None}:
						continue
				elif not include_chrM and not include_non_canonical_chr:
					if l[2] not in {"chr1": None, "chr10": None, "chr11": None, "chr12": None, "chr13": None, "chr14": None, "chr15": None, "chr16": None, "chr17": None, "chr18": None, "chr19": None, "chr2": None, "chr20": None, "chr21": None, "chr22": None, "chr3": None, "chr4": None, "chr5": None, "chr6": None, "chr7": None, "chr8": None, "chr9": None, "chrX": None, "chrY": None}:
						continue
				
				# packing exons together to be able to number them (need to start from the end if negative strand, if using IGV as standard)
				exonStarts = l[9][:-1].split(",")
				exonEnds = l[10][:-1].split(",")
				exonsToWrite = []
				
				#deal with utr
				# l[6] = cds start ; l[7] = cds end. To remove UTR on both ends, consider only exons between those columns and truncate them if needed.
				if include_utr_5 and include_utr_3:
					for i in range(len(exonStarts)):
						exonsToWrite.append([l[2], exonStarts[i], exonEnds[i], l[12], l[1], l[3]])
				
				elif not include_utr_5 and include_utr_3:
					for i in range(len(exonStarts)):
						if int(exonStarts[i]) < int(l[6]):
							exonsToWrite.append([l[2], l[6], exonEnds[i], l[12], l[1], l[3]])
						else:
							exonsToWrite.append([l[2], exonStarts[i], exonEnds[i], l[12], l[1], l[3]])
				
				elif include_utr_5 and not include_utr_3:
					for i in range(len(exonStarts)):
						if int(exonEnds[i]) > int(l[7]):
							exonsToWrite.append([l[2], exonStarts[i], l[7], l[12], l[1], l[3]])
						else:
							exonsToWrite.append([l[2], exonStarts[i], exonEnds[i], l[12], l[1], l[3]])
				
				elif not include_utr_5 and not include_utr_3:
					for i in range(len(exonStarts)):
						if int(exonStarts[i]) >= int(l[6]) and int(exonEnds[i]) <= int(l[7]):
							exonsToWrite.append([l[2], exonStarts[i], exonEnds[i], l[12], l[1], l[3]])
						elif int(exonStarts[i]) < int(l[6]) and int(l[6]) <= int(exonEnds[i]) <= int(l[7]):
							exonsToWrite.append([l[2], l[6], exonEnds[i], l[12], l[1], l[3]])
						elif int(l[6])<= int(exonStarts[i]) <= int(l[7]) and int(exonEnds[i]) > int(l[7]):
							exonsToWrite.append([l[2], exonStarts[i], l[7], l[12], l[1], l[3]])
						elif int(exonStarts[i]) <= int(l[6]) <= int(exonEnds[i]) and int(exonStarts[i]) <= int(l[7]) <= int(exonEnds[i]):
							exonsToWrite.append([l[2], l[6], l[7], l[12], l[1], l[3]])
				
				if l[3] == "-":
					k = len(exonsToWrite)
					for i in range(len(exonsToWrite)):
						exonsToWrite[i].append("exon"+str(k))
						k-=1
						fo.write("\t".join(exonsToWrite[i])+"\n")
				else:
					k = 1
					for i in range(len(exonsToWrite)):
						exonsToWrite[i].append("exon"+str(k))
						k+=1
						fo.write("\t".join(exonsToWrite[i])+"\n")

def sort_output(output):
	tmp_name = output+"_unsorted_"+str(time.time())
	os.rename(output, tmp_name)
	#http://azaleasays.com/2014/02/21/sort-chromosome-names-with-linux-sort/
	# get_stdout("sort -k1,1 -V -s " + tmp_name + " > " + output)
	get_stdout("sort -k1,1V -k2,2n " + tmp_name + " > " + output) #Mateusz
	os.remove(tmp_name)

def main_download_and_format(args):
	refgene_file = args.refgene_file
	if args.ncbi:
		url = load_config_dic()["NCBI_RefSeqGene_url"]
	else:
		url = load_config_dic()["UCSC_Refgene_url"]
	refgene_tmpname = os.path.dirname(os.path.normpath(args.output)) + os.path.basename(url)
	assert refgene_tmpname.endswith(".gz")
	cmd = "wget --directory-prefix=" + os.path.dirname(os.path.normpath(args.output)) + " " + url
	get_stdout(cmd)
	if refgene_file == "":
		#use default name from UCSC
		refgene_file = refgene_tmpname[:-3]
		if os.path.exists(refgene_tmpname[:-3]):
			print("[WARNING] No refgene_file specified: will overwrite "+refgene_file)
		get_stdout("gunzip -f -c "+refgene_tmpname+" > "+refgene_file)
	else:
		get_stdout("gunzip -f -c "+refgene_tmpname+" > "+refgene_file)
	get_stdout("rm -f "+refgene_tmpname)
	
	reformat_stark18(refgene_file, args.output, args.include_utr_5, args.include_utr_3, args.include_chrM, args.include_non_canonical_chr, args.include_non_coding_transcripts)
	sort_output(args.output)
	
	print("[INFO] refgene_formater.py done.")
	
def main_format_only(args):
	reformat_stark18(args.refgene_file, args.output, args.include_utr_5, args.include_utr_3, args.include_chrM, args.include_non_canonical_chr, args.include_non_coding_transcripts)
	sort_output(args.output)
	print("[INFO] refgene_formater.py done.")

if __name__=="__main__":
	parser = argparse.ArgumentParser(prog='python refgene_formater.py')
	subparsers = parser.add_subparsers(help='')
	
	parser_a = subparsers.add_parser("downloadAndFormat", help="Download a refgene file and format it")
	parser_a.add_argument("-r", "--refgene_file", type=str, default="", help="How the downloaded refgene file will be renamed. If left empty, will keep the default name from UCSC.")
	parser_a.add_argument("-o", "--output", type=str, default="formated_refgene.tsv", help="Output file. The refgene file will be downloaded to the same directory. Default: ./formated_refgene.tsv")
	parser_a.add_argument("-n", "--ncbi", default=False, action='store_true', help="By default, the UCSC Refgene is used. With this flag, the NCBI version of RefSeqGene is used instead. The main difference is that in the UCSC version, transcripts were realigned and have multiple possible positions. Transcripts have an unique position in the NCBI version.")
	parser_a.add_argument("-u5", "--include_utr_5", default=False, action='store_true')
	parser_a.add_argument("-u3", "--include_utr_3", default=False, action='store_true')
	parser_a.add_argument("-m", "--include_chrM", default=False, action='store_true')
	parser_a.add_argument("-ncc", "--include_non_canonical_chr", default=False, action='store_true', help="Adapted (tested) only with human hg19 chr")
	parser_a.add_argument("-nct", "--include_non_coding_transcripts", default=False, action='store_true', help="Activating this will filter out  all transcripts that are not NM_....")
	parser_a.set_defaults(mode=main_download_and_format)
	
	parser_b = subparsers.add_parser("format", help="Format a refgene file given as input")
	required = parser_b.add_argument_group('required arguments')
	required.add_argument("-r", "--refgene_file", type=str, default="", help="Name of the refGene file that will be used as input for the formatting.", required=True)
	parser_b.add_argument("-o", "--output", type=str, default="formated_refgene.tsv", help="Output file. Default: ./formated_refgene.tsv")
	parser_b.add_argument("-u5", "--include_utr_5", default=False, action='store_true')
	parser_b.add_argument("-u3", "--include_utr_3", default=False, action='store_true')
	parser_b.add_argument("-m", "--include_chrM", default=False, action='store_true')
	parser_b.add_argument("-ncc", "--include_non_canonical_chr", default=False, action='store_true', help="Adapted (tested) only with human hg19 chr")
	parser_b.add_argument("-nct", "--include_non_coding_transcripts", default=False, action='store_true', help="Activating this will filter out  all transcripts that are not NM_....")
	parser_b.set_defaults(mode=main_format_only)
	
	args = parser.parse_args()
	args.mode(args)
	
