#!/usr/bin/env python
# encoding: utf-8
# Title: "Electric Fancy Fox . Python"
# written by sebHTML (S. Boisvert)
# Created 2012-06-11
# Updated 2012-10-17 by plpla (P-L. Plante)
# Updated 2014-07-16 by Antony Le Béchec
# Updated 2015-07-16 by Antony Le Béchec
# Updated 2015-10-14 by Antony Le Béchec
# Updated 2015-10-21 by Antony Le Béchec
# Updated 2016-01-19 by Antony Le Béchec
# Updated 2023-05-12 by Antony Le Béchec
# reason for creating this software:
#  the last Illumina(R) MiSeq(R) run failed to convert the data on board
# because of a faulty cluster density.

"""

The Illumina(R) MiSeq(TM) from Illumina Inc.
produces a sample sheet file on board, but it is not compatible
with Illumina(R) CASAVA.




Legal stuff.

ILLUMINA® is a registered trademark of Illumina, Inc.


<This is a Illumina(R) MiSeq(R) SampleSheet>

[Header]
IEMFileVersion,4
Investigator Name,LR
Project Name,LSPQ
Experiment Name,20120606
Date,6/6/2012
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,TruSeq DNA/RNA
Description,
Chemistry,Default

[Reads]
151
151

[Settings]

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Sample_Project,index,I7_Index_ID,Description,GenomeFolder
90377-1,,LSPQ20120606,A01,,ATCACG,A001,,PhiX\Illumina\RTA\Sequence\WholeGenomeFASTA
90443-1,,LSPQ20120606,A02,,CGATGT,A002,,PhiX\Illumina\RTA\Sequence\WholeGenomeFASTA
90941-1,,LSPQ20120606,A03,,TTAGGC,A003,,PhiX\Illumina\RTA\Sequence\WholeGenomeFASTA
91375-1,,LSPQ20120606,A04,,TGACCA,A004,,PhiX\Illumina\RTA\Sequence\WholeGenomeFASTA
91378-1,,LSPQ20120606,A05,,ACAGTG,A005,,PhiX\Illumina\RTA\Sequence\WholeGenomeFASTA
91655-1,,LSPQ20120606,A06,,GCCAAT,A006,,PhiX\Illumina\RTA\Sequence\WholeGenomeFASTA
91664-1,,LSPQ20120606,A07,,CAGATC,A007,,PhiX\Illumina\RTA\Sequence\WholeGenomeFASTA
91845-1,,LSPQ20120606,A08,,CTTGTA,A012,,PhiX\Illumina\RTA\Sequence\WholeGenomeFASTA
92375-1,,LSPQ20120606,A09,,AGTCAA,A013,,PhiX\Illumina\RTA\Sequence\WholeGenomeFASTA
93302-1,,LSPQ20120606,A10,,AGTTCC,A014,,PhiX\Illumina\RTA\Sequence\WholeGenomeFASTA
93363-1,,LSPQ20120606,A11,,ATGTCA,A015,,PhiX\Illumina\RTA\Sequence\WholeGenomeFASTA
90494-1,,LSPQ20120606,A12,,CCGTCC,A016,,PhiX\Illumina\RTA\Sequence\WholeGenomeFASTA
90505-1,,LSPQ20120606,B01,,GTCCGC,A018,,PhiX\Illumina\RTA\Sequence\WholeGenomeFASTA
93302-2,,LSPQ20120606,B02,,GTGAAA,A019,,PhiX\Illumina\RTA\Sequence\WholeGenomeFASTA

</This is a Illumina(R) MiSeq(R) SampleSheet>

for dual indexes:

Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description

<This is a Illumina(R) CASAVA/BCL2FASTQ(R) SampleSheet>

FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject
000000000-A0VER,1,R6_WT,,GCCAAT,,N,PE_indexing,FR,RNA-seq_spr0058_20120604
000000000-A0VER,1,R6_KOspr0058,,CTTGTA,,N,PE_indexing,FR,RNA-seq_spr0058_20120604

</This is a Illumina(R) CASAVA/BCL2FASTQ(R) SampleSheet>

"""

import sys

programName=sys.argv[0]

usage="# Usage:\n\
# "+programName+" <inputMiSeqSampleSheet> <output> <outputType> <flowcellIdentifier>\n\
#   <inputMiSeqSampleSheet>: MiSeq SampleSheet (usually 'SampleSheet.csv' at run's root )\n\
#   <output>: Output file to generate\n\
#   <outputType>: Type of output file to generate (default: 'bcl2fastqSampleSheet')\n\
#      'casavaSampleSheet' (CASAVA SampleSheet)\n\
#      'bcl2fastqSampleSheet' (BCL2FASTQ SampleSheet)\n\
#      'readsLength' (File with the reads length by line)\n\
#      'mask' (File mask option for CASAVA/BCL2FASTQ demultiplexing)\n\
#      'adapters' (File with the adapters in fasta format)\n\
#   <flowcellIdentifier>: flow cell Identifier, can be the name/foldername of the run\n\
"

if len(sys.argv)<3:
	print(usage)
	sys.exit(1)

# input SampleSheet file
miseqSampleSheet=sys.argv[1]
if miseqSampleSheet.strip()=="":
	print("# Error: No inputMiSeqSampleSheet file defined  (e.g. 'SampleSheet.csv')")
	print(usage)
	sys.exit(1)
try:
	open(miseqSampleSheet)
except IOError:
	print("# Error: '"+miseqSampleSheet+"' file does not exist!")
	print(usage)
	sys.exit(1)

# output file
output=sys.argv[2]
if output.strip()=="":
	print("# Error: No output file defined (e.g. 'SampleSheet.casava.csv')")
	print(usage)
	sys.exit(1)

# output type
outputTypeDefault="bcl2fastqSampleSheet"
try:
	outputType=sys.argv[3]
except IndexError:
	print("# Warning: No output type defined, '"+outputTypeDefault+"' will be used ")
	outputType=outputTypeDefault
if outputType.strip()=="":
	print("# Warning: No output type defined, '"+outputTypeDefault+"' will be used ")
	outputType=outputTypeDefault

# flowcellID
flowcellIdentifierDefault="Unknown"
try:
	flowcellIdentifier=sys.argv[4]
except IndexError:
	flowcellIdentifier=flowcellIdentifierDefault
if flowcellIdentifier.strip()=="":
	flowcellIdentifier=flowcellIdentifierDefault


# Variables

lane="1"
machineOperator="Unknown"
projectName="Unknown"
InvestigatorName="Unknown"

operationCode_Header="[Header]"
operationCode_Manifests="[Manifests]"
operationCode_Reads="[Reads]"
operationCode_Settings="[Settings]"
operationCode_Data="[Data]"


print("# Welcome to this tool.")

print("# I will process the MiSeq sample sheet <"+miseqSampleSheet+"> today.")


cavasaSampleSheet_content="FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\n"
bcl2fastqSampleSheet_content="FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\n"

manifests_content=""
manifests_content_list=""

nb_reads=0
reads_length_array=[]
readsLength_content=""

mask_content=""
maskRead1=""
maskRead2=""
maskIndex=""

adapters_content=""

gotHeader=False

operationCode=""

for line in open(miseqSampleSheet):
	
	linestrip=line.strip()

	if linestrip=="":
		continue

	linestripcols=linestrip.split(",")
	linestripcol1=linestripcols[0]
	if (linestripcol1==operationCode_Header
	or linestripcol1==operationCode_Manifests
	or linestripcol1==operationCode_Reads
	or linestripcol1==operationCode_Settings
	or linestripcol1==operationCode_Data):
		#operationCode=line.strip()
		operationCode=linestripcol1
		print("# "+operationCode+" processing...")
		continue

	#print"# "+operationCode+" and "+linestripcol1+" ??? "

	# [Header]
	if operationCode==operationCode_Header:
		#print"# "+operationCode+" == HEADER "
		if line.find("Project Name,")>=0:
			projectName=line.replace("Project Name,","").strip()
			print("#    Project Name: "+projectName)
		if line.find("Investigator Name,")>=0:
			InvestigatorName=line.replace("Investigator Name,","").strip()
			print("#    Investigator Name: "+InvestigatorName)
		if line.find("Assay,")>=0:
			Assay=line.replace("Assay,","").strip()
			print("#    Assay: "+Assay)
	
	# [Manifests]
	if operationCode==operationCode_Manifests:
		#print"######### "+line
		if line.strip()!="":
			tokens=line.strip().split(",")
			if tokens[1].strip()!="":
				if manifests_content.strip()=="":
					manifests_content=tokens[1]
				else:
					manifests_content=manifests_content+","+tokens[1]
				# Manifest2
				manifests_content_list=manifests_content_list+""+tokens[0]+","+tokens[1]+"\n"

	# [Reads]
	if operationCode==operationCode_Reads:
		if line.strip()!="":
			
			# Sequencing Method
			nb_reads=nb_reads+1
			reads_length_split=line.strip().split(",")
			reads_length=reads_length_split[0]
			
			if reads_length!="":
				reads_length_array.append(reads_length)
				print("#    Reads "+str(nb_reads)+": "+str(reads_length_array[(nb_reads-1)])+"pb")
				readsLength_content=readsLength_content+reads_length+"\n"
			
				sequencing_method="Unknown"
				if nb_reads==1:
					sequencing_method="SE"
					maskRead1="Y"+reads_length_array[0]
				else:
					if nb_reads==2:
						sequencing_method="PE"
						maskRead1="Y"+reads_length_array[0]
						maskRead2="Y"+reads_length_array[1]
	
	# [Settings]
	if operationCode==operationCode_Settings:
		if line.find("Adapter")>=0:
			tokens=line.strip().split(",")
			adapters_content=adapters_content+">"+tokens[0]+"\n"+tokens[1]+"\n"
	
	# [Data]
	if operationCode==operationCode_Data:
		#printline
		tokens=line.split(",")
		
		if not gotHeader: # [Data] Header First line !
			
			# main information
			sampleProject="Unknown"
			if projectName!="" and projectName!="Unknown":
				sampleProject=projectName
			else:
				if InvestigatorName!="" and InvestigatorName!="Unknown":
					sampleProject=InvestigatorName
					
			machineOperator="Unknown"
			if InvestigatorName!="":
				machineOperator=InvestigatorName

			# Index determination
			I5_Index_ID=0
			I7_Index_ID=0
			for i in range(len(tokens)):
				#if tokens[i]=="I5_Index_ID":
				if tokens[i]=="index2":
					I5_Index_ID=i
				if tokens[i]=="index":
					I7_Index_ID=i

			# Dual index determination
			dual=False
			for i in range(len(tokens)):
				if tokens[i]=="I5_Index_ID":
					dual=True
			#if Assay=="HaloPlex":
			#	dual=False

			if dual:
				print("#    Dual index detected")
			else:
				print("#    Single index detected")

			gotHeader=True

			continue

		#tokens=line.split(",")
	
		# MiSeq format
		# Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Sample_Project,index,I7_Index_ID,Description,GenomeFolder
		# 90377-1,,LSPQ20120606,A01,,ATCACG,A001,,PhiX\Illumina\RTA\Sequence\WholeGenomeFASTA

		sampleName=tokens[1].strip()

		# take the identifier if the name is empty.
		if sampleName=="":
			print("#    Warning: Sample_Name is empty, will use Sample_ID for the sample name")
			sampleName=tokens[0].strip()

		if sampleName=="":
			print("#    Error: sample name is empty.")
			sys.exit(1)

		# replace space " " by underscore "_"
		#print"# SampleName0="+sampleName
		#sampleName=sampleName.replace(' ', '_')
		#print"# SampleName1="+sampleName
		
		# Skip name with space
		if " " in sampleName:
			print("#    Warning: sample name 'sampleName' contains a space. This sample will be skipped.")
			continue
		
		# index
		
		
		#if len(tokens)==10:
		if dual:
			if Assay=="HaloPlex":
				index=tokens[5].strip() #+"-"+tokens[7].strip()
				maskIndex="I"+str(len(tokens[I7_Index_ID].strip()))+","+"Y"+str(len(tokens[I5_Index_ID].strip()))
			else:
				index=tokens[5].strip()+"-"+tokens[7].strip()
				maskIndex="I"+str(len(tokens[I7_Index_ID].strip()))+","+"I"+str(len(tokens[I5_Index_ID].strip()))
		else:
			index=tokens[5].strip()
			maskIndex="I"+len(tokens[5].strip())
		

		# CASAVA/BCL2FASTQ format
		#FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject
		#000000000-A0VER,1,R6_WT,,GCCAAT,,N,PE_indexing,FR,RNA-seq_spr0058_20120604

		cavasaSampleSheet_content=cavasaSampleSheet_content+flowcellIdentifier+","	# FCID
		cavasaSampleSheet_content=cavasaSampleSheet_content+lane+","			# Lane
		cavasaSampleSheet_content=cavasaSampleSheet_content+sampleName+","		# SampleID
		cavasaSampleSheet_content=cavasaSampleSheet_content+","				# SampleRef
		cavasaSampleSheet_content=cavasaSampleSheet_content+index+","			# Index
		cavasaSampleSheet_content=cavasaSampleSheet_content+","
		cavasaSampleSheet_content=cavasaSampleSheet_content+"N,"
		cavasaSampleSheet_content=cavasaSampleSheet_content+"PE_indexing,"
		cavasaSampleSheet_content=cavasaSampleSheet_content+machineOperator+","		# Operator
		cavasaSampleSheet_content=cavasaSampleSheet_content+sampleProject
		cavasaSampleSheet_content=cavasaSampleSheet_content+"\n"
		
		bcl2fastqSampleSheet_content=bcl2fastqSampleSheet_content+flowcellIdentifier+","	# FCID
		bcl2fastqSampleSheet_content=bcl2fastqSampleSheet_content+lane+","			# Lane
		bcl2fastqSampleSheet_content=bcl2fastqSampleSheet_content+sampleName+","		# SampleID
		bcl2fastqSampleSheet_content=bcl2fastqSampleSheet_content+","				# SampleRef
		bcl2fastqSampleSheet_content=bcl2fastqSampleSheet_content+index+","			# Index
		bcl2fastqSampleSheet_content=bcl2fastqSampleSheet_content+","
		bcl2fastqSampleSheet_content=bcl2fastqSampleSheet_content+"N,"
		bcl2fastqSampleSheet_content=bcl2fastqSampleSheet_content+"PE_indexing,"
		bcl2fastqSampleSheet_content=bcl2fastqSampleSheet_content+machineOperator+","		# Operator
		bcl2fastqSampleSheet_content=bcl2fastqSampleSheet_content+sampleProject
		bcl2fastqSampleSheet_content=bcl2fastqSampleSheet_content+"\n"

# Mask
mask=maskRead1+","+maskIndex
if sequencing_method=="PE":
	mask=mask+","+maskRead2
mask_content=mask

	
# Write file
stream=open(output,"w")
if outputType=="bcl2fastqSampleSheet":
	stream.write(bcl2fastqSampleSheet_content)
if outputType=="casavaSampleSheet":
	stream.write(cavasaSampleSheet_content)
if outputType=="readsLength":
	stream.write(readsLength_content)
if outputType=="mask":
	stream.write(mask_content)
if outputType=="adapters":
	stream.write(adapters_content)
if outputType=="manifests":
	stream.write(manifests_content)
if outputType=="manifests_list":
	stream.write(manifests_content_list)
stream.close()



print("# Output file <"+output+"> ("+outputType+") is now on the disk")
print("# Have a nice day.")
