#!/usr/bin/perl
############################
# DB data intégration      #
# Author: Antony Le Béchec #
# Copyright: IRC           #
############################

## Main Information
#####################

our %information = ( #
	'release'	=>  	"0.9.2beta",	# Release
	#'beta'		=>  	"beta",		# Man parameter
	'date'		=>  	"20150212",	# Release parameter
	'author'	=>  	"ALB",		# Debug parameter
	'copyright'	=>  	"IRC-GNUGPL",	# Verbose parameter
);


## Modules
############

use Getopt::Long;		# Catch Options
use Pod::Usage;			# Pod
use Time::localtime;		# Time
use Data::Dumper;		# Data
use File::Basename;		# File
use Switch;			# Switch
#use File::Temp qw/ tempfile tempdir tmpnam /;
use File::Temp qw/ tempfile tempdir tmpnam /;
use lib dirname (__FILE__);	# Add lib in the same folder
use Scalar::Util qw(looks_like_number);
use Digest::MD5 qw(md5 md5_hex md5_base64);
use POSIX qw/floor/;

require "functions.inc.pl";	# Common functions


## HELP/MAN
#############

=head1 NAME

DBintegration.pl - Integrate data on DB

=head1 DESCRIPTION

Description

=head1 BUGS

Bugs...

=head1 ACKNOWLEDGEMENTS

Thank U!

=head1 COPYRIGHT

IRC - GNU GPL License

=head1 AUTHOR

ALB

=head1 USAGE

$ARGV[0] [options] --type=<STRING> --input_data=<STRING>

=head1 OPTIONS

=head2 MAIN options

=over 2

=item B<--help|h|?>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=back

=head2 CONFIGURATION

=over 2

=item B<--config|config_file=<file>>

Configuration file for main parameters (default 'config.ini')

=item B<--config_annotation|config_annotation_file=<file>>

Configuration file for annotation parameters (default 'config.annotation.ini').

=back

=head2 INPUT

=over 2

=item B<--type=<string>>

Type of data to integrate (e.g. 'run', 'sample', 'vcf')

=item B<--data=<string>>

Data to integrate, as a list of variable/value separated by a comma. The association variable/value separater is either ':' or '='.

Format example: var1=val1,var2=val2,var3:val3

=back

=head2 OUTPUT

=over 2

=item B<--file_id=<file>>

File with the id of the generated object

=back

=cut


## Parameters
###############

## Parameters default values

#$annotation_default="PZScore,PZFlag,Symbol,hgvs,location,outcome,AlleleFrequency,AD,DP,AF,GQ,Ensembl";
$annotation_default="ALL";
#$sort_by_default="PZFlag,PZScore";
#$order_by_default="DESC,DESC";

our %parameters = ( #
	# Main options 
	'help'		=>  	0,	# Help parameter
	'man'		=>  	0,	# Man parameter
	'release'	=>  	0,	# Release parameter
	'debug'		=>  	0,	# Debug parameter
	'verbose'	=>  	0,	# Verbose parameter
	# Configuration
	'config'		=>	'config.ini',			# Configuration file
	'config_annotation'    	=>  	"config.annotation.ini",	# Configuration annotation file
	#'config_filter'    	=>  	"config.filter.ini",		# Configuration filter file
	# Input
	'type'		=>	undef,	# Type of data to integrate
	'data'		=>	undef,	# Data to integrate (format depend on type)
	# Output
	'file_id'	=>	undef,	# Output file with id
);

## Parameters definition

our @options=(
	# Main options
	'help|h|?',		# Help
	'man',			# Man
	'release',		# Release
	'debug',		# Debug
	'verbose',		# Verbose
	# Configuration
	'config|config_file=s',				# Configuration file
	'config_annotation|config_annotation_file=s',	# Configuration annotation file
	#'config_filter|config_filter_file=s',		# Configuration filter file
	# Input
	'type=s',	# Type of data to integrate
	'data=s',	# Data to integrate
	# output
	'file_id=s',	# Output file with id
);

## Catch Options and put into parameters
GetOptions (\%parameters,@options,@common_options)  or pod2usage();

## Main parameters
$date=timestamp();
$basename = dirname (__FILE__);

## Header
$header="##\n";
$header.="## Script: ".$information{"script"}." (release ".$information{"release"}."/".$information{"date"}.")\n";
$header.="## Excecution Date: ".$date."\n";
$header.="##\n";


## Parameters
###############

require "parameters.inc.pl";	# Parameters

## PrePorcessing
##################


# TYPE
my $type=$parameters{"type"};
if (!defined $type) {
	warn "#[ERROR] No type defined\n";
	pod2usage(2);
};#if

# DATA
my $input_data=$parameters{"data"};
if (!defined $input_data) {
	warn "#[ERROR] No data defined\n";
	pod2usage(2);
};#if

# Data input process
%data = map{split /[=:]/, $_, 2}(split /,/, $input_data);

# FILE_ID
my $output_file_id=$parameters{"file_id"};


## DEBUG
##########

#print Dumper(\%data) if $DEBUG;
#print Dumper(\%config_annotation) if $DEBUG;
#print Dumper(\%config_filter) if $DEBUG;
#print Dumper(\%annotation_filters) if $DEBUG;


## MAIN
#########

# Variables
$output="";


my %variants;
my %annotations;

if ($VERBOSE) {
	$verbose.="# Type: $type\n";
	$verbose.="# Data:\n";	
	while ((my $var, $val) = each(%data)){
		$verbose.="# \t$var \t$val\n";	
	};#while
	$verbose.="# \n";	

};#if


sub insert_group {
# Insert a new RUN if not exists, and return the ID

	# INPUT
	my $data=shift;		# Data list
	#print Dumper(\%$data) if $DEBUG;

	# Variables
	my $id;
	my $query;
	my $type="group";
	warn "# $type Insertion\n" if $DEBUG;

	# CHECK input ID
	$id = check_id("SELECT id FROM `$type` WHERE id='".$data{"$type.id"}."'");
	return return_message($id,"","# $type ID '$id' in input","") if defined $id;

	# CHECK input DATA

	my $ref=$data{"$type.ref"};
	return return_message(undef,"","","#[ERROR] $type REF not defined in input") if (!defined $ref || trim($ref) eq "");
	
	# SEARCH ID
	$query="SELECT id 
		FROM `$type`
		WHERE ref='$ref'
		";
	$id = check_id($query);
	return return_message($id,"","# $type ID '$id' found in DB","") if defined $id;
	
	# INSERT
	$query="INSERT INTO `$type` (`ref`)
		   VALUES ('$ref')
		";
	$id = insert_id($query);
	return return_message($id,"","# $type ID '$id' inserted","") if defined $id;

	# ERROR insertion
	return return_message(undef,"","","#[ERROR] $type insertion failed!");

}

sub insert_project {
# Insert a new RUN if not exists, and return the ID

	# INPUT
	my $data=shift;		# Data list
	#print Dumper(\%$data) if $DEBUG;

	# Variables
	my $id;
	my $query;
	my $type="project";
	warn "# $type Insertion\n" if $DEBUG;

	# CHECK input ID
	$id = check_id("SELECT id FROM `$type` WHERE id='".$data{"$type.id"}."'");
	return return_message($id,"","# $type ID '$id' in input","") if defined $id;

	# CHECK input DATA

	my $ref=$data{"$type.ref"};
	return return_message(undef,"","","#[ERROR] $type REF not defined in input") if (!defined $ref || trim($ref) eq "");

	# Project 'project.id'
	my $group_id=insert_group(\%data);
	return return_message(undef,"","","#[ERROR] $type insertion failed! project ID failed") if !defined $group_id;
	
	# SEARCH ID
	$query="SELECT id 
		FROM `$type`
		WHERE ref='$ref'
		  AND group_id=$group_id
		";
	#print "$query\n" if $DEBUG;
	$id = check_id($query);
	return return_message($id,"","# $type ID '$id' found in DB","") if defined $id;
	
	# INSERT
	$query="INSERT INTO `$type` (`ref`,`group_id`)
		   VALUES ('$ref',$group_id)
		";
	$id = insert_id($query);
	return return_message($id,"","# $type ID '$id' inserted","") if defined $id;

	# ERROR insertion
	return return_message(undef,"","","#[ERROR] $type insertion failed!");

}



sub insert_manifest {
# Insert a new RUN if not exists, and return the ID

	# INPUT
	my $data=shift;		# Data list
	#print Dumper(\%$data) if $DEBUG;

	# Variables
	my $id;
	my $query;
	my $type="manifest";
	warn "# $type Insertion\n" if $DEBUG;

	# CHECK input ID
	$id = check_id("SELECT id FROM `$type` WHERE id='".$data{"$type.id"}."'");
	return return_message($id,"","# $type ID '$id' in input","") if defined $id;

	# CHECK input DATA

	# Manifest file 'manifest'
	my $manifest;
	my $manifest_filename=$data{"manifest"};
	if (defined $manifest_filename && -e $manifest_filename) {
		#$manifest_checked=1;
		open my $f, '<', $manifest_filename or die "#[ERROR] Opening File '$!' ";
		$manifest=do{local $/; <$f>}
	} else {
		#warn "# No Manifest 'manifest' '$manifest''\n" if $DEBUG;
		#return undef;
	};#if

	# Reference 'ref'
	my $ref=$data{"manifest.ref"};
	$ref=basename($manifest_filename) if (!defined $ref || trim($ref) eq ""); # AUTO REF
	
	# SEARCH ID
	$query="SELECT id
		FROM `manifest`
		WHERE manifest_md5=md5(?) AND ref='$ref'
		UNION
		SELECT id
		FROM `manifest`
		WHERE manifest_md5=md5(?)
		UNION
		SELECT id
		FROM `manifest`
		WHERE ref='$ref'
		";
	my @values=($manifest,$manifest);
	$id = check_id($query,\@values);
	return return_message($id,"","# $type ID '$id' found in DB","") if defined $id;
	
	# INSERT check
	return return_message(undef,"","","#[ERROR] $type insertion failed! Input data error") if (trim($manifest) eq "" || trim($ref) eq "");

	# INSERT
	@values=($manifest,$manifest);
	$query="INSERT INTO `manifest` (`ref`,`manifest`,`manifest_md5`)
		   VALUES ('$ref',?,md5(?))
		";
	$id = insert_id($query,\@values);
	return return_message($id,"","# $type ID '$id' inserted","") if defined $id;

	# ERROR Insertion
	return return_message(undef,"","","#[ERROR] $type insertion failed!")
	#warn "#[ERROR] $type Insertion failed!\n" if $DEBUG;
	#return undef;

}

sub insert_snapshot {
# Insert a new SNAPSHOT if not exists, and return the ID

	# INPUT
	my $data=shift;		# Data list
	print Dumper(\%$data) if $DEBUG;

	# Variables
	my $id;
	my $query;
	my $type="snapshot";
	warn "# $type Insertion\n" if $DEBUG;

	# CHECK input ID
	$id = check_id("SELECT id FROM `$type` WHERE id='".$data{"$type.id"}."'");
	return return_message($id,"","# $type ID '$id' in input","") if defined $id;

	# CHECK input DATA

	# SnapShot file 'snapshot'
	my $snapshot;
	my $snapshot_filename=$data{"snapshot"};
	if (defined $snapshot_filename && -e $snapshot_filename) {
		#$snapshot_checked=1;
		open my $f, '<', $snapshot_filename or die "#[ERROR] Opening File '$!' ";
		$snapshot=do{local $/; <$f>}
		#open(my $fh, $snapshot_filename) or die $!;
		#read($fh, $snapshot, "-s $fh");
	} else {
		#warn "# No Manifest 'manifest' '$manifest''\n" if $DEBUG;
		#return undef;
	};#if

	# VCF 'vcf_id'
	my $vcf_id=insert_vcf(\%data);
	return return_message(undef,"","","#[ERROR] $type insertion failed! run ID failed") if !defined $vcf_id;
	insert_snapshot_vcf($vcf_id,$id) if defined $vcf_id && defined $id;

	# Reference 'ref'
	my $ref=$data{"$type.ref"};
	$ref=basename($snapshot_filename) if (!defined $ref || trim($ref) eq ""); # AUTO REF

	# CHROM 'chrom'
	my $chrom=$data{"$type.chrom"};
	$chrom=~s/^chr//; 
	#$ref=basename($snapshot_filename) if (!defined $ref || trim($ref) eq ""); # AUTO REF

	# POS 'pos'
	my $pos=$data{"$type.pos"};
	#$ref=basename($snapshot_filename) if (!defined $ref || trim($ref) eq ""); # AUTO REF

	# START 'start'
	my $start=$data{"$type.start"};
	#$ref=basename($snapshot_filename) if (!defined $ref || trim($ref) eq ""); # AUTO REF

	# STOP 'stop'
	my $stop=$data{"$type.stop"};
	#$ref=basename($snapshot_filename) if (!defined $ref || trim($ref) eq ""); # AUTO REF

	# BAM 'bam'
	my $bam=$data{"$type.bam"};
	#$ref=basename($snapshot_filename) if (!defined $ref || trim($ref) eq ""); # AUTO REF
	
	# SEARCH ID
	$query="SELECT id
		FROM `$type`
		WHERE snapshot_md5=md5(?) AND ref='$ref'
		UNION
		SELECT id
		FROM `$type`
		WHERE snapshot_md5=md5(?)
		UNION
		SELECT id
		FROM `$type`
		WHERE ref='$ref'
		";
	my @values=($snapshot,$snapshot);
	$id = check_id($query,\@values);
	insert_snapshot_vcf($vcf_id,$id) if defined $vcf_id && defined $id;
	return return_message($id,"","# $type ID '$id' found in DB???","") if defined $id;
	
	# INSERT check
	return return_message(undef,"","","#[ERROR] $type insertion failed! Input data error") if (trim($snapshot) eq "" || trim($ref) eq "");

	# INSERT
	@values=($snapshot,$snapshot);
	$query="INSERT IGNORE INTO `snapshot` (`ref`,`snapshot`,`snapshot_md5`,`chrom`,`pos`,`start`,`stop`,`bam`)
		   VALUES ('$ref',?,md5(?),$chrom,$pos,$start,$stop,'$bam')
		";
	$id = insert_id($query,\@values);
	# Association
	insert_snapshot_vcf($vcf_id,$id);
	#if (defined $id) {
	#	@values=();
	#	$query="INSERT INTO `snapshot_vcf` (`vcf_id`,`snapshot_id`)
	#		   VALUES ($vcf_id,$id)
	#		";
	#	$id_association = insert_id($query,\@values);
	#};#if
	return return_message($id,"","# $type ID '$id' inserted","") if defined $id;

	# ERROR Insertion
	return return_message(undef,"","","#[ERROR] $type insertion failed!")
	#warn "#[ERROR] $type Insertion failed!\n" if $DEBUG;
	#return undef;

}

sub insert_snapshot_vcf {
# Insert a snapshot_vcf between a VCF and a SnapShot
	
	# INPUT
	my $vcf_id=$_[0];	# VCF ID
	my $snapshot_id=$_[1];	# SnapShot ID

	if (defined $vcf_id && trim($vcf_id) ne "" && defined $snapshot_id && trim($snapshot_id) ne "") {
		my $query="INSERT IGNORE INTO `snapshot_vcf` (`vcf_id`,`snapshot_id`)
			   VALUES ($vcf_id,$snapshot_id)
			";
		$id = insert_id($query);
		return return_message(1,"","# snapshot_vcf '$snapshot_id.$vcf_id' inserted","");
		#print "ID2: '$id2' form $QC_id,$object_id,'$object_type'\n" if $DEBUG;
	} else {
		return return_message(undef,"","","#[ERROR] 'insert_snapshot_vcf' function failed!");
	};#if
	return 1;
}


sub insert_qc {
# Insert a new QC if not exists, and return the ID

	# INPUT
	my $data=shift;		# Data list
	#print Dumper(\%$data) if $DEBUG;

	# Variables
	my $id;
	my $query;
	my $type="QC";
	warn "# $type Insertion\n" if $DEBUG;

	# ASSOCIATION
	# OBJECT
	my $object_id=$data{"QC.object_id"};
	my $object_type=$data{"QC.object_type"};
	if (defined $object_type && trim($object_type) ne "" && (trim($object_type) eq "vcf" || trim($object_type) eq "sample" || trim($object_type) eq "run")) {
		$object_id=insert_vcf(\%data) if ($object_type eq "vcf");
		$object_id=insert_sample(\%data) if ($object_type eq "sample");
		$object_id=insert_run(\%data) if ($object_type eq "run");
	} elsif (!defined $object_id || trim($object_id) eq "" || !defined $object_type || trim($object_type) eq "") {
		$run_id=insert_run(\%data);
		$sample_id=insert_sample(\%data);
		$vcf_id=insert_vcf(\%data);
		if (defined $vcf_id) {
			$object_id=$vcf_id;
			$object_type="vcf";
		} elsif (defined $sample_id) {
			$object_id=$sample_id;
			$object_type="sample";
		} elsif (defined $run_id) {
			$object_id=$run_id;
			$object_type="run";
		};#if
		#$object_id=basename($QC_filename) if (!defined $ref || trim($ref) eq ""); # AUTO REF
	};#if

	# CHECK input ID
	$id = check_id("SELECT id FROM `$type` WHERE id='".$data{"$type.id"}."'");
	my $association=insert_qc_association($id,$object_id,$object_type) if defined $id;
	return return_message($id,"","# $type ID '$id' in input","") if defined $id && defined $association;

	# CHECK input DATA

	# QC file 'QC'
	my $QC="";
	my $QC_filename=$data{"QC"};
	if (defined $QC_filename && -e $QC_filename) {
		#$manifest_checked=1;
		open my $f, '<', $QC_filename or die "#[ERROR] Opening File '$!' ";
		$QC=do{local $/; <$f>}
	} else {
		#warn "# No Manifest 'manifest' '$manifest''\n" if $DEBUG;
		#return undef;
	};#if

	# Reference 'ref'
	my $ref=$data{"QC.ref"};
	$ref=basename(dirname($QC_filename))."/".basename($QC_filename) if (!defined $ref || trim($ref) eq ""); # AUTO REF
	#print "REF: $ref\n" if $DEBUG;

	# Reference 'metrics'
	my $metrics=$data{"QC.metrics"};
	($metrics)= basename($QC_filename) =~ /\.([^.]+)$/ if (!defined $metrics || trim($metrics) eq ""); # AUTO metrics
	#print "metrics: $metrics\n" if $DEBUG;
	
	# Reference 'type'
	my $qc_type=$data{"QC.type"};
	#($metrics)= basename($QC_filename) =~ /\.([^.]+)$/ if (!defined $metrics || trim($metrics) eq ""); # AUTO metrics
	#print "type: $qc_type\n" if $DEBUG;
	
	# SEARCH ID
	if (1) { # To reduce space on DB # BUG!!!
		$query="SELECT QC.id
			FROM QC
			WHERE QC.ref='$ref'
			  AND QC.metrics='$metrics'
			  AND QC.type='$qc_type'
			  AND QC.QC=?
			";
		print "QUERY: $query\n" if $DEBUG;
		my @values=($QC);
		$id = check_id($query,\@values);
		#$id = check_id($query);
		print "QC.ID FOUND: $id\n" if $DEBUG;
		my $association=insert_qc_association($id,$object_id,$object_type) if defined $id;
		print "ASSOCIATION: $association\n" if $DEBUG;
		return return_message($id,"","# $type ID '$id' found in DB","") if defined $id && defined $association;
	};#if

	# INSERT check
	return return_message(undef,"","","#[ERROR] $type insertion failed! Input data error") if (trim($QC) eq "" || trim($ref) eq "" || trim($object_id) eq "" || trim($object_type) eq "");

	# INSERT
	my @values=($QC,$QC);
	$query="INSERT INTO `QC` (`ref`,`metrics`,`type`,`QC`,`QC_md5`)
		   VALUES ('$ref','$metrics','$qc_type',?,md5(?))
		";
	$id = insert_id($query,\@values);

	# INSERT QC Association
	my $association=insert_qc_association($id,$object_id,$object_type) if defined $id;
	return return_message($id,"","# $type ID '$id' inserted","") if defined $id && defined $association;
	#,`object_id`,`object_type` 

	# INSERT QC special
	#if ($metrics eq "depth") {
	#	insert_qc_depth
	#	my $return_file_insertion=insert_qc_depth_file($id,$QC_filename);
	#};#if

	# ERROR Insertion
	return return_message(undef,"","","#[ERROR] $type insertion failed!")
	#warn "#[ERROR] $type Insertion failed!\n" if $DEBUG;
	#return undef;

}

sub insert_qc_association {
# Insert a QC_association between a QC and a object (id+type)
	
	# INPUT
	my $QC_id=$_[0];	# QC ID
	my $object_id=$_[1];	# Object ID
	my $object_type=$_[2];	# Object TYPE

	if (defined $QC_id && trim($QC_id) ne "" && defined $object_id && trim($object_id) ne "" && defined $object_type && trim($object_type) ne "" ) {
		$query2="INSERT IGNORE INTO `QC_association` (`QC_id`,`object_id`,`object_type`)
			   VALUES ($QC_id,$object_id,'$object_type')
			";
		$id2 = insert_id($query2);
		#print "ID2: '$id2' form $QC_id,$object_id,'$object_type'\n" if $DEBUG;
	} else {
		return return_message(undef,"","","#[ERROR] 'insert_qc_association' function failed!");
	};#if
	return 1;
}

sub insert_qc_depth_file {
# Insert a QC depth file, and return 0 (failed) or 1 (ok)

	# INPUT
	my $id=$_[0];		# input
	my $QC_filename=$_[1];	# input

	# Variables
	my $return;
	my $query;
	my $type="QC depth file";
	warn "# $type insertion\n" if $DEBUG;

	# file QC
	# CREATE TMP TABLE with id NULL
	# INSERT LOAD INFILE $QC_filename
	# UPDATE id with $id
	# INSERT TMP lines into QC_depth TABLE
	

}

sub insert_run {
# Insert a new RUN if not exists, and return the ID

	# INPUT
	my $data=shift;		# Data list
	#print Dumper(\%$data) if $DEBUG;

	# Variables
	my $id;
	my $query;
	my $type="run";
	warn "# $type Insertion\n" if $DEBUG;

	# CHECK input ID
	$id = check_id("SELECT id FROM `$type` WHERE id='".$data{"$type.id"}."'");
	return return_message($id,"","# $type ID '$id' in input","") if defined $id;

	# CHECK input DATA

	# Reference 'ref'
	my $ref=$data{"run.ref"};
	return return_message(undef,"","","#[ERROR] $type REF not defined in input") if (!defined $ref || trim($ref) eq "");

	# SampleSheet file 'samplesheet'
	my $samplesheet=$data{"samplesheet"};
	if (defined $samplesheet && -e $samplesheet) {
		open my $f, '<', $samplesheet or die "#[ERROR] Opening File '$!' ";
		$samplesheet=do{local $/; <$f>}
	} else {
		warn "# No SampleSheet 'samplesheet' '$samplesheet''\n" if $DEBUG;
		return return_message(undef,"","","#[ERROR] $type insertion failed! Input data error 'samplesheet'");
	};#if
	
	# SEARCH ID
	$query="SELECT id 
		FROM `run`
		WHERE ref='$ref'
		";
	$id = check_id($query);
	return return_message($id,"","# $type ID '$id' found in DB","") if defined $id;
	
	# INSERT check
	return return_message(undef,"","","#[ERROR] $type insertion failed! Input data error") if (0);

	# INSERT RUN
	$query="INSERT INTO `run` (`ref`,`samplesheet`)
		   VALUES ('$ref','$samplesheet')
		";
	$id = insert_id($query);
	return return_message($id,"","# $type ID '$id' inserted","") if defined $id;
	
	# ERROR Insertion
	return return_message(undef,"","","#[ERROR] $type insertion failed!");

}


sub insert_sample {
# Insert a new RUN if not exists, and return the ID

	# INPUT
	my $data=shift;		# Data list

	# Variables
	my $id;
	my $query;
	my $type="sample";
	warn "# $type insertion\n" if $DEBUG;

	# CHECK input ID
	$id = check_id("SELECT id FROM `$type` WHERE id='".$data{"$type.id"}."'");
	return return_message($id,"","# $type ID '$id' in input","") if defined $id;

	# CHECK input DATA

	# Reference 'ref'
	my $ref=$data{"sample.ref"};
	return return_message(undef,"","","#[ERROR] $type REF not defined in input") if (!defined $ref || trim($ref) eq "");

	# Run 'run_id'
	my $run_id=insert_run(\%data);
	return return_message(undef,"","","#[ERROR] $type insertion failed! run ID failed") if !defined $run_id;
	
	# Manifest 'manifest_id'
	my $manifest_id=insert_manifest(\%data);
	return return_message(undef,"","","#[ERROR] $type insertion failed! manifest ID failed") if !defined $manifest_id;

	# Project 'project.id'
	my $project_id=insert_project(\%data);
	return return_message(undef,"","","#[ERROR] $type insertion failed! project ID failed") if !defined $project_id;

	
	#warn "REF:$ref\tRUN_ID:$run_id\tMANIFEST_ID:$manifest_id\PROJECT_ID:$project_id\n" if $DEBUG;
	
	# SEARCH ID
	$query="SELECT id 
		FROM `sample`
		WHERE ref='$ref'
		  AND run_id=$run_id
		";
	#  AND manifest_id=$manifest_id
	#  AND project_id=$project_id
	$id = check_id($query);
	return return_message($id,"","# $type ID '$id' found in DB","") if defined $id;

	# INSERT check
	return return_message(undef,"","","#[ERROR] $type insertion failed! Input data error") if (0);

	# INSERT SAMPLE
	$query="INSERT INTO `sample` (`ref`,`run_id`,`manifest_id`,`project_id`)
		   VALUES ('$ref',$run_id,$manifest_id,$project_id)
		";
	$id = insert_id($query);
	return return_message($id,"","# $type ID '$id' inserted","") if defined $id;

	# ERROR Insertion
	return return_message(undef,"","","#[ERROR] $type insertion failed!");

}


sub insert_vcf {
# Insert a new RUN if not exists, and return the ID

	# INPUT
	my $data=shift;		# Data list

	# Variables
	my $id;
	my $query;
	my $type="vcf";
	warn "# $type insertion\n" if $DEBUG;

	# CHECK input ID
	$id = check_id("SELECT id FROM `$type` WHERE id='".$data{"$type.id"}."'");
	return return_message($id,"","# $type ID '$id' in input","") if defined $id;

	# CHECK input DATA

	# VCF file 'vcf'
	my $vcf;
	my $vcf_filename=$data{"vcf"};
	if (defined $vcf_filename && -e $vcf_filename) {
		$vcf_checked=1;
		open my $f, '<', $vcf_filename or die "#[ERROR] Opening File '$!' ";
		$vcf=do{local $/; <$f>}
	} else {
		warn "# No VCF 'vcf' '$vcf''\n" if $DEBUG;
		#return undef;
	};#if

	# BAM 'bam'
	my $bam=$data{"$type.bam"};
	
	# PIPELINE 'pipeline'
	my $pipeline=$data{"$type.pipeline"};
	
	# ALIGNER 'aligner'
	my $aligner=$data{"$type.aligner"};
	
	# CALLER 'caller'
	my $caller=$data{"$type.caller"};
	
	# ANNOTATOR 'annotator'
	my $annotator=$data{"$type.annotator"};
	
	# Reference 'ref'
	my $ref=$data{"$type.ref"};
	$ref=basename($vcf_filename) if (!defined $ref || trim($ref) eq ""); # AUTO REF
	return return_message(undef,"","","#[ERROR] $type REF not defined in input") if (!defined $ref || trim($ref) eq "");

	# Run 'sample_id'
	my $sample_id=insert_sample(\%data);
	return return_message(undef,"","","#[ERROR] $type insertion failed! sample ID failed") if !defined $sample_id;

	# SEARCH ID
	$query="SELECT id 
		FROM `vcf`
		WHERE ref='$ref'
		  AND sample_id=$sample_id
		";
	$id = check_id($query);
	return return_message($id,"","# $type ID '$id' found in DB","") if defined $id;

	# INSERT check
	return return_message(undef,"","","#[ERROR] $type insertion failed! Input data error") if (trim($vcf) eq "" || trim($ref) eq "");

	# INSERT
	$query="INSERT IGNORE INTO `vcf` (`ref`,`sample_id`,`vcf`,`bam`,`pipeline`,`aligner`,`caller`,`annotator`)
		   VALUES ('$ref',$sample_id,?,'$bam','$pipeline','$aligner','$caller','$annotator')
		";
	my @values=($vcf);
	$id = insert_id($query,\@values);

	# INSERT file
	my $return_file_insertion=insert_vcf_file($id,$vcf_filename);
	if (!defined $return_file_insertion) {
		#$query="DELETE FROM `vcf` WHERE id=$id";
		#insert_id($query) if !$return_file_insertion;
		return return_message(undef,"","","#[ERROR] vcf File not inserted") if !$return_file_insertion;
	};#if
	# if error insertion file -> remove ID !!!
	

	return return_message($id,"","# $type ID '$id' inserted","") if defined $id;
	
	# ERROR Insertion
	return return_message(undef,"","","#[ERROR] $type insertion failed!");

}

sub insert_vcf_file {
# Insert a VCF file, and return 0 (failed) or 1 (ok)

	# INPUT
	my $id=$_[0];	# input
	my $vcf=$_[1];	# input
	
	# Variables
	my $return;
	my $query;
	my $split_files=0;
	my $type="VCF file";
	warn "# $type insertion\n" if $DEBUG;

	# Read VCF Header file
	my %vcf_header_hash=read_vcf($vcf,"header");
	#print Dumper(\%vcf_header_hash) if $DEBUG;

	# Read VCF file
	my %vcf_hash=read_vcf($vcf);
	#print Dumper(\%vcf_hash) if $DEBUG;

	# File Variants
	my ($f_variant,$tmp_file_variant)=tempfile();
	#print $tmp_file_variant;
	open my $f_variant, '>', $tmp_file_variant or die "#[ERROR] Opening File '$!' ";
	
	# File Annotations
	my ($f_annotation,$tmp_file_annotation)=tempfile();
	#print $tmp_file_annotation;
	open my $f_annotation, '>', $tmp_file_annotation or die "#[ERROR] Opening File '$!' ";

	# Create TAB file
	my $nb_annotation_line=0;
	my $nb_annotation_content=0;
	my @f_annotation_content_list;
	while ( my ($chr, $poss) = each(%vcf_hash) ) {
	while ( my ($pos, $refs) = each(%{$poss}) ) {
	while ( my ($ref, $alts) = each(%{$refs}) ) {
	while ( my ($allalt, $variant_annotations) = each(%{$alts}) ) {

		# Multiallele
		my $rank=0;

		my $alt_original=$alt;
		my $ref_original=$ref;
		my $pos_original=$pos;
		
		foreach $alt (split(",",$allalt)) {

			#my $alt=$alt_original;
			my $ref=$ref_original;
			my $pos=$pos_original;

			# rank
			$rank++;
			#print "VARIANT: $pos $ref $alt | $allalt $rank\n";
			# Repositionning
			if (length($alt)>1 && length($ref)>1) {
				#print "START $pos $ref $alt\n";
				
				# strip off identical suffixes
				while (substr($alt,-1) == substr($ref,-1) && length($alt)>1 && length($ref)>1) {
					$alt=substr($alt,0,-1);
					$ref=substr($ref,0,-1);
					#print "$pos $ref $alt\n";
				};#if
				
				# strip off identical prefixes and increment position
				while(substr($alt,0,1) == substr($ref,0,1) && length($alt)>1 && length($ref)>1) {
					$alt=substr($alt,1);
					$ref=substr($ref,1);
					$pos++;
					#print "$pos $ref $alt\n";
				};#while

				#print "STOP $pos $ref $alt\n";
				

			};#if

			# VARIANT
			$nb_variant++;
			#print Dumper($variant_annotations) if $DEBUG;
			$chr =~ s/chr//g;
			$chr =~ s/chr//g;
			%replacements = ("X" => "-1", "Y" => "-2", "M" => "-3");
			$chr =~ s/(@{[join "|", keys %replacements]})/$replacements{$1}/g;
	
			# File Variants
			my $QUAL=$$variant_annotations{"QUAL"};
			my $FILTER=$$variant_annotations{"FILTER"};
			my $FORMAT=$$variant_annotations{"FORMAT"};
			my @Q=split(/\|/,$$variant_annotations{"SAMPLE_LIST_CONCAT"});
			my $QUALITY=$Q[0]; #split(/;/,$$variant_annotations{"SAMPLE_LIST_CONCAT"})[0];#SAMPLE_LIST_CONCAT
			my @QUALITY_FORMAT=split(/:/,$FORMAT);
			my @QUALITY_VALUES=split(/:/,$QUALITY);
			my %QUALITIES; @QUALITIES{@QUALITY_FORMAT}=@QUALITY_VALUES;
			#print Dumper(\%QUALITIES) if $DEBUG;
			my $GT=(defined $QUALITIES{"GT"} && trim($QUALITIES{"GT"}) ne "")?$QUALITIES{"GT"}:"\\N";
			my $AD=(defined $QUALITIES{"AD"} && trim($QUALITIES{"AD"}) ne "")?$QUALITIES{"AD"}:"\\N";
			my $DP=(defined $QUALITIES{"DP"} && trim($QUALITIES{"DP"}) ne "")?$QUALITIES{"DP"}:"\\N";
			my $GQ=(defined $QUALITIES{"GQ"} && trim($QUALITIES{"GQ"}) ne "")?$QUALITIES{"GQ"}:"\\N";
			my $PL=(defined $QUALITIES{"PL"} && trim($QUALITIES{"PL"}) ne "")?$QUALITIES{"PL"}:"\\N";
			#print "$id\t\\N\t$chr\t$pos\t$ref\t$alt\t$QUAL\t$FILTER\t$FORMAT\t$QUALITY\t$GT\t$AD\t$DP\t$GQ\t$PL\n" if $DEBUG;

			my $genotype=$GT;

			print $f_variant "$id\t\\N\t$chr\t$pos\t$ref\t$alt\t$rank\t$genotype\t$QUAL\t$FILTER\t$FORMAT\t$QUALITY\t$GT\t$AD\t$DP\t$GQ\t$PL\n";
			#print "@Q[0]\n" if $DEBUG;

			# File Annotations
			my $variant_annotations_infos=$$variant_annotations{"INFOS"};
			while ( my ($annotation_name, $annotation_value) = each(%{$variant_annotations_infos}) ) {
				#print "$id\t\\N\t$chr\t$pos\t$ref\t$alt\t\\N\t$annotation_name\t\\N\t$annotation_value\t\\N\n" if $DEBUG;
				$nb_annotation_line++;
				#my $annotation_type=$config_annotation{$annotation_name}{"annotation_type"};
			


				my $annotation_description=$vcf_header_hash{"INFO"}{$annotation_name}{"Description"};
				$annotation_description=~s/^"|"$//g;
				if (trim($annotation_description) eq "") { $annotation_description=""; };#if NULL

				my %annotation_infos_hash;		
				if ($annotation_description =~ m/^(.*)[(.*)]$/) {
					my $annotation_description_text=$1;
					my $annotation_infos=$2;
					foreach my $ann (split(/;/,$annotation_infos)) {
						(my $var, my $val) = split(/=/,$ann);
						$annotation_infos_hash{$var}=$val;
					};#foreach
				};#if

				my $annotation_type=$annotation_infos_hash{"AnnotationType"};
				if (trim($annotation_type) eq "") { $annotation_type=$config_annotation{$annotation_name}{"annotation_type"}; };#if Check on configuration
				if (trim($annotation_type) eq "") { $annotation_type="unknown"; };#if NULL
				my $annotation_release=$annotation_infos_hash{"Release"};
				if (trim($annotation_release) eq "") { $annotation_release="unknown"; };#if NULL
				my $annotation_date=$annotation_infos_hash{"Date"};
				if (trim($annotation_date) eq "") { $annotation_date="00000000"; };#if NULL
			
				print $f_annotation "$id\t\\N\t$chr\t$pos\t$ref\t$alt\t\\N\tINFO\t$annotation_type\t$annotation_name\t$annotation_release\t$annotation_description\t$annotation_date\t$annotation_value\t$annotation_value\n";
				#if (floor($nb_annotation_line/100)) {
				if ($split_files) {
					$f_annotation_content_list[$nb_annotation_content].="$id\t\\N\t$chr\t$pos\t$ref\t$alt\t\\N\t$annotation_type\t$annotation_name\t$annotation_release\t$annotation_description\t$annotation_date\t$annotation_value\t$annotation_value\n";
					if (($nb_annotation_line % 5000000) == 0) {
						$nb_annotation_content++;
					};#if
				};#if
			};#while

		};#foreach
		
	}}}}

	close($f_variant);
	close($f_annotation);

	# Temporary table id
	my $tmp_table_id = int(rand(1000000000));

	my $query="
		CREATE TEMPORARY TABLE IF NOT EXISTS `tmp_vcf_variant_$tmp_table_id` (
		  `vcf_id` INT,
		  `variant_id` INT DEFAULT NULL,
		  `chr` TINYINT NOT NULL,
		  `pos` INT NOT NULL,
		  `ref` VARCHAR(255) NOT NULL,
		  `alt` VARCHAR(255) NOT NULL,
		  `rank` TINYINT NOT NULL,
		  `genotype` VARCHAR(45) NOT NULL,
		  `QUAL` SMALLINT NOT NULL,
		  `FILTER` VARCHAR(45) NOT NULL,
		  `FORMAT` VARCHAR(100) NOT NULL,
		  `QUALITY` VARCHAR(100) NOT NULL,
		  `GT` VARCHAR(45) NULL,
		  `AD` VARCHAR(45) NULL,
		  `DP` SMALLINT NULL,
		  `GQ` SMALLINT NULL,
		  `PL` VARCHAR(45) NULL
		) ENGINE = InnoDB;

		CREATE TEMPORARY TABLE IF NOT EXISTS `tmp_vcf_annotation_$tmp_table_id` (
		  `vcf_id` INT,
		  `variant_id` INT DEFAULT NULL,
		  `chr` TINYINT NOT NULL,
		  `pos` INT NOT NULL,
		  `ref` VARCHAR(255) NOT NULL,
		  `alt` VARCHAR(255) NOT NULL,
		  `annotation_id` INT NOT NULL,
		  `type` VARCHAR(45) NOT NULL,
		  `annotationType` VARCHAR(45) NOT NULL,
		  `source` VARCHAR(45) NOT NULL,
		  `release` VARCHAR(45) NOT NULL,
		  `description` TEXT NOT NULL,
		  `date` VARCHAR(45) NOT NULL,
		  `value` TEXT NOT NULL,
		  `value_float` FLOAT NULL
		) ENGINE = InnoDB;

		# LOAD DATA VARIANT
		LOAD DATA LOCAL INFILE '$tmp_file_variant'
		INTO TABLE tmp_vcf_variant_$tmp_table_id
		FIELDS TERMINATED BY '\\t'
		ENCLOSED BY ''
		ESCAPED BY '\\\\'
		LINES TERMINATED BY '\\n';

		# INSERT VARIANT
		INSERT IGNORE INTO variant (`chr`,`pos`,`ref`,`alt`)
		SELECT `chr`,`pos`,`ref`,`alt` FROM tmp_vcf_variant_$tmp_table_id;

		# UPDATE VARIANT ID in TMP
		UPDATE tmp_vcf_variant_$tmp_table_id as T
		LEFT JOIN variant as V ON (V.chr=T.chr AND V.pos=T.pos AND V.ref=T.ref AND V.alt=T.alt)
		SET T.variant_id=V.id;

		# INSERT VARIANT_VCF
		INSERT IGNORE INTO variant_vcf (`variant_id`,`vcf_id`,`rank`,`genotype`)
		SELECT `variant_id`,`vcf_id`,`rank`,`genotype` FROM tmp_vcf_variant_$tmp_table_id;

		# INSERT VARIANT_VCF_QUALITY
		INSERT IGNORE INTO variant_vcf_quality (`variant_id`,`vcf_id`,`QUAL`,`FILTER`,`FORMAT`,`QUALITY`,`GT`,`AD`,`DP`,`GQ`,`PL`)
		SELECT `variant_id`,`vcf_id`,`QUAL`,`FILTER`,`FORMAT`,`QUALITY`,`GT`,`AD`,`DP`,`GQ`,`PL` FROM tmp_vcf_variant_$tmp_table_id;

		# LOAD DATA ANNOTATION
		LOAD DATA LOCAL INFILE '$tmp_file_annotation'
		INTO TABLE tmp_vcf_annotation_$tmp_table_id
		FIELDS TERMINATED BY '\\t'
		ENCLOSED BY ''
		ESCAPED BY '\\\\'
		LINES TERMINATED BY '\\n';

		# INSERT ANNOTATION DEFINITION
		INSERT IGNORE INTO annotation (`type`,`annotationType`,`source`,`release`,`description`,`date`)
		SELECT `type`,`annotationType`,`source`,`release`,`description`,`date` FROM tmp_vcf_annotation_$tmp_table_id GROUP BY `source`, `release`;

		# UPDATE VARIANT ID in TMP
		UPDATE tmp_vcf_annotation_$tmp_table_id as T
		LEFT JOIN variant as V ON (V.chr=T.chr AND V.pos=T.pos AND V.ref=T.ref AND V.alt=T.alt)
		SET T.variant_id=V.id;

		# UPDATE ANNOTATION ID in TMP
		UPDATE tmp_vcf_annotation_$tmp_table_id as T
		LEFT JOIN annotation as AD ON (AD.source=T.source AND AD.release=T.release)
		SET T.annotation_id=AD.id;

		# INSERT ANNOTATION
		INSERT IGNORE INTO variant_annotation (`variant_id`,`annotation_id`,`value`,`value_float`)
		SELECT `variant_id`,`annotation_id`,`value`,`value_float` FROM tmp_vcf_annotation_$tmp_table_id;

		";



	@queries=split(/;/,$query);
	#print "@queries\n" if $DEBUG;

	my $nb_q=0;
	foreach my $q (@queries) {
		if (trim($q)) {
			$nb_q++;
			my @values;
			my $query_handle = $DBIconnect->prepare($q);
			my $exec=$query_handle->execute(@values);
			if (!$exec) {	
				warn "#[ERROR] error in query #$nb_q (vcf.id:$id) '$q'\n";
				return return_message(undef,"","","#[ERROR] $type insertion failed! error in query '$q'");
			};#if
		};#if
	};#foreach





	if ($split_files) {
		print "NB lines: $nb_annotation_line\n" if $DEBUG;
		my $f_annotation_content_nb=0;
		my @f_annotation_file_list;
		foreach my $f_annotation_content (@f_annotation_content_list) {
			my ($f_annotation2,$tmp_file_annotation2)=tempfile();
			#print $tmp_file_annotation;
			open my $f_annotation2, '>', $tmp_file_annotation2 or die "#[ERROR] Opening File '$!' ";
			print $f_annotation2 $f_annotation_content;
			$f_annotation_file_list[$f_annotation_content_nb]=$tmp_file_annotation2;
			$f_annotation_content_nb++;
		};#foreach
		#print "SCALAR: ".scalar(@f_annotation_file_list)."\n" if $DEBUG;
		#print "@f_annotation_file_list\n" if $DEBUG;
	
		#print "TMP_VARIANT: $tmp_file_variant\n" if $DEBUG;
		#print "TMP_ANNOTATION: $tmp_file_annotation\n" if $DEBUG;

		


		#DROP TABLE IF EXISTS `tmp_vcf`; TEMPORARY
		my $query="
			CREATE TEMPORARY TABLE IF NOT EXISTS `tmp_vcf_variant_$tmp_table_id` (
			  `vcf_id` INT,
			  `variant_id` INT DEFAULT NULL,
			  `chr` TINYINT NOT NULL,
			  `pos` INT NOT NULL,
			  `ref` VARCHAR(255) NOT NULL,
			  `alt` VARCHAR(255) NOT NULL,
			  `QUAL` SMALLINT NOT NULL,
			  `FILTER` VARCHAR(45) NOT NULL,
			  `FORMAT` VARCHAR(100) NOT NULL,
			  `QUALITY` VARCHAR(100) NOT NULL,
			  `GT` VARCHAR(45) NULL,
			  `AD` VARCHAR(45) NULL,
			  `DP` SMALLINT NULL,
			  `GQ` SMALLINT NULL,
			  `PL` VARCHAR(45) NULL
			) ENGINE = InnoDB;

			CREATE TEMPORARY TABLE IF NOT EXISTS `tmp_vcf_annotation_$tmp_table_id` (
			  `vcf_id` INT,
			  `variant_id` INT DEFAULT NULL,
			  `chr` TINYINT NOT NULL,
			  `pos` INT NOT NULL,
			  `ref` VARCHAR(255) NOT NULL,
			  `alt` VARCHAR(255) NOT NULL,
			  `type` VARCHAR(45) NOT NULL,
			  `source` VARCHAR(45) NOT NULL,
			  `release` VARCHAR(45) NOT NULL,
			  `value` TEXT NOT NULL,
			  `value_float` FLOAT NULL
			) ENGINE = InnoDB;

			# LOAD DATA VARIANT
			LOAD DATA LOCAL INFILE '$tmp_file_variant'
			INTO TABLE tmp_vcf_variant_$tmp_table_id
			FIELDS TERMINATED BY '\\t'
			ENCLOSED BY ''
			ESCAPED BY '\\\\'
			LINES TERMINATED BY '\\n';

			# INSERT VARIANT
			INSERT IGNORE INTO variant (`chr`,`pos`,`ref`,`alt`)
			SELECT `chr`,`pos`,`ref`,`alt` FROM tmp_vcf_variant_$tmp_table_id;

			# UPDATE VARIANT ID in TMP
			UPDATE tmp_vcf_variant_$tmp_table_id as T
			LEFT JOIN variant as V ON (V.chr=T.chr AND V.pos=T.pos AND V.ref=T.ref AND V.alt=T.alt)
			SET T.variant_id=V.id;

			# INSERT VARIANT_VCF
			INSERT IGNORE INTO variant_vcf (`variant_id`,`vcf_id`)
			SELECT `variant_id`,`vcf_id` FROM tmp_vcf_variant_$tmp_table_id;

			# INSERT VARIANT_VCF_QUALITY
			INSERT IGNORE INTO variant_vcf_quality (`variant_id`,`vcf_id`,`QUAL`,`FILTER`,`FORMAT`,`QUALITY`,`GT`,`AD`,`DP`,`GQ`,`PL`)
			SELECT `variant_id`,`vcf_id`,`QUAL`,`FILTER`,`FORMAT`,`QUALITY`,`GT`,`AD`,`DP`,`GQ`,`PL` FROM tmp_vcf_variant_$tmp_table_id;
			";

	

		@queries=split(/;/,$query);
		#print "@queries\n" if $DEBUG;

		my $nb_q=0;
		foreach my $q (@queries) {
			if (trim($q)) {
				$nb_q++;
				my @values;
				#warn "#[QUERY] #$nb_q (vcf.id:$id) '$q'\n" if $DEBUG;
				#if ($nb_q == 3) {
				#	@values=($tmp_file_variant);
				#};#if
				#if ($nb_q == 8) {
				#	@values=($tmp_file_annotation);
				#};#if
				my $query_handle = $DBIconnect->prepare($q);
				my $exec=$query_handle->execute(@values);
				#my $exec=$query_handle->do($q);
				if (!$exec) {	
					warn "#[ERROR] error in query #$nb_q (vcf.id:$id) '$q'\n";
					return return_message(undef,"","","#[ERROR] $type insertion failed! error in query '$q'");
					#return return_message(1,"","","#[ERROR] $type insertion failed! error in query '$q'");
				};#if
			};#if
		};#foreach

	

		foreach my $f_annotation_file (@f_annotation_file_list) {
			my $query="
				# LOAD DATA ANNOTATION
				LOAD DATA LOCAL INFILE '$f_annotation_file'
				INTO TABLE tmp_vcf_annotation
				FIELDS TERMINATED BY '\\t'
				ENCLOSED BY ''
				ESCAPED BY '\\\\'
				LINES TERMINATED BY '\\n';
				";
			#my $query_handle = $DBIconnect->do($query) or die "SQL ERROR $DBI::errstr\n";
			my $query_handle = $DBIconnect->prepare($query);
			my $exec=$query_handle->execute(@values);
			if (!$exec) {	
				#warn "#[ERROR] error in query #$nb_q (vcf.id:$id) '$query'\n";
				return return_message(undef,"","","#[ERROR] $type insertion failed! error in query '$query'");
			};#if
		};#foreach

	

		my $query="
			# UPDATE VARIANT ID in TMP
			UPDATE tmp_vcf_annotation as T
			LEFT JOIN variant as V ON (V.chr=T.chr AND V.pos=T.pos AND V.ref=T.ref AND V.alt=T.alt)
			SET T.variant_id=V.id;
			";
		#my $query_handle = $DBIconnect->do($query) or die "SQL ERROR $DBI::errstr\n";
		my $query_handle = $DBIconnect->prepare($query);
		my $exec=$query_handle->execute(@values);
		if (!$exec) {	
			return return_message(undef,"","","#[ERROR] $type insertion failed! error in query '$query'");
		};#if
		my $query="
			# INSERT ANNOTATION
			INSERT IGNORE INTO annotation (`variant_id`,`type`,`source`,`release`,`value`,`value_float`)
			SELECT `variant_id`,`type`,`source`,`release`,`value`,`value_float` FROM tmp_vcf_annotation;
			";
		#my $query_handle = $DBIconnect->do($query) or die "SQL ERROR $DBI::errstr\n";
		my $query_handle = $DBIconnect->prepare($query);
		my $exec=$query_handle->execute(@values);
		if (!$exec) {	
			return return_message(undef,"","","#[ERROR] $type insertion failed! error in query '$query'");
		};#if

	};#if

	#return return_message(undef,"","","#[ERROR] $type insertion failed! UNDER DEVELOMENT!!!");
	return return_message(1,"","# $type inserted","");

}

sub insert_id {
	my $query=$_[0];	# input
	my $values=$_[1];	# input
	my $id;			# output
	#print "@values\n" if $DEBUG;

	my $query_handle = $DBIconnect->prepare($query);
	my $exec;
	if (scalar(@$values)>0 && "@$values" ne "") {
		$exec=$query_handle->execute(@$values);
	} else {
		$exec=$query_handle->execute();
	};#if
	if (!$exec) {	
		warn "#[ERROR] error in query '$query'\n";
		exit 0;
	};#if
	# Inserted ID
	$id=$query_handle->{mysql_insertid};
	if (defined $id) {
		#warn "# $type ID inserted: $id\n" if $DEBUG;
		return $id;
	};#if
	return undef;

}


sub check_id {
	my $query=$_[0];	# input
	my $values=$_[1];	# input
	my $id;			# output

	my $query_handle = $DBIconnect->prepare($query);
	#my $exec=$query_handle->execute(@values);
	#my $exec=$query_handle->execute();
	if (scalar(@$values)>0 && "@$values" ne "") {
		#print "VALUES: ".scalar(@values)." '@values'\n" if $DEBUG;
		#print "VALUES: ".scalar(@$values)." \n" if $DEBUG;
		$exec=$query_handle->execute(@$values);
	} else {
		$exec=$query_handle->execute();
	};#if
	if (!$exec) {	
		warn "#[ERROR] error in query '$query'\n";
		return undef;
	};#if
	$query_handle->bind_columns(undef, \$id);
	
	while($query_handle->fetch()) {
		return $id;
	};#while
	return undef;

}

sub return_message {
	my $return=$_[0]; #input
	my $message_output=$_[1]; #input
	my $message_verbose=$_[2]; #input
	my $message_debug=$_[3]; #input
	
	print "$message_output\n" if defined $message_output && trim($message_output) ne "";
	print "$message_verbose\n" if defined $message_verbose && trim($message_verbose) ne "" && $VERBOSE;
	print "$message_debug\n" if defined $message_debug && trim($message_debug) ne "" && $DEBUG;
	
	return $return;

}

# Global variables
#my $query_variables="SET global max_allowed_packet=135217728";
#$query_handle = $DBIconnect->do($query_variables);
$query_variables="SET autocommit=0";
$query_handle = $DBIconnect->do($query_variables);

# TMP table
if (0) {
$query="DROP TABLE IF EXISTS `tmp_vcf`";
$query_handle = $DBIconnect->prepare($query);
$query_handle->execute();
$query="CREATE TEMPORARY TABLE IF NOT EXISTS `tmp_vcf` (
		  `vcf_id` INT,
		  `variant_id` INT DEFAULT NULL,
		  `chr` TINYINT NOT NULL,
		  `pos` INT NOT NULL,
		  `ref` VARCHAR(255) NOT NULL,
		  `alt` VARCHAR(255) NOT NULL,
		  `QUAL` SMALLINT NOT NULL,
		  `FILTER` VARCHAR(45) NOT NULL,
		  `FORMAT` VARCHAR(100) NOT NULL,
		  `QUALITY` VARCHAR(100) NOT NULL,
		  `GT` VARCHAR(45) NULL,
		  `AD` VARCHAR(45) NULL,
		  `DP` SMALLINT NULL,
		  `GQ` SMALLINT NULL,
		  `PL` VARCHAR(45) NULL
		) ENGINE = MyISAM";
$query_handle = $DBIconnect->prepare($query);
$query_handle->execute();
};#if

#Lock tables
$query_lock = "LOCK TABLES project WRITE, run WRITE, manifest WRITE, sample WRITE, vcf WRITE, variant WRITE, variant_vcf WRITE, variant_vcf_quality WRITE;";
#$query_lock = "FLUSH TABLES WITH WRITE LOCK";
print "#[QUERY] $query_transaction\n" if $DEBUG;
$query_handle = $DBIconnect->prepare($query_lock);
$query_handle->execute();

#Transaction
$query_transaction = "START TRANSACTION;"; #  SET autocommit=0;
print "#[QUERY] $query_transaction\n" if $DEBUG;
$query_handle = $DBIconnect->prepare($query_transaction);
$query_handle->execute();

# Insert
my $id;
switch ($type) {

	# PROJECT
	case "project" {
		$id=insert_project(\%data);
	}#case "project"

	case "manifest" {
		$id=insert_manifest(\%data);
	}#case "project"

	# RUN
	case "run" {
		$id=insert_run(\%data);
	}#case "run"

	# SAMPLE
	case "sample" {
		$id=insert_sample(\%data);
	}#case "sample"

	# VCF
	case "vcf" {
		$id=insert_vcf(\%data);
	}#case "vcf"
	# QC
	case "QC" {
		$id=insert_qc(\%data);
	}#case "QC"
	# snapshot
	case "snapshot" {
		$id=insert_snapshot(\%data);
	}#case "snapshot"

};#type


#UNLock tables
$query_unlock = "UNLOCK TABLES;";
print "#[QUERY] $query_lock\n" if $DEBUG;
$query_handle = $DBIconnect->prepare($query_unlock);
$query_handle->execute();

if (defined $id) {
	$verbose.="# $type ID: $id\n";
	$query_transaction = "COMMIT;";

	if (defined $output_file_id) {
		open(FILE_OUTPUT, ">$output_file_id") || die "Problem to open the file '$output_file_id': $!";
		print FILE_OUTPUT $id;
		close(FILE_OUTPUT);
	};#if

} else {
	$debug.="#[ERROR] '$type' insertion failed\n";
	$output.="#[ERROR] '$type' insertion failed\n";
	$query_transaction = "ROLLBACK;";
};#if

#Transaction
#$query_transaction = "ROLLBACK;";
print "#[QUERY] $query_transaction\n" if $DEBUG;
$query_handle = $DBIconnect->prepare($query_transaction);
$query_handle->execute();

## PostProcess
################

$header.="##\n";


## OUTPUT
###########

# Header
print $header;
print $debug if $DEBUG;
print $verbose if $VERBOSE;
print $output;


__END__

