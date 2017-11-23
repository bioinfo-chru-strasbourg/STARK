#!/usr/bin/perl

=head1 NAME

check_samplemutation_vcf.pl - check for a sample ID (present in the sample set file) if the mutation is in a VCF
Require 3 files:
- functions.pl
- functions.common.pl
- parameters.common.p

=head1 USAGE

 $ARGV[0] --sample_file=<FILE> --sampleID=<STRING> --vcf_file=<FILE> [options]

=head1 REQUIRE

=over 2

=item B<--sample_file>

Sample file listing all sample with the corresponding found mutation

=item B<--sampleID>

Sample identifier present in the sample file

=item B<--VCF_file>

VCF file listing all variants fund by sequencing the sample

=back

=head1 OPTIONS

=over 2

=item B<--help>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=item B<--window>

Look around the mutation (default=10).

=item B<--show_not_found!>

Print 'NOT FOUND' if variant not found for the sample (default FALSE).

=item B<--header>

Print the header (default=false).

=back

=head1 DESCRIPTION

This program is very cool!!!

=cut

# Modules
#use 5.014;
#use strict;
#use warnings;
use Getopt::Long;
use Pod::Usage;
use DBI;
use Time::localtime;
use File::Temp qw/ tempfile tempdir tmpnam /;
use File::Basename;
use lib dirname (__FILE__);
use List::Util qw(sum);

require 'functions.pl';
require 'functions.common.pl';

# Information
my $script_release="0.9.1";
my $script_date="20130905";
my $script_author="Antony Le Bechec";

# Variables
my @chrom_array=("X", "Y", "M");	# Chromosome index
#$date=timestamp();			# Date
$basename = dirname (__FILE__);		# Base name of the folder

## PARAMETERS ##
################

## Options/Parameters default values

our %parameters = ( # 'our' to allow the visibility by imported scripts such as parameters.common.pl)
	# Main options 
	'help'      =>  	0,
	'man'       =>  	0,
	'debug'     =>  	0,
	'verbose'   =>  	0,
	# Main parameters  
	'config'    =>  	"config.ini",		# config file
	'header'    =>  	0,
	
	# Requirements
	'sample_file'   =>  	"",		# sample file
	'sampleID'    	=>  	"",		# sampleID
	'VCF_file'    	=>  	"",		# VCF file
	'window'    	=>  	10,		# window
	'show_not_found'=>	0,		# show_not_found
	'comment'    	=>  	"",		# comment
);

## Options/Parameters definition

our @options=(
	# Main options
	'help|h|?',
	'man',
	'debug',
	'verbose',
	# Main parameters  
	'config|config_file=s',
	'header=s',
	# Requirements
	'sample_file=s',	# sample file
	'sampleID=s',		# sampleID
	'VCF_file=s',		# VCF file
	# Options
	'window=i',		# window
	'show_not_found!',	# show_not_found
	'comment=s',		# comment
);

### COMMON PARAMETERS ###
#########################

require "parameters.common.pl";

### HEADER ###
##############

if ($parameters{"header"}) {
	print "#\n";
	print "# Execution Date: ".$date."\n";
	print "# Script: ".$0." (release $script_release/$script_date - author $script_author)\n";
	#print "# Configuration: Project folder: $dir_projects - Database connexion: $host:$database\n";
	print "#\n";
};#if

## FUNCTIONS ##
###############

sub sample_list { 
   	my $sample_file=$_[0];	# String to print
	my %sample_list;
	if (!-e $sample_file) {
		return 0;
	};#if
    	open(SAMPLEFILE, $sample_file) || die "Problem to open the file: $!";
	my $line=0;
	#read the file
	while(<SAMPLEFILE>) {
		$line++;
		my $content=trim($_);
		if (grep(/^#.*$/,$content)) {
			# skip
		} else {
			my @cols=split("\t",$content);
			#print $cols[0]."".$cols[2]."".$cols[3]."\n";
			$sample_list{trim($cols[0])}{"variant"}=trim($cols[3]);
			$sample_list{trim($cols[0])}{"HGVS"}=trim($cols[2]);
		};#if
	};#while	
	return %sample_list;
};

sub variant_list { 
   	my $variant_file=$_[0];	# String to print
	my %variant_list;
	if (!-e $variant_file) {
		return 0;
	};#if
    	open(VARIANTFILE, $variant_file) || die "Problem to open the file: $!";
	my $line=0;
	#read the file
	while(<VARIANTFILE>) {
		$line++;
		my $content=trim($_);
		if (grep(/^#.*$/,$content)) {
			# skip
		} else {
			my @cols=split("\t",$content);
			#print $cols[0]."\n";
			my $chrom=trim($cols[0]); $chrom =~ s/chr//i;
			my $pos=trim($cols[1]);
			my $ref=trim($cols[3]);
			my $alt=trim($cols[4]);
			$variant_list{$chrom}{$pos}{$ref}{$alt}=1;
		};#if
	};#while	
	return %variant_list;
};


### INPUT ###
#############


## MAIN ##
##########


&main; exit;

sub main {

	# Initialisation
	my %sample_list;
	my $sampleID=trim($parameters{"sampleID"});
	my $stop=0;

	# DEBUG
	if ($DEBUG || $VERBOSE) {
		print "# INPUT PARAMETERS: \n";
		print "# Sample file: ".$parameters{"sample_file"}."\n";
		print "# VCF file: ".$parameters{"VCF_file"}."\n";
		print "# SampleID: ".$parameters{"sampleID"}."\n";
		print "#\n";
	};#if
	
	# Test and load
	if (-e $parameters{"VCF_file"}) {
		%variant_list=variant_list($parameters{"VCF_file"});
		if ($DEBUG) {		
			while ((my $param_name, my $param_val) = each(%variant_list)) {
				print "   $param_name => $param_val\n";
				while ((my $param_name, my $param_val) = each(%$param_val)) {
					print "      $param_name => $param_val\n";
			    	};#while
		    	};#while
		};#if
	} else {
		print "# VCF file missing! See --VCF_file option\n";
		$stop=1;
	};#if
	if (-e $parameters{"sample_file"}) {
		%sample_list=sample_list($parameters{"sample_file"});
		# check if sample ID in sample file
		
		if ($sampleID ne "") {
			#print "@sample_list";
			#if (!grep {$_ eq $sampleID} @sample_list) {
			if (!defined $sample_list{$sampleID}) {
				#print "# Sample ID not in the sample file! See --sampleID and --sample_file options\n";
				print "# Sample ID not in the sample file!\n";
			};#if
		} else {
			print "# Sample ID missing! See --sampleID option\n";
			$stop=1;
		};#if
		
	} else {
		print "# Sample file missing! See --sample_file option\n";
		$stop=1;
	
		
	};#if
	if ($stop) {
		pod2usage(2);
	};#if

	
	# Print variant
	print "# CHECK Mutation(s) of sample '$sampleID' : ".$sample_list{$sampleID}{"variant"}." (HGVS: ".$sample_list{$sampleID}{"HGVS"}.")\n" if $VERBOSE;
	# Mutation(s)
	# clean
	$sample_list{$sampleID}{"variant"} =~ s/"//g;
	# split
	my @mutations=split("[;,]",$sample_list{$sampleID}{"variant"});
	#print "@mutations\n";
	foreach $mutation (@mutations) {
		my @mutation_split=split(":",$mutation);
		my $chrom=$mutation_split[0];
		my $pos=$mutation_split[1];
		my $ref=$mutation_split[3];
		my $alt=$mutation_split[4];
		my $found=0;
		my $found_message=$parameters{"comment"}." ";
	
		#print "@mutation_split : $chrom:$pos:$ref>$alt\n";
		if (defined $variant_list{$chrom}{$pos}{$ref}{$alt}) {
			$found_message.="[level1] $chrom:$pos:$ref>$alt ";
			print "# CHECK: FOUND! $found_message\n" if $VERBOSE;
			$found=1;
		} else {
			print "# CHECK: NOT FOUND!\n" if $VERBOSE;
		};#if

		if (!$found) {
			# +- indel length
			#print (length($ref)+length($alt)); print " \n \n";
			#print ($pos); print " \n \n";
			#print ($pos-(length($ref)+length($alt)-2)); print " \n \n";
			#print ($pos+(length($ref)+length($alt)-2)); print " \n \n";
			$count = $pos-(length($ref)+length($alt)-2);
			if ((length($ref)+length($alt)-2) > 0 #indel
				&& 	(defined $variant_list{$chrom}{$count}
					|| defined $variant_list{$chrom}{$count})) {
				#print "# FOUND! ";
				while ((my $ref, my $alts) = each(%{$variant_list{$chrom}{$count}})) {
					#print "#   $ref => $alts\n";
					while ((my $alt, my $param_val) = each(%$alts)) {
						#print "#      $alt => $param_val\n";
						$found_message.="[level2] ?$chrom:$count:$ref>$alt (distance ".($count-$pos).") ";
						print "# CHECK: FOUND! $found_message\n" if $VERBOSE;
						$found=1;
				    	};#while
			    	};#while
				while ((my $ref, my $alts) = each(%{$variant_list{$chrom}{$count}})) {
					#print "#   $ref => $alts\n";
					while ((my $alt, my $param_val) = each(%$alts)) {
						#print "#      $alt => $param_val\n";
						$found_message.="[level2] ?$chrom:$count:$ref>$alt (distance ".($count-$pos).") ";
						print "# CHECK: FOUND! $found_message\n" if $VERBOSE;
						$found=1;
				    	};#while
			    	};#while
			} else {
				print "# CHECK: NOT FOUND!\n" if $VERBOSE;
			};#if

		};#if

		if (!$found) {
			# Around
			$window=$parameters{"window"}+(length($ref)+length($alt)-2);
			#$window=$parameters{"window"};
			#print "# Look around the mutation +-$window\n"; 
			my $around_found=0;
			for ($count = $pos-$window; $count <= $pos+$window; $count++) {
				if ($count ne $pos && defined $variant_list{$chrom}{$count}) {
					#print "# FOUND Around! ";
					$around_found=1;
					while ((my $ref, my $alts) = each(%{$variant_list{$chrom}{$count}})) {
						#print "#   $ref => $alts\n";
						while ((my $alt, my $param_val) = each(%$alts)) {
							#print "#      $alt => $param_val\n";
							$found_message.="[level3] ?$chrom:$count:$ref>$alt (distance ".($count-$pos).") ";
							print "# CHECK: FOUND Around! $found_message\n" if $VERBOSE;
							$found=1;
					    	};#while
				    	};#while
				};#if
			};#for
			if (!$around_found) {
				print "# CHECK: NOT FOUND Around!\n" if $VERBOSE;
			};#if
		};#if

		if ($found) {
			print "# CHECK: FOUND! $found_message\n";
		};#if
		if (!$found && $parameters{"show_not_found"}) {
			print "# CHECK: NOT FOUND!\n";
		};#if


	};#foreach
	

}

END{
    if ($parameters{"help"}) {
	print "Script informations:\n";
	print "    $0 (release $script_release/$script_date - author $script_author)\n";
	print " \n";
    };#if
};

__END__
