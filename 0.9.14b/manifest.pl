#!/usr/bin/perl

=head1 NAME

manifest.pl - pepare manifest option

=head1 USAGE

 $ARGV[0] [options]

=head1 OPTIONS

=over 2

=item B<--help>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=item B<--manifest|manifest_file=<file>>

manifest file

=back

=head1 DESCRIPTION

This program is very cool

=cut

# Modules
#use 5.014;
#use strict;
#use warnings;
use Getopt::Long;
use Pod::Usage;
#use DBI;
use Time::localtime;
use File::Temp qw/ tempfile tempdir tmpnam /;
use File::Basename;
use lib dirname (__FILE__);
use List::Util qw(sum);

require 'functions.pl';
require 'functions.common.pl';

# Information
my $script_release="0.9.2";
my $script_date="20140612";
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
	'config'    	=>  	"config.ini",			# config file
	# Options
	'manifest'    	=>  	"",			# Manifest file
	'intervals'    	=>  	""			# Intervals file
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
	# Options
	'manifest|manifest_file=s',
	'intervals|intervals_file=s'
);

### COMMON PARAMETERS ###
#########################

require "parameters.common.pl";

### HEADER ###
##############

print "#\n";
print "# Execution Date: ".$date."\n";
print "# Script: ".$0." (release $script_release/$script_date - author $script_author)\n";
#print "# Configuration: Project folder: $dir_projects - Database connexion: $host:$database\n";
print "#\n";

## FUNCTIONS ##
###############

sub sample_list { 
    my $VCFFile=$_[0];	# VCF file
    #print "VCF file $VCFFile\n";
    open(VCFFILE, $VCFFile) || die "Problem to open the file: $!";
    my $line=0;
    my @sample_list;
    #read the file
    while(<VCFFILE>) {
	    $line++;
	    my $content=$_;
	    if (grep(/^#CHROM.*$/,$content)) {
		#print " $line\n";
		#print " $content\n";
		my @headers=split("\t",$content);
		my $format_header_found=0;
		foreach my $header (@headers) {
		    $header=trim($header);
		    if ($format_header_found) {
			#print "$header\n";
			push(@sample_list,$header);
		    };#if
		    if ($header eq "FORMAT") {
			$format_header_found=1;
		    };
		    
		};#foreach
	    };#if
    };#while
    return @sample_list;
};


### INPUT ###
#############


## MAIN ##
##########


&main; exit;

sub main {

	# OUT
    
	my $tag=0;
	my $col1=0;
	my $sep1=":";
	my $sep2="-";
	
	# INPUT
	@ManifestFiles;
	if (defined $parameters{"manifest"} && trim($parameters{"manifest"}) ne "") {
		@ManifestFiles_split=split(",",$parameters{"manifest"});
		foreach my $ManifestFile (@ManifestFiles_split) {
			if (-e $ManifestFile) {
				push(@ManifestFiles,$ManifestFile);
			};#if
		};#foreach
	};#if

	if (scalar(@ManifestFiles) == 0) {
		print "No Manifest file. See --manifest option\n";
    		pod2usage(1);
	};#if

	# OUT
	if (defined $parameters{"intervals"} && $parameters{"intervals"} ne "") {
		$IntervalsFile=$parameters{"intervals"};
	} elsif (defined $parameters{"bed"} && $parameters{"bed"} ne "") {
		$IntervalsFile=$parameters{"bed"};
		my $sep1="\t";
		my $sep2="\t";
	} else {
		print "# Error in output file! \n";
		pod2usage(1);
	};#if
	open(STDOUT, ">$IntervalsFile") || die "Problem to open the file: $!";

	foreach my $ManifestFile (@ManifestFiles) {
		
		#print "File: ".$ManifestFile."\n";
		open(ManifestFILE, $ManifestFile) || die "Problem to open the file: $!";
		my $line=0;
		my @sample_list;
		#read the file
		while(<ManifestFILE>) {
			$line++;
			my $content=trim($_);
			#if (grep(/^TargetA\tTargetB.*$/,$content)) {
			if (grep(/^Target Region Name\tTarget Region ID.*$/,$content)) {
				$tag=1;
				#$col1=3;
				$col1=5;
			} elsif (grep(/^Name\tChromosome.*$/,$content)) {
				$tag=1;
				$col1=1;
			} else {
				if ($tag) {
					#print "$content\n";
					if (trim($content) ne "") {
						my @cols=split("\t",$content);
						#print $cols[3].":".$cols[4]."-".$cols[5]."\n";
						if (trim($cols[$col1]) ne "" && trim($cols[$col1+1]) ne "" && trim($cols[$col1+2]) ne "") {					
							print $cols[$col1].":".$cols[$col1+1]."-".$cols[$col1+2]."\n";
						};#if
					};#if
				};#if
			};#if
			if (grep(/^\[.*$/,$content)) {
				$tag=0;
			};#if
		
		};#while
		close (ManifestFILE);

	};#foreach
}

END{
    if ($parameters{"help"}) {
	print "Script informations:\n";
	print "    $0 (release $script_release/$script_date - author $script_author)\n";
	print "\n";
    };#if
};

__END__
