#!/usr/bin/perl

=head1 NAME

template.pl - template perl script
Require 3 files:
- functions.pl
- functions.common.pl
- parameters.common.pl

=head1 USAGE

 $ARGV[0] [options]

=head1 OPTIONS

=over 2

=item B<--help>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

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
#use DBI;
use Time::localtime;
use File::Temp qw/ tempfile tempdir tmpnam /;
use File::Basename;
use lib dirname (__FILE__);
use List::Util qw(sum);

require 'functions.pl';
require 'functions.common.pl';

# Information
my $script_release="0.9";
my $script_date="20130711";
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
	# Options
	'command'    =>  	"",	# command to launch
	'endfile'    =>  	"",	# file generated at the end of the script
	'logfile'    =>  	""	# global log file
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
	'command=s',
	'endfile=s',
	'logfile=s'
);

### COMMON PARAMETERS ###
#########################

require "parameters.common.pl";

### HEADER ###
##############

#print "#\n";
#print "# Execution Date: ".$date."\n";
#print "# Script: ".$0." (release $script_release/$script_date - author $script_author)\n";
#print "# Configuration: Project folder: $dir_projects - Database connexion: $host:$database\n";
#print "#\n";

## FUNCTIONS ##
###############

sub print_fc { 
    my $to_print=$_[0];	# String to print
    print $to_print;
};

sub command_fc_old {
	my $command=$_[0];	# Command to run
	my $logFile=$_[1];	# Log file
	my $comment=$_[2];	# Comment
	if (trim($logFile) ne "") {
		if (trim($comment) ne "") {
			system "echo \"# $comment\" >> $logFile";
		};#if
		system "echo \"#COMMAND[".timestamp()."]: $command\" >> $logFile";
	};#if
	print "$command\n" if $DEBUG;
	my $start=time;	
	system $command;
	my $end=time;
	if (trim($logFile) ne "") {
		system "echo \"#EXECTIME = $time sec\" >> $logFile";
	};#if
	return ($end-$start);
	
}


### INPUT ###
#############


## MAIN ##
##########


&main; exit;

sub main {

	if ($parameters{"command"} ne "" && $parameters{"endfile"} ne "") {
		my $command=$parameters{"command"};
		my $logFile=$parameters{"logfile"};
		my $endfile=$parameters{"endfile"};
		my $endfileLog=$endfile.".log";
		
		#$time=command_fc($command,$logFile);
		system "echo \"#COMMAND[".timestamp()."]: $command\" > $endfileLog";
		#$time=command_fc($command,$endfileLog);
		$command = "$command  >> $endfileLog 2>&1";
		$time=command_fc($command);
		system "echo \"#EXECTIME = $time sec\" >> $endfileLog";
		#system $command;
		#system "touch ".$parameters{"endfile"};
		system "cat $endfileLog > $endfile";
		if (trim($logFile) ne "") {
			system "cat $endfileLog >> $logFile";
		};#if
		system "rm $endfileLog";
	};#if

}

END{
    if ($parameters{"help"}) {
	print "Script informations:\n";
	print "    $0 (release $script_release/$script_date - author $script_author)\n";
	print "\n";
    };#if
};

__END__
