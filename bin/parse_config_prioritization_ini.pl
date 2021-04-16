#!/usr/bin/perl
############################
# INI prioritization       #
# Author: Antony Le Béchec #
############################


### USE

use Data::Dumper;
use POSIX qw(locale_h strtod ceil);
use Time::localtime;
use Scalar::Util qw/reftype/;
use List::Util qw(sum min max);
use Getopt::Long;		# Catch Options
use Pod::Usage;			# Pod
use Time::localtime;		# Time
use Data::Dumper;		# Data
use File::Basename;		# File
use Switch;			# Switch
use File::Temp qw/ tmpnam /;
use lib dirname (__FILE__);	# Add lib in the same folder
use Scalar::Util qw(looks_like_number);
use Digest::MD5 qw(md5_hex);


### OPTIONS



our %parameters = ( #
	# Main options
	'help'			=>  	0,	# Help parameter
	'man'			=>  	0,	# Man parameter
	'release'		=>  	0,	# Release parameter
	'debug'			=>  	0,	# Debug parameter
	'verbose'		=>  	0,	# Verbose parameter
	'header'		=>  	0,	# Verbose parameter

	'config_prioritization'    	=>  	"config.prioritization.ini",	# Configuration filter file

	'format'	=>	'tsv',			# Output format (tab or dump)

	'applications'	=>	'',				# Section to print

	'show_applications'	=>	0,				# Section to print
	
	'no_header'	=>	0,				# Section to print


);

our @options=(
	# Main options
	'help|h|?',		# Help
	'man',			# Man
	'release',		# Release
	'debug',		# Debug
	'verbose',		# Verbose
	'header',		# Verbose

	'config_prioritization|config_prioritization_file|config_filter|config_filter_file=s',		# Configuration filter file	
	
	'format=s',		# format output

	'applications=s',		# format output

	'show_applications!',		# format output

	'no_header!',		# format output

);



## HELP/MAN
#############

=head1 NAME

parse_ini.pl - parse configuration prioritization file and print

=head1 DESCRIPTION

Description

=head1 BUGS

Bugs...

=head1 ACKNOWLEDGEMENTS

Thank U!

=head1 COPYRIGHT

HUS - GNU AGPL V3

=head1 AUTHOR

ALB

=head1 USAGE

$ARGV[0] [options] --inpconfig_prioritizationut=<INI_FILE> --applications=<SECTION1,SECTION2> [...]

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

=item B<--config_prioritization=<file>>

Configuration file for prioritization parameters (default 'config.prioritization.ini').

=back

=head2 OPTIONS

=over 2

=item B<--applications=<string>>

list of applications to show

=back

=item B<--show_applications>

List of applications available

=back


=cut



## Catch Options and put into parameters
GetOptions (\%parameters,@options,@common_options)  or pod2usage();


## Use example: ./parse_ini.pl --config_prioritization=/tool/config/howard/config.prioritization.ini --applications=HEMATOLOGY,SOLIDTUMOR --no_header | sort -u -f | sort -k1,2 -f | column -s$'\t' -t


## Help
if ($parameters{"help"}) {
	pod2usage(1);
	exit 0;
};#if

## Man
if ($parameters{"man"}) {
	pod2usage(-exitstatus => 0, -verbose => 2);
};#if

## Header
if ($parameters{"header"}) {
	print $header;
	exit 0;
};#if

## Release
if ($parameters{"release"}) {
	print "## Script Information\n";
	while ((my $var, $val) = each(%information)){
		print "# $var: $val\n";
	};#while
	exit 0;
};#if



### FUNCTIONS


#Function trim
#Perl trim function to remove whitespace from the start and end of the string
sub trim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}


#Function in_array
#return true if an element is in an array
sub in_array {
#parameters:
#$_[0]: The array in which you look in
#$_[1]: The element you look in the array
#$_[3]: case sensitivity (default FALSE)
#usage:  if(in_array(\@arr,'Amlan'))...
	my ($arr,$search_for,$case) = @_;
	($case eq "")?$case=1:$case=$case;
	#my %items = map {$_ => 1} @$arr;
	#return (exists($items{$search_for}))?1:0;
	my %items = map {($case)?$_:lc($_) => 1} @$arr;
	#return (exists($items{($case)?$search_for:lc($search_for)}))?1:0;
	return (exists($items{($case)?$search_for:lc($search_for)}))?1:0;
}


sub read_ini {
# Function read_ini
# Read ini file
#$_[0]: ini file
#$_[1]: header case (either NULL, "lower" or "upper")

	# Parameters
	my $ini_file=$_[0];
	my $case=lc($_[1]);
	my %ini_infos;
	my $header="";
	my $description;

     # Initialisation
	if (-e $ini_file) { #if file exist

		# Open the file
		open(FILE_PROJECT_INFOS, $ini_file) || die "Problem to open the file: $!";

		my $nb_val=0;

		# Read the file		
		while(<FILE_PROJECT_INFOS>) {

			# delete \n character
			chomp;

			# ignore blank lines
			next if /^\s*$/;

			# ignore commented lines
			next if /^\s*\;|#/;

			# check for section header
			if ( /^\s*\[([^\]]+)\]/ ) {

				# read header
				$header = $1;
				
				# case
				if ($case ne "") {
					$header=lc($header) if ($case eq "lower");
					$header=uc($header) if ($case eq "upper");
				};#if

				# remove leading and trailing white spaces
				$header =~ s/^\s*|\s*$//g;
				$ini_infos{$header} = {};

				# nb val
				$nb_val=0;

			};#if

			# check for section data
			if (/^([^=;\r\n]+)=([^;\r\n]*)/) {

				# read line
				$line=$_;

				# split line
				my @line_split1=split("",$line);
				my @line_split2=split("=",$line);
				my $var=$1;
				my $val=$2;
				# Remove comment
				@val_split=split(/([\t;])/,$val);
				$val=trim($val_split[0]);

				# case
				if ($case ne "") {
					$var=lc($var) if ($case eq "lower");
					$var=uc($var) if ($case eq "upper");
				};#if

				# add info
				my $variable=$var;
				my $value=$val;

				$var=~s/\[\]$//g;
				$value=~s/"$//g;
				$value=~s/^"//g;

				# nb val
				$nb_val++;
				
				# split val
				my @val_split1=split(":",$value);
				$criterion=$val_split1[0];
				$action=$val_split1[1];
				my $action_type="";
				$comment=$val_split1[2];
				$comment2=$val_split1[3];
				#$comment2=$val_split1[4];

				if (lc($action) eq "f" || lc($action) eq "p" ) {
					$action_type="Flag";
					$action=~s/f/FILTER/g;
					$action=~s/p/PASS/g;
				} else {
					$action_type="Score";
				};#if

				if (lc($var) eq "description") {
					$description{$header}=$value;
				};
				
				if (trim($variable) =~ /.*\[\]$/) { # Value is a hash
					
					#$ini_infos{$header}{$var}{$nb_val}=\@val_split1;
					$ini_infos{$header}{$var}{$criterion}{$action_type}{$action}="$comment";	
					#$ini_infos{$header}{$var}{$criterion}["comment"]=$comment;	
				#} else {
				#	$ini_infos{$header}{$var}=$value;
				};#if

			};#if

		};#while

		# close file
		close(FILE_PROJECT_INFOS);

		# retunr
		return %ini_infos;

	} else {

		# No ini file
		print "File '$ini_file' DOES NOT exist!\n";
		exit 1;

	};#if

}



### HEADER

if ( ! $parameters{"no_header"} ) {
	print "################################\n";
	print "# PRIORITIZATION CONFIGURATION #\n";
	print "################################\n";
	print "#[INFO] APPLICATIONS=".$parameters{"applications"}."\n";
};




### DATA

$config_prioritization_file=$parameters{"config_prioritization"};
%config_prioritization=read_ini($config_prioritization_file,1);
@applications_array=split(",",$parameters{"applications"});
$nb_applications = @applications_array;
$LIST="#APPLICATION\tANNOTATION\tCRITERION\tTYPE\tFLAG/SCORE\tCOMMENT\n";
$SECTIONS="";


# DEBUG
#print "@applications_array\n";
#print Dumper(%config_prioritization);
#print Dumper(\%description);


while ((my $application, my @filters) = each(%config_prioritization)){
	
	$SECTIONS.="$application\t".$description{$application}."\n";

	#print "application=$application\n";
	#print Dumper($config_prioritization{$application});

	if (in_array(\@applications_array,$application) || ($nb_applications == 0) ) {
	#if (in_array(\@applications_array,$application) ) {
		while ((my $field, my $filter_criteria) = each($config_prioritization{$application})){
			while ((my $criterion, my $filter_criterion_value) = each($filter_criteria)){
				while ((my $type, my $filter_criterion_value2) = each($filter_criterion_value)){
					while ((my $action, my $comment) = each($filter_criterion_value2)){
						print "$application\t$field\t$criterion\t$type\t$action\t$comment\n";
						$LIST.="$application\t$field\t$criterion\t$type\t$action\t$comment\n";
					};
				};
			};
		};
	} else {
		#print "application $application NOT found\n";
	};

};

if ($parameters{"show_applications"}) {
	print $SECTIONS;
	exit 0
};

print $LIST;

exit 0;
