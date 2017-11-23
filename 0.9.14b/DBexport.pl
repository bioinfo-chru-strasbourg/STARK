#!/usr/bin/perl
############################
# DB Export                #
# Author: Antony Le Béchec #
# Copyright: IRC           #
############################

## Main Information
#####################

our %information = ( #
	'release'	=>  	"0.9.3beta",	# Release
	#'beta'		=>  	"beta",		# Man parameter
	'date'		=>  	"20150203",	# Release parameter
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
use Digest::MD5 qw(md5_hex);

require "functions.inc.pl";	# Common functions


## HELP/MAN
#############

=head1 NAME

VCFprioritization.pl - Prioritization of variants from annotations

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

$ARGV[0] [options] --input=<VCF> --output=<VCF> [...]

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

=item B<--config_filter|config_filter_file=<file>>

Configuration file for filter parameters (default 'config.filter.ini').

=back

=head2 INPUT/OUPUT

=over 2

=item B<--output=<file>>

VCF Output file

=back

=head2 OPTIONS

=over 2

=item B<--vcf_id=<string>>

List of VCF id in the DB

=item B<--variant_id=<string>>

List of variant id in the DB

=back

=cut


## Parameters
###############

## Parameters default values

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
	'config_filter'    	=>  	"config.filter.ini",		# Configuration filter file
	# Input
	#'input'		=>	undef,		# Input VCF file
	#'filter'	=>	'ALL',	# Filter sources
	# Output
	'output'	=>	undef,		# Output VCF file
	# Options
	'vcf_id'	=>	undef,		# VCF ID
	'variant_id'	=>	undef,		# VARIANT ID
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
	'config_filter|config_filter_file=s',		# Configuration filter file
	# Input
	#'input|input_file=s',	# Input file
	#'filter=s',		# Filter
	# Output
	'output|output_file=s',	# Output file
	# Options
	'vcf_id=s',		# VCF ID
	'variant_id=s',		# VARIANT ID
);

## Catch Options and put into parameters
GetOptions (\%parameters,@options,@common_options)  or pod2usage();

## Main parameters
$date=timestamp();
$basename = dirname (__FILE__);
$scriptname = basename (__FILE__);

my $script_information=$scriptname." (".join(q{, }, map{qq{$_:$information{$_}}} keys %information).")";

## Header
$header="##\n";
#$header.="## Script: ".$information{"script"}." (release ".$information{"release"}."/".$information{"date"}.")\n";
$header.="## Script: $script_information\n";
$header.="## Excecution Date: ".$date."\n";
$header.="##\n";


## Parameters
###############

require "parameters.inc.pl";	# Parameters

## PrePorcessing
##################


## Output file
my $output_file;
#my $output_filename_pattern="prioritized";
if ($parameters{"output"} eq "") {
	print "# ERROR: NO output file\n";
	pod2usage();
	exit 1;
} else {
	$output_file=$parameters{"output"};
	if (trim($output_file) eq "") {
		$output_file=$input_file;
		$output_file =~ s/\.vcf$/\.$output_filename_pattern\.vcf/g; #/\.vcf$/\.output\.vcf/;
	};#if
	$header.="## Output VCF file: $output_file\n";
};#if

## Options

if (!defined $parameters{"vcf_id"} && !defined $parameters{"variant_id"}) {
	print "# ERROR: vcf_id or variant_id option\n";
	pod2usage();
	exit 1;
};#if

if (!defined $parameters{"vcf_id"}) {
	#print "# ERROR: vcf_id option '".$parameters{"output"}."'\n";
	#pod2usage();
	#exit 1;
} else {
	$vcf_id_input="'".join("','",split(",",$parameters{"vcf_id"}))."'";
};#fi
#print "$vcf_id_input\n" if $DEBUG;
#exit 0;

if (!defined $parameters{"variant_id"}) {
	#print "# ERROR: variant_id option '".$parameters{"output"}."'\n";
	#pod2usage();
	#exit 1;
} else {
	$variant_id_input="".join(",",split(",",$parameters{"variant_id"}))."";
};#fi
#print "$vcf_id_input\n" if $DEBUG;
#exit 0;


## DEBUG
##########

#print Dumper(\%config) if $DEBUG;
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
	print "# VCF id: $vcf_id_input\n";
	print "# Variant id: $variant_id_input\n";
	print "# Output file: $output_file\n";
	#print "# Output format: @output_formats\n";
	#print "# Filter: $filter\n";
	#print "# Filters: ".(map { $_ } keys %filter_list)."\n";# ;
	print "#\n";
};#if




# Open
#open(VCF, $parameters{"input"}) || die "Problem to open the file '".$parameters{"input"}."': $!";
my %variant_list;
my %sample_list;
my %annotation_description; # unused
my %vcf_header;

if (1) {

	if (trim($vcf_id_input) ne "") {
		$vcf_id_query=" AND vcf.id IN ($vcf_id_input) ";
	};#if
	if (trim($variant_id_input) ne "") {
		$variant_id_query=" AND variant.id IN ($variant_id_input) ";
	};#if

	my $query="SELECT vcf.id, vcf.ref,
			variant.id, variant.chr, variant.pos, variant.ref, variant.alt,
			variant_vcf.rank,
			variant_vcf_quality.QUAL, variant_vcf_quality.FILTER, variant_vcf_quality.FORMAT, variant_vcf_quality.QUALITY,
			variant_vcf_quality.GT, variant_vcf_quality.AD, variant_vcf_quality.DP, variant_vcf_quality.GQ, variant_vcf_quality.PL,
			variant_annotation.value, annotation.type, annotation.source, annotation.release, annotation.description, annotation.date
			FROM variant_vcf
			INNER JOIN vcf ON (vcf.id=variant_vcf.vcf_id)
			INNER JOIN variant ON (variant.id=variant_vcf.variant_id)
			INNER JOIN variant_annotation ON (variant.id=variant_annotation.variant_id)
			INNER JOIN annotation ON (variant_annotation.annotation_id=annotation.id)
			INNER JOIN variant_vcf_quality ON (variant_vcf_quality.variant_id=variant_vcf.variant_id
							   AND variant_vcf_quality.vcf_id=variant_vcf.vcf_id)
			WHERE annotation.type <> ''
			  $vcf_id_query
			  $variant_id_query
			ORDER BY variant.chr ASC, variant.pos ASC, variant_vcf.rank DESC
			  ";
	# INNER JOIN annotation ON (annotation.variant_id=variant.id)
	# AND vcf.id IN ($vcf_id_input)
	#print "$query\n" if $DEBUG;
	#exit 0;

	my @values;
	
	my $query_handle = $DBIconnect->prepare($query);
	if (scalar(@$values)>0 && "@values" ne "") {
		$exec=$query_handle->execute(@values);
	} else {
		$exec=$query_handle->execute();
	};#if
	if (!$exec) {	
		warn "#[ERROR] error in query '$query'\n";
	};#if
	$query_handle->bind_columns(undef, \$vcf_id, \$vcf_ref,
					\$variant_id, \$variant_chr, \$variant_pos, \$variant_ref, \$variant_alt,
					\$variant_vcf_rank,
					\$variant_vcf_quality_QUAL, \$variant_vcf_quality_FILTER, \$variant_vcf_quality_FORMAT, \$variant_vcf_quality_QUALITY,
					\$variant_vcf_quality_GT, \$variant_vcf_quality_AD, \$variant_vcf_quality_DP, \$variant_vcf_quality_GQ, \$variant_vcf_quality_PL,
					\$variant_annotation_value, \$annotation_type,  \$annotation_source,  \$annotation_release,  \$annotation_description,  \$annotation_date, 
					);
	
	while($query_handle->fetch()) {
		#"CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HT0100"

		#print "$variant_vcf_quality_FORMAT\n" if $DEBUG;
		#print "$variant_vcf_quality_QUALITY\n\n" if $DEBUG;
		
		#print Dumper(\%sample_quality) if $DEBUG;

		my @global_format=split(":",trim($annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"FORMAT"}));

		# global FILTER
		# ADD QUAL and FILTER ?
		# $variant_vcf_quality_FORMAT
		#if (index($variant_vcf_quality_FORMAT,"FILTER") == -1) {
		print "";
		if ($variant_vcf_quality_FORMAT!=~ /FILTER/) {
			$variant_vcf_quality_FORMAT.=":FILTER";
			$variant_vcf_quality_QUALITY.=":$variant_vcf_quality_FILTER";
		};#if
		if ($variant_vcf_quality_FORMAT!=~ /QUAL/) {
			$variant_vcf_quality_FORMAT.=":QUAL";
			$variant_vcf_quality_QUALITY.=":$variant_vcf_quality_QUAL";
		};#if
		#foreach my $Q (keys %sample_quality) { # @global_format) {
		foreach my $Q (split(":",$variant_vcf_quality_FORMAT)) { # @global_format) {
			#print $Q;
			my %global_format_hash=map{$_=>1} @global_format;
			if (!exists($global_format_hash{$Q})) {
				push @global_format, trim($Q);
			};
			#print @global_format;
			#print "\n";
		};#foreach
		#print "GLOBAL_FORMAT".join(":",@global_format); print "\n";
		#print "GLOBAL_QUAL".$variant_vcf_quality_QUALITY; print "\n";
		
		my %sample_quality;
		@sample_quality{split(":",trim($variant_vcf_quality_FORMAT))}=split(":",trim($variant_vcf_quality_QUALITY));

		# Main QUAL FILTER
		# QUAL
		my $main_QUAL=$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"QUAL"};
		#print "$main_QUAL+0>$variant_vcf_quality_QUAL+0\n";
		if ($main_QUAL+0>$variant_vcf_quality_QUAL+0 || $main_QUAL=="") {
			#print "TEST main_QUAL OK\n";
			$main_QUAL=$variant_vcf_quality_QUAL;
		};#if
		my $main_FILTER=$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"FILTER"}."";
		#print "variant_vcf_quality_FILTER=$variant_vcf_quality_FILTER\n";
		#if ($variant_vcf_quality_FILTER ne "PASS" && $variant_vcf_quality_FILTER ne "" && $variant_vcf_quality_FILTER ne "." ) {
			#print "TEST main_FILTER OK\n";
			if (grep /$variant_vcf_quality_FILTER/, split(",",$main_FILTER)) {
				#print "TESTOK";
			} else {
				$main_FILTER.=(($main_FILTER ne "")?",":"").$variant_vcf_quality_FILTER;
			};#if
		#};#if
		#print "main_QUAL:$main_QUAL\n";
		#print "main_FILTER:$main_FILTER\n";


		my $sample_vcf_ref="$vcf_id.$vcf_ref";
		$variant_list{$variant_id}=1;
		$sample_list{$sample_vcf_ref}=$vcf_id;
		$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"CHROM"}=$variant_chr; #"{$variant_chr}{$variant_pos}{$variant_ref}{$variant_alt}";
		$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"POS"}=$variant_pos; #"{$variant_chr}{$variant_pos}{$variant_ref}{$variant_alt}";
		$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"ID"}=$variant_id; #"{$variant_chr}{$variant_pos}{$variant_ref}{$variant_alt}";
		$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"REF"}=$variant_ref; #"{$variant_chr}{$variant_pos}{$variant_ref}{$variant_alt}";
		$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"ALT"}=$variant_alt; #"{$variant_chr}{$variant_pos}{$variant_ref}{$variant_alt}";
		#$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"QUAL"}=$variant_vcf_quality_QUAL; #"{$variant_chr}{$variant_pos}{$variant_ref}{$variant_alt}";
		#$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"FILTER"}=$variant_vcf_quality_FILTER; #"{$variant_chr}{$variant_pos}{$variant_ref}{$variant_alt}";
		$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"QUAL"}=$main_QUAL; #"{$variant_chr}{$variant_pos}{$variant_ref}{$variant_alt}";
		$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"FILTER"}=$main_FILTER; #"{$variant_chr}{$variant_pos}{$variant_ref}{$variant_alt}";
		$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"INFO"}=""; #"{$variant_chr}{$variant_pos}{$variant_ref}{$variant_alt}";
		#$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"FORMAT"}=$variant_vcf_quality_FORMAT; #"{$variant_chr}{$variant_pos}{$variant_ref}{$variant_alt}";
		$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"FORMAT"}=join(":",@global_format); #"{$variant_chr}{$variant_pos}{$variant_ref}{$variant_alt}";
		$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"SAMPLES_FORMAT"}{$sample_vcf_ref}=$variant_vcf_quality_FORMAT; #}{$variant_alt}";
		$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"SAMPLES"}{$sample_vcf_ref}=$variant_vcf_quality_QUALITY; #"{$variant_chr}{$variant_pos}{$variant_ref}{$variant_alt}";
		$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"SAMPLES_HASH"}{$sample_vcf_ref}=\%sample_quality; #"{$variant_chr}{$variant_pos}{$variant_ref}{$variant_alt}";
		$annotation_outputT{$variant_chr}{$variant_pos}{$variant_vcf_rank}{$variant_ref}{$variant_alt}{"INFOS"}{$annotation_source}=$variant_annotation_value; #"{$variant_chr}{$variant_pos}{$variant_ref}{$variant_alt}";

		$vcf_header{"INFO"}{$annotation_source}{"Number"}=".";
		$vcf_header{"INFO"}{$annotation_source}{"Type"}="String";
		$vcf_header{"INFO"}{$annotation_source}{"Description"}="\"".((trim($annotation_description) eq "")?"unknown":trim($annotation_description))."\"";
		$vcf_header{"INFO"}{$annotation_source}{"Release"}="$annotation_release";
		$vcf_header{"INFO"}{$annotation_source}{"Date"}="$annotation_date";
		$vcf_header{"INFO"}{$annotation_source}{"AnnotationType"}="$annotation_type";

		# QUAL calculation TODO
		# FILTER calculation TODO
		#print "$vcf_id\t$vcf_ref\t
		#	$variant_id\t$variant_chr\t$variant_pos\t$variant_ref\t$variant_alt
		#	$variant_vcf_quality_QUAL\t$variant_vcf_quality_FILTER\t$variant_vcf_quality_FORMAT\t$variant_vcf_quality_QUALITY
		#	$variant_vcf_quality_GT\t$variant_vcf_quality_AD\t$variant_vcf_quality_DP\t$variant_vcf_quality_GQ\t$variant_vcf_quality_PL\n" if $DEBUG;
	};#while

};#if


# READ HASH Variants
while ( my ($chr, $poss) = each(%annotation_outputT) ) {
	while ( my ($pos, $ranks) = each(%{$poss}) ) {
		@rank_array=keys %{$ranks};
		#print "@rank_array \n"; print scalar(@rank_array); print "\n";

		# Multiallele
		if (scalar(@rank_array) > 1) {
			my $ref_merged="";
			my $alt_merged="";
			my %annotation_merged;
			foreach my $rank (sort { $ranks{$a} <=> $ranks{$b} or $a cmp $b } keys %$ranks) {
				while ((my $ref, my $alts) = each(%{$annotation_outputT{$chr}{$pos}{$rank}})){
					$ref_merged.=($ref_merged ne "")?",$ref":"$ref";
					while ( my ($alt, $variant_annotations) = each(%{$alts}) ) {
						#print "$chr\t$pos\t$rank\t$ref\t$alt\n";
						$alt_merged.=($alt_merged ne "")?",$alt":"$alt";
						$annotation_merged=$variant_annotations;
						#$annotation_output{$chr}{$pos}{$ref}{$alt}=$variant_annotations;
					};#foreach

				};#foreach
			};#foreach
			#print "ref='$ref_merged' alt='$alt_merged'\n";

			# biallele only
			$previous_ref=(split(",",$ref_merged))[0];
			$variant_ref=(split(",",$ref_merged))[1];
			$previous_alt=(split(",",$alt_merged))[0];
			$variant_alt=(split(",",$alt_merged))[1];

			#print "# $previous_ref>$previous_alt $variant_ref>$variant_alt\n";

			if (length($previous_ref)>length($variant_ref)) {
				while (length($variant_alt)>0 && length($variant_ref)>0) {
					$variant_alt=substr($variant_alt,1);
					$variant_ref=substr($variant_ref,1);
					#print "$pos ref'$variant_ref' alt'$variant_alt'\n";
				};#if	
				$variant_alt=$previous_ref.$variant_alt;
				$variant_ref=$previous_ref.$variant_ref;
				
			} elsif (length($previous_ref)<length($variant_ref)) {
				while (length($previous_alt)>0 && length($previous_ref)>0) {
					$previous_alt=substr($previous_alt,1);
					$previous_ref=substr($previous_ref,1);
					#print "$pos ref'$previous_ref' alt'$previous_alt'\n";
				};#if	
				$previous_alt=$variant_ref.$previous_alt;
				$previous_ref=$variant_ref.$previous_ref;
				$previousN_ref=$previous_ref;
				$previousN_alt=$previous_alt;
				#print "PreviousN\tref $previous_ref\talt $previous_alt\n";
			};# if
			
			#print "# $previous_ref>$previous_alt $variant_ref>$variant_alt\n";
			#print "# NEW $previous_ref>$previous_alt,$variant_alt\n";
			

			$annotation_output{$chr}{$pos}{$previous_ref}{"$previous_alt,$variant_alt"}=$annotation_merged;
			#Dumper(\%{$annotation_merged});



		} else {
			foreach my $rank (sort { $ranks{$a} <=> $ranks{$b} or $a cmp $b } keys %$ranks) {
				while ((my $ref, my $alts) = each(%{$annotation_outputT{$chr}{$pos}{$rank}})){
					while ( my ($alt, $variant_annotations) = each(%{$alts}) ) {
						#print "$chr\t$pos\t$rank\t$ref\t$alt\n";
						$annotation_output{$chr}{$pos}{$ref}{$alt}=$variant_annotations;
					};#foreach
				};#foreach
			};#foreach
		};#if
	};#while
};

#print Dumper(\%annotation_output);
if ($DEBUG) {
	#print Dumper(\%annotation_output) if $DEBUG;
	#print Dumper(\%annotations) if $DEBUG;filter_list
	#print Dumper(\%filter_list) if $DEBUG;
};#DEBUG

#exit 0;



my $line=0;
my $data_header=0;
my $variant_group=0;
my $variant_index=0;
my $header_VCF;
my $header_VCF_samples;
my $nb_sample=0;
my %VCF_header;
my $vcf;

# READ VCF FILE from DB
#my %annotation_output=read_vcf($tmp_file_variant);



if ($DEBUG) {
	print Dumper(\%annotation_output) if $DEBUG;
	#print Dumper(\%annotations) if $DEBUG;filter_list
	#print Dumper(\%filter_list) if $DEBUG;
};#DEBUG

#exit 0;

my $PZHeader;

my $vcf_header_variant_line="#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT";
while ( my ($sample_vcf_ref, $vcf_id) = each(%sample_list) ) {
	my $sep="\t";
	$vcf_header_variant_line.="$sep$sample_vcf_ref";
};#foreach
$vcf_header_variant_line.="\n";

# CHROM replace
my %replacements = ("X" => "-1", "Y" => "-2", "M" => "-3");
my %replacements_reverse = reverse %replacements;


# READ HASH Variants
while ( my ($chr, $poss) = each(%annotation_output) ) {
while ( my ($pos, $refs) = each(%{$poss}) ) {
while ( my ($ref, $alts) = each(%{$refs}) ) {
while ( my ($alt, $variant_annotations) = each(%{$alts}) ) {

	#"CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HT0100"

	# VARIANT
	#print "$chr:$pos|$ref>$alt\n" if $DEBUG;

	# CHROM
	my $chrom="chr$chr";
	#$chrom =~ s/\.vcf$/\.$output_filename_pattern\.vcf/g; #/\.vcf$/\.output\.vcf/;
	$chrom =~ s/(@{[join "|", keys %replacements_reverse]})/$replacements_reverse{$1}/g;

	# QUAL
	my $ID=$$variant_annotations{"ID"};
	my $QUAL=$$variant_annotations{"QUAL"};
	my $FILTER=$$variant_annotations{"FILTER"};
	my $FORMAT=$$variant_annotations{"FORMAT"};
	#print "\t$ID\t$QUAL\t$FILTER\t$FORMAT\n" if $DEBUG;

	# SAMPLES
	#print Dumper(\%variant_annotations_list) if $DEBUG;	
	my $samples=$$variant_annotations{"SAMPLES"};
	my $samples_hash=$$variant_annotations{"SAMPLES_HASH"};
	my $samples_field;
	while ( my ($sample_name, $sample_quality) = each(%$samples) ) {
		#print "\tA $sample_name=$sample_quality\n" if $DEBUG;
	};#foreach
	while ( my ($sample_vcf_ref, $vcf_id) = each(%sample_list) ) {

		my $sample_quality="0/0"; #(exists $$variant_annotations{"SAMPLES"}{$sample_vcf_ref})?$$variant_annotations{"SAMPLES"}{$sample_vcf_ref}:"0/0";

		if ($DEBUG) {
			print "$FORMAT\n";
			print $$variant_annotations{"SAMPLES"}{$sample_vcf_ref};
			print "\n";
			while ( my ($Q, $QV) = each(%{$$samples_hash{$sample_vcf_ref}}) ) {
			#	while ( my ($Q, $QV) = each(%$QV) ) {
					print " $Q=>$QV\n";
			#	};#while
			};#while
		};#if

		#if (exists $$variant_annotations{"SAMPLES_HASH"}{$sample_vcf_ref}) {
		if (exists $$samples_hash{$sample_vcf_ref}) {
			#print "OK $sample_vcf_ref\n";
			#print $$samples_hash{$sample_vcf_ref}{"GQ"};
			my @QV;
			
			for my $Q (split(":",$FORMAT)) {
				print $$samples_hash{$sample_vcf_ref}{$Q} if $DEBUG;
				if (exists $$samples_hash{$sample_vcf_ref}{$Q} && trim($$samples_hash{$sample_vcf_ref}{$Q}) ne "") {
					print "PUSH\t$Q\t".$$samples_hash{$sample_vcf_ref}{$Q}."\t@QV\n" if $DEBUG;
					push @QV, $$samples_hash{$sample_vcf_ref}{$Q};
				} else {
					if ($Q eq "GT") {
						print "PUSH\t$Q\t0/0\t@QV\n" if $DEBUG;
						push @QV, "0/0";
					} else {
						print "PUSH\t$Q\t.\t@QV\n" if $DEBUG;
						push @QV, ".";
					};#if
				};#if
			};#for
			$sample_quality=join(":",@QV);
			#print join(":",@QV); print "\n";
			#while ( my ($Q, $QV) = each(%{$samples_hash{$sample_vcf_ref}}) ) {
			while ( my ($Q, $QV) = each(%{$$samples_hash{$sample_vcf_ref}}) ) {
			#	while ( my ($Q, $QV) = each(%$QV) ) {
					#print " $Q=>$QV\n";
			#	};#while
			};#while
		};#if


		#print "\t$sample_vcf_ref=$sample_quality\n" if $DEBUG;
		#my $sep=(trim($samples_field) ne "")?"\t":"";
		my $sep="\t";
		$samples_field.="$sep$sample_quality";
	};#foreach

	# Main QUAL/FILTER
	#if (scalar(keys(%sample_list))) { #print "NBSAMPLE"; print scalar(keys(%sample_list)); print "|\n";
	#	$QUAL=".";
	#	$FILTER=".";
	#};#if
	if ($QUAL eq "") {
		$QUAL="0";
	};#if
	if ($FILTER eq "") {
		$FILTER=".";
	};#if

	# INFOS
	my $infos=$$variant_annotations{"INFOS"};
	my $infos_field;
	while ( my ($annotation_name, $annotation_value) = each(%$infos) ) {
		#print "\t$annotation_name=$annotation_value\n" if $DEBUG;
		my $sep=(trim($infos_field) ne "")?";":"";
		$infos_field.="$sep$annotation_name=$annotation_value";
	};#while

	#print "SAMPLES: $samples_field\n" if $DEBUG;
	#print "INFOS  : $infos_field\n" if $DEBUG;

	$vcf.="$chrom\t$pos\t$ID\t$ref\t$alt\t$QUAL\t$FILTER\t$infos_field\t$FORMAT$samples_field\n";

}
}
}
}
#print "$vcf\n" if $DEBUG;
#print "$vcf_header_variant_line\n" if $DEBUG;
#print Dumper(\%vcf_header);

# VCF HEADER
my $vcf_header_content_main="##fileformat=VCFv4.1\n";
$vcf_header_content_main.="##fileDate=$date\n"; 
$vcf_header_content_main.="##source=$script_information\n"; 

my $vcf_header_content="";
foreach my $type (sort { $vcf_header{$a} <=> $vcf_header{$b} or $a cmp $b } keys %vcf_header) {
	foreach my $ID (sort { $vcf_header{$type}{$a} <=> $vcf_header{$type}{$b} or $a cmp $b } keys %{$vcf_header{$type}}) {
		my $line="##$type=<ID=$ID";
		while ((my $info_var, my $info_val) = each(%{$vcf_header{$type}{$ID}})){
			if (trim($info_val) ne "") {
				if ($info_var eq "ID" || $info_var eq "Number" || $info_var eq "Type" || $info_var eq "Description") {
				$line.=",$info_var=$info_val";
				};#if
			};#if
		};#while
		$line.=">\n";
		$vcf_header_content.=$line;
	};#foreach
};#foreach

open(FILE_OUTPUT, ">$output_file") || die "Problem to open the file '$output_file': $!";
print FILE_OUTPUT $vcf_header_content_main.$vcf_header_content.$vcf_header_variant_line.$vcf;
close(FILE_OUTPUT);

exit 0;

## Write VCF
#open the files
open(FILE_OUTPUT, ">$output_file") || die "Problem to open the file '$output_file': $!";
open(FILE_INPUT, $input_file) || die "Problem to open the file '$input_file': $!";
$line=0;
$line_info_read=0;
#read the file
while(<FILE_INPUT>) {

	# init
	chomp; #delete \n character
	$line++;
	$line_content=$_;
	@line_content_split=split("\t",$line_content);

	if (substr($line_content,0,1) eq "#") { # VCF Header
		#check the line in order to insert the description of the new info into the VCF INFO column 
		if (substr($line_content,0,7) eq "##INFO=") {
			$line_info_read=1;
		} elsif ($line_info_read==1) {
			#print FILE_INPUT_annotated "$vcf_header";
			$line_info_read=0;
		};#if
		if (substr($line_content,0,6) eq "#CHROM") {
			# Add Header
			print FILE_OUTPUT $PZHeader;
		};#if
		#keep the line
		print FILE_OUTPUT "$line_content\n";
	} else {
		my $chr=$line_content_split[0];
		my $pos=$line_content_split[1];
		my $ref=$line_content_split[3];
		my $alt=$line_content_split[4];
		my $annotation_input_line=$line_content_split[7];
		
		%{$annotation_input{$chr}{$pos}{$ref}{$alt}} = map{split /=/, $_, 2}(split /;/, $annotation_input_line);

		if ($DEBUG && 0) {
			my @test=split /;/, $annotation_input_line;
			foreach my $t (@test) {
				#print "$t\n" if $DEBUG;	
				my @test2=split /=/, $t, 2;
				foreach my $t2 (@test2) {
					#print "   $t2\n" if $DEBUG;	
					
				};#foreach
			};#foreach
			#print "@test\n" if $DEBUG;
			#print Dumper(\%{$annotation_input{$chr}{$pos}{$ref}{$alt}}) if $DEBUG;
		};#if

		while ((my $annotation_name, my $annotation_result) = each(%{$annotation_output{$chr}{$pos}{$ref}{$alt}})){
		#while ((my $annotation_name, my $annotation_result) = each(%{$annotation_input{$chr}{$pos}{$ref}{$alt}})){
			#my $info_field=$annotations{$chrom}{$pos}{$ref}{$alt};
			
			if ($DEBUG && 0) {
				#print "$annotation_name=$annotation_result\n" if $DEBUG;
				if ($annotation_result =~ /;/) {
					#print "ERROR annotation value for '$annotation_name': $annotation_result\n";
				};#
			};#if


			if (!exists $annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name}
				|| $annotation_result ne $annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name}) {
				#$line_content_split[7].=";$annotation_name=$annotation_result";
				my $original_result=trim($annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name});
				my $sep=((trim($annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name}) eq "")?"":",");
				$annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name}="$original_result$sep$annotation_result";
			};#if
		};#while
		#
		#my $str = join(", ", map { "$_ X $hash{$_}" } keys %hash);
		#print Dumper(\%{$annotation_input{$chr}{$pos}{$ref}{$alt}}) if $DEBUG;
		my $variant_annotation=join(";", map { "$_=$annotation_input{$chr}{$pos}{$ref}{$alt}{$_}" } keys %{$annotation_input{$chr}{$pos}{$ref}{$alt}});
		#print "$variant_annotation\n" if $DEBUG;
		$line_content_split[7]=$variant_annotation;
		$variant_line_join=join("\t",@line_content_split);
		print FILE_OUTPUT "$variant_line_join\n";
	};#if

	

};#while
close(FILE_INPUT);
close(FILE_OUTPUT);


## PostProcess
################

$header.="##\n";


## OUTPUT
###########

# Header
print $header;
print $debug;
print $verbose;
print $output;


__END__

