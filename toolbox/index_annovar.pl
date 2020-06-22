#!/usr/bin/perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: ba2461d35c1c0732560fce2c19319479f5232d60 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2012-10-23 23:32:05 -0700 (Tue, 23 Oct 2012) $';


our ($verbose, $help, $man);
our ($dbfile);
our ($filetype, $bin, $outfile, $skipsort, $commentfile);

GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'filetype=s'=>\$filetype, 'bin=i'=>\$bin, 'outfile=s'=>\$outfile,
	'skipsort'=>\$skipsort, 'commentfile=s'=>\$commentfile) or pod2usage ();
	
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error");

($dbfile) = @ARGV;

$filetype ||= 'A';
$filetype =~ m/^[ABC]$/ or pod2usage ("Error in argument: the -filetype argument can be only 'A' or 'B' or 'C'");
$bin ||= 1000;
$outfile ||= "$dbfile.newdb";

print STDERR "NOTICE: the bin size is set as $bin (use -bin to change this)\n";
print STDERR "NOTICE: Two output files will be generated for use by ANNOVAR: $outfile and $outfile.idx (use -outfile to override)\n";


if (not $skipsort) {
	#step 1: generate the new output file
	print STDERR "NOTICE: Running the first step of indexing (generating $outfile) ...\n";
	
	if ($dbfile eq $outfile) {
		die "Error: your -outfile is identical to input file. Use -skipsort if you are sure that inputfile is sorted\n";
	}
	
	my $command;
	#$command = "echo -n > $outfile";	#create a new empty file
	#system ($command);
	
	if (defined $commentfile) {
		$command = qq{grep -P '^#' $commentfile > $outfile};		#keep the comment lines in the output file
		system ($command);
		print STDERR "NOTICE: Adding comments from commentfile by <$command>\n";
	} else {
		$command = qq{grep -P '^#' $dbfile > $outfile};		#keep the comment lines in the output file
		system ($command);
	}
	
	for my $i (1 .. 22, 'X', 'Y', 'M', 'MT') {	
		if ($filetype eq 'A') {
			$command = qq#grep -P '^(chr)?$i\\t\\d+' $dbfile | sort -n -k 2 >> $outfile#;
		} elsif ($filetype eq 'B') {
			$command = qq#grep -P '^\\w+\\t(chr)?$i\\t\\d+' $dbfile | sort -n -k 3 >> $outfile#;
		} elsif ($filetype eq 'C') {
			$command = qq#grep -P '^\\w+\\t\\w+\\t(chr)?$i\\t\\w+\\t\\d+' $dbfile | sort -n -k 5 >> $outfile#;
		}
		$verbose and print STDERR "NOTICE: Running command: $command\n";
		system ($command);
	}
} else {
	my $command;
	#step 1: generate the new output file
	if ($dbfile ne $outfile) {
		if ($commentfile) {
			print STDERR "NOTICE: Running the first step of indexing (combining $commentfile and $dbfile to generate $outfile) ...\n";
			$command = qq{grep -P '^#' $commentfile > $outfile};		#keep the comment lines in the output file
			print STDERR "NOTICE: Running <$command>\n";
			system ($command);
			$command = qq{grep -v -P '^#' $dbfile >> $outfile};		#keep the comment lines in the output file
			print STDERR "NOTICE: Running <$command>\n";
			system ($command);
		} else {
			print STDERR "NOTICE: Running the first step of indexing (copying $dbfile to $outfile) ...\n";
			system ("cp $dbfile $outfile") and die "Error: cannot run system command 'cp $dbfile $outfile'\n";
		}
	} else {
		if (defined $commentfile) {
			pod2usage ("Error in argument: -outfile must be different from input file when --commentfile is specified");
		}
	}
}


#step 2: generate the index file
print STDERR "NOTICE: Running the second step of indexing (generating $outfile.idx) ...\n";
$dbfile = $outfile;			#now the dbfile is the newdb generated in step 1
my $filesize = -s $dbfile;
my %region = ();
my ($offset, $lastregion, $firstoffset, $firstline) = (0, undef, 0, undef);

open (DB, $dbfile) or die "Error: cannot read from dbfile $dbfile: $!\n";
open (IDX, ">$dbfile.idx") or die "Error: cannot write to index file $dbfile.idx: $!\n";

print IDX "#BIN\t$bin\t$filesize\n";

while (<DB>) {
	my ($chr, $start);
	my $length = length ($_);
	s/[\r\n]+$//;
	
	if (m/^#/) {		#comment line is skipped
		$offset += $length;
		next;
	}
	
	if ($filetype eq 'A') {
		($chr, $start) = split (/\t/, $_);
	} elsif ($filetype eq 'B') {
		(undef, $chr, $start) = split (/\t/, $_);
		$start++;		#UCSC use zero-start
	} elsif ($filetype eq 'C') {
		(undef, undef, $chr, undef, $start) = split (/\t/, $_);
		$start++;
	}
	defined $start or die "Error: unable to find start site from input line <$_>\n";
	$start =~ m/^\d+$/ or die "Error: the start site ($start) is not a positive integer in input line <$_>\n";
	
	my $curbin = $start - ( $start % $bin );
	my $region = "$chr\t$curbin";
	$region{ $region }{ 'min' }   = $offset unless defined( $region{ $region }{ 'min' } );
	$region{ $region }{ 'max' }   = $offset + $length;
	
	$offset =~ m/000$/ and print STDERR sprintf("NOTICE: Indexing $dbfile: %d%%\r", int(100*$offset/$filesize));
	$offset += $length;
}

for my $k ( sort {$a cmp $b} keys %region ) {
	print IDX join("\t", $k, $region{ $k }{ 'min' }, $region{ $k }{ 'max' }), "\n";
}

print STDERR "\nDone!\n";
close(DB);
close(IDX);


=head1 SYNOPSIS

 index_annovar.pl [arguments] <db-file>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
            --filetype <A|B|C>		file type (default: A)
            --bin <int>			BIN size (default: 1000)
            --outfile <file>		prefix of output file name
            --skipsort			skip the pre-sorting procedure
            --commentfile <file>	provie comment lines (starting with #) from a comment file

 Function: generate index for ANNOVAR database files. type A start with chr, type B starts with bin
 
 Example: index_annovar.pl tempdb/hg19_cg69.txt -outfile humandb/hg19_cg69.txt
          index_annovar.pl tempdb/hg19_snp131.txt -outfile humandb/hg19_snp131.txt -filetype B
 
 Version: $LastChangedDate: 2012-10-23 23:32:05 -0700 (Tue, 23 Oct 2012) $
 
 WARNING: THIS PROGRAM IS STILL IN DEVELOPMENT PHASE AND MAY CONTAIN BUGS !

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=back

=head1 DESCRIPTION

This program will generate a new database file as well as an index file, given a 
user-specified database.

The file type A, B and C are explained below:

=over 8

A: first two tab-delimited fields are chr and start

B: first three tab-delimited fields are anything, chr and start

C: first five tab-delimited fields are anything, anything, chr, anything and start

=back

=cut