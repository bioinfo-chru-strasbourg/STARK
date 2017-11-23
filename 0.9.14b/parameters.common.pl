################
#
# PARAMETERS COMMON
#
# version: 1.2
# creation date: 30/04/2013
# author: Antony Le Bechec
#
# Process of common parameters
#
################

require "functions.common.pl";
use File::Basename;
use lib dirname (__FILE__);
$basename = dirname (__FILE__);

## Common options

@common_options=('output=s','projects_show','samples_show','family_show','disease_show', 'filters=s', 'format=s');
@common_options_default=('output'=>'','projects_show'=>0,'samples_show'=>0,'family_show'=>0,'disease_show'=>0, 'filters'=>'', 'format'=>'');


GetOptions (\%parameters,
	@options,
	@common_options
)  or pod2usage();

##
## MAIN options (Help, Man, Debug)
##

pod2usage(1) if $parameters{"help"};
pod2usage(-exitstatus => 0, -verbose => 2) if $parameters{"man"};
$DEBUG=1 if $parameters{"debug"};

##
## CONFIGURATION
##

## Check config file

if (-e $parameters{"config"}) {
    $config_file=$parameters{"config"};
} elsif (-e "$basename/".$parameters{"config"}) {
    $config_file="$basename/".$parameters{"config"};
} elsif (-e "$basename/config.ini") {
    $config_file="$basename/config.ini";
} else {
    print "No config.ini file...\n";
    pod2usage(1);
};#if


## Set parameters from config file

if (-e $config_file) {
	%config_ini=read_ini($config_file);
	#Folders
	$dir_projects=$config_ini{"folders"}{"projects"}."/";
	$annovar_databases=$config_ini{"folders"}{"annovar_databases"}."/";
	$ucsc_folder=$config_ini{"folders"}{"ucsc_folder"}."/";
	# annovar_folder in TRAKXS folder by default (if empty or undef, equal to TRAKXS folder)
	if (defined $config_ini{"folders"}{"annovar_folder"} && -d $config_ini{"folders"}{"annovar_folder"}) {
	    $annovar_folder=$config_ini{"folders"}{"annovar_folder"};
	} else {
	    $annovar_folder=$basename;    
	};#if
	if (trim($annovar_folder) ne "") {$annovar_folder.="/"};
	# vcftools_folder in $PATH by default ("" if not defined)
	if ((defined $config_ini{"folders"}{"vcftools_folder"} && -d $config_ini{"folders"}{"vcftools_folder"}) || trim($config_ini{"folders"}{"vcftools_folder"}) eq "") {
	    $vcftools_folder=$config_ini{"folders"}{"vcftools_folder"};    
	} else {
	    $vcftools_folder=$basename;    
	};#if
	if (trim($vcftools_folder) ne "") {$vcftools_folder.="/"};
	# R_folder in $PATH by default ("" if not defined)
	if ((defined $config_ini{"folders"}{"R_folder"} && -d $config_ini{"folders"}{"R_folder"}) || trim($config_ini{"folders"}{"R_folder"}) eq "") {
	    $R_folder=$config_ini{"folders"}{"R_folder"};    
	} else {
	    $R_folder=$basename;    
	};#if
	if (trim($R_folder) ne "") {$R_folder.="/"};
	
	#Database
	$host = $config_ini{"database"}{"host"};
	$driver = $config_ini{"database"}{"driver"};
	$database= $config_ini{"database"}{"database"};
	$user = $config_ini{"database"}{"user"};
	$pw = $config_ini{"database"}{"pw"};
	$port = $config_ini{"database"}{"port"};
	#Project
	$assembly=$config_ini{"project"}{"assembly"};
	$platform=$config_ini{"project"}{"platform"};
};#if

## Database connexion

#$dsn = "dbi:$driver:$database:$host:$port:mysql_local_infile=1:max_allowed_packet=512M";
#$DBIconnect = DBI->connect($dsn, $user, $pw);

## UCSC database connexion

$use_UCSC=1;
if ($use_UCSC) {
	#Database
	$host_ucsc = $config_ini{"UCSC"}{"host"};
	$driver_ucsc = $config_ini{"UCSC"}{"driver"};
	$database_ucsc = $config_ini{"UCSC"}{"database"};
	$user_ucsc = $config_ini{"UCSC"}{"user"};
	$pw_ucsc = $config_ini{"UCSC"}{"pw"};
	$port_ucsc = $config_ini{"UCSC"}{"port"};
	#$DBIconnect_UCSC = DBI->connect("dbi:$driver_ucsc:$database_ucsc:$host_ucsc:$port_ucsc", "$user_ucsc", "$pw_ucsc")
};#if

## Configuration of the annotation
if (trim($parameters{"config_annotation"}) eq "") {
    $parameters{"config_annotation"}="config.annotation.ini";
};#if
if (-e $parameters{"config_annotation"}) {
    $config_annotation_file=$parameters{"config_annotation"};
} elsif (-e "$basename/".$parameters{"config_annotation"}) {
    $config_annotation_file="$basename/".$parameters{"config_annotation"};
} elsif (-e "$basename/config.annotation.ini") {
    $config_annotation_file="$basename/config.annotation.ini";
} else {
    pod2usage(1);
};#if

# Read the config annotation file
%annotation_sources=read_ini($config_annotation_file);

# Construct the annotation type array
while (($source, $source_infos) = each(%annotation_sources)){
    if (!in_array(\@config_annotation_type_array,lc(trim($$source_infos{"annotation_type"})))) {
	    push(@config_annotation_type_array,lc(trim($$source_infos{"annotation_type"})));
    };#if
};#while


#Help parameter
#if (keys(%parameters)==0 || (exists $parameters{"help"} && $parameters{"help"}==1)) {
if ($parameters{"help"}) {
	help();
	exit 0;
};#if

#DEBUG
if (exists $parameters{"debug"}
    && $parameters{"debug"} ne "0"
    && trim(lc($parameters{"debug"})) ne "false") {
	$DEBUG=1;
	print "
	################
	## DEBUG MODE ##
	################
	
";
} else {
	$DEBUG=0;
};#if

## DEBUG
if ($parameters{"debug"}) {
    print "Options/Parameters:\n";
    while ((my $param_name, my $param_val) = each(%parameters)) {
            print "   $param_name => $param_val\n";
    };#while
    print "\n";
};#if


#VERBOSE
if ($parameters{"verbose"}) {
	$VERBOSE=1;
} else {
	$VERBOSE=0;
};#if

## output file
if (exists $parameters{"output"} && !(trim($parameters{"output"}) eq "") && !(trim($parameters{"output"}) eq "STDOUT")) {
	open (STDOUT,">".$parameters{"output"});
	$STDOUT_open=1;
};#if

#Format of the output, either 'std' or 'tab'
if (exists $parameters{"format"}) {
	$format=$parameters{"format"};
} else {
	$format="std";
};#if


## SHOW all projects
if ($parameters{"projects_show"}) {
	print projects_show($DBIconnect,$format);
	exit 0;
};

## SHOW samples for a project (or all if no --project option, but very long)
if ($parameters{"samples_show"}) {
	print samples_show($DBIconnect,$format);
	exit 0;
};

## SHOW diseases for a project (or all if no --project option, but very long)
if ($parameters{"disease_show"}) {
	print disease_show($DBIconnect,$format);
	exit 0;
};


## SHOW families for a project (or all if no --project option, but very long)
if ($parameters{"family_show"}) {
	print family_show($DBIconnect,$format);
	exit 0;
};

## SHOW variants
if ($parameters{"variants_show"}) {
	print variants_show($DBIconnect,$format);
	exit 0;
};

## SHOW variants with genotype
if ($parameters{"variants_genotype_show"}) {
	print variants_show($DBIconnect,$format,1);
	exit 0;
};


# Show available sources
if (exists $parameters{"source_show"} && $parameters{"source_show"}) {
	#print "Source\t\tAnnovar code\t\tAnnotation type\n";
	printf("%-30s%-30s%-30s\n", "SOURCE", "ANNOVAR CODE", "ANNOTATION TYPE");
	while (($source, $source_infos) = each(%annotation_sources)){
		#printf "$source", 50;
		#printf "\t".$$source_infos{"annovar_code"}."\t\t".$$source_infos{"annotation_type"}."\n";
		if (exists $$source_infos{"available"} && ($$source_infos{"available"} eq "true" || $$source_infos{"available"} eq "1")) {
			printf("%-30s%-30s%-30s\n", $source, $$source_infos{"annovar_code"}, $$source_infos{"annotation_type"})
		};#if
	};#while
	exit 0;
};#if

return 1;

## END

our $help_parameters_common="
Other parameters:
--help   		print help message
--output=<FILE>		put the result on a file 
--projects_show		show projects stored in the database
--samples_show		show samples for a project (or all projects if no --project option, but can take a while)
--family_show		show family for a project (or all projects if no --project option, but can take a while)
--disease_show		show disease for a project (or all projects if no --project option, but can take a while)
--format=<enum>		format of the output for 'show' parameters, either 'std' (text) or 'tab' (tab-delimiter format)
--filters=<string>	define filters for 'show' parameters, with the format 'field1:value1,field2:value2,table.field3:value3...'
			(example 'project filter': --filters=platform:GTF,assembly:19 for all projects from the BGI platform with the assembly hg19)
			(example 'sample filter': --filters=family:Family1,status:patient for all patients of the Family1)
			(example 'sample ambigious field ref': --filters=project.ref:Project1,status:patient for all patients of the Project1. )
";
