################
#
# FUNCTIONS COMMON
#
# version: 1.2.1
# creation date: 18/01/2013
# author: Antony Le BŽchec
#
################

#@common_functions=

#Function projects_show
#List projects stored in the database
sub projects_show {
#$_[0]: The database connexion
#$_[1]: Format of the output, either 'std' (text) or 'tab' (tab-delimiter format)

	#database connexion
	my $DBIconnect=$_[0];
	if ($DBIconnect eq "") {
		return 0;
	};#if
	#output format
	my $format=$_[1];
	if ($format eq "") {
		$format="std";
	};#if
	my $result="";
	
	#Filters
	my %filters;
	if ($parameters{"filters"} ne "") { #global parameters
		@filters_list=split(",",trim($parameters{"filters"}));
		foreach $filter (@filters_list) {
			@filter_values=split(":",trim($filter));
			#print $filter_values[0]." = ".$filter_values[1]."\n";
			$filters{trim($filter_values[0])}=trim($filter_values[1]);
		};#foreach
	};#if
	my $query_parameters_filter="";
	
	#fields of the tables for filters and for length of the fields
	my @tables = ("project", "sample");
	foreach $table (@tables) {
		$query="DESCRIBE $table";
		$query_handle = $DBIconnect->prepare($query);
		$query_handle->execute();
		$query_handle->bind_columns(undef, \$Field, \$Type, \$Null, \$Key, \$Default, \$Extra);
		while($query_handle->fetch()) {
			$field_name=$Field;
			$field_name_table="$table.".$Field;
			$Type =~ s/\(/,/;
			@Type_split=split(",",$Type);
			#print "$field_name $field_name_table Type ".$Type_split[0]."\n";
			$field_type=$Type_split[0];
			if (exists $filters{$field_name}) {
				if ($field_type eq "int") {
					$query_parameters_filter.=" AND $field_name_table=".$filters{$field_name}." ";
				} else {
					$query_parameters_filter.=" AND $field_name_table LIKE '".$filters{$field_name}."' ";
				};#if
			};#if
			if (exists $filters{$field_name_table}) {
				if ($field_type eq "int") {
					$query_parameters_filter.=" AND $field_name_table=".$filters{$field_name_table}." ";
				} else {
					$query_parameters_filter.=" AND $field_name_table LIKE '".$filters{$field_name_table}."' ";
				};#if
			};#if
		};#while
	};#foreach
	
	#Max length of fields
	my %size_fields;
	my %max_length;
	$query="SELECT * 
		FROM  `project`, `sample`
		WHERE project.project_id=sample.project_id
		$query_parameters_filter
		PROCEDURE ANALYSE ( )";
	#print $query."\n";
	$query_handle = $DBIconnect->prepare($query);
	$query_handle->execute();
	$query_handle->bind_columns(undef, \$Field_name, \$Min_value, \$Max_value, \$Min_length, \$Max_length, \$Empties_or_zeros, \$Nulls, \$Avg_value_or_avg_length, \$Std, \$Optimal_fieldtype);
	while($query_handle->fetch()) {
		#find the field
		#@Field_name_split=split(/\./,$Field_name);
		#$field_name=$Field_name_split[2];
		#$size_fields{$field_name}=int($Max_length)+1;
		#$max_length{$field_name}=int($Max_length);
		
		#find the field
		@Field_name_split=split(/\./,$Field_name);
		$field_name=$Field_name_split[2];
		$field_name_table=$Field_name_split[1].".".$Field_name_split[2];
		$size_fields{$field_name_table}=int($Max_length)+1;
		$max_length{$field_name_table}=int($Max_length);
		#print "$field_name_table  ".$max_length{$field_name_table}."\n";
	};#while
	
	@fields_sql_array=("project.project_id", "project.ref", "project.assembly", "project.platform", "project.pi", "count(sample.sample_id) as nb_samples");
	$query = "
		SELECT ".join(", ",@fields_sql_array)." #project.project_id, project.ref, project.assembly, project.platform, project.pi, count(sample.sample_id) as nb_samples
		FROM project, sample
		WHERE project.project_id=sample.project_id
		$query_parameters_filter
		GROUP BY project.project_id
	;";
	$query_handle = $DBIconnect->prepare($query);
	$query_handle->execute();
	$query_handle->bind_columns(undef, \$project_id, \$project_ref, \$project_assembly, \$project_platform, \$project_pi, \$nb_samples);
	#print "List of projects in the DB:\n";
	@fields_head_array=("ID", "REFERENCE", "ASSEMBLY", "PLATFORM", "PI", "#SAMPLES");
	#$std_format="%-4s%-30s%-10s%-10s%-30s%-10s\n";
	$std_format="";
	$field_num=0;
	foreach $field_sql (@fields_sql_array) {
		my @array = (int($max_length{$field_sql})+0, length($fields_head_array[$field_num])+0);
		@array = reverse sort { $a <=> $b } @array;
		$std_format.="%-".($array[0]+2)."s";
		$field_num++;
	};#foreach
	$std_format.="\n";
	if ($format eq "tab") {
		$result=join("\t",@fields_head_array)."\n";
	} else {
		$result=sprintf($std_format, @fields_head_array);
	};#if
	#print $std_format."\n";
	while($query_handle->fetch()) {
		#if (trim($project_assembly) eq "") { $project_assembly=$assembly; };#if
		#if (grep(/^[1..9]*$/,$project_assembly)) { $project_assembly="hg".$project_assembly; }#if
		#if (trim($project_plaform) eq "") { $project_plaform=$platform; };#if
		@fields_array=($project_id, $project_ref, $project_assembly, $project_platform, $project_pi, $nb_samples);
		if ($format eq "tab") {
			$line_print.=join("\t",@fields_array); $line_print.="\n";
		} else {
			$line_print.=sprintf($std_format, @fields_array);
		};#if
	};#while
	$result.=$line_print;
	return $result;
	
}

#Function samples_show
#List samples stored in the database, related to a project or not
sub samples_show {
#$_[0]: The database connexion
#$_[1]: Format of the output, either 'std' (text) or 'tab' (tab-delimiter format)
#$_[2]: project

	#database connexion
	my $DBIconnect=$_[0];
	if ($DBIconnect eq "") {
		return 0;
	};#if
	#output format
	my $format=$_[1];
	if ($format eq "") {
		$format="std";
	};#if
	my $parameters=$_[2];
	if ($parameters{"project"} ne "") {
		$query_projectref_filter="  AND project.ref='".$parameters{"project"}."' ";
	};#if
	
	#Filters
	my %filters;
	if ($parameters{"filters"} ne "") {
		@filters_list=split(",",trim($parameters{"filters"}));
		foreach $filter (@filters_list) {
			@filter_values=split(":",trim($filter));
			#print $filter_values[0]." = ".$filter_values[1]."\n";
			$filters{trim($filter_values[0])}=trim($filter_values[1]);
		};#foreach
	};#if
	my $query_parameters_filter="";
	
	#fields of the tables for filters and for length of the fields
	my @tables = ("sample", "project");
	foreach $table (@tables) {
		$query="DESCRIBE $table";
		$query_handle = $DBIconnect->prepare($query);
		$query_handle->execute();
		$query_handle->bind_columns(undef, \$Field, \$Type, \$Null, \$Key, \$Default, \$Extra);
		while($query_handle->fetch()) {
			$field_name=$Field;
			$field_name_table="$table.".$Field;
			$Type =~ s/\(/,/;
			@Type_split=split(",",$Type);
			#print "$field_name $field_name_table Type ".$Type_split[0]."\n";
			$field_type=$Type_split[0];
			if (exists $filters{$field_name}) {
				if ($field_type eq "int") {
					$query_parameters_filter.=" AND $field_name_table=".$filters{$field_name}." ";
				} else {
					$query_parameters_filter.=" AND $field_name_table LIKE '".$filters{$field_name}."' ";
				};#if
			};#if
			if (exists $filters{$field_name_table}) {
				if ($field_type eq "int") {
					$query_parameters_filter.=" AND $field_name_table=".$filters{$field_name_table}." ";
				} else {
					$query_parameters_filter.=" AND $field_name_table LIKE '".$filters{$field_name_table}."' ";
				};#if
			};#if
		};#while
	};#foreach
	
	#Max length of the fields
	my %size_fields;
	my %max_length;
	$query="SELECT * 
		FROM `sample`
		INNER JOIN project ON (sample.project_id=project.project_id)
		WHERE 1=1
		$query_projectref_filter
		$query_parameters_filter
		PROCEDURE ANALYSE ( )";
	$query_handle = $DBIconnect->prepare($query);
	$query_handle->execute();
	$query_handle->bind_columns(undef, \$Field_name, \$Min_value, \$Max_value, \$Min_length, \$Max_length, \$Empties_or_zeros, \$Nulls, \$Avg_value_or_avg_length, \$Std, \$Optimal_fieldtype);
	while($query_handle->fetch()) {
		#find the field
		@Field_name_split=split(/\./,$Field_name);
		$field_name=$Field_name_split[2];
		$field_name_table=$Field_name_split[1].".".$Field_name_split[2];
		$size_fields{$field_name_table}=int($Max_length)+1;
		$max_length{$field_name_table}=int($Max_length);
	};#while
	
	#query
	@fields_sql_array=("sample.sample_id", "project.ref", "sample.run", "sample.ref", "sample.number", "sample.status", "sample.gender", "sample.disease", "sample.diseasesubtype", "sample.inheritance", "sample.family", "sample.motherid", "sample.fatherid", "sample.share", "sample.flag", "count(variant_id) AS nb_variants");
	$query = "
		SELECT ".join(", ",@fields_sql_array)." #sample.sample_id, project.ref, sample.run, sample.ref, sample.number, sample.status, sample.gender, sample.disease, sample.diseasesubtype, sample.inheritance, sample.family, sample.motherid, sample.fatherid, sample.share, sample.flag, count(variant_id) AS nb_variants
		FROM sample
		INNER JOIN variant_sample ON (sample.sample_id=variant_sample.sample_id)
		INNER JOIN project ON (sample.project_id=project.project_id)
		WHERE 1=1
		$query_projectref_filter
		$query_parameters_filter
		GROUP BY sample.sample_id
		ORDER BY sample.project_id, sample.disease, sample.family, sample.status ASC, sample.ref
	;";
	print $query."\n" if ($DEBUG);
	$query_handle = $DBIconnect->prepare($query);
	$query_handle->execute();
	$query_handle->bind_columns(undef, \$sample_id, \$project_ref, \$sample_run, \$sample_ref, \$sample_number, \$sample_status, \$sample_gender, \$sample_disease, \$sample_diseasesubtype, \$sample_inheritance, \$sample_family, \$sample_motherid, \$sample_fatherid, \$sample_share, \$sample_flag, \$nb_variants);
	#print "List of projects in the DB:\n";
	@fields_head_array=("ID", "PROJECT", "RUN", "REFERENCE", "NUMBER", "STATUS", "GENDER", "DISEASE", "SUBDISEASE", "INHERITANCE", "FAMILY", "MOTHER", "FATHER", "SHARE", "FLAG", "#VARIANTS");
	#$std_format="%-4s%-15s%-15s%-10s%-20s%-20s%-20s%-20s%-15s%-15s%-10s\n";
	$std_format="";
	$field_num=0;
	foreach $field_sql (@fields_sql_array) {
		my @array = (int($max_length{$field_sql})+0, length($fields_head_array[$field_num])+0);
		@array = reverse sort { $a <=> $b } @array;
		$std_format.="%-".($array[0]+2)."s";
		$field_num++;
	};#foreach
	$std_format.="\n";
	#$std_format="%-".$size_fields{"sample_id"}."s%-15s%-15s%-10s%-20s%-20s%-20s%-20s%-15s%-15s%-10s\n";
	#print $std_format."\n";
	if ($format eq "tab") {
		$result=join("\t",@fields_head_array)."\n";
	} else {
		$result=sprintf($std_format, @fields_head_array);
	};#if
	while($query_handle->fetch()) {
		#if ($format eq "std" && length($sample_disease) >= 20) { $sample_disease=substr($sample_disease,0,16)."..." }
		#if ($format eq "std" && length($sample_diseasesubtype) >= 20) { $sample_diseasesubtype=substr($sample_diseasesubtype,0,16)."..." }
		@fields_array=($sample_id, $project_ref, $sample_run, $sample_ref, $sample_number, $sample_status, $sample_gender, $sample_disease, $sample_diseasesubtype, $sample_inheritance, $sample_family, $sample_motherid, $sample_fatherid, $sample_share, $sample_flag, $nb_variants);
		if ($format eq "tab") {
			$line_print.=join("\t",@fields_array)."\n";
		} else {
			$line_print.=sprintf($std_format, @fields_array);
		};#if
	};#while
	$result.=$line_print;
	return $result;
	
}

#Function family_show
#List families stored in the database, related to a project
sub family_show {
#database connexion
	my $DBIconnect=$_[0];
	if ($DBIconnect eq "") {
		return 0;
	};#if
	#output format
	my $format=$_[1];
	if ($format eq "") {
		$format="std";
	};#if
	my $parameters=$_[2];
	if ($parameters{"project"} ne "") {
		$query_projectref_filter="  AND project.ref='".$parameters{"project"}."' ";
	};#if
	my $result="";
	
	#Filters
	my %filters;
	if ($parameters{"filters"} ne "") { #global parameters
		@filters_list=split(",",trim($parameters{"filters"}));
		foreach $filter (@filters_list) {
			@filter_values=split(":",trim($filter));
			#print $filter_values[0]." = ".$filter_values[1]."\n";
			$filters{trim($filter_values[0])}=trim($filter_values[1]);
		};#foreach
	};#if
	my $query_parameters_filter="";
	
	#fields of the tables for filters and for length of the fields
	my @tables = ("project", "sample");
	foreach $table (@tables) {
		$query="DESCRIBE $table";
		$query_handle = $DBIconnect->prepare($query);
		$query_handle->execute();
		$query_handle->bind_columns(undef, \$Field, \$Type, \$Null, \$Key, \$Default, \$Extra);
		while($query_handle->fetch()) {
			$field_name=$Field;
			$field_name_table="$table.".$Field;
			$Type =~ s/\(/,/;
			@Type_split=split(",",$Type);
			#print "$field_name $field_name_table Type ".$Type_split[0]."\n";
			$field_type=$Type_split[0];
			if (exists $filters{$field_name}) {
				if ($field_type eq "int") {
					$query_parameters_filter.=" AND $field_name_table=".$filters{$field_name}." ";
				} else {
					$query_parameters_filter.=" AND $field_name_table LIKE '".$filters{$field_name}."' ";
				};#if
			};#if
			if (exists $filters{$field_name_table}) {
				if ($field_type eq "int") {
					$query_parameters_filter.=" AND $field_name_table=".$filters{$field_name_table}." ";
				} else {
					$query_parameters_filter.=" AND $field_name_table LIKE '".$filters{$field_name_table}."' ";
				};#if
			};#if
		};#while
	};#foreach
	
	#Max length of fields
	my %size_fields;
	my %max_length;
	$query="SELECT * 
		FROM  `project`, `sample`
		WHERE project.project_id=sample.project_id
		$query_projectref_filter
		$query_parameters_filter
		PROCEDURE ANALYSE ( )";
	#print $query."\n";
	$query_handle = $DBIconnect->prepare($query);
	$query_handle->execute();
	$query_handle->bind_columns(undef, \$Field_name, \$Min_value, \$Max_value, \$Min_length, \$Max_length, \$Empties_or_zeros, \$Nulls, \$Avg_value_or_avg_length, \$Std, \$Optimal_fieldtype);
	while($query_handle->fetch()) {
		#find the field
		@Field_name_split=split(/\./,$Field_name);
		$field_name=$Field_name_split[2];
		$field_name_table=$Field_name_split[1].".".$Field_name_split[2];
		$size_fields{$field_name_table}=int($Max_length)+1;
		$max_length{$field_name_table}=int($Max_length);
	};#while
	
	@fields_sql_array=("project.project_id", "project.ref", "sample.family", "count(sample.sample_id) as nb_samples");
	$query = "
		SELECT ".join(", ",@fields_sql_array)." #project.project_id, project.ref, sample.family, count(sample.sample_id) as nb_samples
		FROM project, sample
		WHERE project.project_id=sample.project_id
		$query_projectref_filter
		$query_parameters_filter
		GROUP BY project.project_id, sample.family
	;";
	$query_handle = $DBIconnect->prepare($query);
	$query_handle->execute();
	$query_handle->bind_columns(undef, \$project_id, \$project_ref, \$sample_family, \$nb_samples);
	#print "List of projects in the DB:\n";
	@fields_head_array=("ID", "PROJECT", "FAMILY", "#SAMPLES");
	#$std_format="%-4s%-30s%-10s%-10s%-30s%-10s\n";
	$std_format="";
	$field_num=0;
	foreach $field_sql (@fields_sql_array) {
		my @array = (int($max_length{$field_sql})+0, length($fields_head_array[$field_num])+0);
		@array = reverse sort { $a <=> $b } @array;
		$std_format.="%-".($array[0]+2)."s";
		$field_num++;
	};#foreach
	$std_format.="\n";
	if ($format eq "tab") {
		$result=join("\t",@fields_head_array)."\n";
	} else {
		$result=sprintf($std_format, @fields_head_array);
	};#if
	#print $std_format."\n";
	while($query_handle->fetch()) {
		#if (trim($project_assembly) eq "") { $project_assembly=$assembly; };#if
		#if (grep(/^[1..9]*$/,$project_assembly)) { $project_assembly="hg".$project_assembly; }#if
		#if (trim($project_plaform) eq "") { $project_plaform=$platform; };#if
		@fields_array=($project_id, $project_ref, $sample_family, $nb_samples);
		if ($format eq "tab") {
			$line_print.=join("\t",@fields_array); $line_print.="\n";
		} else {
			$line_print.=sprintf($std_format, @fields_array);
		};#if
	};#while
	$result.=$line_print;
	return $result;
	
}


#Function disease_show
#List families stored in the database, related to a project
sub disease_show {
#database connexion
	my $DBIconnect=$_[0];
	if ($DBIconnect eq "") {
		return 0;
	};#if
	#output format
	my $format=$_[1];
	if ($format eq "") {
		$format="std";
	};#if
	my $parameters=$_[2];
	if ($parameters{"project"} ne "") {
		$query_projectref_filter="  AND project.ref='".$parameters{"project"}."' ";
	};#if
	my $result="";
	
	#Filters
	my %filters;
	if ($parameters{"filters"} ne "") { #global parameters
		@filters_list=split(",",trim($parameters{"filters"}));
		foreach $filter (@filters_list) {
			@filter_values=split(":",trim($filter));
			#print $filter_values[0]." = ".$filter_values[1]."\n";
			$filters{trim($filter_values[0])}=trim($filter_values[1]);
		};#foreach
	};#if
	my $query_parameters_filter="";
	
	#fields of the tables for filters and for length of the fields
	my @tables = ("project", "sample");
	foreach $table (@tables) {
		$query="DESCRIBE $table";
		$query_handle = $DBIconnect->prepare($query);
		$query_handle->execute();
		$query_handle->bind_columns(undef, \$Field, \$Type, \$Null, \$Key, \$Default, \$Extra);
		while($query_handle->fetch()) {
			$field_name=$Field;
			$field_name_table="$table.".$Field;
			$Type =~ s/\(/,/;
			@Type_split=split(",",$Type);
			#print "$field_name $field_name_table Type ".$Type_split[0]."\n";
			$field_type=$Type_split[0];
			if (exists $filters{$field_name}) {
				if ($field_type eq "int") {
					$query_parameters_filter.=" AND $field_name_table=".$filters{$field_name}." ";
				} else {
					$query_parameters_filter.=" AND $field_name_table LIKE '".$filters{$field_name}."' ";
				};#if
			};#if
			if (exists $filters{$field_name_table}) {
				if ($field_type eq "int") {
					$query_parameters_filter.=" AND $field_name_table=".$filters{$field_name_table}." ";
				} else {
					$query_parameters_filter.=" AND $field_name_table LIKE '".$filters{$field_name_table}."' ";
				};#if
			};#if
		};#while
	};#foreach
	
	#Max length of fields
	my %size_fields;
	my %max_length;
	$query="SELECT * 
		FROM  `project`, `sample`
		WHERE project.project_id=sample.project_id
		$query_projectref_filter
		$query_parameters_filter
		PROCEDURE ANALYSE ( )";
	#print $query."\n";
	$query_handle = $DBIconnect->prepare($query);
	$query_handle->execute();
	$query_handle->bind_columns(undef, \$Field_name, \$Min_value, \$Max_value, \$Min_length, \$Max_length, \$Empties_or_zeros, \$Nulls, \$Avg_value_or_avg_length, \$Std, \$Optimal_fieldtype);
	while($query_handle->fetch()) {
		#find the field
		@Field_name_split=split(/\./,$Field_name);
		$field_name=$Field_name_split[2];
		$field_name_table=$Field_name_split[1].".".$Field_name_split[2];
		$size_fields{$field_name_table}=int($Max_length)+1;
		$max_length{$field_name_table}=int($Max_length);
	};#while
	
	@fields_sql_array=("project.project_id", "project.ref", "sample.disease", "count(sample.sample_id) as nb_samples");
	$query = "
		SELECT ".join(", ",@fields_sql_array)." #project.project_id, project.ref, project.assembly, project.platform, project.pi, count(sample.sample_id) as nb_samples
		FROM project, sample
		WHERE project.project_id=sample.project_id
		$query_projectref_filter
		$query_parameters_filter
		GROUP BY project.project_id, sample.disease
	;";
	$query_handle = $DBIconnect->prepare($query);
	$query_handle->execute();
	$query_handle->bind_columns(undef, \$project_id, \$project_ref, \$sample_disease, \$nb_samples);
	#print "List of projects in the DB:\n";
	@fields_head_array=("ID", "PROJECT", "DISEASE", "#SAMPLES");
	#$std_format="%-4s%-30s%-10s%-10s%-30s%-10s\n";
	$std_format="";
	$field_num=0;
	foreach $field_sql (@fields_sql_array) {
		my @array = (int($max_length{$field_sql})+0, length($fields_head_array[$field_num])+0);
		@array = reverse sort { $a <=> $b } @array;
		$std_format.="%-".($array[0]+2)."s";
		$field_num++;
	};#foreach
	$std_format.="\n";
	if ($format eq "tab") {
		$result=join("\t",@fields_head_array)."\n";
	} else {
		$result=sprintf($std_format, @fields_head_array);
	};#if
	#print $std_format."\n";
	while($query_handle->fetch()) {
		#if (trim($project_assembly) eq "") { $project_assembly=$assembly; };#if
		#if (grep(/^[1..9]*$/,$project_assembly)) { $project_assembly="hg".$project_assembly; }#if
		#if (trim($project_plaform) eq "") { $project_plaform=$platform; };#if
		@fields_array=($project_id, $project_ref, $sample_disease, $nb_samples);
		if ($format eq "tab") {
			$line_print.=join("\t",@fields_array); $line_print.="\n";
		} else {
			$line_print.=sprintf($std_format, @fields_array);
		};#if
	};#while
	$result.=$line_print;
	return $result;
	
}


#Function variants_show
#List variants stored in the database, related to a project or not
sub variants_show {
#$_[0]: The database connexion
#$_[1]: Format of the output, either 'std' (text) or 'tab' (tab-delimiter format)
#$_[2]: add genotype

	#database connexion
	my $DBIconnect=$_[0];
	if ($DBIconnect eq "") {
		return 0;
	};#if
	#output format
	my $format=$_[1];
	if ($format eq "") {
		$format="std";
	};#if
	#with genotype
	my $with_genotype=$_[2];
	if ($with_genotype eq "false") {
		$with_genotype=0;
	};#if
	
	#Filters
	my %filters;
	if ($parameters{"filters"} ne "") {
		@filters_list=split(",",trim($parameters{"filters"}));
		foreach $filter (@filters_list) {
			@filter_values=split(":",trim($filter));
			#print $filter_values[0]." = ".$filter_values[1]."\n";
			$filters{trim($filter_values[0])}=trim($filter_values[1]);
		};#foreach
	};#if
	my $query_parameters_filter="";
	my $query_parameters_filter_special="";
	
	#fields of the tables for filters and for length of the fields
	my @tables = ("variant","variant_sample");
	foreach $table (@tables) {
		$query="DESCRIBE $table";
		$query_handle = $DBIconnect->prepare($query);
		$query_handle->execute();
		$query_handle->bind_columns(undef, \$Field, \$Type, \$Null, \$Key, \$Default, \$Extra);
		while($query_handle->fetch()) {
			$field_name=$Field;
			$field_name_table="$table.".$Field;
			$Type =~ s/\(/,/;
			@Type_split=split(",",$Type);
			#print "$field_name $field_name_table Type ".$Type_split[0]."\n";
			$field_type=$Type_split[0];
			if (exists $filters{$field_name}) {
				if ($field_type eq "int" || $field_type eq "smallint" || $field_type eq "tinyint" || $field_type eq "mediumint" || $field_type eq "float") {
					if (nifty_number($filters{$field_name})) {
						$query_parameters_filter.=" AND $field_name_table=".$filters{$field_name}." ";
					} else {
						$query_parameters_filter.=" AND $field_name_table ".$filters{$field_name}." ";
					};#if
				} else {
					$query_parameters_filter.=" AND $field_name_table LIKE '".$filters{$field_name}."' ";
				};#if
			};#if
			if (exists $filters{$field_name_table}) {
				if ($field_type eq "int" || $field_type eq "smallint" || $field_type eq "tinyint" || $field_type eq "mediumint" || $field_type eq "float") {
					$query_parameters_filter.=" AND $field_name_table=".$filters{$field_name_table}." ";
				} else {
					$query_parameters_filter.=" AND $field_name_table LIKE '".$filters{$field_name_table}."' ";
				};#if
			};#if
		};#while
	};#foreach
	#Special filter
	if (exists $filters{"nb_samples"}) {
		$query_parameters_filter_special.=" HAVING COUNT(variant_sample.sample_id)".$filters{"nb_samples"}." ";
	};#if
	
	#Max length of the fields
	my %size_fields;
	my %max_length;
	$query="SELECT * 
		FROM `variant`
		INNER JOIN variant_sample ON (variant.variant_id=variant_sample.variant_id)
		WHERE variant.variant_id<=100000
		$query_parameters_filter
		PROCEDURE ANALYSE ( )";
	#print "$query\n";
	$query_handle = $DBIconnect->prepare($query);
	if (!$query_handle->execute()) {	
		print "ERROR: error in query\n";
		exit 0;
	};
	$query_handle->bind_columns(undef, \$Field_name, \$Min_value, \$Max_value, \$Min_length, \$Max_length, \$Empties_or_zeros, \$Nulls, \$Avg_value_or_avg_length, \$Std, \$Optimal_fieldtype);
	while($query_handle->fetch()) {
		#find the field
		@Field_name_split=split(/\./,$Field_name);
		$field_name=$Field_name_split[2];
		$field_name_table=$Field_name_split[1].".".$Field_name_split[2];
		$size_fields{$field_name_table}=int($Max_length)+1;
		$max_length{$field_name_table}=int($Max_length);
	};#while
	#patch
	$max_length{"variant.variant_id"}=10;
	
	#query
	$max_line=100;
	@fields_sql_array=("variant.variant_id", "variant.assembly", "variant.chrom", "variant.pos", "variant.ref", "variant.alt", "variant.type", "variant.location", "variant.outcome", "variant.gene", "variant.snpid", "count(variant_sample.sample_id) AS nb_samples");
	$query = "
		SELECT ".join(", ",@fields_sql_array)." 
		FROM `variant_sample`
		INNER JOIN sample ON (sample.sample_id=`variant_sample`.sample_id)
		INNER JOIN variant ON (variant.variant_id=`variant_sample`.variant_id)
		WHERE 1=1
		$query_projectref_filter
		$query_parameters_filter
		GROUP BY variant.variant_id
		$query_parameters_filter_special
		ORDER BY variant.variant_id
		LIMIT $max_line
	;";
	print $query."\n" if ($DEBUG);
	$query_handle = $DBIconnect->prepare($query);
	$query_handle->execute();
	$query_handle->bind_columns(undef, \$variant_id, \$my_assembly, \$chrom, \$pos, \$ref, \$alt, \$type, \$location, \$outcome, \$gene, \$snpid, \$nb_samples);
	#print "List of projects in the DB:\n";
	@fields_head_array=("ID", "ASSEMBLY", "CHROM", "POS", "REF", "ALT", "TYPE", "LOCATION", "OUTCOME", "GENE", "SNPID", "NB SAMPLES");
	#$std_format="%-4s%-15s%-15s%-10s%-20s%-20s%-20s%-20s%-15s%-15s%-10s\n";
	$std_format="";
	$field_num=0;
	foreach $field_sql (@fields_sql_array) {
		my @array = (int($max_length{$field_sql})+0, length($fields_head_array[$field_num])+0);
		@array = reverse sort { $a <=> $b } @array;
		$std_format.="%-".($array[0]+2)."s";
		$field_num++;
	};#foreach
	if ($with_genotype) {
		$std_format.="%-10s";
		push(@fields_head_array,"GENOTYPES");
	};#if
	$std_format.="\n";
	#$std_format="%-".$size_fields{"sample_id"}."s%-15s%-15s%-10s%-20s%-20s%-20s%-20s%-15s%-15s%-10s\n";
	#print $std_format."\n";
	if ($format eq "tab") {
		$result=join("\t",@fields_head_array)."\n";
	} else {
		$result=sprintf($std_format, @fields_head_array);
	};#if
	$nb_line=0;
	while($query_handle->fetch()) {
		$nb_line++;
		#if ($format eq "std" && length($sample_disease) >= 20) { $sample_disease=substr($sample_disease,0,16)."..." }
		#if ($format eq "std" && length($sample_diseasesubtype) >= 20) { $sample_diseasesubtype=substr($sample_diseasesubtype,0,16)."..." }
		@fields_array=($variant_id, $my_assembly, $chrom, $pos, $ref, $alt, $type, $location, $outcome, $gene, $snpid, $nb_samples);
		
		if ($with_genotype) {
			$genotypes="";
			$query_genotype = "
			SELECT sample.ref, genotype 
			FROM `variant_sample`
			INNER JOIN sample ON (variant_sample.sample_id=sample.sample_id)
			INNER JOIN variant ON (variant.variant_id=`variant_sample`.variant_id)
			WHERE variant_sample.variant_id=$variant_id
			$query_projectref_filter
			$query_parameters_filter
			;";
			#print $query_genotype."\n" if ($DEBUG || 1);
			$query_handle_genotype = $DBIconnect->prepare($query_genotype);
			$query_handle_genotype->execute();
			$query_handle_genotype->bind_columns(undef, \$sample_ref, \$genotype);
			while($query_handle_genotype->fetch()) {
				$genotypes.="$sample_ref:$genotype, ";
			}#while
			#print $genotypes."\n";
			push(@fields_array,$genotypes);
		};#if
		#print "$std_format, @fields_array\n";
		
		if ($format eq "tab") {
			$line_print.=join("\t",@fields_array)."\n";
		} else {
			$line_print.=sprintf($std_format, @fields_array);
		};#if
	};#while
	if ($nb_line == $max_line) {
		$line_print.="...\n";
	};#if
	$result.=$line_print;
	return $result;
	
}




### read_infos ###
#Read an infos tab-delimiter file
sub read_infos {
#$_[0]: info file
#$_[1]: boolean. To test the unicity of the sampleRef
	#parameters
	$samples_infos_file=$_[0];
	$only_integrate=$_[1]; #To test the unicity of the sampleRef
	my %infos;
	#init
	@samples_header=[];
	my %samples_header_inv;
	if (-e $samples_infos_file) { #if file exist		
		#open the file
		open(FILE_SAMPLES_INFOS, $samples_infos_file) || die "Problem to open the file: $!";
		$line=0;
		#read the file
		while(<FILE_SAMPLES_INFOS>) {
			$line++;
			if ($line==1) { #header at the first line
				#split columns
				@samples_header=split("\t",$_);
				$l=0; # first column as index by default
				for $head_value (@samples_header) {
					#clean the value
					$head_value =~ s/_//; $head_value =~ s/ //; $head_value =~ s/-//; $head_value =~ s/#//;
					#find the sampleid
					if (lc($head_value) eq "sampleid" || lc($head_value) eq "sampleref"
					    || lc($head_value) eq "id" || lc($head_value) eq "ref") {
						$sampleid_line=$l;
					};#if
					if (lc($head_value) eq "integrate") {
						$integrate_line=$l;
					};#if	
					if (lc($head_value) eq "runref") {
						$runref_line=$l;
					};#if	
					#put name of the header in the header table
					$samples_header[$l]=trim(lc($head_value));			
					$samples_header_inv{trim(lc($head_value))}=$l;			
					$l++;
				};#if
				
			} else { #Each line correspond to a sample
				#split columns
				@sample_info=split("\t",$_);
				#test if exlusions
				if (!$only_integrate || $sample_info[$integrate_line] eq "1") {
					#prepare the hash table of samples
					my $sampleid=trim($sample_info[$sampleid_line]);
					my $runref=trim($sample_info[$runref_line]);
					if (exists $infos{$runref}{$sampleid}) { # if index already in the hash, create an error
						$infos{$runref}{"error"}{"multipleIndex"}=1;
						#print "ERROR $line $sampleid\n";
					};#if
					$l=0;
					for $sample_value (@sample_info) {
						#Put values on the hash table of samples
						$infos{$runref}{$sampleid}{$samples_header[$l]} = trim($sample_value);
						$l++;
					};#for
				};#if
			}#if			
		};
		close(FILE_SAMPLES_INFOS);
		return %infos;
	} else {
		print "File '$samples_infos_file' DOES NOT exist\n";
		exit 1;
	};#if
}


### write_infos ###
#Write an infos tab-delimiter file
sub write_infos {
#$_[0]: info file
#$_[1]: hash

	#parameters
	my $samples_infos_file=$_[0];
	my $samples_infos_hash=$_[1];
	
	# Create a saveid
	use Time::localtime;
	my $t = localtime;
	my $save_id=sprintf("%04d%02d%02d%02d%02d%02d",$t->year + 1900, $t->mon + 1, $t->mday, $t->hour, $t->min, $t->sec);

	print "save_id:$save_id\n" if $DEBUG;
	
	if (-e $samples_infos_file) {
		$cmd_mv = "mv $samples_infos_file $samples_infos_file.$save_id";
		#print $cmd_mv."\n";
		$result = `$cmd_mv 2>&1`;
	};#if
	
	open(SAMPLES_INFOS, ">$samples_infos_file") || die "Problem: $!";
	print SAMPLES_INFOS "sampleRef\tProjectRef\tRunRef\tsampleNumber\tGender\tFamily\tmotherID\tfatherID\trelationship\tstatus\tdisease\tdiseaseSubtype\tInheritance\tintegrate\tshare\tflag\tlocation\tcomment\n";
	
	#print Dumper($samples_infos_hash);

	# Create index to sort the file
	my %samples_infos_hash_sorted;
	if (1) {
	while ((my $runref, my $sample_list) = each(%$samples_infos_hash)) {
	while ((my $sampleref, my $val) = each(%$sample_list)) {
		if (trim($sampleref) ne "") {
			my $index=$$val{"runref"}."_".$$val{"disease"}."_".$$val{"family"}."_".$sampleref;
			$samples_infos_hash_sorted{$index}{$sampleref}=$$sample_list{$sampleref};
		};#if
	};#while
	};#while
	};#if

	#print Dumper(%samples_infos_hash_sorted);
	
	# Write the file
	foreach my $index (sort keys %samples_infos_hash_sorted) {
		
	while ((my $sampleref, my $val) = each(%{$samples_infos_hash_sorted{$index}})) {
	    #print "   $name\n";
	    print SAMPLES_INFOS $sampleref;
	    print SAMPLES_INFOS "\t".$$val{"projectref"};
	    print SAMPLES_INFOS "\t".$$val{"runref"};
	    print SAMPLES_INFOS "\t".$$val{"samplenumber"};
	    print SAMPLES_INFOS "\t".$$val{"gender"};
	    print SAMPLES_INFOS "\t".$$val{"family"};
	    print SAMPLES_INFOS "\t".$$val{"motherid"};
	    print SAMPLES_INFOS "\t".$$val{"fatherid"};
	    print SAMPLES_INFOS "\t".$$val{"relationship"};
	    print SAMPLES_INFOS "\t".$$val{"status"};
	    print SAMPLES_INFOS "\t".$$val{"disease"};
	    print SAMPLES_INFOS "\t".$$val{"diseasesubtype"};
	    print SAMPLES_INFOS "\t".$$val{"inheritance"};
	    print SAMPLES_INFOS "\t".$$val{"integrate"};
	    print SAMPLES_INFOS "\t".$$val{"share"};
	    print SAMPLES_INFOS "\t".$$val{"flag"};
	    print SAMPLES_INFOS "\t".$$val{"location"};
	    print SAMPLES_INFOS "\t".$$val{"comment"};
	    print SAMPLES_INFOS "\n";
	    
	    #while ((my $name2, my $val2) = each(%$val)) {
		#print "   $name2=$val2\n";
	    #};#while
	    #print "\n";
	};#while
	};#foreach
	#print "\n";

	close(SAMPLES_INFOS);
	
}

#Function create_project_ini
#Read an ini file
sub create_project_ini {
#$_[0]: project reference (usually the folder name)
#$_[1]: path of the file 'project.ini' (usually $dir_projects/$project_id)
#$_[2]: default plateform (usually 'GTF')
#$_[3]: default assembly (usually 'hg19')
#The 'new' project is not integrated by default (manual chek first)
	
	#Parameters
	my $project_ref=$_[0];
	my $path=$_[1];
	my $my_platform=$_[2]; #by default, change if in the 'project.ini' file
	my $my_assembly=$_[3]; #by default, change if in the 'project.ini' file
	$my_assembly =~ s/hg//g;
	my $date=timestamp("date");
	
	#Open the file
	open(PROJECT_INI, ">$path/project.ini") || die "Problem: $!";
	print PROJECT_INI "[Project]\n";
	print PROJECT_INI "Reference	$project_ref\n";
	print PROJECT_INI "Name	$project_ref\n";
	print PROJECT_INI "Description	Project $project_ref defined as a default project\n";
	print PROJECT_INI "Date	$date	; generated by default\n";
	print PROJECT_INI "Assembly	$my_assembly	; generated by default\n";
	print PROJECT_INI "Integrate	0	; change to 1 to integrate the project\n";
	print PROJECT_INI "Platform	$my_platform\n";
	print PROJECT_INI "[Contacts]\n";
	print PROJECT_INI "Recipient	Unknown\n";
	print PROJECT_INI "RecipientEmail	Unknown\n";
	print PROJECT_INI "RecipientAddress	Unknown\n";
	print PROJECT_INI "PI	Unknown\n";
	print PROJECT_INI "PIEmail	Unknown\n";
	print PROJECT_INI "PIAddress	Unknown\n";
	print PROJECT_INI "[Storage]\n";
	print PROJECT_INI "Location	$path\n";
	close(PROJECT_INI);
	return 1;
}


return 1;
