<?php


### INCLUDES
###############

include "config.php";


// Script start time
$rustart = getrusage();



### VARIABLES
###############

$PATH=$_GET["PATH"];
#$HOME="repositories";
$SEARCH=$_GET["search"];


$MAX_REPORTS=2000;

$VERBOSE=0;



### MODULES/CONTAINERS
#############

$MODULES=array();

# SERVER JARVIS CONTAINER NAME
$DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME=$_ENV["DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME"];
if ($DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME!="") {
	$MODULES["JARVIS"]=$DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME;
};

# SERVER IGV CONTAINER NAME
$DOCKER_STARK_SERVICE_IGV_CONTAINER_NAME=$_ENV["DOCKER_STARK_SERVICE_IGV_CONTAINER_NAME"];
if ($DOCKER_STARK_SERVICE_IGV_CONTAINER_NAME!="") {
	$MODULES["IGV"]=$DOCKER_STARK_SERVICE_IGV_CONTAINER_NAME;
};


### ENVIRONMENT
################

$DOCKER_STARK_SERVICE_DATA_SUBFOLDER_REPOSITORIES=$_ENV["DOCKER_STARK_SERVICE_DATA_SUBFOLDER_REPOSITORIES"];
if ($DOCKER_STARK_SERVICE_DATA_SUBFOLDER_REPOSITORIES=="") {
	$DOCKER_STARK_SERVICE_DATA_SUBFOLDER_REPOSITORIES="repositories";
};


### VARIABLES
###############


if ($PATH=="") {
	$PATH=$DOCKER_STARK_SERVICE_DATA_SUBFOLDER_REPOSITORIES;
}



### FUNCTIONS
###############

function path_full($PATH="") {
	$PATH_SPLIT=explode ( "/" , $PATH );
	if (isset($PATH_SPLIT[0]) && $PATH_SPLIT[0]!="") {$ROOT=$PATH_SPLIT[0];} else {$ROOT="*";};
	if (isset($PATH_SPLIT[1]) && $PATH_SPLIT[1]!="") {$REPOSITORY=$PATH_SPLIT[1];} else {$REPOSITORY="*";};
	if (isset($PATH_SPLIT[2]) && $PATH_SPLIT[2]!="") {$GROUP=$PATH_SPLIT[2];} else {$GROUP="*";};
	if (isset($PATH_SPLIT[3]) && $PATH_SPLIT[3]!="") {$PROJECT=$PATH_SPLIT[3];} else {$PROJECT="*";};
	if (isset($PATH_SPLIT[4]) && $PATH_SPLIT[4]!="") {$RUN=$PATH_SPLIT[4];} else {$RUN="*";};
	if (isset($PATH_SPLIT[5]) && $PATH_SPLIT[5]!="") {$SAMPLE=$PATH_SPLIT[5];} else {$SAMPLE="*";};
	return "$ROOT/$REPOSITORY/$GROUP/$PROJECT/$RUN/$SAMPLE/";
};

function path_short($PATH="") {
	return str_replace("//","/",str_replace("*/","",path_full($PATH)));
};

function path_html($PATH="") {
	# <a href='?PATH=$REPOSITORY/$GROUP/$PROJECT/$RUN/$SAMPLE/'>$SAMPLE</a>
	$PATH_SHORT=path_short($PATH);
	$PATH_HREF="";

	#$last_value="";
	foreach (explode("/",$PATH_SHORT) as $key => $value) {

		if ($value!="") {

			$PATH_HREF.="$value/";
			if ($key==0) {
				$PATH_RETURN="
					<big>
					<a class='' href='?PATH=$PATH_HREF'>
						<div class='card-img align-self-center mbr-bold'>
							<span class='mbr-iconfont mbri-home mbr-bold' style='color: rgb(20, 157, 204); fill: rgb(20, 157, 204);'></span>
						</div>
						HOME
					</a>
					</big>
				";
			} else {

				#$PATH_RETURN.=" > <a class='btn btn-primary-outline display-7' href='?PATH=$PATH_HREF'>$value</a>";
				$PATH_RETURN.="
				&nbsp;
				&nbsp;
				<a class='' href='?PATH=$PATH_HREF'>
					<div class='align-self-center bold'>
						<span class='mbr-iconfont mbri-arrow-next mbr-bold' style='color: rgb(20, 157, 204); fill: rgb(20, 157, 204);'></span>
						<span class='mbr-iconfont mbri-folder mbr-bold' style='color: rgb(20, 157, 204); fill: rgb(20, 157, 204);'></span>
					</div>
					&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$value
				</a>
				&nbsp;
				";
			};
			# card-img
		};
		#$last_value=$value;
	};

	if ($key<6) {
		$PATH_RETURN.="

		&nbsp;
		&nbsp;
		<a class='' href='?PATH=$PATH_HREF'>
			<div class='align-self-center bold'>

				<span class='mbr-iconfont mbri-arrow-down mbr-bold' style='color: rgb(20, 157, 204); fill: rgb(20, 157, 204);'></span>

				&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;".path_links($PATH)."
			</div>

		</a>
		&nbsp;

		";
	}

	# mbri-cust-feedback

	return $PATH_RETURN;
};

function path_links($PATH="") {
	$PATH_SHORT=path_short($PATH);
	if (count(explode( "/", $PATH ))<=6) {
		foreach (glob ( $PATH_SHORT."/*", GLOB_ONLYDIR ) as $key => $value) {
			$PATH_LINKS.="
				<a class='' href='?PATH="."$PATH_SHORT".end(explode( "/", $value ))."/'>
					<div class=' align-self-left align-left bold'>
						".end(explode( "/", $value ))."
					</div>
				</a>
			";
		};
	};
	#echo "L $LINKS_HTML L"; btn btn-primary-outline display-8 display-8
	return $PATH_LINKS;
};

function rutime($ru, $rus, $index) {
    return ($ru["ru_$index.tv_sec"]*1000 + intval($ru["ru_$index.tv_usec"]/1000))
     -  ($rus["ru_$index.tv_sec"]*1000 + intval($rus["ru_$index.tv_usec"]/1000));
}

function tags_extract($TAGS="") {
	$return=str_replace(" ","\n",str_replace("!"," ",$TAGS));
	preg_match_all('/.*#.*/i', $return, $matchWords);
	return implode(" ",$matchWords[0]);
    #  $matchWords = $matchWords[0];
    #  return $return;
};


### HEADER
############

echo '<!-- Site made with Mobirise Website Builder v4.10.3, https://mobirise.com -->
<meta charset="UTF-8">
<meta http-equiv="X-UA-Compatible" content="IE=edge">
<meta name="generator" content="Mobirise v4.10.3, mobirise.com">
<meta name="viewport" content="width=device-width, initial-scale=1, minimum-scale=1">
<link rel="shortcut icon" href="assets/favicon.ico" type="image/x-icon">
<meta name="description" content="STARK Report">
<title>STARK Reports</title>
<link rel="stylesheet" href="assets/web/assets/mobirise-icons/mobirise-icons.css">
<link rel="stylesheet" href="assets/tether/tether.min.css">
<link rel="stylesheet" href="assets/bootstrap/css/bootstrap.min.css">
<link rel="stylesheet" href="assets/bootstrap/css/bootstrap-grid.min.css">
<link rel="stylesheet" href="assets/bootstrap/css/bootstrap-reboot.min.css">
<link rel="stylesheet" href="assets/dropdown/css/style.css">
<link rel="stylesheet" href="assets/as-pie-progress/css/progress.min.css">
<link rel="stylesheet" href="assets/datatables/data-tables.bootstrap4.min.css">
<link rel="stylesheet" href="assets/theme/css/style.css">
<link rel="stylesheet" href="assets/gallery/style.css">
<link rel="stylesheet" href="assets/mobirise/css/mbr-additional.css" type="text/css">
<style>
	.div-wrapper {
		overflow: auto;
		max-height:400px;
		}
</style>
';

// <style>
// 	.table-wrapper {
// 		overflow: auto;
// 		max-height:800px;
// 		}
// </style>


echo '
<section class="menu cid-qTkzRZLJNu" once="menu" id="menu1-3">
  <nav class="navbar navbar-expand beta-menu navbar-dropdown align-items-center navbar-fixed-top collapsed bg-color transparent">
	  <button class="navbar-toggler navbar-toggler-right" type="button" data-toggle="collapse" data-target="#navbarSupportedContent" aria-controls="navbarSupportedContent" aria-expanded="false" aria-label="Toggle navigation">
		  <div class="hamburger">
			  <span></span>
			  <span></span>
			  <span></span>
			  <span></span>
		  </div>
	  </button>
	  <div class="menu-logo">
		  <div class="navbar-brand">
			  <span class="navbar-logo">
				  <a href="/">
					   <img src="assets/logo.png" alt="Mobirise" title="" style="height: 6rem;">
				  </a>
			  </span>
			  <span class="navbar-caption-wrap"><a class="navbar-caption text-secondary display-2" href="">Reports</a></span>
		 </div>
	  </div>
	  <div class="collapse navbar-collapse align-center" id="navbarSupportedContent">
	  <!--
		  <p class="mbr-text pb-3 mbr-fonts-style display-5">
			  <img src="assets/logo.png" alt="STARK" title="" style="height: 24rem;">
		  </p>
		  <h1 class="mbr-section-title mbr-bold pb-3 mbr-fonts-style display-2">
			  $ENV_NAME
		  </h1>
		  <p class="mbr-text pb-3 mbr-fonts-style display-5">
			  $ENV_RELEASE - $ENV_DATE
			  <BR>
			  $ENV_DESCRIPTION
			  <BR>
			  © Copyright $ENV_COPYRIGHT - All Rights Reserved
			  <BR>
			  $ENV_LICENCE Licence
			  <BR>
			  <BR>
			  <BR>
		  </p>
		-->
	  </div>
  </nav>
</section>
';



$PATH_HTML=path_html($PATH);
$PATH_LINKS=path_links($PATH);

echo '
	<section class="header1 cid-ru7OEConn1" id="header16-1k">
	    <div class="container">
			<h3 class="mbr-section-subtitle mbr-fonts-style align-center mbr-light display-2">
				&nbsp;
		  	</h3>
	    </div>
		<div class="container div-wrapper">

			<h3 class="mbr-section-subtitle mbr-fonts-style align-left mbr-light display-5">
				<small><small>
				<div class="media">
					'.$PATH_HTML.'
				</div>
				<br>

				<div class="div-wrapper">

				</div>

				</small></small>
		  	</h3>

	    </div>
	</section>

';

#<div class="mbr-section-btn align-center"> <div class="media"> '.$PATH_LINKS.'


### TABLE
###########

### THEAD

$thead='
	<tr class="table-heads">
		<th class="head-item mbr-fonts-style display-7">
			Reports
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Sample
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Analysis/Run
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Project
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Group
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Repository
		</th>
	</tr>
';


### TBODY


$root="repositories";

$tbody="";
$reports=glob ( "".path_full($PATH)."*stark.report.html" );


if (count($reports)<=$MAX_REPORTS) {

	$tags_sample_files=glob ( "".path_full($PATH)."*tag" );
	$tags_analysis_files=glob ( "".path_full($PATH)."*analysis.tag" );

	foreach ($reports as $key => $report_html) {

		# Find infos
		$report_html_split=explode ( "/" , $report_html );
		$root=$report_html_split[0];
		$repository=$report_html_split[1];
		$group=$report_html_split[2];
		$project=$report_html_split[3];
		$run=$report_html_split[4];
		$sample=$report_html_split[5];
		$report_html_file=$report_html_split[6];
		$report_html_file_split=explode ( "." , $report_html_file );
		$report_html_id=$report_html_file_split[1];

		# Files
		$sample_path=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run.'/'.$sample;
		$run_path=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run;
		$project_path=$root.'/'.$repository.'/'.$group.'/'.$project;
		$group_path=$root.'/'.$repository.'/'.$group;
		$report_tsv=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run.'/'.$sample.'/'.$sample.'.final.tsv';
		$report_vcf_gz=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run.'/'.$sample.'/'.$sample.'.final.vcf.gz';
		$report_bed=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run.'/'.$sample.'/'.$sample.'.bed';




		# TAGS
		$tags_file=glob ( "".$sample_path."/".$sample.".tag" );
		$analysis_tags_file=glob ( "".$sample_path."/".$sample.".analysis.tag" );
		if (count($tags_file)>0) {

			$myfile = fopen($tags_file[0], "r") or die("Unable to open file!");
			$tags=tags_extract(trim(fread($myfile,filesize($tags_file[0]))));
			fclose($myfile);
		};
		if (count($analysis_tags_file)>0) {
			$myfile = fopen($analysis_tags_file[0], "r") or die("Unable to open file!");
			$tags_analysis=tags_extract(trim(fread($myfile,filesize($analysis_tags_file[0]))));
			fclose($myfile);
		};
		if ($tags != "") {
			$tags="<br><small><span style='color:gray'>$tags</span></small>";
		};
		if ($tags_analysis != "") {
			$tags_analysis="<br><small><span style='color:gray'>$tags_analysis</span></small>";
		};

		# IGV LINKS
		if (isset($MODULES["IGV"])) {
			$IGV_SAMPLE_LINK="<a target='IGV' href='index.igv.php?SEARCH_LEVEL=0&PROCESS=1&SAMPLE_PATH[]=$sample_path'>IGV</a>";
			$IGV_RUN_LINK="<a target='IGV' href='index.igv.php?SEARCH_LEVEL=1&PROCESS=0&SAMPLE_PATH[]=$run_path'>IGV</a>";
		}

		# JARVIS LINKS
		if (isset($MODULES["JARVIS"])) {
			$JARVIS_SAMPLE_LINK="<a target='JARVIS' href='index.jarvis.php?SEARCH_LEVEL=0&PROCESS=1&SAMPLE_PATH[]=$sample_path'>JARVIS</a>";
			$JARVIS_RUN_LINK="<a target='JARVIS' href='index.jarvis.php?SEARCH_LEVEL=1&SEARCH_EXT=final.vcf.gz,final.bcf&PROCESS=0&SAMPLE_PATH[]=$run_path'>JARVIS</a>";
			$JARVIS_PROJECT_LINK="<a target='JARVIS' href='index.jarvis.php?SEARCH_LEVEL=2&SEARCH_EXT=final.vcf.gz,final.bcf&PROCESS=0&SAMPLE_PATH[]=$project_path'>JARVIS</a>";
			$JARVIS_GROUP_LINK="<a target='JARVIS' href='index.jarvis.php?SEARCH_LEVEL=3&SEARCH_EXT=final.vcf.gz,final.bcf&PROCESS=0&SAMPLE_PATH[]=$group_path'>JARVIS</a>";
		}

		# METRICS LINKS
		$METRICS_RUN_LINK="<a target='METRICS' href='index.metrics.php?PATH[]=$run_path'>METRICS</a>";
		$METRICS_PROJECT_LINK="<a target='METRICS' href='index.metrics.php?SEARCH_LEVEL=1&PATH[]=$project_path'>METRICS</a>";
		$METRICS_GROUP_LINK="<a target='METRICS' href='index.metrics.php?SEARCH_LEVEL=2&PATH[]=$group_path'>METRICS</a>";

		# TBODY
		$tbody=$tbody.'
			<tr class="table-heads">
				<td class="head-item mbr-fonts-style display-7">
					<small>
						<a href="'.$report_html.'"><b>'.$report_html_id.'</b></a>
						<br>
						<a href="'.$report_tsv.'" download>TSV</a>
						<a href="'.$report_vcf_gz.'" download>VCF</a>
						<a href="'.$report_bed.'" download>BED</a>
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						<b>'.$sample.'</b>
						'.$tags.'
						<br>
						'.$IGV_SAMPLE_LINK.'
						'.$JARVIS_SAMPLE_LINK.'
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						'.$run.'
						'.$tags_analysis.'
						<br>
						'.$IGV_RUN_LINK.'
						'.$JARVIS_RUN_LINK.'
						'.$METRICS_RUN_LINK.'
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						'.$project.'
						<br>
						'.$JARVIS_GROUP_LINK.'
						'.$METRICS_PROJECT_LINK.'
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						'.$group.'
						<br>
						'.$JARVIS_GROUP_LINK.'
						'.$METRICS_GROUP_LINK.'
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						'.$repository.'
						<br>
						'.$JARVIS_REPOSITORY_LINK.'
					</small>
				</td>
			</tr>
		';
	};


	### SECTION
	#############

	echo '

	<section class="section-table cid-ruiBoanIwc" id="REPORTS">

	  <div class="container container-table">

		  <div class="table-wrapper">

			<div class="container">
			  <div class="row search">
				<div class="col-md-6"></div>
				<div class="col-md-6">
					<div class="dataTables_filter">
					  <label class="searchInfo mbr-fonts-style display-7">Search:</label>
					  <input class="form-control input-sm" disabled="">
					</div>
				</div>
			  </div>
			</div>

			<div class="container scroll">
				<table class="table isSearch" cellspacing="0">
					<thead>
						'.$thead.'
					</thead>
					<tbody>
						'.$tbody.'
					</tbody>
				</table>
			</div>

			<div class="container table-info-container">
			  <div class="row info">
				<div class="col-md-6">
				  <div class="dataTables_info mbr-fonts-style display-7">
					<span class="infoBefore">Showing</span>
					<span class="inactive infoRows"></span>
					<span class="infoAfter">entries</span>
					<span class="infoFilteredBefore">(filtered from</span>
					<span class="inactive infoRows"></span>
					<span class="infoFilteredAfter"> total entries)</span>
				  </div>
				</div>
				<div class="col-md-6"></div>
			  </div>
			</div>

		  </div>

		</div>

	</section>
	';

} else {

	echo '

	<section class="section-table cid-ruiBoanIwc" id="REPORTS">

	  <p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
	  	Too many reports selected. Please browse folders...
	  </p>

	</section>
	';


};

### SCRIPTS
#############

echo '
	<script src="assets/web/assets/jquery/jquery.min.js"></script>
	<script src="assets/popper/popper.min.js"></script>
	<script src="assets/tether/tether.min.js"></script>
	<script src="assets/bootstrap/js/bootstrap.min.js"></script>
	<script src="assets/smoothscroll/smooth-scroll.js"></script>
	<script src="assets/dropdown/js/script.min.js"></script>
	<script src="assets/touchswipe/jquery.touch-swipe.min.js"></script>
	<script src="assets/viewportchecker/jquery.viewportchecker.js"></script>
	<script src="assets/as-pie-progress/jquery-as-pie-progress.min.js"></script>
	<script src="assets/vimeoplayer/jquery.mb.vimeo_player.js"></script>
	<script src="assets/datatables/jquery.data-tables.min.js"></script>
	<script src="assets/datatables/data-tables.bootstrap4.min.js"></script>
	<script src="assets/masonry/masonry.pkgd.min.js"></script>
	<script src="assets/imagesloaded/imagesloaded.pkgd.min.js"></script>
	<script src="assets/bootstrapcarouselswipe/bootstrap-carousel-swipe.js"></script>
	<script src="assets/mbr-switch-arrow/mbr-switch-arrow.js"></script>
	<script src="assets/theme/js/script.js"></script>
	<script src="assets/gallery/player.min.js"></script>
	<script src="assets/gallery/script.js"></script>
	<script src="assets/slidervideo/script.js"></script>
	<script src="assets/mbr-tabs/mbr-tabs.js"></script>
';


# Script end time

if ($VERBOSE) {

	$ru = getrusage();
	echo "This process used " . rutime($ru, $rustart, "utime") .
	    " ms for its computations\n";
	echo "It spent " . rutime($ru, $rustart, "stime") .
	    " ms in system calls\n";

};




?>
