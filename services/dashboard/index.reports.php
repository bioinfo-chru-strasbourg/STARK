<?php


	### VARIABLES
	###############

	$PATH=$_GET["PATH"];
	$HOME="repositories";
	$SEARCH=$_GET["search"];


	### FUNCTIONS
	###############

	function path_full($PATH="") {
		$PATH_SPLIT=explode ( "/" , $PATH );
		if (isset($PATH_SPLIT[0]) && $PATH_SPLIT[0]!="") {$REPOSITORY=$PATH_SPLIT[0];} else {$REPOSITORY="*";};
		if (isset($PATH_SPLIT[1]) && $PATH_SPLIT[1]!="") {$GROUP=$PATH_SPLIT[1];} else {$GROUP="*";};
		if (isset($PATH_SPLIT[2]) && $PATH_SPLIT[2]!="") {$PROJECT=$PATH_SPLIT[2];} else {$PROJECT="*";};
		if (isset($PATH_SPLIT[3]) && $PATH_SPLIT[3]!="") {$RUN=$PATH_SPLIT[3];} else {$RUN="*";};
		if (isset($PATH_SPLIT[4]) && $PATH_SPLIT[4]!="") {$SAMPLE=$PATH_SPLIT[4];} else {$SAMPLE="*";};
		return "$REPOSITORY/$GROUP/$PROJECT/$RUN/$SAMPLE/";
	};

	function path_short($PATH="") {
		return str_replace("*/","",path_full($PATH));
	};

	function path_html($HOME="repositories",$PATH="") {
		# <a href='?PATH=$REPOSITORY/$GROUP/$PROJECT/$RUN/$SAMPLE/'>$SAMPLE</a>
		$PATH_SHORT=path_short($PATH);
		$PATH_HREF="";
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
		#$last_value="";
		foreach (explode("/",$PATH_SHORT) as $key => $value) {
			if ($value!="") {
				$PATH_HREF.="$value/";
				#$PATH_RETURN.=" > <a class='btn btn-primary-outline display-7' href='?PATH=$PATH_HREF'>$value</a>";
				$PATH_RETURN.="
				&nbsp;
				<!--
				<div class='card-img align-self-center bold'>
					<span class='mbr-iconfont mbri-arrow-next mbr-bold' style='color: rgb(20, 157, 204); fill: rgb(20, 157, 204);'></span>
				</div>
				-->
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
				# card-img
			};
			#$last_value=$value;
		};

		if ($key<5) {
			$PATH_RETURN.="

			&nbsp;
			<!--
			<div class='card-img align-self-center bold'>
				<span class='mbr-iconfont mbri-arrow-next mbr-bold' style='color: rgb(20, 157, 204); fill: rgb(20, 157, 204);'></span>
			</div>
			-->
			&nbsp;
			<a class='' href='?PATH=$PATH_HREF'>
				<div class='align-self-center bold'>

					<span class='mbr-iconfont mbri-arrow-down mbr-bold' style='color: rgb(20, 157, 204); fill: rgb(20, 157, 204);'></span>

					&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;".path_links($HOME,$PATH)."
				</div>

			</a>
			&nbsp;

			";
		}

		# mbri-cust-feedback

		return $PATH_RETURN;
	};

	function path_links($HOME="repositories",$PATH="") {
		$PATH_SHORT=path_short($PATH);
		if (count(explode( "/", $PATH ))<=5) {
			foreach (glob ( $HOME."/".$PATH_SHORT."/*", GLOB_ONLYDIR ) as $key => $value) {
				// $PATH_LINKS.="
				// 	&nbsp;
				// 	<small><small>
				// 	<a class='btn display-8 small' href='?PATH="."$PATH".end(explode( "/", $value ))."/'>
				// 		<div class='card-img align-self-center'>
				// 			<span class='mbr-iconfont mbri-folder mbr-bold' style='color: rgb(20, 157, 204); fill: rgb(20, 157, 204);'></span>
				// 				<br><small>".end(explode( "/", $value ))."</small>
				// 		</div>
				//
				// 	</a>
				// 	</small></small>
				// 	&nbsp;<br>
				// 	";
				$PATH_LINKS.="
					<a class='' href='?PATH="."$PATH".end(explode( "/", $value ))."/'>
						<div class=' align-self-left align-left bold'>
							<!--
							<span class='mbr-iconfont mbri-folder mbr-bold' style='color: rgb(20, 157, 204); fill: rgb(20, 157, 204);'></span>
							&nbsp;&nbsp;
							-->
							".end(explode( "/", $value ))."
							<!--
							&nbsp;&nbsp;
							<span class='mbr-iconfont mbri-cursor-click mbr-bold' style='color: rgb(20, 157, 204); fill: rgb(20, 157, 204);'></span>
							-->

						</div>
					</a>
				";
			};
		};
		#echo "L $LINKS_HTML L"; btn btn-primary-outline display-8 display-8
		return $PATH_LINKS;
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






	$PATH_HTML=path_html($HOME,$PATH);
	$PATH_LINKS=path_links($HOME,$PATH);

	echo '
		<section class="header1 cid-ru7OEConn1" id="header16-1k">
			<!--
			<div class="container">
			   	<div class="row justify-content-center  align-center">
			   		<div class="card p-6 col-12 col-md-12">
			   			<div class="card-box">
			   				<h1 class="mbr-section-title mbr-bold pb-3 mbr-fonts-style display-2">
				   				<img src="assets/logo.png" width="128">
			   				</h1>
			   			</div>
			   		</div>
				</div>
			</div>
			-->
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
		</tr>
	';


	### TBODY

	$tbody="";
	$reports=glob ( "repositories/".path_full($PATH)."*stark.report.html" );
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
		$report_tsv=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run.'/'.$sample.'/'.$sample.'.final.tsv';
		$report_vcf_gz=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run.'/'.$sample.'/'.$sample.'.final.vcf.gz';
		$report_bed=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run.'/'.$sample.'/'.$sample.'.bed';
		# TBODY
		$tbody=$tbody.'
			<tr class="table-heads">
				<td class="head-item mbr-fonts-style display-7">
					<a href="'.$report_html.'">'.$report_html_id.'</a>
					<br>
					<a href="'.$report_tsv.'">TSV</a>
					<a href="'.$report_vcf_gz.'">VCF</a>
					<a href="'.$report_bed.'">BED</a>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					'.$sample.'
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>'.$run.'</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>'.$project.'</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>'.$group.'</small>
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

?>
