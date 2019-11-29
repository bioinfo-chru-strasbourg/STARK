<?php


	### VARIABLES
	###############


	### FUNCTIONS
	###############


	### HEADER
	############

	echo '<!-- Site made with Mobirise Website Builder v4.10.3, https://mobirise.com -->
	<meta charset="UTF-8">
	<meta http-equiv="X-UA-Compatible" content="IE=edge">
	<meta name="generator" content="Mobirise v4.10.3, mobirise.com">
	<meta name="viewport" content="width=device-width, initial-scale=1, minimum-scale=1">
	<link rel="shortcut icon" href="assets/favicon.ico" type="image/x-icon">
	<meta name="description" content="STARK Dashboard">
	<title>STARK Dashboard</title>
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
				  <span class="navbar-caption-wrap"><a class="navbar-caption text-secondary display-2" href="">Dashboard</a></span>
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

		</section>

	';

	#<div class="mbr-section-btn align-center"> <div class="media"> '.$PATH_LINKS.'


	### CONTENT
	#############


	// echo '
	// <pre>
	// '.print_r($_SERVER).'
	// </pre>
	// ';

	$BROWSER_SERVER=$_SERVER["SERVER_NAME"];
	$BROWSER_PORT=$_ENV["DOCKER_STARK_SERVICE_PORT_PATTERN"].$_ENV["DOCKER_STARK_SERVICE_BROWSER_PORT"];
	$BROWSER_URL=$BROWSER_SERVER.":".$BROWSER_PORT;

	# class='btn btn-primary-outline display-7'

	echo '



	<section class="features10 cid-ru7OEDbxhA" id="features10-1l">

		<div class="container ">
			<div class="row justify-content-center">

				<div class="card p-3 col-12 col-md-6 mb-4">
					<a href="http://'.$BROWSER_URL.'" class="navbar-caption text-secondary">
						<div class="media mb-2">
							<div class="card-img align-self-center">
								<span class="mbr-iconfont mbri-browse" style="color: rgb(20, 157, 204); fill: rgb(20, 157, 204);"></span>
							</div>
							<h4 class="card-title media-body py-3 mbr-fonts-style display-5">Data Browser</h4>
						</div>
						<div class="card-box">
							<p class="mbr-text mbr-fonts-style display-7">
								Data Browser shares available data, such as raw sequenced data with configuration files (manifest, bed, genes...),
								analyses results in the repository,
								and archived data.
							</p>
						</div>
					</a>
				</div>

				<div class="card p-3 col-12 col-md-6 mb-4">
					<a href="index.reports.php" class="navbar-caption text-secondary ">
						<div class="media mb-2">
							<div class="card-img align-self-center">
								<span class="mbr-iconfont mbri-cust-feedback" style="color: rgb(20, 157, 204); fill: rgb(20, 157, 204);"></span>
							</div>
							<h4 class="card-title media-body py-3 mbr-fonts-style display-5">Reports Browser</h4>
						</div>
						<div class="card-box">
							<p class="mbr-text mbr-fonts-style display-7">
								Reports Browser provides a direct access to STARK reports (validation report, VCF, TSV, BED...) availabled in reporitory.
							</p>
						</div>
					</a>
				</div>

			</div>
		</div>

		<div class="container">
			<div class="row justify-content-center">

				<div class="card p-3 col-12 col-md-6 mb-4">
					<a href="index.stats.php" class="navbar-caption text-secondary ">
						<div class="media mb-2">
							<div class="card-img align-self-center">
								<span class="mbr-iconfont mbri-growing-chart" style="color: rgb(20, 157, 204); fill: rgb(20, 157, 204);"></span>
							</div>
							<h4 class="card-title media-body py-3 mbr-fonts-style display-5">Statistics</h4>
						</div>
						<div class="card-box">
							<p class="mbr-text mbr-fonts-style display-7">
								Statistics on runs, groups, projects, analyses and samples,
								stored in repository, archives and and sequencing raw sequencing data.
							</p>
						</div>
					</a>
				</div>

				<div class="card p-3 col-12 col-md-6 mb-4">
					<a href="index.reports.php" class="navbar-caption text-secondary ">

					</a>
				</div>

			</div>
		</div>

	</section>






	';

	# mbri-cust-feedback mbri-desktop

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
