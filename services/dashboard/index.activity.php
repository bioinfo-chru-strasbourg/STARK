<?php


	### VARIABLES
	###############

	#$count=" count ";
	$count="  ";

	### FUNCTIONS
	###############


	function pp($arr){
	    $retStr = '<ul>';
	    if (is_array($arr)){
	        foreach ($arr as $key=>$val){
	            if (is_array($val)){
	                $retStr .= '<li>' . $key . ' => ' . pp($val) . '</li>';
	            }else{
	                $retStr .= '<li>' . $key . ' => ' . $val . '</li>';
	            }
	        }
	    }
	    $retStr .= '</ul>';
	    return $retStr;
	}



	### HEADER
	############

	echo '<!-- Site made with Mobirise Website Builder v4.10.3, https://mobirise.com -->
	<meta charset="UTF-8">
	<meta http-equiv="X-UA-Compatible" content="IE=edge">
	<meta name="generator" content="Mobirise v4.10.3, mobirise.com">
	<meta name="viewport" content="width=device-width, initial-scale=1, minimum-scale=1">
	<link rel="shortcut icon" href="assets/favicon.ico" type="image/x-icon">
	<meta name="description" content="STARK Activity">
	<title>STARK Activity</title>
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
				  <span class="navbar-caption-wrap"><a class="navbar-caption text-secondary display-2" href="">Activity</a></span>
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


	#$DOCKER_STARK_SERVICE_LAUNCHER_CONTAINER_NAME:$DOCKER_STARK_SERVICE_LAUNCHER_PORT
	#echo readfile("http://$DOCKER_STARK_SERVICE_LAUNCHER_CONTAINER_NAME:$DOCKER_STARK_SERVICE_LAUNCHER_PORT/queue?format=html");
	#echo "http://".$_ENV["DOCKER_STARK_SERVICE_LAUNCHER_CONTAINER_NAME"].":".$_ENV["DOCKER_STARK_SERVICE_LAUNCHER_PORT_INNER"]."/queue";
	#$truc=readfile("http://$DOCKER_STARK_SERVICE_LAUNCHER_CONTAINER_NAME:$DOCKER_STARK_SERVICE_LAUNCHER_PORT");
	#echo "<pre>";
	$URL_LAUNCHER="http://".$_ENV["DOCKER_STARK_SERVICE_LAUNCHER_CONTAINER_NAME"].":".$_ENV["DOCKER_STARK_SERVICE_LAUNCHER_PORT_INNER"]."/queue";

	$lines_array=file($URL_LAUNCHER);
	print_r($lines_array);
	#$url='https://oscarliang.com';
	// using file() function to get content
	#$lines_array=file($url);
	// turn array into one variable
	#print_r($lines_array);
	#$lines_string=implode('',$lines_array);
	//output, you can also save it locally on the server
	#echo $lines_string;


	$BROWSER_SERVER=$_SERVER["SERVER_NAME"];
	$BROWSER_PORT=$_ENV["DOCKER_STARK_SERVICE_PORT_PATTERN"].$_ENV["DOCKER_STARK_SERVICE_BROWSER_PORT"];
	$BROWSER_URL=$BROWSER_SERVER.":".$BROWSER_PORT;


	# REPOSITORIES
	$REPOSITORY=$_ENV["DOCKER_STARK_INNER_FOLDER_OUTPUT_REPOSITORY"];
	$ARCHIVE=$_ENV["DOCKER_STARK_INNER_FOLDER_OUTPUT_ARCHIVE"];
	#print_r(glob ( $REPOSITORY."/*/*/*/*", GLOB_ONLYDIR ));
	foreach (array($REPOSITORY,$ARCHIVE) as $REPO) {
		foreach (glob ( $REPO."/*/*/*/*", GLOB_ONLYDIR ) as $key => $sample_path) {
			$sample_split=explode ( "/" , $sample_path );
			$sample_path_count=count($sample_split);
			#print $sample_path_count;
			#$repository=$sample_split[0];
			$group=$sample_split[$sample_path_count-4];
			$project=$sample_split[$sample_path_count-3];
			$run=$sample_split[$sample_path_count-2];
			$sample=$sample_split[$sample_path_count-1];
			#echo $sample;
			$hash["groups"][$group]++;
			$hash["projects"][$project]++;
			$hash["runs"][$run]++;
			$hash["samples"][$sample]++;
			$hash["tree"][$group][$project][$run][$sample]++;
		};
	};
	$nb_groups=count($hash["groups"]);
	$nb_projects=count($hash["projects"]);
	$nb_runs=count($hash["runs"]);
	$nb_samples=count($hash["samples"]);


	echo '

	<section class="counters1 counters cid-ru7OEH6r7i" id="SECTION_ID" >

		<div class="container">

		   <h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
			   Global Activity
		   </h2>

		   <p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
			  Activity Statistics
		   </p>

		   <h3 class="mbr-section-subtitle mbr-fonts-style align-center mbr-light display-2">
			test
		   </h3>

		   <div class="container pt-4 mt-2">
			<div class="media-container-row">

				<div class="card p-3 align-center col-12 col-md-6 col-lg-2">
					<div class="panel-item">
						<div class="card-img pb-3">
							<span class="mbr-iconfont mbri-align-center"></span>
						</div>

						<div class="card-text">
							<h3 class="'.$count.' pt-3 pb-3 mbr-fonts-style display-2">
								'.count($runs[$input]).'
							</h3>
							<h4 class="mbr-content-title mbr-bold mbr-fonts-style display-7">
								Runs
							</h4>
						</div>
					</div>
				</div>

				<div class="card p-3 align-center col-12 col-md-6 col-lg-2">
					<div class="panel-item">
						<div class="card-img pb-3">
							<span class="mbr-iconfont mbri-target"></span>
						</div>

						<div class="card-text">
							<h3 class="'.$count.' pt-3 pb-3 mbr-fonts-style display-2">
								'.count($designs[$input]).'
							</h3>
							<h4 class="mbr-content-title mbr-bold mbr-fonts-style display-7">
								Target Designs
							</h4>
						</div>
					</div>
				</div>

				<div class="card p-3 align-center col-12 col-md-6 col-lg-2">
					<div class="panel-item">
						<div class="card-img pb-3">
							<span class="mbr-iconfont mbri-menu"></span>
						</div>

						<div class="card-text">
							<h3 class="'.$count.' pt-3 pb-3 mbr-fonts-style display-2">
								'.count($genes[$input]).'
							</h3>
							<h4 class="mbr-content-title mbr-bold mbr-fonts-style display-7">
								Gene Panels
							</h4>
						</div>
					</div>
				</div>

				<div class="card p-3 align-center col-12 col-md-6 col-lg-2">
					<div class="panel-item">
						<div class="card-img pb-3">
							<span class="mbr-iconfont mbri-bookmark"></span>
						</div>

						<div class="card-text">
							<h3 class="'.$count.' pt-3 pb-3 mbr-fonts-style display-2">
								'.count($transcripts[$input]).'
							</h3>
							<h4 class="mbr-content-title mbr-bold mbr-fonts-style display-7">
								Transcript custom lists
							</h4>
						</div>
					</div>
				</div>

			</div>
		   </div>

		   <p class="mbr-content-text mbr-fonts-style display-7">
			   &nbsp;
		  </p>

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
