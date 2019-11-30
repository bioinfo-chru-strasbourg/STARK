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
	<meta name="description" content="STARK Statistics">
	<title>STARK Statisticss</title>
	<script src="assets/plotly-latest.min.js"></script>
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
				  <span class="navbar-caption-wrap"><a class="navbar-caption text-secondary display-2" href="">Statistics</a></span>
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


	### CONTENT
	#############

	$INPUT="INPUT";
	$input="input";

	# INPUTS
	$inputs_folder="inputs";
	foreach (glob ( $inputs_folder."/*/*/*" ) as $key => $input_path) {
		$run_split=explode ( "/" , $input_path );
		$run_path_count=count($run_split);
		$input=$run_split[$run_path_count-3];
		$input_type=$run_split[$run_path_count-2];
		$input_file=$run_split[$run_path_count-1];
		$input_ext=pathinfo($input_path,PATHINFO_EXTENSION);
		#echo "<br>$input | $input_type | $input_file | $input_ext";
		$input_list[$input]++;
		if ($input_type=="manifests" && is_file($input_path)) {
			if ($input_ext=="genes") {
				$genes[$input][$input_file]++;
			} elseif ($input_ext=="transcripts") {
				$transcripts[$input][$input_file]++;
			} else {
				$designs[$input][$input_file]++;
			};
		} elseif ($input_type=="runs" && is_dir($input_path)) {
			$runs[$input][$input_file]++;
		};
	};


	# REPOSITORIES
	$repositories_folder="repositories";
	foreach (glob ( $repositories_folder."/*/*/*/*/*", GLOB_ONLYDIR ) as $key => $sample_path) {
		$sample_split=explode ( "/" , $sample_path );
		$sample_path_count=count($sample_split);
		#print $sample_path_count;
		$repository=$sample_split[$sample_path_count-5];
		$group=$sample_split[$sample_path_count-4];
		$project=$sample_split[$sample_path_count-3];
		$run=$sample_split[$sample_path_count-2];
		$sample=$sample_split[$sample_path_count-1];
		#echo $sample;
		$global[$repository]["groups"][$group]++;
		$global[$repository]["projects"][$project]++;
		$global[$repository]["runs"][$run]++;
		$global[$repository]["samples"][$sample]++;
		$global[$repository]["tree"][$group][$project][$run][$sample]++;

		$run_by_group[$repository][$group][$run]++;
		$sample_by_group[$repository][$group][$sample]++;
		$run_by_group_project[$repository][$group][$project][$run]++;
		$sample_by_group_project[$repository][$group][$project][$sample]++;

	};

	# NB variants
	#echo shell_exec("cut -f1-5 $repositories_folder/*/*/*/*/*/*final.tsv | sort -u | wc -l");


	echo '
		<section class="counters1 counters cid-ru7OEH6r7i" id="SECTION_ID" >

			<div class="container">

			   <h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
				   Global Statistics
			   </h2>

			   <p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
				   Statistics for inputs and repositories folders (duplication not considered)
			   </p>
	';


	foreach ($input_list as $input=>$hash) {

		echo '

			   <h3 class="mbr-section-subtitle mbr-fonts-style align-center mbr-light display-2">
			   	'.$input.'
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

		';

		# add-submenu align-justify bookmark
		# transcripts : bookmark
		# SAPLE : mbri-code

	};

	foreach ($global as $repository=>$hash) {

		$nb_groups=count($hash["groups"]);
		$nb_projects=count($hash["projects"]);
		$nb_runs=count($hash["runs"]);
		$nb_samples=count($hash["samples"]);

		echo '

			   <h3 class="mbr-section-subtitle mbr-fonts-style align-center mbr-light display-2">
			   	'.$repository.'
			   </h3>

			   <div class="container pt-4 mt-2">
			   	<div class="media-container-row">

			   		<div class="card p-3 align-center col-12 col-md-6 col-lg-2">
			   			<div class="panel-item">
			   				<div class="card-img pb-3">
			   					<span class="mbr-iconfont mbri-folder"></span>
			   				</div>

			   				<div class="card-text">
			   					<h3 class="'.$count.' pt-3 pb-3 mbr-fonts-style display-2">
			   						'.$nb_groups.'
			   					</h3>
			   					<h4 class="mbr-content-title mbr-bold mbr-fonts-style display-7">
			   						Groups
			   					</h4>
			   				</div>
			   			</div>
			   		</div>

					<div class="card p-3 align-center col-12 col-md-6 col-lg-2">
			   			<div class="panel-item">
			   				<div class="card-img pb-3">
			   					<span class="mbr-iconfont mbri-folder"></span>
			   				</div>

			   				<div class="card-text">
			   					<h3 class="'.$count.' pt-3 pb-3 mbr-fonts-style display-2">
			   						'.$nb_projects.'
			   					</h3>
			   					<h4 class="mbr-content-title mbr-bold mbr-fonts-style display-7">
			   						Projects
			   					</h4>
			   				</div>
			   			</div>
			   		</div>

					<div class="card p-3 align-center col-12 col-md-6 col-lg-2">
			   			<div class="panel-item">
			   				<div class="card-img pb-3">
			   					<span class="mbr-iconfont mbri-align-center"></span>
			   				</div>

			   				<div class="card-text">
			   					<h3 class="'.$count.' pt-3 pb-3 mbr-fonts-style display-2">
			   						'.$nb_runs.'
			   					</h3>
			   					<h4 class="mbr-content-title mbr-bold mbr-fonts-style display-7">
			   						Analyses / Runs
			   					</h4>
			   				</div>
			   			</div>
			   		</div>

					<div class="card p-3 align-center col-12 col-md-6 col-lg-2">
			   			<div class="panel-item">
			   				<div class="card-img pb-3">
			   					<span class="mbr-iconfont mbri-code"></span>
			   				</div>

			   				<div class="card-text">
			   					<h3 class="'.$count.' pt-3 pb-3 mbr-fonts-style display-2">
			   						'.$nb_samples.'
			   					</h3>
			   					<h4 class="mbr-content-title mbr-bold mbr-fonts-style display-7">
			   						Samples
			   					</h4>
			   				</div>
			   			</div>
			   		</div>

			   	</div>
			   </div>

			   <p class="mbr-content-text mbr-fonts-style display-7">
			       &nbsp;
			  </p>

		';

	};

	echo '
			</div>
		</section>
	';


	# Statistics by group

	$links_section.='
		<section class=" cid-rtQc0shhyd" id="truc">

			<div class="container">

			   <h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
				   Groups Statistics
			   </h2>

			   <p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
				   Statistics on analyses/runs and samples for each projects in groups (duplication not considered)
			   </p>
	';

	foreach ($global as $repository=>$hash) {

		# Analyses/Runs and Samples by group
		$run_by_group_nb=array();
		$sample_by_group_nb=array();
		$run_by_group_project_nb=array();
		$sample_by_group_project_nb=array();
		foreach ($run_by_group[$repository] as $group=>$hash2) {
			$run_by_group_nb[$group]=count($hash2);
			foreach ($run_by_group_project[$repository][$group] as $project=>$hash3) {
				$run_by_group_project_nb[$group][$project]=count($hash3);
			};
		};
		foreach ($sample_by_group[$repository] as $group=>$hash2) {
			$sample_by_group_nb[$group]=count($hash2);
			foreach ($sample_by_group_project[$repository][$group] as $project=>$hash3) {
				$sample_by_group_project_nb[$group][$project]=count($hash3);
			};
		};
		$group_list="'".implode("','",array_keys($run_by_group_nb))."'";
		$run_by_group_list=implode(",",$run_by_group_nb);
		$sample_by_group_list=implode(",",$sample_by_group_nb);


		$run_by_group_html = "
		<div class='section_div'>
			<div id='".$repository."run_by_group_list'>
				<div id='".$repository."run_by_group_list_figure'>
					<div class='figure' id='".$repository."plot_run_by_group_list' style='height:400px;width:500px;'></div>
				</div>

				<script type='text/javascript'>
					var data = [{
						  values: [".$run_by_group_list."],
						  labels: [".$group_list."],
						  type: 'pie'
						}];
					var layout={title:'Number of Analyses/Runs by group ', xaxis:{title:'duplication level'}, yaxis:{title:'Read percent (%) & GC ratio'}};
					Plotly.newPlot('".$repository."plot_run_by_group_list', data, layout);
				</script>

			</div>
		</div>
		";

		$sample_by_group_html = "
			<div class='section_div'>
				<div id='".$repository."sample_by_group_list'>
					<div id='".$repository."sample_by_group_list_figure'>
						<div class='figure' id='".$repository."plot_sample_by_group_list' style='height:400px;width:500px;'></div>
					</div>

					<script type='text/javascript'>
						var data = [{
							  values: [".$sample_by_group_list."],
							  labels: [".$group_list."],
							  type: 'pie'
							}];
						var layout={title:'Number of Samples by group', xaxis:{title:'duplication level'}, yaxis:{title:'Read percent (%) & GC ratio'}};
						Plotly.newPlot('".$repository."plot_sample_by_group_list', data, layout);
					</script>
				</div>
			</div>
		";

		$project_id=$repository."project_id";

		$project_link='
			<li class="nav-item mbr-fonts-style">
				<a class="nav-link active display-7" role="tab" data-toggle="tab" href="#tab'.$project_id.'">
					ALL Groups
				</a>
			</li>
		';

		$project_content='
			<div id="tab'.$project_id.'" class="tab-pane in active mbr-table" role="tabpanel">
				<div class="row tab-content-row">
					<div class="container pt-0 mt-0">
						<div class="media-container-row">
							'.$run_by_group_html.'
							'.$sample_by_group_html .'
						</div>
					</div>
				</div>
			</div>
		';


		foreach ($run_by_group_project_nb as $group=>$projects) {

			$project_id=str_replace(".","_",$repository.$group."project_id");

			$project_list="'".implode("','",array_keys($projects))."'";
			$run_by_group_project_list=implode(",",$projects);
			$sample_by_group_project_list=implode(",",$sample_by_group_project_nb[$group]);

			$run_by_group_html = "
				<div class='section_div'>
					<div id='".$project_id."run_by_group_list'>
						<div id='".$project_id."run_by_group_list_figure'>
							<div class='figure' id='".$project_id."plot_run_by_group_list' style='height:400px;width:500px;'></div>
						</div>

						<script type='text/javascript'>
							var data = [{
								  values: [".$run_by_group_project_list."],
								  labels: [".$project_list."],
								  type: 'pie'
								}];
							var layout={title:'Number of Analyses/Runs by projects', xaxis:{title:'duplication level'}, yaxis:{title:'Read percent (%) & GC ratio'}};
							Plotly.newPlot('".$project_id."plot_run_by_group_list', data, layout);
						</script>

					</div>
				</div>
			";

			$sample_by_group_html = "
				<div class='section_div'>
					<div id='".$project_id."sample_by_group_list'>
						<div id='".$project_id."sample_by_group_list_figure'>
							<div class='figure' id='".$project_id."plot_sample_by_group_list' style='height:400px;width:500px;'></div>
						</div>

						<script type='text/javascript'>
							var data = [{
								  values: [".$sample_by_group_project_list."],
								  labels: [".$project_list."],
								  type: 'pie'
								}];
							var layout={title:'Number of Samples by projects', xaxis:{title:'duplication level'}, yaxis:{title:'Read percent (%) & GC ratio'}};
							Plotly.newPlot('".$project_id."plot_sample_by_group_list', data, layout);
						</script>
					</div>
				</div>
			";


			# foreach group
			$project_link.='
				<li class="nav-item mbr-fonts-style">
					<a class="nav-link  display-7" role="tab" data-toggle="tab" href="#tab'.$project_id.'">
						'.$group.'
					</a>
				</li>
			';

			$project_content.='
				<div id="tab'.$project_id.'" class="tab-pane in mbr-table" role="tabpanel">
					<div class="row tab-content-row">

						<div class="container pt-0 mt-0">
				 		   	<div class="media-container-row">

								'.$run_by_group_html.'

								'.$sample_by_group_html .'

							</div>
						</div>


					</div>
				</div>
			';

		};

		$links_section.='

			<h3 class="mbr-section-subtitle mbr-fonts-style align-center mbr-light display-2">
			 '.$repository.'
			</h3>

			<div class="container-fluid col-md-12">
				<div class="row tabcont">
					<ul class="nav nav-tabs pt-0 mt-0" role="tablist">
						'.$project_link.'
					</ul>
				</div>
			</div>

			<div class="container">
				<div class="row px-1">
					<div class="tab-content">
						'.$project_content.'
					</div>
				</div>
			</div>

		';

	};

	$links_section.='
			</div>
		</section>
	';

 	echo $links_section;



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
