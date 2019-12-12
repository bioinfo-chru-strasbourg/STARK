<?php


	### VARIABLES
	###############

	$DEBUG=1;



	### FUNCTIONS
	###############




	# ENV VARIABLES
	#################

	# SERVER HOSTNAME/IP
	$DOCKER_STARK_SERVER_NAME=$_SERVER["SERVER_NAME"];
	if ($DOCKER_STARK_SERVER_NAME=="") { $DOCKER_STARK_SERVER_NAME="localhost"; };

	# DOCKER_STARK_SERVICE_PORT_PATTERN
	$DOCKER_STARK_SERVICE_PORT_PATTERN=$_ENV["DOCKER_STARK_SERVICE_PORT_PATTERN"];
	if ($DOCKER_STARK_SERVICE_PORT_PATTERN=="") { $DOCKER_STARK_SERVICE_PORT_PATTERN="42"; };

	# SERVER JARVIS PORT
	$DOCKER_STARK_SERVICE_JARVIS_PORT=$_ENV["DOCKER_STARK_SERVICE_JARVIS_PORT"];
	if ($DOCKER_STARK_SERVICE_JARVIS_PORT=="") { $DOCKER_STARK_SERVICE_JARVIS_PORT="92"; };

	# SERVER DATA CONTAINER NAME
	$DOCKER_STARK_SERVICE_DATA_CONTAINER_NAME=$_ENV["DOCKER_STARK_SERVICE_DATA_CONTAINER_NAME"];
	if ($DOCKER_STARK_SERVICE_DATA_CONTAINER_NAME=="") { $DOCKER_STARK_SERVICE_DATA_CONTAINER_NAME="stark-service-data"; };

	# DOCKER STARK SERVICE DATA PORT INNER
	$DOCKER_STARK_SERVICE_DATA_PORT=$_ENV["DOCKER_STARK_SERVICE_DATA_PORT"];
	if ($DOCKER_STARK_SERVICE_DATA_PORT=="") { $DOCKER_STARK_SERVICE_DATA_PORT="02"; };

	# DOCKER STARK SERVICE DATA PORT INNER
	$DOCKER_STARK_SERVICE_DATA_PORT_INNER=$_ENV["DOCKER_STARK_SERVICE_DATA_PORT_INNER"];
	if ($DOCKER_STARK_SERVICE_DATA_PORT_INNER=="") { $DOCKER_STARK_SERVICE_DATA_PORT_INNER="5000"; };

	# SERVER JARVIS CONTAINER NAME
	$DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME=$_ENV["DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME"];
	if ($DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME=="") { $DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME="stark-service-jarvis"; };

	# SERVER DATA PUBLIC DIR
	$DOCKER_STARK_SERVICE_DATA_PUBLIC_DIR=$_ENV["DOCKER_STARK_SERVICE_DATA_PUBLIC_DIR"];
	if ($DOCKER_STARK_SERVICE_DATA_PUBLIC_DIR=="") { $DOCKER_STARK_SERVICE_DATA_PUBLIC_DIR="/static/data/public"; };

	# DOCKER STARK INNER FOLDER OUTPUT LOG
	$DOCKER_STARK_SERVICE_DATA_INNER_FOLDER_DATA=$_ENV["DOCKER_STARK_SERVICE_DATA_INNER_FOLDER_DATA"];
	if ($DOCKER_STARK_SERVICE_DATA_INNER_FOLDER_DATA=="") { $DOCKER_STARK_SERVICE_DATA_INNER_FOLDER_DATA="/STARK/output/log"; };

	# DOCKER STARK INNER FOLDER OUTPUT LOG
	$DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_JARVIS=$_ENV["DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_JARVIS"];
	if ($DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_JARVIS=="") { $DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_JARVIS="/STARK/output/log"; };

	# DOCKER STARK INNER FOLDER OUTPUT LOG
	$DOCKER_STARK_SERVICE_LAUNCHER_FOLDER_LOG=$_ENV["DOCKER_STARK_SERVICE_LAUNCHER_FOLDER_LOG"];
	if ($DOCKER_STARK_SERVICE_LAUNCHER_FOLDER_LOG=="") { $DOCKER_STARK_SERVICE_LAUNCHER_FOLDER_LOG="analyses/stark-services/launcher"; };

	# DOCKER STARK INNER FOLDER OUTPUT LOG
	$DOCKER_STARK_SERVICE_LAUNCHER_CONTAINER_NAME=$_ENV["DOCKER_STARK_SERVICE_LAUNCHER_CONTAINER_NAME"];
	if ($DOCKER_STARK_SERVICE_LAUNCHER_CONTAINER_NAME=="") { $DOCKER_STARK_SERVICE_LAUNCHER_CONTAINER_NAME="STARK-services-launcher"; };

	# DOCKER STARK INNER FOLDER OUTPUT LOG
	$DOCKER_STARK_SERVICE_LAUNCHER_PORT_INNER=$_ENV["DOCKER_STARK_SERVICE_LAUNCHER_PORT_INNER"];
	if ($DOCKER_STARK_SERVICE_LAUNCHER_PORT_INNER=="") { $DOCKER_STARK_SERVICE_LAUNCHER_PORT_INNER="05"; };




	$DATA_URL_EXTERNAL=preg_replace("/\/$/", "", str_replace("//","/","$DOCKER_STARK_SERVER_NAME:$DOCKER_STARK_SERVICE_PORT_PATTERN$DOCKER_STARK_SERVICE_DATA_PORT/$DOCKER_STARK_SERVICE_DATA_PUBLIC_DIR/$DOCKER_STARK_SERVICE_LAUNCHER_FOLDER_LOG"));


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


	### CONTENT
	#############

	$URL_LAUNCHER="http://$DOCKER_STARK_SERVICE_LAUNCHER_CONTAINER_NAME:$DOCKER_STARK_SERVICE_LAUNCHER_PORT_INNER/queue";

	$TS_LIST=file($URL_LAUNCHER);

	// echo "<pre>";
	// print_r($TS_LIST);
	// echo "</pre>";

	### TABLE
	###########

	### THEAD




	foreach ($TS_LIST as $TASK_KEY => $ONE_TASK) {

		# SAMPLE NAME
		$ONE_TASK_split=explode( " ", preg_replace('/\s+/', ' ',($ONE_TASK)) ); # preg_replace('/\s+/', ' ',$row['message']);


		if ($TASK_KEY==0) {
			$thead='
				<tr class="table-heads">
					<th class="head-item mbr-fonts-style display-7">
						'.$ONE_TASK_split[0].'
					</th>
					<th class="head-item mbr-fonts-style display-7">
						Status
					</th>
					<th class="head-item mbr-fonts-style display-7">
						Time
					</th>
					<th class="head-item mbr-fonts-style display-7">
						Analysis/Run & Command ID
					</th>
				</tr>
				';
		} else {

			$ONE_TASK_ID=$ONE_TASK_split[0];
			$ONE_TASK_STATE=$ONE_TASK_split[1];
			$ONE_TASK_OUTPUT=$ONE_TASK_split[2];
			if ($ONE_TASK_STATE=="running" || $ONE_TASK_STATE=="queued" ) {
				$ONE_TASK_E_LEVEL="";
				$ONE_TASK_TIMES="";
				$ONE_TASK_COMMAND=$ONE_TASK_split[3];
			} else {
				$ONE_TASK_E_LEVEL=$ONE_TASK_split[3];
				$ONE_TASK_TIMES=$ONE_TASK_split[4];
				$ONE_TASK_COMMAND=$ONE_TASK_split[5];
			}

			# INFO
			$ONE_TASK_INFO_ARRAY=file("$URL_LAUNCHER?info=$ONE_TASK_ID");
			$ONE_TASK_INFO_CONTENT=implode("\n",$ONE_TASK_INFO_ARRAY);
			$ONE_TASK_INFO_CONTENT_HTML=implode("<br><br>",$ONE_TASK_INFO_ARRAY);

			if ($ONE_TASK_STATE=="running") {
				$ONE_TASK_TIME_RUNNING=str_replace('s','',explode(":",trim($ONE_TASK_INFO_ARRAY[4]))[1]);
			} else {
				$ONE_TASK_TIME_RUNNING="";
			}

			# LOG
			$ONE_TASK_LOG_ARRAY=file("$URL_LAUNCHER?log=$ONE_TASK_ID");
			$ONE_TASK_LOG_CONTENT=implode("",$ONE_TASK_LOG_ARRAY);

			$ONE_TASK_BASENAME=basename($ONE_TASK_OUTPUT);

			preg_match('#\[(.*?)\]#', $ONE_TASK_COMMAND, $match);
			$ONE_TASK_COMMAND_ID=$match[1];

			#$ONE_TASK_RUN=explode('.',$ONE_TASK_COMMAND_ID)[2];
			$ONE_TASK_RUN_TYPE=explode('.',$ONE_TASK_COMMAND_ID)[0];
			$ONE_TASK_RUN=explode('-',explode('.',$ONE_TASK_COMMAND_ID)[2])[3];


			if ($ONE_TASK_OUTPUT!="(file)") {
				$ONE_TASK_LOG_NAME="output";
				$ONE_TASK_LOG_URL="http://$DATA_URL_EXTERNAL/$ONE_TASK_BASENAME";
			} else {
				$ONE_TASK_LOG_NAME="";
				$ONE_TASK_LOG_URL="";
			};

			if ($ONE_TASK_E_LEVEL!="" && $ONE_TASK_E_LEVEL!="0" ) {
				$ONE_TASK_E_LEVEL_HTML="<br>Error '$ONE_TASK_E_LEVEL'";
			} else {
				$ONE_TASK_E_LEVEL_HTML="";
			};

			# Times
			$ONE_TASK_TIME="";
			if ($ONE_TASK_TIMES!="") {
				$ONE_TASK_TIME=str_replace(array("/",":"),array("h","m"),gmdate("H/i:s", explode("/",$ONE_TASK_TIMES)[0]))."s";
			} else if ($ONE_TASK_TIME_RUNNING!="") {
				$ONE_TASK_TIME=str_replace(array("/",":"),array("h","m"),gmdate("H/i:s", explode("/",$ONE_TASK_TIME_RUNNING)[0]))."s";
			}

			# state finished/green queued/black running/orange other/red
			#echo "<BR>$ONE_TASK_STATE && $ONE_TASK_E_LEVEL";

			if ($ONE_TASK_STATE=="finished" && $ONE_TASK_E_LEVEL==0) {
				$ONE_TASK_STATE_COLOR="green";
				$ONE_TASK_TIME_COLOR="";
			} else if ($ONE_TASK_STATE=="finished" && $ONE_TASK_E_LEVEL!=0) {
				$ONE_TASK_STATE_COLOR="red";
				$ONE_TASK_TIME_COLOR="red";
			} else if ($ONE_TASK_STATE=="queued") {
				$ONE_TASK_STATE_COLOR="silver";
				$ONE_TASK_TIME_COLOR="";
			} else if ($ONE_TASK_STATE=="running") {
				$ONE_TASK_STATE_COLOR="orange";
				$ONE_TASK_TIME_COLOR="orange";
			} else {
				$ONE_TASK_STATE_COLOR="red";
				$ONE_TASK_TIME_COLOR="";
			};

			# Last task (by ID)
			if ($ONE_TASK_RUN_TYPE=="STARK") {
				if (!isset($run_task[$ONE_TASK_RUN]["task"]["id"]) || $ONE_TASK_ID>$run_task[$ONE_TASK_RUN]["task"]["id"] ) {
					$run_task[$ONE_TASK_RUN]["task"]["status"]=$ONE_TASK_STATE;
					$run_task[$ONE_TASK_RUN]["task"]["color"]=$ONE_TASK_STATE_COLOR;
					$run_task[$ONE_TASK_RUN]["task"]["time"]=$ONE_TASK_TIME;
					$run_task[$ONE_TASK_RUN]["task"]["time_color"]=$ONE_TASK_TIME_COLOR;
					$run_task[$ONE_TASK_RUN]["task"]["error_level"]=$ONE_TASK_E_LEVEL;
					$run_task[$ONE_TASK_RUN]["task"]["id"]=$ONE_TASK_ID;
				}
			}

			$ONE_TASK_RUN_HEAD="[$ONE_TASK_RUN_TYPE] $ONE_TASK_RUN";

			$tbody=$tbody.'
				<tr class="table-heads">
					<td class="head-item mbr-fonts-style display-7">
						<small>'.$ONE_TASK_ID.'</small>
					</td>
					<td class="head-item mbr-fonts-style display-7">
						<small style="color:'.$ONE_TASK_STATE_COLOR.'">
							'.$ONE_TASK_STATE.'
							'.$ONE_TASK_E_LEVEL_HTML.'
							<br>
							<a href="'.$ONE_TASK_LOG_URL.'" download>'.$ONE_TASK_LOG_NAME.'</a>
						</small>
					</td>
					<td class="head-item mbr-fonts-style display-7">
						<small style="color:'.$ONE_TASK_TIME_COLOR.'">'.$ONE_TASK_TIME.'</small>
					</td>
					<td class="head-item mbr-fonts-style display-7 ">
						'.$ONE_TASK_RUN_HEAD.'
						<br>
						<a>
						<small onclick=\'if(document.getElementById("info_'.$TASK_KEY.'").style.display=="none"){document.getElementById("info_'.$TASK_KEY.'").style.display="block"}else{document.getElementById("info_'.$TASK_KEY.'").style.display="none"};\'>'.$ONE_TASK_COMMAND_ID.'</small>
						</a>
						<div id="info_'.$TASK_KEY.'" style="display:none">
							<small><br>'.$ONE_TASK_INFO_CONTENT_HTML.'</small>
						</div>
					</td>

				</tr>

				';


		}



	};

	### SECTION
	#############

	echo '

	<section class="section-table cid-ruiBoanIwc" id="METRICS">

	  <div class="container container-table">

		  <h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
			  Task Spooler
		  </h2>

		  <p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
			  Task Spooler queue up STARK analyses
		  </p>

		  <div class="div-wrapper" style="max-height:600px;">

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

	<br>

	';


if (1) {


	# INPUTS
	##########

	$runs_infos=array();


	# RUN folders
	$runs_inputs=glob("inputs/*/runs/*",GLOB_ONLYDIR|GLOB_BRACE);
	foreach ($runs_inputs as $runs_inputs_key=>$run_path) {
		$run_path_split=explode("/",$run_path);
		$run=$run_path_split[count($run_path_split)-1];
		$runs_infos[$run]["inputs"]["run_path"][$run_path]=$run_path;
	};


	# Runs specific files
	$runs_inputs_files=glob("inputs/*/*/*/{RTAComplete.txt,SampleSheet.csv}",GLOB_BRACE);
	foreach ($runs_inputs_files as $runs_inputs_key=>$run_file_path) {
		$run_path=dirname($run_file_path);
		$run_file_path_split=explode("/",$run_file_path);
		$run=$run_file_path_split[count($run_file_path_split)-2];
		$run_file=$run_file_path_split[count($run_file_path_split)-1];
		$runs_infos[$run]["inputs"]["run_files"][$run_file][$run_path]=$run_file_path;
		#echo $run_file;
		if ($run_file=="SampleSheet.csv") {
			#echo "<br>SampleSheet";
			$run_file_path_array=file($run_file_path);

			# find INFOS within log file
			$run_file_path_array=file($run_file_path);

			# Explode listener file infos
			#$listener_infos=array();
			$ini_file_content="";
			foreach ($run_file_path_array as $run_file_path_array_key=>$run_file_path_array_info) {

				if (substr(trim($run_file_path_array_info),0,1)!="[") {
					$ini_file_content.=rand(1,10000)."=\"".trim($run_file_path_array_info)."\"\n";
				} else {
					$ini_file_content.=trim($run_file_path_array_info)."\n";
				}

				# ADD
				#$listener_infos[$runs_listener_log_file_array_info_matches[1]]=$runs_listener_log_file_array_info_matches[2];
			};
			#echo "<pre>$ini_file_content</pre>";
			$ini_filename="/tmp/".rand(1,100000);
			$ini_file = fopen($ini_filename, "w") or die("Unable to open file!");
			fwrite($ini_file, $ini_file_content);
			fclose($ini_file);
			$ini_array = parse_ini_file($ini_filename, true);
			unlink($ini_file);
			$runs_infos[$run]["inputs"]["run_files"][$run_file][$run_path]=$ini_array;
			#echo "<pre>";
			#print_r($ini_array["Data"]);

			foreach ($ini_array["Data"] as $ini_array_key=>$ini_array_sample) {
				$ini_array_sample_split=explode(",",$ini_array_sample);
				#print_r($ini_array_sample_split);
				if ($ini_array_sample_split[0] != "Sample_ID") {
					if ($ini_array_sample_split[1]!="") {
						$runs_infos[$run]["inputs"]["run_files"][$run_file]["samples"][$ini_array_sample_split[1]]=$ini_array_sample_split[1];
						$runs_infos[$run]["inputs"]["run_files"][$run_file]["samples_infos"][$ini_array_sample_split[1]]=$ini_array_sample_split;
					} else if ($ini_array_sample_split[0]!="") {
						$runs_infos[$run]["inputs"]["run_files"][$run_file]["samples"][$ini_array_sample_split[0]]=$ini_array_sample_split[0];
						$runs_infos[$run]["inputs"]["run_files"][$run_file]["samples_infos"][$ini_array_sample_split[0]]=$ini_array_sample_split;
					}
				}
			}
			#echo "</pre>";
		}
	};

	#echo "<pre>"; print_r($runs_infos); echo "</pre>";

	# REPOSITORIES
	$runs_repositories=glob("repositories/*/*/*/*/*",GLOB_ONLYDIR);
	foreach ($runs_repositories as $runs_inputs_key=>$sample_path) {
		$sample_path_split=explode("/",$sample_path);
		$sample=basename($sample_path);
		#$run=$sample_path_split[count($sample_path_split)-2];
		$run=basename(dirname($sample_path));
		$run_path=dirname($sample_path);
		#echo "<br>run:$run - sample=$sample";
		$runs_infos[$run]["repositories"]["run_path"][$run_path]=$run_path;
		$runs_infos[$run]["repositories"]["samples"][$sample]=$sample;
	};


	# LISTENER LOG
	$LISTENER_LOG_PATTERN="{ID-*-NAME-*.log}";
	$runs_listener_log=glob("analyses/stark-services/listener/$LISTENER_LOG_PATTERN",GLOB_BRACE);
	foreach ($runs_listener_log as $runs_listener_log_key=>$runs_listener_log_file) {
		$runs_listener_log_file_split=explode("/",$runs_listener_log_file);
		$run_listener_log=$runs_listener_log_file_split[count($runs_listener_log_file_split)-1];

		# find INFOS within log file
		$runs_listener_log_file_array=file($runs_listener_log_file);

		# Explode listener file infos
		$listener_infos=array();
		foreach ($runs_listener_log_file_array as $runs_listener_log_file_array_key=>$runs_listener_log_file_array_info) {
			#echo "<br>$runs_listener_log_file_array_info";
			preg_match("/(.*): (.*)/i", $runs_listener_log_file_array_info,$runs_listener_log_file_array_info_matches);
			#print_r($runs_listener_log_file_array_info_matches);
			# ADD
			$listener_infos[$runs_listener_log_file_array_info_matches[1]]=$runs_listener_log_file_array_info_matches[2];
		};

		# RUN
		$run=$listener_infos["RUN"];

		# ADD
		$runs_infos[$run]["listener"][$run_listener_log]["log"]=$run_listener_log;
		$runs_infos[$run]["listener"][$run_listener_log]["infos"]=$listener_infos;

	};

	# LAUNCHER LOG
	#analyses/stark-services/igv/20191204-080425.432691.json
	#$launcher_log=glob ( "igvjs/static/data/public/stark-services/launcher/*", GLOB_ONLYDIR );
	$LAUNCHER_LOG_PATTERN="{*.json,*.log,*.err,*.info,ts-out.*}";
	$runs_launcher_log=glob("analyses/stark-services/launcher/$LAUNCHER_LOG_PATTERN",GLOB_BRACE);
	foreach ($runs_launcher_log as $runs_launcher_log_key=>$runs_launcher_log_file) {
		#echo "<br><br>";
		#print_r($runs_launcher_log_file);
		$runs_listener_log_file_split=explode("/",$runs_launcher_log_file);
		$run_listener_log=$runs_listener_log_file_split[count($runs_listener_log_file_split)-1];
		$ext = pathinfo($runs_launcher_log_file, PATHINFO_EXTENSION);
		$runs_infos[$run]["launcher"][$runs_launcher_log_file][$ext]=$runs_launcher_log_file;
		if ($ext=="json") {
			#echo "<pre>";
			#print_r(file($runs_launcher_log_file));
			$runs_launcher_log_file_json = json_decode(implode("",file($runs_launcher_log_file)));
			$run_path=$runs_launcher_log_file_json->run;
			$run=basename($run_path);
			if ($run!="") {
				$runs_infos[$run]["launcher"][$runs_launcher_log_file]["run"]=$run;
			};
			#echo "</pre>";
		} else {
			# TODO - inegrate log files !!!
			#$runs_infos[$run]["launcher"][$runs_launcher_log_file]["run"]=$run;
		}
	};




	# TABLE
	# RUN / Sequencing / Analysis / Results


	$thead='
		<tr class="table-heads">
			<th class="head-item mbr-fonts-style display-7">
				Run
			</th>
			<th class="head-item mbr-fonts-style display-7">
				Sequencing
			</th>
			<th class="head-item mbr-fonts-style display-7">
				Analysis
			</th>
			<th class="head-item mbr-fonts-style display-7">
				Repository
			</th>
		</tr>
		';

	$tbody="";

	foreach ($runs_infos as $run=>$run_infos) {

		if ($DEBUG && 0) {
			echo "<br><br>";
			echo "run: $run";
		};

		$run_progress[$run]=array();

		$run_status="unknown";

		# Sequencing
		############

		$sequencing_status_code=isset($run_infos["inputs"])+isset($run_infos["inputs"]["run_files"]["RTAComplete.txt"]);

		switch ($sequencing_status_code) {
	    case 0:
			$sequencing_status="unavailable";
			$sequencing_color="gray";
			$sequencing_message="No run sequencing folder";
			break;
	    case 1:
	        $sequencing_status="in progress";
			$sequencing_color="orange";
			$sequencing_message="Run sequencing seems to be in progress";
	        break;
	    case 2:
	        #$sequencing_status="complete";
			if (!isset($run_infos["inputs"]["run_files"]["SampleSheet.csv"])) {
				$sequencing_status="ready";
				$sequencing_color="dodgerblue";
				$sequencing_message="Ready for analysis (SampleSheet missing)";
			} else {
				$sequencing_status="complete";
				$sequencing_color="green";
				$sequencing_message="Complete for analysis ";
			}
	        break;
		default:
	        $sequencing_status="unknown";
			$sequencing_color="gray";
			$sequencing_message="No information";
	        break;
		}

		$sequencing_message_plus="";
		if (count($runs_infos[$run]["inputs"]["run_files"]["SampleSheet.csv"]["samples"])) {
			$sequencing_message_plus="<br>".implode(" ",array_keys($runs_infos[$run]["inputs"]["run_files"]["SampleSheet.csv"]["samples"]));
		}

		$run_progress[$run]["sequencing"]["status"]=$sequencing_status;
		$run_progress[$run]["sequencing"]["color"]=$sequencing_color;
		$run_progress[$run]["sequencing"]["message"]=$sequencing_message;
		$run_progress[$run]["sequencing"]["message_plus"]=$sequencing_message_plus;

		if ($sequencing_status!="unavailable") {
			$run_status="Sequencing ".$sequencing_status;
			$run_color=$sequencing_color;
		}



		# Analysis
		#############

		// echo "<pre>";
		// print_r($run_task[$run]);
		// echo "</pre>";

		$analysis_status_code=isset($run_infos["listener"])+isset($run_infos["launcher"]);
		switch ($analysis_status_code) {
		case 0:
			$analysis_status="unavailable";
			$analysis_color="gray";
			$analysis_message="No analysis information available";
			break;
		case 1:
			$analysis_status="detected";
			$analysis_color="orange";
			$analysis_message="Run sequencing ready detected";
			break;
		case 2:
			$analysis_status="launched";
			$analysis_color="orange";
			$analysis_message="Run analysis launched";
			#$run_task[$ONE_TASK_RUN]["task"]["time"]=$ONE_TASK_STATE;
			#echo "<br><br>"; print_r($run_task[$ONE_TASK_RUN]);
			if (isset($run_task[$run]["task"])) {
				$analysis_status=$run_task[$run]["task"]["status"];
				$analysis_color=$run_task[$run]["task"]["color"];
				$analysis_message="";
				if ($run_task[$run]["task"]["error_level"]!=0) {
					$analysis_color="red";
					$run_status="error";
					$analysis_message.="<span style='color:red'>Error '".$run_task[$run]["task"]["error_level"]."' </span><br>";
				}
				$analysis_message.="<span style='color:".$run_task[$run]["task"]["time_color"]."'>".$run_task[$run]["task"]["time"]."</span>";

			} else {
				# TODO - Check on log file !!!
				$analysis_status="no information";
				$analysis_color="red";
				$analysis_message="Analysis launched but no Task Spooler information";
			};
			break;
		default:
			$analysis_status="unknown";
			$analysis_color="gray";
			$analysis_message="No information";
			break;
		}

		#echo "<br>$analysis_status";

		$run_progress[$run]["analysis"]["status"]=$analysis_status;
		$run_progress[$run]["analysis"]["color"]=$analysis_color;
		$run_progress[$run]["analysis"]["message"]=$analysis_message;
		$run_progress[$run]["analysis"]["message_plus"]="";

		if ($analysis_status!="unavailable") {
			$run_status="Analysis ".$analysis_status;
			$run_color=$analysis_color;
		}


		# Repository
		#############

		$repository_status_code=count($run_infos["repositories"]["run_path"]);
		if ($repository_status_code==0) {
			$repository_status="unavailable";
			$repository_color="gray";
			$repository_message="No results available";
			$repository_message_plus="";

			# change analysis status...
			// if ($run_progress[$run]["analysis"]["status"]=="finished") {
			// 	$run_progress[$run]["analysis"]["status"]="finished - error";
			// 	$run_progress[$run]["analysis"]["color"]="red";
			// 	$run_progress[$run]["analysis"]["message"].="<br>(Uncomplete and NOT available in repository)<br>Probably Errors in analysis (Please check log files)";
			// }

		} else {
			$repository_status="available";
			$repository_color="green";
			$repository_message="Run results available in folders";
			$repository_message_plus="";

			# change analysis status...
			// if ($run_progress[$run]["analysis"]["status"]=="no information") {
			// 	$run_progress[$run]["analysis"]["color"]="green";
			// 	$run_progress[$run]["analysis"]["message"].="<br>(Complete and available in repository)";
			// }


			# Check number of samples !
			if (isset($runs_infos[$run]["inputs"]["run_files"]["SampleSheet.csv"]) && isset($runs_infos[$run]["inputs"]["run_files"]["SampleSheet.csv"]["samples"])) {
				#$runs_infos[$run]["inputs"]["run_files"]["SampleSheet.csv"]["samples"]["truc"]="truc";
				$result1=array_diff($runs_infos[$run]["repositories"]["samples"],$runs_infos[$run]["inputs"]["run_files"]["SampleSheet.csv"]["samples"]);
				$result_missing=array_diff($runs_infos[$run]["inputs"]["run_files"]["SampleSheet.csv"]["samples"],$runs_infos[$run]["repositories"]["samples"]);
				$result=array_merge($result1,$result_missing);
				if (count($result_missing)) {
					if ($run_progress[$run]["analysis"]["color"]=="green") {
						$missing_samples_color="red";
					} else {
						$missing_samples_color=$run_task[$run]["task"]["color"];
					}
					$repository_message_plus.="<br><span style='color:$missing_samples_color'>Missing Samples (".implode(" ",$result_missing).")</span>";
				}
			}


			foreach ($run_infos["repositories"]["run_path"] as $repository_path=>$repository_path_infos) {
				# http://localhost:42002/index.reports.php?PATH=repositories/Archive/SOMATIC/HEMATOLOGY/RUN_TEST_TAG/
				#$REPORTS_RUN_LINK="<a target='METRICS' href='index.reports.php?PATH=$run_path'>REPORT</a>";
				$repository_path_split=explode("/",$repository_path);
				$repository_message_plus.="<br><a target='REPORT' href='index.reports.php?PATH=$repository_path'>".$repository_path_split[1]."/".$repository_path_split[2]."/".$repository_path_split[3]."</a>";
			}



			// $run_progress[$run]["analysis"]["message_plus"]="";
		};

		$run_progress[$run]["repository"]["status"]=$repository_status;
		$run_progress[$run]["repository"]["color"]=$repository_color;
		$run_progress[$run]["repository"]["message"]=$repository_message;
		$run_progress[$run]["repository"]["message_plus"]=$repository_message_plus;

		if ($repository_status!="unavailable") {
			$run_status="Results ".$repository_status;
			$run_color=$repository_color;
		}

		# <b>'.$run_status.'</b>

		$tbody=$tbody.'
			<tr class="table-heads">
				<td class="head-item mbr-fonts-style display-7">
					<small><b>'.$run.'</b></small>
					<br>
					<small style="color:'.$run_color.'"></small>


				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small style="color:'.$run_progress[$run]["sequencing"]["color"].'"><b>'.$run_progress[$run]["sequencing"]["status"].'</b></small>
					<br>
					<a>
					<small onclick=\'if(document.getElementById("info_sequencing_'.$run.'").style.display=="none"){document.getElementById("info_sequencing_'.$run.'").style.display="block"}else{document.getElementById("info_sequencing_'.$run.'").style.display="none"};\'>'.$run_progress[$run]["sequencing"]["message"].'</small>
					</a>
					<div id="info_sequencing_'.$run.'" style="display:none">
						<small>'.$run_progress[$run]["sequencing"]["message_plus"].'</small>
					</div>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small style="color:'.$run_progress[$run]["analysis"]["color"].'"><b>'.$run_progress[$run]["analysis"]["status"].'</b></small>
					<br>
					<a>
					<small onclick=\'if(document.getElementById("info_analysis_'.$run.'").style.display=="none"){document.getElementById("info_analysis_'.$run.'").style.display="block"}else{document.getElementById("info_analysis_'.$run.'").style.display="none"};\'>'.$run_progress[$run]["analysis"]["message"].'</small>
					</a>
					<div id="info_analysis_'.$run.'" style="display:none">
						<small>'.$run_progress[$run]["analysis"]["message_plus"].'</small>
					</div>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small style="color:'.$run_progress[$run]["repository"]["color"].'"><b>'.$run_progress[$run]["repository"]["status"].'</b></small>
					<br>
					<a>
					<small onclick=\'if(document.getElementById("info_repository_'.$run.'").style.display=="none"){document.getElementById("info_repository_'.$run.'").style.display="block"}else{document.getElementById("info_repository_'.$run.'").style.display="none"};\'>'.$run_progress[$run]["repository"]["message"].'</small>
					</a>
					<div id="info_repository_'.$run.'" style="display:none">
						<small>'.$run_progress[$run]["repository"]["message_plus"].'</small>
					</div>
				</td>

			</tr>

			';

	}


	### SECTION
	#############

	echo '

	<section class="section-table cid-ruiBoanIwc" id="METRICS">

	  <div class="container container-table">

		  <h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
			  Analyses/Runs Progress
		  </h2>

		  <p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
			  Run process activity (sequencing, analysis, repository)
		  </p>

		  <div class="">

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

	<br><br><br>

	';





	if ($DEBUG && 0) {
		echo "<br><br><br>RUNS_INFOS<br>";
		echo "<pre>";
		print_r($run_progress);
		print_r($runs_infos);
		echo "</pre>";
	};


}





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
