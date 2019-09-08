<?php

################
# Task Spooler #
################

$TS_SOCKET="STARK";
$TS_SAVELIST="/tmp/STARKTSLIST";
$TS_ENV=" export TS_SOCKET=$TS_SOCKET && export TS_SAVELIST=$TS_SAVELIST ";
$TS_BIN="ts";
$TS="$TS_ENV && $TS_BIN ";


#########
# INPUT #
#########

// INPUT JSON from GET POST or RAW POST
$analysis_json=$_GET["analysis"].$_POST["analysis"].file_get_contents("php://input");

// INPUT JSON from GET POST or RAW POST for DOCKER RUN
$analysis_json_docker=$_GET["docker"].$_POST["docker"].file_get_contents("php://input");

// INIT
$jsonfile="";
$jsonfile_docker="";


// Create JSON file
if (is_file($analysis_json)) {
	#echo "JSON file: $analysis_json<br>";
	$jsonfile=$analysis_jsonfile;
} else if($analysis_json != "") {
	#echo "JSON string: $analysis_json<br>";
	$STARKAnalysis_dir="/tmp/STARKAnalysis.".rand();
	mkdir($STARKAnalysis_dir);
	$jsonfile="$STARKAnalysis_dir/analysis.json";
	$file = fopen($jsonfile, "w"); fwrite($file , $analysis_json); fclose($file);
	if (!is_file($jsonfile)) {
		echo "JSON file write failed";
		exit(0);
	};
};


// Create JSON file
if (is_file($analysis_json_docker)) {
	#echo "JSON file: $analysis_json<br>";
	$jsonfile_docker=$analysis_jsonfile;
} else if($analysis_json_docker != "") {
	#echo "JSON string: $analysis_json<br>";
	$STARKAnalysis_dir="/tmp/STARKAnalysis.".rand();
	mkdir($STARKAnalysis_dir);
	$jsonfile_docker="$STARKAnalysis_dir/analysis.json";
	$file = fopen($jsonfile, "w"); fwrite($file , $analysis_json_docker); fclose($file);
	if (!is_file($jsonfile_docker)) {
		echo "JSON file write failed";
		exit(0);
	};
};


print_r(getenv());
exit(0)

##########
# ACTION #
##########

$out = '';

// Operation on ID
if(isset($_GET['id']) & isset($_GET['op'])) {
	$id = $_GET['id'];
	$op = $_GET['op'];
	if($op == 'state') {
		$out = exec($TS.' -s '.$id);
	} else if( $op == 'log') {
		$out = exec($TS.' -o '.$id);
		$out=file_get_contents($out);
	} else if($op == 'info') {
		$out = exec($TS.' -i '.$id);
	}
// QUEUE list
} else if(isset($_GET['queue'])) {
	exec($TS.' -l', $out);
	#print_r($out);
	$out = implode(' <br />', $out);
// COMMAND specific
} else if(isset($_GET['command'])) {
	$command = $_GET['command'];
	#$command = base64_decode($command);
	$id = $TS.' '.$command;
	$out = exec($id);
// COMMAND STARK
} else if($jsonfile != "") {
	$command = 'bash -c "STARK --analysis='.$jsonfile.'"';
	#echo "STARK Command $command";
	#$command = base64_decode($command);
	$id = $TS.' '.$command;
	$out = exec($id);
// COMMAND STARK for docker
} else if($jsonfile != "") {
	$command = 'bash -c "docker run stark --analysis='.$jsonfile.'"';
	#echo "STARK Command $command";
	#$command = base64_decode($command);
	$id = $TS.' '.$command;
	$out = exec($id);
}


##########
# OUTPUT #
##########


if ( $_GET['format'] == "html" ) {
	echo "<H1>STARK Service</H1>";
	echo "<pre>";
	print($out);
	echo "</pre>";
} else {
	print($out);


}
