<?php


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
	';

	// <style>
	// 	.table-wrapper {
	// 		overflow: auto;
	// 		max-height:800px;
	// 		}
	// </style>


	echo '
		<section class="header1 cid-ru7OEConn1" id="header16-1k">
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
		    <div class="container">
				<h3 class="mbr-section-subtitle mbr-fonts-style align-center mbr-light display-3">
					Reports
			  	</h3>
		    </div>
		</section>

	';


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


	$tbody="";
	#$reports=glob ( "repositories/*/*/*/*/*/*stark.report.html" );
	$reports=glob ( "repositories/*/*/*/*/*/*stark.report.html" );
	#$reports=glob ( "repositories/*" );
	#print_r($reports);
	foreach ($reports as $key => $report_html) {
		$report_html_split=explode ( "/" , $report_html );
		#echo "<br><a href='$report_html'>".$report_html_split[0].$report_html_split[5]."</a>";
		$root=$report_html_split[0];
		$repository=$report_html_split[1];
		$group=$report_html_split[2];
		$project=$report_html_split[3];
		$run=$report_html_split[4];
		$sample=$report_html_split[5];
		$report_html_file=$report_html_split[6];
		$report_html_file_split=explode ( "." , $report_html_file );
		$report_html_id=$report_html_file_split[1];
		$report_tsv=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run.'/'.$sample.'/'.$sample.'.final.tsv';
		$report_vcf_gz=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run.'/'.$sample.'/'.$sample.'.final.vcf.gz';
		$report_bed=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run.'/'.$sample.'/'.$sample.'.bed';
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



# DER_BIN/csv_to_html.awk
#-v tag_row_head_plus="class='table-heads'"
#-v tag_col_head_plus=" class=\"head-item mbr-fonts-style display-7\""
#-v tag_row_body_plus=""
#-v tag_col_body_plus=" class=\"body-item mbr-fonts-style display-7\"")

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
	'

?>
