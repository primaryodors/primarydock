<?php

// activate.php
//
// Activate a protein by bending its backbone at 6.40:N-CA until the optimum 3.39-6.40 bond occurrs, then bending 5.47:N-CA to match.
//
// Example call syntax:
// php -f predict/activate.php prot=OR1A1
//

chdir(__DIR__);
$prots = json_decode(file_get_contents("../data/receptor.json"));

foreach (@$argv as $a)
{
	$a = explode('=',$a,2);
	$_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}

$protid = @$_REQUEST['prot'] ?: "OR1A1";


