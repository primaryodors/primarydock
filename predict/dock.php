<?php

// dock.php
//
// Performs a dock of an odorant in a receptor using both the inactive and active PDB files.
//
// Example call syntax:
// php -f predict/dock.php prot=OR1A1 lig=geraniol
//

require("protutils.php");

$dock_retries = 5;

$dock_results = [];

foreach (@$argv as $a)
{
	$a = explode('=',$a,2);
	$_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}

$protid = @$_REQUEST['prot'] ?: "OR1A1";
$fam = family_from_protid($protid);

$ligname = @$_REQUEST['lig'] ?: "geraniol";

$cenr3_33 = resno_from_bw($protid, "3.33");
$cenr3_36 = $cenr3_33 + 3;
$cenr5_43 = resno_from_bw($protid, "5.43");
$cenr5_47 = $cenr5_43 + 4;
$cenr6_40 = resno_from_bw($protid, "6.40");
$cenr6_47 = $cenr6_40 + 7;

$tmr2start = $prots[$protid]['region']['TMR2']['start'];
$cyt1end = $tmr2start - 1;
$tmr4end = $prots[$protid]['region']['TMR4']['end'];
$tmr5start = $prots[$protid]['region']['TMR5']['start'];
$exr2start = $tmr4end + 1;
$exr2end = $tmr5start - 1;


$configf = <<<heredoc

PROT pdbs/$fam/$protid.rotated.pdb

LIG sdf/$ligname.sdf

CEN RES $cenr3_33 $cenr3_36 $cenr5_43 $cenr5_47 $cenr6_40 $cenr6_47
PATH 1 REL 0 0 0

NODEPDB 1 pdbs/$fam/$protid.active.pdb

SIZE 6.0 7.5 5.5

EXCL 1 $cyt1end		# Head, TMR1, and CYT1.
EXCL $exr2start $exr2end	# EXR2 between TMR4 and TMR5.

POSE 10
ITER 30

DIFF
ELIM 100

OUT output/$protid-$ligname.pred.dock

heredoc;

$f = fopen("tmp/prediction.config", "wb");
if (!$f) die("File write error. Check tmp folder is write enabled.\n");

fwrite($f, $configf);
fclose($f);

$outlines = [];
if (@$_REQUEST['saved'])
{
	$outlines = explode("\n", file_get_contents($_REQUEST['saved']));
}
else
{
	for ($try = 0; $try < $dock_retries; $try++)
	{
		set_time_limit(300);
		exec("bin/primarydock tmp/prediction.config", $outlines);
		if (count($outlines) >= 100) break;
	}
}
// echo implode("\n", $outlines) . "\n\n";

if (count($outlines) < 100) die("Docking FAILED.\n");

$benerg = [];
$pose = false;
$node = -1;

foreach ($outlines as $ln)
{
	if (substr($ln, 0, 6) == "Pose: ") $pose = intval(explode(" ", $ln)[1]);
	if (substr($ln, 0, 6) == "Node: ") $node = intval(explode(" ", $ln)[1]);
	if (trim($ln) == "TER")
	{
		$pose = false;
		$node = -1;
	}
	
	if ($pose && $node>=0 && substr($ln, 0, 7) == "Total: ") $benerg[$pose][$node] = floatval(explode(" ", $ln)[1]);
}

// print_r($benerg);

$sum = [];
$count = [];
foreach ($benerg as $pose => $data)
{
	foreach ($data as $node => $value)
	{
		if (!isset($sum[$node]  )) $sum[$node] = 0.0;
		if (!isset($count[$node])) $count[$node] = 0;
		
		$sum[$node] += $value;
		$count[$node]++;
	}
}

$average = [];

foreach ($sum as $node => $value)
{
	$average[$node] = $value / (@$count[$node] ?: 1);
}

if (@$_REQUEST['saved'])
	print_r($average);
else
{
	$dock_results[$protid][$ligname] = $average;
	
	$f = fopen("predict/dock_results.json", "wb");
	if (!f) die("File write FAILED.");
	fwrite($f, json_encode_pretty($dock_results));
}






