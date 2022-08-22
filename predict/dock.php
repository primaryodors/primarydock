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
$max_simultaneous_docks = 4;	// If running this script as a cron, we recommend setting this to half the number of physical cores.

$dock_results = [];
$json_file = "predict/dock_results.json";

$odors = json_decode(file_get_contents("data/odorant.json"), true);


if (file_exists($json_file))
{
	$dock_results = json_decode(file_get_contents($json_file), true);
}

foreach (@$argv as $a)
{
	$a = explode('=',$a,2);
	$_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}

if (@$_REQUEST['next'])
{	
	$cmd = "ps -ef | grep ':[0-9][0-9] bin/primarydock' | grep -v grep";
	exec($cmd, $results);
	if (!@$_REQUEST['force'] && trim(@$results[$max_simultaneous_docks-1])) die("Already running.\n".print_r($results, 1));
	$skip = count($results);
	
	$protid = @$_REQUEST['prot'] ?: false;
	$ligname = @$_REQUEST['lig'] ?: false;
	
	foreach ($odors as $o)
	{
		$full_name = str_replace(" ", "_", $o['full_name']);
		if ($ligname && $ligname != $full_name) continue;
		
		if (@$o['activity']) foreach ($o['activity'] as $ref => $acv)
		{
			foreach ($acv as $rcpid => $data)
			{
				if (!isset($dock_results[$rcpid][$full_name]))
				{
					if ($protid && $protid != $rcpid) continue;
					
					if ($skip)
					{
						$skip--;
						continue;
					}
					
					$protid = $rcpid;
					$ligname = $full_name;
					
					if (!file_exists("sdf/$ligname.sdf"))
					{
						$f = fopen("sdf/$ligname.sdf", "wb");
						if (!$f) die("Unable to create sdf/$ligname.sdf, please ensure write access.\n");
						$sdfdat = file_get_contents("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{$o['smiles']}/SDF?record_type=3d");
						fwrite($f, $sdfdat);
						fclose($f);
					}
					
					goto found_next_pair;
				}
			}
		}
	}
	die("All done!");
	
	found_next_pair:
	;
}
else
{
	$protid = @$_REQUEST['prot'] ?: "OR1A1";
	
	$ligname = @$_REQUEST['lig'] ?: "geraniol";
}

echo "Beginning dock of $ligname in $protid...\n\n";
$fam = family_from_protid($protid);

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

# NODEPDB 1 pdbs/$fam/$protid.active.pdb

# Activation Matrix.
# ACVMX region N.x N.y N.z C.x C.y C.z
ACVMX TMR1 -1.432411 -2.746519 2.509232 2.363341 -3.109768 -3.402768
ACVMX TMR2 4.639589 0.024734 -2.108768 -1.090161 -2.695269 1.679482
ACVMX TMR3 0.167589 -2.024517 3.395732 4.231589 1.114231 -4.413017
ACVMX TMR4 1.722089 0.166982 -1.611769 -2.093911 1.330481 2.336232
ACVMX TMR5 -2.133661 3.529982 1.032483 -3.523911 3.119732 2.745732
ACVMX TMR6 -2.557410 1.492481 -5.381518 -2.909160 3.353733 2.692482
ACVMX TMR7 -0.405411 0.009982 2.163982 3.021839 -3.566268 -1.637518

ACVNODE 1

SIZE 6.0 7.5 5.5

EXCL 1 $cyt1end		# Head, TMR1, and CYT1.
EXCL $exr2start $exr2end	# EXR2 between TMR4 and TMR5.

POSE 10
ITER 30

DIFF
ELIM 100

OUT output/$protid-$ligname.pred.dock



heredoc;

if (!file_exists("tmp")) mkdir("tmp");
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
	$average[$node?"Active":"Inactive"] = $value / (@$count[$node] ?: 1);
}

if (min($average) > 0)
{
	$prediction = "Non-Agonist";
}
else
{
	if ($average["Active"] <= 1.333 * $average["Inactive"])
		$prediction = "Agonist";
	else if ($average["Inactive"] > 1.333 * $average["Active"])
		$prediction = "Inverse Agonist";
	else $prediction = "Non-Agonist";
}

$average["Prediction"] = $prediction;

echo "Predicted ligand activity: $prediction\n";

if (@$_REQUEST['saved'])
	print_r($average);
else
{
	$dock_results[$protid][$ligname] = $average;
	
	$f = fopen($json_file, "wb");
	if (!$f) die("File write FAILED. Make sure have access to write $json_file.");
	fwrite($f, json_encode_pretty($dock_results));
	fclose($f);
}

















