<?php

// method_exr_capt.php
//
// Performs a dock of an odorant in a receptor using both the inactive and active PDB files.
//
// Example call syntax:
// php -f predict/method_exr_capt.php prot=TAAR5 lig=trimethylamine
//

// Includes
require("protutils.php");
require("odorutils.php");
require("matrices.php");

// Configurable variables
$dock_retries = 5;
$max_simultaneous_docks = 4;	// If running this script as a cron, we recommend setting this to half the number of physical cores.


// Load data
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
	$protid = @$_REQUEST['prot'] ?: "TAAR5";
	
	$ligname = @$_REQUEST['lig'] ?: "trimethylamine";
}

echo "Beginning dock of $ligname in $protid...\n\n";
$fam = family_from_protid($protid);


$tmr2start = $prots[$protid]['region']['TMR2']['start'];
$cyt1end = $tmr2start - 1;
$tmr4end = $prots[$protid]['region']['TMR4']['end'];
$tmr5start = $prots[$protid]['region']['TMR5']['start'];
$exr2start = $tmr4end + 1;
$exr2end = $tmr5start - 1;



$configf = <<<heredoc

PROT pdbs/$fam/$protid.rotated.pdb

LIG sdf/$ligname.sdf

CEN RES 184 185
PATH 1 RES 185 191 193
PATH 2 RES 94 195
PATH 3 RES 114 265

SIZE 7.0 7.5 7.0

POSE 10
ITER 50

# DIFF
ELIM 30

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

if (@$_REQUEST['echo']) echo implode("\n", $outlines) . "\n\n";

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
	$average["Node $node"] = round($value / (@$count[$node] ?: 1), 3);
}

$prediction = "(insufficient data)";

$actual = best_empirical_pair($protid, $ligname);
$actual = ($actual > 0) ? "Agonist" : ($actual < 0 ? "Inverse Agonist" : "Non-Agonist");

$average["Prediction"] = $prediction;
$average["Actual"] = $actual;

echo "Predicted ligand activity: $prediction\n";
echo "Empirical ligand activity: $actual\n";

if (@$_REQUEST['saved'])
	print_r($average);
else
{
	// Reload to prevent overwriting another process' output.
	if (file_exists($json_file)) $dock_results = json_decode(file_get_contents($json_file), true);
	$dock_results[$protid][$ligname] = $average;
	
	$keys = array_keys($dock_results);
	natsort($keys);
	$out_results = [];
	foreach ($keys as $key) $out_results[$key] = $dock_results[$key];
	
	$f = fopen($json_file, "wb");
	if (!$f) die("File write FAILED. Make sure have access to write $json_file.");
	fwrite($f, json_encode_pretty($out_results));
	fclose($f);
}

















