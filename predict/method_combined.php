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
require("matrix9.php");

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

ensure_sdf_exists($ligname);

echo "Beginning dock of $ligname in $protid...\n\n";
$fam = family_from_protid($protid);

switch ($fam)
{
	case "TAAR":	
	$capturer = resno_from_bw($protid, "45.50");
	$shuttler = resno_from_bw($protid, "45.51");

	$shelf1   = resno_from_bw($protid, "5.29");
	$shelf2   = resno_from_bw($protid, "5.31");

	$flank    = resno_from_bw($protid, "2.64");

	$acid     = resno_from_bw($protid, "3.32");
	$bind     = ($protid == "TAAR6" || $protid == "TAAR8") ? resno_from_bw($protid, "5.40") : "";
	$toggle   = resno_from_bw($protid, "6.48");

	$res644   = resno_from_bw($protid, "6.44");
	$res343   = resno_from_bw($protid, "3.43");

	$acvbrots = <<<heredoc
ACVBROT $toggle CB CG 180
ACVBROT $res644 CA CB 60
ACVBROT $res343 CA CB -60
heredoc;

	if ($protid == "TAAR1" || $protid == "TAAR2")							// TAAR2 and TAAR1 have an insertion in TMR3.
	{
		$shelf1--;
		$shelf2--;
	}
	
	$cenres = "CEN RES $capturer $shuttler";
	$path1 = "PATH 1 RES $shuttler $shelf1 $shelf2";
	$path2 = "PATH 2 RES $flank $shelf2";
	$path3 = "PATH 3 RES $acid $bind $toggle";
	break;
	
	default:
	die("OR predictions will not work until the EXR2 HxxCE motif has been helixed and metal bound.\n");
	
	$captmtl1 = resno_from_bw($protid, "45.47");
	$captmtl2 = resno_from_bw($protid, "45.50");
	$captmtl3 = resno_from_bw($protid, "45.51");
	
	$cenres = "CEN RES $captmtl1 $captmtl2 $captmtl3";
	$path1 = "";
	$path2 = "";
	$path3 = "";
}


$tmr2start = $prots[$protid]['region']['TMR2']['start'];
$cyt1end = $tmr2start - 1;
$tmr4end = $prots[$protid]['region']['TMR4']['end'];
$tmr5start = $prots[$protid]['region']['TMR5']['start'];
$exr2start = $tmr4end + 1;
$exr2end = $tmr5start - 1;

$acv_matrix = "";
foreach ($matrix as $region => $values) $acv_matrix .= "ACVMX $region " . implode(" ", $values)."\n";

chdir(__DIR__);
chdir("..");
if (!file_exists("output")) mkdir("output");
if (!file_exists("output/$fam")) mkdir("output/$fam");
if (!file_exists("output/$fam/$protid")) mkdir("output/$fam/$protid");
if (!file_exists("output/$fam/$protid")) die("Failed to create output folder.\n");

$configf = <<<heredoc

PROT pdbs/$fam/$protid.upright.pdb

LIG sdf/$ligname.sdf

$cenres
$path1
$path2
$path3
PATH 4 REL 0 0 0

ACVNODE 4

$acv_matrix

$acvbrots

SIZE 7.0 7.5 7.0

POSE 10
ITER 50

# DIFF
ELIM 20

OUT output/$fam/$protid/$protid-$ligname.pred.dock

ECHO


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
		$outlines = [];
		passthru("bin/primarydock tmp/prediction.config");
		$outlines = explode("\n", file_get_contents("output/$protid-$ligname.pred.dock"));
		if (count($outlines) >= 100) break;
	}
}

if (@$_REQUEST['echo']) echo implode("\n", $outlines) . "\n\n";

if (count($outlines) < 100) die("Docking FAILED.\n");

$benerg = [];
$lprox = [];
$pose = false;
$node = -1;
$poses_found = 0;
$full_poses = 0;

foreach ($outlines as $ln)
{
	if (substr($ln, 0, 6) == "Pose: ") $pose = intval(explode(" ", $ln)[1]);
	if (substr($ln, 0, 6) == "Node: ") $node = intval(explode(" ", $ln)[1]);
	if (trim($ln) == "TER")
	{
		$pose = false;
		$node = -1;
	}
	
	if ($pose && $node>=0 && substr($ln, 0,  7) == "Total: "    )
	{
		$benerg[$pose][$node] = floatval(explode(" ", $ln)[1]);
		if ($pose && $node==3 && $benerg[$pose][$node] < 0) $full_poses++;
	}
	if ($pose && $node>=0 && substr($ln, 0, 11) == "Proximity: ") $lprox[ $pose][$node] = floatval(explode(" ", $ln)[1]);
	if (false !== strpos($ln, "pose(s) found")) $poses_found = intval($ln);
}

// print_r($benerg);

$sum = [];
$sump = [];
$count = [];
foreach ($benerg as $pose => $data)
{
	foreach ($data as $node => $value)
	{
		if (!isset($sum[$node]  )) $sum[$node] = 0.0;
		if (!isset($sump[$node] )) $sump[$node] = 0.0;
		if (!isset($count[$node])) $count[$node] = 0;
		
		$sum[$node] += $value;
		$count[$node]++;

		$sump[$node] += $lprox[$pose][$node];
	}
}

$average = [];

$average['Poses'] = $poses_found;
$average['Full poses'] = $full_poses;

foreach ($sum as $node => $value)
{
	$average["Node $node"] = round($value / (@$count[$node] ?: 1), 3);
	// $average["Proximity $node"] = round($sump[$node] / (@$count[$node] ?: 1), 3);
}

$capture = max(-$average["Node 0"], -$average["Node 1"], -$average["Node 2"]);
$completion = floatval($full_poses) / $poses_found;
$completion *= (max(0.001, @-$average["Node 3"]) / $capture);
// $heldback = max($average["Proximity 0"], $average["Proximity 1"], $average["Proximity 2"], $average["Proximity 3"]);
$acvratio = -$average["Node 4"] / (-$average["Node 3"] ?: 0.001);

$prediction = "Non-Agonist";
if ($capture >= 20 && $completion >= 0.75)
{
	if ($acvratio < 0.75) $prediction = "Inverse Agonist";
	else $prediction = "Agonist";
}

$actual = best_empirical_pair($protid, $ligname);
if ($actual > $sepyt["?"]) $actual = ($actual > 0) ? "Agonist" : ($actual < 0 ? "Inverse Agonist" : "Non-Agonist");
else $actual = "(unknown)";

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

	ksort($dock_results[$protid]);
	
	$keys = array_keys($dock_results);
	natsort($keys);
	$out_results = [];
	foreach ($keys as $key) $out_results[$key] = $dock_results[$key];
	
	$f = fopen($json_file, "wb");
	if (!$f) die("File write FAILED. Make sure have access to write $json_file.");
	fwrite($f, json_encode_pretty($out_results));
	fclose($f);
}

















