<?php

// dock.php
//
// Performs a dock of an odorant in a receptor using both the inactive and active PDB files.
//
// Example call syntax:
// php -f predict/dock.php prot=OR1A1 lig=geraniol
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
	$protid = @$_REQUEST['prot'] ?: "OR1A1";
	
	$ligname = @$_REQUEST['lig'] ?: "geraniol";
}

echo "Beginning dock of $ligname in $protid...\n\n";
$fam = family_from_protid($protid);

switch ($fam)
{
	case "TAAR":
	$matrix = $matrixADRB2;
	
	$cenr3_32 = resno_from_bw($protid, "3.32");
	$cenr5_43 = resno_from_bw($protid, "5.43");
	$cenr6_48 = resno_from_bw($protid, "6.48");
	$cenr6_52 = resno_from_bw($protid, "6.52");
	$cenr7_43 = resno_from_bw($protid, "7.43");
	
	$cenres = "CEN RES $cenr3_32 $cenr5_43 $cenr6_52 $cenr7_43 $cenr6_48";
	break;
	
	default:
	// $matrix = $matrixbtRho;
	$matrix = $matrixADRB2;
	
	$cenr3_33 = resno_from_bw($protid, "3.33");
	$cenr3_36 = $cenr3_33 + 3;
	$cenr5_43 = resno_from_bw($protid, "5.43");
	$cenr5_47 = $cenr5_43 + 4;
	$cenr6_40 = resno_from_bw($protid, "6.40");
	$cenr6_47 = $cenr6_40 + 7;
	
	$cenres = "CEN RES $cenr3_33 $cenr3_36 $cenr5_43 $cenr5_47 $cenr6_40 $cenr6_47";
}

// If the matrix has a global skew, neutralize it.
$up = [0.0,0.0,0.0];
$dn = [0.0,0.0,0.0];
$upq = $dnq = 0;
foreach ($matrix as $region => $coords)
{
	if (substr($region, 0, 3) == "TMR")
	{
		$tmrno = intval(substr($region, 3));
		if ($tmrno & 1)			// Odd numbered TMRs descend, so "up" is N-terminus and "down" is C-terminus.
		{
			for ($j=0; $j<2; $j++) $up[$j] += $coords[$j];
			$upq++;
			
			if ($tmrno != 5)
			{
				for ($j=0; $j<2; $j++) $dn[$j] += $coords[$j+3];
				$dnq++;
			}
		}
		else					// Even numbered TMRs go the other way.
		{
			for ($j=0; $j<2; $j++) $up[$j] += $coords[$j+3];
			$upq++;
			
			if ($tmrno != 6)
			{
				for ($j=0; $j<2; $j++) $dn[$j] += $coords[$j];
				$dnq++;
			}
		}
	}
}

if ($upq) for ($j=0; $j<2; $j++) $up[$j] /= $upq;
if ($dnq) for ($j=0; $j<2; $j++) $dn[$j] /= $dnq;

foreach ($matrix as $region => $coords)
{
	if (substr($region, 0, 3) == "TMR")
	{
		$tmrno = intval(substr($region, 3));
		if ($tmrno & 1)
		{
			for ($j=0; $j<2; $j++) $matrix[$region][$j  ] -= $up[$j];
			for ($j=0; $j<2; $j++) $matrix[$region][$j+3] -= $dn[$j];
		}
		else
		{
			for ($j=0; $j<2; $j++) $matrix[$region][$j  ] -= $dn[$j];
			for ($j=0; $j<2; $j++) $matrix[$region][$j+3] -= $up[$j];
		}
		
		// Flatten the motion to the XZ plane.
		//$matrix[$region][1] = $matrix[$region][4] = 0;
	}
	
	foreach ($coords as $i => $c) $matrix[$region][$i] = round($matrix[$region][$i], 3);
}


$tmr2start = $prots[$protid]['region']['TMR2']['start'];
$cyt1end = $tmr2start - 1;
$tmr4end = $prots[$protid]['region']['TMR4']['end'];
$tmr5start = $prots[$protid]['region']['TMR5']['start'];
$exr2start = $tmr4end + 1;
$exr2end = $tmr5start - 1;


$acv_matrix = "";
foreach ($matrix as $region => $values) $acv_matrix .= "ACVMX $region " . implode(" ", $values)."\n";

$configf = <<<heredoc

PROT pdbs/$fam/$protid.rotated.pdb

LIG sdf/$ligname.sdf

$cenres
PATH 1 REL 0 0 0

# NODEPDB 1 pdbs/$fam/$protid.active.pdb

# Activation Matrix.
# ACVMX region N.x N.y N.z C.x C.y C.z
$acv_matrix

ACVNODE 1

SIZE 8.0 8.5 8.0

EXCL 1 $cyt1end		# Head, TMR1, and CYT1.
EXCL $exr2start $exr2end	# EXR2 between TMR4 and TMR5.

POSE 10
ITER 50

# DIFF
# ELIM 100

OUT output/$protid/$protid-$ligname.pred.dock



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
	$average[$node?"Active":"Inactive"] = round($value / (@$count[$node] ?: 1), 3);
}

if (min($average) > 0 || !isset($average["Active"]) || !isset($average["Inactive"]))
{
	$prediction = "Non-Agonist";
}
else
{
	if ($average["Active"] <= 1.05 * $average["Inactive"])
		$prediction = "Agonist";
	else if ($average["Inactive"] <= 1.1 * $average["Active"])
		$prediction = "Inverse Agonist";
	else $prediction = "Non-Agonist";
}

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
	
	$f = fopen($json_file, "wb");
	if (!$f) die("File write FAILED. Make sure have access to write $json_file.");
	fwrite($f, json_encode_pretty($dock_results));
	fclose($f);
}

















