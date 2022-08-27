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

// Configurable variables
$dock_retries = 5;
$max_simultaneous_docks = 4;	// If running this script as a cron, we recommend setting this to half the number of physical cores.

// Set the matrix to immutable constant values derived from the btRho activation.
$matrix =
[
	"TMR1" => [ -1.432411, -2.746519,  2.509232,  2.363341, -3.109768, -3.402768 ],
	"TMR2" => [  4.639589,  0.024734, -2.108768, -1.090161, -2.695269,  1.679482 ],
	"TMR3" => [  0.167589, -2.024517,  3.395732,  4.231589,  1.114231, -4.413017 ],
	"TMR4" => [  1.722089,  0.166982, -1.611769, -2.093911,  1.330481,  2.336232 ],
	"TMR5" => [ -2.133661,  3.529982,  1.032483, -3.523911,  3.119732,  2.745732 ],
	"TMR6" => [ -2.557410,  1.492481, -5.381518, -2.909160,  3.353733,  2.692482 ],
	"TMR7" => [ -0.405411,  0.009982,  2.163982,  3.021839, -3.566268, -1.637518 ],
];


// btRho matrix adjustments
// Adjust TMR6 to not clash with TMR5.
$matrix["TMR6"][2] += 3;
$matrix["TMR6"][5] += 1;

// Equalize Y displacements.
foreach ($matrix as $region => &$values)
{
	$y = ($values[1] + $values[4])/2;
	$values[1] = $values[4] = $y;
}

// No increase in distance between TMR5 and TMR3 at the extracellular end.
$matrix["TMR3"][0] -= 3;
$matrix["TMR3"][2] -= 1;
$matrix["TMR5"][0] += 3;
$matrix["TMR5"][2] -= 1;

// Prevent clash of TMR7 with TMR2 and TMR1. Bring TMR7 closer to TMR3.
$matrix["TMR7"][3] -= 3;
$matrix["TMR7"][5] -= 2;


# Activation matrix from ADRB2.pepd:
$matrix =
[
	'TMR1' => [  1.928375, -0.231194,  1.514458, -1.745764,  3.297985, -0.673183 ],
	'TMR2' => [ -1.637074,  2.562085, -0.780735,  6.004886, -0.640827,  1.471912 ],
	'TMR3' => [  3.991817,  1.599894,  0.494513, -4.396341, -0.118191, -2.558083 ],
	'TMR4' => [ -3.715700,  4.363474, -0.580147,  4.709310,  0.645814,  0.620629 ],
	'TMR5' => [  3.459619, -0.976376,  0.259171, -6.087701, -1.455206, -1.813480 ],
	'TMR6' => [ -8.180727, -1.155262, -2.628376,  3.851160, -5.411479,  0.728764 ],
	'TMR7' => [  4.153688, -3.598364,  1.721104, -2.335531,  1.117605,  2.223419 ],
];


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
	$cenr3_32 = resno_from_bw($protid, "3.32");
	$cenr5_43 = resno_from_bw($protid, "5.43");
	$cenr6_48 = resno_from_bw($protid, "6.48");
	$cenr6_52 = resno_from_bw($protid, "6.52");
	$cenr7_43 = resno_from_bw($protid, "7.43");
	
	$cenres = "CEN RES $cenr3_32 $cenr5_43 $cenr6_52 $cenr7_43 $cenr6_48";
	break;
	
	default:
	$cenr3_33 = resno_from_bw($protid, "3.33");
	$cenr3_36 = $cenr3_33 + 3;
	$cenr5_43 = resno_from_bw($protid, "5.43");
	$cenr5_47 = $cenr5_43 + 4;
	$cenr6_40 = resno_from_bw($protid, "6.40");
	$cenr6_47 = $cenr6_40 + 7;
	
	$cenres = "CEN RES $cenr3_33 $cenr3_36 $cenr5_43 $cenr5_47 $cenr6_40 $cenr6_47";
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

SIZE 7.0 7.5 7.0

EXCL 1 $cyt1end		# Head, TMR1, and CYT1.
EXCL $exr2start $exr2end	# EXR2 between TMR4 and TMR5.

POSE 10
ITER 50

# DIFF
# ELIM 100

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

















