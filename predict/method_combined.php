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
require("dock_eval.php");

// Configurable variables
$dock_retries = 5;
$max_simultaneous_docks = 2;	// If running this script as a cron, we recommend setting this to no more than half the number of physical cores.
$dock_metals = false;

// Load data
$dock_results = [];
$json_file = "predict/dock_results.json";

chdir(__DIR__);
chdir("..");
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
	
	foreach (array_keys($prots) as $rcpid)
	{
		if ($protid && $protid != $rcpid) continue;
		$odorids = array_keys(all_empirical_pairs_for_receptor($rcpid));

		foreach ($odorids as $oid)
		{
			$full_name = $odors[$oid]['full_name'];
			$fnnospace = str_replace(" ", "_", $full_name);
			if (!isset($dock_results[$rcpid][$full_name]) && !isset($dock_results[$rcpid][$fnnospace]))
			{
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
	/*
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
	*/
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
	$cmd = "bin/pepteditor predict/acidfinder.pepd pdbs/$fam/$protid.upright.pdb";

	$captures = [];
	$pathres = [];
	$binding = [];
	$outlines = [];

	exec($cmd, $outlines);

	foreach ($outlines as $ol)
	{
		$ol = explode(" ", $ol);
		switch ($ol[0])
		{
			case 'CAPTURE':
			$captures[] = intval(substr(@$ol[2], 1));
			break;

			case 'PATH':
			$pathres[] = intval(substr(@$ol[2], 1));
			break;

			case 'BIND':
			$binding[] = intval(substr(@$ol[2], 1));
			break;

			default:
			;
		}
	}

	$toggle   = resno_from_bw($protid, "6.48");
	$res644   = resno_from_bw($protid, "6.44");
	$res343   = resno_from_bw($protid, "3.43");

	$acvbrots = <<<heredoc
ACVBROT $toggle CB CG 180
ACVBROT $res644 CA CB 60
ACVBROT $res343 CA CB -60
heredoc;

	if ($protid == "TAAR1" || $protid == "TAAR2")							// TAAR2 and TAAR1 have an insertion in TMR5.
	{
		$shelf1--;
		$shelf2--;
	}
	
	$nodeno = 0;
	$paths = [];
	$cenres = "CEN RES ".implode(" ", $captures);

	$pathwas = $captures[count($captures)-1];
	foreach ($pathres as $pr)
	{
		if ($pathwas)
		{
			$nodeno++; $paths[] = "PATH $nodeno RES $pathwas $pr ";
		}
		$pathwas = $pr;
	}
	$nodeno++; $paths[] = "PATH $nodeno RES ".implode(" ", $binding);

	$pocketnode = $nodeno;
	break;
	
	default:
	$acvbrots = "";

	$captmtl1 = resno_from_bw($protid, "45.47");
	$captmtl2 = resno_from_bw($protid, "45.50");
	$captmtl3 = resno_from_bw($protid, "45.51");
	$captbal1 = resno_from_bw($protid, "3.25");
	$captbal2 = resno_from_bw($protid, "7.35");

	extract(binding_site($protid));
	
	$nodeno = 0;
	$paths = [];	
	$cenres = "CEN RES $captmtl1 $captmtl2 $captmtl3 $captbal1 $captbal2";

	$nodeno++; $paths[] = "PATH $nodeno RES $bsr2a $bsr3a $bsr3b $bsr3c $bsr3d $bsr3e $bsr3f $bsr3g $bsr4a $bsr4b $bsr4c $bsr5a $bsr5b $bsr5c $bsr5d $bsr6a $bsr6b $bsr7a $bsr7b $bsr7c";
	$nodeno++; $paths[] = "PATH $nodeno RES $bsr3a $bsr3b $bsr3c $bsr3d $bsr3e $bsr3f $bsr3g $bsr5a $bsr5b $bsr5c $bsr5d";
	$pocketnode = $nodeno;
	$activenode = $pocketnode + 1;
	$paths[] = "PATH $activenode REL 0 0 0";

	// TODO: Go deeper into the protein for e.g. OR2T11.
}

$paths = implode("\n", $paths);

$tmr2start = $prots[$protid]['region']['TMR2']['start'];
$cyt1end = $tmr2start - 1;
$tmr4end = $prots[$protid]['region']['TMR4']['end'];
$tmr5start = $prots[$protid]['region']['TMR5']['start'];
$exr2start = $tmr4end + 1;
$exr2end = $tmr5start - 1;

$acv_matrix = "";
// foreach ($matrix as $region => $values) $acv_matrix .= "ACVMX $region " . implode(" ", $values)."\n";
foreach ($rotations as $region => $values)
{
	$lregion = substr($region,0,4);
	$sr = $prots[$protid]['region'][$lregion]['start'];
	$er = $prots[$protid]['region'][$lregion]['end'];

	switch($lregion)
	{
		case 'TMR2': $mr = resno_from_bw($protid, "2.53"); break;
		case 'TMR3': $mr = resno_from_bw($protid, "3.36"); break;
		case 'TMR4': $mr = resno_from_bw($protid, "4.57"); break;
		case 'TMR5': $mr = resno_from_bw($protid, "5.43"); break;
		case 'TMR6': $mr = resno_from_bw($protid, "6.48"); break;
		case 'TMR7': $mr = resno_from_bw($protid, "7.46"); break;
		
		default:
		$mr = intval(($sr + $er) / 2);
	}

	if (substr($region, -1) == 'n') $er = $mr;
	if (substr($region, -1) == 'c') $sr = $mr;
	$acv_matrix .= "ACVHXR $region $sr $er $mr " . implode(" ", $values)."\n";
}

chdir(__DIR__);
chdir("..");
if (!file_exists("output")) mkdir("output");
if (!file_exists("output/$fam")) mkdir("output/$fam");
if (!file_exists("output/$fam/$protid")) mkdir("output/$fam/$protid");
if (!file_exists("output/$fam/$protid")) die("Failed to create output folder.\n");

$ligname = str_replace(" ", "_", $ligname);
$suffix = ($dock_metals) ? "metal" : "upright";
$pdbfname = "pdbs/$fam/$protid.$suffix.pdb";
if ($dock_metals && !file_exists($pdbfname)) $pdbfname = "pdbs/$fam/$protid.upright.pdb";
if (!file_exists($pdbfname)) die("Missing PDB.\n");
$outfname = "output/$fam/$protid/$protid-$ligname.pred.dock";

$configf = <<<heredoc

PROT $pdbfname

LIG sdf/$ligname.sdf

$cenres
$paths

ACVNODE $activenode

$acv_matrix

$acvbrots

SIZE 7.0 7.5 7.0

POSE 10
ITER 50

# DIFF
ELIM 50

OUT $outfname

ECHO


heredoc;

if (!file_exists("tmp")) mkdir("tmp");
$lignospace = str_replace(" ", "", $ligname);
$cnfname = "tmp/prediction.$protid.$lignospace.config";
$f = fopen($cnfname, "wb");
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
		passthru("bin/primarydock \"$cnfname\"");
		$outlines = explode("\n", file_get_contents($outfname));
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
		if ($pose && $node==$pocketnode && $benerg[$pose][$node] < 0) $full_poses++;
	}
	if ($pose && $node>=0 && substr($ln, 0, 11) == "Proximity: ") $lprox[ $pose][$node] = floatval(explode(" ", $ln)[1]);
	if (false !== strpos($ln, "pose(s) found")) $poses_found = intval($ln);
}

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

$prediction = evaluate_result($average);

$actual = best_empirical_pair($protid, $ligname);
if ($actual > $sepyt["?"]) $actual = ($actual > 0) ? "Agonist" : ($actual < 0 ? "Inverse Agonist" : "Non-Agonist");
else $actual = "(unknown)";

$average["Active node"] = $activenode;
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

















