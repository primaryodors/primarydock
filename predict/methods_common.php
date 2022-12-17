<?php

global $dock_metals, $bias_by_energy, $dock_results, $pdbfname;

// Includes
require("protutils.php");
require("odorutils.php");
require("dock_eval.php");

echo date("Y-m-d H:i:s.u\n");

$method = explode("_", $_SERVER['PHP_SELF'])[1];
$method = explode(".", $method)[0];


// Configurable variables
$dock_retries = 5;
$max_simultaneous_docks = 2;	// If running this script as a cron, we recommend setting this to no more than half the number of physical cores.

// Load data
$dock_results = [];
$json_file = "predict/dock_results_$method.json";

// Version
chdir(__DIR__);
$version = filemtime("method_$method.php");

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

if (@$_REQUEST['simul']) $max_simultaneous_docks = intval($_REQUEST['simul']) ?: 2;

if (@$_REQUEST['next'])
{
	$cmd = "ps -ef | grep ':[0-9][0-9] bin/primarydock' | grep -v grep";
	exec($cmd, $results);
	if (!@$_REQUEST['force'] && trim(@$results[$max_simultaneous_docks-1])) die("Already running.\n".print_r($results, 1));
	$already = implode("\n", $results);
	
	$protid = @$_REQUEST['prot'] ?: false;
	$ligname = @$_REQUEST['lig'] ?: false;
	
	foreach (array_keys($prots) as $rcpid)
	{
		if ($protid && $protid != $rcpid) continue;
		$odorids = array_keys(all_empirical_pairs_for_receptor($rcpid));
		// shuffle($odorids);

		foreach ($odorids as $oid)
		{
			$full_name = $odors[$oid]['full_name'];
			$fnnospace = str_replace(" ", "_", $full_name);
			$lignospace = str_replace(" ", "_", $full_name);
			if ((!isset($dock_results[$rcpid][$full_name]) && !isset($dock_results[$rcpid][$fnnospace]))
                ||
                max(@$dock_results[$rcpid][$full_name]['version'], @$dock_results[$rcpid][$fnnospace]['version']) < $version
                )
			{
				if (false!==strpos($already, $lignospace))
				{
					continue;
				}
				
				$protid = $rcpid;
				$ligname = $full_name;
				
				goto found_next_pair;
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


ensure_sdf_exists($ligname);

echo "Beginning dock of $ligname in $protid...\n\n";
$fam = family_from_protid($protid);


extract(binding_site($protid));

$outfname = "output.dock";

function prepare_outputs()
{
    global $ligname, $dock_metals, $protid, $fam, $outfname, $pdbfname;

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
    $outfname = "output/$fam/$protid/$protid-$ligname.dock";
}

function process_dock()
{
    global $ligname, $protid, $configf, $dock_retries, $outfname, $bias_by_energy, $version, $sepyt, $json_file;
    if (!file_exists("tmp")) mkdir("tmp");
    $lignospace = str_replace(" ", "", $ligname);
    $cnfname = "tmp/prediction.$protid.$lignospace.config";
    $f = fopen($cnfname, "wb");
    if (!$f) die("File write error. Check tmp folder is write enabled.\n");
    
    // echo $configf;

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
    $scenerg = [];
    $vdwrpl = [];
    $dosce = false;
    $dovdw = false;
    $pose = false;
    $node = -1;
    $poses_found = 0;
    $ca_loc = [];
    $sc_loc = [];
    $sc_qty = [];
    foreach ($outlines as $ln)
    {
        if (!strlen(trim($ln))) $dosce = false;
        if (substr($ln, 0, 6) == "Pose: ")
        {
            $pose = intval(explode(" ", $ln)[1]);
            $dovdw = false;
        }
        if (substr($ln, 0, 6) == "Node: ") $node = intval(explode(" ", $ln)[1]);
        if (trim($ln) == "TER")
        {
            $pose = false;
            $node = -1;
        }
    
        if ($dosce)
        {
            $pettia = explode(": ", $ln);
            $resno = intval(substr($pettia[0], 3));
            $e = floatval($pettia[1]);
            $scenerg[$pose][$node][$resno] = $e;
            continue;
        }

        if ($dovdw)
        {
            $pettia = explode(": ", $ln);
            $resno = intval(substr($pettia[0], 3));
            $v = floatval($pettia[1]);
            $vdwrpl[$pose][$node][$resno] = $v;
            continue;
        }
    
        if (substr($ln, 0, 7) == 'BENERG:') $dosce = true;
        if (substr($ln, 0, 7) == 'vdWRPL:') $dovdw = true;
        
        if ($pose && $node>=0 && substr($ln, 0,  7) == "Total: "    )
        {
            $benerg[$pose][$node] = floatval(explode(" ", $ln)[1]);
            $dosce = false;
            continue;
        }
        
        $bias = $bias_by_energy ? max(-$benerg[$pose][$node], 1) : 1;
    
        if (false !== strpos($ln, "pose(s) found")) $poses_found = intval($ln);
    
        if (substr($ln, 0, 5) == 'ATOM ')
        {
            $ln = preg_replace("/\\s+/", " ", trim($ln));
            $ln = explode(" ", $ln);
    
            $aname = $ln[2];
            $resno = intval($ln[4]);
            $x = floatval($ln[5]);
            $y = floatval($ln[6]);
            $z = floatval($ln[7]);
    
            switch ($aname)
            {
                case 'N': case 'H': case 'HN': case 'HA': case 'HB': case 'C': case 'O':
                break;
    
                case 'CB': case 'HB': case 'HB1': case 'HB2': case 'HB3':
                break;
    
                case 'CA':
                $ca_loc[$resno] = [$x,$y,$z];
                if (!isset($sc_loc[$resno]))
                {
                    $sc_loc[$resno] = [0.0,0.0,0.0];
                    $sc_qty[$resno] = 0.0;
                }
                break;
    
                default:
                $sc_loc[$resno][0] += ($x - $ca_loc[$resno][0]) * $bias;
                $sc_loc[$resno][1] += ($y - $ca_loc[$resno][1]) * $bias;
                $sc_loc[$resno][2] += ($z - $ca_loc[$resno][2]) * $bias;
                $sc_qty[$resno] += $bias;
            }
        }
    }
    
    $sum = [];
    $count = [];
    $ssce = [];
    $svdw = [];
    $rsum = [];
    $rsumv = [];
    foreach ($benerg as $pose => $data)
    {
        foreach ($data as $node => $value)
        {
            if (!isset($sum[$node]  )) $sum[$node] = 0.0;
            if (!isset($count[$node])) $count[$node] = 0;
            
            $sum[$node] += $value;
            $count[$node]++;
    
            foreach ($scenerg[$pose][$node] as $resno => $e)
            {
                if (!isset($ssce[$resno])) $ssce[$resno] = floatval($rsum[$resno] = 0);
                $ssce[$resno] += $e;
                $rsum[$resno]++;
            }
    
            foreach ($vdwrpl[$pose][$node] as $resno => $v)
            {
                if (!isset($svdw[$resno])) $svdw[$resno] = floatval($rsumv[$resno] = 0);
                $svdw[$resno] += $v;
                $rsumv[$resno]++;
            }
        }
    }
    
    $sc_avg = [];
    
    $sce = [];
    $vdw = [];
    foreach ($ca_loc as $resno => $a)
    {
        $bw = bw_from_resno($protid, $resno);
    
        if ($sc_qty[$resno]) $sc_avg[$bw] =
        [
            round($sc_loc[$resno][0] / $sc_qty[$resno], 3),
            round($sc_loc[$resno][1] / $sc_qty[$resno], 3),
            round($sc_loc[$resno][2] / $sc_qty[$resno], 3)
        ];
    
        if (@$ssce[$resno] && $rsum[$resno ]) $sce[$bw] = $ssce[$resno] / $rsum[$resno];
        if (@$svdw[$resno] && $rsumv[$resno]) $vdw[$bw] = $svdw[$resno] / $rsumv[$resno];
    }

    
    $average = [];
    $average['version'] = $version;
    $average['Poses'] = $poses_found;
    
    foreach ($sce as $bw => $e)
    {
        $average["BEnerg $bw"] = $e;
    }
    
    foreach ($vdw as $bw => $v)
    {
        $average["vdWrpl $bw"] = $v;
    }
    
    foreach ($sum as $node => $value)
    {
        $average["Node $node"] = round($value / (@$count[$node] ?: 1), 3);
    }
    
    ksort($sc_avg);
    
    foreach ($sc_avg as $bw => $xyz)
    {
        $average["SCW $bw"] = $xyz;
    }
    
    $actual = best_empirical_pair($protid, $ligname);
    if ($actual > $sepyt["?"]) $actual = ($actual > 0) ? "Agonist" : ($actual < 0 ? "Inverse Agonist" : "Non-Agonist");
    else $actual = "(unknown)";
    
    $average["Actual"] = $actual;
    
    
    // Reload to prevent overwriting another process' output.
    if (file_exists($json_file)) $dock_results = json_decode(file_get_contents($json_file), true);
    $dock_results[$protid][$ligname] = $average;
    
    echo "Loaded $json_file with ".count($dock_results)." records.\n";
    
    ksort($dock_results[$protid]);
    
    $keys = array_keys($dock_results);
    natsort($keys);
    $out_results = [];
    foreach ($keys as $key) $out_results[$key] = $dock_results[$key];
    
    echo "Writing $json_file with ".count($dock_results)." records.\n";
    
    $f = fopen($json_file, "wb");
    if (!$f) die("File write FAILED. Make sure have access to write $json_file.");
    fwrite($f, json_encode_pretty($out_results));
    fclose($f);
    
}