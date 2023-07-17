<?php

global $dock_metals, $bias_by_energy, $dock_results, $pdbfname, $fam, $do_scwhere, $metrics_to_process;

// Includes
require("cputemp.php");
require("protutils.php");
require("odorutils.php");
require("dock_eval.php");
require_once("statistics.php");

die_if_too_hot();

echo date("Y-m-d H:i:s.u\n");

$method = explode("_", $_SERVER['PHP_SELF'])[1];
$method = explode(".", $method)[0];

$result = $output = false;
chdir(__DIR__);
chdir("..");
exec("make primarydock", $output, $result);
if ($result) die("Build fail.\n".print_r($output, true));

// Configurable variables
$dock_retries = 5;
$max_simultaneous_docks = 2;	// If running this script as a cron, we recommend setting this to no more than half the number of physical cores.
if (!isset($do_scwhere)) $do_scwhere = false;
$metrics_to_process =
[
  "BENERG" => "BindingEnergy",
  "BENERG.rgn" => "BindingEnergy.rgn",
];

// Load data
$dock_results = [];
$json_file = "predict/dock_results_$method.json";

// Version
chdir(__DIR__);
$version = max(filemtime("method_$method.php"), filemtime("methods_common.php"), filemtime("../bin/primarydock"));

chdir("..");
if (file_exists($json_file))
{
	$dock_results = json_decode(file_get_contents($json_file), true);
}

$prioritize_pairs = [];
$json_file_tmp = "predict/predict_test_pairs.json";
if (file_exists($json_file_tmp))
{
	$prioritize_pairs = json_decode(file_get_contents($json_file_tmp), true);
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
	
	foreach ($prioritize_pairs as $rcpid => $liglist)
	{
	    if ($protid && $protid != $rcpid && !preg_match("/^$protid$/", $rcpid) ) continue;
	    foreach ($liglist as $lig)
	    {
	        $odor = find_odorant($lig);
			$full_name = $odor['full_name'];
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
	
	foreach (array_keys($prots) as $rcpid)
	{
		if ($protid && $protid != $rcpid && !preg_match("/^$protid$/", $rcpid) ) continue;
		// if (!$protid && $rcpid == 'OR1A1') continue;
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

function process_dock($metrics_prefix = "", $noclobber = false)
{
    global $ligname, $protid, $configf, $dock_retries, $outfname, $metrics_to_process, $bias_by_energy, $version, $sepyt, $json_file, $do_scwhere;
    if (!file_exists("tmp")) mkdir("tmp");
    $lignospace = str_replace(" ", "", $ligname);
    $cnfname = "tmp/prediction.$protid.$lignospace.config";
    $f = fopen($cnfname, "wb");
    if (!$f) die("File write error. Check tmp folder is write enabled.\n");

    // echo $configf;

    fwrite($f, $configf);
    fclose($f);

    $retvar = 0;

    $outlines = [];
    if (@$_REQUEST['saved'])
    {
        $outlines = explode("\n", file_get_contents($outfname)); // $_REQUEST['saved']));
    }
    else
    {
    	$elim = 0;
        for ($try = 0; $try < $dock_retries; $try++)
        {            
            set_time_limit(300);
            $outlines = [];
            $cmd = "bin/primarydock \"$cnfname\"";
            if ($elim) $cmd .= " --elim $elim";
            echo "$cmd\n";
            passthru($cmd, $retvar);
            $outlines = explode("\n", file_get_contents($outfname));
            if (count($outlines) >= 200) break;
            if (!$elim) $elim = 99;
            else $elim *= 1.333;
        }
    }

    if (@$_REQUEST['echo']) echo implode("\n", $outlines) . "\n\n";

    if ($retvar || (count($outlines) < 100)) die("Docking FAILED.\n");
    
    unlink($cnfname);

    $mode = "";
    $pose = 0;
    $node = -1;
    $outdqty = [];
    
    if ($noclobber)
    {
        if (file_exists($json_file)) $dock_results = json_decode(file_get_contents($json_file), true);
        $outdata = @$dock_results[$protid][$ligname] ?: [];
    }
    else $outdata = [];

    $weight = [];
    $pose_node_has_weight = [];

    foreach ($outlines as $ln)
    {
        $coldiv = explode(":", $ln);

        if (trim($ln) == "TER")
        {
            $pose = false;
            $node = -1;
            continue;
        }

        if (count($coldiv) == 2)
        {
            if ($coldiv[0] == "Pose")
            {
                $pose = intval($coldiv[1]);
                $node = -1;
                continue;
            }
            else if ($coldiv[0] == "Node")
            {
                $node = intval($coldiv[1]);
                continue;
            }
            else if ($coldiv[0] == 'Total' && !$node)
            {
                $e = -floatval($coldiv[1]);
                if ($e < 1) $e = 1.0 / abs(log(abs($e)));
                $weight[$pose] = $e;
            }
        }
    }
    
    foreach ($outlines as $ln)
    {
        $coldiv = explode(":", $ln);

        if (trim($ln) == "TER")
        {
            $pose = false;
            $node = -1;
            continue;
        }

        if (count($coldiv) == 2)
        {
            if ($coldiv[0] == "Pose")
            {
                $pose = intval($coldiv[1]);
                $node = -1;
                continue;
            }
            else if ($coldiv[0] == "Node")
            {
                $node = intval($coldiv[1]);
                continue;
            }
            else if ($coldiv[0] == 'BENERG' || $coldiv[0] == 'vdWRPL')
            {
                $mode = $coldiv[0];
                continue;
            }
            else if ($pose && $node>=0)
            {
                if (preg_match("/[A-Z][a-z]{2}[0-9]{1,4}/", $coldiv[0]))
                {
                    $resno = intval(substr($coldiv[0], 3));
                    $region = intval(bw_from_resno($protid, $resno));
                    if (isset($metrics_to_process["$mode.rgn"]))
                    {
                        $wmode = str_replace(".rgn", ".$region", $metrics_to_process["$mode.rgn"]);
                        if (!isset($outdata[$wmode])) $outdata[$wmode] = 0.0;
                        $outdata[$wmode] += floatval($coldiv[1]) * $weight[$pose];
                        if (!@$pose_node_has_weight[$wmode][$pose][$node])
                        {
                            if (!isset($outdqty[$wmode])) $outdqty[$wmode] = $weight[$pose];
                            else $outdqty[$wmode] += $weight[$pose];
                            $pose_node_has_weight[$wmode][$pose][$node] = true;
                        }
                    }
                    continue;
                }
                if ($coldiv[0] == "Total")
                {
                    if (isset($metrics_to_process[$mode]))
                    {
                        $wmode = $metrics_to_process[$mode];
                        if (!isset($outdata[$wmode])) $outdata[$wmode] = 0.0;
                        $outdata[$wmode] += floatval($coldiv[1]);
                        if (!isset($outdqty[$wmode])) $outdqty[$wmode] = 1;
                        else $outdqty[$wmode]++;
                    }
                    continue;
                }
                else if ($coldiv[0] == "Ligand polar satisfaction")
                {
                    $mode = "POLSAT";
                    if (isset($metrics_to_process[$mode]))
                    {
                        $wmode = $metrics_to_process[$mode];
                        if (!isset($outdata[$wmode])) $outdata[$wmode] = 0.0;
                        $outdata[$wmode] += floatval($coldiv[1]) * $weight[$pose];
                        if (!isset($outdqty[$wmode])) $outdqty[$wmode] = $weight[$pose];
                        else $outdqty[$wmode] += $weight[$pose];
                    }
                    continue;
                }
                else if ($coldiv[0] == "Protein clashes")
                {
                    $mode = "PCLASH";
                    if (isset($metrics_to_process[$mode]))
                    {
                        $wmode = $metrics_to_process[$mode];
                        if (!isset($outdata[$wmode])) $outdata[$wmode] = 0.0;
                        $outdata[$wmode] += floatval($coldiv[1]) * $weight[$pose];
                        if (!isset($outdqty[$wmode])) $outdqty[$wmode] = $weight[$pose];
                        else $outdqty[$wmode] += $weight[$pose];
                    }
                    continue;
                }
                else if (strpos($coldiv[0], " active theta: "))
                {
                    $morceaux = explode(' ', $coldiv[0]);
                    $mode = "ACVTH.{$morceaux[0]}";
                    if (isset($metrics_to_process[$mode]))
                    {
                        $wmode = $metrics_to_process[$mode];
                        if (!isset($outdata[$wmode])) $outdata[$wmode] = 0.0;
                        $outdata[$wmode] += floatval($coldiv[1]) * $weight[$pose];
                        if (!isset($outdqty[$wmode])) $outdqty[$wmode] = $weight[$pose];
                        else $outdqty[$wmode] += $weight[$pose];
                    }
                    continue;
                }
            }
        }

        if (false !== strpos($ln, "pose(s) found"))
        {
            $mode = "POSES";
            if (isset($metrics_to_process[$mode]))
            {
                $wmode = $metrics_to_process[$mode];
                $outdata[$wmode] = intval($ln);
            }
            continue;
        }

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
                $ca_loc = [$x,$y,$z];
                $region = intval(bw_from_resno($protid, $resno));
                if (isset($metrics_to_process["CALOC.rgn"]))
                {
                    $wmode = str_replace(".rgn", ".$region", $metrics_to_process["CALOC.rgn"]);

                    if (!isset($outdata[$wmode])) $outdata[$wmode] = [0,0,0];
                    for ($x=0; $x<3; $x++) $outdata[$wmode][$x] += floatval(substr($ln, 29+8*x, 8)) * $weight[$pose];

                    if (!isset($outdqty[$wmode])) $outdqty[$wmode] = $weight[$pose];
                    else $outdqty[$wmode] += $weight[$pose];
                }
                break;
    
                default:
                ;
            }
        }
    }
    
    $outdata['version'] = $version;

    $tme = [];
    foreach ($outdata as $k => $v)
    {
        $div = floatval(@$outdqty[$k]);
        if ($div)
        {
            if (is_array($v))
                foreach (array_keys($v) as $k1) $outdata[$k][$k1] /= $div;
            else
                $outdata[$k] = $v / $div;
        }
    }

    $actual = best_empirical_pair($protid, $ligname);
    if ($actual > $sepyt["?"]) $actual = ($actual > 0) ? "Agonist" : ($actual < 0 ? "Inverse Agonist" : "Non-Agonist");
    else $actual = "(unknown)";

    if (function_exists("make_prediction")) $outdata = make_prediction($outdata);
    $outdata["Actual"] = $actual;


    // Reload to prevent overwriting another process' output.
    if (file_exists($json_file)) $dock_results = json_decode(file_get_contents($json_file), true);
    $dock_results[$protid][$ligname] = $outdata;
    
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
