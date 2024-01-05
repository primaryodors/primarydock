<?php

global $dock_metals, $bias_by_energy, $dock_results, $pdbfname, $fam, $do_scwhere, $metrics_to_process, $clashcomp, $best_energy;

// Includes
chdir(__DIR__);
chdir("..");
require("predict/cputemp.php");
require("data/protutils.php");
require("data/odorutils.php");
require("predict/dock_eval.php");
require_once("predict/statistics.php");

die_if_too_hot();

echo date("Y-m-d H:i:s.u\n");

$clashcomp = [];

$method = explode("_", $_SERVER['PHP_SELF'])[1];
$method = explode(".", $method)[0];
echo "Method is $method.\n";

$result = $output = false;
chdir(__DIR__);
chdir("..");
exec("make primarydock", $output, $result);
if ($result) die("Build fail.\n".print_r($output, true));

// Configurable variables
$dock_retries = 1;
$max_simultaneous_docks = 2;	// If running this script as a cron, we recommend setting this to no more than half the number of physical cores.
if (!isset($do_scwhere)) $do_scwhere = false;
$metrics_to_process =
[
  "BENERG" => "BindingEnergy",
  "BENERG.rgn" => "BindingEnergy.rgn",
];

// Load data
$dock_results = [];
$json_file = "predict/dock_results.json";

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
    $cmd = "ps -ef | grep -E ':[0-9][0-9] bin/obabel' | grep -v grep";
	exec($cmd, $results);
    if (count($results)) exit;

	$cmd = "ps -ef | grep -E ':[0-9][0-9] bin/(ic|fyg)_activate_or' | grep -v grep";
	exec($cmd, $results);
    if (count($results)) exit;

	$cmd = "ps -ef | grep -E ':[0-9][0-9] (bin/primarydock|bin/pepteditor|bin/ic|bin/fyg|obabel)' | grep -v grep";
	exec($cmd, $results);
	if (!@$_REQUEST['force'] && trim(@$results[$max_simultaneous_docks-1])) die("Already running.\n".print_r($results, 1));
	$already = implode("\n", $results);
	
	$protid = @$_REQUEST['prot'] ?: false;
	$ligname = @$_REQUEST['lig'] ?: false;
	
	foreach (array_keys($prots) as $rcpid)
	{
		if ($protid
            && $protid != $rcpid
            && !preg_match("/^$protid$/", $rcpid)
            && (substr($protid, -1) != '*' || substr($rcpid, 0, strlen($protid)-1) != substr($protid, 0, -1) )
            ) continue;
		$odorids = array_keys(all_empirical_pairs_for_receptor($rcpid));

		foreach ($odorids as $oid)
		{
			$full_name = $odors[$oid]['full_name'];
			$fnnospace = str_replace(" ", "_", $full_name);
			$lignospace = str_replace(" ", "_", $full_name);
			if ((!isset($dock_results[$rcpid][$full_name]) && !isset($dock_results[$rcpid][$fnnospace]))
                ||
                (   (max(@$dock_results[$rcpid][$full_name]['version'], @$dock_results[$rcpid][$fnnospace]['version']) < $version)
                )
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

echo "Beginning prediction of $ligname in $protid...\n\n";
$fam = family_from_protid($protid);


extract(binding_site($protid));

$outfname = "output.dock";

function dock_failed()
{
    echo "Docking FAILED.\n";

    if (stream_isatty(STDOUT) && file_exists("predict/soundalert"))
    {
        $play_sound = true;
        $sa = explode(",",file_get_contents("predict/soundalert"));
        if (count($sa) == 2)
        {
            $hfrom = intval($sa[0]);
            $hto = intval($sa[1]);
            $hnow = intval(date("H"));

            if ($hnow < $hfrom || $hnow > $hto) $play_sound = false;
        }

        $hassox = [];
        exec("which sox", $hassox);
        if ($play_sound && count($hassox))
        {
            exec("play fail.mp3 &");
        }
    }

    exit;
}

function prepare_outputs()
{
    global $ligname, $dock_metals, $protid, $fam, $outfname, $pdbfname;
    global $binding_pockets, $cenres_active, $cenres_inactive, $size, $search, $num_std_devs;
    global $atomto, $stcr, $flxr, $mcoord, $mbp, $astcr, $istcr, $aflxr, $iflxr;

    chdir(__DIR__);
    chdir("..");
    if (!file_exists("output")) mkdir("output");
    if (!file_exists("output/$fam")) mkdir("output/$fam");
    if (!file_exists("output/$fam/$protid")) mkdir("output/$fam/$protid");
    if (!file_exists("output/$fam/$protid")) die("Failed to create output folder.\n");

    $binding_pockets = json_decode(file_get_contents("data/binding_pocket.json"), true);
    
    $ligname = str_replace(" ", "_", $ligname);
    $suffix = ($dock_metals) ? "metal" : "upright";
    $pdbfname = "pdbs/$fam/$protid.$suffix.pdb";
    if ($dock_metals && !file_exists($pdbfname)) $pdbfname = "pdbs/$fam/$protid.upright.pdb";
    if (!file_exists($pdbfname)) die("Missing PDB.\n");
    $outfname = "output/$fam/$protid/$protid-$ligname.dock";

    $size = "7.5 7.5 7.5";
    $search = "BB";
    $atomto = [];
    $stcr = "";
    $flxr = "";
    $astcr = "";
    $aflxr = "";
    $istcr = "";
    $iflxr = "";
    $mcoord = "";

    $mbp = false;                       // Matched binding pocket.
    
    if (isset($binding_pockets[$protid])) $mbp = $binding_pockets[$protid];
    else foreach ($binding_pockets as $pocketid => $pocket)
    {
        if (substr($pocketid, -1) == '*' && substr($pocketid, 0, -1) == substr($protid, 0, strlen($pocketid)-1))
        {
            $mbp = $pocket;
            echo "Matched $pocketid via wildcard.\n";
            break;
        }
        else if (preg_match("/^$pocketid\$/", $protid))
        {
            $mbp = $pocket;
            echo "Matched $pocketid via regex.\n";
            break;
        }
    }
    
    if ($mbp)
    {
        foreach ($mbp as $k => $v)
        {
            if (substr($k,0,1)=='#') unset($mbp[$k]);
        }

        if (isset($mbp['odorophores']))
        {
            $sdfname = "sdf/".(str_replace(' ', '_', $ligname)).".sdf";
            foreach ($mbp['odorophores'] as $moiety => $params)
            {
                $result = [];
                exec("test/moiety_test \"$sdfname\" \"$moiety\"", $result);
                foreach ($result as $line)
                {
                    if (preg_match("/ occurs [1-9][0-9]* times/", $line))
                    {
                        echo "$line\n";
                        foreach ($params as $key => $value)
                        {
                            if ($key == "mcoord")
                            {
                                if (is_array(@$mbp['mcoord'])) $mbp['mcoord'][] = $value;
                                else if (@$mbp['mcoord']) $mbp['mcoord'] = [$mbp['mcoord'], $value];
                                else $mbp['mcoord'] = $value;
                            }
                            else $mbp[$key] = $value;
                        }
                        break;
                    }
                }
            }
        }
    
        if (isset($mbp["size"])) $size = $mbp["size"];
        if (isset($mbp["search"])) $search = $mbp["search"];
        if (isset($mbp["mcoord"]))
        {
            if (!is_array($mbp['mcoord'])) $mbp['mcoord'] = [$mbp['mcoord']];
            foreach ($mbp['mcoord'] as $mc) $mcoord .= "MCOORD $mc\n";
        }
        if (isset($mbp["stcr" ])) $stcr  = "STCR {$mbp["stcr"]}";
        if (isset($mbp["astcr"])) $astcr = "STCR {$mbp["astcr"]}";
        if (isset($mbp["istcr"])) $istcr = "STCR {$mbp["istcr"]}";
        if (isset($mbp["flxr" ])) $flxr  = "FLXR {$mbp["flxr"]}";
        if (isset($mbp["aflxr"])) $aflxr = "FLXR {$mbp["aflxr"]}";
        if (isset($mbp["iflxr"])) $iflxr = "FLXR {$mbp["iflxr"]}";

        if (isset($mbp["stddev"])) $num_std_devs = floatval($mbp["stddev"]);
    
        if (isset($mbp["atomto"]))
        {
            foreach ($mbp["atomto"] as $a2)
            {
                $atomto[] = "ATOMTO $a2";
            }
        }
    }
    
    if ($mbp && isset($mbp["pocket"]))
    {
        $cenres_active = $cenres_inactive = "CEN RES {$mbp["pocket"]}";
    
        if ($mbp && isset($mbp["active_pocket"]))
        {
            $cenres_active = "CEN RES {$mbp["active_pocket"]}";
        }
        if ($mbp && isset($mbp["inactive_pocket"]))
        {
            $cenres_inactive = "CEN RES {$mbp["inactive_pocket"]}";
        }
    }
    else
    {
        if (substr($fam, 0, 2) == "OR")
        {
            $cenres_active = $cenres_inactive = "CEN RES 3.37 5.47 6.55 7.41";
        }
        else if (substr($fam, 0, 4) == "TAAR")
        {
            $cenres_active = $cenres_inactive = "CEN RES 3.32 3.37 5.43 6.48 7.43";
        }
        else die("Unsupported receptor family.\n");
    }
    
    $atomto = implode("\n", $atomto);    
}

$multicall = 0;
function process_dock($metrics_prefix = "", $noclobber = false, $no_sound_if_clashing = false)
{
    global $ligname, $protid, $configf, $dock_retries, $pdbfname, $outfname, $metrics_to_process, $bias_by_energy, $version;
    global $sepyt, $json_file, $do_scwhere, $multicall, $method, $clashcomp, $best_energy;
    $multicall++;
    if ($multicall > 1) $noclobber = true;

    if (!file_exists("tmp")) mkdir("tmp");

    $retvar = 0;
    $best_energy = false;

    $outlines = [];
    if (@$_REQUEST['saved'])
    {
        $outlines = explode("\n", file_get_contents($outfname)); // $_REQUEST['saved']));
    }
    else
    {
        if (!file_exists("bin/vina"))
        {
            exec("wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 -O bin/vina");
            exec("chmod +x bin/vina");
        }
        set_time_limit(300);
        $outlines = [];
        $cmd = "bin/vina --receptor tmp/prot.pdbqt --flex tmp/flex.pdbqt --ligand tmp/lig.pdbqt --center_x 0 --center_y 15 --center_z 0 --size_x 20 --size_y 20 --size_z 20 --exhaustiveness 20 --cpu 1";
        echo "$cmd\n";
        exec($cmd, $outlines, $retvar);
    }

    if (@$_REQUEST['echo']) echo implode("\n", $outlines) . "\n\n";

    if ($retvar)
    {
        dock_failed();
        return 0;
    }
    

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
    $num_poses = 0;

    $outpdb = file_get_contents($pdbfname);

    $posesln = false;
    $affinities = [];
    $mustattribute = true;
    $attribution = "";
    foreach ($outlines as $lno => $ln)
    {
        echo "$ln\n";

        if ($mustattribute && $lno > 2)
        {
            if (false !== strpos($ln, "####")) $mustattribute = false;
            $ln = trim(str_replace('#', '', $ln));
            $attribution .= "REMARK   1 $ln\n";
            continue;
        }

        if (false !== strpos($ln, "affinity"))
        {
            $posesln = true;
            continue;
        }

        if (intval($ln))
        {
            $ln = trim($ln);
            $ln = preg_replace("/\\s+/", " ", $ln);
            $ln = explode(" ", $ln);
            $affinities[] = floatval($ln[1]) * 4.184;
            $num_poses++;
        }
    }

    $outpdb = $attribution.$outpdb;

    if (!$posesln)
    {
        dock_failed();
        return 0;
    }

    $outpdb = str_replace("TER\n", "", $outpdb);
    $outpdb = str_replace("END\n", "", $outpdb);

    $ol = explode("\n", file_get_contents("tmp/lig_out.pdbqt"));
    $fatno = 0;
    foreach($ol as $ln)
    {
        if (substr($ln, 0, 5) != "ATOM ") continue;
        $atno = intval(substr($ln, 7, 4));
        if ($atno < $fatno) break;
        $esym = trim(substr($ln, 13, 2));
        $ln = "HETATM ".(9000+$atno)."  ".str_pad($esym.$atno, 4)."LIG".substr($ln, 20,57).$esym;
        $outpdb .= "$ln\n";
        $fatno = $atno;
    }

    $outpdb .= "\nTER\nEND\n";

    $fp = fopen($outfname, "w");
    if (!$fp) die("Failed to write to $outfname.\n");
    fwrite($fp, $outpdb);
    fclose($fp);

    if ($metrics_prefix && substr($metrics_prefix, -1) != '_') $metrics_prefix .= '_';

    $outdata[$metrics_prefix."POSES"] = $num_poses;
    // $outdata[$metrics_prefix."BENERG"] = $affinities[0];
    $scoring = [];
    $cmd = "bin/score_pdb $outfname | grep \"Total: \"";
    exec($cmd, $scoring);
    echo "$cmd\n";
    $outdata[$metrics_prefix."BENERG"] = floatval(substr($scoring[0], 7));
    echo $scoring[0]."\n";
    
    $outdata['version'] = $version;
    $outdata['method'] = $method;

    // print_r($outdata);

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

    if (function_exists("make_prediction"))
    {
        $oc = count($outdata);
        $prediction = make_prediction($outdata);
        $outdata = array_merge($outdata, $prediction);

        if (stream_isatty(STDOUT) && isset($outdata['Predicted']) && file_exists("predict/soundalert"))
        {
            $play_sound = true;
            if ($no_sound_if_clashing && count($clashcomp)) $play_sound = false;
            $sa = explode(",",file_get_contents("predict/soundalert"));
            if (count($sa) == 2)
            {
                $hfrom = intval($sa[0]);
                $hto = intval($sa[1]);
                $hnow = intval(date("H"));

                if ($hnow < $hfrom || $hnow > $hto) $play_sound = false;
            }

            $hassox = [];
            exec("which sox", $hassox);
            if ($play_sound && count($hassox))
            {
                if (strtolower($outdata['Predicted']) == 'agonist')
                {
                    if (strtolower($outdata['Predicted']) == strtolower($actual)) exec("play success.mp3 &");
                    else if (strtolower($actual) == "(unknown)") exec("play agonist.mp3 &");
                    else exec("play fail.mp3 &");
                }
                else
                {
                    if (strtolower($actual) == 'agonist') exec("play fail.mp3 &");
                    else if (strtolower($actual) == "(unknown)") exec("play non-agonist.mp3 &");
                    else exec("play success.mp3 &");
                }
            }
        }
    }
    $outdata["Actual"] = $actual;

    // Reload to prevent overwriting another process' output.
    if (file_exists($json_file)) $dock_results = json_decode(file_get_contents($json_file), true);
    if ($noclobber)
    {
        $newdata = $outdata;
        $outdata = $dock_results[$protid][$ligname];
        foreach ($newdata as $k => $v) $outdata[$k] = $v;
    }
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

    return $num_poses;
}
