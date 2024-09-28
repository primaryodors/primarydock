<?php

global $docker, $dock_metals, $bias_by_energy, $dock_results, $pdbfname, $fam, $do_scwhere, $metrics_to_process, $clashcomp, $best_energy;

// Includes
chdir(__DIR__);
chdir("..");
require("predict/cputemp.php");
require("data/protutils.php");
require("data/odorutils.php");
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
$elim = 0;
$max_simultaneous_docks = 2;	// If running this script as a cron, we recommend setting this to no more than half the number of physical cores.
if (!isset($do_scwhere)) $do_scwhere = false;
$metrics_to_process =
[
  "BENERG" => "BindingEnergy",
  "BENERG.rgn" => "BindingEnergy.rgn",
  "POSES" => "POSES",
  "WCLASH" => "WCLASH"
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

$docker = "pd";
if (@$_REQUEST["docker"]) $docker = $_REQUEST["docker"];
echo "Docker is $docker.\n";

if (@$_REQUEST['simul']) $max_simultaneous_docks = intval($_REQUEST['simul']) ?: 2;

if (@$_REQUEST['next'])
{
    $cmd = "ps -ef | grep -E ':[0-9][0-9] bin/obabel' | grep -v grep";
	exec($cmd, $results);
    if (count($results)) exit;

	$cmd = "ps -ef | grep -E ':[0-9][0-9] (bin/primarydock|bin/pepteditor|bin/ic|obabel)' | grep -v grep";
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

$odor = find_odorant($ligname);
$ligname = $odor['full_name'];


ensure_sdf_exists($ligname);

echo "Beginning prediction of $ligname in $protid...\n\n";
$fam = family_from_protid($protid);
$isomers = check_isomers($ligname);


extract(binding_site($protid));

$outfname = "output.dock";

function bad_docker_code()
{
    die("Unrecognized docker. The choices are pd or vina.\n");
}

function dock_failed($reason = "")
{
    if ($reason) echo "Docking FAILED: $reason.\n";
    else echo "Docking FAILED.\n";

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
            exec("play crash.mp3 &");
        }
    }

    exit;
}

function prepare_receptor($pdbinfname, $lflxr)
{
    global $docker, $protid, $mcoord;
    switch(strtolower($docker))
    {
        case "pd":
        break;

        case "vina":
        if ($mcoord)
        {
            $pieces = explode(" ", $mcoord);
            foreach ($pieces as $k => $pc)
            {
                if (preg_match("/^[0-9]{1,2}[.x][0-9]{1,2}$/", $pc)) $pieces[$k] = "%$pc";
            }
            $mcoord = implode(" ", $pieces);
            $pid = getmypid();
            $lpdbinfname = "tmp/mtl$pid.pdb";
            file_put_contents("tmp/mtl$pid.pepd", <<<heredoc
LOAD $pdbinfname
$mcoord
SAVE $lpdbinfname


heredoc
            );
            $cmd = "bin/pepteditor tmp/mtl$pid.pepd";
            echo "$cmd\n";
            exec($cmd);

            if (!file_exists($lpdbinfname))
            {
                dock_failed();
                exit;
            }
        }
        else
        {
            $lpdbinfname = $pdbinfname;
        }

        $lines = explode("\n", file_get_contents($lpdbinfname));
        $rf = split_pdb_to_rigid_and_flex($protid, $lines, explode(" ", "$lflxr"));
        $fp = fopen("tmp/prot.pdb", "w");
        if (!$fp) die("Failed to write to tmp/prot.pdb.\n");
        fwrite($fp, implode("\n",$rf[0]));
        fclose($fp);
        $fp = fopen("tmp/flex.pdb", "w");
        if (!$fp) die("Failed to write to tmp/flex.pdb.\n");
        fwrite($fp, implode("\n",$rf[1]));
        fclose($fp);

        if ($mcoord)
        {
            unlink($lpdbinfname);
            unlink("tmp/mtl$pid.pepd");
        }

        // Convert to PDBQT format.
        exec("obabel -i pdb tmp/prot.pdb -xr -o pdbqt -O tmp/prot.pdbqt");
        exec("obabel -i pdb tmp/flex.pdb -xs -o pdbqt -O tmp/flex.pdbqt");
        break;

        default:
        bad_docker_code();
    }
}

function prepare_ligand($lig_name)
{
    global $docker;
    switch(strtolower($docker))
    {
        case "pd":
        break;

        case "vina":
        exec("obabel -i sdf \"sdf/$lig_name.sdf\" -o pdbqt -O tmp/lig.pdbqt");
        break;

        default:
        bad_docker_code();
    }
}

function prepare_outputs()
{
    global $ligname, $dock_metals, $protid, $fam, $outfname, $pdbfname, $docker;
    global $binding_pockets, $cenres_active, $cenres_inactive, $size, $search, $num_std_devs;
    global $atomto, $excl, $nodel, $stcr, $flxr, $mcoord, $mbp, $astcr, $istcr, $aflxr, $iflxr;

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
    $odor = find_odorant($ligname);

    $size = "7.5 7.5 7.5";
    $search = "CS";
    $atomto = [];
    $excl = "";
    $nodel = "";
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
                $mol = ($odor && $odor['smiles']) ? $odor['smiles'] : $sdfname;
                $cmd = "test/moiety_test \"$mol\" \"$moiety\"";
                // echo "$cmd\n";
                exec("$cmd", $result);
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
                        break 2;
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
        if (isset($mbp["excl"]))
        {
            $se = explode(" ", trim($mbp["excl"]));
            if (false !== strpos($se[0], ".")) $sr = resno_from_bw($protid, $se[0]);
            else $sr = intval($se[0]);
            if (false !== strpos($se[1], ".")) $er = resno_from_bw($protid, $se[1]);
            else $er = intval($se[1]);
            $excl = "EXCL $sr $er\n";
        }
        if (isset($mbp["nodel"])) $nodel = "NODEL {$mbp["nodel"]}";
    
        if (isset($mbp["atomto"]))
        {
            if (!is_array($mbp["atomto"])) $mbp["atomto"] = [$mbp["atomto"]];
            foreach ($mbp["atomto"] as $a2)
            {
                $atomto[] = "ATOMTO $a2";
            }
        }
    }
    
    if ($mbp && (isset($mbp["pocket"]) || (isset($mbp["active_pocket"]) && isset($mbp["inactive_pocket"]))))
    {
        if (isset($mbp["pocket"]))
        {
            if (is_array($mbp["pocket"]))
            {
                $cenres_inactive = "";
                foreach ($mbp["pocket"] as $p) $cenres_inactive .= "CEN RES $p\n";
                $cenres_active = $cenres_inactive;
            }
            else $cenres_active = $cenres_inactive = "CEN RES {$mbp["pocket"]}";
        }
    
        if ($mbp && isset($mbp["active_pocket"]))
        {
            if (is_array($mbp["active_pocket"]))
            {
                $cenres_active = "";
                foreach ($mbp["active_pocket"] as $p) $cenres_active .= "CEN RES $p\n";
            }
            else $cenres_active = "CEN RES {$mbp["active_pocket"]}";
        }
        if ($mbp && isset($mbp["inactive_pocket"]))
        {
            if (is_array($mbp["inactive_pocket"]))
            {
                $cenres_inactive = "";
                foreach ($mbp["inactive_pocket"] as $p) $cenres_inactive .= "CEN RES $p\n";
            }
            else $cenres_inactive = "CEN RES {$mbp["inactive_pocket"]}";
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
        else if (substr($fam, 0, 4) == "VN1R")
        {
            $cenres_active = $cenres_inactive = "CEN RES 3.32 3.36 3.40 7.38";
        }
        else die("Unsupported receptor family.\n");
    }

    $atomto = implode("\n", $atomto);
}

$multicall = 0;
function process_dock($metrics_prefix = "", $noclobber = false, $no_sound_if_clashing = false)
{
    global $ligname, $isomers, $protid, $configf, $pdbfname, $outfname, $metrics_to_process, $bias_by_energy, $version;
    global $sepyt, $json_file, $do_scwhere, $multicall, $method, $clashcomp, $best_energy, $_REQUEST, $nodel;
    global $cenres, $size, $docker, $mcoord, $atomto, $excl, $stcr, $flxr, $search, $pose, $elim, $flex_constraints, $iter, $flex;
    $multicall++;
    if ($multicall > 1) $noclobber = true;

    if (!file_exists("tmp")) mkdir("tmp");
    $modelfname = preg_replace("/.dock$/", ".model%o.pdb", $outfname);
    $retvar = 0;

    switch(strtolower($docker))
    {
        case "vina":

        if ($metrics_prefix && substr($metrics_prefix, -1) != '_') $metrics_prefix .= '_';
        $cenresno = [];
        $center = "--center_x 0 --center_y 15 --center_z 0";
        $size = "--size_x 20 --size_y 20 --size_z 20";
        if ($cenres)
        {
            foreach (explode(" ", explode("\n", $cenres)[0]) as $bw)
            {
                if (!preg_match("/[A-Z]*[0-9]+([.][0-9]+)?!?/", $bw)) continue;
                $resno = resno_from_bw($protid, $bw);
                if ($resno) $cenresno[] = $resno;
            }

            $cenresno = implode(" ", $cenresno);
            $cmd = "bin/pepteditor predict/center.pepd $pdbfname $cenresno";
            echo "$cmd\n";
            $censz = [];
            exec($cmd, $censz);

            foreach ($censz as $ln)
            {
                if (substr($ln, 0, 5) == "CEN: ")
                {
                    $ln = substr($ln, 5);
                    $pcen = explode(",",str_replace('[','',str_replace(']','',$ln)));
                    $center = "--center_x {$pcen[0]} --center_y {$pcen[1]} --center_z {$pcen[2]}";
                }
                if (substr($ln, 0, 4) == "SZ: ")
                {
                    $ln = substr($ln, 4);
                    $psiz = explode(",",str_replace('[','',str_replace(']','',$ln)));
                    $sx = floatval($psiz[0]) + 5;
                    $sy = floatval($psiz[1]) + 5;
                    $sz = floatval($psiz[2]) + 5;
                    $size = "--size_x $sx --size_y $sy --size_z $sz";
                }
            }
        }

        $best_energy = false;

        if (@$_REQUEST['size']) $size = "--size_x {$_REQUEST['size']} --size_y {$_REQUEST['size']} --size_z {$_REQUEST['size']}";

        $outlines = [];
        if (@$_REQUEST['saved'])
        {
            $outlines = explode("\n", file_get_contents($outfname));
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
            $cmd = "bin/vina --receptor tmp/prot.pdbqt --flex tmp/flex.pdbqt --ligand tmp/lig.pdbqt $center $size --exhaustiveness 20 --cpu 1";
            echo "$cmd\n";
            passthru("$cmd | tee tmp/vina.out", $retvar);
            $outlines = explode("\n", file_get_contents("tmp/vina.out"));
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

        echo "Loading $pdbfname...\n";
        $origpdb = $outpdb = file_get_contents($pdbfname);

        $posesln = false;
        $affinities = [];
        $mustattribute = true;
        $attribution = "";
        foreach ($outlines as $lno => $ln)
        {
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
                $affinities[] = floatval(@$ln[1]) * 4.184;
                $num_poses++;
            }
        }

        $outpdb = $attribution.$outpdb;

        if (!$posesln)
        {
            dock_failed("no affinities");
            return 0;
        }

        $outpdb = str_replace("TER\n", "", $outpdb);
        $outpdb = str_replace("END\n", "", $outpdb);

        $dockpdb = [];
        $ol = explode("\n", file_get_contents("tmp/lig_out.pdbqt"));
        $fatno = 0;
        $lnposeno = 1;
        $dockpdb[$lnposeno] = "";
        foreach($ol as $ln)
        {
            if (substr($ln, 0, 5) != "ATOM ") continue;
            $atno = intval(substr($ln, 7, 4));
            if ($atno < $fatno)
            {
                $lnposeno++;
                $dockpdb[$lnposeno] = "";
            }
            $esym = trim(substr($ln, 13, 2));
            $ln = "HETATM ".(9000+$atno)."  ".str_pad($esym.$atno, 4)."LIG".substr($ln, 20,57).$esym;
            $dockpdb[$lnposeno] .= "$ln\n";
            $fatno = $atno;
        }

        $dockout = <<<heredoc
PDB file: $pdbfname
Ligand: sdf/$ligname.sdf

External dock using third party docker.

heredoc;
        
        foreach ($dockpdb as $poseno => $lpdb)
        {
            $loutpdb = $outpdb."\n".$lpdb;
            $loutpdb .= "\nTER\nEND\n";

            $lmodelfname = str_replace("%o", $poseno, $modelfname);
            $fp = fopen($lmodelfname, "w");
            if (!$fp) die("Failed to write to $lmodelfname.\n");
            fwrite($fp, $loutpdb);
            fclose($fp);

            $outdata[$metrics_prefix."POSES"] = $num_poses;
            $scoring = [];
            $cmd = "bin/score_pdb -n \"$lmodelfname\"";
            exec($cmd, $scoring);
            // echo "$cmd\n";
            $scoring = implode("\n", $scoring);

            $dockout .= <<<heredoc

Pose: $poseno
Node: 0
$scoring

# PDB Data
PDBDAT:
$lpdb
TER
END


heredoc;

        }

        $dockout .= <<<heredoc

$num_poses pose(s) found.

heredoc;

        echo "$dockout\n";

        $dockout .= <<<heredoc

Original PDB:
$origpdb

heredoc;
        $fp = fopen($outfname, "w");
        if (!$fp) die("Failed to write to $outfname.\n");
        fwrite($fp, $dockout);
        fclose($fp);
        $outlines = explode("\n", file_get_contents($outfname));

        break;

        case "pd":
        if ($cenres)
        {
            foreach (explode(" ", explode("\n", $cenres)[0]) as $bw)
            {
                if (!preg_match("/[A-Z]*[0-9]+([.][0-9]+)?!?/", $bw)) continue;
                $resno = resno_from_bw($protid, $bw);
                if ($resno) $cenresno[] = $resno;
            }

            $cavsr = resno_from_bw($protid, "1.28");
            $caver = resno_from_bw($protid, "7.56");

            $fam = family_from_protid($protid);
            $cavfname = "pdbs/$fam/$protid.upright.cvty";
            if (!file_exists($cavfname) || filemtime($cavfname) < filemtime("bin/cavity_search"))
            {
                $cmd = "bin/cavity_search -p pdbs/$fam/$protid.upright.pdb -o $cavfname --ymin 5 --ymax 23 --xzrlim 15 --sr $cavsr --er $caver";
                echo "$cmd\n";
                passthru($cmd);
            }
            $cavfname = "pdbs/$fam/$protid.active.cvty";
            if (!file_exists($cavfname) || filemtime($cavfname) < filemtime("bin/cavity_search"))
            {
                $cmd = "bin/cavity_search -p pdbs/$fam/$protid.active.pdb -o $cavfname --ymin 5 --ymax 23 --xzrlim 15 --sr $cavsr --er $caver";
                echo "$cmd\n";
                passthru($cmd);
            }

            $cenresno = implode(" ", $cenresno);
            $cmd = "bin/pepteditor predict/center.pepd $pdbfname $cenresno";
            echo "$cmd\n";
            $censz = [];
            exec($cmd, $censz);

            switch (family_from_protid($protid))
            {
                case "OR":
                $softness = "1.0";
                break;

                case "TAAR":
                $softness = "0.0";
                break;

                default:
                $softness = "1.0";
            }

            foreach ($censz as $ln)
            {
                if (substr($ln, 0, 4) == "SZ: ")
                {
                    $ln = substr($ln, 4);
                    $psiz = explode(",",str_replace('[','',str_replace(']','',$ln)));
                    $sx = floatval($psiz[0])+2;
                    $sy = floatval($psiz[1])+2;
                    $sz = floatval($psiz[2])+2;
                    $size = "$sx $sy $sz";
                }
            }
        }

        if ($fam == "TAAR") $nodel .= "\nNODEL 45.60 5.38";
        else $nodel .= "\nNODEL 45.51 45.59";

        if ($isomers)
        {
            $iso = [];
            foreach ($isomers as $k => $v) $iso[] = "ISO sdf/".str_replace(' ','_',$v).".sdf";
            $iso = implode("\n", $iso);
            if ($pose < 4*count($isomers)) $pose = 4*count($isomers);
        }
        else $iso = "";

        $excl1 = resno_from_bw($protid, "2.37");
        $soft = /*(($metrics_prefix != "i" && $metrics_prefix != "i_") && $softness) ? "SOFT $softness 1 2 3 4 45 5 6 7" :*/ "";

        $configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf
$iso

$cenres
SIZE $size
# H2O 5
$mcoord
$atomto
$stcr
$flxr

EXCL 1 $excl1		# Head, TMR1, and CYT1.
$excl
ATOMTO 3.39 EXTENT 7.53
BRIDGE 3.39 7.49
STCR 3.39 7.49

SEARCH $search
POSE $pose
ELIM $elim
$flex_constraints
ITERS $iter
PROGRESS
# OUTMC 1
# OUTVDWR 1
# MOVIE

FLEX 1
$soft
$nodel

OUT $outfname
OUTPDB 1 $modelfname
OPEND


heredoc;
        
        chdir(__DIR__);
        chdir("..");
        if (!file_exists("output/$fam")) mkdir("output/$fam");
        if (!file_exists("output/$fam/$protid")) mkdir("output/$fam/$protid");

        $lignospace = str_replace(" ", "", $ligname);
        $prefixfn = $metrics_prefix ? "_$metrics_prefix" : "";
        $cnfname = "tmp/prediction.$protid$prefixfn.$lignospace.config";
        $f = fopen($cnfname, "wb");
        if (!$f) die("File write error. Check tmp folder is write enabled.\n");

        if ($metrics_prefix && substr($metrics_prefix, -1) != '_') $metrics_prefix .= '_';

        fwrite($f, $configf);
        fclose($f);

        $retvar = 0;
        $best_energy = false;

        $outlines = [];
        if (@$_REQUEST['saved'])
        {
            $outlines = explode("\n", file_get_contents($outfname));
        }
        else
        {
            set_time_limit(300);
            $outlines = [];
            $cmd = "bin/primarydock \"$cnfname\"";
            echo "$cmd\n";
            passthru($cmd, $retvar);
            $outlines = explode("\n", file_get_contents($outfname));
        }

        if (@$_REQUEST['echo']) echo implode("\n", $outlines) . "\n\n";
        $num_poses = 0;
        foreach ($outlines as $ln)
        {
            if (false !== strpos($ln, "pose(s) found") || false !== strpos($ln, "poses found") || false !== strpos($ln, "pose found"))
            {
                $num_poses = intval($ln);
            }
        }

        if (!$retvar && $num_poses && !file_exists("tmp/nodelete")) unlink($cnfname);
        else echo "WARNING: Not deleting a temporary config file in tmp/.\n"
            ."Recommend periodically checking the tmp/ folder and manually deleting old files that are no longer necessary.\n";

        break;

        default:
        bad_docker_code();
    }

    if ($retvar)
    {
        dock_failed("return value");
        return 0;
    }    

    $mode = "";
    $lpose = 0;
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

    $posesln = false;
    foreach ($outlines as $ln)
    {
        if (preg_match("/TMR[1-7][.](nseg|center|cseg)[.]clashdir = /", $ln))
        {
            $tmrno = intval(substr($ln, 3, 1));

            $segment = false;
            switch (@explode('.', $ln)[1])
            {
                case "nseg":
                $segment = ($tmrno & 1) ? "exr" : "cyt";
                break;

                case "center":
                $segment = "x.50";
                break;

                case "cseg":
                $segment = ($tmrno & 1) ? "cyt" : "exr";
                break;

                default:
                ;
            }

            if ($segment)
            {
                $coords = explode(" = ", $ln)[1];
                $coords = str_replace('[', '', $coords);
                $coords = str_replace(']', '', $coords);
                $coords = explode(',', $coords);
                foreach ($coords as $k => $c) $coords[$k] = floatval($c);

                $clashcomp[$tmrno][$segment] = $coords;
                continue;
            }
        }

        if (false !== strpos($ln, "pose(s) found") || false !== strpos($ln, "poses found") || false !== strpos($ln, "pose found"))
        {
            $num_poses = intval($ln);
            $posesln = true;
        }
    }

    if (!$posesln)
    {
        dock_failed("no poses");
        return 0;
    }

    if (!$num_poses) $outdata[$metrics_prefix."POSES"] = $num_poses;
    else foreach ($outlines as $ln)
    {
        $coldiv = explode(":", $ln);

        if (trim($ln) == "TER")
        {
            $lpose = false;
            $node = -1;
            continue;
        }

        if (count($coldiv) == 2)
        {
            if ($coldiv[0] == "Pose")
            {
                $lpose = intval($coldiv[1]);
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
                if ($e < 1) $e = 0; // 1.0 / abs(log(abs($e)));
                $weight[$lpose] = $e;
            }
        }
    }

    foreach ($outlines as $lno => $ln)
    {
        $coldiv = explode(":", $ln);

        if (trim($ln) == "TER")
        {
            $lpose = false;
            $node = -1;
            continue;
        }

        if (count($coldiv) == 2)
        {
            if ($coldiv[0] == "Pose")
            {
                $lpose = intval($coldiv[1]);
                $node = -1;
                continue;
            }
            else if ($coldiv[0] == "Node")
            {
                $node = intval($coldiv[1]);
                continue;
            }
            else if ($coldiv[0] == 'BENERG' || $coldiv[0] == 'MC' || $coldiv[0] == 'vdWRPL')
            {
                $mode = $coldiv[0];
                continue;
            }
            else if ($lpose && $node>=0)
            {
                if (preg_match("/[A-Z][a-z]{2}[0-9]{1,4}/", $coldiv[0]))
                {
                    $resno = intval(substr($coldiv[0], 3));
                    $region = intval(bw_from_resno($protid, $resno));
                    if (isset($metrics_to_process["$mode.rgn"]))
                    {
                        $wmode = str_replace(".rgn", ".$region", $metrics_to_process["$mode.rgn"]);
                        if (!isset($outdata[$metrics_prefix.$wmode])) $outdata[$metrics_prefix.$wmode] = 0.0;
                        $outdata[$metrics_prefix.$wmode] += floatval($coldiv[1]) * $weight[$lpose];
                        if (!@$pose_node_has_weight[$wmode][$lpose][$node])
                        {
                            if (!isset($outdqty[$metrics_prefix.$wmode])) $outdqty[$metrics_prefix.$wmode] = $weight[$lpose];
                            else $outdqty[$metrics_prefix.$wmode] += $weight[$lpose];
                            $pose_node_has_weight[$wmode][$lpose][$node] = true;
                        }
                    }
                    continue;
                }
                if ($coldiv[0] == "Total")
                {
                    if (isset($metrics_to_process[$mode]) && floatval($coldiv[1]) < 0)
                    {
                        $wmode = $metrics_to_process[$mode];
                        if (!isset($outdata[$metrics_prefix.$wmode])) $outdata[$metrics_prefix.$wmode] = 0.0;
                        $outdata[$metrics_prefix.$wmode] += floatval($coldiv[1]);
                        if (!isset($outdqty[$metrics_prefix.$wmode])) $outdqty[$metrics_prefix.$wmode] = 1;
                        else $outdqty[$metrics_prefix.$wmode]++;
                    }

                    if ($mode == "BENERG")
                    {
                        $benerg = floatval($coldiv[1]);
                        if (false===$best_energy || $benerg < $best_energy) $best_energy = $benerg;
                    }

                    if ($mode == "BENERG" && isset($metrics_to_process["BEST"]))
                    {
                        $wmode = $metrics_to_process["BEST"];
                        if (!isset($outdata[$metrics_prefix.$wmode])) $outdata[$metrics_prefix.$wmode] = floatval($coldiv[1]);
                    }
                    continue;
                }
                else if ($coldiv[0] == "Ligand polar satisfaction")
                {
                    $mode = "POLSAT";
                    if (isset($metrics_to_process[$mode]))
                    {
                        $wmode = $metrics_to_process[$mode];
                        if (!isset($outdata[$metrics_prefix.$wmode])) $outdata[$metrics_prefix.$wmode] = 0.0;
                        $outdata[$metrics_prefix.$wmode] += floatval($coldiv[1]) * $weight[$lpose];
                        if (!isset($outdqty[$metrics_prefix.$wmode])) $outdqty[$metrics_prefix.$wmode] = $weight[$lpose];
                        else $outdqty[$metrics_prefix.$wmode] += $weight[$lpose];
                    }
                    continue;
                }
                else if ($coldiv[0] == "A100 score")
                {
                    $mode = "A100";
                    if (isset($metrics_to_process[$mode]))
                    {
                        $wmode = $metrics_to_process[$mode];
                        if (!isset($outdata[$metrics_prefix.$wmode])) $outdata[$metrics_prefix.$wmode] = 0.0;
                        $outdata[$metrics_prefix.$wmode] += floatval($coldiv[1]) * $weight[$lpose];
                        if (!isset($outdqty[$metrics_prefix.$wmode])) $outdqty[$metrics_prefix.$wmode] = $weight[$lpose];
                        else $outdqty[$metrics_prefix.$wmode] += $weight[$lpose];
                    }
                    continue;
                }
                else if ($coldiv[0] == "Ligand occlusion")
                {
                    $mode = "occlusion";
                    if (isset($metrics_to_process[$mode]))
                    {
                        $wmode = $metrics_to_process[$mode];
                        if (!isset($outdata[$metrics_prefix.$wmode])) $outdata[$metrics_prefix.$wmode] = 0.0;
                        $outdata[$metrics_prefix.$wmode] += floatval($coldiv[1]) * $weight[$lpose];
                        if (!isset($outdqty[$metrics_prefix.$wmode])) $outdqty[$metrics_prefix.$wmode] = $weight[$lpose];
                        else $outdqty[$metrics_prefix.$wmode] += $weight[$lpose];
                    }
                    continue;
                }
                else if ($coldiv[0] == "Protein clashes")
                {
                    $mode = "PCLASH";
                    if (isset($metrics_to_process[$mode]))
                    {
                        $wmode = $metrics_to_process[$mode];
                        if (!isset($outdata[$metrics_prefix.$wmode])) $outdata[$metrics_prefix.$wmode] = 0.0;
                        $outdata[$metrics_prefix.$wmode] += floatval($coldiv[1]) * $weight[$lpose];
                        if (!isset($outdqty[$metrics_prefix.$wmode])) $outdqty[$metrics_prefix.$wmode] = $weight[$lpose];
                        else $outdqty[$metrics_prefix.$wmode] += $weight[$lpose];
                    }
                    continue;
                }
                else if (strpos($coldiv[0], " active theta"))
                {
                    $morceaux = explode(' ', $coldiv[0]);
                    $mode = "ACVTH.{$morceaux[0]}";
                    if (isset($metrics_to_process[$mode]))
                    {
                        $wmode = $metrics_to_process[$mode];
                        if (!isset($outdata[$metrics_prefix.$wmode])) $outdata[$metrics_prefix.$wmode] = 0.0;
                        $outdata[$metrics_prefix.$wmode] += floatval($coldiv[1]) * $weight[$lpose];
                        if (!isset($outdqty[$metrics_prefix.$wmode])) $outdqty[$metrics_prefix.$wmode] = $weight[$lpose];
                        else $outdqty[$metrics_prefix.$wmode] += $weight[$lpose];
                    }
                    continue;
                }
            }
        }

        if (false !== strpos($ln, "pose(s) found") || false !== strpos($ln, "poses found") || false !== strpos($ln, "pose found"))
        {
            $mode = "POSES";
            if (isset($metrics_to_process[$mode]))
            {
                $wmode = $metrics_to_process[$mode];
                $outdata[$metrics_prefix.$wmode] = intval($ln);
            }
            continue;
        }
        else if (false !== strpos($ln, "Best worst clash"))
        {
            $mode = "WCLASH";
            if (isset($metrics_to_process[$mode]))
            {
                $wmode = $metrics_to_process[$mode];
                $pettia = explode(':',$ln);
                $outdata[$metrics_prefix.$wmode] = floatval(@$pettia[1]);
            }
            continue;
        }
        else if (false !== strpos($ln, "Repeatability"))
        {
            $mode = "Repeatability";
            if (isset($metrics_to_process[$mode]))
            {
                $wmode = $metrics_to_process[$mode];
                $outdata[$metrics_prefix.$wmode] = intval($ln);
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

                    if (!isset($outdata[$metrics_prefix.$wmode])) $outdata[$metrics_prefix.$wmode] = [0,0,0];
                    for ($x=0; $x<3; $x++) $outdata[$metrics_prefix.$wmode][$x] += floatval(substr($ln, 29+8*x, 8)) * $weight[$lpose];

                    if (!isset($outdqty[$metrics_prefix.$wmode])) $outdqty[$metrics_prefix.$wmode] = $weight[$lpose];
                    else $outdqty[$metrics_prefix.$wmode] += $weight[$lpose];
                }
                break;
    
                default:
                ;
            }
        }
    }

    $outdata['version'] = $version;
    $outdata['method'] = $method;
    $outdata['docker'] = $docker;

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
                if (strtolower($outdata['Predicted']) == strtolower($actual)) exec("play success.mp3 &");
                else if (strtolower($actual) == "non-agonist" && strtolower($outdata['Predicted']) == 'inverse agonist') exec("play success.mp3 &");
                else if (strtolower($actual) == "(unknown)")
                {
                    if (strtolower($outdata['Predicted']) == 'agonist') exec("play agonist.mp3 &");
                    else exec("play non-agonist.mp3 &");
                }
                else
                {
                    $d = dir(getcwd());
                    $failsound = [];
                    while (false !== ($entry = $d->read()))
                    {
                        if (substr(strtolower($entry), 0, 4) == "fail")
                        {
                            $ext = substr($entry, strrpos($entry, '.'));
                            if ($ext == '.mp3' || $ext == '.wav') $failsound[] = $entry;
                        }
                    }
                    $fn = count($failsound) ? $failsound[rand(0, count($failsound))] : "fail.mp3";
                    exec("play $fn &");
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
