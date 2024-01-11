<?php

// method_fygactive.php
//
// Performs a dock of an odorant inside inactive and FYG-motif active conformer PDB files.
//
// Example call syntax:
// php -f predict/method_fygactive.php prot=OR1A1 lig=d-limonene
//

chdir(__DIR__);
require_once("methods_common.php");
chdir(__DIR__);
require_once("../data/protutils.php");
chdir(__DIR__);
require_once("template.php");
chdir(__DIR__);

// Configurable variables
$flex = 1;                      // Flexion (0 or 1).
$pose = 15;
$iter = 30;
$elim = 1e3;                    // Energy limit for poses. (Not the tailor/spy from the space station.)
$num_std_devs = 2.0;            // How many standard deviations to move the helices for active clash compensation.

prepare_outputs();

$metrics_to_process["BEST"] = "Pose1";

function make_prediction($data)
{
    global $protid, $ligname;

    if (isset($data["a_BENERG"]) || isset($data["a_BindingEnergy"]))
    {
        // TODO: For ascore and iscore each, require contact between ligand and TMR6 and
        // contact between ligand and any of TMR3, TMR4, EXR2, TMR5, TMR7. If either condition
        // not met, zero out that score. For FYG activation receptors with no rock6, also require
        // that the ligand extend far enough in the -Y direction impinge on the vdW space of at
        // least one side chain displaced by FYG activation.
        $ascore = min(0, floatval(@$data['a_BENERG']), floatval(@$data['a_BindingEnergy']));
        $iscore = min(0, floatval(@$data['i_BENERG']), floatval(@$data['i_BindingEnergy']));

        if ($ascore < 0 && $ascore < $iscore)
        {
            $data['Predicted'] = 'Agonist';
            $data['DockScore'] = (min($iscore, 0) - $ascore);
        }
        else if ($iscore < 0 && $iscore < $ascore)
        {
            $data['Predicted'] = 'Inverse Agonist';
            $data['DockScore'] = (min($iscore, 0) - $ascore);
        }
        else
        {
            $data['Predicted'] = 'Non-agonist';
            $data['DockScore'] = 0.0;
        }

        echo "\nProtein: $protid\nLigand: $ligname";
        echo "\nResult: " . print_r($data, true) . "\n";
    }
    else if (isset($data["a_POSES"]) && !$data["a_POSES"])
    {
        $data['Predicted'] = 'Non-agonist';
        $data['DockScore'] = 0.0;

        echo "\nProtein: $protid\nLigand: $ligname";
        echo "\nResult: " . print_r($data, true) . "\n";
    }

    return $data;
}

chdir(__DIR__);
chdir("..");

$pdbfname_active = str_replace(".upright.pdb", ".active.pdb", $pdbfname);
$paramfname = str_replace(".upright.pdb", ".params", $pdbfname);
$template = [];

if (!file_exists($pdbfname_active) || filemtime($pdbfname_active) < filemtime("bin/fyg_activate_or"))
{
    if (check_already_fyg_activating($protid))
    {
        echo "Waiting on activation...";
        sleep(30);
        while (check_already_fyg_activating($protid))
        {
            echo ".";
            sleep(10);
        }
        echo "\n";

        if (!file_exists($pdbfname_active)) dock_failed();
    }
    else do_templated_activation();
}

$flex_constraints = "";
if (file_exists($paramfname)) $flex_constraints = file_get_contents($paramfname);

$fam = family_from_protid($protid);
$outfname = "output/$fam/$protid/$protid.$ligname.inactive.dock";
$cenres = substr($cenres_inactive, 8);

prepare_receptor($pdbfname, "$flxr $iflxr");

// Convert ligand as well.
prepare_ligand($ligname);

chdir(__DIR__);
chdir("..");
if (!file_exists("output/$fam")) mkdir("output/$fam");
if (!file_exists("output/$fam/$protid")) mkdir("output/$fam/$protid");

if (!@$_REQUEST["acvonly"]) process_dock("i");


$pdbfname = $pdbfname_active;
$outfname = "output/$fam/$protid/$protid.$ligname.active.dock";
$cenres = substr($cenres_active, 8);

prepare_receptor($pdbfname, "$flxr $aflxr");

$poses = process_dock("a");

exit;     // TODO: Bring back all that stuff in methods_common get from scorpdbee what used to get from primarysuck.

if ((!$poses || $best_energy >= 0) && count($clashcomp) && $num_std_devs)
{
    if (!count($template)) build_template();

    // Ensure template contains standard deviations. If not, average the ones from the templates that do.
    foreach ($template as $hxno => $metrics)
    {
        foreach ($metrics as $metric => $xyz)
        {
            if (!isset($xyz["sigma"]))
            {
                $averages = [];
                $counts = [];

                foreach ($cryoem as $lrecep => $lhelixes)
                {
                    foreach ($lhelixes as $lhelix => $lmetrics)
                    {
                        foreach ($lmetrics as $lmetric => $lxyz)
                        {
                            if (isset($lxyz["sigma"]))
                            {
                                if (!isset($averages[$lhelix][$lmetric])) $averages[$lhelix][$lmetric] = floatval($lxyz["sigma"]);
                                else $averages[$lhelix][$lmetric] += floatval($lxyz["sigma"]);

                                if (!isset($counts[$lhelix][$lmetric])) $counts[$lhelix][$lmetric] = 1;
                                else $counts[$lhelix][$lmetric]++;
                            }
                        }
                    }
                }

                foreach ($averages as $lhelix => $lmetrics)
                {
                    foreach ($lmetrics as $lmetric => $total)
                    {
                        if (!isset($template[$lhelix][$lmetric]["sigma"]))
                        {
                            if ($counts[$lhelix][$lmetric]) $total /= $counts[$lhelix][$lmetric];
                            $template[$lhelix][$lmetric]["sigma"] = $total;
                        }
                    }
                }

                break;
            }
        }
    }

    foreach ($clashcomp as $tmrno => $cc)
    {
        foreach ($cc as $segment => $xyz)
        {
            // Check each clash compensation is within the SD limit. If not, reduce the magnitude to the deviation limit.
            $sigma = floatval(@$template[$tmrno][$segment]["sigma"]);
            $rlimit = $sigma*$num_std_devs;
            if (!$sigma) continue;
            $tmparr = []; foreach (['x','y','z'] as $idx => $var) $tmparr[$var] = floatval($xyz[$idx]); extract($tmparr);
            $r = sqrt($x*$x+$y*$y+$z*$z);

            if ($r > $rlimit)
            {
                $divisor = $rlimit / $r;
                $x *= $divisor; $y *= $divisor; $z *= $divisor;
            }

            // Apply the corrected compensation.
            foreach (['x','y','z'] as $var)
            {
                $$var = round($$var, 3);
                $template[$tmrno][$segment][$var] = floatval($template[$tmrno][$segment][$var]) + $$var;
                echo "Compensating TMR$tmrno.$segment.$var to {$template[$tmrno][$segment][$var]}...\n";
            }
        }
    }

    // Create a temporary custom output PDB and retry the active dock.
    $args = "$protid";
    foreach ($template as $hxno => $metrics)
    {
        foreach ($metrics as $metric => $dimensions)
        {
            if ($hxno == 6)
            {
                if ($metric == "cyt" && ($has_fyg || $has_rock6)) continue;
                else if ($metric == "exr" && $has_rock6) continue;
            }

            $cmdarg = "--" . substr($metric, 0, 1) . $hxno;
            $args .= " $cmdarg {$dimensions['x']} {$dimensions['y']} {$dimensions['z']}";
        }
    }

    $tmpoutpdb = "tmp/$protid.".getmypid().".pdb";
    $args .= " -o $tmpoutpdb";

    $cmd = "bin/fyg_activate_or $args";
    echo "$cmd\n";
    passthru($cmd);


    $configf = str_replace($pdbfname, $tmpoutpdb, $configf);
    $poses = process_dock("ad");

    // Delete the tmp PDB.
    unlink($tmpoutpdb);
}

































