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
    global $protid, $ligname, $pose;

    if (isset($data["a_Pose1"]) || isset($data["a_BindingEnergy"]))
    {
        $ascore = min(0, floatval(@$data['a_Pose1']));
        $iscore = min(0, floatval(@$data['i_Pose1']));

        if ($ascore >= $iscore)
        {
            $ascore = min(0, floatval(@$data['a_BindingEnergy']));
            $iscore = min(0, floatval(@$data['i_BindingEnergy']));
        }

        $aeff = isset($data["a_POSES"]) ? (floatval($data["a_POSES"]) / $pose) : 1;
        $ieff = isset($data["i_POSES"]) ? (floatval($data["i_POSES"]) / $pose) : 1;

        if (floatval(@$data['a_Pose1']) < floatval(@$data['i_Pose1']))
        {
            $ascore = floatval(@$data['a_Pose1']);
            $iscore = floatval(@$data['i_Pose1']);
        }

        // For ascore and iscore each, require contact between ligand and TMR6 and
        // contact between ligand and any of TMR3, TMR4, EXR2, TMR5, TMR7. If either condition
        // not met, zero out that score.
        if (!@$data["a_BindingEnergy.6"]) $ascore = 0;
        if (!@$data["a_BindingEnergy.3"] && !@$data["a_BindingEnergy.4"] && !@$data["a_BindingEnergy.45"]
            && !@$data["a_BindingEnergy.5"] && !@$data["a_BindingEnergy.7"]) $iscore = 0;
        if (!@$data["i_BindingEnergy.6"]) $iscore = 0;
        if (!@$data["i_BindingEnergy.3"] && !@$data["i_BindingEnergy.4"] && !@$data["i_BindingEnergy.45"]
            && !@$data["i_BindingEnergy.5"] && !@$data["i_BindingEnergy.7"]) $iscore = 0;

        // TODO: For FYG activation receptors with no rock6, also require
        // that the ligand extend far enough in the -Y direction impinge on the vdW space of at
        // least one side chain displaced by FYG activation.

        if ($ascore < 0 && $ascore < $iscore)
        {
            $data['Predicted'] = 'Agonist';
            $data['DockScore'] = (min($iscore, 0) - $ascore)*$aeff;
        }
        else if ($iscore < 0 && $iscore < $ascore)
        {
            $data['Predicted'] = 'Inverse Agonist';
            $data['DockScore'] = (min($iscore, 0) - $ascore)*$ieff;
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

if (false && (!$poses || $best_energy >= 0) && count($clashcomp) && $num_std_devs)
{
    dynamic_clash_compensation();

    $pdbfname = $tmpoutpdb;
    $outfname = "output/$fam/$protid/$protid.$ligname.dynamic.dock";
    prepare_receptor($pdbfname, "$flxr $aflxr");
    $poses = process_dock("ad");

    // Delete the tmp PDB.
    unlink($tmpoutpdb);
}
