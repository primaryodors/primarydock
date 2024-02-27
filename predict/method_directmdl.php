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
$elim = 0;                      // Energy limit for poses.
$num_std_devs = 2.0;            // How many standard deviations to move the helices for active clash compensation.

prepare_outputs();

$metrics_to_process["BEST"] = "Pose1";

$aposes = 0;
function make_prediction($data)
{
    global $protid, $ligname, $pose, $aposes;

    if (isset($data["a_Pose1"]) || isset($data["a_BindingEnergy"]))
    {
        $ascore = min(0, floatval(@$data['a_BindingEnergy']));
        $iscore = min(0, floatval(@$data['i_BindingEnergy']));
        $score345 = floatval(@$data['a_BindingEnergy.3'])
            + floatval(@$data['a_BindingEnergy.4'])
            + floatval(@$data['a_BindingEnergy.45'])
            + floatval(@$data['a_BindingEnergy.5']);

        $aeff = isset($data["a_POSES"]) ? (floatval($data["a_POSES"]) / $pose) : 1;
        $ieff = isset($data["i_POSES"]) ? (floatval($data["i_POSES"]) / $pose) : 1;

        if ($ascore < 0 && $score345 <= 0)
        {
            $data['Predicted'] = 'Agonist';
            $data['DockScore'] = -$ascore*$aeff;
        }
        else if ($iscore < 0) // && $iscore < $ascore)
        {
            $data['Predicted'] = 'Inverse Agonist';
            $data['DockScore'] = -$ascore*$ieff;
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

$pdbfname_inactive = str_replace(".upright.pdb", ".apo.pdb", $pdbfname);
$pdbfname_active = str_replace(".upright.pdb", ".bound.pdb", $pdbfname);
$paramfname = str_replace(".upright.pdb", ".params", $pdbfname);

if (!file_exists($pdbfname_inactive))              // If no apo model, just use upright.
    $pdbfname_inactive = $pdbfname;

if (!file_exists($pdbfname_active))
{
    if (substr($protid, 0, 4) == "OR51" || substr($protid, 0, 4) == "OR52")
    {
        exec("php -f predict/cryoem_motions.php");
        if (substr($protid, 0, 4) == "OR51") exec("bin/pepteditor data/OR51.pepd");
        if (substr($protid, 0, 4) == "OR52") exec("bin/pepteditor data/OR52.pepd");
    }
    else $pdbfname_active = str_replace(".upright.pdb", ".evolved.pdb", $pdbfname);
}

if (!file_exists($pdbfname_active)) die("No bound model.\n");

$flex_constraints = "";
if (file_exists($paramfname)) $flex_constraints = file_get_contents($paramfname);

$fam = family_from_protid($protid);
$pdbfname = $pdbfname_active;
$outfname = "output/$fam/$protid/$protid.$ligname.active.dock";
$cenres = substr($cenres_active, 8);

prepare_receptor($pdbfname, "$flxr $aflxr");

// Convert ligand as well.
prepare_ligand($ligname);

chdir(__DIR__);
chdir("..");
if (!file_exists("output/$fam")) mkdir("output/$fam");
if (!file_exists("output/$fam/$protid")) mkdir("output/$fam/$protid");

if (!@$_REQUEST["iacvonly"]) process_dock("a");

if (!$aposes)
{
    $pdbfname = $pdbfname_inactive;
    $outfname = "output/$fam/$protid/$protid.$ligname.inactive.dock";
    $cenres = substr($cenres_inactive, 8);

    prepare_receptor($pdbfname, "$flxr $iflxr");

    $poses = process_dock("i");
}
