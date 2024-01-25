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

    if (isset($data["a_BENERG"]) || isset($data["a_BindingEnergy"]))
    {
        $ascore = min(0, floatval(@$data['a_BENERG']));
        $iscore = min(0, floatval(@$data['i_BENERG']));

        if ($ascore > $iscore)
        {
            $ascore = min(0, floatval(@$data['a_BindingEnergy']));
            $iscore = min(0, floatval(@$data['i_BindingEnergy']));
        }

        $aeff = isset($data["a_POSES"]) ? (floatval($data["a_POSES"]) / $pose) : 1;
        $ieff = isset($data["i_POSES"]) ? (floatval($data["i_POSES"]) / $pose) : 1;

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

$pdbfname_inactive = str_replace(".upright.pdb", ".apo.pdb", $pdbfname);
$pdbfname_active = str_replace(".upright.pdb", ".bound.pdb", $pdbfname);
$paramfname = str_replace(".upright.pdb", ".params", $pdbfname);

if (!file_exists($pdbfname_inactive) && file_exists($pdbfname_active))
    $pdbfname_inactive = $pdbfname;

if (!file_exists($pdbfname_active) && (substr($protid, 0, 4) == "OR51" || substr($protid, 0, 4) == "OR52"))
{
    exec("php -f predict/cryoem_motions.php");
    exec("bin/pepteditor data/OR51.pepd");
    exec("bin/pepteditor data/OR52.pepd");
}

if (!file_exists($pdbfname_active)) die("No bound model.\n");

$flex_constraints = "";
if (file_exists($paramfname)) $flex_constraints = file_get_contents($paramfname);

$fam = family_from_protid($protid);
$pdbfname = $pdbfname_inactive;
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

// TODO: Separate dynamic_clash_compensation() into FYG-activation and direct-model editions.
/* if ((!$poses || $best_energy >= 0) && count($clashcomp) && $num_std_devs)
{
    dynamic_clash_compensation();

    $pdbfname = $tmpoutpdb;
    $outfname = "output/$fam/$protid/$protid.$ligname.dynamic.dock";
    prepare_receptor($pdbfname, "$flxr $aflxr");
    $poses = process_dock("ad");

    // Delete the tmp PDB.
    unlink($tmpoutpdb);
} */
