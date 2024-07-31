<?php

// method_directmdl.php
//
// Performs a dock of an odorant inside inactive and active conformer PDB files.
// The active PDB may be an experimental structure from a crystallography or cryo-EM
// measurement, or it may be a homology model based on experimental models.
//
// Example call syntax:
// php -f predict/method_directmdl.php prot=OR1A1 lig=d-limonene
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
$pose = 10;
$iter = 30;
$elim = 1e3;                    // Energy limit for poses. (Not the tailor/spy from the space station.)
$num_std_devs = 2.0;            // How many standard deviations to move the helices for active clash compensation.

prepare_outputs();

$metrics_to_process["BEST"] = "Pose1";
$metrics_to_process["A100"] = "A100";

function make_prediction($data)
{
    global $protid, $ligname, $pose;

    if (isset($data["a_Pose1"]) || isset($data["a_BindingEnergy"]))
    {
        $ascore = min(0, floatval(@$data['a_Pose1']));
        $iscore = min(0, floatval(@$data['i_Pose1']));

        $aa100 = floatval(@$data['a_A100']);
        $ia100 = floatval(@$data['i_A100']);

        $aeff = isset($data["a_POSES"]) ? (floatval($data["a_POSES"]) / $pose) : 1;
        $ieff = isset($data["i_POSES"]) ? (floatval($data["i_POSES"]) / $pose) : 1;

        if ($ascore < 0 && $ascore < $iscore)
        {
            $data['Predicted'] = 'Agonist';
            $dock_score = -$ascore*$aeff;
            if ($aa100 && $ia100) $dock_score *= max(1, ($aa100 - $ia100) / 50);
            $data['DockScore'] = $dock_score;
        }
        else if ($iscore < 0 && $iscore < $ascore)
        {
            $data['Predicted'] = 'Inverse Agonist';
            $dock_score = $iscore*$aeff;
            if ($aa100 && $ia100) $dock_score *= max(1, ($aa100 - $ia100) / 50);
            $data['DockScore'] = $dock_score;
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

$pdbfname_inactive = $pdbfname;
$pdbfname_active = str_replace(".upright.pdb", ".active.pdb", $pdbfname);
$paramfname = str_replace(".upright.pdb", ".params", $pdbfname);

if (!file_exists($pdbfname_active)) die("No active model.\n");

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
