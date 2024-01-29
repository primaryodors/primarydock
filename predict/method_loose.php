<?php

// method_loose.php
//
// Performs a dock of an odorant inside an inactive model, then isolates the ligand binding residues,
// allowing them to move freely, before constructing an active model that approximates the residues'
// freely conformed positions.
//
// Example call syntax:
// php -f predict/method_loose.php prot=OR1A1 lig=geraniol
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

$pdbfname_inactive = $pdbfname;
$pdbfname_active = "tmp/$protid.$ligname.loose.pdb";

$fam = family_from_protid($protid);
$pdbfname = $pdbfname_inactive;
$outfname = "output/$fam/$protid/$protid.$ligname.inactive.dock";
$cenres = substr($cenres_inactive, 8);

chdir(__DIR__);
chdir("..");
if (!file_exists("output/$fam")) mkdir("output/$fam");
if (!file_exists("output/$fam/$protid")) mkdir("output/$fam/$protid");

prepare_receptor($pdbfname, "$flxr $iflxr");
prepare_ligand($ligname);
if (!@$_REQUEST["acvonly"]) process_dock("i");

$cmd = "bin/loose -o '$pdbfname_active' '$outfname'";
echo "$cmd\n";
exec($cmd);

if (!file_exists($pdbfname_active)) dock_failed("The bin/loose tool did not generate an active model.");

$pdbfname = $pdbfname_active;
$outfname = "output/$fam/$protid/$protid.$ligname.active.dock";
$cenres = substr($cenres_active, 8);

prepare_receptor($pdbfname, "$flxr $aflxr");
$poses = process_dock("a");
