<?php

// method_fygactive.php
//
// Performs a dock of an odorant inside inactive and FYG-motif active conformer PDB files.
//
// Example call syntax:
// php -f predict/method_fygactive.php prot=OR1A1 lig=d-limonene
//

// Configurable variables
$flex = 1;                      // Flexion (0 or 1) for active dock.
$flxi = 1;                      // Flexion for inactive dock.
$pose = 20;
$iter = 40;

chdir(__DIR__);
require_once("methods_common.php");
chdir(__DIR__);
require_once("../data/protutils.php");
chdir(__DIR__);

prepare_outputs();

$metrics_to_process =
[
    "BENERG" => "BindingEnergy",
    "BENERG.rgn" => "BindingEnergy.rgn",
    "BEST" => "Pose1"
];

function make_prediction($data)
{
    global $protid, $ligname;

    if (isset($data["a_BENERG"]))
    {
        $ascore = min(0, floatval( $data['a_BENERG']));
        $iscore = min(0, floatval(@$data['i_BENERG']));

        if ($ascore < 0 && $ascore < $iscore)
        {
            $data['Predicted'] = 'Agonist';
            $data['DockScore'] = (min($iscore, 0) - $ascore) / 2;
        }
        else if ($iscore < 0 && $iscore < $ascore)
        {
            $data['Predicted'] = 'Inverse Agonist';
            $data['DockScore'] = (min($i1, 0) - $a1) / 2;
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

if (!file_exists($pdbfname_active))
{
    exec("php -f predict/cryoem_motions.php");
    exec("bin/pepteditor data/OR51.pepd");
    exec("bin/pepteditor data/OR52.pepd");
}

$flex_constraints = "";
if (file_exists($paramfname)) $flex_constraints = file_get_contents($paramfname);

$fam = family_from_protid($protid);
$outfname = "output/$fam/$protid/$protid.$ligname.inactive.dock";

// Filter out everything except ATOM records.
exec("cat $pdbfname_inactive | grep ATOM > tmp/prot.pdb");

// Convert to PDBQT format.
exec("obabel -i pdb tmp/prot.pdb -xr -o pdbqt -O tmp/prot.pdbqt");

// Convert ligand as well.
exec("obabel -i sdf sdf/$ligname.sdf -o pdbqt -O tmp/lig.pdbqt");

chdir(__DIR__);
chdir("..");
if (!file_exists("output/$fam")) mkdir("output/$fam");
if (!file_exists("output/$fam/$protid")) mkdir("output/$fam/$protid");

if (!@$_REQUEST["acvonly"]) process_dock("i");


$pdbfname = $pdbfname_active;
$outfname = "output/$fam/$protid/$protid.$ligname.active.dock";

// Filter out everything except ATOM records.
exec("cat $pdbfname_active | grep ATOM > tmp/prot.pdb");

// Convert to PDBQT format.
exec("obabel -i pdb tmp/prot.pdb -xr -o pdbqt -O tmp/prot.pdbqt");

$poses = process_dock("a");
