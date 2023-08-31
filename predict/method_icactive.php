<?php

// method_icactive.php
//
// Performs a dock of an odorant inside inactive and internal-contact active conformer PDB files.
//
// Example call syntax:
// php -f predict/method_icactive.php prot=OR1A1 lig=d-limonene
//

// Configurable variables
$dock_metals = false;

chdir(__DIR__);
require("methods_common.php");
chdir(__DIR__);
$binding_pockets = json_decode(file_get_contents("../data/binding_pocket.json"), true);

prepare_outputs();

$size = "5.0 6.0 5.0";

if (substr($fam, 0, 2) == "OR")
{
    $sub = intval(substr($fam, 2));

    if (isset($binding_pockets[$protid]["pocket"]))
    {
        $cenres_active = $cenres_inactive = "CEN RES {$binding_pockets[$protid]["pocket"]}";
        if (isset($binding_pockets[$protid]["size"])) $size = $binding_pockets[$protid]["size"];
    }
    /*else if ($sub >= 50)
    {
        // https://doi.org/10.1101/2022.12.20.520951
        $cenres_active = "CEN RES 4.57 4.60 5.39 45.52 6.59";
        $cenres_inactive = "CEN RES 4.57 4.60 5.39 45.52";
    }*/
    else
    {
        $cenres_active = $cenres_inactive = "CEN RES 3.37 5.47 6.55 7.41";
    }
}
else if (substr($fam, 0, 4) == "TAAR")
{
    die("There is not yet an internal contacts activation app for TAARs.\n");
    $cenres_active = $cenres_inactive = "CEN RES 3.32 3.37 5.43 6.48 7.43";
}
else die("Unsupported receptor family.\n");

$metrics_to_process =
[
    "BENERG" => "BindingEnergy",
    "BENERG.rgn" => "BindingEnergy.rgn",
    "BEST" => "Pose1"
];

function make_prediction($data)
{
    if (isset($data["a_Pose1"]))
    {
        $ae = floatval(@$data['a_BindingEnergy']);
        $ie = floatval(@$data["i_BindingEnergy"]);
        $a1 = floatval( $data['a_Pose1']);
        $i1 = floatval(@$data['i_Pose1']);
        if ($a1 < 0 && $a1 < $i1)
        {
            $data['Predicted'] = 'Agonist';
            $data['DockScore'] = (min($i1, 0) - $a1) / 2;
        }
        else if ($ae < 0 && $ae < $ie)
        {
            $data['Predicted'] = 'Agonist';
            $data['DockScore'] = (min($ie, 0) - $ae) / 2;
        }
        else if ($a1 < 0 && $a1 > $i1)
        {
            $data['Predicted'] = 'Inverse Agonist';
            $data['DockScore'] = (min($i1, 0) - $a1) / 2;
        }
        else
        {
            $data['Predicted'] = 'Non-agonist';
            $data['DockScore'] = 0.0;
        }

        echo "\nResult: " . print_r($data, true) . "\n";
    }

    return $data;
}


chdir(__DIR__);
chdir("..");

$pdbfname_active = str_replace(".upright.pdb", ".active.pdb", $pdbfname);
$paramfname = str_replace(".upright.pdb", ".params", $pdbfname);

if (!file_exists($pdbfname_active) || filemtime($pdbfname_active) < filemtime("bin/ic_activate_or"))
{
    passthru("bin/ic_activate_or $protid");
}

$flex_constraints = "";
if (file_exists($paramfname)) $flex_constraints = file_get_contents($paramfname);


$outfname = "output/$fam/$protid/$protid.$ligname.inactive.dock";

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres_inactive
SIZE $size

EXCL 1 56		# Head, TMR1, and CYT1.

SEARCH BB
POSE 10
ELIM 5000
ITERS 50
PROGRESS

FLEX 1
WET

OUT $outfname
OUTPDB 1 output/$fam/$protid/%p.%l.inactive.model%o.pdb


heredoc;

chdir(__DIR__);
chdir("..");
if (!file_exists("output/$fam")) mkdir("output/$fam");
if (!file_exists("output/$fam/$protid")) mkdir("output/$fam/$protid");

if (!@$_REQUEST["acvonly"]) process_dock("i");


$pdbfname = $pdbfname_active;

$outfname = "output/$fam/$protid/$protid.$ligname.active.dock";

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres_active
SIZE $size

EXCL 1 56		# Head, TMR1, and CYT1.

SEARCH BB
POSE 10
ELIM 5000
$flex_constraints
ITERS 50
PROGRESS

FLEX 1
WET

OUT $outfname
OUTPDB 1 output/$fam/$protid/%p.%l.active.model%o.pdb


heredoc;

process_dock("a");
