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

prepare_outputs();

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
    if ($protid == "OR51E2")
    {
        $pepd = <<<heredoc
LOAD pdbs/OR51/OR51E2.8f76.pdb
CENTER
UPRIGHT
HYDRO
SAVE $pdbfname_active


heredoc;
        $fp = fopen("tmp/8f76.pepd", "wb");
        fwrite($fp, $pepd);
        fclose($fp);
        passthru("bin/pepteditor tmp/8f76.pepd");
    }
    else passthru("bin/ic_activate_or $protid");
}

$flex_constraints = "";
if (file_exists($paramfname)) $flex_constraints = file_get_contents($paramfname);


$outfname = "output/$fam/$protid/$protid.$ligname.inactive.dock";

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres_inactive
SIZE $size
# H2O 5
$mcoord
$atomto
$stcr
$flxr

EXCL 1 56		# Head, TMR1, and CYT1.

SEARCH $search
POSE 10
ELIM 5000
$flex_constraints
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
# H2O 5
$mcoord
$atomto
$stcr
$flxr

EXCL 1 56		# Head, TMR1, and CYT1.

# SOFT 4.53 4.64
# SOFT 5.50 5.33
# SOFT 6.48 6.59

SEARCH $search
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

$poses = (!@$_REQUEST["softonly"]) ? process_dock("a") : 0;

// This has not borne fruit with ionones in OR51E2.
if (false) // !$poses && !@$_REQUEST["nosoft"])
{
    $configf = str_replace("# SOFT ", "SOFT ", $configf);
    process_dock("a");
}
