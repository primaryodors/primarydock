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

    if (isset($data["a_Pose1"]))
    {
        $ae = min(0, floatval(@$data['a_BindingEnergy']));
        $ie = min(0, floatval(@$data["i_BindingEnergy"]));
        $a1 = min(0, floatval( $data['a_Pose1']));
        $i1 = min(0, floatval(@$data['i_Pose1']));
        $a6 = min(0, floatval(@$data["a_BindingEnergy.6"]));
        $i6 = min(0, floatval(@$data["i_BindingEnergy.6"]));

        $ascore = -min(0, $ae - $a6) * $a6 + $a1;
        $iscore = -min(0, $ie - $i6) * $i6 + $i1;

        if ($ascore < 0 && $ascore < $iscore)
        {
            $data['Predicted'] = 'Agonist';
            $data['DockScore'] = (min($i1, 0) - $a1) / 2;
        }
        else if ($i1 < 0 && $iscore < $ascore)
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

$configf = <<<heredoc

PROT $pdbfname_inactive
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
POSE $pose
ELIM 10000
$flex_constraints
ITERS $iter
PROGRESS

FLEX $flxi
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

SEARCH $search
POSE $pose
ELIM 10000
$flex_constraints
ITERS $iter
PROGRESS

FLEX $flex
WET

OUT $outfname
OUTPDB 1 output/$fam/$protid/%p.%l.active.model%o.pdb
APPENDPROT


heredoc;

$poses = process_dock("a");
