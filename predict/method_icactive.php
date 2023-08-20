<?php

// method_icactive.php
//
// Performs a dock of an odorant inside a default inactive-conformer PDB file.
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

if (substr($fam, 0, 2) == "OR")
{
    $sub = intval(substr($fam, 2));

    if ($sub >= 50)
    {
        // https://doi.org/10.1101/2022.12.20.520951
        $cenres = "CEN RES 4.57 4.60 5.39 45.52";
    }
    else
    {
        $cenres = "CEN RES 3.37 5.47 6.55 7.41";
    }
}
else if (substr($fam, 0, 4) == "TAAR")
{
    $cenres = "CEN RES 3.32 3.37 5.43 6.48 7.43";
}
else die("Unsupported receptor family.");

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
        if ($ae < 0 && $ae < $ie) $data['Predicted'] = 'Agonist';
        else if ($a1 < 0 && $a1 < $i1) $data['Predicted'] = 'Agonist';
        else if ($a1 < 0 && $a1 > $i1) $data['Predicted'] = 'Inverse Agonist';
        else $data['Predicted'] = 'Non-agonist';

        echo "\nResult: " . print_r($data, true) . "\n";
    }

    return $data;
}

$outfname = "output/$fam/$protid/$protid.$ligname.inactive.dock";

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres
SIZE 7.0 7.0 7.0

EXCL 1 56		# Head, TMR1, and CYT1.

SEARCH TS
POSE 10
ELIM 500
ITERS 50
PROGRESS

FLEX 1
WET

OUT $outfname


heredoc;

chdir(__DIR__);
chdir("..");
if (!file_exists("output/$fam")) mkdir("output/$fam");
if (!file_exists("output/$fam/$protid")) mkdir("output/$fam/$protid");

process_dock("i");


chdir(__DIR__);
chdir("..");

$original_pdb = $pdbfname;
$pdbfname = str_replace(".upright.pdb", ".icactive.pdb", $original_pdb);
$result = [];
exec("php -f predict/icactive.php $protid", $result);

$config_params = [];
$reading_params = false;
$flex_constraints = [];
foreach ($result as $line)
{
    if ($reading_params)
    {
        $config_params[] = $line;
        if (substr($line, 0, 5) == "STCR " || substr($line, 0, 5) == "FLXR ") $flex_constraints[] = $line;
    }
    else if (trim($line) == "PRIMARYDOCK_CONFIG_PARAMS:") $reading_params = true;
}

$config_params = implode("\n", $config_params);
$flex_constraints = implode("\n", $flex_constraints);

if (!file_exists($pdbfname))
{
    $config_params = <<<heredoc

PROT $original_pdb
$config_params
OUTPDB $pdbfname

heredoc;

    $fp = fopen("tmp/icactive.config", "wb") or die("Failed to write to tmp/ folder.\n");
    fwrite($fp, $config_params);
    fclose($fp);

    passthru("bin/ic_active_pdb tmp/icactive.config");
}

$outfname = "output/$fam/$protid/$protid.$ligname.icactive.dock";

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres
SIZE 7.0 7.0 7.0

EXCL 1 56		# Head, TMR1, and CYT1.

SEARCH TS
POSE 10
ELIM 500
$flex_constraints
ITERS 50
PROGRESS

FLEX 1
WET

OUT $outfname
OUTPDB 1 output/$fam/$protid/%p.%l.model%o.pdb


heredoc;

process_dock("a");
