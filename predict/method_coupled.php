<?php

// method_coupled.php
//
// Performs a dock of an odorant inside a receptor that has previously been coupled
// with a G protein using the generate_couple.php script.
//
// Example call syntax:
// php -f predict/method_coupled.php prot=OR1A1 lig=R-carvone
//

// Configurable variables
$dock_metals = false;

// Internal variables
$elima =  5;
$elimi = 25;

chdir(__DIR__);
require("methods_common.php");
chdir(__DIR__);


if (!isset($_REQUEST['force']))
{
    $results = [];
    exec("ps -ef | grep bin/couple | grep -v sh | grep -v grep", $results);
    if ($results) exit;
    $results = [];
    exec("ps -ef | grep bin/pepteditor | grep -v sh | grep -v grep", $results);
    if ($results) exit;
}

prepare_outputs();

$opdbname = $pdbfname;
$ooutname = $outfname;

// $cenres = "CEN RES 2.53 3.29 3.32 3.33 3.36 3.37 3.40 3.41 4.53 4.57 4.60 45.49 45.52 5.39 5.43 5.46 5.47 6.48 6.51 7.38 7.39 7.42";

switch ($fam)
{
    case "OR51":
    case "OR52":
    case "OR56":
    $cenres = "CEN RES 3.33 4.57 4.60 5.39 45.52\nREQSR 0 45.52\nSTCR 3.33 4.57 4.60 5.39";
    break;

    default:
    $cenres = "CEN RES 3.37 4.60 5.47 6.55";
}

function make_prediction($data)
{
    if (!isset($data['inactive']) && !isset($data['hGNAO1'])) return $data;
    if (!isset($data['hGNAL']) && !isset($data['hGNAS2']) && !isset($data['hGNAQ'])) return $data;
    
    $i = floatval(@$data['inactive'] ?: 0);
    $o = floatval(@$data['hGNAO1'] ?: 0);

    $l = floatval(@$data['hGNAL'] ?: 0);
    $s = floatval(@$data['hGNAS2'] ?: 0);
    $q = floatval(@$data['hGNAQ'] ?: 0);

    if ($i >= 13.5) $data['Predicted'] = 'Non-agonist';
    else
    {
        if ($i < 0 && $l >= 0 && $s >= 0 && $q >= 0) $data['Predicted'] = 'Inverse Agonist';
        if ($l < $i && $l < $o && $l < 0) $data['Predicted'] = 'Agonist';
        else if ($s < $i && $s < $o && $s < 0) $data['Predicted'] = 'Agonist';
        else if ($q < $i && $q < $o && $q < 0) $data['Predicted'] = 'Agonist';
        else $data['Predicted'] = 'Non-agonist';
    }

    return $data;
}


// Inactive state

$metrics_to_process =
[
  "BENERG" => "inactive",
  "BENERG.rgn" => "inactive.rgn",
];

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres
SIZE 7.0 7.0 7.0

EXCL 1 56		# Head, TMR1, and CYT1.

SEARCH TS
POSE 4
ELIM $elimi

FLEX 1
# H2O 15
WET

ITERS 100

OUT $outfname
ECHO
OPEND



heredoc;

process_dock("");



foreach (["hGNAL", "hGNAS2", "hGNAQ", "hGNAO1"] as $gpid)
{
    $pdbfname = str_replace("pdbs/", "pdbs/coupled/", $opdbname);
    $pdbfname = str_replace(".upright.", "_$gpid.", $pdbfname);

    $outfname = str_replace("$protid-", "{$protid}_$gpid-", $ooutname);

    if (!file_exists($pdbfname) || filesize($pdbfname) < 100000)
    {
        $fp = fopen($outfname, "wb");
        fwrite($fp, "No protein.");
        fclose($fp);
    }

    $metrics_to_process =
    [
      "BENERG" => "$gpid",
      "BENERG.rgn" => "$gpid.rgn",
    ];

    $configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres
SIZE 7.0 7.0 7.0

EXCL 1 56		# Head, TMR1, and CYT1.

SEARCH TS
POSE 4
ELIM $elima

FLEX 1
# H2O 15
WET

ITERS 40

OUT $outfname
ECHO
OPEND



heredoc;

    process_dock("", true);
}

