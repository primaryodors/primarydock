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
$bias_by_energy = true;
$elima =  5;
$elimi = 20;

chdir(__DIR__);
require("methods_common.php");
chdir(__DIR__);

$cenres = "CEN RES 2.53 3.29 3.32 3.33 3.36 3.37 3.40 3.41 4.53 4.57 4.60 45.49 45.52 5.39 5.43 5.46 5.47 6.48 6.51 7.38 7.39 7.42";

if (!isset($_REQUEST['force']))
{
    $results = [];
    exec("ps -ef | grep bin/couple | grep -v sh | grep -v grep", $results);
    if ($results) exit;
}

prepare_outputs();

$opdbname = $pdbfname;
$ooutname = $outfname;


function make_prediction($data)
{
    if (!isset($data['inactive'])) return $data;
    $i = floatval($data['inactive']);
    if (isset($data['hGNAL']) && floatval($data['hGNAL']) < $i) $data['Predicted'] = 'Agonist';
    else if (isset($data['hGNAS2']) && floatval($data['hGNAS2']) < $i) $data['Predicted'] = 'Agonist';
    else if (isset($data['hGNAL']) && floatval($data['hGNAL']) > $i) $data['Predicted'] = 'Non-agonist';
    else if (isset($data['hGNAS2']) && floatval($data['hGNAS2']) > $i) $data['Predicted'] = 'Non-agonist';
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
POSE 10
ELIM $elimi

FLEX 1
# H2O 15
WET

ITERS 50

OUT $outfname
ECHO



heredoc;

process_dock("");



// Golf a.k.a. hGNAL

$pdbfname = str_replace("pdbs/", "pdbs/coupled/", $opdbname);
$pdbfname = str_replace(".upright.", "_hGNAL.", $pdbfname);

$outfname = str_replace("$protid-", "{$protid}_hGNAL-", $ooutname);

if (!file_exists($pdbfname) || filesize($pdbfname) < 100000)
{
    $fp = fopen($outfname, "wb");
    fwrite($fp, "No protein.");
    fclose($fp);
}

$metrics_to_process =
[
  "BENERG" => "hGNAL",
  "BENERG.rgn" => "hGNAL.rgn",
];

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres
SIZE 7.0 7.0 7.0

EXCL 1 56		# Head, TMR1, and CYT1.

SEARCH TS
POSE 10
ELIM $elima

FLEX 1
# H2O 15
WET

ITERS 50

OUT $outfname
ECHO



heredoc;

process_dock("", true);



// Gαs a.k.a. hGNAS2

$pdbfname = str_replace("pdbs/", "pdbs/coupled/", $opdbname);
$pdbfname = str_replace(".upright.", "_hGNAS2.", $pdbfname);

$outfname = str_replace("$protid-", "{$protid}_hGNAS2-", $ooutname);

if (!file_exists($pdbfname) || filesize($pdbfname) < 100000)
{
    $fp = fopen($outfname, "wb");
    fwrite($fp, "No protein.");
    fclose($fp);
}

$metrics_to_process =
[
  "BENERG" => "hGNAS2",
  "BENERG.rgn" => "hGNAS2.rgn",
];

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres
SIZE 7.0 7.0 7.0

EXCL 1 56		# Head, TMR1, and CYT1.

SEARCH TS
POSE 10
ELIM $elima

FLEX 1
# H2O 15
WET

ITERS 50

OUT $outfname
ECHO



heredoc;

process_dock("", true);




