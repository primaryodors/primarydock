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

chdir(__DIR__);
require("methods_common.php");
chdir(__DIR__);

$cenres = "CEN RES 2.53 3.29 3.32 3.33 3.36 3.37 3.40 3.41 4.53 4.57 4.60 45.49 45.52 5.39 5.43 5.46 5.47 6.48 6.51 7.38 7.39 7.42";

prepare_outputs();

$pdbfname = str_replace("pdbs/", "pdbs/coupled/", $pdbfname);
$pdbfname = str_replace(".upright.", "_hGNAL.", $pdbfname);

$outfname = str_replace("$protid-", "{$protid}_hGNAL-", $outfname);

if (!file_exists($pdbfname) || filesize($pdbfname) < 100000)
{
    $fp = fopen($outfname, "wb");
    fwrite($fp, "No protein.");
    fclose($fp);
}

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres
SIZE 7.0 7.0 7.0

EXCL 1 56		# Head, TMR1, and CYT1.

POSE 10
ELIM 99

FLEX 1
# H2O 15
WET

ITERS 50

OUT $outfname
ECHO



heredoc;

process_dock("");



