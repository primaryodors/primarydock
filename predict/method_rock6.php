<?php

// method_rock6.php
//
// Performs a dock of an odorant in a receptor, using the TM6 rock method
// based on the direct observation by Billesbølle et al (2022).
// https://doi.org/10.1101/2022.12.20.520951
//
// Example call syntax:
// php -f predict/method_rock6.php prot=OR51E2 lig=propionic_acid
//

// Configurable variables
$dock_metals = false;
$bias_by_energy = true;

chdir(__DIR__);
require("methods_common.php");
chdir(__DIR__);

$cenres = "CEN RES 3.33 4.57 4.60 45.52 45.53 5.39 5.43 6.55 6.59";

prepare_outputs();

$outfname = str_replace(".dock", "_rock6a.dock", $outfname);

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres
SIZE 5.0 4.0 5.0

ACVNODE 0
ACVHXR TMR6 6.26 6.60 6.48 0 0 0 -0.5448 0 0.8385 13

EXCL 1 56		# Head, TMR1, and CYT1.

POSE 10
ELIM 99

FLEX 1
H2O 15

ITERS 50

OUT $outfname
ECHO



heredoc;

process_dock("active_");

$outfnamei = str_replace("_rock6a.dock", "_rock6i.dock", $outfname);
$configf = str_replace("ACVNODE 0", "ACVNODE 9999", $configf);
$configf = str_replace("OUT $outfname", "OUT $outfnamei", $configf);
$outfname = $outfnamei;
process_dock("inactive_", true);
