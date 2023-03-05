<?php

// method_dyn6.php
//
// Performs a path-based dock of an odorant in a receptor, using the TM6 rock method
// based on the direct observation by Billesbølle et al (2022).
// https://doi.org/10.1101/2022.12.20.520951
//
// Example call syntax:
// php -f predict/method_dyn6.php prot=OR51E2 lig=propionic_acid
//

// Configurable variables
$dock_metals = false;
$bias_by_energy = true;

chdir(__DIR__);
require("methods_common.php");
chdir(__DIR__);

$cenres = "CEN RES 6.59";
$path = [];
$path[] = "PATH ".(count($path)+1)." RES 5.40 5.43 5.44 5.47 6.51";
$path[] = "PATH ".(count($path)+1)." RES 3.33 45.52"; $acvnode = count($path);
$path[] = "PATH ".(count($path)+1)." RES 4.57 4.60 5.43";
$path = implode("\n", $path);

prepare_outputs();

$outfname = str_replace(".dock", "_dyn6.dock", $outfname);

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres
$path
SIZE 8.0 4.0 8.0

ACVNODE $acvnode
ACVHXR TMR6 6.26 6.60 6.48 0 0 0 -0.5448 0 0.8385 13

EXCL 1 56		# Head, TMR1, and CYT1.

POSE 10
ELIM 99

FLEX 1
# H2O 15

ITERS 50

OUT $outfname
ECHO



heredoc;

process_dock();

/*
process_dock("active_");

$outfnamei = str_replace("_rock6a.dock", "_rock6i.dock", $outfname);
$configf = str_replace("ACVNODE 0", "ACVNODE 9999", $configf);
$configf = str_replace("OUT $outfname", "OUT $outfnamei", $configf);
$outfname = $outfnamei;
process_dock("inactive_", true);*/
