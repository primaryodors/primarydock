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

$cenres = "CEN RES 4.57 4.60 5.39";

prepare_outputs();

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres
SIZE 8.0 5.0 8.0

ACVNODE 0
ACVHXR TMR6 6.26 6.60 6.48 0 0 0 0.5015 0 0.8651 10?

EXCL 1 56		# Head, TMR1, and CYT1.

POSE 10
ELIM 99

FLEX 1
# H2O 15

ITERS 50

OUT $outfname
ECHO



heredoc;

process_dock("");
