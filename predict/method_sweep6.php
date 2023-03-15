<?php

// method_sweep6.php
//
// Performs a path-based dock of an odorant in a receptor, using the TM6 rock method
// based on the direct observation by Billesbølle et al (2022).
// https://doi.org/10.1101/2022.12.20.520951
//
// Example call syntax:
// php -f predict/method_sweep6.php prot=OR51E2 lig=propionic_acid
//

// Configurable variables
$dock_metals = false;
$bias_by_energy = true;

chdir(__DIR__);
require("methods_common.php");
chdir(__DIR__);

prepare_outputs();

$outfname = str_replace(".dock", "_sweep6.dock", $outfname);

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

CEN RES 4.57 4.60 5.39 6.59
PATH 1 REL 0 0 0
PATH 2 REL 0 0 0
PATH 3 REL 0 0 0
PATH 4 REL 0 0 0
PATH 5 REL 0 0 0

SIZE 8.0 4.0 8.0

HXR 0 TMR6 6.26 6.60 6.48 0 2 0 0.5015 0 0.8651 0
HXR 1 TMR6 6.26 6.60 6.48 0 0 0 0.5015 0 0.8651 3
HXR 2 TMR6 6.26 6.60 6.48 0 0 0 0.5015 0 0.8651 3
HXR 3 TMR6 6.26 6.60 6.48 0 0 0 0.5015 0 0.8651 3
HXR 4 TMR6 6.26 6.60 6.48 0 0 0 0.5015 0 0.8651 3
HXR 5 TMR6 6.26 6.60 6.48 0 0 0 0.5015 0 0.8651 3

EXCL 1 56		# Head, TMR1, and CYT1.

POSE 5
ELIM 99

FLEX 1
WET

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
