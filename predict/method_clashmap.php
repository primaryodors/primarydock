<?php

// method_clashmap.php
//
// Performs a strongly clashing best-binding dock of an odorant in a receptor to search
// for regions that could move to accommodate the ligand in the protein's active state.
//
// Example call syntax:
// php -f predict/method_clashmap.php prot=OR2AT4 lig=sandalore
//


// Configurable variables
$dock_metals = false;
$bias_by_energy = true;

require("methods_common.php");

$nodeno = 0;
$paths = [];	
$cenres = "CEN RES $bsr2a $bsr3a $bsr3b $bsr3c $bsr3d $bsr3e $bsr3f $bsr3g $bsr4a $bsr4b $bsr4c $bsr5a $bsr5b $bsr5c $bsr5d $bsr6a $bsr6b $bsr7a $bsr7b $bsr7c";

$paths = implode("\n", $paths);

prepare_outputs();

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres
$paths
SIZE 7.0 8.5 7.0

POSE 5
ITER 8
ELIM 99999
SEARCH BB

OUT $outfname
ECHO


heredoc;


process_dock();

