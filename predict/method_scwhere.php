<?php

// method_scwhere.php
//
// Perform a dock of an odorant in a receptor and record the relative positions of side chain atoms.
//
// Example call syntax:
// php -f predict/method_scwhere.php prot=OR1D2 lig=floralozone
//

// Configurable variables
$dock_metals = true;
$bias_by_energy = true;

require("methods_common.php");

$nodeno = 0;
$paths = [];	
$cenres = "CEN RES $bsr2a $bsr3a $bsr3b $bsr3c $bsr3d $bsr3e $bsr3f $bsr3g $bsr4a $bsr4b $bsr4c $bsr5a $bsr5b $bsr5c $bsr5d $bsr6a $bsr6b $bsr7a $bsr7b $bsr7c";

$paths = implode("\n", $paths);

$tmr2start = $prots[$protid]['region']['TMR2']['start'];
$cyt1end = $tmr2start - 1;
$tmr4end = $prots[$protid]['region']['TMR4']['end'];
$tmr5start = $prots[$protid]['region']['TMR5']['start'];
$exr2start = $tmr4end + 1;
$exr2end = $tmr5start - 1;

prepare_outputs();

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres
SIZE 7.0 7.5 7.0
EXCL $tmr4end $tmr5start
# H2O 1

POSE 10
ITER 50
ELIM -0.001

OUT $outfname
ECHO


heredoc;


process_dock();

