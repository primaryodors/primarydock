<?php

// method_acv1st.php
//
// Performs a dock of an odorant in a receptor using first the active and then the inactive configurations.
//
// Example call syntax:
// php -f predict/method_acv1st.php prot=OR1A1 lig=geraniol
//


// Configurable variables
$dock_metals = true;
$bias_by_energy = true;

require("methods_common.php");
require("matrix9.php");

$nodeno = 0;
$paths = [];	
$cenres = "CEN RES $bsr2a $bsr3a $bsr3b $bsr3c $bsr3d $bsr3e $bsr3f $bsr3g $bsr4a $bsr4b $bsr4c $bsr5a $bsr5b $bsr5c $bsr5d $bsr6a $bsr6b $bsr7a $bsr7b $bsr7c";

$paths[] = "PATH 1 REL 0 0 0";

$paths = implode("\n", $paths);

$tmr2start = $prots[$protid]['region']['TMR2']['start'];
$cyt1end = $tmr2start - 1;
$tmr4end = $prots[$protid]['region']['TMR4']['end'];
$tmr5start = $prots[$protid]['region']['TMR5']['start'];
$exr2start = $tmr4end + 1;
$exr2end = $tmr5start - 1;

$acv_matrix = "";
foreach ($rotations as $region => $values) $acv_matrix .= "ACVMR $region " . implode(" ", $values)."\n";

prepare_outputs();

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres
$paths
SIZE 7.0 7.5 7.0
EXCL $tmr4end $tmr5start
# H2O 1

# Activation Matrix.
$acv_matrix
ACVNODE 0

NODEPDB 1 pdbs/$fam/$protid.active.pdb

POSE 10
ITER 50
ELIM -0.001

OUT $outfname
ECHO


heredoc;


process_dock();

