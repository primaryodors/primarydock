<?php

// method_acv1st.php
//
// Performs a dock of an odorant in a receptor using first the active and then the inactive configurations.
//
// Example call syntax:
// php -f predict/method_acv1st.php prot=OR1A1 lig=geraniol
//


// Configurable variables
$dock_metals = false;
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
foreach ($rotations as $region => $values)
{
	$lregion = substr($region,0,4);
	$sr = $prots[$protid]['region'][$lregion]['start'];
	$er = $prots[$protid]['region'][$lregion]['end'];

	switch($lregion)
	{
		case 'TMR2': $mr = resno_from_bw($protid, "2.53"); break;
		case 'TMR3': $mr = resno_from_bw($protid, "3.36"); break;
		case 'TMR4': $mr = resno_from_bw($protid, "4.57"); break;
		case 'TMR5': $mr = resno_from_bw($protid, "5.43"); break;
		case 'TMR6': $mr = resno_from_bw($protid, "6.48"); break;
		case 'TMR7': $mr = resno_from_bw($protid, "7.46"); break;

		default:
		$mr = intval(($sr + $er) / 2);
	}

	if (substr($region, -1) == 'n') $er = $mr;
	if (substr($region, -1) == 'c') $sr = $mr;
	$acv_matrix .= "ACVHXR $region $sr $er $mr " . implode(" ", $values)."\n";
}


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

DEACVNODE 1

POSE 10
ITER 50
ELIM 50

OUT $outfname
ECHO


heredoc;


process_dock();

