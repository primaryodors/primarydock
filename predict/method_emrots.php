<?php

// method_emrots.php
//
// Performs a path-based dock of an odorant in a receptor, using TMR rotations
// to approximate the PDB file from direct observation by Billesbølle et al (2022).
// https://doi.org/10.1101/2022.12.20.520951
//
// Example call syntax:
// php -f predict/method_emrots.php prot=OR51E2 lig=propionic_acid
//

// Configurable variables
$dock_metals = false;
$bias_by_energy = true;

chdir(__DIR__);
require("methods_common.php");
chdir(__DIR__);

prepare_outputs();

# Active rotations from A51E2.pepd:
$rotations =
[
	'TMR1' => [  0.046653, -0.683896,  0.988195, -0.403776, -0.078198, -0.911510,  5.233483 ],
	'TMR2' => [ -0.075481, -0.803532,  0.298527,  0.461675, -0.314792, -0.829314,  1.933021 ],
	'TMR3' => [  0.056110, -0.062618, -0.694922,  0.804089,  0.342890,  0.485663,  2.996447 ],
	'TMR4' => [ -0.018695, -0.978119,  0.699819,  0.903594, -0.170428,  0.393030,  1.552474 ],
	'TMR5' => [  0.037672,  0.319781, -2.055663,  0.973715,  0.223967,  0.041455,  7.376118 ],
	'TMR6' => [  0.405011,  1.857748, -0.752629,  0.951442, -0.048377,  0.304004, 14.500003 ],
	'TMR7' => [ -0.395103, -0.028852,  2.120022,  0.287893, -0.032532, -0.957110,  6.944020 ],
];


$acv_matrix = "";
foreach ($rotations as $region => $values)
{
	$lregion = substr($region,0,4);
	$sr = $prots[$protid]['region'][$lregion]['start'];
	$er = $prots[$protid]['region'][$lregion]['end'];

    $values[0] = $values[1] = $values[2] = 0;

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

	$acv_matrix .= "ACVHXR $region $sr $er $mr " . implode(" ", $values)."\n";
}


$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

CEN RES 4.57 4.60 5.39
PATH 1 RES 4.57 4.60 5.39 6.59

REQSR ALL 4.57 4.60 5.39

SIZE 7.5 6.0 7.5

# Correct for slight misalignment in AlphaFold model
HXR 1 TMR6 6.49 6.60 6.49 0 0 0 helical -20

# Activation Matrix.
$acv_matrix
ACVNODE 1

EXCL 1 56		# Head, TMR1, and CYT1.

POSE 10
ELIM 20

FLEX 1
FLXR 6.59 45.53
STCR 4.56 4.57 4.60
WET
# MOVIE

ITERS 50

OUT $outfname
ECHO



heredoc;

process_dock();


