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

// $outfname = str_replace(".dock", "_emrot.dock", $outfname);

# Active rotations from A51E2.pepd:
$rotations =
[
	'TMR1.n' => [ -1.433261, -0.495408, -0.108732,  0.127823, -0.141131, -0.981704,  9.203959 ],
	'TMR1.c' => [ -1.433261, -0.495408, -0.108732, -0.888311,  0.399802, -0.225970,  5.120095 ],

	'TMR2.n' => [ -0.369847, -0.596285,  0.435198,  0.881323,  0.026436, -0.471774,  1.794605 ],
	'TMR2.c' => [ -0.369847, -0.596285,  0.435198,  0.062762, -0.431871, -0.899749,  2.726275 ],

	'TMR3.n' => [ -0.460842, -0.310634,  0.145399,  0.931975,  0.350756,  0.091609,  0.550008 ],
	'TMR3.c' => [ -0.460842, -0.310634,  0.145399,  0.790443,  0.325104,  0.519141,  5.437398 ],

	'TMR4.n' => [  0.155436, -0.635536,  0.519279,  0.950716, -0.149542,  0.271617,  0.755370 ],
	'TMR4.c' => [  0.155436, -0.635536,  0.519279,  0.882852, -0.158706,  0.442023,  1.967520 ],

	'TMR5.n' => [  0.025806, -0.135105, -1.168056,  0.996638,  0.056400, -0.059425,  5.066973 ],
	'TMR5.c' => [  0.025806, -0.135105, -1.168056,  0.917328,  0.376384,  0.129786,  9.479465 ],

	'TMR6.n' => [  2.101592,  1.852155, -1.640432,  0.930834,  0.343676,  0.124235, 11.108210 ],
	'TMR6.c' => [  2.101592,  1.852155, -1.640432,  0.831236, -0.417401,  0.367183, 17.869820 ],

	'TMR7.n' => [ -0.131216,  1.079784,  0.610644,  0.459033, -0.414481, -0.785808, 11.035672 ],
	'TMR7.c' => [ -0.131216,  1.079784,  0.610644, -0.436363,  0.148272, -0.887470,  6.800398 ],

];

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


$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

CEN RES 4.57 4.60 5.39
PATH 1 RES 4.57 4.60 5.39 6.59

REQSR ALL 4.57 4.60 5.39

SIZE 8.0 6.0 8.0

# Correct for slight misalignment in AlphaFold model
HXR 1 TMR6 6.49 6.60 6.49 0 0 0 helical -20

# Activation Matrix.
$acv_matrix
ACVNODE 1

EXCL 1 56		# Head, TMR1, and CYT1.

POSE 1
ELIM 199

FLEX 1
WET
# MOVIE

ITERS 50

OUT $outfname
ECHO



heredoc;

process_dock();


