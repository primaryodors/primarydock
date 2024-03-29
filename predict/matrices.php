<?php

global $matrixbtRho, $matrixADRB2;

// Immutable constant matrix values from btRho.pepd:
$matrixbtRho =
[
	"TMR1" => [ -1.432411, -2.746519,  2.509232,  2.363341, -3.109768, -3.402768 ],
	"TMR2" => [  4.639589,  0.024734, -2.108768, -1.090161, -2.695269,  1.679482 ],
	"TMR3" => [  0.167589, -2.024517,  3.395732,  4.231589,  1.114231, -4.413017 ],
	"TMR4" => [  1.722089,  0.166982, -1.611769, -2.093911,  1.330481,  2.336232 ],
	"TMR5" => [ -2.133661,  3.529982,  1.032483, -3.523911,  3.119732,  2.745732 ],
	"TMR6" => [ -2.557410,  1.492481, -5.381518, -2.909160,  3.353733,  2.692482 ],
	"TMR7" => [ -0.405411,  0.009982,  2.163982,  3.021839, -3.566268, -1.637518 ],
];

// btRho matrix adjustments
// Adjust TMR6 to not clash with TMR5.
$matrixbtRho["TMR6"][2] += 3;
$matrixbtRho["TMR6"][5] += 1;

// Equalize Y displacements.
foreach ($matrixbtRho as $region => $values)
{
	$y = ($values[1] + $values[4])/2;
	$matrixbtRho[$region][1] = $matrixbtRho[$region][4] = $y;
}

// No increase in distance between TMR5 and TMR3 at the extracellular end.
$matrixbtRho["TMR3"][0] -= 3;
$matrixbtRho["TMR3"][2] -= 1;
$matrixbtRho["TMR5"][0] += 3;
$matrixbtRho["TMR5"][2] -= 1;

// Prevent clash of TMR7 with TMR2 and TMR1. Bring TMR7 closer to TMR3.
$matrixbtRho["TMR7"][3] -= 3;
$matrixbtRho["TMR7"][5] -= 2;


# Immutable constant matrix values from ADRB2.pepd:
$matrixADRB2 =
[
	'TMR1' => [ -2.085748, 0.999161, 0.541523, 1.020319, 1.053395, -0.199438 ],
	'TMR2' => [ 1.129009, 0.317495, -0.306989, 1.990763, 0.589527, 0.498976 ],
	'TMR3' => [ -0.022306, 2.830248, -0.478422, -1.630260, -2.362781, -2.084338 ],
	'TMR4' => [ -0.949617, 2.118884, -0.106401, 0.695187, 1.876168, -0.352305 ],
	'TMR5' => [ -0.554504, 0.253978, -0.713764, -3.321619, -3.699796, -1.339734 ],
	'TMR6' => [ -5.414645, -3.399852, -2.154630, -0.162964, -4.181125, -0.244171 ],
	'TMR7' => [ 0.139565, -2.368010, 0.748168, 0.430551, -1.126985, 2.697166 ],
];
/*[
	'TMR1' => [  1.928375, -0.231194,  1.514458, -1.745764,  3.297985, -0.673183 ],
	'TMR2' => [ -1.637074,  2.562085, -0.780735,  6.004886, -0.640827,  1.471912 ],
	'TMR3' => [  3.991817,  1.599894,  0.494513, -4.396341, -0.118191, -2.558083 ],
	'TMR4' => [ -3.715700,  4.363474, -0.580147,  4.709310,  0.645814,  0.620629 ],
	'TMR5' => [  3.459619, -0.976376,  0.259171, -6.087701, -1.455206, -1.813480 ],
	'TMR6' => [ -8.180727, -1.155262, -2.628376,  3.851160, -5.411479,  0.728764 ],
	'TMR7' => [  4.153688, -3.598364,  1.721104, -2.335531,  1.117605,  2.223419 ],
];*/

// Equalize Y displacements.
foreach ($matrixADRB2 as $region => $values)
{
	$y = ($values[1] + $values[4])/2;
	$matrixADRB2[$region][1] = $matrixADRB2[$region][4] = $y;
}

// No increase in distance between TMR5 and TMR3 at the extracellular end
// or TMR7 and TMR3 overall.
// $matrixADRB2["TMR3"][0] -= 2;
$matrixADRB2["TMR3"][2] -= 1;
$matrixADRB2["TMR5"][0] += 1;
$matrixADRB2["TMR5"][2] -= 1;
$matrixADRB2["TMR7"][0] -= 2;
$matrixADRB2["TMR7"][2] -= 1;
$matrixADRB2["TMR7"][3] -= 2;
$matrixADRB2["TMR7"][5] -= 1;





