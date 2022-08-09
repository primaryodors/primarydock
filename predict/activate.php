<?php

// activate.php
//
// Activate a protein by bending its backbone at 6.40:N-CA until the optimum 3.39-6.40 bond occurrs, then bending 5.47:N-CA to match.
//
// Example call syntax:
// php -f predict/activate.php prot=OR1A1
//

chdir(__DIR__);
$prots = json_decode(file_get_contents("../data/receptor.json"));

foreach (@$argv as $a)
{
	$a = explode('=',$a,2);
	$_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}

$protid = @$_REQUEST['prot'] ?: "OR1A1";

$pdisdat = <<<heredoc

LOAD "pdbs/OR1/$protid.rotated.pdb"

LET &bendamt = 0.5
LET &mbendamt = 0.0 - &bendamt
LET &curr_angle = 0
LET &best_angle = 0
LET &best_energy = 0		# The final best energy should be a best binding and not a least clash.

# Begin main loop.
_loop:

BRIDGE 251 111 100
BENERG 111 251 &e
ECHO &curr_angle " degrees, " &e " energy...     " ~

LET \$file = "output/" + \$PROTEIN + ".acv" + &curr_angle + ".pdb"
SAVE \$file

IF &curr_angle < 10 GOTO _Ive_had_better
IF &e > &best_energy GOTO _Ive_had_better
LET &best_energy = &e
LET &best_angle = &curr_angle
_Ive_had_better:

BEND 206 229 "N-CA" &mbendamt
BEND 254 230 "C-CA" &mbendamt
# CONNECT 218 230 50

LET &curr_angle += &bendamt

IF &curr_angle < 35 GOTO _loop

ECHO "Best binding energy " &best_energy " found at " &best_angle " degrees."



heredoc;



