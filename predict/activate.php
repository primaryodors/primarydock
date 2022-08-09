<?php

// activate.php
//
// Activate a protein by bending its backbone at 6.40:N-CA until the optimum 3.39-6.40 bond occurrs, then bending 5.47:N-CA to match.
//
// Example call syntax:
// php -f predict/activate.php prot=OR1A1
//

function resno_from_bw($protid, $bw)
{
	global $prots;
	if (!isset($prots[$protid])) die("Protein not found: $protid.\n");
	
	$pettia = explode(".", $bw);
	$tmrno = intval($pettia[0]);
	$offset = intval($pettia[1]);
	
	$res50 = intval(@$prots[$protid]["bw"]["$tmrno.50"]) or die("Unknown Ballesteros-Weinstein number: $bw.\n");
	
	return $res50 + $offset - 50;
}

chdir(__DIR__);
$prots = json_decode(file_get_contents("../data/receptor.json"));

if (!file_exists("temp")) mkdir("temp");

foreach (@$argv as $a)
{
	$a = explode('=',$a,2);
	$_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}

$protid = @$_REQUEST['prot'] ?: "OR1A1";

$lockarom = resno_from_bw($protid, "6.40");
$lockhphl = resno_from_bw($protid, "3.39");
$pivot5 = resno_from_bw($protid, "5.47");
$pivot6 = $lockarom + 3;

$endof5 = intval($prots[$protid]["region"]["TMR5"]["end"]);
$start6 = intval($prots[$protid]["region"]["TMR6"]["start"]);

$midway = intval(($endof5+$start6)/2) + 2;
$midwa1 = $midway+1;

$pdisdat = <<<heredoc

LOAD "pdbs/OR1/$protid.rotated.pdb"

LET &bendamt = 0.5
LET &mbendamt = 0.0 - &bendamt
LET &curr_angle = 0
LET &best_angle = 0
LET &best_energy = 0		# The final best energy should be a best binding and not a least clash.

# Begin main loop.
_loop:

BRIDGE $lockarom $lockhphl1 100
BENERG $lockhphl1 $lockarom &e
ECHO &curr_angle " degrees, " &e " energy...     " ~

LET \$file = "temp/" + \$PROTEIN + ".acv" + &curr_angle + ".pdb"
SAVE \$file

IF &curr_angle < 10 GOTO _Ive_had_better
IF &e > &best_energy GOTO _Ive_had_better
LET &best_energy = &e
LET &best_angle = &curr_angle
_Ive_had_better:

BEND $pivot5 $midway "N-CA" &mbendamt
BEND $pivot6 $midwa1 "C-CA" &mbendamt
# CONNECT $endof5 $midwa1 50

LET &curr_angle += &bendamt

IF &curr_angle < 35 GOTO _loop

ECHO "Best binding energy " &best_energy " found at " &best_angle " degrees."

heredoc;

$f = fopen("temp/activate.pdis", "wb");
if (!$f) die("File write error. Check temp folder is write enabled.\n");

fwrite($f, $pdisdat);
fclose($f);

$outlines = [];
exec("bin/interpreter temp/activate.pdis", $outlines);
















