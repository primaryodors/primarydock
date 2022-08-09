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
chdir("..");
$prots = json_decode(file_get_contents("data/receptor.json"), true);

if (!file_exists("tmp")) mkdir("tmp");

foreach (@$argv as $a)
{
	$a = explode('=',$a,2);
	$_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}

$protid = @$_REQUEST['prot'] ?: "OR1A1";

if (substr($protid, 0, 2) == "OR") $fam = "OR".intval(substr($protid, 2, 2));
else $fam = substr($protid, 0, 4);

$lockarom = resno_from_bw($protid, "6.40");
$lockhphl = resno_from_bw($protid, "3.39");
$pivot5 = resno_from_bw($protid, "5.47");
$pivot6 = $lockarom + 3;

$endof5 = intval($prots[$protid]["region"]["TMR5"]["end"]);
$start6 = intval($prots[$protid]["region"]["TMR6"]["start"]);

$midway = intval(($endof5+$start6)/2) + 2;
$midwa1 = $midway+1;
$startconn = $midway-5;

$pdisdat = <<<heredoc

LOAD "pdbs/$fam/$protid.rotated.pdb"

LET &bendamt = 0.5
LET &mbendamt = 0.0 - &bendamt
LET &curr_angle = 0
LET &best_angle = 0
LET &best_energy = 0		# The final best energy should be a best binding and not a least clash.

# Begin main loop.
_loop:

BRIDGE $lockarom $lockhphl 100
BENERG $lockhphl $lockarom &e
ECHO &curr_angle " degrees, " &e " energy...     " ~

LET \$file = "tmp/" + \$PROTEIN + ".acv" + &curr_angle + ".pdb"
SAVE \$file

IF &curr_angle < 10 GOTO _Ive_had_better
IF &e > &best_energy GOTO _Ive_had_better
LET &best_energy = &e
LET &best_angle = &curr_angle
_Ive_had_better:

BEND $pivot5 $midway "N-CA" &mbendamt
BEND $pivot6 $midwa1 "C-CA" &mbendamt

LET &curr_angle += &bendamt

IF &curr_angle < 35 GOTO _loop

ECHO "Best binding energy " &best_energy " found at " &best_angle " degrees."

LET \$file = "tmp/" + \$PROTEIN + ".acv" + &best_angle + ".pdb"
LOAD \$file
# CONNECT $startconn $midwa1 20
SAVE "pdbs/$fam/$protid.active.pdb" QUIT

heredoc;

$f = fopen("tmp/activate.pdis", "wb");
if (!$f) die("File write error. Check tmp folder is write enabled.\n");

fwrite($f, $pdisdat);
fclose($f);

$outlines = [];
exec("bin/interpreter tmp/activate.pdis", $outlines);

foreach ($outlines as $ln)
{
	if (false!==strpos($ln, "Best binding energy")) echo "$ln\n";
}














