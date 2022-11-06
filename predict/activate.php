<?php

// activate.php
//
// Activates a GPCR by bending its backbone at a specific TMR6 residue until the optimum TMR6-TMR3 bridge occurrs, then bending TMR5 to match.
//
// Example call syntax:
// php -f predict/activate.php prot=OR1A1
//
// May be run as a cron to generate active PDBs for all olfactory receptors, e.g.
// php -f predict/activate.php next
//
// May be run on the entire set of receptors in one call, e.g.
// php -f predict/activate.php all
//

require("protutils.php");

if (!file_exists("tmp")) mkdir("tmp");

foreach (@$argv as $a)
{
	$a = explode('=',$a,2);
	$_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}

_loop4all:
$protid = @$_REQUEST['prot'] ?: "OR1A1";

if ($protid == "next" || $protid == "all")
{
	foreach ($prots as $k => $v)
	{
		$protid = $k;
		$fam = family_from_protid($protid);
		if (!file_exists("pdbs/$fam/$protid.active.pdb")) goto _found_next;
	}
	die("All finished!\n");
	
	_found_next:
	echo "Processing $protid... ";
}

if (substr($protid, 0, 2) != "OR") die("Only ORs are currently supported.\n");
$fam = family_from_protid($protid);

$lockarom = resno_from_bw($protid, "6.40");
$lockhphl = resno_from_bw($protid, "3.39");
$pivot5 = resno_from_bw($protid, "5.47");
$pivot6 = $lockarom + 3;

$endof5 = intval($prots[$protid]["region"]["TMR5"]["end"]);
$start6 = intval($prots[$protid]["region"]["TMR6"]["start"]);

$midway = intval(($endof5+$start6)/2) + 2;
$midwa1 = $midway+1;
$startconn = $midway-5;

$pepddat = <<<heredoc

LOAD "pdbs/$fam/$protid.upright.pdb"

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

IF &curr_angle < 20 GOTO _loop

# IF &best_energy < 10 ECHO "Failed to find best binding energy. " &best_energy
# IF &best_energy < 10 EXIT

ECHO "$protid" " best binding energy " &best_energy " found at " &best_angle " degrees."

LET \$file = "tmp/" + \$PROTEIN + ".acv" + &best_angle + ".pdb"
LOAD \$file
# CONNECT $startconn $midwa1 20
SAVE "pdbs/$fam/$protid.active.pdb" QUIT

heredoc;

$f = fopen("tmp/activate.pepd", "wb");
if (!$f) die("File write error. Check tmp folder is write enabled.\n");

fwrite($f, $pepddat);
fclose($f);

$outlines = [];
set_time_limit(300);
exec("bin/peptiditor tmp/activate.pepd", $outlines);

foreach ($outlines as $ln)
{
	if (false!==strpos($ln, "est binding energy")) echo "$ln\n";
}

if (@$_REQUEST['prot'] == "all") goto _loop4all;












