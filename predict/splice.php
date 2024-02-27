<?php

// Necessary to add the AlphaFold version of the EXR2 loop to pdbs/OR5K1.1015.pdb
// since the latter doesn't have the right shape in this part of the protein.

require_once("data/protutils.php");

$lines1 = explode("\n", file_get_contents("pdbs/OR5K1.1015.pdb"));
$lines2 = explode("\n", file_get_contents("pdbs/OR5/OR5K1.upright.pdb"));

$start = resno_from_bw("OR5K1", "4.67");
$end   = resno_from_bw("OR5K1", "5.30");

$remarks = [];
$before_atoms = [];
$during_atoms = [];
$after_atoms = [];

foreach ($lines1 as $ln)
{
    if (substr($ln, 0, 7) == "REMARK ") $remarks[] = $ln;
    else if (substr($ln, 0, 5) == "ATOM ")
    {
        $resno = intval(substr($ln, 22, 4));
        if ($resno < $start) $before_atoms[] = $ln;
        else if ($resno > $end) $after_atoms[] = $ln;
    }
}

foreach ($lines2 as $ln)
{
    if (substr($ln, 0, 7) == "REMARK ") $remarks[] = $ln;
    else if (substr($ln, 0, 5) == "ATOM ")
    {
        $resno = intval(substr($ln, 22, 4));
        if ($resno >= $start && $resno <= $end) $during_atoms[] = $ln;
    }
}

$out = [];

foreach ($remarks as $ln) $out[] = $ln;

$atno = 1;
foreach ($before_atoms as $ln)
{
    $ln = substr($ln, 0, 7).str_pad($atno, 4, " ", STR_PAD_LEFT).substr($ln, 11);
    $out[] = $ln;
    $atno++;
}
foreach ($during_atoms as $ln)
{
    $ln = substr($ln, 0, 7).str_pad($atno, 4, " ", STR_PAD_LEFT).substr($ln, 11);
    $out[] = $ln;
    $atno++;
}
foreach ($after_atoms as $ln)
{
    $ln = substr($ln, 0, 7).str_pad($atno, 4, " ", STR_PAD_LEFT).substr($ln, 11);
    $out[] = $ln;
    $atno++;
}

$fp = fopen("pdbs/OR5K1.1015.spliced.pdb", "w");
if (!$fp) die("FAILED.\n");
fwrite($fp, implode("\n", $out));
fclose($fp);
