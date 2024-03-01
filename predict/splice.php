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

$n3x25 = resno_from_bw("OR5K1", "3.25");
$n45x50 = resno_from_bw("OR5K1", "45.50");

foreach ($lines1 as $ln)
{
    if (substr($ln, 0, 7) == "REMARK ") $remarks[] = $ln;
    else if (substr($ln, 0, 5) == "ATOM ")
    {
        $resno = intval(substr($ln, 22, 4));
        if ($resno < $start) $before_atoms[] = $ln;
        else if ($resno > $end) $after_atoms[] = $ln;

        if (trim(substr($ln, 12, 4)) == "SG" && $resno == $n3x25)
        {
            $x325 = floatval(substr($ln, 30, 8));
            $y325 = floatval(substr($ln, 38, 8));
            $z325 = floatval(substr($ln, 46, 8));
            echo substr($ln,30,24)." $x325 $y325 $z325\n";
        }
    }
}

foreach ($lines2 as $ln)
{
    if (substr($ln, 0, 10) == "REMARK   1") $remarks[] = $ln;
    else if (substr($ln, 0, 5) == "ATOM ")
    {
        $resno = intval(substr($ln, 22, 4));
        if ($resno >= $start && $resno <= $end) $during_atoms[] = $ln;

        if (trim(substr($ln, 12, 4)) == "SG" && $resno == $n45x50)
        {
            $x45x50 = floatval(substr($ln, 30, 8));
            $y45x50 = floatval(substr($ln, 38, 8));
            $z45x50 = floatval(substr($ln, 46, 8));
            echo substr($ln,30,24)." $x45x50 $y45x50 $z45x50\n";
        }
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

if (isset($x325) && isset($x45x50))
{
    $x = $x325 - $x45x50;
    $y = $y325 - $y45x50;
    $z = $z325 - $z45x50;

    $x += 1.414;
    $z -= 1.414;

    /*$r = sqrt($x*$x + $y*$y + $z*$z);
    $mult = ($r-2.07)/$r;
    echo "r = $r; mult = $mult\n";

    $x *= $mult;
    $y *= $mult;
    $z *= $mult;*/
}
else
{
    echo "Warning: no Cys3.25-Cys45.50 cross link.\n";
    $x = $z = 0;
    $y = 7;
}

$pepd = <<<heredoc
LOAD pdbs/OR5K1.1015.spliced.pdb
LET @xlat45 = [$x,$y,$z]
MOVEREL $start $end @xlat45
SAVE pdbs/OR5K1.1015.spliced.pdb
SAVE pdbs/OR5/OR5K1.active.pdb
heredoc;

echo "$pepd\n";

$fn = "tmp/exr2.pepd";
$fp = fopen($fn, "wb");
fwrite($fp, $pepd);
fclose($fp);
exec("bin/pepteditor $fn");
unlink($fn);
