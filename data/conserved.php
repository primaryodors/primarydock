<?php

chdir(__DIR__);
chdir("..");

require("predict/statistics.php");

foreach (@$argv as $a)
{
	$a = explode('=',$a,2);
	$_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}

$bw = @$_REQUEST["bw"] ?: "6.48";
$pieces = explode(".", $bw);
$rgno = intval($pieces[0]);
$ofst = intval(@$pieces[1]) ?: 50;

$c = file_get_contents("data/sequences_aligned.txt");
$lines = explode("\n", $c);

$ttpd = [];
$fish = [];
$taar = [];
$vn1r = [];

$ttpd_has = [];
$fish_has = [];
$taar_has = [];
$vn1r_has = [];

foreach ($lines as $ln)
{
    if (preg_match("/TMR[0-9]-*/", $ln))
    {
        $lookfor = "TMR" . substr($rgno, 0, 1);
        $j = strpos($ln, $lookfor);
        if ($j < 15) die("Unknown region $rgno.\n");

        if ($rgno < 10)
        {
            $col50 = strpos($ln, "|", $j);
        }
        else
        {
            $col50 = strpos($ln, " |", $j) + 1;
        }

        $col = false;
    }
    else if (substr($ln, 0, 2) == "% ")
    {
        $rel = $ofst - 50;
        $j = 0;
        for ($i=0; $j!=$rel; $i+=sgn($rel))
        {
            $i1 = $col50+$i;
            if ($i1<15 || $i1>=strlen($ln)) continue 2;
            $c = substr($ln, $i1, 1);
            if ($c != " ") $j += sgn($rel);
        }

        $col = $col50 + $j;
    }
    else
    {
        if (@$col < 15) continue;
        $orid = trim(substr($ln, 0, 15));
        if (!$orid) continue;
        $orid = explode(" ", $orid)[0];

        $c = substr($ln, $col, 1);
        if ($c == ' ') $c = '-';

        $var = false;
        if (substr($orid, 0, 2) == "OR")
        {
            $fam = intval(substr($orid, 2, 2));
            $var = ($fam < 50) ? "ttpd" : "fish";
        }
        else if (substr($orid, 0, 4) == "TAAR") $var = "taar";
        else if (substr($orid, 0, 4) == "VN1R") $var = "vn1r";

        if (!$var) continue;
        
        if (!isset($$var[$c])) $$var[$c] = 1;
        else $$var[$c]++;

        $varhas = $var."_has";
        $$varhas[$c][] = $orid;
    }
}


foreach (["ttpd", "fish", "taar", "vn1r"] as $k => $var)
{
    if (!count($$var)) continue;
    arsort($$var);
    $varhas = $var."_has";

    echo ["Tetrapod-like", "Fish-like", "Trace amine", "Vomeronasal"][$k] . ":\n";

    $n = count($$var);
    $ttl = array_sum($$var);
    $thr = 0.8 * $ttl;
    $fnd80 = $sum = 0;
    for ($i=0; $i<$n; $i++)
    {
        $sum += array_values($$var)[$i];
        if ($sum >= $thr)
        {
            $fnd80 = $i;
            break;
        }
    }

    for ($i=0; $i<$n; $i++)
    {
        $c = array_keys($$var)[$i];
        echo "$c ";
        $pcnt = round(floatval(array_values($$var)[$i]) / $ttl * 100, 3);
        echo "$pcnt%";

        if ($i > $fnd80) echo " " . implode(" ", $$varhas[$c]);

        echo "\n";
    }

    echo "\n";
}
