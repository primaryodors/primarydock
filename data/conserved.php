<?php

chdir(__DIR__);
chdir("..");

require("predict/statistics.php");
require("data/odorutils.php");

$colors =
[
    "PAILV" => "\033[38;5;247m",
    "G" => "\033[38;5;243m",
    "M" => "\033[38;5;185m",
    "C" => "\033[38;5;220m",
    "STNQ" => "\033[38;5;49m",
    "DE" => "\033[38;5;196m",
    "KR" => "\033[38;5;27m",
    "H" => "\033[38;5;99m",
    "FWY" => "\033[38;5;171m"
];

$orphan = "\033[2m\033[3m";
$deorphan = ""; // "\033[1m";
$reset = "\033[22m\033[23m\033[24m";

echo "Legend: {$deorphan}Deorphaned receptor$reset {$orphan}Orphan receptor$reset\n\n";

foreach (@$argv as $a)
{
	$a = explode('=',$a,2);
	$_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}

$bw = @$_REQUEST["bw"];

if (!$bw) die("Error no BW number.\n");

$list = @$_REQUEST["list"];
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
        if (!$j)
        {
            $lookfor = "HXR" . substr($rgno, 0, 1);
            $j = strpos($ln, $lookfor);
        }
        if ($j < 15) continue; // die("Unknown region $rgno.\n");

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
        for ($i=0; ; $i+=sgn($rel))
        {
            $i1 = $col50+$i;
            if ($i1<15 || $i1>=strlen($ln)) continue 2;
            $c = substr($ln, $i1, 1);
            if ($c != " ")
            {
                if ($j == $rel)
                {
                    $col = $col50 + $i;
                    break;
                }
                $j += sgn($rel);
            }
        }
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
        $ligands = count(all_empirical_pairs_for_receptor($orid, true, true));
        if ($ligands) $$varhas[$c][] = "$deorphan$orid$reset";
        else $$varhas[$c][] = "$orphan$orid$reset";
    }
}


foreach (["ttpd", "fish", "taar", "vn1r"] as $k => $var)
{
    if (!count($$var)) continue;
    arsort($$var);
    $varhas = $var."_has";

    echo ["Class II", "Class I", "Trace amine", "Vomeronasal"][$k] . ":\n";

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
        if (!$c) continue;

        foreach ($colors as $kc => $esc) if (false !== strpos($kc, $c)) echo $esc;

        echo "$c ";
        $pcnt = round(floatval(array_values($$var)[$i]) / $ttl * 100, 3);
        echo "$pcnt%";

        if ($i > $fnd80 || $c == $list) echo " " . implode(" ", $$varhas[$c]);

        echo "\033[0m\n";
    }

    echo "\n";
}
