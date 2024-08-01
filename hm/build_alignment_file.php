<?php

chdir(__DIR__);
require_once("../data/protutils.php");
chdir(__DIR__);

$exp = file_get_contents("experimental.ali");

foreach (explode("\n", $exp) as $ln)
{
    if (substr($ln, 0, 4) == ">P1;")
    {
        $rcsbid = substr($ln, 4);
        $fn = "$rcsbid.pdb";
        if (file_exists($fn)) continue;
        $rcsbidu = strtoupper($rcsbid);
        $url = "https://files.rcsb.org/download/$rcsbidu.pdb";
        file_put_contents($fn, file_get_contents($url));
        if (file_exists($fn)) echo "Downloaded $rcsbid.\n";
    }
}

$fp = fopen("allgpcr.ali", "w");
if (!$fp) die("FAIL; check folder permissions.\n");
fwrite($fp, $exp);
fwrite($fp, "\n\n");

foreach ($prots as $rcpid => $p)
{
    if (!isset($p['aligned'])) continue;

    $p1row = ">P1;$rcpid";
    fwrite($fp, "$p1row\n");

    $deets = "sequence:$rcpid:1     :A:".strlen($p['sequence'])."  :A:";
    $fam = family_from_protid($rcpid);
    $mem = member_from_protid($rcpid);
    switch ($fam)
    {
        case "TAAR":
        $deets .= "Trace amine-associated receptor $mem";
        break;

        case "VN1R":
        $deets .= "Vomeronasal type 1 receptor number $mem";
        break;

        case "MS4A":
        $deets .= "Membrane-spanning 4A receptor $mem";
        break;

        default:
        $famn = substr($fam, 2);
        $sub = subfamily_from_protid($rcpid);
        $deets .= "Olfactory receptor family $famn subfamily $sub number $mem";
    }
    $deets .= ":Homo sapiens: 1.90: 0.19";
    fwrite($fp, "$deets\n");

    fwrite($fp, "{$p['aligned']}*\n\n");
}

fclose($fp);
echo "Wrote alignments file.\n";

