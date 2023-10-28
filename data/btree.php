<?php

chdir(__DIR__);
include_once("protutils.php");

$tree = file_get_contents("btree");
$lines = explode("\n", $tree);

foreach ($prots as $rcpid => $pdata) unset($prots[$rcpid]['btree']);

foreach ($lines as $lno => $ln)
{
    if (substr($ln, 0, 1) == '#') continue;
    $pieces = explode(" ", $ln);
    if (count($pieces) != 2) continue;

    $binary = $pieces[0];
    $receptor = $pieces[1];

    if (!preg_match("/^[01]+$/", $binary)) continue;

    foreach (explode("/", $receptor) as $rcpid)
    {
        if (isset($prots[$rcpid]))
        {
            $prots[$rcpid]['btree'] = $binary;
            echo "$rcpid\n";
        }
    }
}

foreach ($prots as $rcpid => $pdata)
{
    if (substr($rcpid, 0, 2) != "OR") continue;
    if (!isset($pdata['btree']))
    {
        $seq1 = $pdata['sequence'];
        $best_match = false;
        $best_btree = false;
        $best_distance = 10000;
        foreach ($prots as $rcp2 => $p2)
        {
            if (!isset($p2['btree'])) continue;

            $seq2 = $p2['sequence'];
            $lev = levenshtein(substr($seq1, 0, 200), substr($seq2, 0, 200)) + levenshtein(substr($seq1, 200), substr($seq2, 200));
            $lev += levenshtein($rcpid, $rcp2);
            if ($lev < $best_distance)
            {
                $best_match = $rcp2;
                $best_btree = $p2['btree'];
                $best_distance = $lev;
            }
        }

        if ($best_match)
        {
            $prots[$best_match]['btree'] .= '0';
            $prots[$rcpid]['btree'] = $best_btree.'1';
            echo "$rcpid -- $best_match\n";
        }
        else die("Nobody likes $rcpid.\n");
    }
}

$fp = fopen("../data/receptor.json", "w");
if ($fp)
{
	fwrite($fp, json_encode_pretty($prots));
	fclose($fp);
}
