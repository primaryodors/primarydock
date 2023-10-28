<?php

chdir(__DIR__);
include_once("protutils.php");

function long_levenshtein($str1, $str2)
{
    return levenshtein(substr($str1, 0, 200), substr($str2, 0, 200)) + levenshtein(substr($str1, 200), substr($str2, 200));
}

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
            $lev = long_levenshtein($seq1, $seq2);
            // $lev += levenshtein($rcpid, $rcp2);
            if ($lev < $best_distance)
            {
                $best_match = $rcp2;
                $best_btree = $p2['btree'];
                $best_distance = $lev;
            }
        }

        if ($best_match)
        {
            /*
            $prots[$best_match]['btree'] .= '0';
            $prots[$rcpid]['btree'] = $best_btree.'1';
            */
            // echo "$rcpid -- $best_match\n";

            $seq2 = $prots[$best_match]['sequence'];
            $mrca = $best_btree;
            $n = strlen($mrca);
            foreach ($prots as $p3)
            {
                if (!isset($p3['btree'])) continue;
                $seq3 = $p3['sequence'];
                $lev = long_levenshtein($seq2, $seq3);
                if ($lev < $best_distance)
                {
                    foreach (str_split($p3['btree']) as $i => $c)
                    {
                        if ($i >= $n) break;
                        if ($c != substr($mrca, $i, 1))
                        {
                            $mrca = substr($mrca, 0, $i);
                            $n = strlen($mrca);
                            break;
                        }
                    }
                }
            }

            // echo "MRCA is $mrca\n";
            $n = strlen($mrca);
            foreach ($prots as $l => $p3)
            {
                if (!isset($p3['btree'])) continue;
                if (substr($p3['btree'], 0, $n) === $mrca)
                {
                    // echo "\t$l btree was {$p3['btree']}";
                    $prots[$l]['btree'] = $mrca.'0'.substr($p3['btree'], $n);
                    // echo " now {$prots[$l]['btree']}\n";
                }
            }

            $prots[$rcpid]['btree'] = $mrca.'1';
            echo "$rcpid btree is now {$prots[$rcpid]['btree']}\n";
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
