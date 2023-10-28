<?php

chdir(__DIR__);
include_once("protutils.php");

define("inhouse", false);

function long_levenshtein($str1, $str2)
{
    return levenshtein(substr($str1, 0, 200), substr($str2, 0, 200)) + levenshtein(substr($str1, 200), substr($str2, 200));
}

$tree = file_get_contents("btree");
$lines = explode("\n", $tree);

foreach ($prots as $rcpid => $pdata) unset($prots[$rcpid]['btree']);


$prots["TAAR5"]["btree"] = "001";
$prots["VN1R1"]["btree"] = "01";
// $prots["MS4A1"]["btree"] = "1";

if (inhouse)
{
    $prots["OR2W1"]["btree"] = "0000";
    $prots["OR52D1"]["btree"] = "0001";
}
else foreach ($lines as $lno => $ln)
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
            $prots[$rcpid]['btree'] = "000$binary";
            echo "$rcpid\n";
        }
    }
}

foreach ($prots as $rcpid => $pdata)
{
    // if (!inhouse && substr($rcpid, 0, 2) != "OR") continue;
    if (!isset($pdata["region"]["TMR3"]["start"])) continue;
    if (!isset($pdata['btree']))
    {
        $len = intval($pdata["region"]["TMR7"]["end"]) - intval($pdata["region"]["TMR3"]["start"]);
        $seq1 = substr($pdata['sequence'], intval($pdata["region"]["TMR3"]["start"]), $len);
        $best_match = false;
        $best_btree = false;
        $best_distance = 10000;
        foreach ($prots as $rcp2 => $p2)
        {
            if (!isset($p2['btree'])) continue;

            $len = intval($p2["region"]["TMR7"]["end"]) - intval($p2["region"]["TMR3"]["start"]);
            $seq2 = substr($p2['sequence'], intval($p2["region"]["TMR3"]["start"]), $len);
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
                $len = intval($p3["region"]["TMR7"]["end"]) - intval($p3["region"]["TMR3"]["start"]);
                $seq3 = substr($p3['sequence'], intval($p3["region"]["TMR3"]["start"]), $len);
                $lev = long_levenshtein($seq2, $seq3);
                if ($lev < $best_distance)
                {
                    foreach (str_split($p3['btree']) as $i => $c)
                    {
                        if ($i >= $n) break;
                        if ($c != substr($mrca, $i, 1))
                        {
                            $mrca = substr($mrca, 0, $i-1);
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
