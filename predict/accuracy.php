<?php

$file = "dock_results_coupled.json";

chdir(__DIR__);
$data = json_decode(file_get_contents($file), true);

$nc = [];
$fp = [];
$fn = [];

foreach ($data as $orid => $pairs)
{
    foreach ($pairs as $odor => $r)
    {
        if (!isset($r['Actual']) || !isset($r['Predicted'])) continue;
        
        if (strtolower($r['Actual']) == "(unknown)") continue;
        
        if (strtolower($r['Actual']) == strtolower($r['Predicted']))
        {
            if (!isset($nc['All'])) $nc['All'] = 1;
            else $nc['All']++;
            if (!isset($nc[$orid])) $nc[$orid] = 1;
            else $nc[$orid]++;
        }
        else if (strtolower($r['Predicted']) == "agonist")
        {
            if (!isset($fp['All'])) $fp['All'] = 1;
            else $fp['All']++;
            if (!isset($fp[$orid])) $fp[$orid] = 1;
            else $fp[$orid]++;
        }
        else
        {
            if (!isset($fn['All'])) $fn['All'] = 1;
            else $fn['All']++;
            if (!isset($fn[$orid])) $fn[$orid] = 1;
            else $fn[$orid]++;
        }
    }
}


foreach (array_keys($data) as $orid)
{
    echo "Results for $orid:\n";
    
    $lnc = (@$nc[$orid] ?: 0);
    $lfp = (@$fp[$orid] ?: 0);
    $lfn = (@$fn[$orid] ?: 0);
    
    $pcnt = round(floatval($lnc) / ($lnc + $lfp + $lfn) * 100, 2);
    
    echo "$lnc correct predictions ($pcnt%).\n";
    echo "  $lfp false positives.\n";
    echo "  $lfn false negatives.\n";
    echo "\n";
}


echo "Total:\n";
    
$lnc = (@$nc["All"] ?: 0);
$lfp = (@$fp["All"] ?: 0);
$lfn = (@$fn["All"] ?: 0);

$pcnt = round(floatval($lnc) / ($lnc + $lfp + $lfn) * 100, 2);

echo "$lnc correct predictions ($pcnt%).\n";
echo "  $lfp false positives.\n";
echo "  $lfn false negatives.\n";
echo "\n";













