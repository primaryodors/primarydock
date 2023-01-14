<?php

$inpfile = "dock_results_clashmap.json";
$rplstrength = 100;

chdir(__DIR__);

$data = json_decode(file_get_contents($inpfile), true);

$array = [];

foreach ($data as $rcp => $rdat)
{
    if (!isset($array[$rcp])) $array[$rcp] = [];

    foreach ($rdat as $ligand => $ldat)
        foreach ($ldat as $metric => $mdat)
        {
            $valmult = 0;

            if (substr($metric, 0, 7) == "BEnerg ") $valmult = 1;
            if (substr($metric, 0, 7) == "vdWrpl ") $valmult = $rplstrength;

            if ($valmult)
            {
                $pettia = explode(" ", $metric);
                if (count($pettia) < 2) continue;
                $morceaux = explode(".", $pettia[1]);
                if (count($morceaux) < 2) continue;

                $tmrno = intval($morceaux[0]);
                $bw = intval($morceaux[1]);
                if (!$tmrno || !$bw) continue;

                if (!isset($array[$rcp][$tmrno][$bw])) $array[$rcp][$tmrno][$bw] = floatval($mdat) * $valmult;
                else $array[$rcp][$tmrno][$bw] += floatval($mdat) * $valmult;
            }
        }
}

foreach ($array as $rcp => $rdat) foreach ($rdat as $tmrno => $tdat) ksort($array[$rcp][$tmrno]);

echo "<pre>".print_r($array, true);

