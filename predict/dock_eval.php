<?php

function evaluate_result($array)
{
    $activenode = 0;
    foreach (array_keys($array) as $key) if (strtolower(substr($key, 0, 5)) == "node ")
    {
        $nodeno = intval(explode(" ", $key)[1]);
        if ($nodeno > $activenode) $activenode = $nodeno;
    }
    $pocketnode = $activenode - 1;
    if ($pocketnode < 0) return "(unknown)";

    $capture = 0.0;
    foreach ($array as $key => $value) if (strtolower(substr($key, 0, 5)) == "node ")
    {
        $nodeno = intval(explode(" ", $key)[1]);
        if ($nodeno < $pocketnode) $capture -= floatval($value);
    }

    $completion = floatval($array["Full poses"]) / $array["Poses"];
    $cpl_weighted = $completion * (max(0.001, @-$array["Node $pocketnode"], -$array["Node $activenode"]) / $capture);
    $acvratio = -$array["Node $activenode"] / (-$array["Node $pocketnode"] ?: 0.001);

    $prediction = "Non-Agonist";
    if ($capture >= 20 && (/*$completion >= 0.75 ||*/ $cpl_weighted >= 0.75))
    {
        if ($acvratio < 0.75) $prediction = "Inverse Agonist";
        else $prediction = "Agonist";
    }

    return $prediction;
}