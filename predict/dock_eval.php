<?php

function evaluate_result($array, $pocketnode)
{
    $activenode = $pocketnode + 1;
    $capture = max(-$array["Node 0"], -$array["Node 1"]/*, -$average["Node 2"]*/);
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