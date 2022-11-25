<?php

$pwd = getcwd();
chdir(__DIR__);
require_once("statistics.php");
chdir($pwd);

function get_active_node($array)
{
    $activenode = 0;
    if (@$array['Active node']) $activenode = intval($array['Active node']);

    if (!$activenode) foreach (array_keys($array) as $key) if (strtolower(substr($key, 0, 5)) == "node ")
    {
        $nodeno = intval(explode(" ", $key)[1]);
        if ($nodeno > $activenode) $activenode = $nodeno;
    }

    return $activenode;
}

function path_dock_probability($array)
{
    $activenode = get_active_node($array);
    $pocketnode = $activenode - 1;

    $probability = 0;
    for ($n = 1; $n < $pocketnode && isset($array["Node $n"]); $n++)
    {
        if ($n == 1) $probability = 1.0;
        $m = $n - 1;
        $ratio = -floatval($array["Node $n"]) / -floatval($array["Node $m"]);
        $s = sigmoid(log($ratio));
        $probability *= $s;
    }

    return $probability;
}

function active_binding_ratio($array)
{
    $activenode = get_active_node($array);
    if ($activenode)
    {
        $pocketnode = $activenode - 1;
        $ratio = -floatval($array["Node $activenode"]) / -floatval($array["Node $pocketnode"]);
        return max($ratio, 0);
    }
    else return 0;
}

function evaluate_result($array)
{
    $activenode = get_active_node($array);
    
    $pocketnode = $activenode - 1;
    if ($pocketnode <= 0)
    {
        $pocketnode = 0;
        $activenode = $pocketnode + 1;

        if (-$array["Node $pocketnode"] > 2*@-$array["Node $activenode"])
        {
            $array['Prediction'] = "Inverse Agonist";
            return $array;
        }
        else if (@-$array["Node $activenode"] > 0.95*@-$array["Node $pocketnode"])
        {
            $array['Prediction'] =  "Agonist";
            return $array;
        }
        else
        {
            $array['Prediction'] = "Non-Agonist";
            return $array;
        }
    }

    $probability = path_dock_probability($array);
    $acvbndratio = active_binding_ratio($array);

    $capture = 0.0;
    $capqty = 0;
    foreach ($array as $key => $value) if (strtolower(substr($key, 0, 5)) == "node ")
    {
        $nodeno = intval(explode(" ", $key)[1]);
        if ($nodeno < $pocketnode)
        {
            $capture -= floatval($value);
            $capqty++;
        }
    }
    if ($capqty) $capture /= $capqty;

    if (!isset($array["Node $activenode"])) return "Non-Agonist";
    $completion = $array["Poses"] ? (floatval($array["Full poses"]) / $array["Poses"]) : 0;
    $cpl_weighted = $completion * (max(0.001, @-$array["Node $pocketnode"], -$array["Node $activenode"]) / $capture);
    $acvratio = -$array["Node $activenode"] / (-$array["Node $pocketnode"] ?: 0.001);

    // echo "$capture, $completion, $cpl_weighted\n";

    $prediction = "Non-Agonist";
    if ($capture >= 15 && (/*$completion >= 0.75 ||*/ $cpl_weighted >= 0.5))
    {
        if ($acvratio < 0.9) $prediction = "Inverse Agonist";
        else $prediction = "Agonist";
    }

    $array['Prediction'] = $prediction;
    $array['Probability'] = $probability;
    $array['Activation Binding Ratio'] = $acvbndratio;

    return $array;
}
