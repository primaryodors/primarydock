<?php

chdir(__DIR__);
require_once("statistics.php");

$json_file = "dock_results_scwhere.json";
$scw_data = [];
if (file_exists($json_file))
{
    $scw_data = json_decode(file_get_contents($json_file), true);
}

$xvals = [];
$yvals = [];

foreach ($scw_data as $rcp => $ligs)
{
    foreach ($ligs as $ligand => $pair)
    {
        $idx = "$rcp|$ligand";

        $kJmol = @floatval($pair['Node 0']) ?: 0;

        switch (@$pair['Actual'])
        {
            case 'Agonist':
            $xvals[$idx] = 1;
            break;

            case 'Non-Agonist':
            $xvals[$idx] = 0;
            break;

            case 'Inverse Agonist':
            $xvals[$idx] = -1;
            break;

            default:
            goto _skip_pair;
        }

        foreach ($pair as $k => $v)
        {
            if (substr($k, 0, 4) == "SCW ")
            {
                $bw = substr($k, 4);
                $yvals["$bw.x"][$idx] = $v[0] * $kJmol;
                $yvals["$bw.y"][$idx] = $v[1] * $kJmol;
                $yvals["$bw.z"][$idx] = $v[2] * $kJmol;
            }
        }

        _skip_pair:
        ;
    }
}

set_time_limit(600);

$corrs = [];
foreach ($yvals as $metric => $ly)
{
    $x = [];
    $y = [];
    foreach ($xvals as $idx => $lx)
    {
        if (isset($ly[$idx]))
        {
            $x[] = $lx;
            $y[] = $ly[$idx];
        }
    }

    if (count($x) >= 10 && count($y) >= 10)
    {
        $corr = correlationCoefficient($x, $y);
        $p = calculate_p($x, $y, $corr, 100);
        if ($p <= 0.2 && abs($corr) > 0.5) $corrs[$metric] = $corr;
    }
}

function corrrsort($a, $b)
{
    if (abs($a) == abs($b)) return 0;
    return (abs($a) > abs($b)) ? -1 : 1;
}

uasort($corrs, 'corrrsort');

print_r($corrs);