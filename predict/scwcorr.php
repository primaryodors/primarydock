<?php

chdir(__DIR__);
require_once("statistics.php");
$version = filemtime("method_scwhere.php");

$json_file = "dock_results_scwhere.json";
$scw_data = [];
if (file_exists($json_file))
{
    $scw_data = json_decode(file_get_contents($json_file), true);
}

$utd = 0;
if (count($scw_data))
{
    foreach ($scw_data as $d) if (@$d['version'] == $version) $utd++;
    echo round(100.0*$utd/count($scw_data), 2) . "% up to date.\n\n";
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
        if ($p <= 0.05 && abs($corr) > 0.25) $corrs[$metric] = round($corr, 3);
    }
}

function corrrsort($a, $b)
{
    if (abs($a) == abs($b)) return 0;
    return (abs($a) > abs($b)) ? -1 : 1;
}

uasort($corrs, 'corrrsort');

$maxnatc = max(max($corrs), -min($corrs));
$cc = count($corrs);

for ($bits=0; $bits < 4096; $bits++)
{
    if ($bits >= pow(2, $cc)) break;
    $x = [];
    $y = [];
    $metstr = "";
    for ($i=0; $i<$cc; $i++)
    {
        $sci = sgn(array_values($corrs)[$i]);
        $pi = pow(2, $i);
        if ($pi > $bits) break;
        if ($bits & $pi)
        {
            $metric = array_keys($corrs)[$i];
            $ly = $yvals[$metric];
            foreach ($xvals as $idx => $lx)
            {
                if (isset($ly[$idx]))
                {
                    if (!isset($y[$idx])) $y[$idx] = 0.0;
                    $y[$idx] += $ly[$idx] * $sci;
                }
            }
            if ($metstr) $metstr .= ($sci < 0 ? " - " : " + ");
            $metstr .= $metric;
        }
    }

    foreach ($y as $idx => $ly) $x[$idx] = $xvals[$idx];

    if (count($x) >= 10 && count($y) >= 10)
    {
        $corr = correlationCoefficient($x, $y);
        // $p = calculate_p($x, $y, $corr, 100);
        if (/* $p <= 0.1 &&*/ abs($corr) > $maxnatc) $corrs[$metstr] = round($corr, 3);
    }
}

echo "Correlations: ";
print_r($corrs);
echo "\n\n";

passthru("ps -ef | grep ':[0-9][0-9] bin/primarydock' | grep -v grep");