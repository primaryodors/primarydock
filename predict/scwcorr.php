<?php

chdir(__DIR__);
require_once("protutils.php");
require_once("statistics.php");
$version = filemtime("method_scwhere.php");

$json_file = "dock_results_scwhere.json";
$scw_data = [];
if (file_exists($json_file))
{
    $scw_data = json_decode(file_get_contents($json_file), true);
}

$utd = 0;
$ttl = 0;
if ($scw_data && count($scw_data))
{
    foreach ($scw_data as $scw) foreach ($scw as $d)
    {
        $ttl++;
        if (@$d['version'] == $version) $utd++;
    }
    if (!$ttl) $ttl = 1;
    echo round(100.0*$utd/$ttl, 2) . "% ($utd of $ttl) up to date with version $version.\n\n";
}
else goto _nodata;

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

$cxv = count($xvals);
$threshold = 0.75 * $cxv;
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

    if (count($x) >= $threshold && count($y) >= $threshold && (max($y) - min($y)) > 0.1 )
    {
        // if ($metric == "7.46.y") print_r($y);
        $corr = correlationCoefficient($x, $y);
        $p = (count($x) >= 20) ? calculate_p($x, $y, $corr, 100) : 0;
        if ($p <= 0.05 && abs($corr) > 0.5) $corrs[$metric] = round($corr, 3);
    }
}

function corrrsort($a, $b)
{
    if (abs($a) == abs($b)) return 0;
    return (abs($a) > abs($b)) ? -1 : 1;
}

uasort($corrs, 'corrrsort');

$cc = count($corrs);
$maxnatc = $cc ? (max(max($corrs), -min($corrs))) : 0.5;

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
        if (/* $p <= 0.1 &&*/ abs($corr) > max($corrs)) $corrs[$metstr] = round($corr, 3);

        if (count($corrs) > 40) break;
    }
}

echo "Correlations (residues of $rcp):\n";
// print_r($corrs);
foreach ($corrs as $metric => $corr)
{
    echo str_pad($metric, 10);

    if (false===strpos($metric, " + ") && false===strpos($metric, " - "))
    {
        $pettia = explode('.', $metric);
        $resno = resno_from_bw($rcp, "{$pettia[0]}.{$pettia[1]}");
        $aa = substr($prots[$rcp]['sequence'], $resno-1, 1);
        echo str_pad("$aa$resno", 7);
    }

    if ($corr >= 0) echo ' ';
    echo $corr;
    echo "\n";
}
echo "\n\n";

_nodata:

passthru("ps -ef | grep ':[0-9][0-9] bin/primarydock' | grep -v grep");