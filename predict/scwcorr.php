<?php

chdir(__DIR__);
require_once("protutils.php");
require_once("odorutils.php");
require_once("statistics.php");

function set_color($r, $g, $b)
{
    $lr =intval( max(0,min(5,$r*6)));
    $lg = intval(max(0,min(5,$g*6)));
    $lb = intval(max(0,min(5,$b*6)));

    $ccode = intval(16 + $lb + 6*$lg + 36*$lr);

    echo "\x1b[48;5;{$ccode}m";
}

function clear_color()
{
    echo "\x1b[49m";
}

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
            $emp = best_empirical_pair($rcp, $ligand, true);
            // print_r($emp); exit;
            if ($emp)
            {
                if (isset($emp['adjusted_curve_top']) && isset($emp['ec50']))
                    $x = (floatval($emp['adjusted_curve_top']) - min(10, -2.0*floatval($emp['ec50']))) / 2;
                else if (isset($emp['ec50']))
                    $x = min(10, -2.0*floatval($emp['ec50']));
                else if (isset($emp['adjusted_curve_top']))
                    $x = floatval($emp['adjusted_curve_top']);
                else $x = 1;
            }
            else $x = 1;
            $xvals[$rcp][$idx] = $x;
            break;

            case 'Non-Agonist':
            $xvals[$rcp][$idx] = 0;
            break;

            case 'Inverse Agonist':
            $xvals[$rcp][$idx] = -1;
            break;

            default:
            goto _skip_pair;
        }

        // print_r($pair);

        foreach ($pair as $k => $v)
        {
            if (substr($k, 0, 7) == "BEnerg ")
            {
                $bw = substr($k, 7);
                $yvals[$rcp]["$bw.e"][$idx] = $v;
            }
            if (substr($k, 0, 7) == "vdWrpl ")
            {
                $bw = substr($k, 7);
                $yvals[$rcp]["$bw.v"][$idx] = $v;
            }
            if (substr($k, 0, 4) == "SCW ")
            {
                $bw = substr($k, 4);
                $yvals[$rcp]["$bw.x"][$idx] = $v[0] * $kJmol;
                $yvals[$rcp]["$bw.y"][$idx] = $v[1] * $kJmol;
                $yvals[$rcp]["$bw.z"][$idx] = $v[2] * $kJmol;
            }
        }

        _skip_pair:
        ;
    }
}

// echo "yvals: "; print_r($yvals);
// echo "xvals: "; print_r($xvals);

set_time_limit(600);

$corrs = [];
foreach ($yvals as $rcp => $yv)
{   
    $cxv = count($xvals[$rcp]);
    $threshold = 0.75 * $cxv;

    foreach ($yv as $metric => $ly)
    {
        $x = [];
        $y = [];
        foreach ($xvals[$rcp] as $idx => $lx)
        {
            if (isset($ly[$idx]))
            {
                $x[] = $lx;
                $y[] = $ly[$idx];
            }
        }

        // echo "Count x: " . count($x) . "; count y: " . count($y) . "\n";
        if (count($x) >= $threshold && count($y) >= $threshold && (max($y) - min($y)) > 0.1 )
        {
            // if ($metric == "7.46.y") print_r($y);
            $corr = correlationCoefficient($x, $y);
            echo "Correlation: $corr\n";
            $p = (count($x) >= 20) ? calculate_p($x, $y, $corr, 100) : 0;
            if ($p <= 0.05 && abs($corr) > 0.25) $corrs[$rcp][$metric] = round($corr, 3);
        }
    }
}

function corrrsort($a, $b)
{
    if (abs($a) == abs($b)) return 0;
    return (abs($a) > abs($b)) ? -1 : 1;
}

if (!count($corrs)) echo "Insufficient data for correlation processing.\n\n";
else
foreach ($corrs as $rcp => $c)
{
    $cxv = count($xvals[$rcp]);
    if ($cxv < 5) continue;

    uasort($c, 'corrrsort');

    $cc = count($c);
    $maxnatc = $cc ? (max(max($c), -min($c))) : 0.5;

    for ($bits=0; $bits < 16384; $bits++)
    {
        if ($bits >= pow(2, $cc)) break;
        $x = [];
        $y = [];
        $metstr = "";
        for ($i=0; $i<$cc; $i++)
        {
            $sci = sgn(array_values($c)[$i]);
            $pi = pow(2, $i);
            if ($pi > $bits) break;
            if ($bits & $pi)
            {
                $metric = array_keys($c)[$i];
                $ly = $yvals[$rcp][$metric];
                foreach ($xvals[$rcp] as $idx => $lx)
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

        foreach ($y as $idx => $ly) $x[$idx] = $xvals[$rcp][$idx];

        if (count($x) >= 10 && count($y) >= 10)
        {
            $corr = correlationCoefficient($x, $y);
            // $p = calculate_p($x, $y, $corr, 100);
            if (/* $p <= 0.1 &&*/ abs($corr) > max($c)) $c[$metstr] = round($corr, 3);

            if (count($c) > 40) break;
        }
    }

    uasort($c, 'corrrsort');
    $cc = count($c);
    $maxallc = $cc ? (max(max($c), -min($c))) : 0.5;

    echo "Correlations (residues of $rcp, $cxv ligands):\n";
    $j = 0;
    foreach ($c as $metric => $corr)
    {
        if ($corr < 0.5*$maxallc) continue;
        $j++;
        if ($j > 5) break;
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

    $idx = array_keys($c)[0];
    if (!@$yvals[$rcp][$idx]) continue;

    $xk = array_unique($xvals[$rcp]);
    rsort($xk);

    $maxyv = max( 
        max($yvals[$rcp][$idx]),
        -min($yvals[$rcp][$idx])
                );
    $matrix = [];
    /*foreach ($xvals[$rcp] as $oid => $xv)
    {
        foreach ($xk as $xv) $matrix[$xv][$oid] = 0;
    }*/
    for ($y=0; $y<=10; $y++) foreach ($xk as $xv) $matrix[$xv][$y] = 0;
    foreach ($xvals[$rcp] as $oid => $xv)
    {
        if (isset($yvals[$rcp][$idx][$oid]))
        {
            // $matrix[$xv][$oid] = $yvals[$rcp][$idx][$oid];
            $y = intval(abs($yvals[$rcp][$idx][$oid]) * 10.0/$maxyv);
            $matrix[$xv][$y]++;
        }
    }

    $maxyv = 0;
    // for ($xv = -1; $xv <= 1; $xv++)
    foreach ($xk as $xv)
    {
        if (@$matrix[$xv])
        {
            $m = max(max($matrix[$xv]), -min($matrix[$xv]));
            if ($m > $maxyv) $maxyv = $m;
        }
    }

    // for ($xv = 1; $xv >= -1; $xv--)
    foreach ($xk as $xv)
    {
        $j = 0;
        if (@$matrix[$xv])
        foreach (array_values($matrix[$xv]) as $i)
        {
            $k = abs($i);
            if ($maxyv) $k /= $maxyv;
            // echo intval($k*9.999)."";
            $red = $green = $blue = $k;

            $k *= M_PI*2;
            $blue  = max(0, sin($k+0.17));
            $red   = max(0, sin($k-1));
            $green = max(0, sin($k-2));
            $red = max($red, $green);

            set_color($red, $green, $blue);
            echo " ";
            clear_color();
            $j++;
            if ($j > 99) break;
        }
        echo "\n";
    }
    echo "\n";
}

_nodata:
;

echo "Processes running:\n";
passthru("ps -ef | grep ':[0-9][0-9] bin/primarydock' | grep -v grep");