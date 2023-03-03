<?php

chdir(__DIR__);
require_once("protutils.php");
require_once("odorutils.php");
require_once("statistics.php");

define('_EC50MUL', 1.3);

function set_color($r, $g, $b)
{
    $lr =intval( max(0,min(5,$r*6)));
    $lg = intval(max(0,min(5,$g*6)));
    $lb = intval(max(0,min(5,$b*6)));

    $ccode = intval(16 + $lb + 6*$lg + 36*$lr);

    echo "\x1b[48;5;{$ccode}m";
}

function fire_color($f)
{
    $k = $f * 4.2;
    $blue  = max(0, sin(($k+0.17)*1.15));
    $red   = pow(max(0, sin($k-1.0)), 0.666);
    $green = pow(max(0, sin($k-2.6)), 0.666);
    $red = max($red, $green);

    if ($k > 6) $red = $green = $blue = 1;
    else if ($k > 3.5)
    {
        $w     = max(0, sin(($k-3.5)*2));
        $red   = max($red, $w);
        $green = max($green, $w);
        $blue  = max($blue, $w);
    }

    set_color($red, $green, $blue);
}

function clear_color()
{
    echo "\x1b[49m";
}

function corrrsort($a, $b)
{
    if (abs($a) == abs($b)) return 0;
    return (abs($a) > abs($b)) ? -1 : 1;
}

foreach (@$argv as $a)
{
	$a = explode('=',$a,2);
	$_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}



$version = filemtime("method_rock6.php");

$json_file = "dock_results_rock6.json";
$dock_data = [];
if (file_exists($json_file))
{
    $dock_data = json_decode(file_get_contents($json_file), true);
}

$plotsz = 20; $plot2 = $plotsz / 2;


$utd = 0;
$ttl = 0;
if ($dock_data && count($dock_data))
{
    foreach ($dock_data as $dd) foreach ($dd as $d)
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

foreach ($dock_data as $rcp => $ligs)
{
    foreach ($ligs as $ligand => $pair)
    {
        if (@$pair['version'] < $version) continue;
        $idx = "$rcp|$ligand";
        $gidx = "All|$ligand";

        switch (@$pair['Actual'])
        {
            case 'Agonist':
            $emp = best_empirical_pair($rcp, $ligand, true);
            $xvals["is_ag"][$rcp][$idx] = 1;
            if (isset($emp['ec50'])) $xvals["ec50"][$rcp][$idx] = floatval($emp['ec50']);
            if (isset($emp['adjusted_curve_top'])) $xvals["top"][$rcp][$idx] = floatval($emp['adjusted_curve_top']);
            $xvals["is_ag"]["All"][$idx] = 1;
            if (isset($emp['ec50'])) $xvals["ec50"]["All"][$idx] = floatval($emp['ec50']);
            if (isset($emp['adjusted_curve_top'])) $xvals["top"]["All"][$idx] = floatval($emp['adjusted_curve_top']);
            break;

            case 'Non-Agonist':
            $xvals["is_ag"][$rcp][$idx] = 0;
            $xvals["top"][$rcp][$idx] = 0;
            $xvals["is_ag"]["All"][$idx] = 0;
            $xvals["top"]["All"][$idx] = 0;
            break;

            case 'Inverse Agonist':
            $xvals["is_ag"][$rcp][$idx] = 0; // -1;
            $xvals["top"][$rcp][$idx] = floatval($emp['adjusted_curve_top']) ?: 0;
            if (isset($emp['ec50'])) $xvals["ec50"][$rcp][$idx] = floatval($emp['ec50']);
            $xvals["is_ag"]["All"][$idx] = -1;
            $xvals["top"]["All"][$idx] = floatval($emp['adjusted_curve_top']) ?: 0;
            if (isset($emp['ec50'])) $xvals["ec50"]["All"][$idx] = floatval($emp['ec50']);
            break;

            default:
            goto _skip_pair;
        }

        foreach ($pair as $k => $v)
        {
            $jk = explode("_", $k);
            if (count($jk) > 1)
            {
                $ia = (strtolower($jk[0]) == 'active') ? "a" : "i";
                $k1 = $jk[1];

                /*if (substr($k1, 0, 7) == "BEnerg ")
                {
                    $bw = substr($k1, 7);
                    $yvals[$rcp]["$ia.$bw.e"][$idx] = $v;
                }*/
            
                if (substr($k1, 0, 5) == "Node ")
                {
                    $yvals[$rcp]["$ia.$k1"][$idx] = $v;
                }
            }
            
            if (substr($k, 0, 7) == "active_")
            {
                $k1 = "in$k";
                if (!isset($pair[$k1])) continue;
                $v1 = floatval($v) + floatval($pair[$k1]);
                $v2 = floatval($v) - floatval($pair[$k1]);
                $k2 = substr($k, 7);
            }
            else continue;

            /*if (substr($k2, 0, 7) == "BEnerg ")
            {
                $bw = substr($k2, 7);
                $yvals[$rcp]["$bw.e"][$idx] = $v2;
            }*/
            
            /*if (substr($k2, 0, 7) == "vdWrpl ")
            {
                $bw = substr($k2, 7);
                $yvals[$rcp]["$bw.v"][$idx] = $v2;
            }*/
            
            /*if (substr($k2, 0, 3) == "TMR")
            {
                $bw = substr($k2, 7);
                $yvals[$rcp][$k2][$idx] = $v;
                $yvals["All"][$k2][$idx] = $v2;
            }*/
            
            if (substr($k2, 0, 5) == "Node ")
            {
                $yvals[$rcp]["s.$k2"][$idx] = $v1;
                $yvals[$rcp]["d.$k2"][$idx] = $v2;
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
$bestcorr = 0.0;
foreach ($yvals as $rcp => $yv)
{
    foreach ($xvals as $xmet => $lxv)
    {
        $cxv = count($lxv[$rcp]);
        $threshold = max(0.5 * $cxv, 5);
        if ($rcp == "All") $threshold = 0.1;

        foreach ($yv as $metric => $ly)
        {
            $x = [];
            $y = [];
            foreach ($lxv[$rcp] as $idx => $lx)
            {
                if (isset($ly[$idx]))
                {
                    $x[] = $lx;
                    $y[] = $ly[$idx];
                }
            }

            // echo "Count x: " . count($x) . "; count y: " . count($y) . "\n";
            if ((count($x) >= $threshold && count($y) >= $threshold && (max($y) - min($y)) > 0.1 ) || $metric == 'acv.d')
            {
                // if ($metric == "7.46.y") print_r($y);
                $corr = round(correlationCoefficient($x, $y), 3);
                echo "$rcp: $xmet - $metric correlation: $corr\n";
                if (abs($corr) > abs($bestcorr)) $bestcorr = $corr;
                $p = (count($x) >= 20) ? calculate_p($x, $y, $corr, 100) : 0;
                if (($p <= 0.05 && abs($corr) > 0.25) || $metric == 'acv.d') $corrs[$rcp][$metric] = round($corr, 3);
            }
        }
    }
}

// print_r($corrs);

if (count($yvals))  echo "Best single-metric correlation: $bestcorr.\n\n";
else                echo "Insufficient data for correlation processing.\n\n";

exit;

if (count($corrs))
foreach ($corrs as $rcp => $c)
{
    foreach ($xvals as $xmet => $lxv)
    {
        $cxv = count($lxv[$rcp]);
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
                    foreach ($lxv[$rcp] as $idx => $lx)
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

            foreach ($y as $idx => $ly)
            {
                $x[$idx] = $lxv[$rcp][$idx];
                $yvals[$rcp][$metstr][$idx] = $y[$idx];
            }

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
            if (($corr < 0.5*$maxallc) && $metric != 'acv.d') continue;
            $j++;
            if (($j > 10) && $metric != 'acv.d') continue;
            echo str_pad($metric, 10);

            if (false===strpos($metric, " + ")
                &&
                false===strpos($metric, " - ")
                &&
                preg_match("/^[0-9]+[.][0-9]+/", $metric)
               )
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

        $xk = array_unique($lxv[$rcp]);
        rsort($xk);

        $maxyv = max( 
             max($yvals[$rcp][$idx]),
            -min($yvals[$rcp][$idx])
                    );
        $matrix = [];
        for ($yv=0; $yv<=$plotsz; $yv++) for ($xv=0; $xv<=$plotsz; $xv++) $matrix[$yv][$xv] = 0;

        $xmin = min($lxv[$rcp]);
        $xmax = max($lxv[$rcp]);
        $xscl = floatval($plotsz) / (($xmax-$xmin)?:1);
        $yscl = floatval($plot2) / ($maxyv ?: 1);

        foreach ($lxv[$rcp] as $oid => $xv)
        {
            $x = ($xv-$xmin)*$xscl;
            if (isset($yvals[$rcp][$idx][$oid]))
            {
                $y = $yvals[$rcp][$idx][$oid]*$yscl+$plot2;
                $matrix[$y][$x]++;
            }
        }

        echo "Plot:\n";
        foreach ($matrix as $y => $ln)
        {
            foreach ($ln as $x => $k)
            {
                fire_color($k);            
                echo ($y==$plot2?'-':" ");
                clear_color();
            }
            echo "\n";
        }
        echo "\n";
    }
}

_nodata:
;

echo "Processes running:\n";
passthru("ps -ef | grep ':[0-9][0-9] bin/primarydock' | grep -v grep");
