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

function is_priority($metric)
{
    if ($metric == "4<->6") return true;
    if (substr($metric, 0, 8) == "AcvTheta") return true;
    // if (substr($metric, 0, 3) == "TMR") return true;
    return false;
}

$version = filemtime("method_softrock.php");

$json_file = "dock_results_softrock.json";
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

$a = []; $na = []; $va = [];

foreach ($dock_data as $rcp => $ligs)
{
    $a[$rcp] = $na[$rcp] = $va[$rcp] = 0;
    
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
            
            $a[$rcp]++;
            break;

            case 'Non-Agonist':
            $xvals["is_ag"][$rcp][$idx] = 0;
            $xvals["top"][$rcp][$idx] = 0;
            $xvals["is_ag"]["All"][$idx] = 0;
            $xvals["top"]["All"][$idx] = 0;
            $na[$rcp]++;
            break;

            case 'Inverse Agonist':
            $xvals["is_ag"][$rcp][$idx] = 0; // -1;
            $xvals["top"][$rcp][$idx] = floatval($emp['adjusted_curve_top']) ?: 0;
            if (isset($emp['ec50'])) $xvals["ec50"][$rcp][$idx] = floatval($emp['ec50']);
            $xvals["is_ag"]["All"][$idx] = -1;
            $xvals["top"]["All"][$idx] = floatval($emp['adjusted_curve_top']) ?: 0;
            if (isset($emp['ec50'])) $xvals["ec50"]["All"][$idx] = floatval($emp['ec50']);
            $va[$rcp]++;
            break;


            default:
            goto _skip_pair;
        }

        foreach ($pair as $k => $v)
        {
            if (substr($k, 0, 6) == "BEnerg")
            {
                $bw = substr($k, 7);
                if (false !== strpos($bw, '.'))
                {
                    // $yvals[$rcp]["$bw.e"][$idx] = $v;
                
                    $tmr = intval($bw);
                    if (!isset($yvals[$rcp]["TMR$tmr.e"][$idx])) $yvals[$rcp]["TMR$tmr.e"][$idx] = $v;
                    else $yvals[$rcp]["TMR$tmr.e"][$idx] += $v;
                }
            }
            else
            if (substr($k, 0, 6) == "vdWrpl")
            {
                $bw = substr($k, 7);
                // $yvals[$rcp]["$bw.v"][$idx] = $v;
            }
            else
            if (substr($k, 0, 4) == "Node "
                ||
                substr($k, 0, 6) == "PolSat"
                ||
                substr($k, 0, 8) == "AcvTheta"
               )
            {
                $bw = substr($k, 7);
                $yvals[$rcp][$k][$idx] = $v;
                $yvals["All"][$k][$idx] = $v;
            }
        }

        _skip_pair:
        ;
    }
    
    foreach ($ligs as $ligand => $pair)
    {
        if (@$pair['version'] < $version) continue;
        $idx = "$rcp|$ligand";
        $yvals[$rcp]["4<->6"][$idx] = (@$yvals[$rcp]["TMR4.e"][$idx] ?: 0) * (@$yvals[$rcp]["TMR6.e"][$idx] ?: 0);
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
                else
                {
                    $x[] = $lx;
                    $y[] = 0;
                }
            }

            // echo "Count x: " . count($x) . "; count y: " . count($y) . "\n";
            if ((count($x) >= $threshold && count($y) >= $threshold && (max($y) - min($y)) > 0.1 ) || is_priority($metric))
            {
                // if ($metric == "7.46.y") print_r($y);
                $corr = round(correlationCoefficient($x, $y), 3);
                // echo "Correlation: $corr\n";
                if (abs($corr) > abs($bestcorr)) $bestcorr = $corr;
                $p = (count($x) >= 20) ? calculate_p($x, $y, $corr, 100) : 0;
                if (($p <= 0.05 && abs($corr) > 0.25) || is_priority($metric)) $corrs[$rcp][$metric] = round($corr, 3);
            }
        }
    }
}

// print_r($corrs);

if (count($yvals))  echo "Best single-metric correlation: $bestcorr.\n\n";
else                echo "Insufficient data for correlation processing.\n\n";

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

        echo "Correlations ($rcp:$xmet, $cxv ligands):\n";
        $j = 0;
        foreach ($c as $metric => $corr)
        {
            if ((abs($corr) < 0.5*$maxallc) && !is_priority($metric)) continue;
            $j++;
            if (($j > 10) && !is_priority($metric)) continue;
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
        
        // print_r($c);

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
        $mxmtrx = 0.0;

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
                if ($matrix[$y][$x] > $mxmtrx) $mxmtrx = $matrix[$y][$x];
            }
        }

        echo "Plot ($rcp: $xmet vs. $idx):\n";
        foreach ($matrix as $y => $ln)
        {
            foreach ($ln as $x => $k)
            {
                fire_color(floatval($k) / $mxmtrx);            
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
