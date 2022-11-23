<?php

function sgn($num)
{   
    if ($num < 0) return -1;
    if ($num > 0) return 1;
    return 0;
}

// http://www.sjsu.edu/faculty/gerstman/StatPrimer/correlation.pdf
function correlationCoefficient($xarr, $yarr)
{   
    if (!is_array($xarr) || !is_array($yarr))   throw new Exception("correlationCoefficient requires two equal length arrays.");
    if (count($xarr) != count($yarr))           throw new Exception("correlationCoefficient requires two equal length arrays.");
    if (!count($xarr) || !count($yarr))         throw new Exception("correlationCoefficient passed empty arrays.");
    
    $xarr = array_values($xarr);
    $yarr = array_values($yarr);
    
    foreach ($xarr as $i => $x)
    {   
        $y = $yarr[$i];
        if (!is_numeric($x)) $xarr[$i] = (trim($x) == '' ? -500 : 500);
        if (!is_numeric($y)) $yarr[$i] = (trim($y) == '' ? -500 : 500);
    }
    
    $xavg = array_sum($xarr) / count($xarr); 
    $yavg = array_sum($yarr) / count($yarr); 
    
    $ssxx = 0.0; $ssyy = 0.0; $ssxy = 0.0;
    
    foreach ($xarr as $i => $x)
    {   
        $y = $yarr[$i];
        /*if (!is_numeric($x)) $x = (trim($x) == '' ? -500 : 500);
        if (!is_numeric($y)) $y = (trim($y) == '' ? -500 : 500);*/
        
        $ssxx += pow($x - $xavg, 2);
        $ssyy += pow($y - $yavg, 2);
        $ssxy += ($x - $xavg)*($y - $yavg);
    }
    
    if ($ssxx == 0 || $ssyy == 0) return '(insufficient data)';
    
    $r = $ssxy/pow($ssxx*$ssyy,0.5);
    return $r;
}

function array_shuffle($arr)
{   
    $arr1 = [];
    foreach (array_values($arr) as $v) 
    {   
        do
        {   
            $k = rand(0,16777215);
        } while (isset($arr1[$k]));
        $arr1[$k] = $v;
    }
    ksort($arr1);
    return array_values($arr1);
}


function hill($concn, $emax, $ec50, $nH = 2.0)
{   
    return $emax / (1.0 + pow($ec50 / $concn, $nH));
}

function scatter_to_best_fit_hill($arr, $iters = 50)       // Index = concentration; value = observed response.
{   
    $emin = min($arr);
    $emax = max($arr);
    
    foreach ($arr as &$a) $a -= $emin;
    $emax -= $emin;
    
    $ec50 = array_sum(array_keys($arr))/count($arr);
    $nH = 2.0;
    $delta = 9999999999.99999;
    
    if (!$emax) $ec50 = 0;
    else
    {
        $half = $iters/2;
        for ($i=0; $i<$iters; $i++)
        {   
            $delta = 0.0;
            foreach ($arr as $concn => $obs) $delta += abs($obs - hill(floatval($concn), $emax, $ec50, $nH));
            
            $ecold = $ec50;
            $ec50 *= pow(10, 0.01*rand(-100, 100));
            
            $deltan= 0.0;
            foreach ($arr as $concn => $obs) $deltan+= abs($obs - hill(floatval($concn), $emax, $ec50, $nH));
            
            if ($deltan > $delta) $ec50 = $ecold;
            else $delta = $deltan;
            
            if ($i >= $half && $nH >= 0.1 && $nH <= 10)
            {   $nHold = $nH;
                $nH *= pow(10, 0.0001*rand(-100, 100));
                
                $deltan= 0.0;
                foreach ($arr as $concn => $obs) $deltan+= abs($obs - hill(floatval($concn), $emax, $ec50, $nH));
                
                if ($deltan > $delta) $nH = $nHold;
                else $delta = $deltan;
            }
        }
    }
    
    return [$emax, $emin, $ec50, $nH, $delta];
}

function median($arr)
{   
    if (!is_array($arr)) return $arr;
    
    sort($arr);
    $half = 0.5 * count($arr);
    
    if ($half == floor($half))
    {   
        return 0.5*$arr[$half] + 0.5*$arr[$half-1];
    }
    else
    {   
        return $arr[floor($half)];
    }
}

function find_angle($dx, $dy)
{
    $d = sqrt(pow($dx,2)+pow($dy,2)); if ($d == 0) return 0;		// undefined
	$dx /= $d;
	$dy /= $d;
	$as = asin($dx);
	$ac = acos(-$dy);
	if ($as > 0) return $ac;
	else return pi()*2-$ac;
}

function sigmoid($x)
{
    // https://en.wikipedia.org/wiki/Sigmoid_function
    $emx = pow(M_E, -$x);
    return 1.0 / (1+$emx);
}























