<?php

chdir(__DIR__);
require_once("statistics.php");
require_once("protutils.php");

$inpfile = "dock_results_clashmap.json";
$rplstrength = 3;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


$data = json_decode(file_get_contents($inpfile), true);

$array = [];

foreach ($data as $rcp => $rdat)
{
    if (!isset($array[$rcp])) $array[$rcp] = [];

    foreach ($rdat as $ligand => $ldat)
    {
    	$sign = (@$ldat['Actual'] == 'Agonist') ? 1 : -1;
        foreach ($ldat as $metric => $mdat)
        {
            $valmult = 0;

            if (substr($metric, 0, 7) == "BEnerg ")
            {
            	$valmult = 1;
            	$mdat = floatval($mdat);
            	if ($mdat > 0) $mdat = pow($mdat, 1.0/3);
        	}
            if (substr($metric, 0, 7) == "vdWrpl ") $valmult = $rplstrength;

            if ($valmult)
            {
                $pettia = explode(" ", $metric);
                if (count($pettia) < 2) continue;
                $morceaux = explode(".", $pettia[1]);
                if (count($morceaux) < 2) continue;

                $tmrno = intval($morceaux[0]);
                $bw = intval($morceaux[1]);
                if (!$tmrno || !$bw) continue;

                if (!isset($array[$rcp][$tmrno][$bw])) $array[$rcp][$tmrno][$bw] = floatval($mdat) * $valmult * $sign;
                else $array[$rcp][$tmrno][$bw] += floatval($mdat) * $valmult * $sign;

                if (!isset($array['all'][$tmrno][$bw])) $array['all'][$tmrno][$bw] = floatval($mdat) * $valmult * $sign;
                else $array['all'][$tmrno][$bw] += floatval($mdat) * $valmult * $sign;
            }
        }
    }
}

foreach ($array as $rcp => $rdat) foreach ($rdat as $tmrno => $tdat) ksort($array[$rcp][$tmrno]);

echo "<pre>".print_r($array, true);


// Get xyz coordinates for OR1A1 to use as a model.

$orid = "OR1A1";
$fam = family_from_protid($orid);
$c = file_get_contents("../pdbs/$fam/$orid.upright.pdb");
$lines = explode("\n", $c);
$resxyz = [];

foreach ($lines as $ln)
{
    if (substr($ln, 0, 7) == "ATOM   ")
    {
        $ln = preg_replace("/\\s+/", " ", $ln);
        $ln = explode(" ", $ln);
        if ($ln[2] == "CA")
        {
            $resno = intval($ln[4]);
            $x = floatval($ln[5]);
            $y = floatval($ln[6]);
            $z = floatval($ln[7]);

            $resxyz[$orid][$resno] = [$x,$y,$z];
        }
    }
}


// Set up an output image.

$pad = 25;
$w = 360 * 3 + $pad * 2;
$h = $w / 1.618;
$im = imagecreatetruecolor($w, $h);

$black = imagecolorallocate($im, 0,0,0);
$tm1col = imagecolorallocate($im, 64,64,64);
$tm2col = imagecolorallocate($im, 0,128,255);
$tm3col = imagecolorallocate($im, 0,255,128);
$tm4col = imagecolorallocate($im, 128,192,0);
$tm5col = imagecolorallocate($im, 255,192,0);
$tm6col = imagecolorallocate($im, 255,104,0);
$tm7col = imagecolorallocate($im, 225,0,32);

imagefilledrectangle($im, 0,0, $w-1,$h-1, $black);


// Each bw number, get the OR1A1 resno, get its xyz coords, plot its azimuth and y value,
// and draw a filled circle of the TMR color with a size proportional to the perturbation.
$maxptbn = 0.0;
foreach ($array['all'] as $tmrno => $tdat)
	foreach ($tdat as $bw => $perturbation)
		if ($perturbation > $maxptbn) $maxptbn = $perturbation;

$ptbm = 53.81 / sqrt($maxptbn);
foreach ($array['all'] as $tmrno => $tdat)
{
    if ($tmrno < 1 || $tmrno > 7) continue;
    $ltmcol = "tm{$tmrno}col";
    $ltmcol = $$ltmcol;

    $vals = [];
    $avgpert = array_sum($tdat) / (count($tdat) ?: 1);

    $yavg = 0.0;
    $ydiv = 0;
    $azc = [];
    $yc = [];
    foreach ($tdat as $bw => $perturbation)
    {
        $resno = resno_from_bw($orid, "$tmrno.$bw");
        $xyz = $resxyz[$orid][$resno];
        $yavg += $xyz[1];
        $ydiv++;

        $azimuth = find_angle($xyz[0], -$xyz[2]) + 0.3 - pi();
        if ($tmrno != 7 && $azimuth < 0) $azimuth += pi()*2;
        $azc[$bw] = $azimuth;
        $yc[$bw] = $xyz[1];
    }

    $lrgr = regression_line($yc, $azc);

    if ($ydiv) $yavg /= $ydiv;

    foreach ($tdat as $bw => $perturbation)
    {
        $resno = resno_from_bw($orid, "$tmrno.$bw");
        $xyz = $resxyz[$orid][$resno];
        $azimuth = find_angle($xyz[0], -$xyz[2]) + 0.3 - pi();
        if ($azimuth < 0) $azimuth += pi()*2;
        $adjaz = $azimuth - (/*$lrgr[1] +*/ $lrgr[0]*$xyz[1]);

        $vals['az'][$bw] = -$adjaz;
        $vals['y'][$bw] = $xyz[1];
        $vals['azy'][$bw] = $adjaz * ($xyz[1] - $yavg);

        $x = intval(3 * (180.0 / pi()) * $azimuth + $pad);
        $y = ($h/2) + $xyz[1]*20 + $pad;
        
        if ($perturbation < 0)
        {
        	$r = sqrt(-$perturbation)*$ptbm;
        	imageellipse($im, $x,$y, $r,$r, $ltmcol);
        }
        else
        {
        	$r = sqrt($perturbation)*$ptbm;
        	imagefilledellipse($im, $x,$y, $r,$r, $ltmcol);
    	}
    }

    $rxform = $avgpert;
    $thxform = correlationCoefficient($tdat, $vals['az']);
    $throt8 = correlationCoefficient($tdat, $vals['y']);
    $rrot8 = correlationCoefficient($tdat, $vals['azy']);

    echo "TMR$tmrno radial transformation: $rxform\n";
    echo "TMR$tmrno angular transformation: $thxform\n";
    echo "TMR$tmrno centrifugal rocking: $throt8\n";
    echo "TMR$tmrno torsional rocking: $rrot8\n";
    echo "\n";
}

imagepng($im, "../tmp/clashmap.png");

