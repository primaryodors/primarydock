<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");

function lum($r, $g, $b)
{
    $rcontrib = 0.299;
    $bcontrib = 0.114;
    $gcontrib = 1.0 - $rcontrib - $bcontrib;
    return $rcontrib*$r + $gcontrib*$g + $bcontrib*$b;
}

$bkcolor = [0x15, 0x1a, 0x37];

$bkcol_tetrapod = $bkcolor;
$bkcol_tetrapod[1] += 5;

$bkcol_fishlike = $bkcolor;
$bkcol_fishlike[2] += 13;

$bkcol_taar = $bkcolor;
$bkcol_taar[0] += 10;

$bklum = lum($bkcolor[0], $bkcolor[1], $bkcolor[2]);

$odor = find_odorant(@$_REQUEST['o']);
if (!$odor)
{
    header("Location: odorants.php");
    exit;
}

$mode = 'e';            // Empirical.
if (isset($_REQUEST["m"])) $mode = $_REQUEST["m"];

$t = [];
$e = [];
$p = [];

switch ($mode)
{
    case 'e':
    foreach ($odor['activity'] as $ref => $acv)
    {
        foreach ($acv as $rcpid => $a)
        {
            $act = @$a['adjusted_curve_top'] ?: false;
            if (!isset($t[$rcpid])) $t[$rcpid] = $act;
            else
            {
                if ($t[$rcpid] < $act) $t[$rcpid] = $act;
            }

            $ec50 = @$a['ec50'] ?: 0;
            if (!isset($e[$rcpid])) $e[$rcpid] = $ec50;
            else
            {
                if ($ec50)
                {
                    if (!$e[$rcpid] || $ec50 < $e[$rcpid]) $e[$rcpid] = $ec50;
                }
            }
        }
    }
    break;

    case 'p':
    $predictions = [];
    chdir(__DIR__);
    $dock_results = json_decode(file_get_contents("../predict/dock_results.json"), true);
    $odorname_under = str_replace(' ', '_', $odor["full_name"]);
    foreach ($dock_results as $protid => $dr)
    {
        if (isset($dr[$odorname_under]))
        {
            $p[$protid] = floatval($dr[$odorname_under]["DockScore"]);
        }
    }
    break;

    default:
    die("Unsupported mode $mode.\n");
}

$res = 3;
$xbuf = 80;
$ybuf = 80;

$w = count($prots)*$res + $xbuf;
$h = 300;

$im = imagecreatetruecolor($w, $h);
imagefilledrectangle($im, 0,0, $w,$h, imagecolorallocate($im,$bkcolor[0],$bkcolor[1],$bkcolor[2]));

$last4 = "";
$xmax = count($prots);
$dxmax = $xmax * $res + $xbuf/2;
foreach (array_keys($prots) as $x => $orid)
{
    $dx  = $x * $res + $xbuf/2;

    $first4 = substr($orid, 0, 4);

    if ($first4 != $last4)
    {
        switch ($first4)
        {
            case "OR1A":
            imagefilledrectangle($im, $dx,0, $dxmax,$h, imagecolorallocate($im,$bkcol_tetrapod[0],$bkcol_tetrapod[1],$bkcol_tetrapod[2]));
            break;

            case "OR51":
            imagefilledrectangle($im, $dx,0, $dxmax,$h, imagecolorallocate($im,$bkcol_fishlike[0],$bkcol_fishlike[1],$bkcol_fishlike[2]));
            break;

            case "TAAR":
            imagefilledrectangle($im, $dx,0, $dxmax,$h, imagecolorallocate($im,$bkcol_taar[0],$bkcol_taar[1],$bkcol_taar[2]));
            break;

            default:
            ;
        }
    }

    $last4 = $first4;
}

$maxt = count($t) ? ( @max($t) ?: 1 ) : 1;
$maxe = count($e) ? ( @max($e) ?: 0 ) : 0;
$mine = count($e) ? ( @min($e) ?: -6 ) : -6;
if ($maxe) $maxe += 0.5;
$mine -= 0.5;
$maxp = count($p) ? ( @max($p) ?: 1 ) : 1;

if ($maxt < 1) $maxt = 1;
if ($maxp < 1) $maxp = 1;
if ($maxp > 50) $maxp = 50;

if ($maxe <= $mine+2) { $maxe += 1; $mine -= 1; }

if ($mine <= -3 && $maxt < 2) $maxt = 2;

$tscale = floatval($h-$ybuf) / $maxt;
$escale = floatval($h-$ybuf) / ($maxe-$mine);
$pscale = floatval($h-$ybuf) / $maxp;

$red   = imagecolorallocate($im,255,96,80);
$wine  = imagecolorallocate($im,128,32,48);
$green = imagecolorallocate($im,64,144,96);
$brown = imagecolorallocate($im,96,80,64);
$yellow= imagecolorallocate($im,192,160,96);
$blue  = imagecolorallocate($im,128,160,255);
$cyan  = imagecolorallocate($im,32,80,104);
$pink  = imagecolorallocate($im,192,176,218);
$white = imagecolorallocate($im,240,240,240);
$azure = imagecolorallocate($im,32,96,255);
$sapphire = imagecolorallocate($im,32,16,224);

for ($fam=1; $fam<=16; $fam++)
{
    $i = $fam-1;
    $var = "or$fam";

    $i1 = boolval($i & 0x1);
    $i2 = boolval($i & 0x2);
    $i4 = boolval($i & 0x4);
    $i8 = boolval($i & 0x8);

    if (!($i & 0x8))
    {
        $r = ($i1 xor $i2 xor $i4) ? 255 : 0;
        $g =  $i1 ? 255 : 0;
        $b = ($i1 xor $i4) ? 255 : 0;
    }
    else
    {
        $r  = ((!$i1 xor $i4) && ($i != 0xd)) ? 0xc0 : 0;
        $r += 0x30;
        $r += !($i1 and $i2) ? 0x0f : 0;

        $g  = ($i1 and $i4) ? 0xc0 : 0;
        $g += ($i != 9) ? 0x30 : 0;
        $g += !($i1 and $i2) ? 0x0f : 0;

        $b  = ($i2 and $i4) ? 0xc0 : 0;
        $b += ($i2 or $i4) ? 0x30 : 0;
        $b += (!$i1 or $i1 == 0xd) ? 0x0c : 0;
        $b += !($i1 and $i2) ? 0x03 : 0;
    }

    $$var = imagecolorallocate($im,$r,$g,$b);
}

$or1    = imagecolorallocate($im,255,255,224);          // off-white, representing waxy, aldehydic, citronella.
$or2    = imagecolorallocate($im,250,192,  0);          // yellow representing sulfur, pineapple, Allium.
$or3    = imagecolorallocate($im,128,160, 96);          // earthy olive-green.
$or9    = imagecolorallocate($im, 64, 32,  0);          // dark brown representing leather, horse stable.
$or10   = imagecolorallocate($im,224,128, 32);          // orange-brown representing clove, resin, vanilla.
$or11   = imagecolorallocate($im,128,160,192);          // aqua-gray representing water, rock, mold, sweat.
$or51   = imagecolorallocate($im,255,255,192);
$or52   = imagecolorallocate($im,192,  0,224);
$or56   = imagecolorallocate($im,  0,128,  0);
$taar   = imagecolorallocate($im,  0,  0,160);
$vn1r   = imagecolorallocate($im,255,240,224);
$ms4a   = imagecolorallocate($im, 64,128,255);

$base = $h-$ybuf/2;
$bsht = 8;

imageline($im, 0,$base, $w,$base, $blue);

if (count($t))
{
    imagestring($im, 3, $w-29, 0, "Rel.", $red);
    imagestring($im, 3, $w-29,15, "Top" , $red);
}
else if (count($p))
{
    imagestring($im, 3, $w-29,0, "Dock", $azure);
    imagestring($im, 3, $w-36,15, "Score" , $azure);
}

if (count($t) || count($e))
{
    // Right labels first, so that left lines take precedence.
    for ($top = 1; $top <= floor($maxt); $top += 1)
    {   
        $dy = intval($base-1 - $tscale*$top);
        
        if (!($top & 1)) imageline($im, $xbuf/3,$dy, $w-$xbuf/3,$dy, $wine );
        imagestring($im, 3, $w-$xbuf/6,$dy-8, $top, $red);
    }

    // Left labels.
    imagestring($im, 3, 2,0, "log10", $green);
    imagestring($im, 3, 2,15, "EC50", $green);

    for ($ec = floor($maxe); $ec >= ceil($mine); $ec -= 1)
    {   
        $dy = intval($base-1 - $escale*($maxe-$ec));
        
        imageline($im, $xbuf/3,$dy, $w-$xbuf/3,$dy, $cyan);
        imagestring($im, 3, 2,$dy-8, $ec, $green);
    }
}

if (count($p))
{
    for ($score = floor($maxp); $score > 1; $score -= 7)
    {   
        $dy = intval($base-1 - $pscale*$score);

        imageline($im, $xbuf/3,$dy, $w-$xbuf/3,$dy, $sapphire);
        imagestring($im, 3, $w-$xbuf/6,$dy-8, $score, $azure);
    }
}

$bytree = [];
foreach ($prots as $rcpid => $pp)
{
    if (!isset($pp['btree'])) continue;
    $bytree[$pp['btree']] = $rcpid;
}

ksort($bytree, SORT_STRING);

foreach ($prots as $rcpid => $pp)
{
    if (!isset($pp['btree'])) $bytree[$rcpid] = $rcpid;
}

// foreach (array_keys($prots) as $x => $orid)
foreach (array_values($bytree) as $x => $orid)
{   
    $rcpcol = $white;
    switch (substr($orid,0,2))
    {   case 'OR':
        $bcol = 'or'.intval(preg_replace("/[^0-9]/","",substr($orid,2,2)));
        $bcol = $$bcol;
        break;
        
        case 'TA':
        $bcol = $taar;
        break;
        
        case 'VN':
        $bcol = $vn1r;
        break;
        
        case 'MS':
        $bcol = $ms4a;
        break;
        
        default:
        $bcol = $blue;
    }
    
    $dx  = $x * $res + $xbuf/2;
    $dyt = isset($t[$orid])
         ? intval($base-1 - $tscale*$t[$orid])
         : false;
    $dye = (isset($e[$orid]) && $e[$orid])
         ? intval($base-1 - $escale*($maxe-$e[$orid]))
         : false;
    $dyp = isset($p[$orid])
        ? intval($base-1 - $pscale*$p[$orid])
        : false;
    
    $base1 = $base2 = $base;
    if ($dyt >= $base) { $dyt += $bsht; $base1 += $bsht; }
    if ($dye >= $base) { $dye += $bsht; $base2 += $bsht; }

    $dy = $base;
    if (false!==$dyt && false===$dye)   imagefilledrectangle($im, $dx,$base1, $dx+$res-2,$dy=$dyt, $red);
    if (false!==$dye && false===$dyt)   imagefilledrectangle($im, $dx,$base2, $dx+$res-2,$dy=$dye, $green);
    if (false!==$dye && false!==$dyt)
    {   
        $dy = min($dye,$dyt);
        if ($dye < $dyt)
        {   if ($t[$orid] >= 0) imagefilledrectangle($im, $dx,$base2, $dx+$res-3,$dye, $green);
            imagefilledrectangle($im, $dx,$base1, $dx+$res-2,$dyt, $yellow);
        }
        else
        {   imagefilledrectangle($im, $dx+1,$base1, $dx+$res-2,$dyt, $red);
            imagefilledrectangle($im, $dx,$base2, $dx+$res-2,$dye, $yellow);
        }
    }

    if (false!==$dyp) imagefilledrectangle($im, $dx+2,$base2, $dx+$res-1,$dy=$dyp, $azure);
    
    if ($dy < $h/7) $dy = $h/7;
    if ($dy < $h/1.4) $texts[] = [$dx, $dy-5, $orid];
    
    imagefilledrectangle($im, $dx,$base, $dx+$res-1,$base+$bsht, $bcol);

}

if (!file_exists($fontfile = "assets/Montserrat.ttf"))
{
    $wofffile = str_replace(".ttf", ".woff2", $fontfile);
    $c = file_get_contents("https://fonts.gstatic.com/s/montserrat/v25/JTUHjIg1_i6t8kCHKm4532VJOt5-QNFgpCtr6Hw5aXo.woff2");
    $fp = fopen($wofffile, "wb");
    if (!$fp) die("Failed to open $wofffile for writing!");
    fwrite($fp, $c);
    fclose($fp);

    exec("woff2_decompress $wofffile");
}


$x = 53*$res;
$x1 = $x + 10*strlen($odor['full_name']);

foreach ($texts as $txt)
{   
    imagettftext($im, 9, 35, $txt[0], $txt[1], $pink, $fontfile, $txt[2]);

    if ($txt[0] >= $x && $txt[0] <= $x1) $x = $txt[0] + 50;
}

if ((!count($t) || !max($t)) && (!count($e) || !min($e)) && !count($p))
{   
    imagettftext($im, 28, 0, $w*0.42, $h*0.44, $pink, $fontfile, "(no data)");
}

// Odorant Name
imagestring($im, 5, $x, 2, $odor['full_name'], $blue);

header('Content-Type: image/png');
header('Content-Disposition: inline; filename="barchart.png"');
imagepng($im);
