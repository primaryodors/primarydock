<?php
chdir(__DIR__);
require_once("../predict/protutils.php");
require_once("../predict/odorutils.php");

$odor = find_odorant(@$_REQUEST['o']);
if (!$odor)
{
    header("Location: odorants.php");
    exit;
}


$md5 = md5($odor['smiles']);
if (!file_exists($imgfname = "assets/pngs/$md5.png"))
{
    $oimage = file_get_contents($url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{$odor['smiles']}/PNG");
    if (!$oimage) die("Empty response from $url.");

    $im = imagecreatefromstring($oimage);
    $sx = imagesx($im);
    $sy = imagesy($im);
    if (!imageistruecolor($im)) imagepalettetotruecolor($im);

    for ($y=0; $y<$sy; $y++)
        for ($x=0; $x<$sx; $x++)
        {
            $c = imagecolorat($im, $x, $y);
            $arr['red']   = ($c & 0xff0000) >> 16;
            $arr['green'] = ($c &   0xff00) >> 8;
            $arr['blue']  = ($c &     0xff);

            $cmax = max($arr['red'], $arr['green'], $arr['blue']);
            $cmin = min($arr['red'], $arr['green'], $arr['blue']);
            $cavg = ($cmax + $cmin) / 2;
            $span = $cmax - $cmin;

            $newavg = 255 - $cavg;

            $red   = max(0, min(255, $newavg + ($arr['red']   - $cavg)));
            $green = max(0, min(255, $newavg + ($arr['green'] - $cavg)));
            $blue  = max(0, min(255, $newavg + ($arr['blue']  - $cavg)));

            $c = imagecolorallocate($im, $red, $green, $blue);
            imagesetpixel($im, $x, $y, $c);
        }

    if (!file_exists("assets/pngs")) mkdir("assets/pngs");
    $fp = fopen($imgfname, "wb");
    if ($fp)
    {
        imagepng($im, $fp);
        fclose($fp);
    }
    else die("Failed to open assets folder for writing!");
}

include("header.php");

?>
<h1><?php echo $odor['full_name']; ?></h1>

<img src="<?php echo $imgfname; ?>">