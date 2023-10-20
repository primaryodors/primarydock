<?php
require_once("../data/odorutils.php");

$oid = $_REQUEST['oid'];
$imgfname = "assets/pngs/$oid.png";
$odor = $odors[$oid];

if ( @$_REQUEST['refresh'] || !file_exists($imgfname))
{
    $smilesn = str_replace("[O-]", "O", $odor['smiles']);
    $smilesu = urlencode($smilesn);
    $url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/PNG";

    $ch = curl_init( $url );
    curl_setopt( $ch, CURLOPT_POST, 1);
    curl_setopt( $ch, CURLOPT_POSTFIELDS, "smiles=$smilesu");
    curl_setopt( $ch, CURLOPT_FOLLOWLOCATION, 1);
    curl_setopt( $ch, CURLOPT_HEADER, 0);
    curl_setopt( $ch, CURLOPT_RETURNTRANSFER, 1);

    $oimage = curl_exec( $ch );

    if (!$oimage)
    {
        $im = imagecreatetruecolor(300, 300);
        imagestring($im, 4, 37, 130, "Empty response from PubChem.", imagecolorallocate($im, 255,255,255));
        $imgfname = "assets/pngs/noimg.png";
        error_log("Empty response from $url smiles=$smilesu");
    }
    else
    {
        $im = imagecreatefromstring($oimage);
        if (!$im) die("Bad image data: \n\n$oimage\n\nSMILES: {$odor['smiles']}");
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
    }

    chdir(__DIR__);
    if (!file_exists("assets/pngs")) mkdir("assets/pngs");
    $fp = fopen($imgfname, "wb");
    if ($fp)
    {
        imagepng($im, $fp);
        fclose($fp);
    }
    else die("Failed to open $imgfname for writing!");
}

echo "<img src=\"$imgfname\">";