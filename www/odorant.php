<?php
chdir(__DIR__);
require_once("../predict/protutils.php");
require_once("../predict/odorutils.php");
require_once("receptor_notes.php");

$odor = find_odorant(@$_REQUEST['o']);
if (!$odor)
{
    header("Location: odorants.php");
    exit;
}

if (!file_exists($pqcorr = "../data/receptor_pq.json"))
{
    $correlations = correlate_receptors_aromanotes();
    // die("<pre>".print_r($correlations, true));

    $f = fopen($pqcorr, "wb");
    if (!$f) die("FAILED to open $pqcorr for writing.");
    fwrite($f, preg_replace("/([ \t]*)([^\\s]*) ([\\[{])\n/", "\$1\$2\n\$1\$3\n", json_encode($correlations, JSON_PRETTY_PRINT)));
    fclose($f);
}
else
{
    chdir(__DIR__);
    $correlations = json_decode(file_get_contents($pqcorr), true);
}

$md5 = md5($odor['smiles']);
if (!file_exists($imgfname = "assets/pngs/$md5.png"))
{
    $smilesu = urlencode($odor['smiles']);
    $oimage = file_get_contents($url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/$smilesu/PNG");
    if (!$oimage)
    {
        $im = imagecreatetruecolor(300, 300);
        imagestring($im, 4, 37, 130, "Empty response from PubChem.", imagecolorallocate($im, 255,255,255));
        $imgfname = "assets/pngs/noimg.png";
    }
    else
    {
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

$page_title = $odor['full_name'];

include("header.php");

?>
<h1><?php echo $odor['full_name']; ?></h1>

<div class="scrollw">
    <div>
        <img class="skeletal" src="<?php echo $imgfname; ?>">
        <img class="barchart" src="barchart.php?o=<?php echo $odor['full_name']; ?>">
    </div>
</div>

<p class="aromainfo">
    <strong>Aroma Description:</strong>
    <br>
    <?php 
    $refno = 1;
    foreach ($odor['aroma'] as $refurl => $notes)
    {
        $comma = false;
        foreach ($notes as $note)
        {
            if ($comma) echo ", ";
            echo "$note";
            $comma = true;
        }
        echo "<sup><a href=\"$refurl\">$refno</a></sup><br>";
        $refno++;
    }
    ?>
</p>

<div class="scrollh">
<table class="rcplist">

<tr>
    <th>Receptor</th>
    <th>log10 ec<sub>50</sub></th>
    <th>Adjusted Top</th>
    <th>Associated Perceptual Qualities</th>
</tr>

<?php

$sorted = [];
$tbltops = [];
$tblec50 = [];
foreach ($odor['activity'] as $refurl => $acv)
{
    foreach ($acv as $rcpid => $a)
    {
        if (!isset($sorted[$rcpid])) $sorted[$rcpid] = 0.0;
        $ssamples = 0;
        if (isset($a['adjusted_curve_top']))
        {
            if (!isset($tbltops[$rcpid])) $tbltops[$rcpid] = "";
            else $tbltops[$rcpid] .= ", ";

            $tbltops[$rcpid] .= round($a['adjusted_curve_top'], 4) . "<sup><a href=\"$refurl\">$refno</a></sup>";
            $sorted[$rcpid] += $a['adjusted_curve_top'];
            $ssamples++;
        }
        if (isset($a['ec50']))
        {
            if (!isset($tblec50[$rcpid])) $tblec50[$rcpid] = "";
            else $tblec50[$rcpid] .= ", ";

            $tblec50[$rcpid] .= round($a['ec50'], 4) . "<sup><a href=\"$refurl\">$refno</a></sup>";
            $sorted[$rcpid] -= $a['ec50'];
            $ssamples++;
        }
        if ($ssamples) $sorted[$rcpid] /= $ssamples;
    }
    $refno++;
}

arsort($sorted);


foreach (array_keys($sorted) as $rcpid)
{
    echo "<tr>\n";
    echo "<td><a href=\"receptor.php?r=$rcpid\">$rcpid</a></td>\n";

    echo "<td>" . $dispec50 = (@$tblec50[$rcpid] ?: "-") . "</td>\n";
    echo "<td>" . $disptop = (@$tbltops[$rcpid] ?: "-") . "</td>\n";

    if (floatval($disptop) > 0 || (floatval($dispec50) && floatval($disptop) >= 0))
        echo "<td>" . substr(get_notes_for_receptor($rcpid, $correlations), 0, 123) . "</td>\n";
    else
        echo "<td>&nbsp;</td>\n";

    echo "</tr>\n";
}

?>
</table>
</div>
