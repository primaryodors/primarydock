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
$imgfname = "assets/pngs/$md5.png";
if ( @$_REQUEST['refresh'] || !file_exists($imgfname))
{
    $smilesu = urlencode($odor['smiles']);
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
        error_log("Empty response from $url");
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

$page_title = $odor['full_name'];
$extra_js = ['js/tabs.js'];
$extra_css = ['assets/tabs.css'];

include("header.php");

?>
<script>
var viewer_loaded = false;
function load_viewer(obj)
{
    openTab(obj, 'Structure');
    if (!viewer_loaded)
    {
        window.setTimeout( function()
        {
            $('#viewer').on('load', function()
            {
                var embdd = $('#viewer')[0];
                $("[type=file]", embdd.contentDocument).hide();
                var filediv = $("#filediv", embdd.contentDocument)[0];
                
                filediv.innerText = "<?php echo @$odrow['full_name']; ?>";
            });
            $('#viewer')[0].src = '<?php echo "viewer.php?url=sdf.php&mol={$odor['oid']}"; ?>'; 
        }, 259); 
        viewer_loaded = true;
    }
}

window.setTimeout( function()
{
    var boundary = parseInt($(".tab")[0].getClientRects()[0].bottom);
    $("#tabAroma").click();
}, 123);

</script>
<div class="tab" style="display: inline-block; margin-top: 30px;">
    <button class="tabstatic" id="tabFullName"><?php echo $odor['full_name']; ?></button>
	<button class="tablinks" id="tabAroma" onclick="openTab(this, 'Aroma');">Notes & Receptors</button>
	<button	class="tablinks"
			id="tabStructure"
			onclick="load_viewer(this);"
			>3D Structure</button>
</div>

<div id="Aroma" class="tabcontent">

<div class="scrollw">
    <div>
        <img class="skeletal" src="<?php echo $imgfname; ?>">
        <img class="barchart" src="barchart.php?o=<?php echo urlencode($odor['smiles']); ?>">
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
$agonist = [];
foreach ($odor['activity'] as $refurl => $acv)
{
    $maxcurvtop = [];
    $minec50 = [];
    foreach ($acv as $rcpid => $a)
    {
        $maxcurvtop[$rcpid] = false;
        $minec50[$rcpid] = false;

        if (!isset($sorted[$rcpid])) $sorted[$rcpid] = 0.0;
        $ssamples = 0;
        if (isset($a['adjusted_curve_top']))
        {
            if (!isset($tbltops[$rcpid])) $tbltops[$rcpid] = "";
            else $tbltops[$rcpid] .= ", ";

            $tbltops[$rcpid] .= round($a['adjusted_curve_top'], 4) . " <sup><a href=\"$refurl\">$refno</a></sup>";
            $sorted[$rcpid] += $a['adjusted_curve_top'];
            $ssamples++;

            if (false===$maxcurvtop[$rcpid] || $a['adjusted_curve_top'] > $maxcurvtop[$rcpid]) $maxcurvtop[$rcpid] = $a['adjusted_curve_top'];
        }
        if (isset($a['ec50']))
        {
            if (!isset($tblec50[$rcpid])) $tblec50[$rcpid] = "";
            else $tblec50[$rcpid] .= ", ";

            $tblec50[$rcpid] .= round($a['ec50'], 4) . " <sup><a href=\"$refurl\">$refno</a></sup>";
            $sorted[$rcpid] -= $a['ec50']*1.666;
            $ssamples++;

            if (false===$minec50[$rcpid] || $a['ec50'] < $minec50[$rcpid]) $minec50[$rcpid] = $a['ec50'];
        }
        if ($ssamples) $sorted[$rcpid] /= $ssamples;
    }
    $refno++;

    foreach (array_keys($maxcurvtop) as $rcpid)
    {
        if ($maxcurvtop[$rcpid] > 0) $agonist[$rcpid] = true;
        if ($maxcurvtop[$rcpid] >= 0 && $minec50[$rcpid] < 0) $agonist[$rcpid] = true;
    }
}

arsort($sorted);


foreach (array_keys($sorted) as $rcpid)
{
    echo "<tr>\n";
    echo "<td><a href=\"receptor.php?r=$rcpid\">$rcpid</a></td>\n";

    echo "<td style=\"white-space: nowrap;\">" . $dispec50 = (@$tblec50[$rcpid] ?: "-") . "</td>\n";
    echo "<td style=\"white-space: nowrap;\">" . $disptop = (@$tbltops[$rcpid] ?: "-") . "</td>\n";

    if (@$agonist[$rcpid])
        echo "<td style=\"white-space: nowrap;\">" . substr(get_notes_for_receptor($rcpid, $correlations), 0, 123) . "</td>\n";
    else
        echo "<td>&nbsp;</td>\n";

    echo "</tr>\n";
}

?>
</table>
</div>

</div>

<div id="Structure" class="tabcontent">
    <iframe id="viewer"></iframe>
</div>
