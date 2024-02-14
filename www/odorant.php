<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");
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

$predictions = [];
chdir(__DIR__);
$dock_results = json_decode(file_get_contents("../predict/dock_results.json"), true);
$odorname_under = str_replace(' ', '_', $odor["full_name"]);
foreach ($dock_results as $protid => $dr)
{
    if (isset($dr[$odorname_under]))
    {
        $predictions[$protid] = $dr[$odorname_under]["DockScore"];
    }
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
                $("#posey", embdd.contentDocument).css("position", "absolute").css("left", "-262144px");

                <?php if (@$odor['isomers'])
                {
                    echo "filediv.innerHTML = \"".@$odor['full_name'];
                    foreach (array_keys($odor['isomers']) as $iso)
                    {
                        $isou = urlencode($iso);
                        echo " &#xb7; <small><a href=\\\"#\\\" onclick=\\\"load_remote_sdf('sdf.php?mol={$odor['oid']}&iso=$isou');\\\">$iso</a></small>";
                    }
                    echo "\";";
                }
                else
                {
                    ?>filediv.innerText = "<?php echo @$odor['full_name']; ?>";
                <?php } ?>
            });
            $('#viewer')[0].src = '<?php echo "viewer.php?url=sdf.php&mol={$odor['oid']}"; ?>'; 
        }, 259); 
        viewer_loaded = true;
    }
}

window.setTimeout( function()
{
    $('#skeletal')[0].innerHTML = svg_from_smiles('<?php echo str_replace("\\", "\\\\", $odor['smiles']); ?>', 300, 300);
    var boundary = parseInt($(".tab")[0].getClientRects()[0].bottom);
    $("#tabAroma").click();
}, 123);


</script>
<div class="tab" style="display: inline-block; margin-top: 30px;">
    <button class="tabstatic" id="tabFullName"><?php echo $odor['full_name']; ?></button>
	<button class="tablinks" id="tabAroma" onclick="openTab(this, 'Aroma');">Notes & Receptors</button>
    <?php if (count($predictions)) { ?>
	<button class="tablinks" id="tabPred" onclick="openTab(this, 'Predict');">Predictions</button>
    <?php } ?>
    <button class="tablinks" id="tabRefs" onclick="openTab(this, 'Refs');">References</button>
	<button	class="tablinks"
			id="tabStructure"
			onclick="load_viewer(this);"
			>3D Structure</button>
</div>

<div id="Aroma" class="tabcontent">

<div class="scrollw">
    <div>
        <div id="skeletal">&nbsp;</div>
        <img class="barchart" src="barchart.php?o=<?php echo urlencode($odor['smiles']); ?>">
    </div>
</div>

<p class="aromainfo">
    <?php
    for ($i=1; isset($odor["name$i"]); $i++)
    {
        if ($i > 1) echo ", ";
        echo $odor["name$i"];
    }
    if (isset($odor["iupac"]))
    {
        if ($i > 1) echo ", ";
        echo $odor["iupac"];
        if ($i == 1) echo "<br><br>";
    }
    if ($i > 1) echo "<br><br>";
    ?>
	<strong>SMILES:</strong><br>
	<?php echo $odor['smiles']; ?><br>
	<br>
    <strong>Aroma Description:</strong>
    <br>
    <?php 
    $refno = 1;
    $lrefs = [];
    foreach ($odor['aroma'] as $refurl => $notes)
    {
        $comma = false;
        foreach ($notes as $note)
        {
            if ($comma) echo ", ";
            echo "$note";
            $comma = true;
        }
        echo "<sup><a href=\"#\" onclick=\"openTab($('#tabRefs')[0], 'Refs');\">$refno</a></sup><br>";
        $lrefs[] = $refurl;
        $refno++;
    }

    if (@$odor['comment']) echo "<br>{$odor['comment']}<br>";
    ?>
</p>

<div class="box" style="height: auto!important;">
<div class="row content scrollh">
<table class="rcplist">

<tr>
    <th>Receptor</th>
    <th>log10 EC<sub>50</sub></th>
    <th>Adj. Top</th>
    <th>Antagonist?</th>
    <th>Correlated Perceptual Qualities</th>
</tr>

<?php

$sorted = [];
$tbltops = [];
$tblec50 = [];
$agonist = [];
$antagonist = [];
if (@$odor['activity']) foreach ($odor['activity'] as $refurl => $acv)
{
    $maxcurvtop = [];
    $minec50 = [];
    $lrefs[] = $refurl;
    foreach ($acv as $rcpid => $a)
    {
        $maxcurvtop[$rcpid] = false;
        $minec50[$rcpid] = false;

        if (@$a['antagonist']) $antagonist[$rcpid] = "Y";

        if (!isset($sorted[$rcpid])) $sorted[$rcpid] = 0.0;
        $ssamples = 0;
        if (!isset($a['adjusted_curve_top']) && @$a['type'] == 'na') $a['adjusted_curve_top'] = 0;
        if (isset($a['adjusted_curve_top']))
        {
            if (!isset($tbltops[$rcpid])) $tbltops[$rcpid] = "";
            else $tbltops[$rcpid] .= ", ";

            $tbltops[$rcpid] .= round($a['adjusted_curve_top'], 4) . " <sup><a href=\"#\" onclick=\"openTab($('#tabRefs')[0], 'Refs');\">$refno</a></sup>";
            $sorted[$rcpid] += $a['adjusted_curve_top'];
            $ssamples++;

            if (false===$maxcurvtop[$rcpid] || $a['adjusted_curve_top'] > $maxcurvtop[$rcpid]) $maxcurvtop[$rcpid] = $a['adjusted_curve_top'];
        }
        if (isset($a['ec50']))
        {
            if (!isset($tblec50[$rcpid])) $tblec50[$rcpid] = "";
            else $tblec50[$rcpid] .= ", ";

            $tblec50[$rcpid] .= round($a['ec50'], 4) . " <sup><a href=\"#\" onclick=\"openTab($('#tabRefs')[0], 'Refs');\">$refno</a></sup>";
            if (!isset($a['adjusted_curve_top']) || floatval($a['adjusted_curve_top']) > 0)
            {
              $sorted[$rcpid] -= $a['ec50']*2;
              $ssamples++;
            }

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
    echo "<td style=\"white-space: nowrap;\">" . $disptop =
    (
        (floatval(@$tbltops[$rcpid]) || !floatval(@$tblec50[$rcpid]))
        ? @$tbltops[$rcpid]
        : "-"
    ) . "</td>\n";

    if (@$antagonist[$rcpid]) echo "<td>Y</td>";
    else echo "<td>&nbsp;</td>";

    if (@$agonist[$rcpid])
    {
        $notes = substr(get_notes_for_receptor($rcpid, $correlations), 0, 123);
        $notes = make_clickable_notes(array_slice(explode(", ", $notes), 0, 10));
        $notes = implode(", ", $notes);
        if (substr($notes, 0, 1) == '(')
        {
            $notes = "<i class=\"dim\">$notes</i>";
        }
        echo "<td style=\"white-space: nowrap;\">$notes</td>\n";
    }
    else
    {
        echo "<td>&nbsp;</td>\n";
    }    

    echo "</tr>\n";
}

?>
</table>
</div>
</div>

</div>
<div id="Predict" class="tabcontent">

<div class="scrollw">
    <div>
        <img class="skeletal" src="<?php echo $imgfname; ?>">
        <img class="barchart" src="barchart.php?m=p&o=<?php echo urlencode($odor['smiles']); ?>&t=<?php echo time(); ?>">
    </div>
</div>

<p class="aromainfo">
    <?php
    for ($i=1; isset($odor["name$i"]); $i++)
    {
        if ($i > 1) echo ", ";
        echo $odor["name$i"];
    }
    if (isset($odor["iupac"]))
    {
        if ($i > 1) echo ", ";
        echo $odor["iupac"];
        if ($i == 1) echo "<br><br>";
    }
    if ($i > 1) echo "<br><br>";
    ?>
	<strong>SMILES:</strong><br>
	<?php echo $odor['smiles']; ?><br>
	<br>
    <strong>Aroma Description:</strong>
    <br>
    <?php 
    foreach ($odor['aroma'] as $refurl => $notes)
    {
        $comma = false;
        foreach ($notes as $note)
        {
            if ($comma) echo ", ";
            echo "$note";
            $comma = true;
        }
    }

    if (@$odor['comment']) echo "<br>{$odor['comment']}<br>";
    ?>
</p>

<div class="box" style="height: auto!important;">
<div class="row content scrollh">
<table class="rcplist">

<tr>
    <th>Receptor</th>
    <th>Docking Score</th>
    <th>Correlated Perceptual Qualities</th>
</tr>

<?php

arsort($predictions);

foreach ($predictions as $rcpid => $score)
{
    echo "<tr>\n";
    echo "<td><a href=\"receptor.php?r=$rcpid\">$rcpid</a></td>\n";

    echo "<td style=\"white-space: nowrap;\">$score</td>\n";

    if ($score > 0)
    {
        $notes = substr(get_notes_for_receptor($rcpid, $correlations), 0, 123);
        $notes = make_clickable_notes(explode(", ", $notes));
        $notes = implode(", ", $notes);
        if (substr($notes, 0, 1) == '(')
        {
            $notes = "<i class=\"dim\">$notes</i>";
        }
        echo "<td style=\"white-space: nowrap;\">$notes</td>\n";
    }
    else
    {
        echo "<td>&nbsp;</td>\n";
    }    

    echo "</tr>\n";
}

?>
</table>
</div>
</div>

</div>

<div id="Refs" class="tabcontent">
<div class="box">
<div class="row content scrollh">
<?php
foreach ($lrefs as $idx => $refurl)
{
    echo "<a href=\"$refurl\"><p>\n";
    $idx1 = $idx + 1;
    echo "$idx1.) ";
    echo $refs[$refurl]['citation'];
    echo "</p></a>\n";
}
?>
<!-- p>References for aroma perceptual qualities should not be taken to indicate that the authors of outside studies necessarily
assigned aroma notes to the neurons that receive input from any given receptor.
Rather, the findings of outside studies often constitute the information on which we base our own perceptual quality assignments.
</p -->
</div>
</div>
</div>

<div id="Structure" class="tabcontent">
    <iframe id="viewer"></iframe>
</div>
