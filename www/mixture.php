<?php
chdir(__DIR__);
require_once("../predict/protutils.php");
require_once("../predict/odorutils.php");
require_once("receptor_notes.php");

$mixttl = preg_replace("/[^A-Za-z0-9 +_-]/", "", @$_REQUEST['mixttl']) ?: "Mixture";

$mixture = [];
foreach ($_REQUEST as $key => $value)
{
    if (floatval($value))
    {
        $o = find_odorant($key);
        if ($o) $mixture[$o['oid']] = floatval($value);
    }
}

if (@$_REQUEST['add'])
{
    $o = find_odorant($_REQUEST['add']);
    if ($o) $mixture[$o['oid']] = floatval(@$_REQUEST['addwt'] ?: 1);
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


$page_title = $mixttl;

include("header.php");

// echo get_mol_wt("civetone");

?>
<form action="" method="GET">
<br>Mixture name:
<input type="text" name="mixttl" value="<?php echo $mixttl; ?>">
<table class="rcplist">
<tr>
    <th>Compound</th>
    <th>Rel. Amt. w/w</th>
    <th>Mol. wt.</th>
    <th>% molar</th>
<?php
$mixtbl = [];
$ttlmol = 0.0;
foreach ($mixture as $oid => $amtww)
{
    $odor = $odors[$oid];
    $r = [$odor['full_name']];

    $r[] = "<input size=\"4\" type=\"text\" name=\"$oid\" value=\"$amtww\" style=\"text-align: right;\"><br>\n";

    $r[] = $mw = get_mol_wt($odor['full_name']);
    $ttlmol += $r[] = $mw * $amtww;

    $mixtbl[$oid] = $r;
}

$moldiv = $ttlmol ? (100.0 / $ttlmol) : 1;
foreach ($mixtbl as $k => $r) $mixtbl[$k][3] = round( $r[3] * $moldiv, 2 );

foreach ($mixtbl as $r) echo "<tr><td>".implode("</td><td>", $r)."</td></tr>\n";
?>
<tr>
    <td>Add: <select name="add">
        <option value="">-----</option>
        <?php foreach ($odors as $oid => $o) echo "<option value=\"{$oid}\">{$o['full_name']}</option>\n"; ?>
        </select>
    </td>
    <td colspan="2">
        <input size="4" type="text" name="addwt" value="" style="text-align: right;">
    </td>
    <td>
        <input type="submit" value="Update">
    </td>
</table>
</form>

<div class="scrollh" style="heaight: 300px;">
<?php
$dbgec50 = 0;
$agby = [];
$antby = [];
$rcptbl = [];
foreach ($mixtbl as $oid => $r)
{
    $odor = $odors[$oid];
    $amtww = $mixture[$oid];
    $molarp = $r[3] / 100;

    $lacv = [];
    if (@$odor['activity']) foreach ($odor['activity'] as $ref => $acv)
    {
        // echo "$ref\n";
        foreach ($acv as $rcpid => $a)
        {
            if (!isset($lacv[$rcpid])) $lacv[$rcpid] = $a;
            else
            {
                if (@$a['adjusted_curve_top'] > $lacv[$rcpid]['adjusted_curve_top']) $lacv[$rcpid]['adjusted_curve_top'] = $a['adjusted_curve_top'];
                if (@$a['ec50'] && $a['ec50'] < $lacv[$rcpid]['ec50']) $lacv[$rcpid]['ec50'] = $a['ec50'];
                if (@$a['antagonist']) $lacv[$rcpid]['antagonist'] = 1;
            }
        }
    }

    foreach ($lacv as $rcpid => $a)
    {
        $p = @$rcptbl[$rcpid] ?: [$rcpid, 0, 0, "", substr(get_notes_for_receptor($rcpid, $correlations), 0, 123)];

        if (!isset($agby[$rcpid])) $agby[$rcpid] = [];
        if (!isset($antby[$rcpid])) $antby[$rcpid] = [];

        $top = $ec50 = 0;
        if (isset($a['adjusted_curve_top']) && isset($a['ec50']))
        {
            $top = floatval($a['adjusted_curve_top']);
            $ec50 = floatval($a['ec50']);
        }
        else if (isset($a['ec50']))
        {
            $f = floatval($a['ec50']);
            $ec50 = $f;
            $top = -($f+3.7)*5.3;
        }
        else if (isset($a['adjusted_curve_top']))
        {
            $top = floatval($a['adjusted_curve_top']);
            $ec50 = -abs($top/5.3 - 3.7);
        }

        if ($top > 0 || $ec50 < 0)
        {
            if ($p[1])
            {
                $f = pow(10, $ec50) / $molarp;
                $e = pow(10, $p[1]);
                if ($dbgec50) echo "$rcpid already {$p[1]} $e {$odor['full_name']} $ec50 $f ";
                $p[1] = min(log10( 1.0 / (1.0/$e + 1.0/$f)), 0);
                if ($dbgec50) echo "{$p[1]} <br>\n";
            }
            else
            {
                $f = pow(10, $ec50) / $molarp;
                if ($dbgec50) echo "$rcpid new {$odor['full_name']} $ec50 $molarp $f ";
                $p[1] = min(log10( pow(10, $ec50) / $molarp ), 0);
                if ($dbgec50) echo "{$p[1]} <br>\n";
            }

            if ($p[2])
            {
                $p[2] += max($top * $molarp, 0);
            }
            else $p[2] = max($top * $molarp, 0);
        }

        if ($top < 0 || @$a['antagonist'])
        {
            $p[3] = "Y";
            $antby[$rcpid][] = $odor['full_name'];
        }
        else if ($top > 0 || $ec50 < 0) $agby[$rcpid][] = $odor['full_name'];

        if ($top || $ec50 || @$a['antagonist']) $rcptbl[$rcpid] = $p;
    }
}

function rcpsort($a, $b)
{
    $a = $a[2]-1.666*$a[1];
    $b = $b[2]-1.666*$b[1];
    if ($a == $b) return 0;
    else return ($a < $b) ? 1 : -1;
}

uasort($rcptbl, 'rcpsort');

?>

<table class="rcplist">
<tr>
    <th>Receptor</th>
    <th>Est. EC<sub>50</sub></th>
    <th>Est. Top</th>
    <th>Antagonized?</th>
    <th>Correlated Notes</th>
</tr>
<?php
foreach ($rcptbl as $rcpid => $r)
{
    $r[1] = round($r[1], 3);
    $r[2] = round($r[2], 3);
    $bk = $r[3] ? "background-color: #2b2622;" : "";
    echo "<tr onclick=\"$('.hidable_$rcpid').show();\" style=\"$bk\"><td>".implode("</td><td>", $r)."</td></tr>\n";

    if (count($agby[$rcpid]))
    {
        $ag = implode(", ", array_unique($agby[$rcpid]));
        echo "<tr class=\"hidable_$rcpid\" style=\"background-color: #2b4e42; display: none;\" onclick=\"$(this).hide();\"><td colspan=\"5\">Agonized by: $ag</td></tr>\n";
    }

    if (count($antby[$rcpid]))
    {
        $ant = implode(", ", array_unique($antby[$rcpid]));
        echo "<tr class=\"hidable_$rcpid\" style=\"background-color: #2b1c26; display: none;\" onclick=\"$(this).hide();\"><td colspan=\"5\">Antagonized by: $ant</td></tr>\n";
    }
}

