<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");
require_once("../predict/statistics.php");

$note = @$_REQUEST['n'];
if (!$note)
{
    header("Location: odorants.php");
    exit;
}

if (!file_exists($pqcorr = "../data/receptor_pq.json"))
{
    $correlations = correlate_receptors_aromanotes();

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

$orlist = [];
foreach ($correlations as $orid => $ornotes)
{
    if (!isset($ornotes[$note])) continue;

    $orlist[$orid] = floatval($ornotes[$note]);
}

$max = max($orlist);
arsort($orlist);

$page_title = $note;
include("header.php");

echo "<h1>Receptors that correlate to $note</h1>";
if (!count($orlist)) die("Insufficient data.");
echo "<table width=\"50%\">\n";
foreach ($orlist as $orid => $strength)
{
    $poid = find_poid($orid);
    if ($strength < 0.25 * $max) break;
    $width = $strength * 100.0 / $max;
    echo "<tr><th width=\"5%\"><a href=\"receptor.php?r=$poid\">$poid</a></th>
        <td><div style=\"display: inline-block; background-color: #ffc; width: $width%;\">&nbsp;</div></td></tr>\n";
}
