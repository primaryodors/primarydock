<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");

include("header.php");

?>
<h1>Odorants List</h1>

<div class="box">
<div class="row content scrollh">
<?php

if (!file_exists($pqcorr = "../data/receptor_pq.json"))
{
    echo "Note: The first time an odorant's page is loaded, the app will one-time run receptor-aroma-note correlations.";
    echo " This calculation takes a few moments.<br><br>";
}

echo "<div class=\"odorlist\">\n";
for ($l = 1; $l < 28; $l++)
{
    if ($l <= 26) echo "<h2>".chr(64+$l)."</h2>";
    else echo "<h2>Others</h2>";
    foreach (array_keys($odors) as $oid)
    {
        $odor = $odors[$oid];
        $C = strtoupper(substr(trim_prefixes($odor['full_name']), 0, 1));
        if ($l <= 26 && $C != chr(64+$l)) continue;
        if ($l > 26 && ($C >= 'A' && $C <= 'Z')) continue;
        $fnu = urlencode($odor['full_name']);
        echo "<a href=\"odorant.php?o=$oid\">".hellenicize($odor['full_name'])."</a> ";
        // echo trim_prefixes($odor['full_name']);
        echo "<br>\n";
    }
    echo "<hr>";
}
echo "</div>\n";
