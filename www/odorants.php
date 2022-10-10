<?php
chdir(__DIR__);
require_once("../predict/protutils.php");
require_once("../predict/odorutils.php");

include("header.php");

?>
<h1>Odorants List</h1>

<div class="scrollh" style="height: 780px;">
<?php

if (!file_exists($pqcorr = "../data/receptor_pq.json"))
{
    echo "Note: The first time an odorant's page is loaded, the app will one-time run receptor-aroma-note correlations.";
    echo " This calculation takes a few moments.<br><br>";
}

foreach ($odors as $oid => $odor)
{
    $fnu = urlencode($odor['full_name']);
    echo "<a href=\"odorant.php?o=$oid\">{$odor['full_name']}</a><br>\n";
}