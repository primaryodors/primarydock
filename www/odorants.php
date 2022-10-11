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

$sortable = [];
foreach ($odors as $oid => $odor)
{
    $fns = $fns1 = strtolower($odor['full_name']); 
    $fns = str_replace("alpha-",   "", $fns);
    $fns = str_replace("beta-",    "", $fns);
    $fns = str_replace("gamma-",   "", $fns);
    $fns = str_replace("delta-",   "", $fns);
    $fns = str_replace("epsilon-", "", $fns);
    $fns = preg_replace("/([(][A-Za-z0-9+-][)],?)+-/", "", $fns);
    $fns = preg_replace("/[(]([0-9]+[rs],?)+[)]-/", "", $fns);
    $fns = preg_replace("/[A-Za-z](,[A-Za-z])*-/", "", $fns);
    $fns = preg_replace("/^[a-zA-Z]-/", "", $fns);
    $fns = preg_replace("/[^a-zA-Z0-9][a-zA-Z]-/", "", $fns);
    $fns = preg_replace("/([0-9]+,?)+-/", "", $fns);
    $fns = preg_replace("/[^A-Za-z]/", "", $fns);   
    
    $fns1 = str_replace("alpha-",   "&#x3B1;-", $fns1);
    $fns1 = str_replace("beta-",    "&#x3B2;-", $fns1);
    $fns1 = str_replace("gamma-",   "&#x3B3;-", $fns1);
    $fns1 = str_replace("delta-",   "&#x3B4;-", $fns1);
    $fns1 = str_replace("epsilon-", "&#x3B5;-", $fns1);

    if (!$fns) $fns = $fns1;

    $sortable[$oid] = $fns . " ~ " . $fns1 . " ~ " . $odor['full_name'];
}

asort($sortable);

$limit = count($sortable) / 3;
$i = 0;

echo "<div class=\"odorlist\">\n";
foreach (array_keys($sortable) as $oid)
{
    $odor = $odors[$oid];
    $fnu = urlencode($odor['full_name']);
    // echo "<!-- {$sortable[$oid]} -->\n";
    echo "<a href=\"odorant.php?o=$oid\">{$odor['full_name']}</a><br>\n";

    $i++;
    if ($i > $limit)
    {
        echo "</div>\n";
        echo "<div class=\"odorlist\">\n";
        $i = 0;
    }
}
echo "</div>\n";
