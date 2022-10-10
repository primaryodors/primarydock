<?php
chdir(__DIR__);
require_once("../predict/protutils.php");
require_once("../predict/odorutils.php");

$rcpid = @$_REQUEST['r'];
if (!$rcpid)
{
    header("Location: receptors.php");
    exit;
}

$receptor = @$prots[$rcpid];
if (!$receptor)
{
    header("Location: receptors.php");
    exit;
}

$page_title = $rcpid;

include("header.php");

echo "<h1>$rcpid</h1>\n";

?>
<div class="scrollh" style="height: 780px;">
<table class="liglist">
    <tr>
        <th>Odorant</th>
        <th>EC<sub>50</sub></th>
        <th>Adjusted Top</th>
        <th>Reference</th>
        <th>Aroma Notes</th>
    </tr>

<?php 

foreach ($odors as $oid => $odor)
{
    $pq = [];
    foreach ($odor['aroma'] as $refurl => $notes) $pq = array_merge($pq, $notes);
    $pq = array_unique($pq);

    // TODO: This approach is VERY inefficient.
    $pair = best_empirical_pair($rcpid, $oid, true);
    if ($pair)
    {
        $ufn = urlencode($odor['full_name']);
        echo "<tr>\n";
        echo "<td><a href=\"odorant.php?o=$ufn\">{$odor['full_name']}</a></td>\n";

        echo "<td>" . $dispec50 = (@$pair['ec50'] ?: "-") . "</td>\n";
        echo "<td>" . $disptop = (round(@$pair['adjusted_curve_top'], 4) ?: "-") . "</td>\n";
        echo "<td><a href=\"{$pair['ref']}\">&#x29c9;</a></td>\n";

        echo "<td>" . implode(", ",$pq) . "</td>\n";
        echo "</tr>\n";
    }
}