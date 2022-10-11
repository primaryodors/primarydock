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

$pairs = all_empirical_pairs_for_receptor($rcpid);
// die("<pre>".print_r($pairs,true));

$refs = [];
foreach ($pairs as $oid => $pair)
{
    $odor = $odors[$oid];

    if (@$pair['ec50_ref'])
    {
        $refno_ec50 = array_search($pair['ec50_ref'], $refs);
        if (false === $refno_ec50)
        {
            $refno_ec50 = count($refs)+1;
            $refs[$refno_ec50] = $pair['ec50_ref'];
        }
    }
    if (@$pair['top_ref'])
    {
        $refno_top = array_search($pair['top_ref'], $refs);
        if (false === $refno_top)
        {
            $refno_top = count($refs)+1;
            $refs[$refno_top] = $pair['top_ref'];
        }
    }

    $pq = [];
    foreach ($odor['aroma'] as $refurl => $notes) $pq = array_merge($pq, $notes);
    $pq = array_unique($pq);

    $ufn = urlencode($odor['full_name']);
    echo "<tr>\n";
    echo "<td><a href=\"odorant.php?o=$oid\">{$odor['full_name']}</a></td>\n";

    echo "<td>" . $dispec50 = (@$pair['ec50'] ? ("{$pair['ec50']} <sup><a href=\"{$pair['ec50_ref']}\">$refno_ec50</a></sup>") : "-") . "</td>\n";
    echo "<td>" . $disptop = (@$pair['adjusted_curve_top'] ? (round(@$pair['adjusted_curve_top'], 4) . " <sup><a href=\"{$pair['top_ref']}\">$refno_top</a>") : "-") . "</sup></td>\n";

    echo "<td>" . implode(", ",$pq) . "</td>\n";
    echo "</tr>\n";
}