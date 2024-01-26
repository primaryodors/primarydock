<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");

include("header.php");

chdir(__DIR__);
$dock_results = json_decode(file_get_contents("../predict/dock_results.json"), true);
$predictions = [];

foreach ($dock_results as $rcpid => $results)
{
    if (isset($_REQUEST['r']))
    {
        if ($rcpid != $_REQUEST['r']
            &&
            (   substr($_REQUEST['r'], -1) != '*'
                ||
                substr($rcpid, 0, strlen($_REQUEST['r'])-1) != substr($_REQUEST['r'], 0, -1)
            )
            &&
            !preg_match("/^{$_REQUEST['r']}$/", $rcpid)
        )
        continue;
    }

    foreach ($results as $lig => $dr)
    {
        $odor = find_odorant($lig);
        $oid = $odor['oid'];

        if (isset($_REQUEST['o']) && $_REQUEST['o'] != $oid) continue;

        $newp['r'] = $rcpid;
        $newp['o'] = $oid;
        $newp['p'] = $dr['DockScore'];
        $ep = best_empirical_pair($rcpid, $oid);
        $newp['a'] = $ep;
        $predictions[] = $newp;
    }
}

?>
<h1>Predictions List</h1>

<div class="box">
<div class="row content scrollh">

<?php if (isset($_REQUEST['r']) || isset($_REQUEST['o'])) { ?>
<a href="predlist.php">Clear filters</a>
<?php } ?>

<table class="liglist">
    <tr><th>Receptor</th>
        <th>Odorant</th>
        <th>Empirical</th>
        <th>Prediction</th>
    </tr>
<?php

foreach ($predictions as $p)
{
    echo "<tr>\n";

    echo "<td><a href=\"predlist.php?r={$p['r']}\">{$p['r']}</a></td>\n";

    $o = $odors[$p['o']];
    $fn = $o['full_name'];
    echo "<td><a href=\"predlist.php?o={$p['o']}\">$fn</a></td>\n";

    $type = $typenames[intval($p['a'])];
    echo "<td>$type</td>\n";
    echo "<td>{$p['p']}</td>\n";
}


?></table>
