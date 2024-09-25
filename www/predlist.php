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
        $ep = best_empirical_pair($rcpid, $oid, true);
        $newp['a'] = $ep;
        $newp['s'] = (isset($dr['Actual']) && @$dr['Actual'] != '(unknown)') ? (((@$dr['Predicted'] == 'Agonist') == (@$dr['Actual'] == 'Agonist')) ? "Y" : "N") : "?";
        $newp['d'] = $dr['CalculateDate'];
        $predictions[] = $newp;
    }
}

// print_r($predictions);

$fh = 44;
$fw = 30;
$nih = 25;
$nt = 20;
$noh = $nih + $nt;
$nw = 18;
$filter_svgdat = "m 0,0 $fw,$fh 0,$nih $nw,$nt 0,-$noh $fw,-$fh Z"

?>
<h1>Predictions List</h1>

<?php 
// echo "<svg height=\"81px\" viewBox=\"0 0 80 90\" xmlns=\"http://www.w3.org/2000/svg\"><path fill=\"magenta\" d=\"$filter_svgdat\"></path></svg>";
?>

<div class="box">
<div class="row content scrollh">

<?php if (isset($_REQUEST['r']) || isset($_REQUEST['o'])) { ?>
<a href="predlist.php">Clear filters</a>
<?php } ?>

<table class="liglist">
    <tr><th>Receptor</th>
        <th>Odorant</th>
        <th>Date</th>
        <th>Successful?</td>
    </tr>
<?php

$frcp = false;
$flig = false;
foreach ($predictions as $p)
{
    echo "<tr>\n";

    echo "<td><a href=\"receptor.php?r={$p['r']}\">{$p['r']}</a>";
    if ($frcp != $p['r'])
        echo " <a href=\"predlist.php?r={$p['r']}\"><svg height=\"13px\" viewBox=\"0 0 80 90\" xmlns=\"http://www.w3.org/2000/svg\"><path fill=\"#50cea8\" d=\"$filter_svgdat\"></path></svg></a>";
    echo "</td>\n";
    $frcp = $p['r'];

    $o = $odors[$p['o']];
    $fn = $o['full_name'];
    echo "<td><a href=\"odorant.php?o={$p['o']}\">$fn</a>";
    if ($flig != $p['o'])
        echo " <a href=\"predlist.php?o={$p['o']}\"><svg height=\"13px\" viewBox=\"0 0 80 90\" xmlns=\"http://www.w3.org/2000/svg\"><path fill=\"#50cea8\" d=\"$filter_svgdat\"></path></svg></a>";
    echo "</td>\n";

    echo @"<td>{$p['d']}</td>\n";

    switch ($p['s'])
    {
        case 'Y':
            $ps = "<span style=\"color: #0c0; font-weight: bold;\">Yes</span>";
            break;
        case 'N':
            $ps = "<span style=\"color: #f00; font-weight: bold;\">No</span>";
            break;
        default:
            $ps = "<span style=\"color: #9cf;\">???</span>";
    }
    echo @"<td>$ps</td>\n";
}


?></table>
