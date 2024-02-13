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

echo "<p><a name=\"top\"></a>";
for ($l = 1; $l < 28; $l++)
{
    if ($l <= 26) echo "<a href=\"#letter_$l\">".chr(64+$l)."</a> ";
    else echo "<a href=\"#letter_123\">#</a>";
}
echo "</p>";

echo "<div class=\"odorlist\">\n";
for ($l = 1; $l < 28; $l++)
{
    if ($l <= 26) echo "<a name=\"letter_$l\"><h2>".chr(64+$l)."</h2></a>";
    else echo "<a name=\"letter_123\"><h2>Others</h2></a>";

    $tl = [];
    foreach (array_keys($odors) as $oid)
    {
        $odor = $odors[$oid];
        $C = strtoupper(substr(trim_prefixes($odor['full_name']), 0, 1));
        if ($l <= 26 && $C != chr(64+$l)) continue;
        if ($l > 26 && ($C >= 'A' && $C <= 'Z')) continue;

        if (@$_REQUEST['gsac'])
        {
            $notes = @$odor['aroma']["http://www.thegoodscentscompany.com"];
            if (is_array($notes) && count($notes))
            {
                if (count($notes) < 2) continue;
                natsort($notes);
                if (implode(" ", $notes) != implode(" ", $odor['aroma']["http://www.thegoodscentscompany.com"])) continue;
            }
            else continue;
        }

        $tl[$oid] = hellenicize($odor['full_name']);
    }

    $cols = 5;
    $pcnt = intval(100/$cols);
    $pcol = ceil(count($tl)/$cols);

    echo "<div style=\"display: flex; width: 100%;\">\n";
    echo "<div style=\"width: $pcnt%;\">\n";
    $o = 0;
    $col = 1;
    foreach ($tl as $oid => $name)
    {
        $fnu = urlencode($odor['full_name']);
        echo "<a href=\"odorant.php?o=$oid\">$name</a> ";
        echo "<br>\n";
        $o++;
        if ($o >= $pcol && $col < $cols)
        {
            echo "</div><div style=\"width: $pcnt%;\">\n";
            $o = 0;
            $col++;
        }
    }
    echo "</div>\n</div>\n";
    echo "<p><a href=\"#top\">(back to top)</a></p>";
    echo "<hr>";
}
echo "</div>\n";
