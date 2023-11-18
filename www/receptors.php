<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");

function filter_prot($protid, $filter)
{
    global $fnum, $prots;

    $disp = true;
    switch ($filter)
    {
        case 'a':                   // Has agonists.
        $ep = all_empirical_pairs_for_receptor($protid, true);
        foreach ($ep as $k => $v) if ($v <= 0) unset($ep[$k]);
        if (!count($ep)) $disp = false;
        else if ($fnum && count($ep) < $fnum) $disp = false;
        else if (max($ep) <= 0) $disp = false;
        break;

        case 'i':                   // Has inverse agonists and/or antagonists.
        $ep = all_empirical_pairs_for_receptor($protid, true);
        if (!count($ep)) $disp = false;
        else if (min($ep) >= 0) $disp = false;
        if (has_antagonists($protid)) $disp = true;
        break;

        case 'o':                   // Orphan.
        $ep = all_empirical_pairs_for_receptor($protid, true);
        foreach ($ep as $k => $v) if ($v == 0) unset($ep[$k]);
        if (count($ep)) $disp = false;
        break;

        case 's':
        $ftxt = substr($_REQUEST['f'], 1);
        if (substr($protid, 0, strlen($ftxt)) != $ftxt) $disp = false;
        break;

        case 'q':
        $disp = (@$prots[$protid]['hypothesized']);
        break;

        case '':
        break;

        default:
        $disp = false;
    }

    return $disp;
}


$filter = @$_REQUEST['f'] ?: "";
$fnum = intval(@substr($filter, 1));
$filter = substr($filter, 0, 1);

$extra_js = ['js/tabs.js'];
$extra_css = ['assets/tabs.css'];
include("header.php");

echo "<h1>Receptors</h1>\n";

?>

<div class="box">
<div class="row header">
<p>
<a href="receptors.php" class="<?php if ($filter == '') echo "hilite" ?>">All</a>
|
<a href="receptors.php?f=a55" class="<?php if ($filter == 'a' && $fnum == 55) echo "hilite" ?>">Big Five</a>
|
<a href="receptors.php?f=a" class="<?php if ($filter == 'a' && $fnum <= 1) echo "hilite" ?>">Having Agonists</a>
|
<a href="receptors.php?f=a2" class="<?php if ($filter == 'a' && $fnum == 2) echo "hilite" ?>">Having Multiple Agonists</a>
|
<a href="receptors.php?f=a10" class="<?php if ($filter == 'a' && $fnum == 10) echo "hilite" ?>">Having Many Agonists</a>
|
<a href="receptors.php?f=i" class="<?php if ($filter == 'i') echo "hilite" ?>">Having Antagonists</a>
|
<a href="receptors.php?f=o" class="<?php if ($filter == 'o') echo "hilite" ?>">Orphans</a>
</p>
</div>
<p style="position: absolute; right: 0px"><a href="#" onclick="$('#aboutpo').show();">About PO Numbering</a></p>

<div id="aboutpo" style="display: none;"><p>
The numbering system used on these pages is designed to address some problems with the prevailing OR numbering system.
The current dominant system groups olfactory receptors into families, denoted by a number 1-15 or 51-56;
then a subfamily denoted by one or two letters, with Z being followed by AA, AB, etc;
then a number representing the individual member gene of the subfamily, as well as the protein it encodes.
The families and subfamilies are assigned on the basis of sequence similarity.
However, phylogenetic analysis reveals that the family numbers are unreliable;
OR7 is for example a subset of OR1, while OR8, OR9, and most of OR5 form a single large branch of intermingled
subfamilies and individuals with no pattern to distinguish which parts should be numbered 5, 8, or 9.
Furthermore, several of the higher lettered OR5 subfamilies are actually part of OR14.</p>

<p>To address these issues, a new numbering scheme is proposed, one that strives to keep the existing subfamilies
and individuals as much as possible, but which more faithfully adheres to the phylogeny of the receptors.
</p></div>


<div class="row content scrollh">


<div class="tab" style="display: inline-block; margin-top: 30px;">
<button class="tablinks <?php if (!@$_REQUEST['tree']) echo "default"; ?>" id="tabList" onclick="openTab(this, 'List');">List</button>
<button class="tablinks <?php if ( @$_REQUEST['tree']) echo "default"; ?>" id="tabTree" onclick="openTab(this, 'Tree');">Tree</button>
</div>



<div id="List" class="tabcontent">
<div class="fam">
<?php

$ffam = "OR1";
$echoed = 0;

$byrcpid = [];

foreach ($prots as $poid => $p) $byrcpid[$poid] = $p['id'];
natsort($byrcpid);

echo "<!-- ".print_r($byrcpid,1)." -->\n";

foreach ($byrcpid as $poid => $protid)
{
    $p = $prots[$poid];
    if (!@$p['region']) continue;
    $fam = family_from_protid($protid);
    if ($echoed && $fam != $ffam)
    {
        echo "</div><hr><div class=\"fam\">\n";
        $echoed = 0;
    }

    if (filter_prot($poid, $filter))
    {
        $disp = (substr($poid, 2) == substr($protid, 2)) ? $poid : "$poid ($protid)";
        echo "<a class=\"rcptile\" href=\"receptor.php?r=$poid\">$disp</a>\n";
        $echoed++;
    }
    $ffam = $fam;
}

?>
</div>
</div>


<div id="Tree" class="tabcontent">
<?php
$tree = [];
foreach ($prots as $protid => $p)
{
    if (!isset($p['btree'])) continue;
    if (!filter_prot($protid, $filter)) continue;
    $path = $p['btree'];
    if (!isset($tree[$path])) $tree[$path] = [];
    $tree[$path][] = $protid;
}

ksort($tree, SORT_STRING);

echo "<pre>";
$prev = [];
$path1 = "";
$lno = 0;
$ktree = array_keys($tree);
foreach ($tree as $path => $protids)
{
    $curr = str_split($path);
    // echo str_pad($path, 50);
    foreach ($curr as $i => $c)
    {
        if ($c == 0)
        {
            $sub = substr($path, 0, $i);
            $sub1 = $sub.'1';

            $found1 = false;
            for ($j=$lno+1; $j<count($tree); $j++)
            {
                $lookahead = $ktree[$j];
                if (substr($lookahead, 0, $i) != $sub)
                {
                    $c = '1';
                    break;
                }
                else if (substr($lookahead, 0, $i+1) == $sub1)
                {
                    $found1 = true;
                    break;
                }
            }
            if (!$found1) $c = '1';
        }
        if (isset($prev[$i]))
        {
            if (substr($path1, 0, $i+1) != substr($path, 0, $i+1))
            {
                if ($c === '0') echo "&#x252c;&#x2500;&#x2500;&#x2500;";       // +
                else if (substr($path1, 0, $i) == substr($path, 0, $i)) echo "&#x2514;&#x2500;&#x2500;&#x2500;";                  // `
                else echo "&#x2500;&#x2500;&#x2500;&#x2500;";                      // -
            }
            else
            {
                if ($c === '0') echo "&#x2502;&nbsp;&nbsp;&nbsp;";       // |
                else echo "&nbsp;&nbsp;&nbsp;&nbsp;";                         //
            }
        }
        else
        {
            if ($c == 0) echo "&#x252c;&#x2500;&#x2500;&#x2500;";              // +
            else echo "&#x2500;&#x2500;&#x2500;&#x2500;";                      // -
        }
    }

    echo "&#x2500;&#x2500;&#x2500;";
    foreach ($protids as $r)
    {
        echo " <a href=\"receptor.php?r=$r\">$r</a>";
        if ($filter == 'q') echo " ".$prots[$r]['hypothesized'];
    }
    echo "\n";

    $prev = str_split($path);
    $path1 = $path;
    $lno++;
}

?></div>
