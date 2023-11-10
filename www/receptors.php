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

    if (filter_prot($protid, $filter))
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
