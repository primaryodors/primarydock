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
<p style="position: absolute; right: 0px"><a href="#" onclick="$('#aboutpo').show(); document.body.style.overflowY = 'auto';">About Odr Numbering</a></p>

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
The new numbers begin with Odr, short for &quot;odor receptor&quot;, followed by a trunk number, a branch letter,
and a member number, so that the format stays very similar to the old system.
By design, many receptors have the same number in both systems, such as Odr1G1/OR1G1; others necessarily required a
different number, for example Odr1G2/OR1I1.
</p>

<p>While the new Odr scheme strives to keep all its trunks and branches monophyletic, Odr5 is not.
Therefore, Odr5/8/9 suffers from the same limitation as OR1/7 does.
However, some classification groups are paraphyletic, such as the terms "fish", "crustacean", and "monkey",
and so is the human Odr5 branch.
</p>

<p>Here's an overview of the correspondences that will allow recognition of the majority of Odr numbers:</p>

<ul>
    <li>The following receptors and subfamilies are the same in both systems:<br>
    1A, 1B1, 1D1-5, 1E1-3, 1F1-2, 1G1, 1J, 1L1-8, 1M, 1S, 3A, 4A, 4B1, 4C11, 4C16, 4D1-4, 4D6-11, 4E, 4F, 4K, 4L, 4N1-9, 4Q1-3,
    4S1-2, 5A1-2, 5H1-15, 5F1, 5K, 5M1, 5M10, 5P1-3, 8B1-12, 8D1-4, 8J, 8K, 8U1-9, 10A1-7, 10C1, 10H, 10J1-5, 10K, 10Q1, 10R2, 11, 12,
    14A2, 14I1, 14J1, 51A, 51B1-6, 51E1-2, 51F, 51G1-2, 51H, 51I, 51L, 51M, 51Q, 51S, 51T, 52A1-4, 52B1-6, 52E1-8, 52K, 52L, 52P, 52W, 56.
    <br>&nbsp;
    </li>

    <li>The following are direct equivalents:<br>
    Odr1A = OR7, Odr21G = OR2AG, Odr22D = OR2D, Odr22F = OR2F, Odr22H = OR13H, Odr22J = OR13J, Odr23B1-6 = OR2B1-6, Odr4N11-12 = OR4M1-2,
    Odr23J = OR2J, Odr23W = OR2W, Odr4AC = OR4C (not 11 or 16), Odr5H32 = OR5AC2, Odr5P6-7 = OR5AK2-3, Odr5GK1,2,4,9 = OR9G1,K2,G4,G9,
    Odr16K,N = OR6K,N, Odr17F = OR6F, Odr17J1 = OR6J1, Odr17S = OR6S, Odr17T = OR6T, Odr17V1 = OR6V1, Odr17V2-4 = OR9A2-4, Odr17X = OR6X,
    Odr18P1 = OR6P1, Odr18Q = OR6Q, Odr8B21-25 = OR8G1-5, Odr8M = OR5M, Odr9B = OR5B, Odr9C1 = OR5C1, Odr9D = OR5D, Odr9L = OR5L,
    Odr13A1 = OR10AD1, Odr19D = OR10D, Odr52A5 = OR52N5, Odr52A6 = OR52N4, Odr58A = OR52I.
    <br>&nbsp;
    </li>

    <li>The following are mostly direct equivalents:<br>
    Odr21A = 2A (Odr21A34 = OR2AT4), Odr21L = OR2L (Odr21L21 = OR2AJ1), Odr22C = OR13C (Odr22C10 = OR2S2), Odr23G = OR2G (Odr23G11 = OR2B11),
    Odr23H = OR2H (Odr23H3 = OR2C1), Odr5Q = OR9Q (Odr5Q3 = OR9I1), Odr17C = OR6C (Odr17C80 = OR2AP1), Odr18B = OR6B (Odr18B12 = OR6A2),
    Odr9T = OR5T (Odr9T4 = OR5AP2), Odr19G = OR10G (Odr19G19 = OR10S1).
    <br>&nbsp;
    </li>

    <li>OR2M/T/V/Z have been reassigned into Odr21M, 21T, and 21V according to their phylogeny.</li>
</ul>

    <a href="#" onclick="$('#aboutpo').hide();">(close)</a>
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
uksort($byrcpid, "rcpid_cmp");

// echo "<!-- ".print_r($byrcpid,1)." -->\n";

foreach ($byrcpid as $poid => $protid)
{
    $p = $prots[$poid];
    if (!@$p['region']) continue;
    $fam = family_from_protid($poid);
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
