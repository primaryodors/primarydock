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

<p>Here's an overview of the correspondences that will allow recognition of the majority of Odr numbers:
    <ul>
        <li>Odr1 single-letter generally equals OR1, with a few exceptions;</li>
        <li>Odr1AA through Odr1AG equals OR7 exactly;</li>
        <li>Odr21A generally equals OR2A and OR2AT;</li>
        <li>Odr21G corresponds to OR2AG;</li>
        <li>Odr21L/W/Y equal OR2L/W/Y;</li>
        <li>Odr21M, 21T, and 21V have been regrouped from original OR2M, 2T, 2V, and 2Z;</li>
        <li>Odr22C,H,J equal OR13C,H,J; Odr22K consists of OR13D and OR13F;</li>
        <li>Odr22D and Odr22F equal OR2D and OR2F;</li>
        <li>Odr23B 2,3,6 equals OR2B 2,3,6;</li>
        <li>Odr23G mostly corresponds to OR2G, however Odr23G11 is OR2B11;</li>
        <li>Odr23H/J/K mostly equals OR2H/J/K;</li>
        <li>Odr3A equals OR3A;</li>
        <li>Odr3B equals OR13A and OR13G;</li>
        <li>Odr4 mostly equals OR4, with a few exceptions;</li>
        <li>Some OR4C receptors have been assigned to Odr4AC;</li>
        <li>OR4M1,2 have been assigned to Odr4N11,12, as have OR4X1,2 to Odr4B11,12;</li>
        <li>Odr4S4 is OR4P4;</li>
        <li>Odr5A consists of OR5A and OR5AN;</li>
        <li>Odr5F,H,K correspond to OR5F,H,K, except Odr5F also includes OR8H and OR8I;</li>
        <li>OR5M has been split into Odr5M and Odr8M;</li>
        <li>Odr5P corresponds to OR5P, OR5J, and OR5AK;</li>
        <li>OR5U and most of the two-letter OR5 receptors have been assigned to Odr14;</li>
        <li>Odr5GK consists of OR9G and OR9K; Odr5Q corresponds to OR9Q;</li>
        <li>Odr9B,D,L equal OR5B,D,L;</li>
        <li>Odr9C1-5 consist of OR5C,W,I,AR,AS;</li>
        <li>Odr9T consists of OR5T and OR5AP;</li>
        <li>OR6 has been split up, with almost all members keeping the same letter and member number, into
            Odr18B, 17C, 17F, 17J, 16K, 16N, 18P, 18Q, 17S, 17T, 17V, and 17X;</li>
        <li>Odr18B12 is OR6A2; Odr17J2 is OR6M1; Odr18P2 is OR6Y1;</li>
        <li>Odr8B2-12 corresponds to OR8B2-12, and Odr8B21-25 to OR8G1-5;</li>
        <li>Odr8D equals OR8D, except Odr8D5 is OR8A1;</li>
        <li>Odr8J,K,U equal OR8J,K,U;</li>
        <li>Odr10A mostly equals OR10A; Odr10H equals OR10H; Odr10K equals OR10K;</li>
        <li>Odr10C consists of OR10C and OR10P; Odr10J of OR10J,X,Z; Odr10Q of OR10Q,V,W,AC; Odr10R of OR10R,T;</li>
        <li>Odr10A10 is OR8S1; Odr10AV1 is OR5V1;</li>
        <li>Odr11 and Odr12 equal OR11 and OR12 exactly;</li>
        <li>Odr17V is OR9A;</li>
        <li>Odr19D and Odr19G equal OR10D and OR10G, except Odr19G19 is OR10S1;</li>
        <li>14A2, 14I1, and 14J1 are the same in both systems;</li>
        <li>Odr51 is identical to OR51 except: Odr51E3 is OR51D1; Odr51G3 is OR51J1; Odr51B7 is OR51V1;</li>
        <li>Odr58A equals OR52I;</li>
        <li>Odr52E10 is OR52D1; Odr52B8 is OR52H1; Odr52A6 is OR52N4; Odr52A5 is OR52N5;</li>
        <li>All other Odr52 receptors equal OR52 exactly;</li>
        <li>Odr56 equals OR56 exactly.</li>
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
