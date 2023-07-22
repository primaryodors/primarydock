<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");

$filter = @$_REQUEST['f'] ?: "";
$fnum = intval(@substr($filter, 1));
$filter = substr($filter, 0, 1);

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
<div class="row content scrollh"><div class="fam">
<?php

$ffam = "OR1";
$echoed = 0;
foreach ($prots as $protid => $p)
{
    $fam = family_from_protid($protid);
    if ($echoed && $fam != $ffam)
    {
        echo "</div><hr><div class=\"fam\">\n";
        $echoed = 0;
    }

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

        case '':
        break;

        default:
        $disp = false;
    }

    if ($disp)
    {
        echo "<a class=\"rcptile\" href=\"receptor.php?r=$protid\">$protid</a>\n";
        $echoed++;
    }
    $ffam = $fam;
}
