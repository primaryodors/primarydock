<?php
chdir(__DIR__);
require_once("../predict/protutils.php");
require_once("../predict/odorutils.php");

include("header.php");

echo "<h1>Receptors</h1>\n";

?><div class="scrollw"><div class="fam"><?php

$ffam = "OR1";
foreach ($prots as $protid => $p)
{
    $fam = family_from_protid($protid);
    if ($fam != $ffam) echo "</div><div class=\"fam\">\n";
    echo "<a class=\"rcptile\" href=\"receptor.php?r=$protid\">$protid</a>\n";
    $ffam = $fam;
}