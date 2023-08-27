<?php
chdir(__DIR__);

$c = file_get_contents("../viewer.htm");

if (@$_REQUEST['view'] == "pred")
{
    require("../data/protutils.php");
    $protid = $_REQUEST["prot"];
    $fam = family_from_protid($protid);
    $odor = $_REQUEST["odor"];
    $mode = $_REQUEST["mode"];      // active or inactive.
    $n = @$_REQUEST["n"] ?: 1;

    $path = "../output/$fam/$protid/$protid.$odor.$mode.model$n.pdb";
    if (!file_exists($path)) $path = "../output/$fam/$protid/$protid.$odor.model$n.pdb";
    if (!file_exists($path)) die("Something went wrong. $path");
    $pdb = file_get_contents($path);

    $c = str_replace("var literal_pdb = false;\n", "var literal_pdb = `$pdb`;\n", $c);
    $c = str_replace("var literal_fname = \"\";\n", "var literal_fname = \"$protid.$odor.$mode.model$n.pdb\";\n", $c);
}

echo $c;
