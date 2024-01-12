<?php

chdir(__DIR__);
require_once("../data/protutils.php");
chdir(__DIR__);
require_once("template.php");

chdir(__DIR__);
chdir("..");
foreach ($prots as $protid => $prot)
{
    if (isset($argv[1]) && $protid != $argv[1] && !preg_match("/^{$argv[1]}$/", $protid)) continue;
    $fam = family_from_protid($protid);

    if (substr($fam, 0, 2) == "OR")
    {
        $fno = intval(substr($fam, 2));
        if ($fno == 51 || $fno == 52)
        {
            if (@$argv[1] == $protid) die("OR51 and OR52 receptors should use the direct model method.\n");
            else continue;
        }
    }

    $pdbfname_active = "pdbs/$fam/$protid.active.pdb";
    if (!file_exists($pdbfname_active) || filemtime($pdbfname_active) < filemtime("bin/fyg_activate_or") || $protid == @$argv[1])
    {
        if (check_already_fyg_activating($protid)) continue;
        build_template();
        if (check_already_fyg_activating($protid)) continue;
        else
        {
            do_templated_activation();
            if (file_exists($pdbfname_active) && filemtime($pdbfname_active) > filemtime("bin/fyg_activate_or")) exit;
        }
    }
}
