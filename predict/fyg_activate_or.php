<?php

chdir(__DIR__);
require_once("../data/protutils.php");
chdir(__DIR__);
require_once("template.php");

chdir(__DIR__);
chdir("..");
foreach ($prots as $protid => $prot)
{
    $fam = family_from_protid($protid);
    $pdbfname_active = "pdbs/$fam/$protid.active.pdb";
    if (!file_exists($pdbfname_active) || filemtime($pdbfname_active) < filemtime("bin/fyg_activate_or"))
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
