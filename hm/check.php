<?php

chdir(__DIR__);
chdir("..");

require_once("data/protutils.php");
$dofix = @$argv[1] == "fix";

foreach ($prots as $rcpid => $p)
{
    $fam = family_from_protid($rcpid);
    if (file_exists("pdbs/$fam/$rcpid.upright.pdb") && !file_exists("pdbs/$fam/$rcpid.active.pdb"))
    {
        echo "$rcpid\n";
        if ($dofix) passthru("php -f hm/dohm.php $rcpid");
    }
}
