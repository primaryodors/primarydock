<?php

// metal.php
//
// Creates a zinc-binding EXR2 helix for every OR that has a compatible HXXCD/E motif and saves the metalloprotein PDB.
//

require("protutils.php");

foreach ($prots as $k => $v)
{
    $protid = $k;
    $fam = family_from_protid($protid);
    if (!file_exists("pdbs/$fam/$protid.upright.pdb")) continue;
    if (file_exists("pdbs/$fam/$protid.metal.pdb")) continue;
    
    echo "Processing $protid...\n";
    chdir(__DIR__);
    chdir("..");
    set_time_limit(600);
    passthru("bin/peptiditor predict/exrhelix.pepd $protid");
    echo "\n";
}
