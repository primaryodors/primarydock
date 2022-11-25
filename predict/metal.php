<?php

// metal.php
//
// Creates a zinc-binding EXR2 helix for every OR that has a compatible HXXCD/E motif and saves the metalloprotein PDB.
//
// If a change is made to this file or to exrhelix.pepd, this script can be rerun with the force option:
// php -f predict/metal.php force
//

require("protutils.php");

foreach (@$argv as $a)
{
	$a = explode('=',$a,2);
	$_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}

foreach ($prots as $k => $v)
{
    $protid = $k;
    $fam = family_from_protid($protid);
    if (!file_exists("pdbs/$fam/$protid.upright.pdb")) continue;
    if (!@$_REQUEST['force'] && file_exists("pdbs/$fam/$protid.metal.pdb")) continue;
    
    echo "Processing $protid...\n";
    chdir(__DIR__);
    chdir("..");
    set_time_limit(600);
    // $script = "exrhelix.pepd";
    $script = "nonexrmtl.pepd";
    passthru("bin/pepteditor predict/$script $protid");
    echo "\n";
}
