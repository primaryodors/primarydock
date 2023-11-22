<?php

chdir(__DIR__);
require_once("odorutils.php");
require_once("protutils.php");

function echo_usage()
{
    echo "Usage:\n\nphp -f gensdf.php molecule_name SMILES\n\n";
    exit;
}

$name = $argv[1];
if (!$name) echo_usage();

$smiles = $argv[2];
if (!$smiles) echo_usage();

$hash = md5($smiles);

$odors[$hash]["full_name"] = $name;
$odors[$hash]["smiles"] = $smiles;

ksort($odors);

$fp = fopen("odorant.json", "wb");
fwrite($fp, json_encode_pretty($odors));

chdir("..");
passthru("obabel -:\"$smiles\" --gen3D -osdf -O\"sdf/$name.sdf\"");
