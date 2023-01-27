<?php

// method_soft.php
//
// Performs a soft-binding dock of an odorant in a receptor.
//
// Example call syntax:
// php -f predict/method_soft.php prot=OR1A1 lig=geraniol
//

// Configurable variables
$dock_metals = false;
$bias_by_energy = true;

require("methods_common.php");

$cenres = "CEN RES $bsr2a $bsr3a $bsr3b $bsr3c $bsr3d $bsr3e $bsr3f $bsr3g $bsr4a $bsr4b $bsr4c $bsr5a $bsr5b $bsr5c $bsr5d $bsr6a $bsr6b $bsr7a $bsr7b $bsr7c";

$outfname = str_replace(".dock", "_soft.dock", $outfname);

prepare_outputs();

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres
SIZE 7.0 7.5 7.0

POSE 5
ITER 50
ELIM 99
SOFT TMR2 TMR4 TMR5 TMR6 TMR7

OUT $outfname
ECHO


heredoc;

function dockline_callback($inplines)
{
    $dosoft = false;
    $result = [];
    foreach ($inplines as $ln)
    {
        if ($dosoft)
        {
            $morceaux = explode(": ", $ln);
            if (count($morceaux) == 2) $result[$morceaux[0]] = floatval($morceaux[1]);
        }

        if (trim($ln) == 'Soft transformations:') $dosoft = true;
        if (trim($ln) == '') $dosoft = false;
    }

    return $result;
}

process_dock();

