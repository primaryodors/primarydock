<?php

// method_basic.php
//
// Performs a dock of an odorant inside a default inactive-conformer PDB file.
//
// Example call syntax:
// php -f predict/method_basic.php prot=OR1A1 lig=R-carvone
//

// Configurable variables
$dock_metals = false;                   // This switches whether to use OR####.metal.pdb files when available, or to always use OR####.upright.pdb files.

chdir(__DIR__);
require("methods_common.php");
chdir(__DIR__);

// These residues are based on the binding site identified by Billesbolle et al:
// https://doi.org/10.1101/2022.12.20.520951
$cenres = "CEN RES 4.57 4.60 5.39 45.52";

// This will set some variables for the dock, including handling the "next" parameter.
prepare_outputs();

// Optional custom array of output metrics. If removed, the default will be as shown.
$metrics_to_process =
[
    "BENERG" => "BindingEnergy",
    "BENERG.rgn" => "BindingEnergy.rgn",
];

// Optional callback function to make a prediction directly in the JSON file.
// This function *must* return the $data parameter, including any modifications the
// function may perform.
function make_prediction($data)
{
    if (@$data['BindingEnergy'] <= -5) $data['Predicted'] = 'Agonist';
    else $data['Predicted'] = 'Non-agonist';

    return $data;
}

// The config file contents for bin/primarydock.
// The $pdbfname, $ligname, and $outfname variables have been set by prepare_outputs(),
// however you may modify them if necessary.
$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres
SIZE 7.0 7.0 7.0

EXCL 1 56		# Head, TMR1, and CYT1.

SEARCH TS       # Can use either search method.
POSE 10
ELIM 5          # If no poses are found within this energy limit, the limit will be increased until poses are returned.

FLEX 1          # Enabling flexion increases processing time, but allows better conformer searching.
WET             # Hydrophobic groups stick together as strongly as polar groups do.

ITERS 50        # More iters mean slower processing and better results, albeit with diminishing returns.

OUT $outfname   # Note %p and %l cannot be used for this line because methods_common.php must know the output file name.


heredoc;


// The following will execute the bin/primarydock app with the provided options, then
// save the data specified in $metrics_to_process to a JSON file, in this case named
// predict/dock_results_basic.json reflecting the name of the prediction method.
// If a make_prediction() function has been specified, it will be called before saving
// the dock result so that the changes made by make_prediction() will be included in
// the JSON file.
process_dock();
