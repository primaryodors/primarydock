
# Predicting Receptor Responses to Ligands

At this time, no method yet exists to accurately identify agonist and non-agonist ligands from PrimaryDock output.
To the best of our knowledge, either no one has developed this capability, or they have not released it to the public.
Nevertheless, we continue to strive towards accurate predictions of GPCR agonism by odorants.

# Prediction Methods

PrimaryDock prediction methods are cronnable PHP scripts that perform docks of known ligands in PDB files and aggregate
the resulting data into a JSON file. A prediction method begins with `require("methods_common.php");`, followed by a
call to `prepare_outputs();` which generates the output file names as well as some other custom variables. A string
value named `$configf` is then set, containing the config file to be generated for the `primarydock` app, and then the
`process_dock()` function is called.

The most current prediction method is `method_coupled.php`. To use this prediction method, one must first set up the
required data.

1. First, run `bin/pepteditor predict/gprot.pepd` to download the experimental models from RCSB.
2. Next, the PHP script `predict/generate_couple.php` must be run for every odor receptor. A good way to do this is
    using a cron, e.g.:
    ```
    * * * * * php -f /path/to/primarydock/predict/generate_couple.php next simul=8
    ```
3. Finally, the prediction method can be run for all known receptor-ligand pairs. This also lends itself well to
    running as a cron, e.g.:
    ```
    * * * * * php -f /path/to/primarydock/predict/method_coupled.php next simul=8
    ```

For both scripts, the `next` parameter means find the next GPCR/G-protein pair or GPCR/ligand pair yet to process, and 
process it.

For `generate_couple.php`, the "next" pair to process is the first combination of GPCR/G-protein that does not yet
have a PDB file in the `pdbs/coupled/` folder.

For `method_coupled.php`, the "next" pair to process is the first pair of GPCR + known ligand (whether agonist or not)
that either is not present in `predict/dock_results_coupled.json` or has a `version` value older than either the
prediction method PHP, `methods_common.php`, or the `/bin/primarydock` executable.

The `simul=8` parameter indicates not to run more than 8 concurrent processes. We recommend a value of no more than
the number of physical cores on your machine, so if you have a server with one 8-core processor, then `simul=8` would
be the maximum recommended value. If you have `lm-sensors` installed (`sudo apt-get install lm-sensors`), then both PHP
scripts will automatically limit themselves depending on the CPU temperature.

You can use the `predict/progress.sh` shell script to monitor the progress of both PHP scripts.

# Writing Your Own Prediction Methods

Until an accurate prediction method exists, we encourage you to write your own prediction methods. You can use the
`method_basic.php` script as a template, and then apply any variations you wish. If you create a method that makes
predictions with at least 80% accuracy, we'd be delighted to accept your pull request.

Configurable variables include:
`$dock_metals`          If true, the dock will use PDBs that end in `.metal.pdb` instead of the default `.upright.pdb`.
                        You probably won't ever have to use the `$dock_metals` feature.
`$dock_retries`         The number of times to retry if no poses returned. Each retry increases the energy limit.
                        Default = 5.
`$bias_by_energy`       TODO: This has not been implemented yet.
`$metrics_to_process`   An array of key-value pairs identifying a dock metric and its new key in the output JSON.
                        Valid keys can be found in the table below; the value becomes the metric's key in the JSON.

Available parameters for `$metrics_to_process` include:
`BENERG`                Total binding energy between ligand and protein;
`BENERG.rgn`            Total binding energy between ligand and specific regions of the protein, e.g. TMR6;
`vdWRPL`                Total van der Waals repulsion (i.e. where atoms are < 4Å apart) ligand to protein;
`vdWRPL.rgn`            Total vdW repulsion between ligand and regions of protein;
`POLSAT`                Ligand's sum polar satisfaction, a measure of polar-polar and nonpolar-nonpolar proximity;
`PCLASH`                Clashes between residues of the protein, minus clashes in the original PDB model;

`ACVTH.TMR1` thru `ACVTH.TMR7`
                        Angles of helix rotations performed for a "soft dock", in which the helices are allowed to move.

`POSES`                 Total number of poses found.
`CALOC.rgn`             Average relative motion of alpha carbon atoms within each region as a result of a "soft dock".

For keys ending with `.rgn`, the `rgn` will be replaced with the actual region name. For example, if the method PHP
specifies `BENERG.rgn` => `energy.rgn`, then the JSON will contain metrics labeled `energy.TMR3`, `energy.TMR6`, etc.




