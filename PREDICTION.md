
# Overview

PrimaryDock offers a facility for making predictions about receptor responses to ligands. However, since the exact details
of the activation mechanisms of most olfactory receptors are not yet well understood, there is as yet insufficient basis to
make accurate predictions. Nevertheless, we continue to strive toward this goal.


# Predictions from the Command Line

There is a `run_prediction.sh` script in the project root folder that can be used to generate a prediction with minimal
user input, e.g.:

```
./run_prediction.sh OR51E2 propionic_acid
```

The first argument is the receptor ID, and the second is the name of the ligand. After processing finishes, it will output
a result showing the empirically observed activity type (e.g. agonist, inverse agonist, non-agonist) if known, as well as
the predicted activity type and a dock score representing the estimated relative degree of agonism.

Behind the scenes, `run_prediction.sh` calls a prediction method (see section below) which handles such things as making
sure the code is compiled, making sure the PDB and SDF models exist, and maintaining the prediction results for all
receptor-ligand pairs in a single JSON file.


# Prediction Methods

PrimaryDock prediction methods are cronnable PHP scripts that perform docks of empirically measured (or yet to be measured)
ligands in PDB files and aggregate the resulting data into a JSON file. A prediction method script begins with
`require("methods_common.php");`, followed by a call to `prepare_outputs();` which generates the output file names based on
command line input, as well as some other custom variables. A string value named `$configf` is then set, containing the config
file to be generated for the `primarydock` app, and then the `process_dock()` function is called.

The most current prediction method is `method_icactive.php`. The text "icactive" stands for Internal Contacts Activation,
meaning the active state of the receptor is predicted based on contacts made between the side chains of its amino acids.
An example of the command line for running a prediction might be:

```
php -f predict/method_icactive.php prot=OR1A1 lig=cis-3-hexen-1-ol
```

If the name of the ligand contains spaces, then the spaces should be replaced with underscores, e.g. `lig=ethyl_vanillin`.
If the ligand name contains parentheses, it's a good idea to put the entire lig= parameter in quotes, e.g.:

```
php -f predict/method_icactive.php prot=OR1A1 "lig=(S)-(+)-carvone"
```

The prediction method lends itself well to running as a cron, e.g.:

```
* * * * * php -f /path/to/primarydock/predict/method_icactive.php next simul=8
```

In this example, the `next` parameter means to find the next receptor/ligand pair yet to process, and go ahead with the
calculations for it. The "next" pair to process is the first pair of GPCR + known ligand (whether agonist or not) that
either is not present in `predict/dock_results_optimized.json` or has a `version` timestamp older than either: the
prediction method .php, `methods_common.php`, or the `/bin/primarydock` executable.

The `simul=8` parameter indicates not to run more than 8 concurrent processes. We recommend a value of no more than
the number of processor cores or threads on your machine, so if you have a server with one 8-core 16-thread processor,
then `simul=8` would be the maximum recommended value. We strongly recommend installing `lm-sensors` if your system does
not already have this tool (`sudo apt-get install lm-sensors`), that way the prediction script will automatically limit
itself depending on the CPU temperature.

You can use the `predict/progress.sh` shell script to monitor the progress of the predictions.


# Writing Your Own Prediction Methods

Until an accurate prediction method exists, we encourage you to write your own prediction methods. You can use the
`method_basic.php` script as a template, and then apply any variations you wish. If you create a method that makes
predictions with at least 80% accuracy, we'd be delighted to accept your pull request.

Configurable variables include:
- `$dock_metals`          If true, the dock will use PDBs that end in `.metal.pdb` instead of the default `.upright.pdb`.
                          You probably won't ever have to use the `$dock_metals` feature.
- `$dock_retries`         The number of times to retry if no poses returned. Each retry increases the energy limit.
                          Default = 5.
- `$bias_by_energy`       TODO: This has not been implemented yet.
- `$metrics_to_process`   An array of key-value pairs identifying a dock metric and its new key in the output JSON.
                          Valid keys can be found in the table below; the value becomes the metric's key in the JSON.

Available parameters for `$metrics_to_process` include:
- `BENERG`                Total binding energy between ligand and protein;
- `BENERG.rgn`            Total binding energy between ligand and specific regions of the protein, e.g. TMR6;
- `BEST`                  Total ligand-protein binding energy of the most favorable pose.
- `vdWRPL`                Total van der Waals repulsion (i.e. where atoms are < 4Å apart) ligand to protein;
- `vdWRPL.rgn`            Total vdW repulsion between ligand and regions of protein;
- `POLSAT`                Ligand's sum polar satisfaction, a measure of polar-polar and nonpolar-nonpolar proximity;
- `PCLASH`                Clashes between residues of the protein, minus clashes in the original PDB model;
- `ACVTH.TMR1` thru `ACVTH.TMR7`
                          Angles of helix rotations performed for a "soft dock", in which the helices are allowed to move.
- `POSES`                 Total number of poses found.
- `CALOC.rgn`             Average relative motion of alpha carbon atoms within each region as a result of a "soft dock".

For keys ending with `.rgn`, the `rgn` will be replaced with the actual region name. For example, if the method PHP
specifies `BENERG.rgn` => `energy.rgn`, then the JSON will contain metrics labeled `energy.TMR3`, `energy.TMR6`, etc.




