
# Overview

PrimaryDock offers a facility for making predictions about receptor responses to ligands. Advances in cryogenic electron
microscopy (cryo-EM) have begun to reveal the structural changes by receptors in response to odorants, making it possible
to estimate odorant activity for each receptor by docking the odorant in the protein's active and inactive states. There
is still a long way to go to be able to make accurate predictions for all combinations of receptor and odorant, and we
continue to strive toward that goal.


# Predictions from the Command Line

The current best way to make predictions is to use the `run_prediction.sh` script in the project root folder, e.g.:

```
./run_prediction.sh OR51E2 propionic_acid
```

The first argument is the receptor ID, and the second is the name of the ligand. After processing finishes, it will output
a result showing the empirically observed activity type (e.g. agonist, inverse agonist, non-agonist) if known, as well as
the predicted activity type and a dock score representing the estimated relative degree of agonism.

Behind the scenes, `run_prediction.sh` calls a prediction method (see section below) which handles such things as making
sure the code is compiled, making sure the PDB and SDF models exist, and maintaining the prediction results for all
receptor-ligand pairs in a single JSON file.

It is also possible to choose which molecular docker to use to make the prediction. Currently the options are PrimaryDock
(`pd`) or AutoDock Vina (`vina`). PrimaryDock is not as fast as Vina, but offers support for ionic and metallic bonds which
Vina lacks. This is important for receptors for acids, amines, and thiols. To specify the docker to use, simply add `pd`
or `vina` as an extra command line argument like so:

```
./run_prediction.sh OR51E2 propionic_acid vina
```


# Prediction Methods

PrimaryDock prediction methods are cronnable PHP scripts that perform docks of empirically measured (or yet to be measured)
ligands in PDB files and aggregate the resulting data into a JSON file. Prediction method scripts start by including
`methods_common.php`, followed by several calls to functions that set up necessary local variables for docking. When the
`process_dock()` function is called, a dock is performed and the difference between active and inactive dock results is
used by the `make_prediction()` function to ascertain whether the odorant is likely to be an agonist for the receptor.

Currently the only prediction method in use is `method_directmdl.php`.

The `directmdl` method, by contrast, is used when a cryo-EM model of the active state is available for either the target
receptor, or one that is deemed similar enough for the cryo-EM model to act as a substitute. Currently, only the OR51,
OR52, and TAAR families are eligible for the `directmdl` method. It works by taking the backbone skeleton of the cryo-EM
model and substituting its side chains with those from the target protein, using Ballesteros-Weinstein numbers to align
the two sequences. Any residue which is conserved in both proteins is not replaced. This results in mostly good results for
SCFA predictions in OR51E2, though predictions for other class I ORs and/or non-acidic ligands may be less optimal.

An example of the command line for running a prediction might be:

```
php -f predict/method_directmdl.php prot=OR1A1 lig=cis-3-hexen-1-ol
```

If the name of the ligand contains spaces, then the spaces should be replaced with underscores, e.g. `lig=isoamyl_acetate`.
If the ligand name contains parentheses, it's a good idea to put the entire lig= parameter in quotes, e.g.:

```
php -f predict/method_directmdl.php prot=OR1A1 "lig=(S)-(+)-carvone"
```


# Alert Sounds

If there is a file named `soundalert` in the `predict/` folder, a sound will play on completion of the prediction,
indicating whether the prediction was correct or not, or in the case of an empirically unknown receptor-ligand pair,
whether the algorithm predicted an agonist. Playing the sound requires the `sox` package to be installed.

To enable the alert sounds, simply create a new text file e.g. `echo "y" > predict/soundalert`. You can also specify
a time range to allow the sounds to play, e.g. `echo "7,21" > predict/soundalert` to only allow sounds between 7:00 and
21:59 local time.

To disable the alert sounds, simply delete the `predict/soundalert` file.


# Writing Your Own Prediction Methods

You may wish to write your own prediction methods. You can make a copy of the `method_directmdl.php` script and then
change the copy as you see fit. If you create a method that makes predictions with at least 80% accuracy, we'd be
delighted to accept your pull request.

Configurable variables include:
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




