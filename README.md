# PrimaryDock
PrimaryOdors.org molecular docker.<br>
http://www.primaryodors.org

PrimaryDock is a lightweight stochastic molecular docking software package that offers the following advantages:
- Path-based docking;
- Native support for side-chain flexion;
- Per-residue binding strength output;
- Per-binding-type binding strength output;
- Per-residue van der Waals proximity repulsion output;
- Small self contained codebase with no extraordinary dependencies;
- Does not require CMake, but can be built on any recent *nix system using make, g++, and the C++14 Standard Library;
- Interatomic parameters stored in flat text files that can be edited without recompiling the application.

PrimaryDock comes with Pepteditor, a tool for editing proteins using a scripting language.

PrimaryDock has a prediction feature that attempts to predict receptor responses to odorants.
The prediction feature requires php and <a href="https://openbabel.org">OpenBabel</a> to be installed.
See `PREDICTION.md` for more info.

PrimaryDock also offers a web interface that allows you to run a local copy of the same data explorer pages that power
the PrimaryOdors website. The web interface also requires OpenBabel.

To use PrimaryDock and Pepteditor, please first clone the repository, then enter the `primarydock/` folder and execute
the following command:

```
make
```

The application will require 3D maps of your target receptor(s) in PDB format. If your model contains heavy atoms only,
the accuracy of docking results may be severely compromised. Fortunately, PrimaryDock can hydrogenate PDBs automatically
during docking, or you can use Pepteditor to hydrogenate them before docking. PDBs for human olfactory receptors are provided
in the `pdbs/` folder for olfactory docking. They have been modified from the PDBs available from AlphaFold.

It will also be necessary to obtain 3D models of your ligand(s). Currently, only SDF format is supported.
SDFs can be obtained a few different ways:
<ul>
  <li>Using openbabel, you can generate 3 dimensional SDFs from SMILES input. Example syntax:<br>
    <code>obabel -:'CCO' --gen3D -osdf -Osdf/ethanol.sdf</code><br>
    There is also a PHP script located at `data/gensdf.php` that will make the `obabel` call, as well as add the molecule to
    `odorants.json`, if given a compound name and a SMILES string.
  </li>
  <li>SDFs are available from PubChem at https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/SDF?record_type=3d</li>
  <li>Or at https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{SMILES}/SDF?record_type=3d</li>
</ul>

Please take a look at the `primarydock.config` file as a sample of the format for dock settings. You will want to edit this file,
or create a new one, for each receptor+ligand pair that you wish to dock. There are lines for repointing to your PDB and SDF
model input files, as well as various other options that may be useful to your purposes.

Known issue: please make sure to add an empty line at the end of your config file, or primarydock will throw an exception.

Once your .config file is ready, and the PrimaryDock code is compiled, simply cd to the primarydock folder and run the following command:

```
./bin/primarydock {config file}
```

(...replacing `{config file}` with the actual name of your config file.)

After a little while, depending on your config settings, PrimaryDock will output data about one or more poses, including binding energy
per residue, binding energy per type, total binding energy, PDB data of the ligand, and (if flexing is enabled) PDB data of the flexed
residues and binding residues. This output can be captured and parsed by external code written in your language of choice, for further
omputation, storage in a database, etc.

Most of the time, PrimaryDock will find poses with favorable energy levels because of the Best-Binding algorithm, which seeks to match
binding pocket features with compatible ligand features. However, there may be a few cases when a ligand is a poor fit for a binding
pocket, and PrimaryDock might find zero output poses. If this happens, you can try increasing the energy limit or performing a soft dock.
PrimaryDock is stochastic so that its output will be different each time, and rerunning the application can often catch poses that
previous runs may have missed.

Contributions are always welcome! Please create a branch off of stable, then submit a pull request.
All PRs that change the C++ classes or any of the apps must pass the master unit tests (`test/unit_test_master.sh`) before merge.

Note to developers: if you run PrimaryDock under a memory utility such as Valgrind, you are likely to see a lot of errors saying that
uninitialized variables are being used or that conditional jumps depend on them. Most of these are false positives. Many places in the
code create temporary arrays of pointers and then assign those pointers addresses of objects that persist throughout the entire program
execution. The memory tool "thinks" the objects have not been initialized even when they have. We recommend using the 
`--undef-value-errors=no` option with valgrind or the equivalent switch in your utility of choice.


# Utilities

A few utility apps are also provided in the bin/ dir.


# cavity_search

Scans a protein's 3 dimensional structure looking for places where a ligand may be able to dock. An output file can be specified to receive
the coordinates of collections of spheres that form the shapes of the cavities.

If the output file is in the same directory as a PDB model, and its name is identical except for having a .cav extension instead of .pdb,
then certain 3D views of the web app will recognize the file and allow the user to see the found pockets in the protein's 3D structure.

The prediction feature automatically generates .cav files for the active and inactive models used in the prediction, if no .cav file already
exists when the prediction begins.


# ic

Finds internal contacts within a protein. Simply pass it the pathname of a PDB file and it will output a series of contacts between residues
in that model.


# pepteditor

A tool for editing protein models using a scripting language. See the `SCRIPTING.md` file for more information.


# ramachandran

Allows visualizing Ramachandran plots of protein models. Takes the pathname of a PDB file as its only required argument.

Optionally you can specify the -n parameter to output the plot as numbers instead of colorized blocks.


# ringflip

Performs flips of atoms in flexible rings. The canonical example would be the conversion of cyclohexane rings between chair form and boat form.


# score_pdb

Analyzes a .pdb file containing HETATM records and scores the interactions between ligand atoms and protein atoms. Internally, it uses the same
code as the scoring function of `primarydock`, however this code depends on certain details not stored in the PDB format, so the results may
vary from the output of `primarydock`.


# Web Application

You may optionally host your own PrimaryDock web interface for viewing the contents of the JSON files in the data folder. It is the same web
application that is used for the Primary Odors website.

![Web app screenshot](www/assets/webapp.png?raw=true "Web App")

To enable the web app:
- Either set up a local web server or checkout primarydock in a folder on a web host.
- Make sure your server has the `php`, `php-curl`, `php-gd`, and `openbabel` packages installed.
- After installing `php-curl`, it's important to restart the web service e.g. `sudo apache2ctl -k restart`.
- Then open the `www/symlink.sh` file in a text editor, make sure the destination folder is correct (by default it will show `/var/www/html/`
  which is usually correct for Apache2 installations), make sure you have write permissions in the 
  folder (or use `sudo`), and execute `www/symlink.sh` in a command line.
- The `data` and `www/assets` folders and all contents must also be recursively made writable by the web user.
- If on a local server, you will now have an instance of the web app at http://127.0.0.1/primarydock/ whereas if you are using a web host
  then you may have to configure your hosting to point one of your registered domains or subfolders to the `primarydock/www` folder.

If you get a 403 Forbidden error, please make sure that every containing folder of the `primarydock/www` folder has public execute access.


# Adding Data

To add a new receptor protein to the PrimaryDock database, there is a series of steps and utilities that facilitate this process:

- Add the protein to `data/receptor.json`. It must have an "id": and a "sequence":, and should also have a "uniprot_id":.
- Add the ID and sequence to `data/sequences_aligned.txt`, in alphanumeric order with related proteins, and manually align the sequence with dashes.
- In a command line, run `php -f data/sequence_update.php` and then `php -f data/btree.php`.
- Next, run `php -f www/getpdbs.php` and look out for lines similar to `Wrote pdbs/***/*****.upright.pdb.` in the output.
- Optionally, if you have MODELLER installed, you can generate an active-state PDB for the new protein:
  - Firstly, align the new protein's sequence to the sequences in `hm/experimental.ali`. *Note this is not the same alignment as in `data/sequences_aligned.txt`.*
  - Next, add the new alignment as a node called `"aligned"` in the new protein's record in `data/receptor.json`.
  - Then run `php -f hm/build_alignment_file.php` to generate the alignment file for all included GPCRs.
  - Finally, run `php -f hm/dohm.php RCPID` changing `RCPID` to the ID of the new protein.
- Create a new branch if necessary and `git add -f` each of the `*.upright.pdb` (and `*.active.pdb` if present) files from the getpdbs output.
- Check in the new and updated files and create a pull request.

