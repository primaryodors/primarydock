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

To Use PrimaryDock, please first clone the repository, then execute the following command:

```
make primarydock
```

If you are a developer contributing to the project, you can use `make` to build everything and run the test reports, or 
`make code` to just build the code and run only the amino aldehyde test. This test is critical to the function of PrimaryDock
because any change to the code that causes it to fail, means the docking functionality will be impaired. (If it fails, just
try running `make code` again. Usually it will take a few tries. It is hoped to one day have this test succeed every time.)

The application will require 3D maps of your target receptor(s) in PDB format. Please note that PrimaryDock does not currently
hydrogenate PDB models that do not include hydrogen atoms, so if your model contains heavy atoms only, the accuracy of
docking results may be severely compromised. PDBs for human olfactory receptors are provided in the pdbs folder for olfactory
docking. They have been modified from the PDBs available at the GPCR-I-TASSER website: https://zhanggroup.org/GPCR-I-TASSER/

It will also be necessary to obtain 3D models of your ligand(s). Currently, only SDF format is supported.
SDFs can be obtained a few different ways:
<ul>
  <li>If you have <a href="https://openbabel.org">obabel</a>, you can generate 3 dimensional SDFs from SMILES input. Example syntax:<br>
    obabel -:'CCO' --gen3D -osdf -Oethanol.sdf
  </li>
  <li>SDFs are available from PubChem at https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/SDF?record_type=3d</li>
  <li>Or at https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{SMILES}/SDF?record_type=3d</li>
</ul>

Please take a look at the primarydock.config file as a sample of the format for dock settings. You will be editing this file
(or creating a new one) for each receptor+ligand pair that you wish to dock. There are lines for repointing to your PDB and SDF
model input files, as well as various other options that may be useful to your purposes.

Known issue: please make sure to add an empty line at the end of your config file, or primarydock will throw an exception.

Once your .config file is ready, and the PrimaryDock code is compiled, simply cd to the primarydock folder and run the following command:

./bin/primarydock [config file]

(...replacing "[config file]" with the actual name of your config file.)

After a little while, depending on your config settings, PrimaryDock will output data about one or more poses, including binding energy
per residue, binding energy per type, total binding energy, PDB data of the ligand, and (if flexing is enabled) PDB data of the binding
residues. This output can be captured and parsed by external code, written in your language of choice, for further computation, storage
in a database, etc.

Note if PrimaryDock does not output any poses, please try rerunning it a few times until it gives results. PrimaryDock has pseudo-random
calculations built in so that its output will be different each time, that rerunning the application can catch poses that previous runs
may have missed.

If you would like to contribute to this project:
<ol><li>I would be sooo very grateful for the help!</li>
<li>Please create a branch off of stable, then submit a pull request;</li>
<li>Use whichever { style you prefer; as long as the code is readable and it works, that's all I care about.</li>
<li>Have fun and try not to let the project vex you. (:</li>
</ol>

Note to developers: if you run PrimaryDock under a memory utility such as valgrind, you are likely to see a lot of errors saying that
uninitialized variables are being used or that conditional jumps depend on them. Most of these are false positives. Many places in the
code create temporary arrays of pointers and then assign those pointers addresses of objects that persist throughout the entire program
execution. The memory tool "thinks" the objects have not been initialized even when they have. We recommend using the --undef-value-errors=no
option with valgrind or the equivalent switch in your utility of choice.


# Web Application

You may now host your own PrimaryDock web interface for viewing the contents of the JSON files in the data folder. It is the same web application
as is used for the Primary Odors website.

![Web app screenshot](www/assets/webapp.png?raw=true "Web App")

To enable the web app, either set up a local web server or checkout primarydock in a folder on a web host. Make sure your server has the 
`php` and `php-gd` packages installed. Then open the `www/symlink.sh` file in a text editor, make sure the destination folder is correct (by 
default it will show `/var/www/html/` which is usually correct for Apache2 installations), make sure you have write permissions in the 
folder (or use `sudo`), and execute `www/symlink.sh` in a command line. If on a local server, you will now have an instance of the web app 
at http://127.0.0.1/primarydock/ whereas if you are using a web host then you may have to configure your hosting to point one of your 
registered domains or subfolders to the `primarydock/www` folder.






