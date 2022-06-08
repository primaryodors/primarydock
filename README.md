# podock
PrimaryOdors.org molecular docker.<br>
http://www.primaryodors.org

PODock is a lightweight molecular docking software package that allows path-based docking, side-chain flexion, and binding strength
measurement at the residue level and per binding type. It uses data files so atomic and interatomic parameters can be fully customized
without having to recompile the code.

To Use PODock:

First clone the repository, then execute the following commands:

For Linux:
```
git checkout stable
make
```

The application will require 3D maps of your target receptor(s) in PDB format. Please note that PODock does not currently
hydrogenate PDB models that do not include hydrogen atoms, so if your model contains heavy atoms only, the accuracy of
docking results may be severely compromised.

PDBs for most olfactory receptors are available here:<br>
http://primaryodors.org/cached/

It will also be necessary to obtain 3D models of your ligand(s). Currently, only SDF format is supported.
SDFs can be obtained a few different ways:
<ul>
  <li>If you have <a href="https://openbabel.org">obabel</a>, you can generate 3 dimensional SDFs from SMILES input. Example syntax:<br>
    obabel -:'CCO' --gen3D -osdf -Oethanol.sdf
  </li>
  <li>SDFs are available from PubChem at https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/SDF?record_type=3d</li>
  <li>Or at https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{SMILES}/SDF?record_type=3d</li>
</ul>

Please take a look at the podock.config file as a sample of the format for dock settings. You will be editing this file
(or creating a new one) for each receptor+ligand pair that you wish to dock. There are lines for repointing to your PDB and SDF
model input files, as well as various other options that may be useful to your purposes.

Known issue: please make sure to add an empty line at the end of your config file, or PODock will throw an exception.

Once your .config file is ready, and the PODock code is compiled, simply cd to the podock folder and run the following command:

./bin/podock [config file]

(...replacing "[config file]" with the actual name of the file you edited or created.)

After a little while, depending on your config settings, PODock will output data about one or more poses, including binding energy
per residue, binding energy per type, total binding energy, PDB data of the ligand, and (if flexing is enabled) PDB data of the binding
residues. This output can be captured and parsed by external code, written in your language of choice, for further computation, storage
in a database, etc.

If you would like to contribute to this project:
<ol><li>I would be sooo very grateful for the help!</li>
<li>Please create a branch off of main, then submit a pull request;</li>
<li>Use whichever { style you prefer; as long as the code is readable and it works, that's all I care about.</li>
<li>Have fun and try not to let the project vex you. (:</li>
</ol>

Note to developers: if you run podock under a memory utility such as valgrind, you are likely to see a lot of errors saying that
uninitialized variables are being used or that conditional jumps depend on them. These are false positives. Many places in the code
create temporary arrays of pointers and then assign those pointers addresses of objects that persist throughout the entire program
execution. The memory tool "thinks" the objects are uninitialized even when they are. We recommend using the --undef-value-errors=no
option with valgrind or the equivalent switch in your utility of choice.


