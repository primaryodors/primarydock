
# Pepteditor Documentation

Pepteditor: The Scripted Peptide Editor.

To run a Pepteditor script, after building PrimaryDock, please use the following command:

```
bin/pepteditor path/to/script.pepd
```


# Strands

The PDB format allows each amino acid residue to be given a strand ID, typically a letter between `A` and `Z`.
Pepteditor allows holding up to 26 strands in memory at any given time.
One strand is always the current working strand; the default is `A`.
All commands apply only to the working strand unless otherwise noted.
It is possible to load data from multiple PDBs into multiple strands, move protein data between strands, and delete strands.
When saving an output PDB, all strands are included in the written file.


# Variables

With few exceptions, any variable can stand in for any parameter of any command.
Variable names are case sensitive.

All integers begin with `%`, e.g. `%start_resno`.
All floats begin with `&`, e.g. `&transpose_amt`.
All Cartesian 3D locations begin with `@`, e.g. `@pocket_center`.
All strings begin with `$`, e.g. `$name`.

Cartesians have members .x, .y, and .z that behave as floats.

Casting a float to an integer rounds the value (e.g. 0.4 rounds to zero but 0.5 rounds to 1).

Casting a float to a Cartesian normally sets the .x member to the float value, leaving .y=0 and .z=0.
But a command like `LET @foo = @bar.y` will result in only the value of `@foo.y` being nonzero.

Casting an integer to a Cartesian obtains the location of the CA atom for that residue number, if it exists in the working strand.
But note that e.g. `LET @foo = %motif + 2` will not work. The arithmetic must come first e.g. `LET %foo = %motif + 2` followed by
`LET @foo = %foo`.

Casting a Cartesian back to float or integer obtains the magnitude of the Cartesian, equal to sqrt(x^2 + y^2 + z^2).

The command line arguments are made available to the script as $arg1, $arg2, etc. Normally, $arg1 will be the .pepd script filename.


Array functionality is available, though not implemented as true arrays. Each element is internally stored as its own discrete variable, for example
if we define `LET $array[3] = "three"` then a string variable called `$array[3]` will be assigned the value of `three`. That means any reference to
`$array` will obtain an empty string since that variable technically hasn't been set.

Variables as array indices are possible because of nested variable names. If an integer `%i` is set, and has a value of `3`, then `$array[%i]` and `$array[3]` will be equivalent. Nested variables are limited to integer and string types, so while `$array[$key]` and `@array[$key]` are allowed, 
`$array[@location.x]` is not. If there is any ambiguity in variable names, the longest and earliest matching variable is used.

This works similarly to the `$$` syntax of PHP, where if `$varname = "foo"` then `$$varname` is the same as `$foo`. However, unlike in PHP, the nested
variable name can occur anywhere in a pepteditor variable, there can be more than one nested variable, and variable names can be nested multiple levels
deep, so for example `$array[%i][%array2[%j]]` is legal pepteditor syntax as long as `%i`, `%j`, and `%array2[%j]` have all been set.

This also means that square brackets are not required for array functionality, and one could just as easily define variable `$array[3]` as `$array3`,
`$array(3)`, or `$arrayfoo3bar`, and the syntax of `%i` (or any other int variable) in place of `3` will still work.

Example of how nested variable names work:

```
LET %sub = 3

LET $array[3] = "three"
ECHO $array[%sub]         # outputs three

LET $$array[3] = "Yup!"   # same as LET $three = "Yup!" because $array[3] = "three".
ECHO $$array[%sub]        # outputs Yup!
ECHO $three               # outputs Yup!
```


Some "magic variables" are supplied upon loading a protein. The strand ID is included in the variable name; for the below list, the strand is A:

- `$PDBA` the path and name of the source PDB file.
- `$PROTEINA` the name of the protein, derived from a `REMARK 6` record if present.
- `%SEQLENA` the length of the protein sequence.
- `$SEQUENCEA` the sequence of the protein in standard one-letter amino acid code.
- Region start and end residue numbers from any `REMARK 650 HELIX` records, e.g. `%A.TMR3.s` and `%A.TMR3.e` for the start and end resnos of TMR3.
- Ballesteros-Weinstein n.50 numbers e.g. `%A.1.50`, `%A.2.50` etc., from `REMARK 800 SITE BW` records of the PDB if present.

The following commands are supported:


# ALIGN
Example:
```
ALIGN 174 182 4 @location1 @location2
ALIGN 1 9999 4 174 @location1 182 @location2
```

Repositions a piece of the protein so that the residues near the start resno occur at one location, and the residues near the end resno occur in the
direction of the other location.

The first two numbers are the starting and ending residue numbers. They define the range of residues that will be moved.

The next number is the averaging length from each end of the region. In the example this number is 4, so the locations of the CA atoms of four residues
at each end will be averaged together, in this case 174-177 and 179-182, to determine the alignment. The piece will be moved so that the average of CA
atom locations of residues 174-177 will now equal @location1, and the average of CA locations of residues 179-182 will fall along a straight line
pointing from @location1 to @location2.

The second example aligns residues 174-177 and 179-182 with the same two points, but affects the entire protein instead of just the region between
the indicated residues.


# ATOMTO
Example:
```
ATOMTO %residue_number $atom_name @target
```

Flexes the side chain of `%residue_number` to bring `$atom_name` as close to `@target` as physically possible.


# BEND
Example:
```
BEND 206 218 "N-CA" 15
```

Effects a partial bend of the protein backbone at the start residue (first parameter) by rotating in the bond direction (third parameter) by the rotation
angle (fourth parameter, degrees), moving all subsequent residues up to the end residue (second parameter). The result is a bend in the chain
corresponding to either a phi or psi rotation of the start residue.

Note that the end residue will be disconnected from its neighbor and will require to be reconnected.

Valid bond directions are:
N-CA		Keep the residue's N terminus stationary and rotate its C terminus and subsequent residues about the N-CA (phi) bond;
CA-C		Keep the residue's N terminus stationary and rotate its C terminus and subsequent residues about the CA-C (psi) bond;
CA-N		Keep the residue's C terminus stationary and rotate its N terminus and preceding residues about the CA-N (phi) bond;
C-CA		Keep the residue's C terminus stationary and rotate its N terminus and preceding residues about the C-CA (psi) bond;

Note if bending in the N-CA or CA-C direction, i.e. towards higher numbered residues, then the end residue must be greater than the start residue.
If bending in the C-CA or CA-N direction, the start residue must be greater than the end residue.

See here for an introduction to phi and psi bonds:
https://proteopedia.org/wiki/index.php/Tutorial:Ramachandran_principle_and_phi_psi_angles


# BENERG
Example:
```
BENERG 251 111 &binding
```

Reads the non-covalent energy level (negative for binding, positive for clashes) in kJ/mol between two protein residues, and writes that value to a
specified float variable.


# BRIDGE
Example:
```
BRIDGE 251 111
BRIDGE 251 111 50
```

Iteratively conform the side chains of two residues (parameters 1 and 2) to maximize their non-covalent bonding to each other. If the optional third
parameter is given, it represents the number of iterations (default 50). More iterations give better results, but takes longer to process, and diminishing
returns occur with large iteration values.


# BWCENTER
Example:
```
BWCENTER
```

If the current working strand is a seven-helix protein (7HP), with Ballesteros-Weinstein numbering for all its transmembrane helices,
then `BWCENTER` will obtain the center of all n.50:CA atom locations and recenter the entire protein.
This is useful for comparing various active and inactive states of GPCRs; the 1.50 through 7.50 residues will be assumed to be the most stationary
parts of the protein, and variations in structure can be observed with minimal global transformational anomalies.


# CENTER
Example:
```
CENTER
CENTER [0,0,0]
```

Centers the entire protein at the optional specified coordinates, or [0, 0, 0] if omitted.


# CONNECT
Examples:
```
CONNECT 160 174
CONNECT 160 174 250
```

Flexes a section of protein backbone to try to reconnect a broken strand, for example after a HELIX and/or ALIGN command.

The first parameter is the starting residue number of the section to be flexed. The second parameter is the target residue to reconnect to. In these
examples, the program will flex the range of 160-173 and attempt to reunite 173 with 174. The last optional parameter is the number of iterations, the
default being 50.


# CTNRG
Examples:
```
CTNRG A B &energy
CTNRG A B &energy @direction
CTNRG A %start1 %end1 B %start2 %end2 &energy @direction
CTNRG A B %start2 %end2 &energy @direction
```

Gets the contact energy between two strands. In these examples, the interaction energy between strands A and B will be stored in `&energy`.
Negative values mean favorable interactions, expressed in kJ/mol. Positive values mean atom clashes.
The A and B parameters can also be strings.

The second example introduces a direction of motion to move the second strand (in this case, B strand) to optimize the contacts between proteins.
This is useful if e.g. two proteins are too close together and clashing, or too far away and not making good contact.

The third example limits the function to a specific region of each strand. This is useful to for example optimize only part of the second strand.
The fourth example limits only the second strand to a region. If a region is specified for only one strand, it must be the second one.
This may seem counterintuitive but it makes sense for optimizing only part of the second protein against the whole of the first protein.

The direction parameter is always optional, irrespective of any start and end residue numbers given.


# DELETE
Example:
```
DELETE 1 57     # Remove all residues before number 58.
```

Deletes a range of residues from the model of the protein in memory.


# DISULF
Example:
```
DISULF
DISULF 97 179
```

Creates disulfide bonds between residues that have sulfhydryl side chains, such as cysteine or homocysteine. If called with two parameters,
`DISULF` will attempt to join the two residues and generate a warning if the attempt fails. If called without parameters, the command will try
all combinations of sulfhydryl side chains, and output only the ones it was able to successfully join.

Note that for `DISULF` to work, the two sulfur atoms must already be moved into position, and the side chains must be hydrogenated. No attempt
will be made to flex bonds to bring the sulfurs into proximity, and if the sulfurs are not already bound to hydrogens, the command will fail.


# DOWNLOAD
Example:
```
DOWNLOAD RCSB 8F76 "pdbs/8f76.pdb"
DOWNLOAD AF Q5QD04 "pdbs/mTAAR9.pdb" ONCE
```

Downloads a protein model from one of the download sources registered in the `data/dlsrc.dat` file.
The second parameter is the source's unique identifier for the model, e.g. 8F76 for the RCSB model of propionate-activated hOR51E2 on GNAS2,
or UniProt ID Q5QD04 for mTAAR9 from AlphaFold.
The third parameter is the destination file name.
ONCE is an optional final parameter that tells Pepteditor not to request the source URL if the destination file already exists.


# DUMP
Example:
```
DUMP
```

Useful for debugging, this outputs a list of all variables and their values to stdout. It does not take any parameters.


# ECHO
Examples:
```
ECHO $PROTEIN
ECHO $SEQUENCE " is " %SEQLEN " residues long."
ECHO @location ~
```

Writes information to stdout.

Unlike the `LET` command, `ECHO` does not require a plus sign for concatenation.

If a tilde `~` occurs at the end of the statement, no newline will be output. Otherwise, a newline will be added to whatever was echoed.


# END / EXIT / QUIT
Examples:
```
END
EXIT
QUIT
DIE
EXIT -1
DIE "An error has occurred."
```

Ceases script execution.

Optionally, an integer return value may be passed, or a string to echo, but not both.


# GEN
```
GEN "MAYDRYVAIC"
```

Creates a peptide using the specified sequence.


# GOSUB
Example:
```
GOSUB label
...
label:
...
RETURN
```

Like `GOTO` (see below), except the label begins a subroutine that ends with a `RETURN` command. Upon finishing the subroutine, script execution
resumes on the line following the `GOSUB`.


# GOTO
Example:
```
GOTO label
...
label:
```

Diverts script execution to the specified label. A line must exist in the script that consists of only the label and a colon `:` character, with
no spaces and no other characters.


# HELIX
Examples:
```
HELIX ALPHA 174 182
HELIX -57.8 -47.0 174 182
```

Create a helix of the specified residue range. If the helix type is specified, it cannot be a variable and must be one of: `ALPHA`, `PI`, `3.10`,
`PPRO1` (polyproline I), `PPRO2` (polyproline II), `BETA`, or `STRAIGHT`. Note that while the last two aren't helices, they use the same syntax
because like with helices, the residues' phi and psi bonds are rotated to preset angles.

Phi and psi angles can also be manually specified as in the second example.


# HYDRO
Examples:
```
HYDRO
```

Hydrogenates the protein.


# IF
Examples:
```
IF %var >= 10 LET %var = 0
IF %var THEN ECHO "Yep."
IF &var1 > 0 AND &var1 < &var2 SAVE $name QUIT

IF $match = "HFFCE" ECHO $message
ELSE ECHO "Nope."

LET %iter = 1
loop:
ECHO %iter
LET %iter ++
IF %iter <= 10 GOTO loop

IF EXISTS $filename ECHO "File exists."
```

The `IF` command evaluates a conditional expression and, if true, executes another command. Currently, only simple A = B type expressions are
supported, i.e. straightforward comparison of two values, where either value can be a single variable or a single constant; such expressions
as A > B + 1 do not yet work.

`IF` can be used to test a single variable, as in the `IF %var THEN` statement in the examples. For this syntax, the `THEN` keyword is required,
otherwise the interpreter would treat the following statement (in this case `ECHO`) as an unrecognized operator and throw an error. This also
means `IF %foo AND %bar` will not work; the syntax instead would have to be `IF %foo THEN IF %bar THEN` [statement]. All other times, the `THEN`
keyword is optional.

`AND` and `OR` are supported, however parentheses are not, and while an unlimited number of `AND`s and/or `OR`s can be chained together, their
evaluation will be entirely sequential and the chain will evaluate false as soon as any chain of `OR`ed expressions all evaluate false. Note that
any `ELSE` after such a statement will execute no matter which part of the chain evaluated to false.

The operators available for `IF` are `=` `!=` `>` `<` `>=` `<=`. Note that `=` means comparison, not assignment. The operator `==` is a synonym
of `=` for purposes of `IF` statements. Strings also support the `=*` operator, which means "contains", e.g. `"ABCDE" =* "BCD"` evaluates true.

If the optional `ELSE` subcommand is provided, it occurs on the next line. `ELSE` without `IF`, or `ELSE` separated from its `IF` even if only by
a blank line or a comment, will silently fail to execute so please be careful.

There is no limitation to which commands can be paired with `IF` or `ELSE` (keeping in mind that `ELSE` is not itself a command). It is possible
to write a series of `ELSE IF`s, or to combine `IF`s, e.g. `IF &var1 > 0 IF &var1 < &var2 GOTO label`, in fact the `AND` keyword is implemented
in exactly this way by transparently replacing `AND` with `IF` behind the scenes.

When using `IF` with `GOTO`, it is possible to create loops as seen in the last example above.


# LET
Examples:
```
LET %i = 1
LET $name = "PrimaryDock"
LET $range = $SEQUENCE FROM 174 FOR 20
LET $sub = $name FROM 3 FOR 4
LET &j = %i
LET @loc = [5,3,-8]             # Literal XYZ coordinates. Do not include spaces!
LET @loc = [5,&j,%i]
LET @res174loc = 174            # Gets location of residue 174's alpha carbon.
LET &j += @res174loc.y
LET @loc3 = @loc1 + @loc2
LET $message = "TMR4 ends on residue " + %TMR4.e + " and TMR5 starts on residue " + %TMR5.s + "."
```

Assigns a value to a new or existing variable.

For assigning integer and float variables, the following operators are allowed: `=` `+=` `-=` `*=` `/=`. Integer assignments also support the bitwise
`&=` (and) and `|=` (or) operators.

Integer variables are also allowed `++` `--`, which do not take an r-value (a variable or literal value to apply to the variable assignment).

For algebraic assignments, e.g. `LET &k = %i + &j`, the following operators are supported: `+` `-` `*` `/` `^`, the last being an exponent operator.
Integers also support the bitwise `&` (and) and `|` (or) operators. Note parentheses are not currently supported. Also, the data type of the variable
being assigned determines how all r-value variables are interpreted, so even though `%i` is an integer in the above statement, its value is cast to a 
float (the data type of `&k`) before the addition takes place.

For strings, the operators `=` `+=` and `+` are allowed, with the plus sign indicating concatenation.
Substrings are _one-based_ and are indicated with the word FROM and (optinally) FOR, e.g. `LET $end = $string FROM 10` or `LET $range = $SEQUENCE FROM %start FOR %length`.


# LOAD
Example:
```
LOAD "path/to/protein.pdb"
LOAD "path/to/protein.pdb" A
LOAD "path/to/protein.pdb" A B
```

Loads the specified protein into working memory. Note only one protein is allowed in memory at a given time.
If one strand ID is provided, then that strand of the PDB will be loaded. The default strand is A.
If a second strand ID is provided, then the protein will be loaded into that strand locally and the current working strand will be updated.
The default is to load the protein into the existing working strand ID.


# MCOORD

Example:
```
MCOORD Ag 1 %res1 %res2 %res3 %res4			# Coordinate Ag+ to the residues indicated by the variables.
MCOORD Th8 Cu 2 123 154 203 259 			# Coordinate Cu2+ to the indicated residues, forming thiolates of any cysteines or selenolates of selenocysteines.
MCOORD YAr Zn 2 123 234 356 YO 467			# Coordinate Zn2+, cation-pi bonding any tyrosines except the last listed residue.
```

Coordinate a metal ion to the indicated residues. A minimum of 3 residues is required.

The optional Th8 parameter indicates that the metal ion should be thiolated to any CYS or SEC residues in the list, i.e. the metal will be covalently
bonded to the S or Se, replacing the residue's H. TODO: This feature has not been implemented yet.

Tyrosine is treated as both an aromatic residue and an oxygen coordinating residue. Aromatic residues are capable of forming cation-pi bonds, therefore
phenylalanine and tryptophan will also be cation-pi coordinated to the metal if found. The default coordination method for tyrosine is cation-pi, however
the default can be overridden with YO to force oxygen coordination or YAr to force pi coordination.

Yar or YO before the metal symbol applies globally to the coordination residues. Yar or YO directly before a residue number applies only to that residue,
overriding the default or any global.


# MOVE
# MOVEREL

Example:
```
MOVE %start_res %end_res @newcen
MOVE 160 176 [0,25,0]
MOVEREL 160 176 [0,0,1]
```

`MOVE` Recenters the indicated region at the indicated Cartesian location.
`MOVEREL` adds the specified relative Cartesian value to the region's location.


# PTALIGN

Example:
```
PTALIGN @point @target @origin @normal &angle
```

Finds the rotation axis and angle necessary to rotate `@point` about `@origin` to get as close as possible to `@target`. The `@normal` parameter
is an output parameter; a variable must be used and its value will be set to the rotational axis _relative to the origin_. Note this is not a point
in absolute space. The final parameter `&angle` is an output parameter that receives the necessary rotation angle in degrees.


# PTROTATE

Example:
```
PTROTATE @point @origin @axis &angle @result
```

Rotates `@point` in space about `@origin`, using `@axis` as a _relative_ rotation axis (i.e. add `@axis` to `@origin` to get an actual point in space
along this imaginary line) by `&angle` degrees and stores the result in the output parameter `@result`.


# REGION
Example:
```
REGION TMR5 202 238
REGION "TMR5" 202 238
REGION $region_name %start_res %end_res
```

Defines a region in the protein object, overwriting any existing region with the same name. Also creates the %region.s and %region.e variables.


# RENUMBER
Example:
```
RENUMBER 1058 9999 58   # Remove offset of 1000 from PDB data.
```

Renumbers all residues within the range specified by the first two arguments, so that the range will now start on the number given by the third
argument. This is useful in cases where e.g. a chimeric peptide segment is prepended to a receptor, and the receptor's native residue numbers are offset 
to compensate for the longer chain; after deleting the prepended segment, the remaining sequence can be restored to its original numbering.


# ROTATE
Example:
```
ROTATE @axis &angle
ROTATE @axis &angle %center_resno
ROTATE @axis &angle %start_resno %end_resno
ROTATE @axis &angle %center_resno %start_resno %end_resno
```

Rotates the entire protein about its center, using the first parameter as a _relative_ axis of rotation,
and the second parameter as the rotation angle in degrees.
In other words, the first param is added to the protein center to get the true axis of rotation.

If three arguments are given, the third argument defines a residue whose alpha carbon functions as the center of rotation.

If four arguments are given, the third and fourth define the start and end residue number of the segment to apply the rotation to.

Five arguments means axis, angle, center of rotation, start residue, and end residue.


# SAVE
Example:
```
SAVE "path/to/output.pdb"
SAVE "path/to/output.pdb" QUIT
```

Saves the protein in working memory to the specified file path (can be a string variable).
All strands that contain protein data will be included.

If `QUIT`, `EXIT`, or `END` follows the path parameter, script execution will terminate immediately after saving the protein.


# SCENV
Side Chain ENVironment.

Example:
```
SCENV %resno $result
SCENV 185 4 $result
```

Examines the local environment near the side chain of the indicated residue. Of all residues close enough to interact with the
target residue, if any side chain atom forms a non-covalent bond with any atom of the residue's side chain, then the nearby residue's
letter code and residue number will be added to the output. For example, if a protein has ASP200 and ARG254 forming a salt bridge,
then `SCENV 200 $result` will place the value `R254 ` into the variable `$result`.

If there are two arguments, then the first is the residue number and the second is a string variable in which to place the output
data. If there are three arguments, then the middle argument is a bitmask that determines which types of interatomic binding to look
for, using the following values:

1 van der Waals bonding;
2 hydrogen bonding;
4 ionic bonding and metal coordination;
8 pi stacking and polar-pi interactions.

For example to only search for hydrogen bonds and ionic/metallic bonds on a given side chain, add the numbers 2 + 4 to get the
bitmask value 6.

Note Asn and Gln will register as charged side chains for the purpose of the search, because each amide group has a slight 
zwitterionic charge.


# SEARCH
Example:
```
SEARCH 1 %SEQLEN "MAYDRYVAIC" %olfr_motif
SEARCH 1 %SEQLEN "MAYDRYVAIC" TH 5 %olfr_motif %similarity
```

Performs a fuzzy search on the specified residue range looking for the closest match to the indicated search term, and stores the starting resno 
of the result in the variable given as the output parameter (`%olfr_motif` in the examples). The fuzzy search matches amino acids based on their
similarity using the metrics of hydrophilicity, aromaticity, metal binding capability, and charge. By default the search finds the closest match
**no matter what**, so it will not by itself determine whether a motif exists in the protein, but only find where the search string most closely
matches the sequence.

The optional `TH` keyword indicates a threshold number of residues that must match the search string **exactly**. Note `X` counts as an exact match
for any residue. If no match is found with at least the threshold number of exact match residues, then the output variable will be set to zero.
This is useful for finding whether a motif exists in the specified region.

The optional second output variable receives the total similarity of the best match, or zero if no match was found. This is useful for determining 
which of two or more motifs a region has, for example it is possible to search a region for both `"XHXXCE"` and `"HYXXCE"` and then compare similarity
values and use an `IF` to take a different action depending which is the better match.


# STRAND
Example:
```
STRAND B
STRAND B C

```

Changes the current working strand to the strand ID specified. If two strand IDs are given, the protein in the first strand is moved to the second
strand and the working strand is set to the second ID.


# STRLEN
Example:
```
STRPOS $string %out_var
```

Gets the length of the input string and sets `%out_var` to that number.


# STRPOS
Example:
```
STRPOS $haystack $needle %out_var
```

Searches the `$haystack` string for the `$needle` and, if found, sets the `%out_var` to the 1-based position of the start of the searched string.
If the second string is not found within the first string, `%out_var` is set to zero.


# UPRIGHT
Example:
```
UPRIGHT
UPRIGHT A
UPRIGHT DEFGH
```

Turns the protein "upright", so that the extracellular domain is in the +Y direction and the cytoplasmic domain is in the -Y direction.
The entire protein will also be centered at [0, 0, 0].

Valid only for transmembrane proteins with at least one transmembrane helix. The protein must have TMR1 through TMR{n} defined, where n <= 7.
(Helices beyond 7 will be ignored.)
If the PDB does not contain suitable `REMARK 650 HELIX` records, then the TMRs must be manually defined with the `REGION` command.

If the protein has a TMR4, it will be rotated to the +Z direction from TMR1. This affords a good view of the binding pocket of 7-helix GPCRs, but also
maintains compatibility with MS4A proteins since it does not depend on the existence of a TMR5, TMR6, or TMR7.

`UPRIGHT` takes an optional parameter to indicate which strand(s) are to be moved as a unit with the uprighting of the working strand.
The default is for all strands in memory to be included.


# UNCHAIN
Example:
```
UNCHAIN A
```

Deletes the indicated strand (peptide chain). Note the current working strand cannot be deleted; an error will result.


