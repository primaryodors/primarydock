
# PrimaryDock Interpreted Script Documentation

To run a PrimaryDock interpreted script, after building PrimaryDock, please use the following command:

```
bin/interpreter path/to/script.pdis
```

# Variables

With few exceptions, any variable can stand in for any parameter of any command.

All integers begin with `%`, e.g. `%start_resno`.
All floats begin with `&`, e.g. `&transpose_amt`.
All Cartesian 3D locations begin with `@`, e.g. `@pocket_center`.
All strings begin with `$`, e.g. `$name`.

Cartesians have members .x, .y, and .z that behave as floats.

Casting a float to an integer rounds the value (e.g. 0.4 rounds to zero but 0.5 rounds to 1).

Casting a float to a Cartesian normally sets the .x member to the float value, leaving .y=0 and .z=0. But a command like `LET @foo = @bar.y` will result in only the value of `@foo.y` being nonzero.

Casting an integer to a Cartesian obtains the location of the CA atom for that residue number, if it exists.

Casting a Cartesian back to float or integer obtains the magnitude of the Cartesian, equal to sqrt(x^2 + y^2 + z^2).

The command line arguments are made available to the script as $arg1, $arg2, etc.

The following "magic variables" are supplied upon loading a protein:
- `$PDB` the path and name of the source PDB file.
- `$PROTEIN` the name of the protein, derived from a `REMARK 6` record if present.
- `%SEQLEN` the length of the protein sequence.
- `$SEQUENCE` the sequence of the protein in standard one-letter amino acid code.
- Region start and end residue numbers from any `REMARK 650 HELIX` records, e.g. `%TMR3.s` and `%TMR3.e` for the start and end resnos of TMR3.

The following commands are supported:


# ALIGN
Example:
```
ALIGN 174 182 4 @location1 @location2
```

Repositions a piece of the protein so that the residues near the start resno occur at one location, and the residues near the end resno occur in the
direction of the other location.

The first two numbers are the starting and ending residue numbers. They define the range of residues that will be moved.

The next number is the averaging length from each end of the region. In the example this number is 4, so the locations of the CA atoms of four residues
at each end will be averaged together, in this case 174-177 and 179-182, to determine the alignment. The piece will be moved so that the average of CA
atom locations of residues 174-177 will now equal @location1, and the average of CA locations of residues 179-182 will fall along a straight line pointing
from @location1 to @location2.


# BEND
Example:
```
BEND 206 218 "N-CA" 15
```

Effects a partial bend of the protein backbone at the start residue (first parameter) by rotating in the bond direction (third parameter) by the rotation
angle (fourth parameter, degrees), moving all subsequent residues up to the end residue (second parameter). The result is a bend in the chain corresponding
to either a phi or psi rotation of the start residue.

Note that the end residue will be disconnected from its neighbor and will require to be reconnected.

Valid bond directions are:
N-CA		Keep the residue's N terminus stationary and rotate its C terminus about the N-CA (phi) bond;
CA-C		Keep the residue's N terminus stationary and rotate its C terminus about the CA-C (psi) bond;
CA-N		Keep the residue's C terminus stationary and rotate its N terminus about the CA-N (phi) bond;
C-CA		Keep the residue's C terminus stationary and rotate its N terminus about the C-CA (psi) bond;

Note if bending in the N-CA or CA-C direction, i.e. towards higher numbered residues, then the end residue must be greater than the start residue.
If bending in the C-CA or CA-N direction, the start residue must be greater than the end residue.

See here for an introduction to phi and psi bonds:
https://proteopedia.org/wiki/index.php/Tutorial:Ramachandran_principle_and_phi_psi_angles


# BENERG
Example:
```
BEND 251 111 &binding
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
```

Ceases script execution.


# GEN
```
GEN "MAYDRYVAIC"
```

Creates a peptide using the specified sequence.


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


# IF
Examples:
```
IF %var >= 10 LET %var = 0
IF &var1 > 0 AND &var1 < &var2 SAVE $name QUIT

IF $match = "HFFCE" ECHO $message
ELSE ECHO "Nope."

LET %iter = 1
loop:
ECHO %iter
LET %iter ++
IF %iter <= 10 GOTO loop
```

The `IF` command evaluates a conditional expression and, if true, executes another command. Currently, only simple A = B type expressions are
supported, i.e. straightforward comparison of two values, where either value can be a single variable or a single constant; such expressions as
A = B + 1 do not yet work.

`AND` and `OR` are supported, however parentheses are not, and while an unlimited number of `AND`s and/or `OR`s can be chained together, their
evaluation will be entirely sequential and the chain will evaluate false as soon as any chain of `OR`ed expressions all evaluate false. Note that
any `ELSE` after such a statement will execute no matter which part of the chain evaluated to false.

The operators available for `IF` are `=` `!=` `>` `<` `>=` `<=`. Note that `=` means comparison, not assignment.

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
LET @res174loc = 174
LET &j += @res174loc.y
LET @loc3 = @loc1 + @loc2
LET $message = "TMR4 ends on residue " + %TMR4.e + " and TMR5 starts on residue " + %TMR5.s + "."
```

Assigns a value to a new or existing variable.

For assigning integer and float variables, the following operators are allowed: `=` `+=` `-=` `*=` `/=`.

Integer variables are also allowed `++` `--`, which do not take an r-value (a variable or literal value to apply to the variable assignment).

For algebraic assignments, e.g. `LET &k = %i + &j`, the following operators are supported: `+` `-` `*` `/` `^`, the last being an exponent operator.
Note parentheses are not currently supported. Also, the data type of the variable being assigned determines how all r-value variables are interpreted,
so even though `%i` is an integer in the above statement, its value is cast to a float before `&j` is added to it.

For strings, the operators `=` `+=` and `+` are allowed, with the plus sign indicating concatenation.
Substrings are indicated with the word FROM and (optinally) FOR, e.g. `LET $end = $string FROM 10` or `LET $range = $SEQUENCE FROM %start FOR %length`.


# LOAD
Example:
```
LOAD "path/to/protein.pdb"
```

Loads the specified protein into working memory. Note only one protein is allowed in memory at a given time.


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


# SAVE
Example:
```
SAVE "path/to/output.pdb"
SAVE "path/to/output.pdb" QUIT
```

Saves the protein in working memory to the specified file path (can be a string variable).

If `QUIT`, `EXIT`, or `END` follows the path parameter, script execution will terminate immediately after saving the protein.


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











