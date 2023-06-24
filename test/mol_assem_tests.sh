
# Note since {AtomName} is a nonstandard feature, the code will always attempt to natively decode a string
# containing { rather than integrate with a third party app.
# The ! char suppresses the warning about native functionality often getting rings wrong.

allmols=("benzene" "toluene" "toluene1" "#cyclohexane" "#cyclopentane" "#cyclobutane" "#cyclopropane" "#cyclobutene" "#cyclobutadiene" "#cyclopropene" \
         "#calicene" "#2_6-xylenol" "#proline" "#total_substitution" "#o-methylfuranium_ion" \
         \
        )
allsmiles=("c{C1}1ccccc1!" "C{C1}c1ccccc1!" "c{C1}1ccccc1!C" "C{C1}1CCCCC1!" "C{C1}1CCCC1!" "C{C1}1CCC1!" "C{C1}1CC1!" "C{C1}1=CCC1!" "C{C1}1=CC=C1!" "C{C1}1=CC1!" \
         "C{C1}1(C=CC=C1!)=C2C=C2" "O{O1}c1c(C)cccc1!C" "N{N1}1[C@H](C(=O)O)CCC1!" "C{C1}c1c([Cl])c(CSC)c(O)c(C(=O)[O-])c1![NH3+]" "c{C1}1cc[o+](C)c1!" \
         \
        )

RED="\033[38;5;226m\033[48;5;196m"
GRN="\033[38;5;46m"
CYN="\033[38;5;68m"
NC="\033[0m"

for i in ${!allmols[@]}; do
    MOLECULE=${allmols[$i]}
    CHAR=${MOLECULE:0:1}

    if [ $CHAR == "#" ]; then
        printf "${CYN}Skipping $MOLECULE known to be failing.${NC}\n"
        continue
    fi

    SMILES=${allsmiles[$i]}
    REPORT="test/mol_assem_test.$MOLECULE.approved.txt"
    test/mol_assem_test "$SMILES" "test/received/$MOLECULE.received.sdf" | sed '/^#/d' > "test/received/mol_assem_test.$MOLECULE.received.txt"
    RESULT=$(diff --unified $REPORT "test/received/mol_assem_test.$MOLECULE.received.txt")
    if [ -z "$RESULT" ]; then
        printf "${GRN}Molecule assembly test succeeded for $MOLECULE.${NC}\n"
    else
        printf "${RED}Protein test FAILED for $MOLECULE.${NC}\n"
        diff --color --unified $REPORT "test/received/mol_assem_test.$MOLECULE.received.txt"
    fi
done


# test/mol_assem_test 'c{C1}1c2ccccc2ccc1!' test/naphthalene.received.sdf
# test/mol_assem_test 'c{C1}1ccccc1!Oc2ccccc2' test/diphenyl_ether.received.sdf
# test/mol_assem_test 'c{C1}1ccccc1!N=Nc2ccccc2' test/azobenzene.received.sdf
# test/mol_assem_test 'n{N1}1c[nH]cc1!' test/imidazole.received.sdf
# test/mol_assem_test 'C{C1}C/C=C\CCO' test/leaf_alcohol.received.sdf
# test/mol_assem_test 'C{C1}([C@@H]1[C@H]([C@@H]([C@H](C(O1!)O)O)O)O)O' test/glucose.received.sdf
# test/mol_assem_test 'C{C1}([C@@H]1[C@H]([C@@H]([C@H](C(S1!)S)O)S)S)O' test/tetrathioglucose.received.sdf
# test/mol_assem_test 'C{C1}([C@@H]1[C@H]([C@@H]([C@H](C(S1!)S)S)S)S)S' test/hexathioglucose.received.sdf

