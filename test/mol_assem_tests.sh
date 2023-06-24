
# Note since {AtomName} is a nonstandard feature, the code will always attempt to natively decode a string
# containing { rather than integrate with a third party app.
# The ! char suppresses the warning about native functionality often getting rings wrong.

allmols=("benzene" "toluene")
allsmiles=("c{C1}1ccccc1!" "C{C1}c1ccccc1!")

for i in ${!allmols[@]}; do
    MOLECULE=${allmols[$i]}
    SMILES=${allsmiles[$i]}
    REPORT="test/mol_assem_test.$MOLECULE.approved.txt"
    test/mol_assem_test "$SMILES" test/benzene.received.sdf | sed '/^#/d' > "test/mol_assem_test.$MOLECULE.received.txt"
    RESULT=$(diff --unified $REPORT "test/mol_assem_test.$MOLECULE.received.txt")
    if [ -z "$RESULT" ]; then
        echo "Molecule assembly test succeeded for $MOLECULE."
    else
        echo "Protein test FAILED for $MOLECULE."
        diff --color --unified $REPORT "test/mol_assem_test.$MOLECULE.received.txt"
    fi
done

# test/mol_assem_test 'C{C1}c1ccccc1!' test/toluene.received.sdf
# test/mol_assem_test 'c{C1}1ccccc1!C' test/toluene1.received.sdf
# test/mol_assem_test 'C{C1}1CCCCC1!' test/cyclohexane.received.sdf
# test/mol_assem_test 'C{C1}1CCCC1!' test/cyclopentane.received.sdf
# test/mol_assem_test 'C{C1}1CCC1!' test/cyclobutane.received.sdf
# test/mol_assem_test 'C{C1}1CC1!' test/cyclopropane.received.sdf
# test/mol_assem_test 'C{C1}1=CCC1!' test/cyclobutene.received.sdf
# test/mol_assem_test 'C{C1}1=CC=C1!' test/cyclobutadiene.received.sdf
# test/mol_assem_test 'C{C1}1=CC1!' test/cyclopropene.received.sdf
# test/mol_assem_test 'C{C1}1(C=CC=C1!)=C2C=C2' test/calicene.received.sdf
# test/mol_assem_test 'O{O1}c1c(C)cccc1!C' test/2_6-xylenol.received.sdf
# test/mol_assem_test 'N{N1}1[C@H](C(=O)O)CCC1!' test/proline.received.sdf
# test/mol_assem_test 'C{C1}c1c([Cl])c(CSC)c(O)c(C(=O)[O-])c1![NH3+]' test/total_substitution.received.sdf
# test/mol_assem_test 'c{C1}1cc[o+](C)c1!' test/o-methylfuranium_ion.received.sdf
# test/mol_assem_test 'c{C1}1c2ccccc2ccc1!' test/naphthalene.received.sdf
# test/mol_assem_test 'c{C1}1ccccc1!Oc2ccccc2' test/diphenyl_ether.received.sdf
# test/mol_assem_test 'c{C1}1ccccc1!N=Nc2ccccc2' test/azobenzene.received.sdf
# test/mol_assem_test 'n{N1}1c[nH]cc1!' test/imidazole.received.sdf
# test/mol_assem_test 'C{C1}C/C=C\CCO' test/leaf_alcohol.received.sdf
# test/mol_assem_test 'C{C1}([C@@H]1[C@H]([C@@H]([C@H](C(O1!)O)O)O)O)O' test/glucose.received.sdf
# test/mol_assem_test 'C{C1}([C@@H]1[C@H]([C@@H]([C@H](C(S1!)S)O)S)S)O' test/tetrathioglucose.received.sdf
# test/mol_assem_test 'C{C1}([C@@H]1[C@H]([C@@H]([C@H](C(S1!)S)S)S)S)S' test/hexathioglucose.received.sdf

