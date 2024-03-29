#!/bin/bash

RED="\033[31;1m"
GRN="\033[32m"
NC="\033[0m"

cd "$(dirname "$0")"
cd ..


make


REPORT="testdata/molecule_test1.approved.txt"
./test/molecule_test 'CCO' 'CCO' | tee temp | sed '/^#/d' >testdata/received/molecule_test1.received.txt
RESULT=$(diff --unified $REPORT testdata/received/molecule_test1.received.txt)
CODE_RESULT=$?
if [ $CODE_RESULT -eq 0 ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Molecule test FAILED.${NC}\n"
    diff --color --unified $REPORT testdata/received/molecule_test1.received.txt
fi


REPORT="testdata/point_test.approved.txt"
./test/point_test | sed '/^#/d' >testdata/received/point_test.received.txt
RESULT=$(diff --unified $REPORT testdata/received/point_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Point test FAILED.${NC}\n"
    diff --color --unified $REPORT testdata/received/point_test.received.txt
fi


REPORT="testdata/atom_test.approved.txt"
./test/atom_test H | sed '/^#/d' >testdata/received/atom_test.received.txt
RESULT=$(diff --unified $REPORT testdata/received/atom_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Atom test FAILED.${NC}\n"
    diff --color --unified $REPORT testdata/received/atom_test.received.txt
fi


REPORT="testdata/aniso_test.approved.txt"
./test/aniso_test --asciiart | sed '/^#/d' >testdata/received/aniso_test.received.txt
RESULT=$(diff --unified $REPORT testdata/received/aniso_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Anisotropy test FAILED.${NC}\n"
    diff --color --unified $REPORT testdata/received/aniso_test.received.txt
fi


REPORT="testdata/amino_test.approved.txt"
./test/amino_test | sed '/^#/d' >testdata/received/amino_test.received.txt
RESULT=$(diff --unified $REPORT testdata/received/amino_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Amino test FAILED.${NC}\n"
    diff --color --unified $REPORT testdata/received/amino_test.received.txt
fi


REPORT="testdata/ifnot_test.approved.txt"
./bin/pepteditor test/ifnot.pepd | sed '/^#/d' >testdata/received/ifnot_test.received.txt
RESULT=$(diff --unified $REPORT testdata/received/ifnot_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}IF NOT test FAILED.${NC}\n"
    diff --color --unified $REPORT testdata/received/ifnot_test.received.txt
fi


REPORT="testdata/measure_test.approved.txt"
./bin/pepteditor test/measure.pepd | sed '/^#/d' >testdata/received/measure_test.received.txt
RESULT=$(diff --unified $REPORT testdata/received/measure_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}MEASURE test FAILED.${NC}\n"
    diff --color --unified $REPORT testdata/received/measure_test.received.txt
fi


REPORT="testdata/bwmagic.approved.txt"
bin/pepteditor test/bwmagic.pepd | sed '/^#/d' >testdata/received/bwmagic.received.txt
RESULT=$(diff --unified $REPORT testdata/received/bwmagic.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}BW magic variables test FAILED.${NC}\n"
    diff --color --unified $REPORT testdata/received/bwmagic.received.txt
fi


MOLECULE="benzaldehyde"
test/bond_rotation_test "c1ccccc1C=O" | sed '/^#/d' > testdata/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt
fi


MOLECULE="acetophenone"
test/bond_rotation_test "c1ccccc1C(=O)C" | sed '/^#/d' > testdata/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt
fi


MOLECULE="cinnamaldehyde"
test/bond_rotation_test "c1ccccc1C=CC=O" | sed '/^#/d' > testdata/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt
fi


MOLECULE="leaf_alcohol"
test/bond_rotation_test "CC\\C=C/CCO" | sed '/^#/d' > testdata/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt
fi


MOLECULE="azobenzene"
test/bond_rotation_test "c1ccccc1N=Nc2ccccc2" | sed '/^#/d' > testdata/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt
fi


MOLECULE="phenol"
test/bond_rotation_test "c1ccccc1O" | sed '/^#/d' > testdata/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt
fi


MOLECULE="ethyl_pyrazine"
test/bond_rotation_test "c1nccnc1CC" | sed '/^#/d' > testdata/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt
fi


MOLECULE="hydrogen_peroxide"
test/bond_rotation_test "OO" | sed '/^#/d' > testdata/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt
fi


MOLECULE="2-butyne"
test/bond_rotation_test "CC#CC" | sed '/^#/d' > testdata/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt
fi


MOLECULE="ethyl_acetate"
test/bond_rotation_test "CCOC(=O)C" | sed '/^#/d' > testdata/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt
fi


MOLECULE="eucalyptol"
test/bond_rotation_test "sdf/eucalyptol.sdf" | sed '/^#/d' > testdata/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt
fi


MOLECULE="dipeptide_GG"
test/bond_rotation_test "NCC(=O)NCC(=O)O" | sed '/^#/d' > testdata/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/brot.$MOLECULE.approved.txt testdata/received/brot.$MOLECULE.received.txt
fi


# group_tests_mol.sh

test/group_test_res pdbs/TAAR/TAAR8.upright.pdb sdf/cadaverine.sdf 3.32 3.37 5.42 6.48 7.43 | grep "atom_group\[ N" | sed '/^#/d' > testdata/received/group_test_res.received.txt
RESULT=$(diff --unified testdata/group_test_res.taar8.approved.txt testdata/received/group_test_res.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Residue group test FAILED for TAAR8.${NC}\n"
    diff --color --unified testdata/group_test_res.taar8.approved.txt testdata/received/group_test_res.received.txt
fi

bin/pepteditor testdata/8f76.pepd 2>&1 > /dev/null

test/group_test_res pdbs/OR51/OR51E2.active.pdb sdf/propionic_acid.sdf 4.57 6.59 | sed '/^#/d' > testdata/received/group_test_res.received.txt
RESULT=$(diff --unified testdata/group_test_res.51e2.approved.txt testdata/received/group_test_res.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Residue group test FAILED for OR51E2.${NC}\n"
    diff --color --unified testdata/group_test_res.51e2.approved.txt testdata/received/group_test_res.received.txt
fi

test/ring_test sdf/histidine.sdf | sed '/^#/d' > testdata/received/ring.histidine.received.txt
RESULT=$(diff --unified testdata/ring.histidine.approved.txt testdata/received/ring.histidine.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Ring test FAILED for histidine.${NC}\n"
    diff --color --unified testdata/ring.histidine.approved.txt testdata/received/ring.histidine.received.txt
fi

test/ring_test "O=C1COc2c(OC1)cc(cc2)C" | sed '/^#/d' > testdata/received/ring.calone.received.txt
RESULT=$(diff --unified testdata/ring.calone.approved.txt testdata/received/ring.calone.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Ring test FAILED for calone.${NC}\n"
    diff --color --unified testdata/ring.calone.approved.txt testdata/received/ring.calone.received.txt
fi

test/ring_test sdf/androstenone.sdf | sed '/^#/d' > testdata/received/ring.androstenone.received.txt
RESULT=$(diff --unified testdata/ring.androstenone.approved.txt testdata/received/ring.androstenone.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Ring test FAILED for androstenone.${NC}\n"
    diff --color --unified testdata/ring.androstenone.approved.txt testdata/received/ring.androstenone.received.txt
fi

test/moiety_test sdf/cinnamaldehyde.sdf "HC(C)=O" | sed '/^#/d' > testdata/received/moiety.aldehyde.received.txt
RESULT=$(diff --unified testdata/moiety.aldehyde.approved.txt testdata/received/moiety.aldehyde.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Moiety test (aldehyde) FAILED.${NC}\n"
    diff --color --unified testdata/moiety.aldehyde.approved.txt testdata/received/moiety.aldehyde.received.txt
fi

test/moiety_test sdf/phenethyl_alcohol.sdf "c1ccccc1" | sed '/^#/d' > testdata/received/moiety.bzring.received.txt
RESULT=$(diff --unified testdata/moiety.bzring.approved.txt testdata/received/moiety.bzring.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Moiety test (benzene ring) FAILED.${NC}\n"
    diff --color --unified testdata/moiety.bzring.approved.txt testdata/received/moiety.bzring.received.txt
fi

test/moiety_test sdf/cis-3-hexen-1-ol.sdf "cc" | sed '/^#/d' > testdata/received/moiety.pibond.received.txt
RESULT=$(diff --unified testdata/moiety.pibond.approved.txt testdata/received/moiety.pibond.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Moiety test (pi bond) FAILED.${NC}\n"
    diff --color --unified testdata/moiety.pibond.approved.txt testdata/received/moiety.pibond.received.txt
fi

test/moiety_test sdf/linalool.sdf "COH" | sed '/^#/d' > testdata/received/moiety.alcohol.received.txt
RESULT=$(diff --unified testdata/moiety.alcohol.approved.txt testdata/received/moiety.alcohol.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Moiety test (alcohol) FAILED.${NC}\n"
    diff --color --unified testdata/moiety.alcohol.approved.txt testdata/received/moiety.alcohol.received.txt
fi

test/moiety_test "[Cl]c1ccc([Cl])cc1" "C[Cl]" | sed '/^#/d' > testdata/received/moiety.halide.received.txt
RESULT=$(diff --unified testdata/moiety.halide.approved.txt testdata/received/moiety.halide.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Moiety test (halide) FAILED.${NC}\n"
    diff --color --unified testdata/moiety.halide.approved.txt testdata/received/moiety.halide.received.txt
fi

test/moiety_test sdf/isovalerate.sdf "oc[o-]" | sed '/^#/d' > testdata/received/moiety.acid.received.txt
RESULT=$(diff --unified testdata/moiety.acid.approved.txt testdata/received/moiety.acid.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Moiety test (acid) FAILED.${NC}\n"
    diff --color --unified testdata/moiety.acid.approved.txt testdata/received/moiety.acid.received.txt
fi

test/moiety_test sdf/cadaverine.sdf "N[H+]" | sed '/^#/d' > testdata/received/moiety.amine.received.txt
RESULT=$(diff --unified testdata/moiety.amine.approved.txt testdata/received/moiety.amine.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}\u2588${NC}"
else
    printf "\n${RED}Moiety test (amine) FAILED.${NC}\n"
    diff --color --unified testdata/moiety.amine.approved.txt testdata/received/moiety.amine.received.txt
fi

printf "\n\n"
