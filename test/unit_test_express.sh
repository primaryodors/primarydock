#!/bin/bash

RED="\033[38;5;226m\033[48;5;196m"
GRN="\033[38;5;46m"
NC="\033[0m"

cd "$(dirname "$0")"
if [ ! -f "received" ]; then mkdir -p "received"; fi
cd ..

clear

make code | sed '/^#/d' > test/received/make_code.txt
CODE_RESULT=$?
if [ $CODE_RESULT -eq 0 ]; then
    printf "${GRN}Make code succeeded.${NC}\n"
else
    printf "${RED}Make code FAILED.${NC}\n"
    echo $CODE_RESULT
fi


REPORT="test/point_test.approved.txt"
./test/point_test | sed '/^#/d' >test/received/point_test.received.txt
RESULT=$(diff --unified $REPORT test/received/point_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Point test succeeded.${NC}\n"
else
    printf "${RED}Point test FAILED.${NC}\n"
    diff --color --unified $REPORT test/received/point_test.received.txt
fi


REPORT="test/atom_test.approved.txt"
./test/atom_test H | sed '/^#/d' >test/received/atom_test.received.txt
RESULT=$(diff --unified $REPORT test/received/atom_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Atom test succeeded.${NC}\n"
else
    printf "${RED}Atom test FAILED.${NC}\n"
    diff --color --unified $REPORT test/received/atom_test.received.txt
fi


REPORT="test/aniso_test.approved.txt"
./test/aniso_test --asciiart | sed '/^#/d' >test/received/aniso_test.received.txt
RESULT=$(diff --unified $REPORT test/received/aniso_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Anisotropy test succeeded.${NC}\n"
else
    printf "${RED}Anisotropy test FAILED.${NC}\n"
    diff --color --unified $REPORT test/received/aniso_test.received.txt
fi


REPORT="test/amino_test.approved.txt"
./test/amino_test | sed '/^#/d' >test/received/amino_test.received.txt
RESULT=$(diff --unified $REPORT test/received/amino_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Amino test succeeded.${NC}\n"
else
    printf "${RED}Amino test FAILED.${NC}\n"
    diff --color --unified $REPORT test/received/amino_test.received.txt
fi


REPORT="test/ifnot_test.approved.txt"
./bin/pepteditor test/ifnot.pepd | sed '/^#/d' >test/received/ifnot_test.received.txt
RESULT=$(diff --unified $REPORT test/received/ifnot_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}IF NOT test succeeded.${NC}\n"
else
    printf "${RED}IF NOT test FAILED.${NC}\n"
    diff --color --unified $REPORT test/received/ifnot_test.received.txt
fi

REPORT="test/bwmagic.approved.txt"
bin/pepteditor test/bwmagic.pepd | sed '/^#/d' >test/received/bwmagic.received.txt
RESULT=$(diff --unified $REPORT test/received/bwmagic.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}BW magic variables test succeeded.${NC}\n"
else
    printf "${RED}BW magic variables test FAILED.${NC}\n"
    diff --color --unified $REPORT test/received/bwmagic.received.txt
fi


MOLECULE="benzaldehyde"
test/bond_rotation_test "c1ccccc1C=O" | sed '/^#/d' > test/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Rotatable bond test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt
fi


MOLECULE="acetophenone"
test/bond_rotation_test "c1ccccc1C(=O)C" | sed '/^#/d' > test/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Rotatable bond test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt
fi


MOLECULE="cinnamaldehyde"
test/bond_rotation_test "c1ccccc1C=CC=O" | sed '/^#/d' > test/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Rotatable bond test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt
fi


MOLECULE="leaf_alcohol"
test/bond_rotation_test "CC\\C=C/CCO" | sed '/^#/d' > test/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Rotatable bond test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt
fi


MOLECULE="azobenzene"
test/bond_rotation_test "c1ccccc1N=Nc2ccccc2" | sed '/^#/d' > test/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Rotatable bond test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt
fi


MOLECULE="phenol"
test/bond_rotation_test "c1ccccc1O" | sed '/^#/d' > test/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Rotatable bond test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt
fi


MOLECULE="ethyl_pyrazine"
test/bond_rotation_test "c1nccnc1CC" | sed '/^#/d' > test/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Rotatable bond test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt
fi


MOLECULE="hydrogen_peroxide"
test/bond_rotation_test "OO" | sed '/^#/d' > test/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Rotatable bond test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt
fi


MOLECULE="2-butyne"
test/bond_rotation_test "CC#CC" | sed '/^#/d' > test/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Rotatable bond test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt
fi


MOLECULE="ethyl_acetate"
test/bond_rotation_test "CCOC(=O)C" | sed '/^#/d' > test/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Rotatable bond test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt
fi


MOLECULE="dipeptide_GG"
test/bond_rotation_test "NCC(=O)NCC(=O)O" | sed '/^#/d' > test/received/brot.$MOLECULE.received.txt
RESULT=$(diff --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Rotatable bond test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Rotatable bond test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/brot.$MOLECULE.approved.txt test/received/brot.$MOLECULE.received.txt
fi


MOLECULE="phenethyl_alcohol"
test/group_test_mol "c1ccccc1CCO" | sed '/^#/d' > test/received/groups.$MOLECULE.received.txt
RESULT=$(diff --unified test/groups.$MOLECULE.approved.txt test/received/groups.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Molecule group test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Molecule group test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/groups.$MOLECULE.approved.txt test/received/groups.$MOLECULE.received.txt
fi

MOLECULE="indole"
test/group_test_mol "C12=C(C=CN2)C=CC=C1" | sed '/^#/d' > test/received/groups.$MOLECULE.received.txt
RESULT=$(diff --unified test/groups.$MOLECULE.approved.txt test/received/groups.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Molecule group test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Molecule group test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/groups.$MOLECULE.approved.txt test/received/groups.$MOLECULE.received.txt
fi

MOLECULE="pyrazine"
test/group_test_mol "n1ccncc1" | sed '/^#/d' > test/received/groups.$MOLECULE.received.txt
RESULT=$(diff --unified test/groups.$MOLECULE.approved.txt test/received/groups.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Molecule group test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Molecule group test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/groups.$MOLECULE.approved.txt test/received/groups.$MOLECULE.received.txt
fi

MOLECULE="cinnamaldehyde"
test/group_test_mol "c1ccccc1C=CC=O" | sed '/^#/d' > test/received/groups.$MOLECULE.received.txt
RESULT=$(diff --unified test/groups.$MOLECULE.approved.txt test/received/groups.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Molecule group test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Molecule group test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/groups.$MOLECULE.approved.txt test/received/groups.$MOLECULE.received.txt
fi

MOLECULE="arabinose"
test/group_test_mol "C1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O" | sed '/^#/d' > test/received/groups.$MOLECULE.received.txt
RESULT=$(diff --unified test/groups.$MOLECULE.approved.txt test/received/groups.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Molecule group test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Molecule group test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/groups.$MOLECULE.approved.txt test/received/groups.$MOLECULE.received.txt
fi

MOLECULE="adenosine"
test/group_test_mol "n2c1c(ncnc1n(c2)[C@@H]3O[C@@H]([C@@H](O)[C@H]3O)CO)N" | sed '/^#/d' > test/received/groups.$MOLECULE.received.txt
RESULT=$(diff --unified test/groups.$MOLECULE.approved.txt test/received/groups.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Molecule group test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Molecule group test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified test/groups.$MOLECULE.approved.txt test/received/groups.$MOLECULE.received.txt
fi

test/group_test_res | sed '/^#/d' > test/received/group_test_res.received.txt
RESULT=$(diff --unified test/group_test_res.approved.txt test/received/group_test_res.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Residue group test succeeded.${NC}\n"
else
    printf "${RED}Residue group test FAILED.${NC}\n"
    diff --color --unified test/group_test_res.approved.txt test/received/group_test_res.received.txt
fi
