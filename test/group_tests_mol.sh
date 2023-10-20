
MOLECULE="phenethyl_alcohol"
test/group_test_mol "sdf/$MOLECULE.sdf" | sed '/^#/d' > testdata/received/groups.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/groups.$MOLECULE.approved.txt testdata/received/groups.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Molecule group test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Molecule group test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/groups.$MOLECULE.approved.txt testdata/received/groups.$MOLECULE.received.txt
fi

MOLECULE="indole"
test/group_test_mol "sdf/$MOLECULE.sdf" | sed '/^#/d' > testdata/received/groups.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/groups.$MOLECULE.approved.txt testdata/received/groups.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Molecule group test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Molecule group test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/groups.$MOLECULE.approved.txt testdata/received/groups.$MOLECULE.received.txt
fi

MOLECULE="pyrazine"
test/group_test_mol "sdf/$MOLECULE.sdf" | sed '/^#/d' > testdata/received/groups.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/groups.$MOLECULE.approved.txt testdata/received/groups.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Molecule group test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Molecule group test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/groups.$MOLECULE.approved.txt testdata/received/groups.$MOLECULE.received.txt
fi

MOLECULE="cinnamaldehyde"
test/group_test_mol "sdf/$MOLECULE.sdf" | sed '/^#/d' > testdata/received/groups.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/groups.$MOLECULE.approved.txt testdata/received/groups.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Molecule group test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Molecule group test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/groups.$MOLECULE.approved.txt testdata/received/groups.$MOLECULE.received.txt
fi

MOLECULE="arabinose"
test/group_test_mol "sdf/l-arabinose.sdf" | sed '/^#/d' > testdata/received/groups.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/groups.$MOLECULE.approved.txt testdata/received/groups.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Molecule group test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Molecule group test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/groups.$MOLECULE.approved.txt testdata/received/groups.$MOLECULE.received.txt
fi

MOLECULE="adenosine"
test/group_test_mol "sdf/$MOLECULE.sdf" | sed '/^#/d' > testdata/received/groups.$MOLECULE.received.txt
RESULT=$(diff --unified testdata/groups.$MOLECULE.approved.txt testdata/received/groups.$MOLECULE.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Molecule group test succeeded for $MOLECULE.${NC}\n"
else
    printf "${RED}Molecule group test FAILED for $MOLECULE.${NC}\n"
    diff --color --unified testdata/groups.$MOLECULE.approved.txt testdata/received/groups.$MOLECULE.received.txt
fi
