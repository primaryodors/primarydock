#!/bin/bash

RED="\033[38;5;226m\033[48;5;196m"
GRN="\033[38;5;46m"
NC="\033[0m"

cd "$(dirname "$0")"
if [ ! -f "received" ]; then mkdir -p "received"; fi
cd ..

clear


./test/unit_test_express.sh

./test/mol_assem_tests.sh


REPORT="test/protein_test.approved.txt"
./test/protein_test AAAAAAAAAA | sed '/^#/d' > test/received/protein_test.received.txt
RESULT=$(diff --unified $REPORT test/received/protein_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Protein test succeeded.${NC}\n"
else
    printf "${RED}Protein test FAILED.${NC}\n"
    diff --color --unified $REPORT test/received/protein_test.received.txt
fi


REPORT="test/coplanar_test.approved.txt"
./bin/pepteditor test/planar_h.pepd | sed '/^#/d' > test/received/coplanar_test.received.txt
RESULT=$(diff --unified $REPORT test/received/coplanar_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Coplanar test succeeded.${NC}\n"
else
    printf "${RED}Coplanar test FAILED.${NC}\n"
    diff --color --unified $REPORT test/received/coplanar_test.received.txt
fi


REPORT="test/delwork_test.approved.txt"
./bin/pepteditor test/delwork.pepd | sed '/^#/d' > test/received/delwork_test.received.txt
RESULT=$(diff --unified $REPORT test/received/delwork_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Delete working strand test succeeded.${NC}\n"
else
    printf "${RED}Delete working strand test FAILED.${NC}\n"
    diff --color --unified $REPORT test/received/delwork_test.received.txt
fi


REPORT="test/motif_test.approved.txt"
./bin/pepteditor test/motif_test.pepd | sed '/^#/d' >test/received/motif_test.received.txt
RESULT=$(diff --unified $REPORT test/received/motif_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Motif test succeeded.${NC}\n"
else
    printf "${RED}Motif test FAILED.${NC}\n"
    diff --color --unified $REPORT test/received/motif_test.received.txt
fi


bin/primarydock test/testTAAR8.config --colorless --pose 5 --iter 50 --congress > test/received/TAAR8_CAD.txt
TAAR_RESULT=$?
if [ "$TAAR_RESULT" -eq "0" ]; then
    POSES=$( cat test/received/TAAR8_CAD.txt | grep "pose(s) found" )
    if [ -z "$POSES" ]; then
        printf "${RED}TAAR test FAILED: no poses.${NC}\n"
    else
        ASP111=$( cat test/received/TAAR8_CAD.txt | grep -m 1 "Asp111: " )
        ASP111="${ASP111/Asp111: /}"
        ASP111="${ASP111/[.][0-9]*/}"
        ASP201=$( cat test/received/TAAR8_CAD.txt | grep -m 1 "Asp201: " )
        ASP201="${ASP201/Asp201: /}"
        ASP201="${ASP201/[.][0-9]*/}"
        TAARPDB=$( cat output/test_TAAR8_cadaverine.dock | grep -m 1 "HETATM 9000  N1  LIG " )
        if [[ $ASP111 -gt "-25"  ]] || [[ $ASP201 -gt "-25"  ]]; then
            printf "${RED}TAAR test FAILED: bad contacts.${NC}\n"
        else
            if [ -z "$TAARPDB" ]; then
                printf "${RED}TAAR test FAILED: bad PDB data.${NC}\n"
            else
                printf "${GRN}TAAR test succeeded.${NC}\n"
            fi
        fi
    fi
else
    printf "${RED}TAAR test FAILED: return value.${NC}\n"
fi


bin/primarydock test/test1A1.config --colorless --pose 5 --iter 50 --congress > test/received/OR1A1_dLIMN.txt
DLIMN_RESULT=$?
if [ "$DLIMN_RESULT" -eq "0" ]; then
    POSES=$( cat test/received/OR1A1_dLIMN.txt | grep "pose(s) found" )
    if [ -z "$POSES" ]; then
        printf "${RED}d-limonene test FAILED: no poses.${NC}\n"
    else
        printf "${GRN}d-limonene test succeeded.${NC}\n"
    fi
else
    printf "${RED}d-limonene test FAILED.${NC}\n"
fi

