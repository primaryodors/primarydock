#!/bin/bash

RED="\033[31;1m"
GRN="\033[32m"
CYN="\033[36m"
NC="\033[0m"

cd "$(dirname "$0")"
cd ..

clear


./test/unit_test_express.sh

./test/mol_assem_tests.sh


REPORT="testdata/protein_test.approved.txt"
./test/protein_test AAAAAAAAAA | sed '/^#/d' > testdata/received/protein_test.received.txt
RESULT=$(diff --unified $REPORT testdata/received/protein_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Protein test succeeded.${NC}\n"
else
    printf "${RED}Protein test FAILED.${NC}\n"
    diff --color --unified $REPORT testdata/received/protein_test.received.txt
fi


REPORT="testdata/coplanar_test.approved.txt"
./bin/pepteditor test/planar_h.pepd | sed '/^#/d' > testdata/received/coplanar_test.received.txt
RESULT=$(diff --unified $REPORT testdata/received/coplanar_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Coplanar test succeeded.${NC}\n"
else
    printf "${RED}Coplanar test FAILED.${NC}\n"
    diff --color --unified $REPORT testdata/received/coplanar_test.received.txt
fi


REPORT="testdata/motif_test.approved.txt"
./bin/pepteditor test/motif_test.pepd | sed '/^#/d' >testdata/received/motif_test.received.txt
RESULT=$(diff --unified $REPORT testdata/received/motif_test.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}Motif test succeeded.${NC}\n"
else
    printf "${RED}Motif test FAILED.${NC}\n"
    diff --color --unified $REPORT testdata/received/motif_test.received.txt
fi



printf "${CYN}Running prediction tests and docking tests; these will take some time. Please wait.${NC}\n"


REPORT="testdata/OR51E2_propionate_pred.approved.txt"
php -f predict/method_directmdl.php prot=OR51E2 lig=propionic_acid | grep '[[]Predicted[]] => ' > testdata/received/OR51E2_propionate_pred.received.txt
RESULT=$(diff --unified $REPORT testdata/received/OR51E2_propionate_pred.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}OR51E2 propionate prediction test succeeded.${NC}\n"
else
    printf "${RED}OR51E2 propionate prediction test FAILED.${NC}\n"
    diff --color --unified $REPORT testdata/received/OR51E2_propionate_pred.received.txt
fi


REPORT="testdata/TAAR8_cadaverine_pred.approved.txt"
php -f predict/method_fygactive.php prot=TAAR8 lig=cadaverine | grep '[[]Predicted[]] => ' > testdata/received/TAAR8_cadaverine_pred.received.txt
RESULT=$(diff --unified $REPORT testdata/received/TAAR8_cadaverine_pred.received.txt)
if [ -z "$RESULT" ]; then
    # printf "${GRN}TAAR8 cadaverine prediction test succeeded.${NC}\n"
    ASP111=$( cat output/TAAR/TAAR8/TAAR8.cadaverine.active.dock | grep -m 1 "Asp111(3.32): " )
    ASP111="${ASP111/Asp111(3.32): /}"
    ASP111="${ASP111/[.][0-9]*/}"
    ASP201=$( cat output/TAAR/TAAR8/TAAR8.cadaverine.active.dock | grep -m 1 "Asp201(5.43): " )
    ASP201="${ASP201/Asp201(5.43): /}"
    ASP201="${ASP201/[.][0-9]*/}"
    if [[ $ASP111 -gt "-30"  ]] || [[ $ASP201 -gt "-30"  ]]; then
        printf "${RED}TAAR test FAILED: bad contacts.${NC}\n"
    else
        printf "${GRN}TAAR test succeeded.${NC}\n"
    fi
else
    printf "${RED}TAAR8 cadaverine prediction test FAILED to predict agonist.${NC}\n"
    diff --color --unified $REPORT testdata/received/OR51E2_propionate_pred.received.txt
fi
