#!/bin/bash

RED="\033[38;5;226m\033[48;5;196m"
GRN="\033[38;5;46m"
NC="\033[0m"

cd "$(dirname "$0")"
if [ ! -f "../testdata/received" ]; then mkdir -p "../testdata/received"; fi
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


echo "Running prediction tests and docking tests; these will take some time. Please wait."

REPORT="testdata/OR51E2_propionate_pred.approved.txt"
# php -f predict/method_icactive.php prot=OR51E2 lig=propionic_acid | tee >( grep '[[]Predicted[]] => ' > testdata/received/OR51E2_propionate_pred.received.txt)
php -f predict/method_icactive.php prot=OR51E2 lig=propionic_acid nosoft=1 | grep '[[]Predicted[]] => ' > testdata/received/OR51E2_propionate_pred.received.txt
RESULT=$(diff --unified $REPORT testdata/received/OR51E2_propionate_pred.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}OR51E2 propionate prediction test succeeded.${NC}\n"
else
    printf "${RED}OR51E2 propionate prediction test FAILED.${NC}\n"
    diff --color --unified $REPORT testdata/received/OR51E2_propionate_pred.received.txt
fi


REPORT="testdata/OR1A1_hexenol_pred.approved.txt"
php -f predict/method_icactive.php prot=OR1A1 lig=cis-3-hexen-1-ol nosoft=1 | grep '[[]Predicted[]] => ' > testdata/received/OR1A1_hexenol_pred.received.txt
RESULT=$(diff --unified $REPORT testdata/received/OR1A1_hexenol_pred.received.txt)
if [ -z "$RESULT" ]; then
    printf "${GRN}OR1A1 hexenol prediction test succeeded.${NC}\n"
else
    printf "${RED}OR1A1 hexenol prediction test FAILED.${NC}\n"
    diff --color --unified $REPORT testdata/received/OR1A1_hexenol_pred.received.txt
fi


# php -f predict/method_icactive.php prot=OR1A1 lig=d-limonene | tee >( grep '[[]Predicted[]] => ' > testdata/received/OR1A1_d-limonene_pred.received.txt)


bin/primarydock testdata/test_TAAR8.config --colorless --iter 50 --congress > testdata/received/TAAR8_CAD.txt
TAAR_RESULT=$?
if [ "$TAAR_RESULT" -eq "0" ]; then
    POSES=$( cat testdata/received/TAAR8_CAD.txt | grep "pose(s) found" )
    if [ -z "$POSES" ]; then
        printf "${RED}TAAR test FAILED: no poses.${NC}\n"
    else
        ASP111=$( cat testdata/received/TAAR8_CAD.txt | grep -m 1 "Asp111: " )
        ASP111="${ASP111/Asp111: /}"
        ASP111="${ASP111/[.][0-9]*/}"
        ASP201=$( cat testdata/received/TAAR8_CAD.txt | grep -m 1 "Asp201: " )
        ASP201="${ASP201/Asp201: /}"
        ASP201="${ASP201/[.][0-9]*/}"
        if [[ $ASP111 -gt "-15"  ]] || [[ $ASP201 -gt "-10"  ]]; then
            printf "${RED}TAAR test FAILED: bad contacts.${NC}\n"
        else
            printf "${GRN}TAAR test succeeded.${NC}\n"
        fi
    fi
else
    printf "${RED}TAAR test FAILED: return value.${NC}\n"
fi


# bin/primarydock testdata/test_1A1.config --colorless --congress > testdata/received/OR1A1_dLIMN.txt
# DLIMN_RESULT=$?
# if [ "$DLIMN_RESULT" -eq "0" ]; then
#     POSES=$( cat testdata/received/OR1A1_dLIMN.txt | grep "pose(s) found" | sed 's/[^0-9]//g' )
#     if [ "$POSES" -eq "0" ]; then
#         printf "${RED}d-limonene test FAILED: no poses.${NC}\n"
#     else
#         printf "${GRN}d-limonene test succeeded.${NC}\n"
#     fi
# else
#     printf "${RED}d-limonene test FAILED: return value.${NC}\n"
# fi

