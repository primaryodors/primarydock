#!/bin/bash
cd "$(dirname "$0")"
if [ ! -f "received" ]; then mkdir -p "received"; fi
cd ..

clear

make code | sed '/^#/d' > test/received/make_code.txt
CODE_RESULT=$?
if [ $CODE_RESULT -eq 0 ]; then
    echo "Make code succeeded."
else
    echo "Make code FAILED."
    echo $CODE_RESULT
fi


REPORT="test/point_test.approved.txt"
./test/point_test | sed '/^#/d' >test/received/point_test.received.txt
RESULT=$(diff --unified $REPORT test/received/point_test.received.txt)
if [ -z "$RESULT" ]; then
    echo "Point test succeeded."
else
    echo "Point test FAILED."
    diff --color --unified $REPORT test/received/point_test.received.txt
fi


REPORT="test/atom_test.approved.txt"
./test/atom_test H | sed '/^#/d' >test/received/atom_test.received.txt
RESULT=$(diff --unified $REPORT test/received/atom_test.received.txt)
if [ -z "$RESULT" ]; then
    echo "Atom test succeeded."
else
    echo "Atom test FAILED."
    diff --color --unified $REPORT test/received/atom_test.received.txt
fi


REPORT="test/aniso_test.approved.txt"
./test/aniso_test --asciiart | sed '/^#/d' >test/received/aniso_test.received.txt
RESULT=$(diff --unified $REPORT test/received/aniso_test.received.txt)
if [ -z "$RESULT" ]; then
    echo "Anisotropy test succeeded."
else
    echo "Anisotropy test FAILED."
    diff --color --unified $REPORT test/received/aniso_test.received.txt
fi


./test/mol_assem_tests.sh


REPORT="test/amino_test.approved.txt"
./test/amino_test | sed '/^#/d' >test/received/amino_test.received.txt
RESULT=$(diff --unified $REPORT test/received/amino_test.received.txt)
if [ -z "$RESULT" ]; then
    echo "Amino test succeeded."
else
    echo "Amino test FAILED."
    diff --color --unified $REPORT test/received/amino_test.received.txt
fi


REPORT="test/protein_test.approved.txt"
./test/protein_test AAAAAAAAAA | sed '/^#/d' > test/received/protein_test.received.txt
RESULT=$(diff --unified $REPORT test/received/protein_test.received.txt)
if [ -z "$RESULT" ]; then
    echo "Protein test succeeded."
else
    echo "Protein test FAILED."
    diff --color --unified $REPORT test/received/protein_test.received.txt
fi


REPORT="test/motif_test.approved.txt"
./bin/pepteditor test/motif_test.pepd | sed '/^#/d' >test/received/motif_test.received.txt
RESULT=$(diff --unified $REPORT test/received/motif_test.received.txt)
if [ -z "$RESULT" ]; then
    echo "Motif test succeeded."
else
    echo "Motif test FAILED."
    diff --color --unified $REPORT test/received/motif_test.received.txt
fi


bin/primarydock test/testTAAR8.config --colorless --pose 5 --iter 50 --congress > test/received/TAAR8_CAD.txt
TAAR_RESULT=$?
if [ "$TAAR_RESULT" -eq "0" ]; then
    POSES=$( cat test/received/TAAR8_CAD.txt | grep "pose(s) found" )
    if [ -z "$POSES" ]; then
        echo "TAAR test FAILED: no poses."
    else
        ASP111=$( cat test/received/TAAR8_CAD.txt | grep -m 1 "Asp111: " )
        ASP111="${ASP111/Asp111: /}"
        ASP111="${ASP111/[.][0-9]*/}"
        ASP201=$( cat test/received/TAAR8_CAD.txt | grep -m 1 "Asp201: " )
        ASP201="${ASP201/Asp201: /}"
        ASP201="${ASP201/[.][0-9]*/}"
        if [[ $ASP111 -gt "-30"  ]] || [[ $ASP201 -gt "-30"  ]]; then
            echo "TAAR test FAILED: bad contacts."
        else
            echo "TAAR test succeeded."
        fi
    fi
else
    echo "TAAR test FAILED: return value."
fi


bin/primarydock test/test1A1.config --colorless --pose 5 --iter 50 --congress > test/received/OR1A1_dLIMN.txt
DLIMN_RESULT=$?
if [ "$DLIMN_RESULT" -eq "0" ]; then
    POSES=$( cat test/received/OR1A1_dLIMN.txt | grep "pose(s) found" )
    if [ -z "$POSES" ]; then
        echo "d-limonene test FAILED: no poses."
    else
        echo "d-limonene test succeeded."
    fi
else
    echo "d-limonene test FAILED."
fi

