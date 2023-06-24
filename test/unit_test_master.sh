#!/bin/bash
cd "$(dirname "$0")"
if [ ! -f "received" ]; then mkdir -p "received"; fi
cd ..

clear

make code > test/received/make_code.txt
CODE_RESULT=$?
if [ $CODE_RESULT -eq 0 ]; then
    echo "Make code succeeded."
else
    echo "Make code FAILED."
fi


REPORT="test/point_test.approved.txt"
./test/point_test >test/point_test.received.txt
RESULT=$(diff --color --unified $REPORT test/point_test.received.txt)
if [ -z "$RESULT" ]; then
    echo "Point test succeeded."
else
    echo "Point test FAILED."
fi


REPORT="test/atom_test.approved.txt"
./test/atom_test H >test/atom_test.received.txt
RESULT=$(diff --unified $REPORT test/atom_test.received.txt)
if [ -z "$RESULT" ]; then
    echo "Atom test succeeded."
else
    echo "Atom test FAILED."
    diff --color --unified $REPORT test/atom_test.received.txt
fi


REPORT="test/aniso_test.approved.txt"
./test/aniso_test --asciiart >test/aniso_test.received.txt
RESULT=$(diff --unified $REPORT test/aniso_test.received.txt)
if [ -z "$RESULT" ]; then
    echo "Anisotropy test succeeded."
else
    echo "Anisotropy test FAILED."
    diff --color --unified $REPORT test/aniso_test.received.txt
fi


# TODO: This is where mol_assem_test would go, but currently its output must be checked visually.


REPORT="test/amino_test.approved.txt"
./test/amino_test >test/amino_test.received.txt
RESULT=$(diff --color --unified $REPORT test/amino_test.received.txt)
if [ -z "$RESULT" ]; then
    echo "Amino test succeeded."
else
    echo "Amino test FAILED."
fi


# TODO: A way to check these automatically.
REPORT="test/protein_test.approved.txt"
./test/protein_test AAAAAAAAAA > test/protein_test.received.txt
# Straight Strand PDB.
echo "Content of test.pdb:" >> $REPORT
cat test.pdb >> $REPORT
# Alpha helix PDB.
echo "Content of test_alpha.pdb:" >> $REPORT
cat test_alpha.pdb >> $REPORT
# Beta pleated PDB.
echo "Content of test_beta.pdb:" >> $REPORT
cat test_beta.pdb >> $REPORT
# 3.10 helix PDB.
echo "Content of test_310.pdb:" >> $REPORT
cat test_310.pdb >> $REPORT
# Pi helix PDB.
echo "Content of test_pi.pdb:" >> $REPORT
cat test_pi.pdb >> $REPORT
# SDF of most recent PDB.
echo "Content of test2.sdf:" >> $REPORT
sed '2d' test2.sdf >> $REPORT


REPORT="test/motif_test.approved.txt"
./bin/pepteditor test/motif_test.pepd >test/motif_test.received.txt
RESULT=$(diff --color --unified $REPORT test/motif_test.received.txt)
if [ -z "$RESULT" ]; then
    echo "Motif test succeeded."
else
    echo "Motif test FAILED."
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

