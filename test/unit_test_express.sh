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



