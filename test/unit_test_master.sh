#!/bin/bash
cd "$(dirname "$0")"
if [ ! -f "received" ]; then mkdir -p "received"; fi
cd ..

make code > test/received/make_code.txt
CODE_RESULT=$?
if [ $CODE_RESULT -eq 0 ]; then
    echo "Make code succeeded."
else
    echo "Make code FAILED."
fi

# TODO: Some way to evaluate whether dock results exist and have good enough contacts on Asp111 and Asp201.
bin/primarydock test/testTAAR8.config --colorless --pose 2 --iter 25 --congress > test/received/TAAR8_CAD.txt
TAAR_RESULT=$?
if [ "$TAAR_RESULT" -eq "0" ]; then
    echo "TAAR test succeeded."
else
    echo "TAAR test FAILED."
fi

# TODO: Some way to evaluate whether at least one dock result was returned.
bin/primarydock test/test1A1.config --colorless --pose 3 --iter 30 --congress > test/received/TAAR8_CAD.txt
DLIMN_RESULT=$?
if [ "$DLIMN_RESULT" -eq "0" ]; then
    echo "d-limonene test succeeded."
else
    echo "d-limonene test FAILED."
fi

