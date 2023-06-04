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
        if [[ $ASP111 -gt "-35"  ]] || [[ $ASP201 -gt "-35"  ]]; then
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

