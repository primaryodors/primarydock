#!/bin/bash
clear

PROT="$1"
LIG="$2"

if [ -z "$PROT" ]; then
    PROT="OR51E2"
fi

if [ -z "$LIG" ]; then
    LIG="propionic_acid"
fi

php -f predict/method_fygactive.php "prot=$PROT" "lig=$LIG"
