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

FIRST4=${PROT:0:4}

if [[ "$FIRST4" = "OR51" ]]; then
    php -f predict/method_directmdl.php "prot=$PROT" "lig=$LIG"
elif [[ "$FIRST4" = "OR52" ]]; then
    php -f predict/method_directmdl.php "prot=$PROT" "lig=$LIG"
else
    php -f predict/method_fygactive.php "prot=$PROT" "lig=$LIG"
fi
