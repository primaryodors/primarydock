#!/bin/bash
clear

PROT="$1"
LIG="$2"

if [ -z "$PROT" ]; then
    PROT="PO51E2"
fi

if [ -z "$LIG" ]; then
    LIG="propionic_acid"
fi

php -f predict/method_icactive.php "prot=$PROT" "lig=$LIG"
