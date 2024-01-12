#!/bin/bash
clear

PROT="$1"
LIG="$2"
DKR="$3"

if [ -z "$PROT" ]; then
    PROT="OR51E2"
fi

if [[ "$LIG" = "next" ]]; then
    LIG="next=1"
elif [ -z "$LIG" ]; then
    LIG="lig=propionic_acid"
else
    LIG="lig=$LIG"
fi

DKRARG=""
if [[ "$DKR" = "pd" ]]; then
    DKRARG="docker=$DKR"
elif [[ "$DKR" == "vina" ]]; then
    DKRARG="docker=$DKR"
fi

FIRST4=${PROT:0:4}

if [[ "$FIRST4" = "OR51" ]]; then
    php -f predict/method_directmdl.php "prot=$PROT" "$LIG" "$DKRARG"
elif [[ "$FIRST4" = "OR52" ]]; then
    php -f predict/method_directmdl.php "prot=$PROT" "$LIG" "$DKRARG"
else
    php -f predict/method_fygactive.php "prot=$PROT" "$LIG" "$DKRARG"
fi
