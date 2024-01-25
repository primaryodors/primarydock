#!/bin/bash
clear

PROT="$1"
LIG="$2"
DKR="$3"

if [ -z "$4" ]; then
    SMILES="$3"
    DKR="$4"
fi

if [ -z "$PROT" ]; then
    PROT="OR51E2"
fi

if [[ "$LIG" = "next" ]]; then
    LIGARG="next=1"
elif [ -z "$LIG" ]; then
    LIGARG="lig=propionic_acid"
else
    LIGF="sdf/$LIG.sdf"
    LIGARG="lig=$LIG"
    if ! test -f "$LIGF"; then
        if [ -z "$SMILES" ]; then
            echo "Odorant not found; please specify SMILES string for third argument."
            exit
        else
            php -f "data/gensdf.php" "$LIG" "$SMILES"
        fi
    fi
fi

DKRARG=""
if [[ "$DKR" = "pd" ]]; then
    DKRARG="docker=$DKR"
elif [[ "$DKR" == "vina" ]]; then
    DKRARG="docker=$DKR"
fi

FIRST4=${PROT:0:4}

if [[ "$FIRST4" = "OR51" ]]; then
    php -f predict/method_directmdl.php "prot=$PROT" "$LIGARG" "$DKRARG"
elif [[ "$FIRST4" = "OR52" ]]; then
    php -f predict/method_directmdl.php "prot=$PROT" "$LIGARG" "$DKRARG"
else
    php -f predict/method_fygactive.php "prot=$PROT" "$LIGARG" "$DKRARG"
fi
