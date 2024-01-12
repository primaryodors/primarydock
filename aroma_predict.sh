#!/bin/bash

LIG="$1"
DKR="$2"

DKRARG=""
if [[ "$DKR" = "pd" ]]; then
    DKRARG="$DKR"
elif [[ "$DKR" == "vina" ]]; then
    DKRARG="$DKR"
fi

while read ln; do
  PROT=${ln%%[[:space:]]*}
  if [[ "$PROT" =~ ^(OR[0-9]{1,2}[A-Z]{1,2}[0-9]{1,2}|TAAR[0-9]) ]]; then
    echo "$PROT"
    bash ./run_prediction.sh "$PROT" "$LIG" "$DKRARG"
  fi
done <data/sequences_aligned.txt
