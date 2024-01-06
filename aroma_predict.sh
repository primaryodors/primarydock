#!/bin/bash

LIG="$1"

while read ln; do
  PROT=${ln%%[[:space:]]*}
  if [[ "$PROT" =~ ^(OR[0-9]{1,2}[A-Z]{1,2}[0-9]{1,2}|TAAR[0-9]) ]]; then
    echo "$PROT"
    bash ./run_prediction.sh "$PROT" "$LIG"
  fi
done <data/sequences_aligned.txt
