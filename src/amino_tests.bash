AMINOS=$1
for (( i=0; i<${#AMINOS}; i++ )); do
  LETTER=${AMINOS:$i:1}
  echo "Running amino $LETTER test"
  REPORT=test/amino_test.$LETTER.approved.txt
  RECEIVED=test/amino_test.$LETTER.received.txt
  test/amino_test $LETTER | sed '/^#/d' >$RECEIVED
  # Ignoring sdf file because it is expected to vary and must be checked visually.
  #  echo "Content of test.sdf:" >>$RECEIVED
  #  sed "2d" test.sdf >>$RECEIVED
  diff --color --unified $REPORT $RECEIVED
done
