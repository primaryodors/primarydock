AMINOS=$1
for (( i=0; i<${#AMINOS}; i++ )); do
  LETTER=${AMINOS:$i:1}
  echo "Running amino $LETTER test"
  REPORT=test/amino_test.$LETTER.approved.txt
  test/amino_test $LETTER >$REPORT
  # Ignoring sdf file for now - varies
  #  echo "Content of test.sdf:" >>$REPORT
  #  sed "2d" test.sdf >>$REPORT
done
