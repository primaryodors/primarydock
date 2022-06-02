AMINOS=$1
for (( i=0; i<${#AMINOS}; i++ )); do
  LETTER=${AMINOS:$i:1}
  echo "Running amino $LETTER test"
  REPORT=test/amino_test.$LETTER.approved.txt
  ./amino_test $LETTER >$REPORT
  echo "Content of test.pdb:" >>$REPORT
  cat test.pdb >>$REPORT
  echo "Content of test.sdf:" >>$REPORT
  sed "2d" test.sdf >>$REPORT
done
