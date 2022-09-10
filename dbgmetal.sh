
RCPID="OR1A1"

clear
make code

wget -O 2bmetald.dat http://www.primaryodors.org/recepinfo.php?recep=$RCPID
wget -O.upright.pdb http://www.primaryodors.org/orpdb.php?recep=$RCPID

#valgrind ./bin/metal -p.upright.pdb -o output/$RCPID.metal.pdb -e Zn -c 2 -i 2bmetald.dat
./bin/metal -p.upright.pdb -o output/$RCPID.metal.pdb -e Zn -c 2 -i 2bmetald.dat

