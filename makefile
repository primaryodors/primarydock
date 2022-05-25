all: point.o atom.o intera.o molecule.o aminoacid.o protein.o point_test atom_test mol_test mol_assem_test amino_test aniso_test protest bktest metal podock

# TODO: https://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

CC=g++
CFLAGS=-g -fpermissive -fextended-identifiers -std=c++11

clean:
	rm *.o

point.o: point.h point.cpp constants.h
	$(CC) -c point.cpp -o point.o $(CFLAGS)

atom.o: atom.h atom.cpp
	$(CC) -c atom.cpp -o atom.o $(CFLAGS)

intera.o: intera.h intera.cpp
	$(CC) -c intera.cpp -o intera.o $(CFLAGS)

molecule.o: molecule.h molecule.cpp
	$(CC) -c molecule.cpp -o molecule.o $(CFLAGS)

aminoacid.o: aminoacid.h aminoacid.cpp
	$(CC) -c aminoacid.cpp -o aminoacid.o $(CFLAGS)

protein.o: protein.h protein.cpp
	$(CC) -c protein.cpp -o protein.o $(CFLAGS)

point_test: point_test.cpp point.o
	$(CC) point_test.cpp point.o -o point_test $(CFLAGS)

atom_test: atom_test.cpp point.o atom.o
	$(CC) atom_test.cpp atom.o point.o -o atom_test $(CFLAGS)

mol_test: mol_test.cpp point.o atom.o molecule.o intera.o
	$(CC) mol_test.cpp atom.o point.o intera.o molecule.o -o mol_test $(CFLAGS)

aniso_test: aniso_test.cpp point.o atom.o molecule.o intera.o
	$(CC) aniso_test.cpp atom.o point.o intera.o molecule.o -o aniso_test $(CFLAGS)

mol_assem_test: mol_assem_test.cpp point.o atom.o molecule.o intera.o
	$(CC) mol_assem_test.cpp atom.o point.o intera.o molecule.o -o mol_assem_test $(CFLAGS)

amino_test: amino_test.cpp point.o atom.o molecule.o intera.o aminoacid.o
	$(CC) amino_test.cpp atom.o point.o intera.o molecule.o aminoacid.o -o amino_test $(CFLAGS)

metal: metal.cpp point.o atom.o molecule.o intera.o aminoacid.o protein.o
	$(CC) metal.cpp atom.o point.o intera.o molecule.o aminoacid.o protein.o -o metal $(CFLAGS)

protest: protest.cpp point.cpp atom.o molecule.o intera.o aminoacid.o protein.o
	$(CC) protest.cpp atom.o point.o intera.o molecule.o aminoacid.o protein.o -o protest $(CFLAGS)

bktest: backbone_test.cpp atom.o point.o intera.o molecule.o aminoacid.o protein.o
	$(CC) backbone_test.cpp atom.o point.o intera.o molecule.o aminoacid.o protein.o -o bktest $(CFLAGS)

podock: podock.cpp point.cpp atom.o molecule.o intera.o aminoacid.o protein.o
	$(CC) podock.cpp atom.o point.o intera.o molecule.o aminoacid.o protein.o -o podock $(CFLAGS)




