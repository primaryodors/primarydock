all: point.o atom.o intera.o molecule.o aminoacid.o protein.o point_test atom_test mol_test mol_assem_test amino_test aniso_test protest bktest metal podock

# TODO: https://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

CC=g++

# Default CFLAG - no code coverage
CFLAGS=-g -fpermissive -Wwrite-strings -fextended-identifiers -std=c++14

# For code coverage instrumentation, switch to these CFLAGS (slower performance):
#CFLAGS=-g -fpermissive -Wwrite-strings -fextended-identifiers -std=c++14 -fprofile-arcs -ftest-coverage

clean:
	rm *.o *.gcov *.gcno *.gcda

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

protein_test: protein_test.cpp point.o atom.o molecule.o intera.o aminoacid.o protein.o
	$(CC) protein_test.cpp atom.o point.o intera.o molecule.o aminoacid.o protein.o -o protein_test $(CFLAGS)

bktest: backbone_test.cpp atom.o point.o intera.o molecule.o aminoacid.o protein.o
	$(CC) backbone_test.cpp atom.o point.o intera.o molecule.o aminoacid.o protein.o -o bktest $(CFLAGS)

podock: podock.cpp point.o atom.o molecule.o intera.o aminoacid.o protein.o
	$(CC) podock.cpp atom.o point.o intera.o molecule.o aminoacid.o protein.o -o podock $(CFLAGS)

performance_test: podock testdata/test_TAAR8.config testdata/TAAR8.rotated.pdb testdata/CAD_ion.sdf
	./podock testdata/test_TAAR8.config

# low-tooling regression tests below
amino_report: REPORT="amino_test.approved.txt"
amino_report: amino_test
	bash amino_tests.bash ARNDCEQGHILKMFPUSTWYV

atom_report: REPORT="atom_test.approved.txt"
atom_report: atom_test
	./atom_test H >$(REPORT)

aniso_report: REPORT="aniso_test.approved.txt"
aniso_report: aniso_test
	./aniso_test >$(REPORT)

point_report: REPORT="point_test.approved.txt"
point_report: point_test
	./point_test >$(REPORT)

mol_report: REPORT="mol_test.approved.txt"
mol_report: mol_test
	./mol_test >$(REPORT)
	echo "Content of output.sdf:" >> $(REPORT)
	sed '2d' output.sdf >> $(REPORT)

mol_assem_report: REPORT="mol_assem_test.approved.txt"
mol_assem_report: mol_assem_test
	./mol_assem_test >$(REPORT)
	echo "Content of test.sdf:" >> $(REPORT)
	sed '2d' test.sdf >> $(REPORT)  # remove line 2 (date stamp)

#ARNDCEQGHILKMFPUSTWYV

protein_report: REPORT="protein_test.approved.txt"
protein_report: protein_test
	./protein_test AAAAAAAAAA >$(REPORT)
	echo "Content of test.pdb:" >> $(REPORT)
	cat test.pdb >> $(REPORT)
	echo "Content of test1.pdb:" >> $(REPORT)
	cat test1.pdb >> $(REPORT)
	echo "Content of test2.sdf:" >> $(REPORT)
	sed '2d' test2.sdf >> $(REPORT)

reports: amino_report atom_test aniso_report point_report mol_report mol_assem_report protein_report
