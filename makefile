all: obj/point.o obj/atom.o obj/intera.o obj/molecule.o obj/aminoacid.o obj/protein.o test/point_test test/atom_test test/mol_test test/mol_assem_test test/amino_test test/aniso_test test/protein_test test/bktest bin/metal bin/podock

# TODO: https://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

CC=g++

# Default CFLAG - no code coverage
CFLAGS=-g -fpermissive -Wwrite-strings -fextended-identifiers -std=c++14

# For code coverage instrumentation, switch to these CFLAGS (slower performance):
#CFLAGS=-g -fpermissive -Wwrite-strings -fextended-identifiers -std=c++14 -fprofile-arcs -ftest-coverage

clean:
	rm obj/*.o *.gcov *.gcno *.gcda

obj/point.o: src/classes/point.h src/classes/point.cpp src/classes/constants.h
	$(CC) -c src/classes/point.cpp -o obj/point.o $(CFLAGS)

obj/atom.o: src/classes/atom.h src/classes/atom.cpp
	$(CC) -c src/classes/atom.cpp -o obj/atom.o $(CFLAGS)

obj/intera.o: src/classes/intera.h src/classes/intera.cpp
	$(CC) -c src/classes/intera.cpp -o obj/intera.o $(CFLAGS)

obj/molecule.o: src/classes/molecule.h src/classes/molecule.cpp
	$(CC) -c src/classes/molecule.cpp -o obj/molecule.o $(CFLAGS)

obj/aminoacid.o: src/classes/aminoacid.h src/classes/aminoacid.cpp
	$(CC) -c src/classes/aminoacid.cpp -o obj/aminoacid.o $(CFLAGS)

obj/protein.o: src/classes/protein.h src/classes/protein.cpp
	$(CC) -c src/classes/protein.cpp -o obj/protein.o $(CFLAGS)

test/point_test: src/point_test.cpp obj/point.o
	$(CC) src/point_test.cpp obj/point.o -o test/point_test $(CFLAGS)

test/atom_test: src/atom_test.cpp obj/point.o obj/atom.o
	$(CC) src/atom_test.cpp obj/atom.o obj/point.o -o test/atom_test $(CFLAGS)

test/mol_test: src/mol_test.cpp obj/point.o obj/atom.o obj/molecule.o obj/intera.o
	$(CC) src/mol_test.cpp obj/atom.o obj/point.o obj/intera.o obj/molecule.o -o test/mol_test $(CFLAGS)

test/aniso_test: src/aniso_test.cpp obj/point.o obj/atom.o obj/molecule.o obj/intera.o
	$(CC) src/aniso_test.cpp obj/atom.o obj/point.o obj/intera.o obj/molecule.o -o test/aniso_test $(CFLAGS)

test/mol_assem_test: src/mol_assem_test.cpp obj/point.o obj/atom.o obj/molecule.o obj/intera.o
	$(CC) src/mol_assem_test.cpp obj/atom.o obj/point.o obj/intera.o obj/molecule.o -o test/mol_assem_test $(CFLAGS)

test/amino_test: src/amino_test.cpp obj/point.o obj/atom.o obj/molecule.o obj/intera.o obj/aminoacid.o
	$(CC) src/amino_test.cpp obj/atom.o obj/point.o obj/intera.o obj/molecule.o obj/aminoacid.o -o test/amino_test $(CFLAGS)

bin/metal: src/metal.cpp obj/point.o obj/atom.o obj/molecule.o obj/intera.o obj/aminoacid.o obj/protein.o
	$(CC) src/metal.cpp obj/atom.o obj/point.o obj/intera.o obj/molecule.o obj/aminoacid.o obj/protein.o -o bin/metal $(CFLAGS)

test/protein_test: src/protein_test.cpp obj/point.o obj/atom.o obj/molecule.o obj/intera.o obj/aminoacid.o obj/protein.o
	$(CC) src/protein_test.cpp obj/atom.o obj/point.o obj/intera.o obj/molecule.o obj/aminoacid.o obj/protein.o -o test/protein_test $(CFLAGS)

test/bktest: src/backbone_test.cpp obj/atom.o obj/point.o obj/intera.o obj/molecule.o obj/aminoacid.o obj/protein.o
	$(CC) src/backbone_test.cpp obj/atom.o obj/point.o obj/intera.o obj/molecule.o obj/aminoacid.o obj/protein.o -o test/bktest $(CFLAGS)

bin/podock: src/podock.cpp obj/point.o obj/atom.o obj/molecule.o obj/intera.o obj/aminoacid.o obj/protein.o
	$(CC) src/podock.cpp obj/atom.o obj/point.o obj/intera.o obj/molecule.o obj/aminoacid.o obj/protein.o -o bin/podock $(CFLAGS)

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
