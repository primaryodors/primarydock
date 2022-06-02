OBJDIR=obj
BINDIR=bin

all: $(OBJDIR) $(BINDIR) $(OBJDIR)/point.o $(OBJDIR)/atom.o $(OBJDIR)/intera.o $(OBJDIR)/molecule.o $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o test/point_test test/atom_test test/mol_test test/mol_assem_test test/amino_test test/aniso_test test/protein_test test/bktest $(BINDIR)/metal $(BINDIR)/podock

# TODO: https://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

CC=g++

# Default CFLAG - no code coverage
CFLAGS=-g -fpermissive -Wwrite-strings -fextended-identifiers -std=c++14

# For code coverage instrumentation, switch to these CFLAGS (slower performance):
#CFLAGS=-g -fpermissive -Wwrite-strings -fextended-identifiers -std=c++14 -fprofile-arcs -ftest-coverage

clean:
	rm $(OBJDIR)/*.o *.gcov *.gcno *.gcda

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

$(OBJDIR)/point.o: src/classes/point.h src/classes/point.cpp src/classes/constants.h
	$(CC) -c src/classes/point.cpp -o $(OBJDIR)/point.o $(CFLAGS)

$(OBJDIR)/atom.o: src/classes/atom.h src/classes/atom.cpp
	$(CC) -c src/classes/atom.cpp -o $(OBJDIR)/atom.o $(CFLAGS)

$(OBJDIR)/intera.o: src/classes/intera.h src/classes/intera.cpp
	$(CC) -c src/classes/intera.cpp -o $(OBJDIR)/intera.o $(CFLAGS)

$(OBJDIR)/molecule.o: src/classes/molecule.h src/classes/molecule.cpp
	$(CC) -c src/classes/molecule.cpp -o $(OBJDIR)/molecule.o $(CFLAGS)

$(OBJDIR)/aminoacid.o: src/classes/aminoacid.h src/classes/aminoacid.cpp
	$(CC) -c src/classes/aminoacid.cpp -o $(OBJDIR)/aminoacid.o $(CFLAGS)

$(OBJDIR)/protein.o: src/classes/protein.h src/classes/protein.cpp
	$(CC) -c src/classes/protein.cpp -o $(OBJDIR)/protein.o $(CFLAGS)

test/point_test: src/point_test.cpp $(OBJDIR)/point.o
	$(CC) src/point_test.cpp $(OBJDIR)/point.o -o test/point_test $(CFLAGS)

test/atom_test: src/atom_test.cpp $(OBJDIR)/point.o $(OBJDIR)/atom.o
	$(CC) src/atom_test.cpp $(OBJDIR)/atom.o $(OBJDIR)/point.o -o test/atom_test $(CFLAGS)

test/mol_test: src/mol_test.cpp $(OBJDIR)/point.o $(OBJDIR)/atom.o $(OBJDIR)/molecule.o $(OBJDIR)/intera.o
	$(CC) src/mol_test.cpp $(OBJDIR)/atom.o $(OBJDIR)/point.o $(OBJDIR)/intera.o $(OBJDIR)/molecule.o -o test/mol_test $(CFLAGS)

test/aniso_test: src/aniso_test.cpp $(OBJDIR)/point.o $(OBJDIR)/atom.o $(OBJDIR)/molecule.o $(OBJDIR)/intera.o
	$(CC) src/aniso_test.cpp $(OBJDIR)/atom.o $(OBJDIR)/point.o $(OBJDIR)/intera.o $(OBJDIR)/molecule.o -o test/aniso_test $(CFLAGS)

test/mol_assem_test: src/mol_assem_test.cpp $(OBJDIR)/point.o $(OBJDIR)/atom.o $(OBJDIR)/molecule.o $(OBJDIR)/intera.o
	$(CC) src/mol_assem_test.cpp $(OBJDIR)/atom.o $(OBJDIR)/point.o $(OBJDIR)/intera.o $(OBJDIR)/molecule.o -o test/mol_assem_test $(CFLAGS)

test/amino_test: src/amino_test.cpp $(OBJDIR)/point.o $(OBJDIR)/atom.o $(OBJDIR)/molecule.o $(OBJDIR)/intera.o $(OBJDIR)/aminoacid.o
	$(CC) src/amino_test.cpp $(OBJDIR)/atom.o $(OBJDIR)/point.o $(OBJDIR)/intera.o $(OBJDIR)/molecule.o $(OBJDIR)/aminoacid.o -o test/amino_test $(CFLAGS)

$(BINDIR)/metal: src/metal.cpp $(OBJDIR)/point.o $(OBJDIR)/atom.o $(OBJDIR)/molecule.o $(OBJDIR)/intera.o $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o
	$(CC) src/metal.cpp $(OBJDIR)/atom.o $(OBJDIR)/point.o $(OBJDIR)/intera.o $(OBJDIR)/molecule.o $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o -o $(BINDIR)/metal $(CFLAGS)

test/protein_test: src/protein_test.cpp $(OBJDIR)/point.o $(OBJDIR)/atom.o $(OBJDIR)/molecule.o $(OBJDIR)/intera.o $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o
	$(CC) src/protein_test.cpp $(OBJDIR)/atom.o $(OBJDIR)/point.o $(OBJDIR)/intera.o $(OBJDIR)/molecule.o $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o -o test/protein_test $(CFLAGS)

test/bktest: src/backbone_test.cpp $(OBJDIR)/atom.o $(OBJDIR)/point.o $(OBJDIR)/intera.o $(OBJDIR)/molecule.o $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o
	$(CC) src/backbone_test.cpp $(OBJDIR)/atom.o $(OBJDIR)/point.o $(OBJDIR)/intera.o $(OBJDIR)/molecule.o $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o -o test/bktest $(CFLAGS)

$(BINDIR)/podock: src/podock.cpp $(OBJDIR)/point.o $(OBJDIR)/atom.o $(OBJDIR)/molecule.o $(OBJDIR)/intera.o $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o
	$(CC) src/podock.cpp $(OBJDIR)/atom.o $(OBJDIR)/point.o $(OBJDIR)/intera.o $(OBJDIR)/molecule.o $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o -o $(BINDIR)/podock $(CFLAGS)

performance_test: $(BINDIR)/podock testdata/test_TAAR8.config testdata/TAAR8.rotated.pdb testdata/CAD_ion.sdf
	./$(BINDIR)/podock testdata/test_TAAR8.config

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
