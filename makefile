OBJDIR=obj
BINDIR=bin
OUTDIR=output

DIRS=$(OBJDIR) $(BINDIR) $(OUTDIR)
OBJS=$(OBJDIR)/point.o $(OBJDIR)/atom.o $(OBJDIR)/intera.o $(OBJDIR)/molecule.o $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o
TESTS=test/point_test test/atom_test test/molecule_test test/mol_assem_test test/amino_test test/aniso_test test/protein_test test/bktest
APPS=$(BINDIR)/metal $(BINDIR)/podock

all: $(DIRS) \
	 $(OBJS) \
	 $(TESTS) \
	 $(APPS)

# TODO: https://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

CC=g++

# Default CFLAG - no code coverage
CFLAGS=-g -fpermissive -Wwrite-strings -fextended-identifiers -std=c++14


# For code coverage instrumentation, switch to these CFLAGS (slower performance):
#CFLAGS=-g -fpermissive -Wwrite-strings -fextended-identifiers -std=c++14 -fprofile-arcs -ftest-coverage

clean:
	rm $(OBJDIR)/*.o *.gcov *.gcno *.gcda

$(OBJDIR):
	if [ ! -f $(OBJDIR) ]; then mkdir -p $(OBJDIR); fi

$(BINDIR):
	if [ ! -f $(BINDIR) ]; then mkdir -p $(BINDIR); fi

$(OUTDIR):
	if [ ! -f $(OUTDIR) ]; then mkdir -p $(OUTDIR); fi

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

test/molecule_test: src/molecule_test.cpp $(OBJDIR)/point.o $(OBJDIR)/atom.o $(OBJDIR)/molecule.o $(OBJDIR)/intera.o
	$(CC) src/molecule_test.cpp $(OBJDIR)/atom.o $(OBJDIR)/point.o $(OBJDIR)/intera.o $(OBJDIR)/molecule.o -o test/molecule_test $(CFLAGS)

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
amino_report: REPORT="test/amino_test.approved.txt"
amino_report: test/amino_test
	bash src/amino_tests.bash ARNDCEQGHILKMFPUSTWYV

atom_report: REPORT="test/atom_test.approved.txt"
atom_report: test/atom_test
	./test/atom_test H >$(REPORT)

aniso_report: REPORT="test/aniso_test.approved.txt"
aniso_report: test/aniso_test
	./test/aniso_test >$(REPORT)

point_report: REPORT="test/point_test.approved.txt"
point_report: test/point_test
	./test/point_test >$(REPORT)

molecule_report: REPORT="test/molecule_test.approved.txt"
molecule_report: test/molecule_test
	./test/molecule_test 'CC(=O)[O-]' 'C[NH+](C)C' | sed '/^#/d' >$(REPORT)  # ignore lines starting with #

mol_assem_report: REPORT="test/mol_assem_test.approved.txt"
mol_assem_report: test/mol_assem_test
	./test/mol_assem_test >$(REPORT)
	echo "Content of test.sdf:" >> $(REPORT)
	sed '2d' test.sdf >> $(REPORT)  # remove line 2 (date stamp)

protein_report: REPORT="test/protein_test.approved.txt"
protein_report: test/protein_test
	./test/protein_test AAAAAAAAAA >$(REPORT)
	# Straight Strand PDB.
	echo "Content of test.pdb:" >> $(REPORT)
	cat test.pdb >> $(REPORT)
	# Alpha helix PDB.
	echo "Content of test_alpha.pdb:" >> $(REPORT)
	cat test_alpha.pdb >> $(REPORT)
	# Beta pleated PDB.
	echo "Content of test_beta.pdb:" >> $(REPORT)
	cat test_beta.pdb >> $(REPORT)
	# 3.10 helix PDB.
	echo "Content of test_310.pdb:" >> $(REPORT)
	cat test_310.pdb >> $(REPORT)
	# Pi helix PDB.
	echo "Content of test_pi.pdb:" >> $(REPORT)
	cat test_pi.pdb >> $(REPORT)
	# SDF of most recent PDB.
	echo "Content of test2.sdf:" >> $(REPORT)
	sed '2d' test2.sdf >> $(REPORT)

reports: amino_report atom_report aniso_report point_report molecule_report mol_assem_report protein_report
