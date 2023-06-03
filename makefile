OBJDIR=obj
BINDIR=bin
OUTDIR=output
SDFDIR=sdf

DIRS=$(OBJDIR) $(BINDIR) $(OUTDIR) $(SDFDIR)
OBJS=$(OBJDIR)/misc.o $(OBJDIR)/point.o $(OBJDIR)/atom.o $(OBJDIR)/intera.o $(OBJDIR)/molecule.o $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o $(OBJDIR)/group.o
TESTS=test/point_test test/atom_test test/molecule_test test/pi_stack_test test/mol_assem_test test/aniso_test \
	  test/group_test_mol test/protein_test test/backbone_test
APPS=$(BINDIR)/primarydock $(BINDIR)/pepteditor
REPORTS=amino_report atom_report aniso_report point_report molecule_report mol_assem_report protein_report
all: $(DIRS) \
	 $(OBJS) \
	 $(TESTS) \
	 $(APPS) \
	 $(REPORTS)
code: $(OBJS) $(TESTS) molecule_report $(APPS)
primarydock: $(DIRS) $(OBJS) $(BINDIR)/primarydock
pepteditor: $(DIRS) $(OBJS) $(BINDIR)/pepteditor

CC=g++

# Default CFLAG - no code coverage
CFLAGS=-g -Wwrite-strings -fextended-identifiers -std=c++14

# For code coverage instrumentation, switch to these CFLAGS (slower performance):
#CFLAGS=-g -Wwrite-strings -fextended-identifiers -std=c++14 -fprofile-arcs -ftest-coverage

clean:
	rm $(OBJDIR)/*.o *.gcov *.gcno *.gcda

$(OBJDIR):
	if [ ! -f $(OBJDIR) ]; then mkdir -p $(OBJDIR); fi

$(BINDIR):
	if [ ! -f $(BINDIR) ]; then mkdir -p $(BINDIR); fi

$(OUTDIR):
	if [ ! -f $(OUTDIR) ]; then mkdir -p $(OUTDIR); fi

$(SDFDIR):
	if [ ! -f $(SDFDIR) ]; then mkdir -p $(SDFDIR); fi

$(OBJDIR)/misc.o: src/classes/misc.h src/classes/misc.cpp src/classes/constants.h
	$(CC) -c src/classes/misc.cpp -o $(OBJDIR)/misc.o $(CFLAGS)

$(OBJDIR)/point.o: src/classes/point.h src/classes/point.cpp $(OBJDIR)/misc.o src/classes/constants.h
	$(CC) -c src/classes/point.cpp -o $(OBJDIR)/point.o $(CFLAGS)

$(OBJDIR)/atom.o: src/classes/atom.h src/classes/atom.cpp $(OBJDIR)/point.o
	$(CC) -c src/classes/atom.cpp -o $(OBJDIR)/atom.o $(CFLAGS)

$(OBJDIR)/intera.o: src/classes/intera.h src/classes/intera.cpp $(OBJDIR)/atom.o
	$(CC) -c src/classes/intera.cpp -o $(OBJDIR)/intera.o $(CFLAGS)

$(OBJDIR)/molecule.o: src/classes/molecule.h src/classes/molecule.cpp $(OBJDIR)/intera.o
	$(CC) -c src/classes/molecule.cpp -o $(OBJDIR)/molecule.o $(CFLAGS)

$(OBJDIR)/aminoacid.o: src/classes/aminoacid.h src/classes/aminoacid.cpp $(OBJDIR)/molecule.o
	$(CC) -c src/classes/aminoacid.cpp -o $(OBJDIR)/aminoacid.o $(CFLAGS)

$(OBJDIR)/protein.o: src/classes/protein.h src/classes/protein.cpp $(OBJDIR)/aminoacid.o
	$(CC) -c src/classes/protein.cpp -o $(OBJDIR)/protein.o $(CFLAGS)

$(OBJDIR)/group.o: src/classes/group.h src/classes/group.cpp $(OBJDIR)/protein.o
	$(CC) -c src/classes/group.cpp -o $(OBJDIR)/group.o $(CFLAGS)

test/point_test: src/point_test.cpp $(OBJDIR)/point.o
	$(CC) src/point_test.cpp $(OBJDIR)/point.o $(OBJDIR)/misc.o -o test/point_test $(CFLAGS)

test/atom_test: src/atom_test.cpp $(OBJDIR)/point.o $(OBJDIR)/atom.o
	$(CC) src/atom_test.cpp $(OBJDIR)/misc.o $(OBJDIR)/atom.o $(OBJDIR)/point.o -o test/atom_test $(CFLAGS)

test/molecule_test: src/molecule_test.cpp $(OBJS)
	$(CC) src/molecule_test.cpp $(OBJS) -o test/molecule_test $(CFLAGS)

test/pi_stack_test: src/pi_stack_test.cpp $(OBJS)
	$(CC) src/pi_stack_test.cpp $(OBJS) -o test/pi_stack_test $(CFLAGS)

test/aniso_test: src/aniso_test.cpp $(OBJS)
	$(CC) src/aniso_test.cpp $(OBJS) -o test/aniso_test $(CFLAGS)

test/group_test_mol: src/group_test_mol.cpp $(OBJS)
	$(CC) src/group_test_mol.cpp $(OBJS) -o test/group_test_mol $(CFLAGS)

test/mol_assem_test: src/mol_assem_test.cpp $(OBJS)
	$(CC) src/mol_assem_test.cpp $(OBJS) -o test/mol_assem_test $(CFLAGS)

test/amino_test: src/amino_test.cpp $(OBJS) $(OBJDIR)/aminoacid.o
	$(CC) src/amino_test.cpp $(OBJS) -o test/amino_test $(CFLAGS)

test/protein_test: src/protein_test.cpp $(OBJS) $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o
	$(CC) src/protein_test.cpp $(OBJS) -o test/protein_test $(CFLAGS)

test/backbone_test: src/backbone_test.cpp $(OBJS) $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o
	$(CC) src/backbone_test.cpp $(OBJS) -o test/backbone_test $(CFLAGS)

$(BINDIR)/primarydock: src/primarydock.cpp $(OBJS) $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o $(OBJDIR)/group.o
	$(CC) src/primarydock.cpp $(OBJS) -o $(BINDIR)/primarydock $(CFLAGS)

$(BINDIR)/pepteditor: src/interpreter.cpp $(OBJS) $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o $(OBJDIR)/group.o
	$(CC) src/interpreter.cpp $(OBJS) -o $(BINDIR)/pepteditor $(CFLAGS)

performance_test: $(BINDIR)/primarydock testdata/test_TAAR8.config testdata/TAAR8.upright.pdb testdata/CAD_ion.sdf
	./$(BINDIR)/primarydock testdata/test_TAAR8.config

# low-tooling regression tests below
amino_report: REPORT="test/amino_test.approved.txt"
amino_report: test/amino_test
	./test/amino_test >test/amino_test.received.txt
	diff --color --unified $(REPORT) test/amino_test.received.txt

atom_report: REPORT="test/atom_test.approved.txt"
atom_report: test/atom_test
	./test/atom_test H >test/atom_test.received.txt
	diff --color --unified $(REPORT) test/atom_test.received.txt

aniso_report: REPORT="test/aniso_test.approved.txt"
aniso_report: test/aniso_test
	./test/aniso_test >test/aniso_test.received.txt
	diff --color --unified $(REPORT) test/aniso_test.received.txt

point_report: REPORT="test/point_test.approved.txt"
point_report: test/point_test
	./test/point_test >test/point_test.received.txt
	diff --color --unified $(REPORT) test/point_test.received.txt

molecule_report: REPORT="test/molecule_test1.approved.txt"
molecule_report: test/molecule_test
	./test/molecule_test 'CCO' 'CCO' | tee temp | sed '/^#/d' >test/molecule_test1.received.txt; cat temp # ignore lines starting with #
	diff --color --unified $(REPORT) test/molecule_test1.received.txt

mol_assem_report: REPORT="test/mol_assem_test.approved.txt"
mol_assem_report: test/mol_assem_test
	bash ./test/mol_assem_tests.sh	# these must be checked visually.
	./test/mol_assem_test >test/molecule_test.received.txt
	echo "Content of test.sdf:" >> test/molecule_test.received.txt
	sed '2d' test.sdf >> test/molecule_test.received.txt  # remove line 2 (date stamp)
	diff --color --unified $(REPORT) test/molecule_test.received.txt

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

reports: $(REPORTS)
