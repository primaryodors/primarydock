OBJDIR=obj
BINDIR=bin
OUTDIR=output
SDFDIR=sdf
TMPDIR=tmp

DIRS=$(OBJDIR) $(BINDIR) $(OUTDIR) $(SDFDIR) $(TMPDIR)
OBJS=$(OBJDIR)/misc.o $(OBJDIR)/point.o $(OBJDIR)/atom.o $(OBJDIR)/intera.o $(OBJDIR)/molecule.o $(OBJDIR)/aminoacid.o \
	$(OBJDIR)/protein.o $(OBJDIR)/group.o $(OBJDIR)/dynamic.o $(OBJDIR)/moiety.o $(OBJDIR)/scoring.o $(OBJDIR)/conj.o \
	$(OBJDIR)/search.o
TESTS=test/point_test test/atom_test test/molecule_test test/pi_stack_test test/mol_assem_test test/aniso_test test/amino_test \
	  test/group_test_mol test/group_test_res test/protein_test test/backbone_test test/bond_rotation_test test/moiety_test \
	  test/flexion_test test/histidine_test test/ring_test test/eclipsing_test test/cs_test test/mcoord_test
APPS=$(BINDIR)/primarydock $(BINDIR)/pepteditor $(BINDIR)/ic \
	 $(BINDIR)/score_pdb $(BINDIR)/ramachandran $(BINDIR)/ringflip
REPORTS=amino_report atom_report aniso_report point_report molecule_report mol_assem_report protein_report motif_report
all: $(DIRS) \
	 $(OBJS) \
	 $(TESTS) \
	 $(APPS)
code: $(DIRS) $(OBJS) $(TESTS) $(APPS)
primarydock: $(DIRS) $(OBJS) $(BINDIR)/primarydock
pepteditor: $(DIRS) $(OBJS) $(BINDIR)/pepteditor
ic: $(DIRS) $(OBJS) $(BINDIR)/ic

CC=g++

# Default CFLAGS - release mode
# CFLAGS=-O -Wwrite-strings -fextended-identifiers -std=c++14

# Debug CFLAGS - allows gdb, valgrind
CFLAGS=-g -Wwrite-strings -fextended-identifiers -std=c++14

# For code coverage instrumentation, switch to these CFLAGS (slower performance):
# CFLAGS=-g -Wwrite-strings -fextended-identifiers -std=c++14 -fprofile-arcs -ftest-coverage

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

$(TMPDIR):
	if [ ! -f $(TMPDIR) ]; then mkdir -p $(TMPDIR); fi

$(OBJDIR)/misc.o: src/classes/misc.h src/classes/misc.cpp src/classes/constants.h
	$(CC) -c src/classes/misc.cpp -o $(OBJDIR)/misc.o $(CFLAGS)

$(OBJDIR)/point.o: src/classes/point.h src/classes/point.cpp $(OBJDIR)/misc.o src/classes/constants.h
	$(CC) -c src/classes/point.cpp -o $(OBJDIR)/point.o $(CFLAGS)

$(OBJDIR)/atom.o: src/classes/atom.h src/classes/atom.cpp $(OBJDIR)/point.o
	$(CC) -c src/classes/atom.cpp -o $(OBJDIR)/atom.o $(CFLAGS)

$(OBJDIR)/conj.o: src/classes/conj.h src/classes/conj.cpp $(OBJDIR)/atom.o
	$(CC) -c src/classes/conj.cpp -o $(OBJDIR)/conj.o $(CFLAGS)

$(OBJDIR)/intera.o: src/classes/intera.h src/classes/intera.cpp $(OBJDIR)/conj.o
	$(CC) -c src/classes/intera.cpp -o $(OBJDIR)/intera.o $(CFLAGS)

$(OBJDIR)/molecule.o: src/classes/molecule.h src/classes/molecule.cpp $(OBJDIR)/intera.o
	$(CC) -c src/classes/molecule.cpp -o $(OBJDIR)/molecule.o $(CFLAGS)

$(OBJDIR)/aminoacid.o: src/classes/aminoacid.h src/classes/aminoacid.cpp $(OBJDIR)/molecule.o
	$(CC) -c src/classes/aminoacid.cpp -o $(OBJDIR)/aminoacid.o $(CFLAGS)

$(OBJDIR)/protein.o: src/classes/protein.h src/classes/protein.cpp $(OBJDIR)/aminoacid.o
	$(CC) -c src/classes/protein.cpp -o $(OBJDIR)/protein.o $(CFLAGS)

$(OBJDIR)/group.o: src/classes/group.h src/classes/group.cpp $(OBJDIR)/protein.o
	$(CC) -c src/classes/group.cpp -o $(OBJDIR)/group.o $(CFLAGS)

$(OBJDIR)/search.o: src/classes/search.h src/classes/search.cpp $(OBJDIR)/group.o
	$(CC) -c src/classes/search.cpp -o $(OBJDIR)/search.o $(CFLAGS)

$(OBJDIR)/dynamic.o: src/classes/dynamic.h src/classes/dynamic.cpp $(OBJDIR)/protein.o
	$(CC) -c src/classes/dynamic.cpp -o $(OBJDIR)/dynamic.o $(CFLAGS)

$(OBJDIR)/moiety.o: src/classes/moiety.h src/classes/moiety.cpp $(OBJDIR)/molecule.o
	$(CC) -c src/classes/moiety.cpp -o $(OBJDIR)/moiety.o $(CFLAGS)

$(OBJDIR)/scoring.o: src/classes/scoring.h src/classes/scoring.cpp $(OBJDIR)/protein.o
	$(CC) -c src/classes/scoring.cpp -o $(OBJDIR)/scoring.o $(CFLAGS)

test/point_test: src/test/point_test.cpp $(OBJDIR)/point.o
	$(CC) src/test/point_test.cpp $(OBJDIR)/point.o $(OBJDIR)/misc.o -o test/point_test $(CFLAGS)

test/atom_test: src/test/atom_test.cpp $(OBJDIR)/point.o $(OBJDIR)/atom.o
	$(CC) src/test/atom_test.cpp $(OBJDIR)/misc.o $(OBJDIR)/atom.o $(OBJDIR)/point.o -o test/atom_test $(CFLAGS)

test/molecule_test: src/test/molecule_test.cpp $(OBJS)
	$(CC) src/test/molecule_test.cpp $(OBJS) -o test/molecule_test $(CFLAGS)

test/ring_test: src/test/ring_test.cpp $(OBJDIR)/molecule.o
	$(CC) src/test/ring_test.cpp $(OBJS) -o test/ring_test $(CFLAGS)

test/pi_stack_test: src/test/pi_stack_test.cpp $(OBJS)
	$(CC) src/test/pi_stack_test.cpp $(OBJS) -o test/pi_stack_test $(CFLAGS)

test/aniso_test: src/test/aniso_test.cpp $(OBJS)
	$(CC) src/test/aniso_test.cpp $(OBJS) -o test/aniso_test $(CFLAGS)

test/group_test_mol: src/test/group_test_mol.cpp $(OBJS)
	$(CC) src/test/group_test_mol.cpp $(OBJS) -o test/group_test_mol $(CFLAGS)

test/group_test_res: src/test/group_test_res.cpp $(OBJS)
	$(CC) src/test/group_test_res.cpp $(OBJS) -o test/group_test_res $(CFLAGS)

test/moiety_test: src/test/moiety_test.cpp $(OBJS)
	$(CC) src/test/moiety_test.cpp $(OBJS) -o test/moiety_test $(CFLAGS)

test/mol_assem_test: src/test/mol_assem_test.cpp $(OBJS)
	$(CC) src/test/mol_assem_test.cpp $(OBJS) -o test/mol_assem_test $(CFLAGS)

test/eclipsing_test: src/test/eclipsing_test.cpp $(OBJS)
	$(CC) src/test/eclipsing_test.cpp $(OBJS) -o test/eclipsing_test $(CFLAGS)

test/amino_test: src/test/amino_test.cpp $(OBJS) $(OBJDIR)/aminoacid.o
	$(CC) src/test/amino_test.cpp $(OBJS) -o test/amino_test $(CFLAGS)

test/protein_test: src/test/protein_test.cpp $(OBJS) $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o
	$(CC) src/test/protein_test.cpp $(OBJS) -o test/protein_test $(CFLAGS)

test/backbone_test: src/test/backbone_test.cpp $(OBJS) $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o
	$(CC) src/test/backbone_test.cpp $(OBJS) -o test/backbone_test $(CFLAGS)

test/bond_rotation_test: src/test/bond_rotation_test.cpp $(OBJS) $(OBJDIR)/molecule.o
	$(CC) src/test/bond_rotation_test.cpp $(OBJS) -o test/bond_rotation_test $(CFLAGS)

test/flexion_test: src/test/flexion_test.cpp $(OBJS)
	$(CC) src/test/flexion_test.cpp $(OBJS) -o test/flexion_test $(CFLAGS)

test/histidine_test: src/test/histidine_test.cpp $(OBJS)
	$(CC) src/test/histidine_test.cpp $(OBJS) -o test/histidine_test $(CFLAGS)

test/cs_test: src/test/cs_test.cpp $(OBJS)
	$(CC) src/test/cs_test.cpp $(OBJS) -o test/cs_test $(CFLAGS)

test/mcoord_test: src/test/mcoord_test.cpp $(OBJS)
	$(CC) src/test/mcoord_test.cpp $(OBJS) -o test/mcoord_test $(CFLAGS)

$(BINDIR)/primarydock: src/primarydock.cpp $(OBJS) $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o $(OBJDIR)/group.o $(OBJDIR)/search.o $(OBJDIR)/scoring.o
	$(CC) src/primarydock.cpp $(OBJS) -o $(BINDIR)/primarydock $(CFLAGS)

$(BINDIR)/pepteditor: src/interpreter.cpp $(OBJS) $(OBJDIR)/aminoacid.o $(OBJDIR)/protein.o $(OBJDIR)/group.o
	$(CC) src/interpreter.cpp $(OBJS) -o $(BINDIR)/pepteditor $(CFLAGS)

$(BINDIR)/ic: src/ic.cpp $(OBJS) $(OBJDIR)/protein.o
	$(CC) src/ic.cpp $(OBJS) -o $(BINDIR)/ic $(CFLAGS)

$(BINDIR)/score_pdb: src/score_pdb.cpp $(OBJS) $(OBJDIR)/protein.o $(OBJDIR)/scoring.o
	$(CC) src/score_pdb.cpp $(OBJS) -o $(BINDIR)/score_pdb $(CFLAGS)

$(BINDIR)/ramachandran: src/ramachandran.cpp $(OBJS) $(OBJDIR)/protein.o
	$(CC) src/ramachandran.cpp $(OBJS) -o $(BINDIR)/ramachandran $(CFLAGS)

$(BINDIR)/ringflip: src/ringflip.cpp $(OBJS) $(OBJDIR)/molecule.o
	$(CC) src/ringflip.cpp $(OBJS) -o $(BINDIR)/ringflip $(CFLAGS)

performance_test: $(BINDIR)/primarydock testdata/test_TAAR8.config testdata/TAAR8.upright.pdb testdata/CAD_ion.sdf
	./$(BINDIR)/primarydock testdata/test_TAAR8.config

# low-tooling regression tests below
amino_report: REPORT="testdata/amino_test.approved.txt"
amino_report: test/amino_test
	./test/amino_test >testdata/received/amino_test.received.txt
	diff --color --unified $(REPORT) testdata/received/amino_test.received.txt

atom_report: REPORT="testdata/atom_test.approved.txt"
atom_report: test/atom_test
	./test/atom_test H >testdata/received/atom_test.received.txt
	diff --color --unified $(REPORT) testdata/received/atom_test.received.txt

aniso_report: REPORT="testdata/aniso_test.approved.txt"
aniso_report: test/aniso_test
	./test/aniso_test --asciiart >testdata/received/aniso_test.received.txt
	diff --color --unified $(REPORT) testdata/received/aniso_test.received.txt

point_report: REPORT="testdata/point_test.approved.txt"
point_report: test/point_test
	./test/point_test >testdata/received/point_test.received.txt
	diff --color --unified $(REPORT) testdata/received/point_test.received.txt

molecule_report: REPORT="testdata/molecule_test1.approved.txt"
molecule_report: test/molecule_test
	./test/molecule_test 'CCO' 'CCO' | tee temp | sed '/^#/d' >testdata/received/molecule_test1.received.txt; cat temp # ignore lines starting with #
	diff --color --unified $(REPORT) testdata/received/molecule_test1.received.txt

mol_assem_report: REPORT="testdata/mol_assem_test.approved.txt"
mol_assem_report: test/mol_assem_test
	bash ./test/mol_assem_tests.sh	# these must be checked visually.
	./test/mol_assem_test >testdata/received/molecule_test.received.txt
	echo "Content of test.sdf:" >> testdata/received/molecule_test.received.txt
	sed '2d' test.sdf >> testdata/received/molecule_test.received.txt  # remove line 2 (date stamp)
	diff --color --unified $(REPORT) testdata/received/molecule_test.received.txt

protein_report: REPORT="testdata/protein_test.approved.txt"
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

motif_report: REPORT="testdata/motif_test.approved.txt"
motif_report: bin/pepteditor test/motif_test.pepd
	./bin/pepteditor test/motif_test.pepd > testdata/received/motif_test.received.txt
	diff --color --unified $(REPORT) testdata/received/motif_test.received.txt

reports: $(REPORTS)
