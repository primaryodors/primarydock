
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "protein.h"

using namespace std;

int xition_no = 1;

// This is used only for debugging. Once the utility is perfected, it can be removed.
void save_transitional_pdb(Protein* p)
{	char filename[64];
	sprintf(filename, "temp.%d.pdb", xition_no++);
	FILE* pf = fopen(filename, "wb");
	p->save_pdb(pf);
	p->end_pdb(pf);
	fclose(pf);
}

int main(int argc, char** argv)
{
	if (argc < 2)
	{	cout << "No input file." << endl;
		return -1;
	}
	
	char* inpfile = NULL;
	char* outfile = NULL;
	int mcoordres[8]{};
	char* mcoorda[8]{};
	char* esym = NULL;
	int i, j, charge=2;
	Point pocketcen;
	bool pocketset = false;
	bool dohelix = true;
	
	j=0;
	for (i=1; i<argc; i++)
	{	if 		(!strcmp(argv[i], "-p")) inpfile = argv[++i];
		else if (!strcmp(argv[i], "-r"))
		{	mcoordres[j] = atoi(argv[++i]);
			mcoorda[j++] = strstr(argv[i], ":")+1;
		}
		else if (!strcmp(argv[i], "-o")) outfile = argv[++i];
		/*else if (!strcmp(argv[i], "-b"))
		{	pocketcen.x = atof(argv[++i]);
			pocketcen.y = atof(argv[++i]);
			pocketcen.z = atof(argv[++i]);
			pocketset = true;
		}*/
		else if (!strcmp(argv[i], "-c")) charge = atoi(argv[++i]);
		// else if (!strcmp(argv[i], "--nh")) dohelix = false;
		else if (!strcmp(argv[i], "-e")) esym = argv[++i];
	}
	
	int numres = j;
	if (!inpfile || !numres || !mcoordres[0] || (mcoorda[0] < argv[1]) || !esym || !outfile)
	{	cout << "Usage:" << endl;
		cout << "metal -p {path/to/input.pdb} -r {resno}:{aname} -r {resno}:{aname} -r {resno}:{aname} -e {elem_sym} [-b {pocket.x} {pocket.y} {pocket.z}] -o {path/to/output.pdb}" << endl;
		cout << endl << "Example:" << endl;
		cout << "metal -p ";
		cout << (inpfile ? inpfile : "../pdbs/OR1A1/OR1A1.rotated.pdb");
		cout << " -r 176:OD1 -r 179:SG -r 180:OD1 -o output/OR1A1.metal.pdb -b -4.93883 5.08067 -4.91533 -e Zn" << endl;
		cout << endl << "Optionally, you may indicate the metal's charge, e.g. -c 1 for a charge of +1. The default charge, if -c is omitted, is +2." << endl;
		cout << endl << "By default, a helix will be created and moved above the binding pocket. If this is not a helix-bound metal ion,"
			 << " or if its helix already exists in the protein, use the --nh option to leave the protein structure intact." << endl;
		cout << endl;
		return -1;
	}
	
	Protein p(inpfile);
	FILE* pf = fopen(inpfile, "r");
	if (!pf)
	{	cout << "Error trying to read " << inpfile << endl;
		return 0xbadf12e;
	}
	p.load_pdb(pf);
	fclose(pf);
	
	cout << "TMR1 angle: " << (p.get_helix_orientation( 21,  49)*fiftyseven-360) << endl;
	cout << "TMR2 angle: " << (p.get_helix_orientation( 58,  81)*fiftyseven) << endl;
	cout << "TMR3 angle: " << (p.get_helix_orientation( 94, 123)*fiftyseven-360) << endl;
	cout << "TMR4 angle: " << (p.get_helix_orientation(143, 163)*fiftyseven) << endl;
	cout << "TMR5 angle: " << (p.get_helix_orientation(196, 218)*fiftyseven-360) << endl;
	cout << "TMR6 angle: " << (p.get_helix_orientation(236, 261)*fiftyseven) << endl;
	cout << "TMR7 angle: " << (p.get_helix_orientation(268, 289)*fiftyseven-360) << endl;
	
	int startres = mcoordres[0]-1, endres = mcoordres[numres-1]+1, stopat = mcoordres[numres-1]+15;
	Point mtgt = pocketcen, wayuphigh = pocketcen;
	mtgt.y = 5000;
	wayuphigh.y += 4000;
	AminoAcid* aa = p.get_residue(stopat);
	Atom *Cend = aa->get_atom("C"), *Oend = aa->get_atom("O");
	Point Cpt = Cend->get_location(), Opt = Oend->get_location();
	if (dohelix)
	{	cout << "Extending strand." << endl;
		p.conform_backbone(startres-15, stopat, Cend, wayuphigh, 25);
		
		save_transitional_pdb(&p);
		
		cout << "Generating an alpha helix." << endl;
		p.make_helix(startres, endres, stopat, ALPHA_PHI, ALPHA_PSI);
		save_transitional_pdb(&p);
	}
	
	cout << "Coordinating metal ion." << endl;
	Atom* m = new Atom(esym, &mtgt);
	m->name = esym;
	strcpy(m->aa3let, "MTL");
	m->increment_charge(charge);
	p.coordinate_metal(m, numres, mcoordres, mcoorda);
	save_transitional_pdb(&p);
	mtgt.y = pocketcen.y + 15;
	
	if (dohelix)
	{
		/*cout << "Extending strand." << endl;
		p.conform_backbone(mcoordres[0]-12, stopat, m, wayuphigh, 25);
		save_transitional_pdb(&p);
		cout << "Extending strand further." << endl;
		p.conform_backbone(endres, stopat, Cend, wayuphigh, 25);
		save_transitional_pdb(&p);*/
		cout << "Bringing the metal over the pocket center." << endl;
		p.conform_backbone(mcoordres[0]-15, stopat, m, mtgt, 50);
		save_transitional_pdb(&p);
		cout << "Leveling the helix." << endl;
		cout << "Old helix angle: " << (p.get_helix_orientation(startres, endres)*fiftyseven) << endl;
		cout << "New helix angle: " << (p.orient_helix(startres, endres, stopat, 0, 15)*fiftyseven) << endl;
		save_transitional_pdb(&p);
		return 0;
		cout << "Extending post-helix strand." << endl;
		p.conform_backbone(endres, stopat, Cend, wayuphigh, 25);
		save_transitional_pdb(&p);
		
		cout << "Reconnecting broken pieces." << endl;
		p.conform_backbone(endres, stopat, Cend, Cpt, Oend, Opt, 50);
		save_transitional_pdb(&p);
		
		// Reconform the head segment.
		// p.conform_backbone(20, 1, 50, 15);
		save_transitional_pdb(&p);
	}
	
	pf = fopen(outfile, "wb");
	
	p.save_pdb(pf);
	p.end_pdb(pf);
	fclose(pf);
	
	return 0;
}
