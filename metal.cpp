
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
	char* rcpfile = NULL;
	char* rcpinfo = NULL;
	
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
		else if (!strcmp(argv[i], "-i"))
		{	if (argc < i+2 || !argv[i+1]) throw 0xbadf12e;
			rcpfile = argv[i+1];
		}
	}
	
	int numres = j;
	if (!inpfile 
		||	(	(!numres || !mcoordres[0] || (mcoorda[0] < argv[1]) || !esym )
				&& !rcpfile
			)
		|| !outfile)
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
	
	if (rcpfile)
	{	pf = fopen(rcpfile, "rb");
		if (!pf) throw 0xbadf12e;
		while (!feof(pf))
		{	char buffer[16384];
			fgets(buffer, 16361, pf);
			if (buffer[0] >= ' ')
			{	char** fields = chop_spaced_fields(buffer, '|');
				
				if (!strcmp(fields[0], "REGION"))
				{	p.set_region(fields[1], atoi(fields[2]), atoi(fields[3]));
				}
				else if (!strcmp(fields[0], "METAL"))
				{	esym = new char[strlen(fields[1])+4];
					strcpy(esym, fields[1]);
					
					for (i=2; fields[i]; i++)
					{	if (strstr(fields[i], ":"))
						{	mcoordres[j] = atoi(&fields[i][3]);
							mcoorda[j++] = strstr(fields[i], ":")+1;
						}
						else if (strstr(fields[i], "HELIX="))
						{	dohelix = (fields[i][6] - '0');
						}
					}
				}
				else if (!strcmp(fields[0], "POCKET"))
				{	if (!strcmp(fields[1], "type=A"))
					{	char* pktf = strstr(fields[2], "=")+1;
						pocketcen.x = atof(pktf);
						pktf = strstr(fields[3], "=")+1;
						pocketcen.y = atof(pktf);
						pktf = strstr(fields[4], "=")+1;
						pocketcen.z = atof(pktf);
					}
				}
			}
		}
		fclose(pf);
		numres = j;
		
		for (i=1; i<=7; i++)
		{	char rgname[10]{};
			sprintf(rgname, "TMR%d", i);	
			Region rgn = p.get_region(rgname);
			float turn = (i & 1) ? -360 : 0;
			if (rgn.start) cout << rgname << " angle: " << (p.get_helix_orientation(rgn.start, rgn.end)*fiftyseven+turn) << endl;
		}
	}
	
	
	
	
	int startres, endres, stopat;
	
	Point mtgt = pocketcen, wayuphigh = pocketcen;
	mtgt.y = 5000;
	wayuphigh.y += 4000;
	AminoAcid* aa = p.get_residue(stopat);
	Atom *Cend = aa->get_atom("C"), *Oend = aa->get_atom("O");
	Point Cpt = Cend->get_location(), Opt = Oend->get_location();
	if (dohelix)
	{	startres = mcoordres[0]-1;
		endres = mcoordres[numres-1]+1;
		Region rgn = p.get_region("TMR5");
		stopat = rgn.start-1;
		
		/*cout << "Extending strand." << endl;
		p.conform_backbone(startres-15, stopat, Cend, wayuphigh, 25);
		
		save_transitional_pdb(&p);
		*/
		
		cout << "Generating an alpha helix." << endl;
		p.make_helix(startres, endres, stopat, ALPHA_PHI, ALPHA_PSI);
		save_transitional_pdb(&p);
		
		// Move the start of the helix to above TMR5.
		
		// Move the end of the helix to about 2/3 of the way to a point above TMR3.
		
		// Find a point along the region between the helix and TMR5 and move it between the helix and pocket center.
		
		// Just from the point found in the previous step, reconnect the loose end to TMR5.
		
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
		save_transitional_pdb(&p);
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
		save_transitional_pdb(&p);*/
		
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
