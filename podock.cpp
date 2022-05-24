#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include "protein.h"

using namespace std;

#define _DBG_STEPBYSTEP true
#define _DORESPHRES false

struct DockResult
{	int pose;
	float kJmol;
	char** metric;
	float* mkJmol;
	std::string pdbdat;
	float bytype[_INTER_TYPES_LIMIT];
};

char* get_file_ext(char* filename)
{	int i = strlen(filename)-1;
	for (; i>=0; i--)
	{	if (filename[i] == '.') return &filename[i+1];
	}
	return 0;
}

Protein* protein;
int seql = 0;
int mcoord_resno[256] = {};
Molecule* ligand;
Point ligcen_target;
Point size(10,10,10);
Vector path[256] = {};
int pathnodes=1;		// The pocketcen is the initial node.
int poses = 10;
int iters = 50;
bool flex=true;
float kJmol_cutoff = 0.01;
float drift = 0.333;
Molecule** gcfmols = NULL;

void iteration_callback(int iter)
{
	// if (kJmol_cutoff > 0 && ligand->lastbind >= kJmol_cutoff) iter = (iters-1);
	if (iter == (iters-1)) return;
	
	Point bary = ligand->get_barycenter();
					
	if (bary.get_3d_distance(ligcen_target) > size.magnitude())
	{	//cout << "Wrangle! " << bary << ": " << bary.get_3d_distance(ligcen_target) << " vs. " << size.magnitude() << endl;
		bary = ligcen_target;
		// ligand->reset_conformer_momenta();
	}
	else
	{	bary.x += (ligcen_target.x - bary.x) * drift;
		bary.y += (ligcen_target.y - bary.y) * drift;
		bary.z += (ligcen_target.z - bary.z) * drift;
	}
	
	ligand->recenter(bary);
	
	drift *= (1.0 - 0.5/iters);
	
	if (gcfmols && seql)
	{	Star discrete[SPHREACH_MAX+4];
		discrete[0].pmol = gcfmols[0];
		discrete[1].pmol = gcfmols[1];
		
		int i;
		AminoAcid* resphres[SPHREACH_MAX+4] = {};
		int sphres = protein->get_residues_can_collide_ligand(resphres, ligand, bary, size, mcoord_resno);
		for (i=0; i<sphres; i++)
		{	discrete[i+2].paa = resphres[i];
		}
		discrete[sphres+2].n = 0;
		
		sphres += 2;
		for (i=0; i<sphres; i++) gcfmols[i] = discrete[i].pmol;
		gcfmols[sphres] = NULL;
	}
}

int main(int argc, char** argv)
{	char configfname[256] = {};
	char protfname[256] = {};
	char ligfname[256] = {};
	Point pocketcen;
	std::ofstream *output = NULL;
	
	bool configset=false, protset=false, ligset=false, pktset=false;
	
	time_t began = time(NULL);
	
	strcpy(configfname, "podock.config");
	
	if (argc > 1)
	{	strcpy(configfname, argv[1]);
		configset = true;
	}
	
	FILE* pf = fopen(configfname, "r");
	if (!pf)
	{	cout << "Config file not found: " << configfname << ", exiting." << endl;
		return 0xbadf12e;
	}
	
	while (!feof(pf))
	{	char buffer[1024];
		fgets(buffer, 1015, pf);
		if (buffer[0] >= ' ' && buffer[0] != '#')
		{	char** fields = chop_spaced_fields(buffer);
			
			if (!strcmp(fields[0], "PROT"))
			{	strcpy(protfname, fields[1]);
				protset = true;
			}
			else if (!strcmp(fields[0], "LIG"))
			{	strcpy(ligfname, fields[1]);
				ligset = true;
			}
			else if (!strcmp(fields[0], "CEN"))
			{	pocketcen.x = atof(fields[1]);
				pocketcen.y = atof(fields[2]);
				pocketcen.z = atof(fields[3]);
				pktset = true;
			}
			else if (!strcmp(fields[0], "PATH"))
			{	int nodeno = atoi(fields[1]);
				if (nodeno > 255)
				{	cout << "Binding path is limited to 255 nodes." << endl;
					return 0xbad90de;
				}
				if (nodeno)
				{	Point pt(atof(fields[2]),
							 atof(fields[3]),
							 atof(fields[4])
							);
					Vector v(&pt);
					path[nodeno] = v;
					if ((nodeno) > pathnodes) pathnodes = nodeno;
				}
			}
			else if (!strcmp(fields[0], "SIZE"))
			{	size.x = atof(fields[1]);
				if (fields[2])
				{	size.y = atof(fields[2]);
					size.z = atof(fields[3]);
				}
				else size.z = size.y = size.x;
				if (!size.x || !size.y || !size.z)
				{	cout << "Pocket size cannot be zero in any dimension." << endl;
					return 0xbad512e;
				}
			}
			else if (!strcmp(fields[0], "POSE"))
			{	poses = atoi(fields[1]);
			}
			else if (!strcmp(fields[0], "EMIN"))
			{	kJmol_cutoff = atof(fields[1]);
			}
			else if (!strcmp(fields[0], "FLEX"))
			{	flex = (atoi(fields[1]) != 0);
			}
			else if (!strcmp(fields[0], "ITER") || !strcmp(fields[0], "ITERS"))
			{	iters = atoi(fields[1]);
			}
			else if (!strcmp(fields[0], "MCOORD"))
			{	int i, j=0;
				for (i=1; fields[i]; i++)
				{	mcoord_resno[j++] = atoi(fields[i]);
				}
				mcoord_resno[j] = 0;
			}
			else if (!strcmp(fields[0], "DEBUG"))
			{	if (!fields[1])
				{	cout << "Missing debug file name; check config file." << endl;
					throw 0xbadf12e;
				}
				#if _DBG_STEPBYSTEP
				cout << "Starting a debug outstream." << endl;
				#endif
				debug = new std::ofstream(fields[1], std::ofstream::out);
			}
			else if (!strcmp(fields[0], "OUT"))
			{	if (!fields[1])
				{	cout << "Missing output file name; check config file." << endl;
					throw 0xbadf12e;
				}
				#if _DBG_STEPBYSTEP
				cout << "Starting a file outstream." << endl;
				#endif
				output = new std::ofstream(fields[1], std::ofstream::out);
			}
			
			delete fields;
		}
	}
	fclose(pf);
	#if _DBG_STEPBYSTEP
	if (debug) *debug << "Loaded config file." << endl;
	#endif
	
	drift = 1.0 / (iters/25+1);
	
	// Load the protein or return an error.
	Protein p(protfname);
	protein = &p;
	pf = fopen(protfname, "r");
	if (!pf)
	{	cout << "Error trying to read " << protfname << endl;
		return 0xbadf12e;
	}
	p.load_pdb(pf);
	fclose(pf);
	#if _DBG_STEPBYSTEP
	if (debug) *debug << "Loaded protein." << endl;
	#endif
	
	p.mcoord_resnos = mcoord_resno;
	
	// Load the ligand or return an error.
	Molecule m(ligfname);
	ligand = &m;
	char* ext = get_file_ext(ligfname);
	if (!ext)
	{	cout << "Ligand file is missing its extension! " << ligfname << endl;
		return 0xbadf12e;
	}
	
	char buffer[65536] = {};
	switch (ext[0])
	{	case 's':
		case 'S':
		// SDF
		pf = fopen(ligfname, "r");
		fread(buffer, 1, 65535, pf);
		fclose(pf);
		m.from_sdf(buffer);
		break;
		
		case 'p':
		case 'P':
		pf = fopen(ligfname, "r");
		m.from_pdb(pf);
		fclose(pf);
		break;
		
		default:
		cout << "Unrecognized ligand file extension: " << ext << endl;
		return 0xbadf12e;
	}
	#if _DBG_STEPBYSTEP
	if (debug) *debug << "Loaded ligand." << endl;
	#endif
	
	Point box = m.get_bounding_box();
	
	if (debug) *debug << "Ligand bounding box corner (centered at zero): " << box.printable() << endl;
	#if _DBG_STEPBYSTEP
	if (debug) *debug << "Ligand bounding box." << endl;
	#endif
	
	// Identify the ligand atom with the greatest potential binding.
	Atom* ligbb = m.get_most_bindable();
	Atom* ligbbh = NULL;
	intera_type lig_inter_typ = vdW;
	
	if (ligbb->get_Z() == 1)
	{	ligbbh = ligbb;
		ligbb = ligbbh->get_bond_by_idx(0)->btom;
	}	else
	{	ligbbh = ligbb->is_bonded_to("H");
	}
	
	if (fabs(ligbb->get_charge()) >= 1
		||
		ligbb->get_acidbase()
		||
		(ligbbh && fabs(ligbbh->get_charge()) >= 1)
		||
		(ligbbh && ligbbh->get_acidbase())
	   )
		lig_inter_typ = ionic;
	else if (fabs(ligbb->is_polar()) >= 1
			 ||
			 (ligbbh && fabs(ligbbh->is_polar()) >= 1)
			)
		lig_inter_typ = hbond;
	else if (ligbb->is_pi())
		lig_inter_typ = pi;

	#if _DBG_STEPBYSTEP
	if (debug) *debug << "Identified best binding ligand atoms." << endl;
	#endif
	
	int pose, nodeno, iter;
	int i, j, k, l, n;
	Point nodecen = pocketcen;
	seql = p.get_seq_length();
	int rstart = p.get_start_resno();
	AminoAcid* reaches_spheroid[pathnodes+2][SPHREACH_MAX] = {};
	int sphres = 0;
	
	// Filter residues according to which ones are close enough to the spheroid to "reach" it.
	nodecen = pocketcen;
	
	// When docking with a metalloprotein, use this temporary Molecule for interactions the same as
	// we use AminoAcid objects, except don't attempt to flex the metals object.
	Molecule* met = p.metals_as_molecule();
	#if _DBG_STEPBYSTEP
	if (debug) *debug << "Created metals molecule." << endl;
	#endif
	
	float bcoll = 0;
	
	for (l=0; l<seql; l++)
	{	AminoAcid* laa = p.get_residue(l+rstart);
		if (!laa) continue;
		for (n=l+1; n<seql; n++)
		{	AminoAcid* naa = p.get_residue(n+rstart);
			if (!naa) continue;
			bcoll += laa->get_intermol_collisions(naa);
		}
		
		bcoll += m.get_intermol_collisions(laa);
		
		if (met) bcoll += laa->get_intermol_collisions(met);
	}
	if (met) bcoll += m.get_intermol_collisions(met);
	if (debug) *debug << "Initial collisions: " << bcoll << endl;
	
	
	// TODO: Output some basic stats: receptor, ligand, etc.
	
	if (ligbb) 
	{	cout << "# Best binding heavy atom of ligand" << endl << "LBBA: " << ligbb->name << endl;
		if (output) *output << "# Best binding heavy atom of ligand" << endl << "LBBA: " << ligbb->name << endl;
	}
	if (ligbbh)
	{	cout << "# Best binding hydrogen of ligand" << endl << "LBBH: " << ligbbh->name << endl;
		if (output) *output << "# Best binding hydrogen of ligand" << endl << "LBBH: " << ligbbh->name << endl;
	}
	
	DockResult dr[poses+2][pathnodes+2];
	for (i=0; i<poses; i++) dr[i][0].kJmol = 0;
	int drcount = 0;
	
	srand(0xb0ad1cea);
	// srand(time(NULL));
	for (pose=1; pose<=poses; pose++)
	{	
		#if _DBG_STEPBYSTEP
		if (debug) *debug << "Pose " << pose << endl;
		#endif
		nodecen = pocketcen;
		nodecen.weight = 1;
		
		for (nodeno=0; nodeno<pathnodes; nodeno++)
		{	
			#if _DBG_STEPBYSTEP
			if (debug) *debug << "Pose " << pose << endl << "Node " << nodeno << endl;
			#endif
			if (nodeno)
			{	nodecen = nodecen.add(&path[nodeno]);
				#if _DBG_STEPBYSTEP
				if (debug) *debug << "Added whatever points together." << endl;
				#endif
			}
			Point lastnodecen = nodecen;
			ligcen_target = nodecen;
			
			#if _DBG_STEPBYSTEP
			if (debug) *debug << "Saved last nodecen." << endl;
			#endif
			
			// Move the ligand to the new node center.
			if (!nodeno) m.recenter(nodecen);
			#if _DBG_STEPBYSTEP
			if (debug) *debug << "Molecule recenter (or not)." << endl;
			#endif
			m.reset_conformer_momenta();
			#if _DBG_STEPBYSTEP
			if (debug) *debug << "Conformer momenta reset." << endl;
			#endif
			
			
			sphres = p.get_residues_can_collide_ligand(reaches_spheroid[nodeno], &m, nodecen, size, mcoord_resno);
			for (i=sphres; i<SPHREACH_MAX; i++) reaches_spheroid[nodeno][i] = NULL;
			
			/*cout << "Dock residues for node " << nodeno << ": " << endl;
			if (output) *output << "Dock residues for node " << nodeno << ": " << endl;
			for (i=0; i<sphres; i++)
			{	cout << *reaches_spheroid[nodeno][i] << " ";
				if (output) *output << *reaches_spheroid[nodeno][i] << " ";
			}
			cout << endl << endl;
			if (output) *output << endl << endl;*/
			
			
			if (!nodeno)
			{	Molecule* alignment_aa=0;
				#if _DBG_STEPBYSTEP
				if (debug) *debug << "Initialize null AA pointer." << endl;
				#endif
				
				// Rotate the ligand in space so that its strongest binding atom points to a binding
				// pocket feature with a strong potential binding.
				for (i=0; reaches_spheroid[nodeno][i]; i++)
				{	if (!p.aa_ptr_in_range(reaches_spheroid[nodeno][i]))
					{	reaches_spheroid[nodeno][i] = NULL;
						continue;
					}
				
					#if _DBG_STEPBYSTEP
					if (debug)
					{	*debug << "Check capable of inter (" << i << ") ";
						*debug << lig_inter_typ;
						*debug << flush;
						Star s;
						s.paa = reaches_spheroid[nodeno][i];
						*debug << std::hex << s.n << std::dec << " " << flush;
						*debug << *reaches_spheroid[nodeno][i];
						*debug << endl;
					}
					#endif
					if (reaches_spheroid[nodeno][i]->capable_of_inter(lig_inter_typ)
						&&
						(	!alignment_aa 
							||
							(!(rand() % sphres))
						)
					   )
					alignment_aa = reaches_spheroid[nodeno][i];
					#if _DBG_STEPBYSTEP
					if (debug) *debug << "Candidate alignment AA." << endl;
					#endif
				}
				#if _DBG_STEPBYSTEP
				if (debug) *debug << "Selected an alignment AA." << endl;
				#endif
				
				if (met) alignment_aa = met;
				#if _DBG_STEPBYSTEP
				if (debug) *debug << "Alignment AA." << endl;
				#endif
				
				if (alignment_aa)
				{	Atom* alca;
					if (alignment_aa == met)
						alca = alignment_aa->get_nearest_atom(ligbb->get_location());
					else
						alca = alignment_aa->get_atom("CA");
					#if _DBG_STEPBYSTEP
					if (debug) *debug << "Got alignment atom." << endl;
					#endif
					
					if (alca)
					{	Point pt, al, cen;
						pt	= ligbb->get_location();
						al	= alca->get_location();
						cen	= m.get_barycenter();
						
						Rotation rot = align_points_3d(&pt, &al, &cen);
						m.rotate(&rot.v, rot.a);
						
						// Preemptively minimize intermol collisions.
						m.intermol_conform_norecen(alignment_aa, iters, reaches_spheroid[nodeno]);
						if (debug) *debug << "Alignment atom is " << alignment_aa->get_name() << ":" << alca->name << " Z " << alca->get_Z() << endl;
					}
				}
				#if _DBG_STEPBYSTEP
				if (debug) *debug << "Aligned ligand to AA." << endl;
				#endif
			}
			
			// float driftamt = 1.0 / (iters/25+1);
			// cout << pose << ":" << nodeno << " drift " << driftamt << endl;
			int iters_div = iters*0.259;
			
			Molecule* cfmols[SPHREACH_MAX+4] = {};
			gcfmols = cfmols;
			i=0;
			m.movability = MOV_ALL;
			cfmols[i++] = &m;
			if (met)
			{	met->movability = MOV_NONE;
				cfmols[i++] = met;
			}
			for (j=0; j<sphres; j++)
			{	if (reaches_spheroid[nodeno][j]->movability >= MOV_FLEXONLY) reaches_spheroid[nodeno][j]->movability = MOV_FLEXONLY;
				cfmols[i++] = reaches_spheroid[nodeno][j];
			}
			for (; i<SPHREACH_MAX; i++)
				cfmols[i] = NULL;
			
			// time_t preiter = time(NULL);
			Molecule::multimol_conform(
				cfmols,
				iters,
				&iteration_callback
				);
			/*time_t jlgsux = time(NULL);
			cout << "\nIterations took: " << (jlgsux-preiter) << " seconds." << endl;*/
			
			// Add the current pose/path sequentially to the dr[][] array.
			// If the path node # is zero:
			// If it is the first (zeroth) entry, set the pose number to 1.
			// Otherwise, go through all of the preceding entries and:
			// Any entry with a smaller kJ/mol, increment its pose# but remember the smallest pre-increment pose #
			// from the lot of them;
			// Claim that new smallest pose# (which might be 1) as your own.
			#if _DBG_STEPBYSTEP
			if (debug) *debug << "Preparing output." << endl;
			#endif
			
			char metrics[p.get_seq_length()+8][10] = {};
			float mkJmol[p.get_seq_length()+8] = {};
			int metcount=0;
			float btot=0;
			
			for (i=0; i<_INTER_TYPES_LIMIT; i++) total_binding_by_type[i] = 0;
			
			if (debug) *debug << "Pose " << pose << " pathnode " << nodeno /*<< " collisions " << coll*/ << endl;
			
			m.clear_atom_binding_energies();
			
			drcount = pose-1;
			if (met)
			{	met->clear_atom_binding_energies();
				float lb = m.get_intermol_binding(met);
				strcpy(metrics[metcount], "Metals");
				mkJmol[metcount++] = lb;
				btot += lb;
			}
			
			sphres = p.get_residues_can_collide_ligand(reaches_spheroid[nodeno], &m, m.get_barycenter(), size, mcoord_resno);
			// cout << "sphres " << sphres << endl;
			for (i=0; i<sphres; i++)
			{	if (!reaches_spheroid[nodeno][i]) continue;
				if (!p.aa_ptr_in_range(reaches_spheroid[nodeno][i])) continue;
				reaches_spheroid[nodeno][i]->clear_atom_binding_energies();
				float lb = m.get_intermol_binding(reaches_spheroid[nodeno][i]);
				if (lb > 90) lb = 0;
				sprintf(metrics[metcount], "%s%d", reaches_spheroid[nodeno][i]->get_3letter(), reaches_spheroid[nodeno][i]->get_residue_no());
				// cout << metrics[metcount] << ": " << lb << " . ";
				mkJmol[metcount++] = lb;
				btot += lb;
			}
			// cout << btot << endl;
			
			if (btot > 5*m.get_atom_count()) btot = 0;
			
			#if _DBG_STEPBYSTEP
			if (debug) *debug << "Prepared metrics." << endl;
			#endif
			
			dr[drcount][nodeno].kJmol = btot;
			dr[drcount][nodeno].metric = new char*[metcount+2]{};
			dr[drcount][nodeno].mkJmol = new float[metcount]{};
			#if _DBG_STEPBYSTEP
			if (debug) *debug << "Allocated memory." << endl;
			#endif
			
			for (i=0; i<_INTER_TYPES_LIMIT; i++)
			{	dr[drcount][nodeno].bytype[i] = total_binding_by_type[i];
			}
			#if _DBG_STEPBYSTEP
			if (debug) *debug << "Filled btypes." << endl;
			#endif
			
			for (i=0; i<metcount; i++)
			{	dr[drcount][nodeno].metric[i] = new char[max(8,(int)strlen(metrics[i])+4)]{};
				strcpy(dr[drcount][nodeno].metric[i], metrics[i]);
				dr[drcount][nodeno].mkJmol[i] = mkJmol[i];
				// cout << "*" << dr[drcount][nodeno].metric[i] << ": " << dr[drcount][nodeno].mkJmol[i] << endl;
			}
			dr[drcount][nodeno].metric[i] = new char[1]{};
			dr[drcount][nodeno].metric[i][0] = 0;
			#if _DBG_STEPBYSTEP
			if (debug) *debug << "More metrics or something idfk." << endl;
			#endif
			
			std::ostringstream pdbdat;
			
			// Prepare a partial PDB of the ligand atoms and all involved residue sidechains.
			n = m.get_atom_count();
			for (l=0; l<n; l++) m.get_atom(l)->stream_pdb_line(pdbdat, 9000+l);
			#if _DBG_STEPBYSTEP
			if (debug) *debug << "Prepared ligand PDB." << endl;
			#endif
			
			if (flex)
			{	for (k=0; reaches_spheroid[nodeno][k]; k++)
				{	if (!p.aa_ptr_in_range(reaches_spheroid[nodeno][k])) continue;
					n = reaches_spheroid[nodeno][k]->get_atom_count();
					for (l=0; l<n; l++)
					{	reaches_spheroid[nodeno][k]->get_atom(l)->stream_pdb_line(
							pdbdat,
							reaches_spheroid[nodeno][k]->atno_offset+l
						);
					}
				}
				#if _DBG_STEPBYSTEP
				if (debug) *debug << "Prepared flex PDBs." << endl;
				#endif
			}
			
			dr[drcount][nodeno].pdbdat = pdbdat.str();
			if (debug) *debug << "Prepared the PDB strings." << endl;
			
			if (!nodeno)
			{	if (pose==1) dr[drcount][nodeno].pose = pose;
				else
				{	int bestpose = pose;
					for (i=0; i<drcount; i++)
					{	if (dr[i][0].kJmol < btot)
						{	if (dr[i][0].pose < bestpose || bestpose < 0) bestpose = dr[i][0].pose;
							dr[i][0].pose++;
						}
					}
					dr[drcount][nodeno].pose = bestpose;
					// cout << "Around the posie: "; for (i=0; i<=drcount; i++) cout << dr[i][nodeno].pose << " "; cout << endl;
				}
				#if _DBG_STEPBYSTEP
				if (debug) *debug << "Added pose to output array." << endl;
				#endif
			}
			
			drcount = pose;
		}
	}
	#if _DBG_STEPBYSTEP
	if (debug) *debug << "Finished poses." << endl;
	#endif
	
	// Output the dr[][] array in order of increasing pose number.
	cout << endl;
	if (output) *output << endl;
	
	for (i=1; i<=poses; i++)
	{	for (j=0; j<poses; j++)
		{	if (dr[j][0].pose == i)
			{	if (dr[j][0].kJmol >= kJmol_cutoff)
				{	for (k=0; k<pathnodes; k++)
					{	cout << "Pose: " << i << endl << "Node: " << k << endl;
						if (output) *output << "Pose: " << i << endl << "Node: " << k << endl;
						
						cout << "# Binding energies" << endl << "BENERG:" << endl;
						if (output) *output << "# Binding energies" << endl << "BENERG:" << endl;
						for (l=0;
							 dr[j][k].metric
							 	&& dr[j][k].metric[l]
							 	&& dr[j][k].metric[l][0];
							 l++)
						{	cout << dr[j][k].metric[l] << ": " << dr[j][k].mkJmol[l] << endl;
							if (output && dr[j][k].metric[l]) *output << dr[j][k].metric[l] << ": " << dr[j][k].mkJmol[l] << endl;
						}
						
						for (l=0; l<_INTER_TYPES_LIMIT; l++)
						{	char lbtyp[64] = {};
							switch (l+covalent)
							{	case covalent:		continue; /*strcpy(lbtyp, "Total covalent: ");		break;*/
								case ionic:			strcpy(lbtyp, "Total ionic: ");						break;
								case hbond:			strcpy(lbtyp, "Total H-bond: ");					break;
								case pi:			strcpy(lbtyp, "Total pi stack: ");					break;
								case polarpi:		strcpy(lbtyp, "Total polar-pi and cation-pi: ");	break;
								case mcoord:		strcpy(lbtyp, "Total metal coordination: ");		break;
								case vdW:			strcpy(lbtyp, "Total van der Waals: ");				break;
								default:			goto _btyp_unassigned;
							}
							cout << lbtyp << dr[j][k].bytype[l] << endl;
							if (output) *output << lbtyp << dr[j][k].bytype[l] << endl;
						}
						_btyp_unassigned:
						
						if (output) *output << "Total: " << dr[j][k].kJmol << endl << endl;
						cout << "Total: " << dr[j][k].kJmol << endl << endl;
						
						if (!dr[j][k].pdbdat.length())
						{	cout << "Uh-oh!" << endl;
							if (output) *output << "(Missing PDB data.)" << endl;
						}
						else
						{	cout << "# PDB Data" << endl << "PDBDAT:" << endl;
							if (output) *output << "# PDB Data" << endl << "PDBDAT:" << endl;
						
							if (output) *output << dr[j][k].pdbdat << endl;
							cout << dr[j][k].pdbdat << endl;
							
							cout << "TER" << endl << "END" << endl << endl << endl;
							if (output) *output << "TER" << endl << "END" << endl << endl << endl;
						}
					}
				}	else
				{	if (i == 1)
					{	cout << "No poses found within kJ/mol limit." << endl;
						if (output) *output << "No poses found within kJ/mol limit." << endl;
					}
					cout << "Exiting." << endl;
					goto _exitposes;
				}
				
				break;
			}
		}
	}
	_exitposes:
	cout << (i-1) << " pose(s) found." << endl;
	if (output) *output << (i-1) << " pose(s) found." << endl;
	if (debug) *debug << (i-1) << " pose(s) found." << endl;
	
	if (met) delete met;
	
	time_t finished = time(NULL);
	cout << "\nCalculation took: " << (finished-began) << " seconds." << endl;
	if (output) *output << "\nCalculation took: " << (finished-began) << " seconds." << endl;
	if (debug) *debug << "\nCalculation took: " << (finished-began) << " seconds." << endl;
	
	if (debug) debug->close();
	if (output) output->close();
	
	return 0;
}


















