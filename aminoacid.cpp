
#include <math.h>
#include <stdio.h>
#include <iostream> 
#include <string.h>
#include <algorithm>
#include <ctime>
#include <stdint.h>
#include <cassert>
#include "molecule.h"
#include "aminoacid.h"

using namespace std;

AADef aa_defs[26];
char* override_aminos_dat=0;

AminoAcid::AminoAcid(FILE* instream)
{	if (!aa_defs[0]._1let) AminoAcid::load_aa_defs();
	immobile = true;
	movability = MOV_FLEXONLY;
	from_pdb(instream);
	mincoll = get_internal_collisions();
}

AminoAcid::AminoAcid(const char letter, AminoAcid* prevaa)
{	if (!aa_defs[0]._1let) AminoAcid::load_aa_defs();
	immobile = true;
	movability = MOV_FLEXONLY;

	if (!prevaa) residue_no = 1;
	else residue_no = prevaa->residue_no + 1;
	
	int idx = (letter & 0x5f) - 'A';
	aadef = &aa_defs[idx];
	if (!aa_defs[idx].aabonds)
	{	cout << "Cannot load " << letter << " please make sure amino acid exists in aminos.dat and all atoms are defined." << endl;
		throw 0xbadac1d;			// Also known as p!ɔepeq.
	}
	
	name = aa_defs[idx].name;
	
	int i, j, k;
	
	for (j=0; aa_defs[idx].aabonds[j]; j++)
	{	char elemsym[5];
		
		// If atom name starts with a number, then the element symbol begins on the next char.
		if (aa_defs[idx].aabonds[j]->aname[0] < 'A')
			strcpy(elemsym, &aa_defs[idx].aabonds[j]->aname[1]);
		else
			strcpy(elemsym, aa_defs[idx].aabonds[j]->aname);
		
		for (k=4; k; k--)
		{	if (elemsym[k] < 'A') elemsym[k] = 0;		// Remove any numbers from the end of the symbol.
			else
			{	elemsym[k] = 0;							// Remove the Greek letter if present and count it complete.
				break;
			}
		}
		
		if (elemsym[1]) elemsym[1] |= 0x20;
		
		if (!Atom::Z_from_esym(elemsym))
		{	cout << "Bad atom name in aminos.dat " << aa_defs[idx].aabonds[j]->aname << " for " << aa_defs[idx].name << endl;
			throw 0xbadac1d;
		}
		
		Atom* bt = 0;
		bool peptide_bond = false;
		if (aa_defs[idx].aabonds[j]->bname[0] != '<')
		{	bt = get_atom(aa_defs[idx].aabonds[j]->bname);
			if (!bt)
			{	cout << aa_defs[idx].aabonds[j]->bname << " not found or out of sequence in " << aa_defs[idx].name << endl;
				throw 0xbadac1d;
			}
		}
		else
		{	if (prevaa)
			{	bt = prevaa->get_atom(&aa_defs[idx].aabonds[j]->bname[1]);
				prev_aa = prevaa;
				prevaa->next_aa = this;
				peptide_bond = true;
			}
		}
		
		Atom* a;
		if (a = get_atom(aa_defs[idx].aabonds[j]->aname))
		{	// TODO: call close_loop().
			int m=255;
			Atom* ptr;
			Atom* path[256];
			
			bool allarom = (aa_defs[idx].aabonds[j]->cardinality > 1 && aa_defs[idx].aabonds[j]->cardinality < 2);
			//cout << "Existing atom is " << a->name << " and bond-to is " << bt->name << endl;
			path[m--] = 0;
			for (ptr=a; m; m--)
			{	path[m] = ptr;
				//cout << "Path for ring contains " << ptr->name << endl;
				if (ptr == bt) break;
				Bond** ptb = ptr->get_bonds();
				//cout << ptr->name << " is bonded to " << (ptb[0]->btom ? ptb[0]->btom->name : "null") << endl;
				if (!ptb[0]->btom) break;
				ptr = ptb[0]->btom;
				// ptb[0]->can_rotate = true;
			}
			// close_loop(&path[m], aa_defs[idx].aabonds[j]->cardinality);
			
			a->bond_to(bt, aa_defs[idx].aabonds[j]->cardinality);
			make_coplanar_ring(&path[m]);
		}	else
		{	a = add_atom(elemsym, aa_defs[idx].aabonds[j]->aname, bt, aa_defs[idx].aabonds[j]->cardinality);
			// cout << "New atom is " << a->name << " and bond-to is " << (bt ? bt->name : "null") << endl;
			if (!strcmp("CA", aa_defs[idx].aabonds[j]->aname)) a->swap_chirality = true;
			a->aaletter = letter;
			strcpy(a->aa3let, aa_defs[idx]._3let);
			a->residue = residue_no;
			if (   !strcmp(aa_defs[idx].aabonds[j]->aname, "N")
			 	|| !strcmp(aa_defs[idx].aabonds[j]->aname, "HN")
			 	|| !strcmp(aa_defs[idx].aabonds[j]->aname, "CA")
			 	|| !strcmp(aa_defs[idx].aabonds[j]->aname, "C")
			 	|| !strcmp(aa_defs[idx].aabonds[j]->aname, "O")
			   )
				a->is_backbone = true;
			else a->is_backbone = false;
		}
		
		if (peptide_bond)
		{	// a->aromatize();
			a->mirror_geo = 0;
			a->flip_mirror = true;
		}
	}
	
	// hydrogenate();
	minimize_internal_collisions();
}

void AminoAcid::save_pdb(FILE* os, int atomno_offset)
{	int i;
	
	for (i=0; atoms[i]; i++)
	{	atoms[i]->save_pdb_line(os, i+1+atomno_offset);
	}
}

int AminoAcid::from_pdb(FILE* is)
{	/*
ATOM     55  SG  CYS     4       6.721  -8.103   4.542  1.00001.00           S
*/
	char buffer[1024] = {}, res3let[5] = {};
	int added=0, lasttell=0;
	residue_no=0;

	while (!feof(is))
	{	lasttell = ftell(is);
		fgets(buffer, 1003, is);
		char** fields = chop_spaced_fields(buffer);

		if (fields)
		{	if (!strcmp(fields[0], "ATOM")
				||
				!strcmp(fields[0], "HETATM")
			   )
			{	try
				{	// cout << "Resno " << fields[4] << " vs old " << resno << endl;
					if (!residue_no) residue_no = atoi(fields[4]);
					if (!res3let[0])
					{	strcpy(res3let, fields[3]);
						/*int i;
						// cout << res3let << "|";
						for (i=0; i<26; i++)
						{	// cout << aa_defs[i]._3let << "|";
							if (!strcmp(res3let, aa_defs[i]._3let))
							{	cout << res3let << " matches " << aa_defs[i]._3let << ". AA is " << aa_defs[i]._1let << endl;
								aadef = &aa_defs[i];
								break;
							}
						}
						if (!aadef)
						{	aadef = new AADef();
							aadef->_1let = '?';
							cout << res3let << " does not match anything." << endl;
						}*/
					}
					
					if (!atno_offset) atno_offset = atoi(fields[1]);
					
					if (strcmp(res3let, fields[3])
						||
						residue_no != atoi(fields[4])
					   )
					{	fseek(is, lasttell, SEEK_SET);
						goto _return_added;
					}
					
					char esym[7] = {0,0,0,0,0,0,0};
					if (fields[2][0] >= '0' && fields[2][0] <= '9')
						strcpy(esym, &fields[2][1]);
					else
						strcpy(esym, fields[2]);
					
					int i;
					for (i=1; i<6; i++)
					{	if (!esym[i+1]) esym[i] = 0;
						if (esym[i+1] >= '0' && esym[i+1] <= '9') esym[i]=0;
						if (i>1) esym[i] = 0;
						if (!esym[i]) break;
					}
					esym[1] &= 0x5f;
					
					Point aloc(atof(fields[5]), atof(fields[6]),atof(fields[7]));
					
					Atom* a = add_atom(esym, fields[2], &aloc, 0, 0);
					added++;
					
					if (   !strcmp(a->name, "N")
					 	|| !strcmp(a->name, "HN")
					 	|| !strcmp(a->name, "CA")
					 	|| !strcmp(a->name, "C")
					 	|| !strcmp(a->name, "O")
					   )
						a->is_backbone = true;
					else a->is_backbone = false;
					
					a->residue = atoi(fields[4]);
					strcpy(a->aa3let, fields[3]);
					AADef* aaa=0;
					
					name=0;
					for (i=0; i<26; i++)
					{	if (aa_defs[i]._1let && !strcmp(aa_defs[i]._3let, fields[3]))
						{	a->aaletter = aa_defs[i]._1let;
							aaa = &aa_defs[i];
							name = new char[10]{}; // aa_defs[i].name;
							sprintf(name, "%s%d", aa_defs[i]._3let, atoi(fields[4]));
							break; 
						}
					}
					
					if (!aaa)
					{	fseek(is, lasttell, SEEK_SET);
						delete[] fields;
						throw ATOM_NOT_OF_AMINO_ACID;
					}
					else aadef = aaa;
					
					if (aaa && aaa->aabonds)
						for (i=0; aaa->aabonds[i]; i++)
						{	AABondDef* aab = aaa->aabonds[i];
							if (!strcmp(aab->aname, a->name))
							{	Atom* btom = get_atom(aab->bname);
								if (btom)
								{	a->bond_to(btom, aab->cardinality);
									if (aab->cardinality == 1 && !aab->can_rotate)
									{	Bond* b = a->get_bond_between(btom);
										if (b) b->can_rotate = false;
										b = btom->get_bond_between(a);
										if (b) b->can_rotate = false;
									}
								}
							}
						}
					
				}	catch (int ex)
				{	if (ex == ATOM_NOT_OF_AMINO_ACID) throw ex;
				}
			}
		}
		
		delete[] fields;
	}
	
	_return_added:
	
	return added;
}

void AminoAcid::copy_loaded_to_object(char letter, int tbdctr, AABondDef** tmpbdefs)
{	int lidx = (letter & 0x5f) - 'A';
	if (lidx<0 || lidx>26) return;
	aa_defs[lidx].aabonds = new AABondDef*[tbdctr+1]{};
	int j;
	for (j=0; j<tbdctr; j++)
	{	aa_defs[lidx].aabonds[j] = tmpbdefs[j];
	}
	aa_defs[lidx].aabonds[tbdctr] = 0;
}

Bond** AminoAcid::get_rotatable_bonds()
{	// Return ONLY side chain bonds, from lower to higher Greek. E.g. CA-CB but NOT CB-CA.
	// Exclude CA-N and CA-C as these will be managed by the Protein class.
	if (!atoms) return 0;
	Bond* btemp[65536];
	
	int i,j, bonds=0;
	for (i=0; atoms[i]; i++)
	{	Bond** lb = atoms[i]->get_bonds();
		int g = atoms[i]->get_geometry();
		for (j=0; j<g; j++)
		{	if (lb[j]->can_rotate
				&&
				lb[j]->atom && lb[j]->btom
				&&
				(!lb[j]->atom->is_backbone || !strcmp(lb[j]->atom->name, "CA"))
				&&
				!lb[j]->btom->is_backbone
				&&
				greek_from_aname(lb[j]->atom->name) == (greek_from_aname(lb[j]->btom->name)-1)
				&&
				lb[j]->btom->get_Z() > 1
				)
			{	btemp[bonds++] = lb[j];
				btemp[bonds] = 0;
			}
			else lb[j]->can_rotate = false;
		}
		if (lb) delete[] lb;
	}
	
	Bond** retval = new Bond*[bonds+1]{};
	for (i=0; i<=bonds; i++) retval[i] = btemp[i];
	
	return retval;
}

void AminoAcid::load_aa_defs()
{	FILE* pf = fopen(override_aminos_dat ?: "aminos.dat", "rb");
	if (!pf)
		cout << "ERROR failed to open aminos.dat, please verify file exists and you have permissions." << endl;
	else
	{	int i, j;
		Star lastgrk;
		const Star Hellenic = { .cpsz = Greek };
		char buffer[1024];
		char* lastfields[16];
		lastfields[0] = 0;
		AABondDef** tmpbdefs=0;
		int tbdctr=0;
		char lastletter = '\0';
		bool isbb = false;
		while (!feof(pf))
		{	fgets(buffer, 1011, pf);
			if (buffer[0] != '#' && buffer[0] != '\n')
			{	char** fields = chop_spaced_fields(buffer);

				try
				{	
					for (i=0; fields[i]; i++) if (fields[i][0] == '"')
					{	fields[i] = lastfields[i];
					}
					
					int idx = (fields[0][0] & 0x5f) - 'A';
					
					if (!lastletter || fields[0][0] != lastletter)
					{	copy_loaded_to_object(lastletter, tbdctr, tmpbdefs);
						tbdctr = 0;
						
						if (tmpbdefs) delete tmpbdefs;
						tmpbdefs = new AABondDef*[256]{};
					}
					
					// i = (fields[0][0] & 0x5f) - 'A';
					// cout << aa_defs[idx]._1let;
					if (!aa_defs[idx]._1let)
					{	aa_defs[idx]._1let = 'A'+idx;
						// cout << aa_defs[idx]._1let;
						strcpy(aa_defs[idx]._3let, fields[1]);
						strcpy(aa_defs[idx].name, fields[2]);
						aa_defs[idx].reach = 1.09;
					}
					
					tmpbdefs[tbdctr] = new AABondDef();
					strcpy(tmpbdefs[tbdctr]->aname, fields[3]);
					
					isbb = !strcmp(fields[3], "N")
						|| !strcmp(fields[3], "HN")
						|| !strcmp(fields[3], "CA")
						|| !strcmp(fields[3], "HA")
						|| !strcmp(fields[3], "C")
						|| !strcmp(fields[3], "O")
						;
					
					if (!isbb)
					{	int gki = strlen(fields[3])-1;
						if (fields[3][gki] >= '0' && fields[3][gki] <= '9') gki--;
						if (gki && fields[3][gki-1] >= 'A')
						{	char gk = fields[3][gki];
							lastgrk.psz = strchr(Hellenic.psz, gk);
							if (lastgrk.n) lastgrk.n -= Hellenic.n;
							float lreach = (1.54*2*cos((tetrahedral_angle-M_PI/2)*2))/2 * lastgrk.n + 1.09;
							if (lreach > aa_defs[idx].reach)
							{	aa_defs[idx].reach = lreach;
								//cout << aa_defs[idx]._3let << " has a reach of " << lreach << " because lastgreek " << Greek[lastgrk.n] << endl;
							}
						}
					}
					
					tmpbdefs[tbdctr]->cardinality = atof(fields[5]);
					tmpbdefs[tbdctr]->can_rotate
						= tmpbdefs[tbdctr]->cardinality == 1
						  &&
						  !strchr(fields[5], '!')
						;
							
					if (fields[6][0] == '+') fields[6][0] = ' ';
					tmpbdefs[tbdctr]->acharge = atof(fields[6]);
					
					char* comma = strstr(fields[4],",");
					if (comma)
					{	char* part2 = comma+1;
						*comma = 0;
						strcpy(tmpbdefs[tbdctr]->bname, fields[4]);
						// if (i==5) cout << "Fields[4] = " << fields[4] << "; part2 = " << part2 << endl;
						
						tbdctr++;
						tmpbdefs[tbdctr] = new AABondDef();
						strcpy(tmpbdefs[tbdctr]->aname, fields[3]);
						strcpy(tmpbdefs[tbdctr]->bname, part2);
						comma = strstr(fields[5],",");
						part2 = comma+1;
						*comma = 0;
						tmpbdefs[tbdctr-1]->cardinality = atof(fields[5]);
						tmpbdefs[tbdctr-1]->can_rotate
							= tmpbdefs[tbdctr-1]->cardinality == 1
							  &&
							  !strchr(fields[5], '!')
							;
						strcpy(tmpbdefs[tbdctr-1]->bname, fields[4]);
						tmpbdefs[tbdctr]->cardinality = atof(part2);
						tmpbdefs[tbdctr]->can_rotate
							= tmpbdefs[tbdctr]->cardinality == 1
							  &&
							  !strchr(part2, '!')
							;
						tmpbdefs[tbdctr]->acharge = atof(fields[6]);
					}	else
					{	strcpy(tmpbdefs[tbdctr]->bname, fields[4]);
					}
					
					if (fields[7]) strcpy(aa_defs[idx].smiles, fields[7]);
					
					tbdctr++;
					
					lastletter = fields[0][0];
				}
				catch (int e)
				{	cout << "Error while reading aminos.dat, please verify file is in the correct format." << endl;
					throw 0xbadaadef;
				}
				
				if (!lastfields[0])
					for (i=0; i<16; i++)
						lastfields[i] = new char[256]{};
					
				for (i=0; fields[i]; i++)
				{	if (lastfields[i] != fields[i]) strcpy(lastfields[i], fields[i]);
				}
				
				delete[] fields;
			}
		}
		copy_loaded_to_object(lastletter, tbdctr, tmpbdefs);
		fclose(pf);
	}
}

bool AminoAcid::can_reach(AminoAcid* other) const
{	Atom* ca1, *ca2;
	float r;
	
	ca1 = get_atom("CA");
	ca2 = other->get_atom("CA");
	
	if (!ca1 || !ca2)
	{	cout << "Warning: Could not determine reach of " << *this << " - " << *other << " to avoid possibility of collision." << endl;
		return false;
	}
	
	//cout << aadef->_3let << residue_no << " vs " << other->aadef->_3let << other->residue_no;
	
	r = ca1->get_location().get_3d_distance(ca2->get_location());
	
	//cout << ": " << r << " vs. " << get_reach() << " + " << other->get_reach() << endl;
	
	if (r <= 1.15 * (get_reach() + other->get_reach())) return true;
	else return false;
}


std::ostream& operator<<(std::ostream& os, const AminoAcid& aa)
{	if (!&aa) return os;
	try
    {	AADef* raa = aa.get_aa_definition();
    	if (raa) os << raa->_3let;
    }	catch (int ex)
    {	;
    }
    os << aa.get_residue_no();
    return os;
}


bool AminoAcid::capable_of_inter(intera_type inter)
{	if (!atoms) return false;
	
	int i, j;
	for (i=0; atoms[i]; i++)
	{	if (atoms[i]->is_backbone) continue;
		InteratomicForce** iff = InteratomicForce::get_applicable(atoms[i], atoms[i]);
		if (!iff) continue;
		for (j=0; iff[j]; j++)
			if (iff[j]->get_type() == inter)
			{	// cout << *this << " is capable of " << inter << " binding because of atom " << atoms[i]->name << endl;
				delete iff;
				return true;
			}
	}
	
	return false;
}

void AminoAcid::rotate(LocatedVector vector, float theta)
{	if (!atoms) return;
	
	int i;
	for (i=0; i<atcount; i++)
	{	Point loc = atoms[i]->get_location();
		Point nl  = rotate3D(&loc, &vector.origin, &vector, theta);
		atoms[i]->move(&nl);
	}
	
	// If you have a metal coordination, AND YOU ARE THE FIRST COORDINATING RESIDUE OF THE METAL, move the metal with you.
	if (m_mcoord && m_mcoord->coord_res && m_mcoord->coord_res[0] == this)
	{	Point loc = m_mcoord->metal->get_location();
		Point nl  = rotate3D(&loc, &vector.origin, &vector, theta);
		m_mcoord->metal->move(&nl);
	}
}

void AminoAcid::delete_sidechain()
{	Atom* latoms[get_atom_count()+4];
	int i, j=0;
	for (i=0; atoms[i]; i++)
	{	if (atoms[i]->is_backbone)
		{	latoms[j++] = atoms[i];
			latoms[j] = NULL;
		}
		else
		{	delete atoms[i];
		}
	}
	
	for (i=0; latoms[i]; i++)
		atoms[i] = latoms[i];
	
	atoms[i] = NULL;
}

Atom* AminoAcid::previous_residue_C()
{	Atom* N = get_atom("N");
	if (!N) return NULL;
	int i, n = N->get_bonded_atoms_count();
	if (!n) return NULL;
	for (i=0; i<n; i++)
	{	Atom* btom = N->get_bond_by_idx(i)->btom;
		if (!btom) continue;
		if (btom->residue == N->residue-1) return btom;
	}
	return NULL;
}

Atom* AminoAcid::next_residue_N()
{	Atom* C = get_atom("C");
	if (!C) return NULL;
	int i, n = C->get_bonded_atoms_count();
	if (!n) return NULL;
	for (i=0; i<n; i++)
	{	Atom* btom = C->get_bond_by_idx(i)->btom;
		if (!btom) continue;
		if (btom->residue == C->residue+1) return btom;
	}
	return NULL;
}

LocRotation* AminoAcid::flatten()
{	LocRotation* retval = new LocRotation[5]{};
	if (m_mcoord) return retval;	// NO.

	Atom* prevC = previous_residue_C(); if (!prevC) return retval;
	Atom* prevO, *prevCA;
	int i, j, n = prevC->get_bonded_atoms_count();
	if (!n) return retval;
	for (i=0; i<n; i++)
	{	Atom* btom = prevC->get_bond_by_idx(i)->btom;
		if (!btom) continue;
		if (!strcmp(btom->name, "CA")) prevCA = btom;
		if (!strcmp(btom->name, "O" )) prevO = btom;
	}
	Atom* localN  = get_atom("N"); if (!localN) return retval;
	Atom* localHN = get_atom("HN"); if (!localHN) localHN = get_atom("H"); if (!localHN) return retval;
	Atom* localCA  = get_atom("CA"); if (!localCA) return retval;
	Atom* localC  = get_atom("C"); if (!localC) return retval;
	Atom* localO  = get_atom("O"); if (!localC) return retval;
	
	float ad[5] = { 0.1, 0.1, 0.1, 0.1, 0.1 };
	
	for (j=0; j<5; j++)
	{	bool proline = false;
		switch (j)
		{	case 0:
			// Correction about C=N axis.
			retval[0].v = v_from_pt_sub(localN->get_location(), prevC->get_location());
			retval[0].origin = localN->get_location();
			retval[0].a = 0;
			break;
			
			case 1:
			// Correction about N-HN axis.
			retval[1].v = v_from_pt_sub(localHN->get_location(), localN->get_location());
			retval[1].origin = localN->get_location();
			retval[1].a = 0;
			break;
			
			case 2:			
			// Correction about C-O axis.
			retval[2].v = v_from_pt_sub(prevO->get_location(), prevC->get_location());
			retval[2].origin = prevC->get_location();
			retval[2].a = 0;
			break;
			
			case 3:
			// Phi correction.
			retval[3].v = v_from_pt_sub(localCA->get_location(), localN->get_location());
			retval[3].origin = localCA->get_location();
			retval[3].a = 0;
			if (!localN->get_bond_between(localCA)->can_rotate) proline = true;
			break;
			
			case 4:
			// Psi correction.
			retval[4].v = v_from_pt_sub(localC->get_location(), localCA->get_location());
			retval[4].origin = localCA->get_location();
			retval[4].a = 0;
			if (!localCA->get_bond_between(localC)->can_rotate) proline = true;
			break;
			
			default:
			;
		}
		
		if (proline) continue;
		
		float planar = 9999, r = 9999;
		for (i=0; i<50; i++)
		{	
			if (i == 49 && j >= 3) ad[j] = M_PI;
			
			switch(j)
			{	case 3:
				rotate_backbone(N_asc, ad[j]);
				break;
				
				case 4:
				rotate_backbone(CA_asc, ad[j]);
				break;
				
				default:
				;
			}
			
			Point pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8;
			if (j < 3)
			{
				pt1 = prevC->get_location();
				pt2 = prevO->get_location();
				pt3 = rotate3D(localN->get_location(), retval[j].origin, retval[j].v, retval[j].a+ad[j]);
				pt4 = rotate3D(localCA->get_location(), retval[j].origin, retval[j].v, retval[j].a+ad[j]);
				pt5 = rotate3D(localHN->get_location(), retval[j].origin, retval[j].v, retval[j].a+ad[j]);
				pt6 = rotate3D(localC->get_location(), retval[j].origin, retval[j].v, retval[j].a+ad[j]);
				pt7 = rotate3D(localO->get_location(), retval[j].origin, retval[j].v, retval[j].a+ad[j]);
				pt8 = prevCA->get_location();
			}
			else
			{
				pt1 = prevC->get_location();
				pt2 = prevO->get_location();
				pt3 = localN->get_location();
				pt4 = localCA->get_location();
				pt5 = localHN->get_location();
				pt6 = localC->get_location();
				pt7 = localO->get_location();
			}
			
			float lplanar, lr;
			
			switch (j)
			{	case 0:
				lplanar = are_points_planar(pt1, pt2, pt3, pt5);
				lr = 4;
				break;
				
				case 1:
				lplanar = are_points_planar(pt1, pt2, pt3, pt4);
				lr = 4;
				break;
				
				case 2:
				lplanar = are_points_planar(pt1, pt2, pt3, pt8);
				lr = 4;
				break;
				
				case 3:
				lplanar = are_points_planar(pt3, pt4, pt5, pt6);
				lr = pt5.get_3d_distance(pt7);
				break;
				
				case 4:
				lplanar = are_points_planar(pt4, pt5, pt6, pt7);
				lr = pt5.get_3d_distance(pt7);
				break;
				
				default:
				;
			}
			
			if ( (i == 49 && j >= 3) ? (lr <= r) : (lplanar <= planar) )
			{	retval[j].a += ad[j];
				if (fabs(ad[j]) < 0.5) ad[j] *= 1.1;
				planar = lplanar;
				r = lr;
			}
			else
			{	switch(j)
				{	case 3:
					rotate_backbone(N_asc, -ad[j]);
					break;
					
					case 4:
					rotate_backbone(CA_asc, -ad[j]);
					break;
					
					default:
					;
				}
				ad[j] *= -0.5;
			}
		}
		
		LocatedVector lv = retval[j].get_lv();
		switch(j)
		{	case 0:
			case 1:
			case 2:
			rotate(lv, retval[j].a);
			break;
			
			default:
			;
		}
	}
	
	return retval;
}
	
LocRotation AminoAcid::rotate_backbone_abs(bb_rot_dir dir, float angle)
{	// For the N-CA bond, of a residue not at the N-terminus, there's one torsion angle that will place
	// the local C atom as far as possible from the previous residue's C atom.
	// For the CA-C bond, of a residue not at the C-terminus, there is a torsion angle that will place
	// the local N atom as far as possible from the next residue's N atom.
	LocRotation retval;
	if (m_mcoord) return retval;	// NO.

	// Get the location of the previous C, and the location of the local C.
	Atom *atom, *btom;
	
	if (dir == N_asc || dir == CA_desc)
		btom = previous_residue_C();
	else
		btom = next_residue_N();
	
	retval.v.theta = 0.1;
	if (!btom) return retval;
	Point prevC = btom->get_location();
	
	if (dir == N_asc || dir == CA_desc)
		btom = get_atom("C");
	else
		btom = get_atom("N");
	
	retval.v.theta = 0.2;
	if (!btom) return retval;
	Point locC = btom->get_location();
	if (dir == CA_desc || dir == CA_asc)
	{	Point swap = prevC;
		prevC = locC;
		locC = swap;
	}
	
	// Get the N-CA vector and CA as the center of rotation.
	switch (dir)
	{	case N_asc:
		atom = get_atom("N");
		btom = get_atom("CA");
		break;
		
		case CA_asc:
		atom = get_atom("CA");
		btom = get_atom("C");
		break;
		
		case CA_desc:
		atom = get_atom("CA");
		btom = get_atom("N");
		break;
		
		case C_desc:
		atom = get_atom("C");
		btom = get_atom("CA");
		break;
		
		default:
		; //return retval;
	}
	
	retval.v.theta = 0.3;
	if (!atom || !btom) return retval;
	Bond* b = atom->get_bond_between(btom);
	if (!b->can_rotate) return retval;			// Proline.
	Point origin = atom->get_location();
	Point pt = btom->get_location().subtract(origin);
	Vector axis(pt);
	
	// Use the vector as the axis and make an imaginary circle, finding the angle that separates the C atoms most.
	int i, step=10, besti=0;
	float bestrad=0, bestr = 0.00001;
	
	Atom *Cprev, *Hlocal, *Olocal, *Nnext;
	Cprev = previous_residue_C();
	Hlocal = get_atom("NH");
	if (!Hlocal) Hlocal = get_atom("H");
	Olocal = get_atom("O");
	Nnext = next_residue_N();

	float cpl = 1000, delta = 0.5, rad = bestrad;
	int k;
	
	// To this maximum stretch angle, add the input angle and do the backbone rotation.
	LocatedVector retlv = rotate_backbone(dir, bestrad+angle);
	retval.v.r = retlv.r;
	retval.v.theta = retlv.theta;
	retval.v.phi = retlv.phi;
	retval.a = bestrad+angle;
	retval.origin = retlv.origin;
	return retval;
}

LocatedVector AminoAcid::rotate_backbone(bb_rot_dir direction, float angle)
{	Atom *rotcen, *btom;
	LocatedVector retval;
	Point rel, bloc;
	char aname[5], bname[5];
	if (!atoms) return retval;
	if (m_mcoord) return retval;	// NO.
	
	// Determine the center of rotation and rotation vector.
	switch (direction)
	{	case N_asc:
		strcpy(aname, "N");
		strcpy(bname, "CA");
		break;
		
		case CA_asc:
		strcpy(aname, "CA");
		strcpy(bname, "C");
		break;
		
		case CA_desc:
		strcpy(aname, "CA");
		strcpy(bname, "N");
		break;
		
		case C_desc:
		strcpy(aname, "C");
		strcpy(bname, "CA");
		break;
		
		default:
		return retval;
	}
	rotcen = get_atom(aname);
	if (!rotcen) return retval;
	retval.origin = rotcen->get_location();
	btom = get_atom(bname);
	if (!btom) return retval;
	
	Bond* b = rotcen->get_bond_between(btom);
	if (!b->can_rotate)
	{	// cout << "FORBIDDEN cannot rotate " << *this << ":" << rotcen->name << "-" << btom->name << endl;
		retval.r = 0;		// Fail condition.
		return retval;
	}
	
	rel = btom->get_location().subtract(retval.origin);
	Vector v(rel);
	retval.phi = v.phi; retval.theta = v.theta; retval.r = v.r;	
	
	// Rotate the local atoms that would be affected.
	if (direction == C_desc)
	{	btom = get_atom("N");		if (btom) btom->move(rotate3D(btom->get_location(), retval.origin, retval, angle));
	}
	if (direction == CA_desc || direction == C_desc)
	{	btom = get_atom("HN");		if (btom) btom->move(rotate3D(btom->get_location(), retval.origin, retval, angle));
		btom = get_atom("H");		if (btom) btom->move(rotate3D(btom->get_location(), retval.origin, retval, angle));
	}
	if (direction == N_asc)
	{	btom = get_atom("C");		if (btom) btom->move(rotate3D(btom->get_location(), retval.origin, retval, angle));
	}
	if (direction == N_asc || direction == CA_asc || direction == C_desc)
	{	btom = get_atom("O");		if (btom) btom->move(rotate3D(btom->get_location(), retval.origin, retval, angle));
	}
	
	// Include side chain atoms.
	if (direction == N_asc || direction == C_desc)
	{	int i;
		for (i=0; atoms[i]; i++)
		{	if (!atoms[i]->is_backbone)
			{	atoms[i]->move(rotate3D(atoms[i]->get_location(), retval.origin, retval, angle));
			}
		}
	}
	
	return retval;
}

float AminoAcid::get_intermol_binding(AminoAcid** neighbs, bool backbone_atoms_only)
{	if (!neighbs) return 0;
	if (!neighbs[0]) return 0;
	float retval = 0;
	int i;
	if (!backbone_atoms_only)
		for (i=0; neighbs[i]; i++)
		{	Star s;
			s.paa = neighbs[i];
			retval += Molecule::get_intermol_binding(s.pmol);
		}
	else
	{	int j, k;
		for (j=0; j<atcount; j++)
			atoms[j]->last_bind_energy = 0;
		
		for (i=0; neighbs[i]; i++)
		{	for (j=0; j<atcount; j++)
			{	if (!atoms[j]->is_backbone) continue;
				Point aloc = atoms[j]->get_location();
				for (k=0; k<neighbs[i]->atcount; k++)
				{	if (!neighbs[i]->atoms[k]->is_backbone) continue;
					float r = neighbs[i]->atoms[k]->get_location().get_3d_distance(&aloc);
					float abind = InteratomicForce::total_binding(atoms[j], neighbs[i]->atoms[k]);
					if (abind && !isnan(abind) && !isinf(abind))
					{	retval += abind;
						atoms[j]->last_bind_energy += abind;
					}
				}
			}
		}
	}
	return retval;
}

Point MetalCoord::coord_atom_avg_loc()
{	Point pt(0,0,0);
	if (!coord_atoms) return pt;
	if (!coord_atoms[0]) return pt;
	int i;
	for (i=0; coord_atoms[i]; i++)
	{	pt = pt.add(coord_atoms[i]->get_location());
	}
	pt.x /= i; pt.y /= i; pt.z /= i;
	return pt;
}















