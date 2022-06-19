
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
{
    if (!aa_defs[0]._1let) AminoAcid::load_aa_defs();
    immobile = false; // true;
    movability = MOV_FLEXONLY;
    from_pdb(instream);
    minclash = get_internal_clashes();
    mol_typ = MOLTYP_AMINOACID;
}

#define _ALGORITHMIC_GREEK 1
AminoAcid::AminoAcid(const char letter, AminoAcid* prevaa)
{
    if (!aa_defs[0]._1let) AminoAcid::load_aa_defs();
    immobile = false; // true;
    movability = MOV_FLEXONLY;
    mol_typ = MOLTYP_AMINOACID;

    if (!prevaa) residue_no = 1;
    else residue_no = prevaa->residue_no + 1;

    int idx = (letter & 0x5f) - 'A';
    aadef = &aa_defs[idx];
    if (!aa_defs[idx].aabonds)
    {
        cout << "Cannot load " << letter << " please make sure amino acid exists in aminos.dat and all atoms are defined." << endl;
        throw 0xbadac1d;
    }

    name = aa_defs[idx].name;

    int i, j, k, l, n;
    
    if (aa_defs[idx].SMILES.length())
    {
    	from_smiles(aa_defs[idx].SMILES.c_str());
    	
    	#ifdef _ALGORITHMIC_GREEK
    	int ac4 = get_atom_count()+4;
    	int atom_Greek[ac4];
    	int atom_append[ac4];
    	int atom_prepend[ac4];
    	bool atom_isheavy[ac4];
    	int atom_prev[ac4];
    	int numgrk[29];			// 24 letters, one based, plus extra space.
    	
    	Atom *N = nullptr, *HN = nullptr, *C = nullptr, *O = nullptr;
    	
    	for (i=0; i<ac4; i++) atom_Greek[i] = atom_append[i] = atom_prepend[i] = 0;
    	
    	for (i=1; i<=24; i++) numgrk[i] = 0;
    	
    	for (i=0; atoms[i]; i++)
    	{
    		atoms[i]->residue = residue_no;
    		atoms[i]->aaletter = aa_defs[idx]._1let;
    		strcpy(atoms[i]->aa3let, aa_defs[idx]._3let);
    		atom_isheavy[i] = (atoms[i]->get_Z() != 1);
    		
    		// If the atom is hydrogen, its Greek will be the same as its bond[0]->btom. Otherwise, it will be one more.
    		Bond* b0 = atoms[i]->get_bond_by_idx(0);
    		if (!b0 || !b0->btom)
    		{
    			cout << "Error in definition for " << aa_defs[idx].name << "." << endl;
    			throw 0xbadac1d;
    		}
    		j = atom_idx_from_ptr(b0->btom);
    		if (j < 0)
    		{
    			cout << "Error in definition for " << aa_defs[idx].name << "." << endl;
    			throw 0xbadac1d;
    		}
    		
    		atom_prev[i] = j;
		}
    	
    	// Find the alpha carbon. Name it CA.
    	for (i=0; atoms[i]; i++)
    	{
    		if (atoms[i]->get_family() == TETREL)
    		{
    			N = atoms[i]->is_bonded_to(PNICTOGEN, 1);
    			if (N)
				{
					Bond** ab = atoms[i]->get_bonds();
					for (j=0; ab[j]; j++)
					{
						if (ab[j]->btom
							&&
							ab[j]->btom->get_family() == TETREL
							)
						{
							O = ab[j]->btom->is_bonded_to(CHALCOGEN, 2);
							if (O)
							{
								atoms[i]->name = new char[5];
								strcpy(atoms[i]->name, "CA");
								atoms[i]->is_backbone = true;
								k = atom_idx_from_ptr(atoms[i]);
								if (k >= 0)
								{
									atom_Greek[k] = 1;
									numgrk[atom_Greek[k]]++;
								}
								
								C = ab[j]->btom;
								C->name = new char[5];
								strcpy(C->name, "C");
								C->is_backbone = true;
								k = atom_idx_from_ptr(C);
								if (k >= 0) atom_Greek[k] = -1;
								
								O->name = new char[5];
								strcpy(O->name, "O");
								O->is_backbone = true;
								k = atom_idx_from_ptr(O);
								if (k >= 0) atom_Greek[k] = -1;
								
								N->name = new char[5];
								strcpy(N->name, "N");
								N->is_backbone = true;
								k = atom_idx_from_ptr(N);
								if (k >= 0) atom_Greek[k] = -1;
								
								if (N->num_bonded_to("H") > 1)		// Proline conditional.
								{
									ab = N->get_bonds();
									l = 0;
									int n;
									for (n=0; ab[n]; n++)
									{
										if (ab[n]->btom && ab[n]->btom->get_Z() == 1)
										{
											if (!l)
											{
												HN = ab[n]->btom;
												HN->name = new char[5];
												strcpy(HN->name, "HN");
												HN->is_backbone = true;
												k = atom_idx_from_ptr(HN);
												if (k >= 0) atom_Greek[k] = -1;
											}
											
											l++;
										}
									}
								}
								else HN = nullptr;
								
								goto _found_CA;
							}
						}
					}
				}
    		}
    	}
    	
    	cout << "Please check the definition for " << aa_defs[idx].name << " that its SMILES string describes an alpha amino acid." << endl;
    	throw 0xbadac1d;
    	
    	_found_CA:
    	// Set the Greek values for heavy atoms.
    	int maxgrk = 1;
    	for (i=0; atoms[i]; i++)
    	{
    		if (atom_Greek[i] < 0) continue;
    		if (atom_isheavy[i])
    		{
    			j = atom_prev[i];
    			if (atom_Greek[j] > 0 && !atom_Greek[i]) atom_Greek[i] = atom_Greek[j]+1;
    			
    			Bond** ab = atoms[i]->get_bonds();
    			int ag = atoms[i]->get_geometry();
    			for (j=0; j<ag; j++)
    			{
    				if (ab[j]->btom)
    				{
    					k = atom_idx_from_ptr(ab[j]->btom);
    					if (k >= 0 && !atom_Greek[k])
    					{
							atom_Greek[k] = atom_Greek[i]+1;
    					}
    				}
    			}
			}
		}
				
		// Recheck each heavy atom has the lowest possible Greek.
		for (n=0; n<10; n++)
		{
			for (i=0; atoms[i]; i++)
			{
				if (atom_isheavy[i] && atom_Greek[i] > 1)
				{
					Bond** ab = atoms[i]->get_bonds();
					int ag = atoms[i]->get_geometry();
					for (j=0; j<ag; j++)
					{
						if (ab[j]->btom)
						{
							// cout << atoms[i]->name << " is bonded to " << ab[j]->btom->name << name << endl;
							k = atom_idx_from_ptr(ab[j]->btom);
							// cout << "Greek is " << atom_Greek[i] << " vs. bonded Greek " << atom_Greek[k] << endl;
							if (atom_isheavy[k] && atom_Greek[k] > 0)
							{
								l = atom_Greek[k] + 1;
								if (l < atom_Greek[i] || !atom_Greek[i])
								{
									atom_Greek[i] = l;
									atom_prev[i] = k;
									// cout << "Greek is now " << atom_Greek[i] << endl;
								}
							}
						}
					}
				}
			}
    	}
    	
    	// How many heavy atoms assigned to each Greek letter?
    	for (k=1; k<=24; k++)
		{
			numgrk[k] = 0;
			for (i=0; atoms[i]; i++)
			{
				if (!atom_isheavy[i]) continue;
				if (atom_Greek[i] == k) numgrk[k]++;
				if (atom_Greek[i] > maxgrk) maxgrk = atom_Greek[i];
			}
		}
		
		// If any heavy atoms have the same Greek, suffix them 1, 2, 3....
		// But if a heavy atom is bonded to another heavy atom with a suffix, then the new atom has to have the same or greater suffix.
		for (k=1; k<=maxgrk; k++)
		{
			if (numgrk[k] < 2) continue;
			
			l=0;
			for (i=0; atoms[i]; i++)
			{
				if (!atom_isheavy[i]) continue;
				if (atom_Greek[i] == k)
				{
					if (atom_prev[i]>0
						&&
						atom_append[atom_prev[i]]
						&&
						(atom_append[atom_prev[i]]-1) >= l
						)
						l = atom_append[atom_prev[i]]-1;
					atom_append[i] = ++l;
				}
			}
		}
    	
    	// Set the Greek values for hydrogens, including the suffix if present.
    	for (i=0; atoms[i]; i++)
    	{
    		if (atom_Greek[i] < 0) continue;
    		if (!atom_isheavy[i])
    		{
    			j = atom_prev[i];
    			if (atom_Greek[j] < 0) continue;
    			atom_Greek[i] = atom_Greek[j];
    			atom_append[i] = atom_append[j];
			}
		}
    	
    	// If any hydrogens have the same Greek and zero suffix, suffix them 1, 2, 3....
    	for (i=0; atoms[i]; i++)
    	{
    		if (atom_Greek[i] < 0) continue;
    		if (atom_isheavy[i]) continue;
    		if (atom_append[i]) continue;
    		
    		for (j=i+1; atoms[j]; j++)
    		{
    			if (atom_isheavy[j]) continue;
    			if (atom_Greek[j] == atom_Greek[i] && !atom_append[j])
    			{
    				atom_append[i] = 1;
    				atom_append[j] = 2;
    				l = 2;
    				
    				for (k=j+1; atoms[k]; k++)
    				{
    					if (atom_isheavy[k]) continue;
    					if (atom_Greek[k] == atom_Greek[i] && !atom_append[k])
    						atom_append[k] = ++l;
    				}
    			}
    		}
    		
		}
    	
    	// If any hydrogens have the same Greek and the same nonzero suffix, prefix them 1, 2, 3....
    	for (i=0; atoms[i]; i++)
    	{
    		if (atom_Greek[i] < 0) continue;
    		if (atom_isheavy[i]) continue;
    		if (!atom_append[atom_prev[i]]) continue;
    		if (atom_prepend[i]) continue;
    		
    		for (j=i+1; atoms[j]; j++)
    		{
    			if (atom_isheavy[j]) continue;
    			if (atom_Greek[j] == atom_Greek[i] && atom_append[j] == atom_append[i])
    			{
    				atom_prepend[i] = 1;
    				atom_prepend[j] = 2;
    				l = 2;
    				
    				for (k=j+1; atoms[k]; k++)
    				{
    					if (atom_isheavy[k]) continue;
    					if (atom_Greek[k] == atom_Greek[i] && atom_append[k] == atom_append[i])
    						atom_prepend[k] = ++l;
    				}
    				
    				goto _nexti;
    			}
			}
			
			_nexti:
			;
    	}
    	
    	_iagtfksb:
    	// Now name all the Greek atoms.
    	std::string strGrk = Greek;
    	for (i=0; atoms[i]; i++)
    	{
    		/*cout << i
    			 << "\t" << atom_prev[i]
    			 << "\t" << atoms[i]->get_elem_sym()
    			 << "\t" << atoms[i]->name
    			 << "\t" << atom_prepend[i]
    			 << "\t" << atom_Greek[i]
    			 << "\t" << atom_append[i]
    			 << endl;*/
    		
    		if (atom_Greek[i] < 1) continue;
    		
    		std::string newname = "";
    		if (atom_prepend[i]) newname += std::to_string(atom_prepend[i]);
    		newname += atoms[i]->get_elem_sym();					// TODO: Convert to upper case.
    		newname += strGrk.substr(atom_Greek[i]-1, 1);
    		if (atom_append[i]) newname += std::to_string(atom_append[i]);
    		
    		// cout << atoms[i]->name << " is now " << newname << endl;
    		strcpy(atoms[i]->name, newname.c_str());
    	}
    	
    	// Sort all the backbone and sidechain atoms into a new array.
    	Atom** new_atoms = new Atom*[atcount];
    	for (i=0; i<atcount; i++) new_atoms[i] = nullptr;
    	
    	l=0;
    	if (N) new_atoms[l++] = N;
    	if (HN) new_atoms[l++] = HN;
    	
		for (i=0; atoms[i]; i++)
		{
			if (atom_Greek[i] <= 0) continue;
			if (atom_isheavy[i])
			{
				new_atoms[l++] = atoms[i];
				// cout << l-1 << " " << atoms[i]->name << endl;
				for (j=i+1; atoms[j]; j++)
				{
					if (!atom_isheavy[j] 
						&&
						atom_Greek[j] == atom_Greek[i]
						&&
						(!atom_append[i] || atom_append[j] == atom_append[i])
						)
					{
						new_atoms[l++] = atoms[j];
						// cout << l-1 << " " << atoms[j]->name << endl;
					}
				}
			}
			atoms[i]->clear_geometry_cache();
			atoms[i]->clear_all_moves_cache();
		}
    	
    	if (C) new_atoms[l++] = C;
    	if (O) new_atoms[l++] = O;
    	
    	_finish_anames:
    	
    	// Delete any atoms that did not get Greeked.
    	for (i=atcount-1; i>=0; i--) if (!atom_Greek[i]) delete_atom(atoms[i]);
    	
    	// Check correct bonding.
    	/*for (i=0; atoms[i]; i++)
    	{
    		Bond** bb = atoms[i]->get_bonds();
    		int ag = atoms[i]->get_geometry();
    		for (j=0; j<ag; j++)
    		{
    			if (bb[j]->btom) cout << atoms[i]->residue << ":" << atoms[i]->name
    								  << cardinality_printable(bb[j]->cardinality) 
    								  << bb[j]->btom->residue << ":" << bb[j]->btom->name << endl;
    		}
    	}*/
    	
    	if (prevaa)
    	{
    		Atom* prevC = prevaa->get_atom("C");
    		if (prevC) N->bond_to(prevC, 1);
		}
    	
    	// delete[] atoms;
    	atoms = new_atoms;
    	atcount = l;
    	
		#endif
    	
    	get_rotatable_bonds();
    }
    
	// flatten();
	minimize_internal_clashes();
	
	Atom* CA = get_atom("CA");
	Atom* CB = get_atom("CB");
	
	if (CA && CB)
	{
		Bond* b = CA->get_bond_between(CB);
		if (b) b->rotate(triangular);
	}
	
	if (prevaa)
	{
		// Make sure we have all the right atoms.
		Atom* currN = get_atom("N");
		Atom* prevC = prevaa->get_atom("C");
		if (!currN || !prevC) return;
		currN->aromatize();				// Not actually aromatic, but the pi bond geometry is the same.
		prevC->aromatize();
		Atom* prevCA = prevaa->get_atom("CA");
		Atom* prevO  = prevaa->get_atom("O");
		if (!prevCA || !prevO) return;
		Atom* currHN = get_atom("HN");
		if (!currHN) currHN = get_atom("H");
		Atom* currCA = get_atom("CA");
		
		// Proline fix.
		if (!currHN)
		{
			Bond** bb = currN->get_bonds();
			for (i=0; bb[i]; i++)
			{
				if (bb[i]->btom && bb[i]->btom != currCA)
				{
					currHN = bb[i]->btom;
					break;
				}
			}
		}
		
		if (!currHN || !currCA) return;
		
		movability = MOV_ALL;
		immobile = false;
		
		Point ptdbg(10,0,0);
		SCoord sc;
		aamove(ptdbg);
		
		// Get the curr.N - prev.C relative location and scale it to 1.32A.
		sc = prevC->get_location().subtract(currN->get_location());
		
		// Move the entire current molecule so the N is now at that distance.
		sc.r -= 1.32;
		aamove(sc);
		
		for (i=0; i<10; i++)
		{
			// Get the normal from the prev.CA - prev.C - prev.O plane.
			SCoord axis = compute_normal(prevCA->get_location(), prevC->get_location(), prevO->get_location());
			
			// Get the angle along the normal of prev.O - prev.C - curr.N.
			float theta = find_angle_along_vector(prevO->get_location(), currN->get_location(), prevC->get_location(), axis);
			float plus  = triangular - theta;
			float minus = triangular*2 - theta;		// 240deg is the same as -120deg.
			
			// Rotate the current molecule about the prev.C atom, using the normal as an axis,
			// to get a 120 degree angle not clashing curr.N with prev.CA.
			Point maybeP = rotate3D(currN->get_location(), prevC->get_location(), axis, plus);
			Point maybeM = rotate3D(currN->get_location(), prevC->get_location(), axis, minus);
			LocatedVector lv(axis);
			lv.origin = prevC->get_location();
			if (maybeP.get_3d_distance(prevCA->get_location()) > maybeM.get_3d_distance(prevCA->get_location()))
				Molecule::rotate(lv, plus);
			else
				Molecule::rotate(lv, minus);
			
			// Get the normal from the curr.HN - curr.N - curr.CA plane.
			axis = compute_normal(prevCA->get_location(), prevC->get_location(), prevO->get_location());
			
			// Get the angle along the normal of prev.C - curr.N - curr.HN.
			theta = find_angle_along_vector(prevC->get_location(), currHN->get_location(), currN->get_location(), axis);
			plus  = triangular - theta;
			minus = triangular*2 - theta;
			
			// Rotate the current molecule about the curr.N atom, using the normal as an axis,
			// to get a 120 degree angle not clashing prev.C with curr.CA.
			maybeP = rotate3D(currCA->get_location(), currN->get_location(), axis, plus);
			maybeM = rotate3D(currCA->get_location(), currN->get_location(), axis, minus);
			lv.copy(axis);
			lv.origin = currN->get_location();
			if (maybeP.get_3d_distance(prevC->get_location()) > maybeM.get_3d_distance(prevC->get_location()))
				Molecule::rotate(lv, plus);
			else
				Molecule::rotate(lv, minus);
			// return;		// Debug step.
			
			// Get the angle of curr.HN - prev.O along the prev.C - curr.N axis.
			axis = currN->get_location().subtract(prevC->get_location());
			theta = find_angle_along_vector(currHN->get_location(), prevO->get_location(), currN->get_location(), axis);
			
			// Rotate the current molecule around curr.N to get 180 degrees.
			lv.copy(axis);
			lv.origin = currN->get_location();
			Molecule::rotate(lv, M_PI+theta);
			
			// Rotate around prev.C axis prev.C-prev.O to get 180 degrees about that axis.
			axis = prevO->get_location().subtract(prevC->get_location());
			theta = find_angle_along_vector(prevCA->get_location(), currN->get_location(), prevC->get_location(), axis);
			lv.copy(axis);
			lv.origin = prevC->get_location();
			Molecule::rotate(lv, M_PI-theta);
			
			// Rotate around curr.N axis curr.N-curr.HN to get 180 degrees about that axis.
			axis = currHN->get_location().subtract(currN->get_location());
			theta = find_angle_along_vector(currCA->get_location(), prevC->get_location(), currN->get_location(), axis);
			lv.copy(axis);
			lv.origin = currN->get_location();
			Molecule::rotate(lv, M_PI+theta);
		}
		
		movability = MOV_FLEXONLY;
		immobile = true;
	}
}

LocatedVector AminoAcid::predict_next_NH()
{
	Atom* C = get_atom("C");
	if (!C) return LocatedVector();
	C->aromatize();
	Atom* CA = get_atom("CA");
	Atom* O  = get_atom("O");
	if (!CA || !O) return LocatedVector();
	
	Point nextNloc(1.32,0,0), nextHNloc(0,1,0);
	
	nextNloc = C->get_location().add(nextNloc);
	nextHNloc = nextHNloc.add(nextNloc);
	
	LocatedVector lv;
	int i;
	for (i=0; i<10; i++)
	{
		// Step 1: O-C-N angle.
		SCoord axis = compute_normal(CA->get_location(), C->get_location(), O->get_location());
		float theta = find_angle_along_vector(O->get_location(), nextNloc, C->get_location(), axis);
		float plus  = triangular - theta;
		float minus = triangular*2 - theta;
		
		// Rotate to get 120 degrees not clashing next.N with curr.CA.
		Point maybeP = rotate3D(nextNloc, C->get_location(), axis, plus);
		Point maybeM = rotate3D(nextNloc, C->get_location(), axis, minus);
		lv.copy(axis);
		lv.origin = C->get_location();
		if (maybeP.get_3d_distance(CA->get_location()) > maybeM.get_3d_distance(CA->get_location()))
		{
			nextNloc = maybeP;
			nextHNloc = rotate3D(nextHNloc, lv.origin, axis, plus);
		}
		else
		{
			nextNloc = maybeM;
			nextHNloc = rotate3D(nextHNloc, lv.origin, axis, minus);
		}
		
		// Step 2: C-N-H angle.
		axis = compute_normal(CA->get_location(), C->get_location(), O->get_location());
		theta = find_angle_along_vector(C->get_location(), nextHNloc, nextNloc, axis);
		plus  = triangular - theta;
		minus = triangular*2 - theta;
		
		// Rotate to get 120 degrees each way.
		maybeP = rotate3D(nextHNloc, nextNloc, axis, plus);
		maybeM = rotate3D(nextHNloc, nextNloc, axis, minus);
		
		// Whichever point is farthest from O is the direction of HN (trans configuration).
		lv.copy(axis);
		lv.origin = C->get_location();
		if (maybeP.get_3d_distance(O->get_location()) > maybeM.get_3d_distance(O->get_location()))
		{
			nextHNloc = maybeP;
		}
		else
		{
			nextHNloc = maybeM;
		}
		
		// Step 3: Enforce trans configuration 180 degrees.
		axis = nextNloc.subtract(C->get_location());
		theta = find_angle_along_vector(nextHNloc, O->get_location(), nextNloc, axis);
			
		// Rotate the around next.N to get 180 degrees.
		lv.copy(axis);
		lv.origin = nextNloc;
		nextHNloc = rotate3D(nextHNloc, lv.origin, axis, M_PI+theta);
			
		// Rotate around C axis C-O to get 180 degrees about that axis.
		axis = O->get_location().subtract(C->get_location());
		theta = find_angle_along_vector(CA->get_location(), nextNloc, C->get_location(), axis);
		lv.copy(axis);
		lv.origin = C->get_location();
		nextNloc = rotate3D(nextNloc, lv.origin, axis, M_PI-theta);
		nextHNloc = rotate3D(nextHNloc, lv.origin, axis, M_PI-theta);
			
	}
	
	lv.copy(nextHNloc.subtract(nextNloc));
	lv.origin = nextNloc;
	return lv;
}

Molecule** AminoAcid::aas_to_mols(AminoAcid** aas)
{
	if (!aas) return NULL;
	int i, j;
	for (i=0; aas[i]; i++);		// Get count.
	Molecule** mols = new Molecule*[i+4];
	for (j=0; j<i; j++)
	{
		mols[j] = aas[j];
	}
	mols[i] = NULL;

	return mols;
}

void AminoAcid::save_pdb(FILE* os, int atomno_offset)
{
    int i;

    for (i=0; atoms[i]; i++)
    {
        atoms[i]->save_pdb_line(os, i+1+atomno_offset);
    }
}

int AminoAcid::from_pdb(FILE* is)
{
    /*
    ATOM     55  SG  CYS     4       6.721  -8.103   4.542  1.00001.00           S
    */
    char buffer[1024], origbuf[1024], res3let[5];
    int added=0, lasttell=0;
    residue_no=0;
    res3let[0] = 0;

    while (!feof(is))
    {
        lasttell = ftell(is);
        fgets(buffer, 1003, is);
        int thistell = ftell(is);
        strcpy(origbuf, buffer);
        char** fields = chop_spaced_fields(buffer);

        if (fields)
        {
            if (!strcmp(fields[0], "ATOM")
                /*||
                !strcmp(fields[0], "HETATM")*/
               )
            {
                try
                {
                    // cout << "Resno " << fields[4] << " vs old " << resno << endl;
                    if (!residue_no) residue_no = atoi(fields[4]);
                    if (!res3let[0])
                    {
                        strcpy(res3let, fields[3]);
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
                    {
                    	/*cout << res3let << "/" << fields[3] << " | " << residue_no << "/" << fields[4] << " | "
                    		 << origbuf << endl << flush;*/
                        fseek(is, lasttell, SEEK_SET);
                        goto _return_added;
                    }

                    char esym[7] = {0,0,0,0,0,0,0};
                    if (fields[2][0] >= '0' && fields[2][0] <= '9')
                        strcpy(esym, &fields[2][1]);
                    else
                        strcpy(esym, fields[2]);

                    int i;
                    for (i=1; i<6; i++)
                    {
                        if (!esym[i+1]) esym[i] = 0;
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
                    {
                        if (aa_defs[i]._1let && !strcmp(aa_defs[i]._3let, fields[3]))
                        {
                            a->aaletter = aa_defs[i]._1let;
                            aaa = &aa_defs[i];
                            name = new char[10]; // aa_defs[i].name;
                            sprintf(name, "%s%d", aa_defs[i]._3let, atoi(fields[4]));
                            break;
                        }
                    }

                    if (!aaa)
                    {
                        fseek(is, lasttell, SEEK_SET);
                        delete[] fields;
                        throw ATOM_NOT_OF_AMINO_ACID;
                    }
                    else aadef = aaa;

                    if (aaa && aaa->aabonds)
                    {
                        for (i=0; aaa->aabonds[i]; i++)
                        {
                            AABondDef* aab = aaa->aabonds[i];
                            if (!strcmp(aab->aname, a->name))
                            {
                                if (aaa->aabonds[i] && aaa->aabonds[i]->acharge)
                                {
                                    a->increment_charge(aaa->aabonds[i]->acharge);
                                }
                                Atom* btom = get_atom(aab->bname);
                                if (btom)
                                {
                                    a->bond_to(btom, aab->cardinality);
                                    if (aab->cardinality == 1 && !aab->can_rotate)
                                    {
                                        Bond* b = a->get_bond_between(btom);
                                        if (b) b->can_rotate = false;
                                        b = btom->get_bond_between(a);
                                        if (b) b->can_rotate = false;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                    	Atom* bta = nullptr;		// bond to atom.
                    	float cardinality = 1;		// TODO:
                    	
                    	// If HN, bond to N.
                    	// If CA, bond to N.
                    	if (!strcmp(a->name, "HN") || !strcmp(a->name, "CA") || !strcmp(a->name, "H"))
                    	{
                    		bta = get_atom("N");
                    	}
                    	
                    	// If has Greek > A: if hydrogen, bond to same Greek (mind any suffix!), else bind to earlier Greek.
                    	
                    	// If C, bond to CA.
                    	if (!strcmp(a->name, "C"))
                    	{
                    		bta = get_atom("CA");
                    	}
                    	
                    	// If O, bond to C.
                    	if (!strcmp(a->name, "O"))
                    	{
                    		bta = get_atom("C");
                    		cardinality = 2;
                    	}
                    	
                    	if (bta) a->bond_to(bta, cardinality);
                    	
                    }

                }
                catch (int ex)
                {
                    if (ex == ATOM_NOT_OF_AMINO_ACID) throw ex;
                }
            }
        }
        buffer[0] = 0;

        delete[] fields;
    }

_return_added:

    return added;
}

void AminoAcid::copy_loaded_to_object(char letter, int tbdctr, AABondDef** tmpbdefs, bool proline_like)
{
    int lidx = (letter & 0x5f) - 'A';
    if (lidx<0 || lidx>26) return;
    aa_defs[lidx].aabonds = new AABondDef*[tbdctr+1];
    int j;
    for (j=0; j<tbdctr; j++)
    {
        aa_defs[lidx].aabonds[j] = tmpbdefs[j];
        aa_defs[lidx].proline_like = proline_like;
    }
    aa_defs[lidx].aabonds[tbdctr] = 0;
}



void AminoAcid::load_aa_defs()
{
    FILE* pf = fopen(override_aminos_dat ?: "aminos.dat", "rb");
    if (!pf)
    {
        cout << "ERROR failed to open aminos.dat, please verify file exists and you have permissions." << endl;
    }
    else
    {
        int i, j;
        Star lastgrk;
        const Star Hellenic = { .cpsz = Greek };
        char buffer[1024];
        char* lastfields[16];
        lastfields[0] = 0;
        AABondDef** tmpbdefs=0;
        int tbdctr=0;
        char lastletter = '\0';
        bool isbb = false;
        bool proline_like = false;
        while (!feof(pf))
        {
            fgets(buffer, 1011, pf);
            if (buffer[0] != '#' && buffer[0] != '\n')
            {
                char** fields = chop_spaced_fields(buffer);
                if (!fields) continue;

                try
                {
                    for (i=0; fields[i]; i++) if (fields[i][0] == '"')
                        {
                            fields[i] = lastfields[i];
                        }

                    int idx = (fields[0][0] & 0x5f) - 'A';

                    if (!lastletter || fields[0][0] != lastletter)
                    {
                        copy_loaded_to_object(lastletter, tbdctr, tmpbdefs, proline_like);
                        tbdctr = 0;
                        proline_like = false;

                        if (tmpbdefs) delete[] tmpbdefs;
                        tmpbdefs = new AABondDef*[256];
                    }

                    // i = (fields[0][0] & 0x5f) - 'A';
                    // cout << aa_defs[idx]._1let;
                    if (!aa_defs[idx]._1let)
                    {
                        aa_defs[idx]._1let = 'A'+idx;
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
                    {
                        int gki = strlen(fields[3])-1;
                        if (fields[3][gki] >= '0' && fields[3][gki] <= '9') gki--;
                        if (gki && fields[3][gki-1] >= 'A')
                        {
                            char gk = fields[3][gki];
                            lastgrk.psz = strchr(Hellenic.psz, gk);
                            if (lastgrk.n) lastgrk.n -= Hellenic.n;
                            float lreach = (1.54*2*cos((tetrahedral-M_PI/2)*2))/2 * lastgrk.n + 1.09;
                            if (lreach > aa_defs[idx].reach)
                            {
                                aa_defs[idx].reach = lreach;
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
                    if (isbb && strchr(fields[5], '!')) proline_like = true;

                    if (fields[6][0] == '+') fields[6][0] = ' ';
                    tmpbdefs[tbdctr]->acharge = atof(fields[6]);

                    char* comma = strstr(fields[4],",");
                    if (comma)
                    {
                        char* part2 = comma+1;
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
                        if (isbb && strchr(fields[5], '!')) proline_like = true;
                        strcpy(tmpbdefs[tbdctr-1]->bname, fields[4]);
                        tmpbdefs[tbdctr]->cardinality = atof(part2);
                        tmpbdefs[tbdctr]->can_rotate
                            = tmpbdefs[tbdctr]->cardinality == 1
                              &&
                              !strchr(part2, '!')
                              ;
                        if (isbb && strchr(part2, '!')) proline_like = true;
                        tmpbdefs[tbdctr]->acharge = atof(fields[6]);
                    }
                    else
                    {
                        strcpy(tmpbdefs[tbdctr]->bname, fields[4]);
                    }

                    if (fields[7])
                    {
                    	if (fields[7][0] == '=')
                    	{
                    		tmpbdefs[tbdctr]->SMILES_idx = atoi(&fields[7][1]);
                    	}
                    	else if (fields[7][0] == '@')
                    	{
                    		tmpbdefs[tbdctr]->SMILES_idx = 1000 + atoi(&fields[7][1]);
                    	}
                    	else
                    	{
                    		aa_defs[idx].SMILES = fields[7];
                    		tmpbdefs[tbdctr]->SMILES_idx = 1;
                		}
                	}

                    tbdctr++;

                    lastletter = fields[0][0];
                }
                catch (int e)
                {
                    cout << "Error while reading aminos.dat, please verify file is in the correct format." << endl;
                    throw 0xbadaadef;
                }

                if (!lastfields[0])
                    for (i=0; i<16; i++)
                        lastfields[i] = new char[256];

                for (i=0; fields[i]; i++)
                {
                    if (lastfields[i] != fields[i]) strcpy(lastfields[i], fields[i]);
                }

                delete[] fields;
            }
            buffer[0] = 0;
        }
        copy_loaded_to_object(lastletter, tbdctr, tmpbdefs, proline_like);
        fclose(pf);
    }
}

bool AminoAcid::can_reach(AminoAcid* other) const
{
    Atom* ca1, *ca2;
    float r;

    ca1 = get_atom("CA");
    ca2 = other->get_atom("CA");

    if (!ca1 || !ca2)
    {
        cout << "Warning: Could not determine reach of " << *this << " - " << *other << " to avoid possibility of clash." << endl;
        return false;
    }

    //cout << aadef->_3let << residue_no << " vs " << other->aadef->_3let << other->residue_no;

    r = ca1->get_location().get_3d_distance(ca2->get_location());

    //cout << ": " << r << " vs. " << get_reach() << " + " << other->get_reach() << endl;

    if (r <= 1.15 * (get_reach() + other->get_reach())) return true;
    else return false;
}


std::ostream& operator<<(std::ostream& os, const AminoAcid& aa)
{
    if (!&aa) return os;
    try
    {
        AADef* raa = aa.get_aa_definition();
        if (raa) os << raa->_3let;
    }
    catch (int ex)
    {
        ;
    }
    os << aa.get_residue_no();
    return os;
}


bool AminoAcid::capable_of_inter(intera_type inter)
{
    if (!atoms) return false;

    int i, j;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        InteratomicForce** iff = InteratomicForce::get_applicable(atoms[i], atoms[i]);
        if (!iff) continue;
        for (j=0; iff[j]; j++)
            if (iff[j]->get_type() == inter)
            {
                // cout << *this << " is capable of " << inter << " binding because of atom " << atoms[i]->name << endl;
                delete[] iff;
                return true;
            }
    }

    return false;
}

void AminoAcid::aamove(SCoord move_amt)
{
	if (movability < MOV_ALL) return;
    if (!atoms) return;
    int i;

    for (i=0; atoms[i]; i++)
    {
        Point loc = atoms[i]->get_location();
        loc = loc.add(&move_amt);
        atoms[i]->move(&loc);
    }

    // If you have a metal coordination, AND YOU ARE THE FIRST COORDINATING RESIDUE OF THE METAL, move the metal with you.
    if (m_mcoord && m_mcoord->coord_res && m_mcoord->coord_res[0] == this)
    {
        m_mcoord->metal->move(move_amt);
    }
}

void AminoAcid::rotate(LocatedVector SCoord, float theta)
{
    if (!atoms) return;
    // cout << name << " AminoAcid::rotate()" << endl;

    int i;
    for (i=0; i<atcount; i++)
    {
        Point loc = atoms[i]->get_location();
        Point nl  = rotate3D(&loc, &SCoord.origin, &SCoord, theta);
        atoms[i]->move(&nl);
    }

    // If you have a metal coordination, AND YOU ARE THE FIRST COORDINATING RESIDUE OF THE METAL, move the metal with you.
    if (m_mcoord && m_mcoord->coord_res && m_mcoord->coord_res[0] == this)
    {
        Point loc = m_mcoord->metal->get_location();
        Point nl  = rotate3D(&loc, &SCoord.origin, &SCoord, theta);
        m_mcoord->metal->move(&nl);
    }
}

void AminoAcid::delete_sidechain()
{
    Atom* latoms[get_atom_count()+4];
    int i, j=0;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone)
        {
            latoms[j++] = atoms[i];
            latoms[j] = NULL;
        }
        else
        {
            delete atoms[i];
        }
    }

    for (i=0; latoms[i]; i++)
        atoms[i] = latoms[i];

    atoms[i] = NULL;
}

Atom* AminoAcid::previous_residue_C()
{
    Atom* N = get_atom("N");
    if (!N) return NULL;
    int i, n = N->get_bonded_atoms_count();
    if (!n) return NULL;
    for (i=0; i<n; i++)
    {
        Atom* btom = N->get_bond_by_idx(i)->btom;
        if (!btom) continue;
        if (btom->residue == N->residue-1) return btom;
    }
    return NULL;
}

Atom* AminoAcid::next_residue_N()
{
    Atom* C = get_atom("C");
    if (!C) return NULL;
    int i, n = C->get_bonded_atoms_count();
    if (!n) return NULL;
    for (i=0; i<n; i++)
    {
        Atom* btom = C->get_bond_by_idx(i)->btom;
        if (!btom) continue;
        if (btom->residue == C->residue+1) return btom;
    }
    return NULL;
}

LocRotation* AminoAcid::flatten()
{
	Bond* b;
    LocRotation* retval = new LocRotation[5];
    if (m_mcoord) return retval;	// NO.
    
    retval[0] = rotate_backbone_abs(N_asc, M_PI);
    retval[1] = rotate_backbone_abs(CA_asc, M_PI);
    return retval;

    Atom* prevC = previous_residue_C();
    if (!prevC) return retval;
    Atom* prevO, *prevCA;
    int i, j, n = prevC->get_bonded_atoms_count();
    if (!n) return retval;
    for (i=0; i<n; i++)
    {
        Atom* btom = prevC->get_bond_by_idx(i)->btom;
        if (!btom) continue;
        if (!strcmp(btom->name, "CA")) prevCA = btom;
        if (!strcmp(btom->name, "O" )) prevO = btom;
    }
    Atom* localN  = get_atom("N");
    if (!localN) return retval;
    Atom* localHN = get_atom("HN");
    if (!localHN) localHN = get_atom("H");
    if (!localHN) return retval;
    Atom* localCA  = get_atom("CA");
    if (!localCA) return retval;
    Atom* localC  = get_atom("C");
    if (!localC) return retval;
    Atom* localO  = get_atom("O");
    if (!localC) return retval;

    float ad[5] = { 0.1, 0.1, 0.1, 0.1, 0.1 };

    for (j=0; j<5; j++)
    {
        bool proline = false;
        switch (j)
        {
        	case 0:
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
            b = localCA->get_bond_between(localC);
            if (b && !b->can_rotate) proline = true;
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
            {
            case 3:
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
            {
            case 0:
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
            {
                retval[j].a += ad[j];
                if (fabs(ad[j]) < 0.5) ad[j] *= 1.1;
                planar = lplanar;
                r = lr;
            }
            else
            {
                switch(j)
                {
                case 3:
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
        {
        case 0:
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
{
    // For the N-CA bond, of a residue not at the N-terminus, there's one torsion angle that will place
    // the local C atom as far as possible from the previous residue's C atom.
    // For the CA-C bond, of a residue not at the C-terminus, there is a torsion angle that will place
    // the local N atom as far as possible from the next residue's N atom.
    LocRotation retval;
    if (m_mcoord) return retval;	// NO.
    //return retval;

    // Get the location of the previous C, and the location of the local C.
    Atom *atom, *btom;
    
    angle = M_PI+angle;
    
    switch (dir)
    {
    	case N_asc:
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
        return retval;
    }
    if (!atom || !btom) return retval;
    Bond* b = atom->get_bond_between(btom);
    if (!b) return retval;
    if (!m_mcoord) b->can_rotate = true;
    
	switch (dir)
    {
    	case N_asc:
        case CA_desc:
        atom = get_atom("HN");
        if (!atom) get_atom("H");
        btom = get_atom("C");
        break;

    	case CA_asc:
        case C_desc:
        atom = get_atom("N");
        btom = get_atom("O");
        break;

    	default:
        return retval;
    }
    
    // Use the SCoord as the axis and make an imaginary circle, finding the angle that brings HN and O closest.
    int i, step=3;
    float bestrad=0, bestr;
    
    if (!atom || !btom) return retval;
    bestr = atom->get_location().get_3d_distance(btom->get_location());
    for (i=0; i<360; i+=step)
    {
    	b->rotate(fiftyseventh*step, true);
    	float r = atom->get_location().get_3d_distance(btom->get_location());
    	if (r < bestr)
    	{
    		bestrad = fiftyseventh*i;
    		bestr = r;
    	}
    }

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
{
    Atom *rotcen, *btom;
    LocatedVector retval;
    Point rel, bloc;
    char aname[5], bname[5];
    if (!atoms) return retval;
    if (m_mcoord) return retval;	// NO.
    // return retval;

    // Determine the center of rotation and rotation SCoord.
    switch (direction)
    {
		case N_asc:
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
    if (!b || !b->can_rotate)
    {
        // cout << "FORBIDDEN cannot rotate " << *this << ":" << rotcen->name << "-" << btom->name << endl;
        retval.r = 0;		// Fail condition.
        return retval;
    }

    rel = btom->get_location().subtract(retval.origin);
    SCoord v(rel);
    retval.phi = v.phi;
    retval.theta = v.theta;
    retval.r = v.r;

    // Rotate the local atoms that would be affected.
    if (direction == C_desc)
    {
        btom = get_atom("N");
        if (btom) btom->move(rotate3D(btom->get_location(), retval.origin, retval, angle));
    }
    if (direction == CA_desc || direction == C_desc)
    {
        btom = get_atom("HN");
        if (btom) btom->move(rotate3D(btom->get_location(), retval.origin, retval, angle));
        btom = get_atom("H");
        if (btom) btom->move(rotate3D(btom->get_location(), retval.origin, retval, angle));
    }
    if (direction == N_asc)
    {
        btom = get_atom("C");
        if (btom) btom->move(rotate3D(btom->get_location(), retval.origin, retval, angle));
    }
    if (direction == N_asc || direction == CA_asc || direction == C_desc)
    {
        btom = get_atom("O");
        if (btom) btom->move(rotate3D(btom->get_location(), retval.origin, retval, angle));
    }

    // Include side chain atoms.
    if (direction == N_asc || direction == C_desc)
    {
        int i;
        for (i=0; atoms[i]; i++)
        {
            if (!atoms[i]->is_backbone)
            {
                atoms[i]->move(rotate3D(atoms[i]->get_location(), retval.origin, retval, angle));
            }
        }
    }

    return retval;
}

float AminoAcid::get_intermol_binding(AminoAcid** neighbs, bool backbone_atoms_only)
{
    if (!neighbs) return 0;
    if (!neighbs[0]) return 0;
    float retval = 0;
    int i;
    if (!backbone_atoms_only)
        for (i=0; neighbs[i]; i++)
        {
            Star s;
            s.paa = neighbs[i];
            retval += Molecule::get_intermol_binding(s.pmol);
        }
    else
    {
        int j, k;
        for (j=0; j<atcount; j++)
            atoms[j]->last_bind_energy = 0;

        for (i=0; neighbs[i]; i++)
        {
            for (j=0; j<atcount; j++)
            {
                if (!atoms[j]->is_backbone) continue;
                Point aloc = atoms[j]->get_location();
                for (k=0; k<neighbs[i]->atcount; k++)
                {
                    if (!neighbs[i]->atoms[k]->is_backbone) continue;
                    float r = neighbs[i]->atoms[k]->get_location().get_3d_distance(&aloc);
                    float abind = InteratomicForce::total_binding(atoms[j], neighbs[i]->atoms[k]);
                    if (abind && !isnan(abind) && !isinf(abind))
                    {
                        retval += abind;
                        atoms[j]->last_bind_energy += abind;
                    }
                }
            }
        }
    }
    return retval;
}

Point MetalCoord::coord_atom_avg_loc()
{
    Point pt(0,0,0);
    if (!coord_atoms) return pt;
    if (!coord_atoms[0]) return pt;
    int i;
    for (i=0; coord_atoms[i]; i++)
    {
        pt = pt.add(coord_atoms[i]->get_location());
    }
    pt.x /= i;
    pt.y /= i;
    pt.z /= i;
    return pt;
}















