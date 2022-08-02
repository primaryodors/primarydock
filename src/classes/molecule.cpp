
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <ctime>
#include <stdint.h>
#include "molecule.h"
#include "aminoacid.h"

#define _DBGCLSLOOP 0

using namespace std;

float potential_distance = 0;

Molecule::Molecule(char const* lname)
{
    name = new char[strlen(lname)+1];
    strcpy(name, lname);

    atoms = nullptr;
    smiles = nullptr;
    rings = nullptr;
    atcount = 0;
    reset_conformer_momenta();
    rotatable_bonds = nullptr;
}

Molecule::Molecule()
{
    atoms = nullptr;
    smiles = nullptr;
    rings = nullptr;
    atcount = 0;
    reset_conformer_momenta();
    rotatable_bonds = 0;
    paren = nullptr; // not sure what a good default is here, but it was not initialized (warning from clang)
}

int length(Atom** array)
{
    int numAtoms;
    for (numAtoms=0; array[numAtoms]; numAtoms++);	// Get count.
    return numAtoms;
}

bool hasAtoms(Atom** array)
{
    if (array == nullptr)
        return false;
    return array[0] != nullptr;
}

bool noAtoms(Atom** array)
{
    return !hasAtoms(array);
}

Molecule::Molecule(char const* lname, Atom** collection)
{
    if (!collection)
    {
        cout << "Temporary molecule creation attempted from nullptr atom pointer array." << endl;
        throw 0xbadca22;
    }

    int numAtoms = length(collection);
    atoms = new Atom*[numAtoms + 1];
    for (int j=0; j < numAtoms; j++)
        atoms[j] = collection[j];
    atoms[numAtoms] = 0;
    atcount = numAtoms;

    name = new char[strlen(lname)+1];
    strcpy(name, lname);

    smiles = nullptr;
    rings = nullptr;
    reset_conformer_momenta();
    rotatable_bonds = 0;
}

Molecule::~Molecule()
{
	;
}


Pose::Pose()
{
	reset();
}

Pose::Pose(Molecule* m)
{
	reset();
	copy_state(m);
}

void Pose::reset()
{
	sz = 0;
	saved_atom_locs = nullptr;
	saved_from = nullptr;
}

void Pose::copy_state(Molecule* m)
{
	if (!saved_atom_locs || saved_from != m)
	{
		if (saved_atom_locs > reinterpret_cast<void*>(0xff)) delete[] saved_atom_locs;
		saved_from = m;
		if (!m || !m->atoms) return;
		
		sz = m->get_atom_count();
		saved_atom_locs = new Point[sz+16];
	}
	
	int i;
	for (i=0; m->atoms[i] && i<sz; i++)
	{
		saved_atom_locs[i] = m->atoms[i]->get_location();
	}
}

void Pose::restore_state(Molecule* m)
{	
	if (!m || !m->atoms || !sz || m != saved_from) return;
	
	int i;
	for (i=0; i<sz && m->atoms[i]; i++)
	{
		m->atoms[i]->move(saved_atom_locs[i]);
	}
}



void Molecule::delete_atom(Atom* a)
{
	if (!a) return;
	// return;
	int i, j;
	
	if (hasAtoms(atoms))
	{
		for (i=0; atoms[i]; i++)
		{
			if (atoms[i] == a)
			{
				a->unbond_all();
				for (j=i+1; atoms[j]; j++) atoms[j-1] = atoms[j];
				atoms[j-1] = nullptr;
				rotatable_bonds = nullptr;
				atcount--;
				return;
			}
		}
	}
	
	cout << "Attempt to delete atom " << a->name << " from a molecule it is not part of." << endl << flush;
	throw 0xbada70b;
}

void Molecule::reset_conformer_momenta()
{
    srand (time(nullptr));

    lmx = _def_lin_momentum * sgn(0.5-(rand()&1));
    lmy = _def_lin_momentum * sgn(0.5-(rand()&1));
    lmz = _def_lin_momentum * sgn(0.5-(rand()&1));
    amx = _def_ang_momentum * sgn(0.5-(rand()&1));
    amy = _def_ang_momentum * sgn(0.5-(rand()&1));
    amz = _def_ang_momentum * sgn(0.5-(rand()&1));

    Bond** b = get_rotatable_bonds();
    int i;

    if (b)
    {
        for (i=0; b[i]; i++)
        {
            b[i]->angular_momentum = _def_bnd_momentum * sgn(0.5-(rand()&1));
        }
    }
}

void Molecule::reallocate()
{
    if (!(atcount % _def_atc))
    {
        int oac = atcount;
        int ac1 = oac + _def_atc;
        Atom** latoms = new Atom*[ac1+10];

        int i;
        if (atoms && oac)
        {
            for (i=0; i<oac; i++) latoms[i] = atoms[i];
        }
        for (i=oac; i<ac1; i++) latoms[i] = nullptr;
        // delete[] atoms;
        atoms = latoms;
    }
    rotatable_bonds = nullptr;
}

Atom* Molecule::add_atom(char const* elemsym, char const* aname, const Point* location, Atom* bond_to, const float bcard)
{
    Atom* a = new Atom(elemsym, location);
    a->name = new char[strlen(aname)+1];
    a->residue = 0;
    strcpy(a->name, aname);
    if (bond_to && bcard) a->bond_to(bond_to, bcard);

    reallocate();
    atoms[atcount++] = a;
    atoms[atcount] = nullptr;

    return a;
}

Atom* Molecule::add_atom(char const* elemsym, char const* aname, Atom* bondto, const float bcard)
{
    if (!bondto || !bcard)
    {
        Point pt;
        return add_atom(elemsym, aname, &pt, bondto, bcard);
    }
    
    // cout << "Add new " << elemsym << " bonded to " << bondto->name << endl;

    SCoord v = bondto->get_next_free_geometry(bcard);
    Atom* a = new Atom(elemsym);
    a->name = new char[strlen(aname)+1];
    a->residue = 0;
    strcpy(a->name, aname);

    try
    {
        v.r = InteratomicForce::covalent_bond_radius(a, bondto, bcard);
        Point pt(&v);
        Point loc = bondto->get_location();
        pt = pt.add(&loc);
        a->move(&pt);

        reallocate();
        atoms[atcount++] = a;
        atoms[atcount] = nullptr;
        a->bond_to(bondto, bcard);

        if ((atcount & 1) && bondto->get_bonded_atoms_count() == 2)
        {
            Bond* b = bondto->get_bond_between(a);
            if (b && b->can_rotate)
            {
                b->clear_moves_with_cache();
                b->rotate(M_PI);
            }
        }

        return a;
    }
    catch (int err)
    {
        if (err == BOND_DEF_NOT_FOUND)
            cout << "Covalent bond parameters not found for " << elemsym
                 << " to " << bondto->get_elem_sym()
                 << " cardinality " << bcard
                 << ". Please check bindings.dat." << endl;
        throw 0xbad6160;
    }
}

void Molecule::clear_all_bond_caches()
{
    Bond** b = get_rotatable_bonds();
    if (b)
    {
        int i;
        for (i=0; b[i]; i++) b[i]->clear_moves_with_cache();
    }
}

void Molecule::hydrogenate(bool steric_only)
{
    if (noAtoms(atoms)) return;
    int i, j;

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_metal()) continue;					// No hydrides.
        if (atoms[i]->dnh) continue;

        float valence = atoms[i]->get_valence();
        if (valence > 4) valence = 8 - valence;
        
        // cout << atoms[i]->name << " has valence " << valence << endl;

        float bcardsum = 0;

        Bond** aib = atoms[i]->get_bonds();
		int db = 0;
        if (aib)
        {
            for (j=0; aib[j]; j++)
            {
                if (aib[j]->btom) bcardsum += aib[j]->cardinality;
                if (aib[j]->cardinality > 1) db++;
            }
        }
        if (!db && steric_only) continue;

        //cout << " minus existing bonds " << bcardsum ;

        bcardsum -= atoms[i]->get_charge();

        //cout << " given charge makes " << bcardsum << endl;

		int h_to_add = round(valence - bcardsum);
        for (j=0; j<h_to_add; j++)
        {
            char hname[5];
            sprintf(hname, "H%d", atcount+1);
            Atom* H = add_atom("H", hname, atoms[i], 1);
            //cout << "Adding " << hname << " to " << atoms[i]->name << " whose valence is " << valence << " and has " << bcardsum << " bonds already." << endl;
            
            /*atoms[i]->clear_geometry_cache();
            SCoord v = atoms[i]->get_next_free_geometry(1);
            v.r = InteratomicForce::covalent_bond_radius(atoms[i], H, 1);
            H->move(atoms[i]->get_location().add(v));*/
            
            if (atoms[i]->get_geometry() == 3)
            {
                Bond* aib = atoms[i]->get_bond_by_idx(1);
                if (1) // aib && aib->total_rotations) // aib->btom == H)
                {
                    /*Bond* b0 = atoms[i]->get_bond_by_idx(0);
                    Bond* b1 = atoms[i]->get_bond_by_idx(1);
                    if (!b1 || !b1->btom) b1 = atoms[i]->get_bond_by_idx(2);*/
                    
                    int k=0;
                    Bond* b0 = atoms[i]->get_bond_by_idx(k++);
                    if (!b0->btom || b0->btom == H) b0 = atoms[i]->get_bond_by_idx(k++);
                    Bond* b1 = atoms[i]->get_bond_by_idx(k++);
                    if (!b1->btom || b1->btom == H) b1 = atoms[i]->get_bond_by_idx(k++);
                    
                    if (b0->btom && b1->btom)
                    {
                        Point source = atoms[i]->get_location();
                        /*Point axis = b0->btom->get_location();
                        Point avoid = b1->btom->get_location();
                        Point movable = H->get_location();
                        SCoord v(axis.subtract(source));

                        Point visibly_moved = rotate3D(movable, source, v, M_PI);
                        float r_is = movable.get_3d_distance(avoid);
                        float r_wouldbe = visibly_moved.get_3d_distance(avoid);

                        if (r_is < r_wouldbe) H->move(visibly_moved);*/
                        
                        Point movable = b1->btom->get_location();
                        /*cout << "Getting normal of " << atoms[i]->name << " - "
                        	 << b0->btom->name << " - " << b1->btom->name << endl;*/
                        SCoord axis = compute_normal(source, b0->btom->get_location(), b1->btom->get_location());
                        Point plus  = rotate3D(movable, source, axis,  triangular);
                        Point minus = rotate3D(movable, source, axis, -triangular);
                        
                        float rp = plus.get_3d_distance(b0->btom->get_location());
                        float rm = minus.get_3d_distance(b0->btom->get_location());
                        
                        Point pt = ((rp > rm) ? plus : minus).subtract(source);
                        pt.scale(InteratomicForce::covalent_bond_radius(atoms[i], H, 1));
                        H->move(pt.add(source));
                    }
                }
            }
        }

    }

    clear_all_bond_caches();
}

char** Molecule::get_atom_names() const
{
    int i;
    char** retval = new char*[atcount+1];

    for (i=0; i<atcount; i++) retval[i] = atoms[i]->name;
    retval[atcount] = 0;

    return retval;
}

Atom* Molecule::get_atom(char const* aname) const
{
    if (noAtoms(atoms)) return 0;

    int i;
    for (i=0; atoms[i]; i++)
        if (!strcmp(atoms[i]->name, aname)) return atoms[i];

    return 0;
}

int Molecule::atom_idx_from_ptr(Atom* a)
{
	if (!atoms) return -1;
	int i;
	for (i=0; atoms[i]; i++) if (atoms[i] == a) return i;
	
	return -1;		// If not found.	
}

Point Molecule::get_atom_location(char const * aname)
{	if (noAtoms(atoms))
	{
		Point pt;
		return pt;
	}
	Atom* a = get_atom(aname);
	if (!a) return get_barycenter();
	return a->get_location();
}

Atom* Molecule::get_nearest_atom(Point loc) const
{
    if (noAtoms(atoms)) return 0;

    int i, j;
    float best, r;
    for (i=0; atoms[i]; i++)
    {
        r = loc.get_3d_distance(atoms[i]->get_location());
        if (!i || r < best)
        {
            j=i;
            best=r;
        }
    }

    return atoms[j];
}

Atom* Molecule::get_nearest_atom(Point loc, intera_type capable_of) const
{
    if (noAtoms(atoms)) return 0;

    int i, j;
    float best, r;
    for (i=0; atoms[i]; i++)
    {
        if (!InteratomicForce::atom_is_capable_of(atoms[i], capable_of)) continue;
        r = loc.get_3d_distance(atoms[i]->get_location());
        if (!i || r < best)
        {
            j=i;
            best=r;
        }
    }

    return atoms[j];
}

int Molecule::from_sdf(char const* sdf_dat)
{
    if (!sdf_dat) return 0;
    char const* lines[8192];
    int i,j=0,lncount;

    immobile = false;
    
    // cout << sdf_dat << endl;

    lines[j++] = sdf_dat;
    for (i=0; sdf_dat[i]; i++)
    {
        if (sdf_dat[i] == '\n')
        {
            lines[j++] = (sdf_dat+i+1);
        }
        if (j > 8190) break;
    }
    lines[lncount = j] = 0;

    int na, nb;
    int added=0;
    char** fields;

    for (j=3; j<lncount; j++)
    {
        char line[1024];
        strncpy(line, lines[j], 1023);
        fields = chop_spaced_fields(line);

        if (!fields || !fields[0] || !fields[1]) break;
        if (!strcmp(fields[1], "END")) break;

        if (!strcmp(fields[1], "CHG"))
        {
            for (i=3; fields[i] && fields[i+1]; i+=2)
            {
                int aidx = atoi(fields[i]);
                atoms[aidx-1]->increment_charge(atof(fields[i+1]));
            }
        }

        if (j == 3)
        {
            na = atoi(fields[0]);
            nb = atoi(fields[1]);

            atoms = new Atom*[na+4];
            // cout << "Allocated " << na << " atoms." << endl;
        }
        else if (added < na)
        {
            Point* loc = new Point(atof(fields[0]), atof(fields[1]), atof(fields[2]));
            if (fields[3][0] >= 'a' && fields[3][0] <= 'z') fields[3][0] -= 0x20;
            Atom* a = new Atom(fields[3], loc);
            delete loc;
            a->name = new char[16];
            sprintf(a->name, "%s%d", fields[3], added+1);
            a->residue = 0;
            atoms[atcount++] = a;
            atoms[atcount] = nullptr;
            added++;
            
            // cout << "Added " << a->name << endl;
        }
        else
        {
            int a1i = atoi(fields[0]);
            int a2i = atoi(fields[1]);

            if (!a1i || !a2i) break;
            atoms[a1i-1]->bond_to(atoms[a2i-1], atoi(fields[2]));
            // cout << "Bonded " << atoms[a1i-1]->name << " to " << atoms[a2i-1]->name << endl;
        }

        if (fields) delete[] fields;
    }
    atoms[atcount] = 0;
    if (fields) delete[] fields;

    identify_rings();
    identify_acidbase();
    return added;
}

int Molecule::from_pdb(FILE* is)
{
    /*
    ATOM     55  SG  CYS     4       6.721  -8.103   4.542  1.00001.00           S
    */
    char buffer[1024];
    int added=0;

    while (!feof(is))
    {
        fgets(buffer, 1003, is);
        char** fields = chop_spaced_fields(buffer);
        
        if (fields)
        {
            if (!strcmp(fields[0], "ATOM")
                    ||
                    !strcmp(fields[0], "HETATM")
               )
            {
                try
                {
                    char esym[7];
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

                    // a->residue = atoi(fields[4]);

                    for (i=0; atoms[i]; i++)
                    {
                        if (atoms[i] == a) continue;
                        float r = aloc.get_3d_distance(atoms[i]->get_location());
                        if (r < 5)
                        {
                            if (r < 1.05* InteratomicForce::covalent_bond_radius(a, atoms[i], 1))
                                a->bond_to(atoms[i], 1);
                        }
                    }

                }
                catch (int ex)
                {
                    ;
                }
            }
        }
        buffer[0] = 0;

        delete[] fields;
    }

    return added;
}

int Molecule::get_bond_count(bool unidirectional) const
{
    int i, j, bc=0;

    if (noAtoms(atoms)) return 0;
    for (i=0; i<atcount && atoms[i]; i++)
    {
        Bond** b = atoms[i]->get_bonds();
        if (!b) continue;

        for (j=0; b[j]; j++)
        {
            if (b[j]->btom > atoms[i] || !unidirectional) bc++;
        }

        delete[] b;
    }

    return bc;
}

int Molecule::aidx(Atom* a)
{
    if (!a) return -1;
    int i;
    if (noAtoms(atoms)) return -1;

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i] == a) return i;
    }
    return -1;
}

bool Molecule::save_sdf(FILE* os)
{
    return save_sdf(os, 0);
}

bool Molecule::save_sdf(FILE* os, Molecule** lig)
{
    if (!os) return false;
    fprintf(os, "%s\n", name);

    time_t now = time(0);
    tm *gmtm = gmtime(&now);
    
    // If we used obabel or another third party app, give due credit.
    if (sdfgen_aboutline.length()) fprintf(os, "%s\n", sdfgen_aboutline.c_str());
    else
    {
		fprintf(os, "  primarydock-%02d%02d%02d%02d%02d%02d3D\n", gmtm->tm_year % 100, gmtm->tm_mon+1, gmtm->tm_mday,
		        gmtm->tm_hour, gmtm->tm_min, gmtm->tm_sec
		       );

    	fprintf(os, "https://github.com/primaryodors/primarydock\n");
	}

    int ac, bc, chargeds=0;
    ac = get_atom_count();
    bc = get_bond_count(true);

    int i, j, k, l;
    Atom* latoms[65536];
    Bond* lbonds[65536];

    if (hasAtoms(atoms))
        for (j=0; atoms[j]; j++)
            latoms[j] = atoms[j];

    Bond** b = get_all_bonds(true);
    if (b)
        for (l=0; b[l]; l++)
            lbonds[l] = b[l];
    if (b) delete[] b;

    if (lig)
    {
        for (i=0; lig[i] && lig[i]->atoms; i++)
        {
            ac += lig[i]->get_atom_count();
            bc += lig[i]->get_bond_count(true);

            for (k=0; lig[i]->atoms[k]; k++)
                latoms[j++] = lig[i]->atoms[k];

            b = lig[i]->get_all_bonds(true);

            for (k=0; b[k]; k++)
                lbonds[l++] = b[k];

            if (b) delete[] b;
        }
    }

    fprintf(os, " %d %d  0     0  0  0  0  0  0999 V2000\n", ac, bc );

    for (i=0; i<ac; i++)
    {
        Point p = latoms[i]->get_location();
        char const* esym = latoms[i]->get_elem_sym();
        if (!esym) continue;

        if (latoms[i]->get_charge()) chargeds++;

        // This avoids some weird "-0.0000" bullcrap that messes up the alignment and corrupts the SDF.
        if (!p.x) p.x=0;
        if (!p.y) p.y=0;
        if (!p.z) p.z=0;

        /*fprintf(os, "   %s%5.4f   %s%5.4f   %s%5.4f %s%s  0  0  0  0  0  0  0  0  0  0  0  0\n",
        			(p.x<0)?"":" ",p.x, (p.y<0)?"":" ",p.y, (p.z<0)?"":" ",p.z, esym, esym[1]?"":" "
        	   );*/
        char buffer[256];
        string str;
        sprintf(buffer, "%5.4f", p.x);
        str = buffer;
        str = str_pad(str, 10, " ", STR_PAD_LEFT);
        fprintf(os, "%s", str.c_str());
        sprintf(buffer, "%5.4f", p.y);
        str = buffer;
        str = str_pad(str, 10, " ", STR_PAD_LEFT);
        fprintf(os, "%s", str.c_str());
        sprintf(buffer, "%5.4f", p.z);
        str = buffer;
        str = str_pad(str, 10, " ", STR_PAD_LEFT);
        fprintf(os, "%s", str.c_str());
        fprintf(os, " %s%s  0  0  0  0  0  0  0  0  0  0  0  0\n", esym, esym[1]?"":" ");
    }


    for (i=0; i<bc; i++)
    {
        int laidx=0, lbidx=0;

        for (j=0; j<ac; j++)
        {
            if (latoms[j] == lbonds[i]->atom) laidx = j+1;
            if (latoms[j] == lbonds[i]->btom) lbidx = j+1;
        }

        // fprintf(os, " %d %d  %d  0  0  0  0\n", laidx, lbidx, (int)lbonds[i]->cardinality);
        char buffer[256];
        string str;
        sprintf(buffer, "%d", laidx);
        str = buffer;
        str = str_pad(str, 3, " ", STR_PAD_LEFT);
        fprintf(os, "%s", str.c_str());
        sprintf(buffer, "%d", lbidx);
        str = buffer;
        str = str_pad(str, 3, " ", STR_PAD_LEFT);
        fprintf(os, "%s", str.c_str());
        sprintf(buffer, " %3.2f", lbonds[i]->cardinality);
        str = buffer;
        str = str_pad(str, 3, " ", STR_PAD_LEFT);
        fprintf(os, "%s", str.c_str());
        fprintf(os, "  0  0  0  0\n");
    }

    if (chargeds)
    {
        int thisline = min(chargeds, 8);
        fprintf(os, "M  CHG  %d  ", thisline);		// TODO: Multiline if chargeds>8.
        k = 0;
        for (i=0; i<ac; i++)
        {
            float chg = latoms[i]->get_charge();
            if (!chg) continue;
            int ichg = (chg<0) ? floor(chg) : ceil(chg);
            fprintf(os, "%d  %d  ", i+1, ichg);
            k++;
            if (k > 7)
            {
                fprintf(os, "\n");
                chargeds -= k;
                if (chargeds <= 0) break;
                thisline = min(chargeds, 8);
                fprintf(os, "M  CHG  %d  ", thisline);
            }
        }
        fprintf(os, "\n");
    }

    fprintf(os, "M  END\n\n");

    return true;
}

void Molecule::save_pdb(FILE* os, int atomno_offset)
{
    int i;

    for (i=0; atoms[i]; i++)
    {
        atoms[i]->save_pdb_line(os, i+1+atomno_offset);
    }

    fprintf(os, "\nTER\nEND\n");
}

int Molecule::add_ring(Atom** atoms)
{
	int i, ringcount;
	
	if (rings)
	{
		for (i=0; rings[i]; i++);	// Get count.
		ringcount = i;
	}
	else
	{
		ringcount=0;
	}
	
	Ring** ringstmp = new Ring*[ringcount+4];
	
	if (rings)
	{
		for (i=0; i<ringcount; i++) ringstmp[i] = rings[i];
		delete[] rings;
	}
	
	ringstmp[ringcount++] = new Ring(atoms);
	ringstmp[ringcount] = nullptr;
	rings = ringstmp;
	
	return ringcount-1;
}

#define DBG_FND_RNGS 0
int Molecule::identify_rings()
{
    Atom** ringstmp[256];
    int ringcount;
    Atom *a;
    int chainlen[256];
    bool is_ring[256];
    int found_rings=0, chains=0, cnvchain, active, i, j, k, l, m, n, p;
    Atom *cnva, *cnvb;
    Atom *ra, *rb;

    ringcount = 0;

    // Start at any atom, mark it "used".
    a = atoms[0];
    a->used = true;
    
    #if DBG_FND_RNGS
    cout << "Identifying rings for " << (name ? name : "(no name molecule)") << "... " << endl;
    #endif
    
    for (i=0; i<256; i++) ringstmp[i] = nullptr;

    // For each "unused" valence>1 atom attached to the start atom, mark it "used" and start a chain.
    // Include the start atom as the first element of the chain.

    Bond** b = a->get_bonds();
    if (!b) return found_rings;
    for (i=0; b[i]; i++)
    {
        if (b[i]->btom)
            if (b[i]->btom->get_valence() > 1)
            {
                ringstmp[chains] = new Atom*[256];
                for (j=0; j<256; j++) ringstmp[chains][j] = nullptr;
                ringstmp[chains][0] = a;
                ringstmp[chains][1] = b[i]->btom;
                ringstmp[chains][2] = 0;
                is_ring[chains] = 0;
                chainlen[chains] = 2;
                b[i]->btom->used = true;
                chains++;
            }
    }

    if (!chains) return found_rings;

    // LOOP
    do
    {
        active = 0;

		#if DBG_FND_RNGS
        for (i=0; i<chains; i++)
            if (ringstmp[i][0])
            {
                cout << "Begin chain " << i << ": "; Atom::dump_array(ringstmp[i]); cout << endl;
            }
        #endif

        for (i=0; i<chains; i++)
        {
            if (!ringstmp[i] || !chainlen[i]) continue;

            cnva = cnvb = 0;
            cnvchain = 0;

            if (!(a = ringstmp[i][chainlen[i]-1])) continue;
            #if DBG_FND_RNGS
            cout << "Chain " << i << " ends in atom " << a->name << endl;
            #endif

            if (!cnva)
            {
                for (l=0; l<chains; l++)
                {
                    #if DBG_FND_RNGS
                    cout << "Preparing to check chain " << i << " against chain " << l << " for bond convergence...\n";
                    #endif
                    
                    if (l != i && chainlen[l] && ringstmp[l][chainlen[l]-1])
                    {
                        #if DBG_FND_RNGS
                        cout << "Checking chain " << i << " against chain " << l << " for bond convergence...\n";
                        #endif
                        
                        // If two chains converge on a single atom, this is your one converging ato
                        if (ringstmp[l][chainlen[l]-1] == a)
                        {
                            cnva = a;
                            cnvb = 0;
                            
                            #if DBG_FND_RNGS
                            cout << "****** Chains " << i << " and " << l << " converge on atom " << cnva->name << " ******" << endl;
                            #endif
                            
                            cnvchain = l;
                            break;
                        }
                        // If two chains converge on two atoms that are bonded to each other, then there are two converging atoms.
                        else if (ringstmp[l][chainlen[l]-1]->is_bonded_to(a))
                        {
                            cnva = a;
                            cnvb = ringstmp[l][chainlen[l]-1];
                            
                            #if DBG_FND_RNGS
                            cout << "***** Chains " << i << " and " << l << " converge on atoms " << cnva->name << " and " << cnvb->name << " *****" << endl;
                            #endif
                            
                            cnvchain = l;
                            break;
                        }
                    }
                }
_done_l:
                ;
            }



            if (cnva && cnvchain != i)
            {
                // Count backwards from one chain until find the diverging atom, then jump to the other chain and count forward
                // until reach the/a converging atom. This is a RING; add it to the list.
                n = chainlen[i] + chainlen[cnvchain];
                Atom* ring_atoms[n];
                for (m=0; m<n; m++) ring_atoms[m] = nullptr;
                n = 0;
                
                #if DBG_FND_RNGS
                cout << "Counting backwards from " << chainlen[i]-1 << endl;
                #endif
                
                for (m = chainlen[i]-1; m >= 0; m--)
                {
                    ring_atoms[n++] = ringstmp[i][m];
                    
                    #if DBG_FND_RNGS
                    cout << "m: " << ringstmp[i][m]->name << " ";
                    #endif
                    
                    l = in_array(reinterpret_cast<void*>(ringstmp[i][m]),
                                 reinterpret_cast<void**>(ringstmp[cnvchain])
                                );
                    if (l >= 0)
                    {
                        l++;
                        
                        #if DBG_FND_RNGS
                        cout << "\nCounting forwards from " << l << endl;
                        #endif
                        
                        for (; ringstmp[cnvchain][l]; l++)
                        {
                            #if DBG_FND_RNGS
                            cout << "l: " << ringstmp[cnvchain][l]->name << " ";
                            #endif
                            
                            ring_atoms[n++] = ringstmp[cnvchain][l];
                            
                            #if DBG_FND_RNGS
                            cout << "Building " << n << " membered ring: "; Atom::dump_array(ring_atoms); cout << endl;
                            #endif
                            
                            if (n >= 3 && (ringstmp[cnvchain][l] == cnva || ringstmp[cnvchain][l] == cnvb))
                            {
                                #if DBG_FND_RNGS
                                cout << "Found " << n << " membered ring: "; Atom::dump_array(ring_atoms); cout << endl;
                                #endif
                                
                                ring_atoms[n] = 0;

                                int nringid = add_ring(ring_atoms);
#if _ALLOW_FLEX_RINGS
                                if (!rings[nringid]->is_coplanar())
#else
                                if (1)
#endif
                                {
                                    for (p=0; (ra = ring_atoms[p]); p++)
                                    {
                                        rb = ring_atoms[p ? (p-1) : (n-1)];
                                        int card = (int)ra->is_bonded_to(rb);

                                        Bond* ab = ra->get_bond_between(rb);
                                        ab->can_rotate = false;
                                        ab = rb->get_bond_between(ra);
                                        ab->can_rotate = false;
                                    }
                                }

                                goto _exit_m;
                            }
                        }
                    }
                }
_exit_m:
                ;

                // TODO: The algorithm can sometimes give incorrect results in polycyclic molecules
                // if a second ring converges onto a used ato
            }
            // if (b) delete[] b;



            // Advance each chain one step further, incorporating only "unused" atoms.
            b = a->get_bonds();
            if (!b) continue;
            k=0;
            for (j=0; b[j]; j++)
            {
                if (b[j] && b[j]->btom && b[j]->btom->get_valence() > 1)
                {
                    if (!b[j]->btom->used)
                    {
                        // If there is more than one "unused" bonded atom of an eligible element, create new chains for the surplus atoms.
                        if (k)
                        {
                            ringstmp[chains] = new Atom*[256];
                            for (l=0; l<chainlen[i]-1; l++)
                                ringstmp[chains][l] = ringstmp[i][l];
                            ringstmp[chains][l] = b[j]->btom;
                            ringstmp[chains][l+1] = 0;
                            is_ring[chains] = 0;
                            chainlen[chains] = chainlen[i];
                            b[j]->btom->used = true;
                            
                            #if DBG_FND_RNGS
                            cout << "Another new chain " << chains << ": "; Atom::dump_array(ringstmp[chains]); cout << endl;
                            #endif
                            
                            chains++;
                            k++;
                            // a = 0;
                            // return 0;
                        }
                        else
                        {
                            ringstmp[i][chainlen[i]++] = a = b[j]->btom;
                            ringstmp[i][chainlen[i]] = 0;
                            b[j]->btom->used = true;
                            k++;
                            
                            #if DBG_FND_RNGS
                            cout << "Chain " << i << ": "; Atom::dump_array(ringstmp[i]); cout << endl;
                            #endif
                        }
                    }
                }
            }

            // If there are no "unused" bonded atoms, delete the chain.
            if (!k) chainlen[i] = 0;
            else active++;

        }		// for i

        //cout << "------" << endl;
    }
    while (active);

    // if (b) delete[] b;
    for (i=0; ringstmp[i]; i++) delete[] ringstmp[i];

    for (i=0; atoms[i]; i++) atoms[i]->used = false;

    // Return the number of rings found.
    return found_rings;
}

void Molecule::identify_acidbase()
{
    if (noAtoms(atoms)) return;

    // For every atom in the molecule:
    int i, j, k;
    Bond** b = 0;

    for (i=0; atoms[i]; i++)
    {
        // If it is a carbon, double-bonded to any atom, not bonded to a pnictogen,
        // Or if it is a heavier tetrel, a pnictogen, a chalcogen, or a halogen,
        // And it is single-bonded to at least one chalcogen,
        // And the single-bonded chalcogen is either bonded to a hydrogen or carries a negative charge:
        // Then all chalcogens bonded to the carbon are acidic.
        int sbOH = 0;
        int fama = atoms[i]->get_family();
        bool carbon = false;
        if (fama == TETREL || fama == PNICTOGEN || fama == CHALCOGEN || fama == HALOGEN)
        {
            carbon = (atoms[i]->get_Z() == 6);
            //cout << "Atom " << atoms[i]->name << " is of family " << fama << endl;
            b = atoms[i]->get_bonds();
            if (!b) goto _not_acidic;
            if (carbon)
            {
                if (!atoms[i]->is_pi()) goto _not_acidic;
                for (j=0; b[j]; j++)
                {
                    if (!b[j]->btom) continue;
                    if (b[j]->cardinality == 2)
                    {
                        int fam = b[j]->btom->get_family();
                        if (fam != CHALCOGEN) goto _not_acidic;
                    }
                }
            }
            for (j=0; b[j]; j++)
            {
                if (!b[j]->btom) continue;
                int fam = b[j]->btom->get_family();
                if (carbon && fam == PNICTOGEN) goto _not_acidic;
                //cout << "Fam: " << fam << endl;
                if (fam == CHALCOGEN && b[j]->cardinality < 2)
                {
                    if (b[j]->btom->get_charge() < 0)
                    {
                        sbOH++;
                        break;
                    }
                    else
                    {
                        Bond** b1 = b[j]->btom->get_bonds();
                        if (!b1) goto _not_acidic;
                        for (k=0; b1[k]; k++)
                        {
                            if (!b1[k]->btom) continue;
                            if (b1[k]->btom->get_Z() == 1)
                            {
                                sbOH++;
                                //cout << "OH" << endl;
                                break;
                            }
                        }
                        if (b1) delete[] b1;
                    }
                }
            }
            delete[] b;
        }
        if (sbOH)
        {
            b = atoms[i]->get_bonds();
            for (j=0; b[j]; j++)
            {
                if (!b[j]->btom) continue;
                int fam = b[j]->btom->get_family();
                if (fam == CHALCOGEN)
                {
                    b[j]->btom->set_acidbase(-1);
                    //cout << "Atom " << b[j]->btom->name << " is acidic." << endl;
                }
            }
            if (b) delete[] b;
        }
_not_acidic:
        //if (b) delete[] b;

        // If it is a pnictogen,
        // And it is not bonded to a chalcogen or a halogen,
        // And it is not part of an amide,
        // Or if it has a positive charge,
        // TODO: Or if any hydrogen bound to it has a positive charge,
        // Then it is basic.
        if (fama == PNICTOGEN)
        {
            Atom* c;
            c = atoms[i]->is_bonded_to("C");
            if (c)
            {
                if (c->is_bonded_to("O")) goto _not_basic;
            }
            if (atoms[i]->get_charge() <= 0)
            {
                if (atoms[i]->is_bonded_to(CHALCOGEN)) goto _not_basic;
                if (atoms[i]->is_bonded_to(HALOGEN)) goto _not_basic;
            }

            atoms[i]->set_acidbase(1);

        }
_not_basic:
        ;
    }
}

Bond** Molecule::get_rotatable_bonds()
{
    if (noAtoms(atoms)) return 0;
    if (rotatable_bonds) return rotatable_bonds;
    if (mol_typ == MOLTYP_AMINOACID)
    {	
    	// TODO: There has to be a better way.
    	Star s;
    	s.pmol = this;
    	rotatable_bonds = s.paa->get_rotatable_bonds();
    	return rotatable_bonds;
    }
    // cout << name << " Molecule::get_rotatable_bonds()" << endl << flush;

    Bond* btemp[65536];
    int mwblimit = atcount/2;						// Prevent rotations that move most of the molecule.

    int i,j, bonds=0;
    if (!immobile)
        for (i=0; atoms[i]; i++)
        {
            Bond** lb = atoms[i]->get_bonds();
            int g = atoms[i]->get_geometry();
            for (j=0; j<g && lb[j]; j++)
            {
                if (lb[j]->count_moves_with_btom() > mwblimit) continue;
                
                if (!lb[j]->atom || !lb[j]->btom) continue;
                
                bool pia = lb[j]->atom->is_pi(),
                	 pib = lb[j]->btom->is_pi();
                
                int fa = lb[j]->atom->get_family(),
                	fb = lb[j]->btom->get_family();
                
                // Generally, a single bond between pi atoms, or a bond from a pi atom to an amino group, cannot rotate.
                if (	(pia && pib)
                		||
                		(pia && (fb == PNICTOGEN))
                		||
                		(pib && (fa == PNICTOGEN))
                   )
                   lb[j]->can_rotate = false;
                
                if (lb[j]->btom
                        &&
                        lb[j]->atom < lb[j]->btom
                        &&
                        lb[j]->can_rotate
                   )
                {
                    btemp[bonds++] = lb[j];
                    btemp[bonds] = 0;
                }
            }
            if (lb) delete[] lb;
        }
    else
        for (i=0; atoms[i]; i++)
        {
            Bond** lb = atoms[i]->get_bonds();
            int g = atoms[i]->get_geometry();
            for (j=0; j<g; j++)
            {
                if (lb[j]->count_moves_with_btom() > mwblimit) continue;
                
                // Generally, a single bond between pi atoms cannot rotate.
                if (lb[j]->atom && lb[j]->atom->is_pi()
                	&&
                	lb[j]->btom && lb[j]->btom->is_pi()
                   )
                   lb[j]->can_rotate = false;

                if (lb[j]->can_rotate
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
                {
                    btemp[bonds++] = lb[j];
                    btemp[bonds] = 0;
                }
            }
            if (lb) delete[] lb;
        }

    // cout << (name ? name : "") << " has " << bonds << " rotatable bond(s)." << endl;
    rotatable_bonds = new Bond*[bonds+1];
    for (i=0; i<=bonds; i++) rotatable_bonds[i] = btemp[i];
    rotatable_bonds[bonds] = 0;

    return rotatable_bonds;
}

void Molecule::crumple(float theta)
{
	Bond** b = get_rotatable_bonds();
	if (!b) return;
	
	float int_clsh = get_internal_clashes();
	
	int i;
	for (i=0; b[i]; i++)
	{
		float ltheta = theta*randsgn();
		b[i]->rotate(ltheta);
		if (get_internal_clashes() > int_clsh*2) b[i]->rotate(-ltheta);
	}
}

// TODO: There has to be a better way.
Bond** AminoAcid::get_rotatable_bonds()
{
    // cout << name << " AminoAcid::get_rotatable_bonds()" << endl << flush;
    // Return ONLY side chain bonds, from lower to higher Greek. E.g. CA-CB but NOT CB-CA.
    // Exclude CA-N and CA-C as these will be managed by the Protein class.
    if (noAtoms(atoms)) return 0;
    if (rotatable_bonds) return rotatable_bonds;
    if (aadef && aadef->proline_like)
    {
    	// cout << "Proline-like! No rotbonds!" << endl;
    	return nullptr;
	}
    Bond* btemp[65536];
    
    int i,j, bonds=0;
    for (i=0; i<65536; i++) btemp[i] = nullptr;
    if (aadef && aadef->aabonds)
    {
    	for (i=0; aadef->aabonds[i]; i++)
    	{
    		// cout << (name ? name : "(no name)") << "." << *(aadef->aabonds[i]) << endl;
    		if (aadef->aabonds[i]->cardinality == 1
    			&&
    			aadef->aabonds[i]->can_rotate
    		   )
    		{
    			Atom* la = get_atom(aadef->aabonds[i]->aname);
    			if (la
    				&& (!la->is_backbone || !strcmp(la->name, "CA"))
    			   )
    			{
    				Bond* lb = la->get_bond_between(aadef->aabonds[i]->bname);
    				if (!lb)
    				{
    					// TODO: Add the missing bond if possible.
    					cout << "Warning: No bond between " << la->residue << ":" << la->name
							 << " and " << aadef->aabonds[i]->bname
							 << endl << flush;
						
						Bond** lbb = la->get_bonds();
						if (lbb)
						{
							cout << la->name << " is bonded to:";
							int o, ag = la->get_geometry();
							for (o=0; o<ag; o++) if (lbb[o]->btom) cout << " " << lbb[o]->btom->name;
							cout << "." << endl;
						}
						
						Atom* lba = get_atom(aadef->aabonds[i]->bname);
						if (lba)
						{
							lbb = lba->get_bonds();
							if (lbb)
							{
								cout << lba->name << " is bonded to:";
								int o, ag = lba->get_geometry();
								for (o=0; o<ag; o++) if (lbb[o]->btom) cout << " " << lbb[o]->btom->name;
								cout << "." << endl;
							}
						}
						else cout << aadef->aabonds[i]->bname << " not found." << endl;
					}
    				else
    				{
    					// cout << (name ? name : "(no name)") << ":" << *(lb) << endl;
    					// Generally, a single bond between pi atoms cannot rotate.
						if (lb->atom->is_pi() && lb->btom && lb->btom->is_pi())
						   aadef->aabonds[i]->can_rotate = false;
    					
						lb->can_rotate = aadef->aabonds[i]->can_rotate;

						if ((!la->is_backbone || !strcmp(la->name, "CA"))
							&&
							la->get_Z() > 1
				            &&
				            (	greek_from_aname(la->name) == (greek_from_aname(lb->btom->name)+1)
				            	||
				            	greek_from_aname(la->name) == (greek_from_aname(lb->btom->name)-1)
				            )
		                   )
		                {
		                	// cout << "Included." << endl;
		                	
		                	if (greek_from_aname(la->name) < greek_from_aname(lb->btom->name))
		                		btemp[bonds++] = la->get_bond_between(lb->btom);
		                	else
						    	btemp[bonds++] = lb->btom->get_bond_between(la);
						    btemp[bonds] = 0;
						    
						    // cout << (name ? name : "(no name)") << ":" << *(btemp[bonds-1]) << endl;
						}
				    }
				}
    		}
    	}
    	
    	goto _found_aadef;
    }

	// cout << name << " ";
    for (i=0; atoms[i]; i++)
    {
        Bond** lb = atoms[i]->get_bonds();
        int g = atoms[i]->get_geometry();
        for (j=0; j<g; j++)
        {
            if (lb[j]->can_rotate
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
            {
            	// cout << *lb[j] << " ";
                btemp[bonds++] = lb[j];
                btemp[bonds] = 0;
            }
            else lb[j]->can_rotate = false;
        }
        if (lb) delete[] lb;
    }
    // cout << endl;

	_found_aadef:
    Bond** retval = new Bond*[bonds+8];
    for (i=0; i<=bonds; i++) retval[i] = btemp[i];
    retval[i] = 0;
    rotatable_bonds = retval;

    return retval;
}

float Molecule::hydrophilicity()
{
	int i, count;
	float total;
	for (i=0; atoms[i]; i++)
	{
		int Z = atoms[i]->get_Z();
		if (Z==1) continue;
		
		total += atoms[i]->hydrophilicity_rule();
		count++;
	}
	return count ? (total / count) : 0;
}

Bond** Molecule::get_all_bonds(bool unidirectional)
{
    if (noAtoms(atoms)) return 0;
    Bond* btemp[65536];

    int i,j, bonds=0;
    for (i=0; atoms[i]; i++)
    {
        Bond** lb = atoms[i]->get_bonds();
        int g = atoms[i]->get_geometry();
        for (j=0; j<g; j++)
        {
            if (lb[j]->atom < lb[j]->btom
                    ||
                    !unidirectional
               )
            {
                btemp[bonds++] = lb[j];
                btemp[bonds] = 0;
            }
        }
        if (lb) delete[] lb;
    }

    Bond** retval = new Bond*[bonds+1];
    for (i=0; i<=bonds; i++) retval[i] = btemp[i];

    return retval;
}


float Molecule::get_internal_clashes()
{
    int i, j;
    float r;
    Atom* a, *b;
    float clash = 0;

    if (noAtoms(atoms)) return 0;

    for (i=0; atoms[i]; i++)
    {
        Point pta = atoms[i]->get_location();
        float avdW = atoms[i]->get_vdW_radius();
        for (j=i+1; atoms[j]; j++)
        {
            if (atoms[i]->is_bonded_to(atoms[j])) continue;
            if (atoms[i]->shares_bonded_with(atoms[j])) continue;

            Point ptb = atoms[j]->get_location();
            float bvdW = atoms[j]->get_vdW_radius();

            r = pta.get_3d_distance(&ptb);
            if (!r) r += 10e-15;
            if (r < avdW + bvdW)
            {
                float lclash = sphere_intersection(avdW, bvdW, r);
                clash += lclash;
                
                if (false && lclash > 3)
                {
                	cout << atoms[i]->name << " clashes with " << atoms[j]->name << " by " << lclash << " cu. A. resulting in " << clash << endl;
                	int g = atoms[i]->get_geometry();
                	cout << "Geometry: " << g << endl;
                	Bond** b = atoms[i]->get_bonds();
                	int k;
                	for (k=0; k<g; k++)
                		cout << atoms[i]->name << " is bonded to " << hex << b[k]->btom << dec << " "
                			 << (b[k]->btom ? b[k]->btom->name : "") << "." << endl;
            	}
            }
        }
    }

    return clash-base_internal_clashes;
}

float Molecule::get_vdW_repulsion(Molecule* ligand)
{
	if (!ligand) return 0;
	if (ligand == this) return 0;
	if (!atoms || !ligand->atoms) return 0;
	
	int i, j;
	float retval = 0;
	
	for (i=0; atoms[i]; i++)
	{
		float achg = atoms[i]->get_charge();
		bool api = atoms[i]->is_pi();
		
		for (j=0; ligand->atoms[j]; j++)
		{
			float bchg = ligand->atoms[j]->get_charge();
			bool bpi = ligand->atoms[j]->is_pi();
			
			if (!achg || !bchg)
			{
				// TODO: Hard coded values get from bindings.dat instead.
				float rlim = 4, kJmol = 0.4;
				if (api && bpi)
				{
					rlim = 3.87;
					kJmol = 2;
				}
				float halflim = rlim/2;
				float asphere = 4.0/3 * M_PI * halflim * halflim * halflim;
				
				float r = atoms[i]->distance_to(ligand->atoms[j]);
				if (r < rlim)
				{
					retval += sphere_intersection(halflim, halflim, r) * kJmol / asphere;
				}
			}
		}
	}
	
	return retval;
}

float Molecule::get_intermol_clashes(Molecule* ligand)
{
    Molecule * ligands[4];
    ligands[0] = ligand;
    ligands[1] = nullptr;
    return get_intermol_clashes(ligands);
}

float Molecule::get_intermol_clashes(Molecule** ligands)
{
    int i, j, l;
    float r;
    Atom* a, *b;
    float clash = 0;

    if (noAtoms(atoms)) return 0;
    if (!ligands) return 0;
    if (!ligands[0]) return 0;

    for (i=0; atoms[i]; i++)
    {
        Point pta = atoms[i]->get_location();
        float avdW = atoms[i]->get_vdW_radius();
        for (l=0; ligands[l]; l++)
        {
            if (!ligands[l]->atoms) continue;
            if (ligands[l] == this) continue;
            for (j=0; ligands[l]->atoms[j]; j++)
            {
                // We still check for bonded to because we treat coordinations as bonds with cardinality <1.
                if (atoms[i]->is_bonded_to(ligands[l]->atoms[j])) continue;
                if (atoms[i]->shares_bonded_with(ligands[l]->atoms[j])) continue;

                Point ptb = ligands[l]->atoms[j]->get_location();
                float bvdW = ligands[l]->atoms[j]->get_vdW_radius();

                r = pta.get_3d_distance(&ptb) + 1e-3;
                if (r < avdW + bvdW)
                {
                	/*float confidence = 2.5;		// TODO: Get this from the PDB.
					float give = 0.5;			// TODO: Compute this from the receptor secondary structure.
					
					float allowable = give + confidence / sqrt(3);
					
					r += allowable;
					if (r > (avdW + bvdW)) r = avdW + bvdW;*/
    		
                    float lclash = sphere_intersection(avdW, bvdW, r);
                    if (lclash > 0)
                    {
                        clash += lclash;
                        /*cout << name << ":" << atoms[i]->name << " clashes with " <<
                        	ligands[l]->name << ":" << ligands[l]->atoms[j]->name << " by " << lclash << " cu. A." << endl;*/
                    }
                }
            }
        }
    }

    return clash*_kJmol_cuA;
}

void Molecule::move(SCoord move_amt)
{
    if (noAtoms(atoms)) return;
    if (immobile)
    {
        cout << "Warning: Attempt to move \"immobile\" molecule " << name << endl;
        return;
    }
    int i;

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->residue) return;					// If you must move a residue of a protein, use AminoAcid::aamove() instead.
        Point loc = atoms[i]->get_location();
        loc = loc.add(&move_amt);
        atoms[i]->move(&loc);
    }
}

void Molecule::move(Point move_amt)
{
    if (noAtoms(atoms)) return;
    if (immobile)
    {
        cout << "Warning: Attempt to move \"immobile\" molecule " << name << endl;
        return;
    }
    int i;

    for (i=0; atoms[i]; i++)
    {
        // cout << atoms[i]->name << " ";
        if (atoms[i]->residue) return;					// If you must move a residue of a protein, use AminoAcid::aamove() instead.
        Point loc = atoms[i]->get_location();
        loc = loc.add(&move_amt);
        atoms[i]->move(&loc);
    }
}

Point Molecule::get_barycenter(bool bond_weighted) const
{
    if (noAtoms(atoms))
    {
        Point pt;
        return pt;
    }

    Point locs[atcount];
    int i;

    for (i=0; i<atcount; i++)
    {
    	locs[i] = atoms[i]->get_location();
		locs[i].weight = atoms[i]->get_atomic_weight();
    	#if allow_tethered_rotations
    	if (bond_weighted)
    	{
    		locs[i].weight += atoms[i]->last_bind_energy;
    	}
    	#endif
	}

    return average_of_points(locs, atcount);
}

float Molecule::get_charge()
{
	int i;
	float charge=0;
	for (i=0; atoms[i]; i++)
	{
		charge += atoms[i]->get_charge();
	}
	return charge;
}

void Molecule::recenter(Point nl)
{
    Point loc = get_barycenter();
    Point rel = nl.subtract(&loc);
    SCoord v(&rel);
    move(v);
}

void Molecule::rotate(SCoord* SCoord, float theta, bool bond_weighted)
{
    if (noAtoms(atoms)) return;
    // cout << name << " Molecule::rotate()" << endl;

    Point cen = get_barycenter(bond_weighted);

    int i;
    for (i=0; i<atcount; i++)
    {
        // if (atoms[i]->residue) return;
        Point loc = atoms[i]->get_location();
        Point nl  = rotate3D(&loc, &cen, SCoord, theta);
        atoms[i]->move(&nl);
    }
}

void Molecule::rotate(LocatedVector SCoord, float theta)
{
    if (noAtoms(atoms)) return;

    int i;
    for (i=0; i<atcount; i++)
    {
        // if (atoms[i]->residue) return;
        Point loc = atoms[i]->get_location();
        Point nl  = rotate3D(&loc, &SCoord.origin, &SCoord, theta);
        atoms[i]->move(&nl);
    }
}

bool Molecule::shielded(Atom* a, Atom* b) const
{
    int i;
    float r = a->distance_to(b);
    float r6 = r*1.26, r125 = 1.25*r;
    if (r < 2) return false;

    Point aloc = a->get_location(), bloc = b->get_location();
    for (i=0; i<atcount; i++)
    {
        Atom* ai = atoms[i];
        if (ai == a || ai == b) continue;
        float rai = ai->distance_to(a);
        if (rai > r6) continue;
        float rbi = ai->distance_to(b);
        if (rbi > r6) continue;
        if ((rai+rbi) > r125) continue;
        Point sloc = ai->get_location();
        float f3da = find_3d_angle(&aloc, &bloc, &sloc);
        if (f3da > _shield_angle)
        {
            if (last_iter && (a->residue == 114 || b->residue == 114) && ((a->residue + b->residue) == 114))
            {
                /*cout << ai->name << " shields "
                	 << a->residue << ":" << a->name << "..."
                	 << b->residue << ":" << b->name
                	 << " angle " << (f3da*fiftyseven)
                	 << endl;*/
                return true;
            }
        }
    }

    return false;
}

float Molecule::get_atom_mol_bind_potential(Atom* a)
{
    if (noAtoms(atoms)) return 0;
    int i, j;
    float retval=0;
    potential_distance = 0;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        InteratomicForce** ifs = InteratomicForce::get_applicable(a, atoms[i]);
        if (!ifs) continue;
        for (j=0; ifs[j]; j++)
        {
            if (ifs[j]->get_type() == ionic)
            {
                if (sgn(a->get_charge()) != -sgn(atoms[i]->get_charge())) continue;
                retval += 60;
            }
            else
                retval += ifs[j]->get_kJmol();
            /*cout << a->name << " can " << ifs[j]->get_type() << " strength " << ifs[j]->get_kJmol()
            	 << " with " << (atoms[i]->aa3let ? atoms[i]->aa3let : "") << atoms[i]->residue << ":"
            	 << atoms[i]->name << endl;*/

            potential_distance += ifs[j]->get_distance();
        }
        delete[] ifs;
    }

    potential_distance /= retval;

    return retval;
}

float Molecule::get_intermol_binding(Molecule* ligand)
{
    Molecule* ligands[4];
    ligands[0] = ligand;
    ligands[1] = nullptr;
    return get_intermol_binding(ligands);
}

void Molecule::clear_atom_binding_energies()
{
    int i;
    for (i=0; i<atcount; i++)
        atoms[i]->last_bind_energy = 0;
}

float Molecule::get_intermol_potential(Molecule* ligand)
{
    Molecule* ligands[4];
    ligands[0] = ligand;
    ligands[1] = nullptr;
    return get_intermol_potential(ligands);
}

float Molecule::get_intermol_potential(Molecule** ligands)
{
    if (!ligands) return 0;
    if (!ligands[0]) return 0;
    int i, j, l, n;
    float kJmol = 0;
    
    for (i=0; i<atcount; i++)
    {
        Point aloc = atoms[i]->get_location();
        for (l=0; ligands[l]; l++)
        {
            if (ligands[l] == this) continue;
            for (j=0; j<ligands[l]->atcount; j++)
            {
                float r = ligands[l]->atoms[j]->get_location().get_3d_distance(&aloc);
                float f = 1.0 / r;		// Regular invert rather than inv square so that actual bonding will take over at short range.
                InteratomicForce** iff = InteratomicForce::get_applicable(atoms[i], ligands[l]->atoms[j]);

                if (iff) for (n=0; iff[n]; n++)
                {
                	if (iff[n]->get_type() == vdW) continue;
                	
                	if (r < iff[n]->get_distance())
                		kJmol += iff[n]->get_kJmol();
                	else
                    	kJmol += iff[n]->get_kJmol()*f;
                }
                delete[] iff;
            }
        }
    }

    return kJmol;
}

float Molecule::get_intermol_binding(Molecule** ligands)
{
    if (!ligands) return 0;
    if (!ligands[0]) return 0;
    int i, j, l;
    float kJmol = 0;
    kJmol -= get_internal_clashes();
    
    // cout << (name ? name : "") << " base internal clashes: " << base_internal_clashes << "; final internal clashes " << -kJmol << endl;

    for (i=0; i<atcount; i++)
        atoms[i]->last_bind_energy = 0;

    for (i=0; i<atcount; i++)
    {
        Point aloc = atoms[i]->get_location();
        for (l=0; ligands[l]; l++)
        {
            if (ligands[l] == this) continue;
            for (j=0; j<ligands[l]->atcount; j++)
            {
                float r = ligands[l]->atoms[j]->get_location().get_3d_distance(&aloc);
                if (r < _INTERA_R_CUTOFF)
                {
                    if (	!shielded(atoms[i], ligands[l]->atoms[j])
                            &&
                            !ligands[l]->shielded(atoms[i], ligands[l]->atoms[j])
                       )
                    {
                        float abind = InteratomicForce::total_binding(atoms[i], ligands[l]->atoms[j]);
                        if (abind && !isnan(abind) && !isinf(abind))
                        {
                            kJmol += abind;
                            atoms[i]->last_bind_energy += abind;
                            // cout << atoms[i]->name << " to " << ligands[l]->atoms[j]->name << ": " << r << " A; " << abind << " kJ/mol." << endl;
                        }
                    }
                }
            }
        }
    }
    // cout << "Total: " << kJmol << endl;
    // cout << endl;

    return kJmol;
}


void Molecule::minimize_internal_clashes()
{
	// cout << (name ? name : "(no name)");
    if (noAtoms(atoms)) return;

    int i, j, iter;
    float clash = get_internal_clashes();

    if (!clash) return;		// Already zero, nothing to decrease to.

    Bond** b = get_rotatable_bonds();
    if (!b || !b[0])
    {
    	base_internal_clashes = clash;
    	return;		// No bonds to rotate.
	}

    int numrb = 0;
    for (i=0; b[i]; i++) numrb = i;

    float angle[numrb];
    for (i=0; i<numrb; i++) angle[i] = M_PI;

    for (iter=0; iter<500; iter++)
    {
        for (i=0; i<numrb; i++)
        {
            b[i]->rotate(angle[i]);
            float clash1 = get_internal_clashes();

            if (clash1 <= clash)
            {
                clash = clash1;
                angle[i] *= 0.99;
            }
            else
            {
                b[i]->rotate(-angle[i]);		// Put it back.
                angle[i] *= -0.5;
            }
        }
    }

    base_internal_clashes = get_internal_clashes();
    // cout << " base internal clashes: " << base_internal_clashes << endl;
}



#define DBG_BONDFLEX 0
#define DBG_FLEXRES 203
#define DBG_FLEXROTB 0

void Molecule::multimol_conform(Molecule** mm, int iters, void (*cb)(int))
{
    if (!mm) return;

    int i, j, k, l, n, inplen, iter;
    float rad, bestfrrad, bestfrb;

    for (i=0; mm[i]; i++)
    {
        mm[i]->reset_conformer_momenta();
    }
    inplen = i;

    float improvement;
    for (iter=0; iter<iters; iter++)
    {
        float bind = 0, bind1;
        improvement=0;
        last_iter = (iter == (iters-1));
        for (i=0; mm[i]; i++)
        {
            bool nearby[inplen+4];
            Point icen = mm[i]->get_barycenter();

            for (j=0; mm[j]; j++)
            {
                /* if (j == i)
                {
                    nearby[j] = false;
                    continue;
                } */

                Point jcen = mm[j]->get_barycenter();
                Atom* ia = mm[i]->get_nearest_atom(jcen);
                Atom* ja = mm[j]->get_nearest_atom(icen);

                if (ia->distance_to(ja) <= 4)
                {
                    nearby[j] = true;
                }
                else
                {
                    nearby[j] = false;
                    continue;
                }

                if (nearby[j])
                {
                    // get_intermol_binding includes clashes, so we don't have to check them here.
                    bind += mm[i]->get_intermol_binding(mm[j]);
                }
            }
            mm[i]->lastbind = bind;

            float reversal = -0.666; // TODO: Make this binding-energy dependent.
            float accel = 1.1;

            /**** Linear Motion ****/
            #if allow_linear_motion
            if (mm[i]->movability >= MOV_ALL && iter >= 10)
            {
            	#if debug_break_on_move
            	mm[i]->set_atoms_break_on_move(false);
            	#endif
            	
                Point pt(mm[i]->lmx, 0, 0);
                mm[i]->move(pt);
                bind1 = 0;
                for (j=0; mm[j]; j++)
                {
                    if (!nearby[j]) continue;
                    bind1 += mm[i]->get_intermol_binding(mm[j]);
                }
                if (bind1 < bind)
                {
                	// cout << bind << " vs " << bind1 << " x" << endl;
                    pt.x = -pt.x;
                    mm[i]->move(pt);
                    mm[i]->lmx *= reversal;
                }
                else
                {
                	// cout << bind << " vs " << bind1 << " +" << endl;
                    improvement += (bind1 - bind);
                    if (fabs(mm[i]->lmx) < 0.5) mm[i]->lmx *= accel;
                    bind = bind1;
                }

                pt.x = 0;
                pt.y = mm[i]->lmy;
                mm[i]->move(pt);
                bind1 = 0;
                for (j=0; mm[j]; j++)
                {
                    if (!nearby[j]) continue;
                    bind1 += mm[i]->get_intermol_binding(mm[j]);
                }
                if (bind1 < bind)
                {
                    pt.y = -pt.y;
                    mm[i]->move(pt);
                    mm[i]->lmy *= reversal;
                }
                else
                {
                    improvement += (bind1 - bind);
                    if (fabs(mm[i]->lmy) < 0.5) mm[i]->lmy *= accel;
                    bind = bind1;
                }

                pt.y = 0;
                pt.z = mm[i]->lmz;
                mm[i]->move(pt);
                bind1 = 0;
                for (j=0; mm[j]; j++)
                {
                    if (!nearby[j]) continue;
                    bind1 += mm[i]->get_intermol_binding(mm[j]);
                }
                if (bind1 < bind)
                {
                    pt.z = -pt.z;
                    mm[i]->move(pt);
                    mm[i]->lmz *= reversal;
                }
                else
                {
                    improvement += (bind1 - bind);
                    if (fabs(mm[i]->lmz) < 0.5) mm[i]->lmz *= accel;
                    bind = bind1;
                }
                
                Point lmpt(mm[i]->lmx, mm[i]->lmy, mm[i]->lmz);
                if (lmpt.magnitude() > 1.5)
                {
                	float lmm = 0.5 / lmpt.magnitude();
                	mm[i]->lmx *= lmm;
                	mm[i]->lmy *= lmm;
                	mm[i]->lmz *= lmm;
                }
                else
                {
                	float lmm = 0.97;
                	mm[i]->lmx *= lmm;
                	mm[i]->lmy *= lmm;
                	mm[i]->lmz *= lmm;
                }

                mm[i]->lastbind = bind;

            	#if debug_break_on_move
            	mm[i]->set_atoms_break_on_move(true);
            	#endif
            }
            /**** End Linear Motion ****/
            #endif

            #if allow_axial_tumble
            /**** Axial Tumble ****/
            if (mm[i]->movability >= MOV_NORECEN)
            {
            	#if debug_break_on_move
            	mm[i]->set_atoms_break_on_move(false);
            	#endif
            	
                Point pt(1,0,0);
                SCoord v(pt);

                rad = 0;
                bestfrb = 0;
                bestfrrad = nanf("No good results.");

                if (allow_mol_fullrot_iter && !(iter % _fullrot_every))
                {
                	// cout << endl;
                    while ((M_PI*2-rad) > 1e-3)
                    {
                        mm[i]->rotate(&v, _fullrot_steprad);
                        rad += _fullrot_steprad;

                        bind1 = 0;
                        for (j=0; mm[j]; j++)
                        {
                            if (!nearby[j]) continue;
                            bind1 += mm[i]->get_intermol_binding(mm[j]) + intermol_ESP * mm[i]->get_intermol_potential(mm[j]);
                        }
                        
                        // cout << "x " << rad*fiftyseven << "deg " << bind1 << endl;

                        if (bind1 > bestfrb)
                        {
                            bestfrb = bind1;
                            bestfrrad = rad;
                        }
                    }

                    if (!isnan(bestfrrad))
                        mm[i]->rotate(&v, bestfrrad);
                }
                else
                {
                	#if monte_carlo_axial
                	Pose putitback(mm[i]);
                	float ra = fabs(mm[i]->amx);
                	ra = frand(-ra, ra);
                	mm[i]->rotate(&v, ra);
                	#else
                    mm[i]->rotate(&v, mm[i]->amx);
                    #endif
                    bind1 = 0;
                    for (j=0; mm[j]; j++)
                    {
                        if (!nearby[j]) continue;
                        bind1 += mm[i]->get_intermol_binding(mm[j]);
                    }
                    if (bind1 < bind)
                    {
                    	#if monte_carlo_axial
                    	putitback.restore_state(mm[i]);
                        mm[i]->amx *= 0.98;
                    	#else
                        mm[i]->rotate(&v, -mm[i]->amx);
                        mm[i]->amx *= reversal;
                        #endif
                    }
                    else
                    {
                        improvement += (bind1 - bind);
                        bind = bind1;
                        if (fabs(mm[i]->amx) < 0.25) mm[i]->amx *= accel;
                    }
                }

                pt.x=0;
                pt.y=1;
                SCoord v1(pt);

                rad = 0;
                bestfrb = 0;
                bestfrrad = nanf("No good results.");

                if (allow_mol_fullrot_iter && !(iter % _fullrot_every))
                {
                    while ((M_PI*2-rad) > 1e-3)
                    {
                        mm[i]->rotate(&v1, _fullrot_steprad);
                        rad += _fullrot_steprad;

                        bind1 = 0;
                        for (j=0; mm[j]; j++)
                        {
                            if (!nearby[j]) continue;
                            bind1 += mm[i]->get_intermol_binding(mm[j]) + intermol_ESP * mm[i]->get_intermol_potential(mm[j]);
                        }
                        
                        // cout << "y " << rad*fiftyseven << "deg " << bind1 << endl;

                        if (bind1 > bestfrb)
                        {
                            bestfrb = bind1;
                            bestfrrad = rad;
                        }
                    }

                    if (!isnan(bestfrrad))
                        mm[i]->rotate(&v1, bestfrrad);
                }
                else
                {
                	#if monte_carlo_axial
                	Pose putitback(mm[i]);
                	float ra = fabs(mm[i]->amy);
                	ra = frand(-ra, ra);
                	mm[i]->rotate(&v, ra);
                	#else
                    mm[i]->rotate(&v1, mm[i]->amy);
                    #endif
                    
                    bind1 = 0;
                    for (j=0; mm[j]; j++)
                    {
                        if (!nearby[j]) continue;
                        bind1 += mm[i]->get_intermol_binding(mm[j]);
                    }
                    if (bind1 < bind)
                    {
                    	#if monte_carlo_axial
                    	putitback.restore_state(mm[i]);
                        mm[i]->amy *= 0.98;
                    	#else
                        mm[i]->rotate(&v1, -mm[i]->amy);
                        mm[i]->amy *= reversal;
                    	#endif
                    }
                    else
                    {
                        improvement += (bind1 - bind);
                        bind = bind1;
                        if (fabs(mm[i]->amy) < 0.25) mm[i]->amy *= accel;
                    }
                }


                pt.y=0;
                pt.z=1;
                SCoord v2(pt);

                rad = 0;
                bestfrb = 0;
                bestfrrad = nanf("No good results.");

                if (allow_mol_fullrot_iter && !(iter % _fullrot_every))
                {
                    while ((M_PI*2-rad) > 1e-3)
                    {
                        mm[i]->rotate(&v2, _fullrot_steprad);
                        rad += _fullrot_steprad;

                        bind1 = 0;
                        for (j=0; mm[j]; j++)
                        {
                            if (!nearby[j]) continue;
                            bind1 += mm[i]->get_intermol_binding(mm[j]) + intermol_ESP * mm[i]->get_intermol_potential(mm[j]);
                        }
                        
                        // cout << "z " << rad*fiftyseven << "deg " << bind1 << endl;

                        if (bind1 > bestfrb)
                        {
                            bestfrb = bind1;
                            bestfrrad = rad;
                        }
                    }

                    if (!isnan(bestfrrad))
                        mm[i]->rotate(&v2, bestfrrad);
                }
                else
                {
                	#if monte_carlo_axial
                	Pose putitback(mm[i]);
                	float ra = fabs(mm[i]->amz);
                	ra = frand(-ra, ra);
                	mm[i]->rotate(&v, ra);
                	#else
                    mm[i]->rotate(&v2, mm[i]->amz);
                    #endif
                    
                    bind1 = 0;
                    for (j=0; mm[j]; j++)
                    {
                        if (!nearby[j]) continue;
                        bind1 += mm[i]->get_intermol_binding(mm[j]);
                    }
                    if (bind1 < bind)
                    {
                    	#if monte_carlo_axial
                    	putitback.restore_state(mm[i]);
                        mm[i]->amz *= 0.98;
                    	#else
                        mm[i]->rotate(&v2, -mm[i]->amz);
                        mm[i]->amz *= reversal;
                        //cout << "x";
                        #endif
                    }
                    else
                    {
                        improvement += (bind1 - bind);
                        bind = bind1;
                        if (fabs(mm[i]->amz) < 0.25) mm[i]->amz *= accel;
                    }
                }

                mm[i]->lastbind = bind;
                
            	#if debug_break_on_move
            	mm[i]->set_atoms_break_on_move(true);
            	#endif
            }
            /**** End Axial Tumble ****/
            #endif

			#if !monte_carlo_flex
            if ((iter % _fullrot_every)) continue;
            #endif

			#if allow_bond_rots
            /**** Bond Flexion ****/
            // cout << mm[i]->name << ": " << mm[i]->movability << endl;
            if (mm[i]->movability >= MOV_FLEXONLY)
            {
            	#if debug_break_on_move
            	mm[i]->set_atoms_break_on_move(false);
            	#endif
            	
                mm[i]->get_rotatable_bonds();
                int residue = 0;
                
				#if DBG_BONDFLEX
                if (mm[i]->rotatable_bonds && mm[i]->rotatable_bonds[0] && mm[i]->rotatable_bonds[0]->atom)
                	residue = mm[i]->rotatable_bonds[0]->atom->residue;
                #endif

                // Don't know why this is renecessary.
                if (mm[i]->movability == MOV_FLEXONLY)
                {
                    bind = 0;
                    for (j=0; mm[j]; j++)
                    {
                        if (!nearby[j]) continue;
                        // cout << ".";
                        bind += mm[i]->get_intermol_binding(mm[j]);
                    }
                }
                mm[i]->lastbind = bind;

				#if DBG_BONDFLEX
                if (DBG_FLEXRES == residue)
                	cout << "Iter" << iter << " " << mm[i]->name;
                #endif
                
                if (!iter && mm[i]->rings)
                {
                	for (l=0; mm[i]->rings[l]; l++)
                	{
                		for (n=0; mm[i]->atoms[n]; n++)
                			mm[i]->atoms[n]->is_in_ring(mm[i]->rings[l]);
                	}
                }
                
                
                if (mm[i]->rotatable_bonds)
                {
                	#if monte_carlo_flex
		            Pose putitback(mm[i]);
		            #endif
                	
                	#if DBG_BONDFLEX
                    if (DBG_FLEXRES == residue)
                    	cout << " has rotbonds ";
                    #endif
                    for (k=0; mm[i]->rotatable_bonds[k]; k++)
                    {
                    	#if DBG_BONDFLEX
                        if (DBG_FLEXRES == residue)
                        {
                        	cout << k << ".) " << *mm[i]->rotatable_bonds[k] << " ";
                        	Atom** mwb = mm[i]->rotatable_bonds[k]->get_moves_with_btom();
                        	if (mwb)
                        	{
                        		cout << "bringing ";
                        		for (j=0; mwb[j]; j++)
                        		{
                        			cout << mwb[j]->name << " ";
                        		}
                        		cout << endl;
                    		}
                		}
                        #endif

                        rad = 0;
                        bestfrb = -10000;
                        bestfrrad = nanf("No good results.");

                        if (!(iter % _fullrot_every))
                        {
                            while ((M_PI*2-rad) > 1e-3)
                            {
                                mm[i]->rotatable_bonds[k]->rotate(_fullrot_steprad);
                                rad += _fullrot_steprad;

                                bind1 = 0;
                                #if DBG_BONDFLEX
                                if (DBG_FLEXROTB == k && DBG_FLEXRES == residue)
                                	cout << endl << (rad*fiftyseven) << ": ";
                        		#endif
                                for (j=0; mm[j]; j++)
                                {
                                    if (!nearby[j]) continue;
                                    float lbind1 = 
                                    
                                    #if allow_ligand_esp
                            			 (mm[i]->mol_typ == MOLTYP_AMINOACID)
                            			 ?
                                    #endif
                            			 mm[i]->get_intermol_binding(mm[j])
                                    #if allow_ligand_esp
                            			 :  
                            			 mm[i]->get_intermol_potential(mm[j]) - 10 * mm[i]->get_internal_clashes()
                                    #endif
                            			 ;
                            		
                                    bind1 += lbind1;
                                    #if DBG_BONDFLEX
                                    if (DBG_FLEXROTB == k && DBG_FLEXRES == residue)
                                    	cout << mm[j]->name << " " << lbind1 << " ";
                        			#endif
                                }

                                if (bind1 > bestfrb)
                                {
                                    bestfrb = bind1;
                                    bestfrrad = rad;
                                }
                            }
                            #if DBG_BONDFLEX
                            if (DBG_FLEXROTB == k && DBG_FLEXRES == residue)
                            	cout << endl << "(" << (bestfrrad*fiftyseven) << "deg) ";
                			#endif

                            if (!isnan(bestfrrad))
                                mm[i]->rotatable_bonds[k]->rotate(bestfrrad);
                        }
                        else
                        {
                            float ra = mm[i]->rotatable_bonds[k]->angular_momentum;
                            
                            #if monte_carlo_flex
                            ra = frand(-fabs(ra), fabs(ra));
                            #endif
                            
                            mm[i]->rotatable_bonds[k]->rotate(ra);

                            bind1 = 0;
                            for (j=0; mm[j]; j++)
                            {
                                if (!nearby[j]) continue;
                                bind1 += mm[i]->get_intermol_binding(mm[j]);
                            }
                            if (bind1 < bind)
                            {
                            	#if monte_carlo_flex
                            	putitback.restore_state(mm[i]);
                            	mm[i]->rotatable_bonds[k]->angular_momentum *= 0.98;
                            	#else
                                mm[i]->rotatable_bonds[k]->rotate(-ra);
                                mm[i]->rotatable_bonds[k]->angular_momentum *= reversal;
                                #endif
                            }
                            else
                            {
                                improvement += (bind1 - bind);
                                bind = bind1;
                                #if monte_carlo_flex
                                putitback.copy_state(mm[i]);
                                #endif
                                if (fabs(mm[i]->rotatable_bonds[k]->angular_momentum) < 0.25)
                                	mm[i]->rotatable_bonds[k]->angular_momentum *= accel;
                            }
                        }
                    }
                }
                #if DBG_BONDFLEX
                if (DBG_FLEXRES == residue)
                	cout << endl;
    			#endif
                mm[i]->lastbind = bind;
                
            	#if debug_break_on_move
            	mm[i]->set_atoms_break_on_move(true);
            	#endif

                if (!(iter % _fullrot_every)) mm[i]->reset_conformer_momenta();
            }
            /**** End Bond Flexion ****/
            #endif

        }	// for i = 0 to iters
        // cout << "Iteration " << iter << " improvement " << improvement << endl;

        if (cb) cb(iter);
    }	// for iter.
}




Atom* numbered[10];
bool ring_warned = false;

bool Molecule::from_smiles(char const * smilesstr)
{
	if (!strchr(smilesstr, '{'))		// {AtomName} is a nonstandard feature and must be handled by POdock code, not a third party app.
	{
		// Check if OpenBabel is installed.
		FILE* pf = popen(CMD_CHECK_INSTALLED_3P_SMILES_PARSER, "r");
		if (pf)
		{
			char buffer[1024];
			fgets(buffer, 1022, pf);
			if (strlen(buffer))			// TODO: Change this to employ a regex.
			{
				fclose(pf);
				std::string sdfdat = "";
				
				// Temporarily reuse buffer as the obabel command.
				sprintf(buffer, CMD_CALL_3P_SMILES_PARSER, smilesstr);
				pf = popen(buffer, "r");
				
				// Resume using buffer with fgets().
				int lno = 0;
				while (buffer[0] != '$')
				{
					fgets(buffer, 1022, pf);
					lno++;
					sdfdat += buffer;
					
					if (lno == 2) sdfgen_aboutline = buffer;
				}
				fclose(pf);
				
				int result = from_sdf(sdfdat.c_str());
				return (result > 0);
			}
			buffer[0] = 0;
			fclose(pf);
		}
	}
	
	if (strchr(smilesstr, '!')) ring_warned = true;

    smlen = strlen(smilesstr);
    paren = new SMILES_Parenthetical[smlen];
    spnum = 0;

    int i;
    for (i=0; i<10; i++) numbered[i] = 0;
    
    bool retval = from_smiles(smilesstr, nullptr);

    for (i=0; i<spnum; i++)
    {
        retval &= from_smiles(paren[i].smilesstr, paren[i].startsfrom);
    }

    delete[] paren;
    // hydrogenate(true);
    float anomaly = correct_structure();
    if (anomaly > 0.1) cout << "ERROR: Structural anomaly = " << anomaly << endl;
    else cout << "# Structural anomaly = " << anomaly << endl;
    hydrogenate(false);

    return retval;
}

bool Molecule::from_smiles(char const * smilesstr, Atom* ipreva)
{
    Atom* stack[256];
    Atom* sequence[65536];
    bool seqarom[65536];
    int sqidx[10];
    int sp = 0;
    bool bracket=false, prevarom=false;
    Atom* bracketed=0;

    immobile = false;

    int i, j=1, k=0, l, atno=get_atom_count()+1;

    Atom* preva = ipreva;
    float card = ipreva?1:0;
    int len = strlen(smilesstr);
    int lastEZ = 0;

    int numdb = 0, dbi = 0;
    for (i=0; i<len; i++) if (smilesstr[i] == '=') numdb++;

    int EZgiven[numdb+4];
    Atom* EZatom0[numdb+4];
    Atom* EZatom1[numdb+4];

    for (i=0; i<len; i++)
    {
    	if (smilesstr[i] == '!') continue;
    	
        if (smilesstr[i] == '.')
        {
            card = 0;
            continue;
        }

        if (smilesstr[i] == '-')
        {
            card = 1;
            continue;
        }

        if (smilesstr[i] == ':')
        {
            card = 1.5;
            continue;
        }

        if (smilesstr[i] == '=')
        {
            card = 2;
            EZgiven[dbi] = 0;
            EZatom0[dbi] = preva;
            continue;
        }

        if (smilesstr[i] == '#')
        {
            card = 3;
            continue;
        }

        if (smilesstr[i] == '$')
        {
            card = 4;
            continue;
        }

        if (smilesstr[i] == '(')
        {
            // stack[sp++] = preva;
            paren[spnum].startsfrom = preva;
            paren[spnum].smilesstr = new char[smlen-i+4];
            strcpy(paren[spnum].smilesstr, &smilesstr[++i]);
            // while (smilesstr[i] != ')') i++;

            int level = 1;
            while (level)
            {
                i++;
                if (smilesstr[i] == '(') level++;
                if (smilesstr[i] == ')') level--;
                if (!smilesstr[i]) throw 0xbade9c0d;
            }

            /*for (l=0; paren[spnum].smilesstr[l]; l++)
            {	if (paren[spnum].smilesstr[l] == ')')
            	{	paren[spnum].smilesstr[l+1] = nullptr;
            		break;
            	}
            }*/
            spnum++;
            continue;
        }

        if (smilesstr[i] == ')')
        {
            if (!sp) return false;
            if (ipreva) return true;
            // preva = stack[--sp];
            continue;
        }

        // Nonstandard feature.
        if (smilesstr[i] == '{')
        {
            char* aname = new char[15];
            i++;
            j=0;
            while (smilesstr[i] != '}' && j<15)
            {
                aname[j] = smilesstr[i];
                i++;
                j++;
            }
            aname[j] = 0;
            preva->name = aname;
            //i++;
            continue;
        }

        if (smilesstr[i] == '\\'
                ||
                smilesstr[i] == '/'
           )
        {
            if (lastEZ)
            {
                int EZ = (smilesstr[i] == '/') ? 1 : -1;
                preva->EZ_flip = EZgiven[dbi-1] = sgn(EZ) != sgn(lastEZ);
                // cout << preva->name << " EZ flip: " << preva->EZ_flip << endl;
                lastEZ = 0;
            }
            else
            {
                lastEZ = (smilesstr[i] == '/') ? 1 : -1;
            }
            continue;
        }

        if (smilesstr[i] >= '0' && smilesstr[i] <= '9')
        {
        	if (!ring_warned)
        	{
        		cout << "WARNING: Native support of SMILES is still in development. ";
        		cout << "Certain molecules containing rings might not render properly." << endl;
        		cout << "If possible, it is recommended to install OpenBabel (e.g. sudo apt-get install openbabel) ";
        		cout << "so the integration feature can be used and SMILES strings converted seamlessly to ";
        		cout << "the corresponding molecular structures." << endl;
        		ring_warned = true;
        	}
        	
            j = smilesstr[i] - 48;
            if (!numbered[j])
            {
                numbered[j] = preva;
                sqidx[j] = k-1;
                //cout << "+" << endl;;
                continue;
            }
            else
            {
            	//cout << " connecting..." << endl;
                // Rotate bonds to close the loop.
                bool allarom = true;
                Atom* aloop[256];
                for (l=sqidx[j]; l<k; l++)
                {
                    aloop[l-sqidx[j]] = sequence[l];
                    if (!seqarom[l]) allarom = false;
                }
                aloop[l-sqidx[j]] = 0;
                int ringsz = l-sqidx[j];
                
                /*
                if (!ring_atoms)
                {	ring_atoms = new Atom**[16];
                	int n; for (n=0; n<16; n++) ring_atoms[n] = nullptr;
                }
                
                if (!ring_aromatic)
                {	ring_aromatic = new bool[16];
                	int n; for (n=0; n<16; n++) ring_aromatic[n] = false;
                }
                
                if (!ring_atoms[ringcount])
                {	ring_atoms[ringcount] = new Atom*[ringsz+2];
                	int n; for (n=0; n<ringsz; n++) ring_atoms[ringcount][n] = nullptr;
                }*/
                
                add_ring(aloop);          

                // Also default all unmarked double bonds to cis.
                /*for (l=0; l<dbi; l++)
                {	if (!EZgiven[l] && EZatom0[l] && EZatom1[l])
                	{	Bond* lb = EZatom0[l]->get_bond_between(EZatom1[l]);
                		lb->can_rotate = true;
                		lb->flip_angle = M_PI;
                		lb = EZatom1[l]->get_bond_between(EZatom0[l]);
                		if (lb)
                		{	lb->can_flip = true;
                			lb->flip_angle = M_PI;
                		}
                	}
                }*/

                // If the ring has fewer than 5 members, or if it's aromatic, turn it into a regular polygon.
                /*if (ringsz<5 || allarom)
                {
                    Point ringcen;

                    SCoord v;
                    sequence[sqidx[j]]->arom_center = &ringcen;
                    sequence[sqidx[j]]->clear_geometry_cache();
                    SCoord* geob = sequence[sqidx[j]]->get_geometry_aligned_to_bonds();
                    int geon     = sequence[sqidx[j]]->get_geometry();
                    int gvi;
                    if (sqidx[j])		// First atom of ring is not first atom of molecule.
                        gvi = 0;		// The previous atom is outside the ring and forms the incoming SCoord.
                    else				// First atom of ring is first atom of molecule.
                        gvi = 2;		// Geometry vectors 0 and 1 will be ring bonds, so SCoord 2 forms the input.

                    if (geon == 3)
                    {
                        v = geob[gvi];
                    }
                    else
                    {
                        Point _4avg[2];
                        _4avg[0] = geob[gvi];
                        _4avg[1] = geob[3];

                        v = average_of_points(_4avg, 2);
                    }

                    if (allarom) card = 1.5;

                    preva->bond_to(numbered[j], card);

                    // TODO: Un-hardcode the bond lengths below.
                    v.r = polygon_radius(allarom?1.40:1.54, ringsz);
                    ringcen = sequence[sqidx[j]]->get_location().subtract(&v);

                    // Get the normal of v and the first two atoms.
                    Point _4norm[3];
                    _4norm[0] = sequence[sqidx[j]]->get_location().add(&v);
                    _4norm[1] = sequence[sqidx[j]]->get_location();
                    _4norm[2] = sequence[sqidx[j]+1]->get_location();

                    SCoord normal = compute_normal(&_4norm[0], &_4norm[1], &_4norm[2]);

                    for (l=0; l<ringsz; l++)
                    {
                        if (l)
                        {
                            Point lpt = rotate3D(&_4norm[1], &ringcen, &normal, M_PI*2/ringsz*l);
                            sequence[sqidx[j]+l]->move(&lpt);
                        }
                        sequence[sqidx[j]+l]->arom_center = &ringcen;
                        sequence[sqidx[j]+l]->clear_geometry_cache();
                    }
					
					ring_aromatic[ringcount-1] = true;
                    numbered[j] = 0;
                    card = 1;

                    continue;
                }*/

                float anomaly = close_loop(aloop, card);

                if (card) preva->bond_to(numbered[j], card);
                card = 1;
                // ring_aromatic[ringcount-1] = ring_is_aromatic(ringcount-1);

                numbered[j] = 0;

                continue;
            }
        }

        if (smilesstr[i] == '[')
        {
            bracket = true;
            continue;
        }

        if (bracket)
        {
            bool aromatic = false;

            if (smilesstr[i] == '@')
            {
                if (bracketed)
                {
                    bracketed->swap_chirality();
                    continue;
                }
                else throw 0xbade9c0d;
            }

            if (	(smilesstr[i] >= 'A' && smilesstr[i] <= 'Z')
                    ||
                    (smilesstr[i] >= 'a' && smilesstr[i] <= 'z')
               )
            {
                char esym[5];
                int ioff=1;

                esym[0] = smilesstr[i];
                esym[1] = smilesstr[i+1];
                esym[2] = 0;

                if (smilesstr[i+1] < 'a'
                        ||
                        smilesstr[i+1] > 'z'
                        ||
                        !Atom::Z_from_esym(esym)
                   )
                {
                    esym[1] = 0;
                    ioff = 0;
                }

                if (esym[0] >= 'a')
                {
                    aromatic = true;
                    esym[0] &= 0x5f;
                    if (prevarom) card=1.5;
                }

                if (Atom::Z_from_esym(esym))
                {
                    char aname[7];
                    if (!bracketed)
                    {
                        sprintf(aname, "%s%d", esym, atno++);
                        bracketed = add_atom(esym, aname, preva, card);
                        if (aromatic) bracketed->aromatize();
                        if (card == 2)
                        {
                            EZatom1[dbi++] = bracketed;
                        }
                        bracketed->dnh = true;
                        i += ioff;
                        ioff=0;
                    }
                    else
                    {
                        i += ioff;
                        ioff=0;
                        int plex = atoi(smilesstr+i+1);
                        if (!plex) plex=1;

                        for (l=0; l<plex; l++)
                        {
                            sprintf(aname, "%s%d", esym, atno++);
                            Atom* a = add_atom(esym, aname, bracketed, 1);
                            a->dnh = true;
                        }
                    }
                }
                else
                {
                    std::string str = smilesstr, estr = esym;
                    cout << "Bad element symbol " << esym << " (" << i << ")." << endl
                         << str.substr(0, i)
                         << "\x1b[4m" << str.substr(i, estr.length())
                         << "\x1b[0m" << str.substr(i+estr.length())
                         << endl;
                    throw 0xbade9c0d;
                }

                continue;
            }

            if (smilesstr[i] == '+')
            {
                if (bracketed)
                {
                    int plex = atoi(smilesstr+i+1);
                    if (!plex) plex=1;
                    bracketed->increment_charge(plex);
                }
                else throw 0xbade9c0d;
            }

            if (smilesstr[i] == '-')
            {
                if (bracketed)
                {
                    int plex = atoi(smilesstr+i+1);
                    if (!plex) plex=1;
                    bracketed->increment_charge(-plex);
                }
                else throw 0xbade9c0d;
            }

            if (smilesstr[i] == ']')
            {
                bracket = false;
                preva = bracketed;
                bracketed = 0;

                seqarom[k] = aromatic;
                sequence[k++] = preva;
                card = 1;
                prevarom = aromatic;

                continue;
            }

            continue;
        }

        char esym[5], aname[5];
        bool aromatic=false;

        sprintf(esym, "%c", smilesstr[i]);
        if (esym[0] >= 'a' && esym[0] <= 'z')
        {
            esym[0] &= 0x5f;
            aromatic = true;
            if (prevarom) card=1.5;
        }

        sprintf(aname, "%c%d", smilesstr[i], atno++);
        // cout << "Bonding new " << esym << " named " << aname << " to " << (preva?preva->name:"nothing") << " cardinality " << card << endl;
        Atom* a = add_atom(esym, aname, card ? preva : 0, card);
        if (aromatic) a->aromatize();
        if (card == 2)
        {
            EZatom1[dbi++] = a;
        }

        seqarom[k] = aromatic;
        sequence[k++] = a;

        preva = a;
        card = 1;
        prevarom = aromatic;
    }
    atoms[atcount]=0;

    return true;
}

/*
void Molecule::make_coplanar_ring(Atom** ring_members, int ringid)
{
    if (!ring_members) return;
    int i, j, l, ringsz;

    bool allarom = true;
    // cout << "Making coplanar ring of ";
    for (i=0; ring_members[i]; i++)
    {
    	// cout << ring_members[i]->name << " ";
        ringsz = i+1;
        Bond** ab = ring_members[i]->get_bonds();
        bool haspi = false;
        for (j=0; ab[j]; j++)
            if (ab[j]->cardinality > 1 && ab[j]->cardinality < 2) haspi = true;
        if (!haspi) allarom = false;
    }
    // cout << endl;

    if (ringsz<3) return;
    if (ringsz>6) return;

    SCoord normal;
    Point ringcen;

	if (!ring_members || !ring_members[0])
	{
		cout << "Notice: empty ring passed to Molecule::make_coplanar_ring()." << endl;
		return;
	}
    Bond** a0b = ring_members[0]->get_bonds();
    if (!a0b[0]->btom || !a0b[1]->btom)
    {
        cout << "Attempted to form coplanar ring with starting atom bonded to fewer than two other atoms." << endl;
        throw 0xbad12196;					// If you use your imagination, 12196 spells "ring".
    }

    Point A, B, C;
    A = ring_members[0]->get_location();
    B = a0b[0]->btom->get_location();
    C = a0b[1]->btom->get_location();

    normal = compute_normal(&A, &B, &C);
    while (!normal.r)
    {
    	C.x = frand(-1, 1);
    	C.y = frand(-1, 1);
    	C.z = frand(-1, 1);
    	normal = compute_normal(&A, &B, &C);
    }
    
    // TODO: Un-hardcode the bond lengths below.
    ringcen = A.subtract(&B);
    ringcen.scale(polygon_radius(allarom?1.40:1.54, ringsz));
    ringcen = A.add(ringcen);

    for (l=0; l<ringsz; l++)
    {
        if (l)
        {
            Point lpt = rotate3D(&A, &ringcen, &normal, M_PI*2/ringsz*l);
            ring_members[l]->move(&lpt);
        }
        ring_members[l]->arom_center = &ringcen;
        ring_members[l]->clear_geometry_cache();
        ring_members[l]->ring_member = max(1, ring_members[l]->ring_member);
        if (Huckel(ringid)) ring_members[l]->arom_ring_member = max(1, ring_members[l]->arom_ring_member);
        if (ring_members[l]->get_valence() != 4) ring_members[l]->aromatize();

        Bond* b2=0;
        int bgeo = ring_members[l]->get_geometry();
    	// cout << bgeo << ", checking " << ring_members[l]->name << "... " << flush;
        for (i=0; i<bgeo; i++)
        {
        	// cout << i << " ";
            b2 = ring_members[l]->get_bond_by_idx(i);
            if (!b2)
            {
            	// cout << "nullptr bond." << endl;
                continue;
            }
            if (!b2->btom)
            {
            	// cout << "nullptr btom." << endl;
                b2=0;
                continue;
            }
            for (j=0; j<ringsz; j++)
            {
            	if (ring_members[j] == b2->btom)
                {
            		// cout << "btom is part of the ring." << endl;
                    b2 = 0;
                    break;
                }
            }
            if (b2) break;
        }
        if (b2 && b2->btom)
        {
        	// cout << "found " << b2->btom->name << endl;
            Point ptnew = ring_members[l]->get_location().subtract(ringcen);
            ptnew.scale(InteratomicForce::covalent_bond_radius(ring_members[l], b2->btom, b2->cardinality));
            ptnew = ring_members[l]->get_location().add(ptnew);
            b2->btom->move_assembly(&ptnew, ring_members[l]);
        }
    }
}*/

float Molecule::fsb_lsb_anomaly(Atom* first, Atom* last, float lcard, float bond_length)
{
    // Last-should-be and first-should-be positions.
    SCoord lsbv = first->get_next_free_geometry(lcard);
    SCoord fsbv = last->get_next_free_geometry(lcard);
    lsbv.r = fsbv.r = bond_length;
    Point  lsb  = first->get_location().add(&lsbv);
    Point  fsb  = last->get_location().add(&fsbv);

    return first->get_location().get_3d_distance(fsb) + last->get_location().get_3d_distance(lsb);
}

float Molecule::close_loop(Atom** path, float lcard)
{
    Bond* rotables[65536];

    if (!path) return 0;

    int i, j, k=0;
    Atom* first, *last;

    first = path[0];
    if (!first) return 0;

    for (i=0; path[i]; i++)
    {
        Bond** b = path[i]->get_bonds();
        if (!b) continue;
        int geo = path[i]->get_geometry();

        for (j=0; j<geo; j++)
        {
            if (b[j]->btom
                &&
                (	b[j]->can_rotate
                    ||
                    b[j]->can_flip
                )
                &&
                strcmp(b[j]->btom->name, "N")
               ) rotables[k++] = b[j];
        }

        last = path[i];
    }
    rotables[k] = 0;

    if (last == first) return 0;
    last->mirror_geo = -1;
    int ringsize = k;

    float bond_length = InteratomicForce::covalent_bond_radius(first, last, lcard);

    if (ringsize < 5)
    {
        // TODO: Make equilateral ring, except accommodating any differences in bond length.
    }

    int iter;
    float anomaly = fsb_lsb_anomaly(first, last, lcard, bond_length);
    float iclash = get_internal_clashes();
    float oclash = iclash;
    float bondrot[ringsize];

    for (i=0; i<ringsize; i++) bondrot[i] = M_PI/2*randsgn();

    float allowance = 2.5;

    for (iter=0; iter<250+(20*ringsize); iter++)
    {
        for (i=0; rotables[i]; i++)
        {
        	float rr = rotables[i]->atom->get_location().get_3d_distance(rotables[i]->btom->get_location());
        	
        	if (fabs(rr-bond_length) > 0.01)
        	{
        		Point aloc = rotables[i]->atom->get_location();
        		Point bloc = rotables[i]->btom->get_location();
        		
        		bloc = bloc.subtract(aloc);
        		bloc.scale(bond_length);
        		bloc = bloc.add(aloc);
        		
        		rotables[i]->btom->move(bloc);
        	}
        	
            if (rotables[i]->rotate(bondrot[i]))
            {
                float newanom = fsb_lsb_anomaly(first, last, lcard, bond_length);
                float nclash = get_internal_clashes();
                if ((	newanom <= anomaly
                        ||
                        (rotables[i]->can_flip && (rand()%100) < 22)
                    )
                        &&
                        nclash <= allowance * iclash
                   )
                {
                    if (_DBGCLSLOOP) cout << "Anomaly was " << anomaly << " now " << newanom << ", keeping." << endl;
                    anomaly = newanom;
                    iclash = nclash;
                }
                else
                {
                    if (_DBGCLSLOOP) cout << "Anomaly was " << anomaly << " now " << newanom << ", reverting." << endl;
                    rotables[i]->rotate(-bondrot[i]);
                    bondrot[i] *= -0.6;
                }
            }
            else
            {
                bondrot[i] *= -0.5;
            }
        }
        allowance = (allowance-1)*.99+1;
    }

    for (i=0; path[i]; i++)
    {
        path[i]->clear_geometry_cache();
        if (_DBGCLSLOOP) cout << "Reset " << path[i]->name << " geometry." << endl;
    }

    for (i=0; rotables[i]; i++)
        rotables[i]->can_rotate = false;

    return anomaly;
}

#define _DBG_MOLBB 0
Atom** Molecule::get_most_bindable(int max_count)
{
    if (noAtoms(atoms)) return 0;

    int i, j=-1, k, l;
    float best[max_count+2];
    Atom** retval = new Atom*[max_count+2];
    
    for (k=0; k<max_count; k++)
    {
    	best[k]=0;
    	retval[k]=0;
	}
    
    for (i=0; atoms[i]; i++)
    {
    	if (atoms[i]->is_backbone) continue;
    	
    	#if _DBG_MOLBB
    	cout << "Probing atom " << atoms[i]->name << endl;
    	#endif
    	
        float score = 0;
        atoms[i]->clear_geometry_cache();
        
        for (k=0; retval[k] && k<max_count; k++)
        {
        	if (retval[k] == atoms[i])
        	{
				#if _DBG_MOLBB
				cout << "Atom is already in return array." << endl;
				#endif
        		goto _resume;
    		}
        	if (retval[k]->is_bonded_to(atoms[i]) && best[k] >= 10)
        	{
				#if _DBG_MOLBB
				cout << "Atom is bonded to " << retval[k]->name << ", already in return array." << endl;
				#endif
        		goto _resume;
    		}
        	if (retval[k]->shares_bonded_with(atoms[i]) && best[k] >= 10)
        	{
				#if _DBG_MOLBB
				cout << "Atom shares a bond with " << retval[k]->name << ", already in return array." << endl;
				#endif
        		goto _resume;
    		}
        }
        
        if (atoms[i]->get_charge() || atoms[i]->get_acidbase())
            score += 1000;

        if (atoms[i]->is_metal())
            score += 500;

        if (atoms[i]->is_thio())
            score += 20;
        else if (atoms[i]->is_polar())
            score += 100 * fabs(atoms[i]->is_polar());

        if (atoms[i]->is_pi())
            score += 50;

        if (!score) score += 5;		// van der Waals.
        
        #if _DBG_MOLBB
    	cout << "Score is " << score << endl;
    	#endif

		for (k=0; k<max_count; k++)
		{
			if (score > best[k])
		    {
		    	#if _DBG_MOLBB
				cout << "Score claims new #" << k << " spot." << endl;
				#endif
				
		    	for (l=max_count; l>k; l--)
		    	{
		    		best[l] = best[l-1];
		    		retval[l] = retval[l-1];
	    		}
		    	
		        best[k] = score;
		        retval[k] = atoms[i];
		        if (!k) j = i;
		        
		        break;
		    }
		    #if _DBG_MOLBB
			cout << "Score does not exceed previous #" << k << " spot of " << best[k] << endl;
			#endif
	    }
	    
	    _resume:
	    ;
    }
    retval[max_count] = 0;

    if (j < 0) return 0;
    else return retval;
}

#define DBG_BINDABLE 0
Atom** Molecule::get_most_bindable(int max_num, Atom* for_atom)
{
	if (!atoms) return nullptr;
	if (!for_atom) return nullptr;
	
	#if DBG_BINDABLE
	cout << "Molecule::get_most_bindable( " << max_num << ", " << for_atom->name << " )" << endl;
	#endif
	
	int mn2 = max_num+2;
	
	int i, j, k;
	float bb[mn2];
	Atom** bba = new Atom*[mn2];
	
	for (i=0; i<mn2; i++)
	{
		bb[i] = 0;
		bba[i] = nullptr;
	}
	
	for (i=0; atoms[i]; i++)
	{
		if (atoms[i]->is_backbone) continue;
		InteratomicForce** iff = InteratomicForce::get_applicable(atoms[i], for_atom);
		if (!iff) continue;
		
		float lbb = 0;
		for (j=0; iff[j]; j++)
		{
			float kj = iff[j]->get_kJmol();
			
			if (atoms[i]->get_charge()) kj *= fabs(atoms[i]->get_charge());
			if (for_atom->get_charge()) kj *= fabs(for_atom->get_charge());
			
			lbb += kj;
		}
		
		#if DBG_BINDABLE
		// Note: A competent programmer would not have to put these debugs all over the place every stinking function.
		cout << "Binding potential for " << atoms[i]->name << " " << lbb << " kJ/mol." << endl;
		#endif
		
		for (j=0; j<max_num; j++)
		{
			if (lbb > bb[j])
			{
				for (k=max_num-1; k>j; k--)
				{
					bba[k] = bba[k-1];
					bb[k] = bb[k-1];
				}
				
				bba[j] = atoms[i];
				bb[j] = lbb;
				
				#if DBG_BINDABLE
				cout << "Potential is better than previous result " << j << endl;
				#endif
				
				break;
			}
		}
	}
	
	#if DBG_BINDABLE
	cout << "Returning:" << endl << flush;
	for (i=0; i<max_num; i++)
	{
		cout << i << ": " << flush << bba[i] << flush << " ";
		if (bba[i]) cout << bba[i]->name << flush;
		cout << endl << flush;
	}
	cout << endl;
	#endif
	
	return bba;
}

Point Molecule::get_bounding_box() const
{
    if (noAtoms(atoms)) return 0;

    int i;
    float xmax=0, ymax=0, zmax=0;

    for (i=0; atoms[i]; i++)
    {
        SCoord v(atoms[i]->get_location());
        float r = atoms[i]->get_vdW_radius();
        Point pt(v);

        pt.x = fabs(pt.x)+r;
        pt.y = fabs(pt.y)+r;
        pt.z = fabs(pt.z)+r;

        if (pt.x > xmax) xmax = pt.x;
        if (pt.y > ymax) ymax = pt.y;
        if (pt.z > zmax) zmax = pt.z;
    }

    Point pt(xmax, ymax, zmax);
    return pt;
}

bool Molecule::ring_is_coplanar(int ringid)
{
	if (!rings) return false;
	return rings[ringid]->is_coplanar();
}

bool Molecule::ring_is_aromatic(int ringid)
{
	if (!rings) return false;
	return rings[ringid]->get_type() == AROMATIC;
}

Point Molecule::get_ring_center(int ringid)
{
	if (!rings) return Point(0,0,0);
	return rings[ringid]->get_center();
}

SCoord Molecule::get_ring_normal(int ringid)
{
	if (!rings) return SCoord(0,0,0);
	return rings[ringid]->get_normal();
}

Atom** Molecule::get_ring_atoms(int ringid)
{
	if (!rings) return nullptr;
	return rings[ringid]->get_atoms();
}

int Molecule::get_ring_num_atoms(int ringid)
{
	if (!rings) return 0;
	return rings[ringid]->get_atom_count();
}


void Molecule::recenter_ring(int ringid, Point new_ring_cen)
{
	if (!rings) return;
	Point old_ring_cen = get_ring_center(ringid);
	SCoord motion = new_ring_cen.subtract(old_ring_cen);
	int i;
	Atom** ring_atoms = rings[ringid]->get_atoms();
	for (i=0; ring_atoms[i]; i++)
		ring_atoms[i]->move_rel(&motion);
	
	delete[] ring_atoms;
}

void Molecule::rotate_ring(int ringid, Rotation rot)
{
	if (!rings) return;
	Point origin = get_ring_center(ringid);
	int i;
	Atom** ring_atoms = rings[ringid]->get_atoms();
	for (i=0; ring_atoms[i]; i++)
	{
		Point aloc = ring_atoms[i]->get_location();
		aloc = rotate3D(&aloc, &origin, &rot);
		ring_atoms[i]->move(aloc);
	}
	
	delete[] ring_atoms;
}

int Molecule::get_num_rings()
{
	if (!rings) return 0;
	int i;
	for (i=0; rings[i]; i++);	// Get count.
	return i;
}

bool Molecule::in_same_ring(Atom* a, Atom* b)
{
	int i;
	Ring** r = a->get_rings();
	if (!r) return false;
	for (i=0; r[i]; i++)
		if (b->is_in_ring(r[i]))
		{
			delete[] r;
			return true;
		}
	
	delete[] r;
	return false;
}

float Molecule::get_atom_error(int i, LocatedVector* best_lv)
{
	int j;
	float error = 0;
	Point bloc;
	Atom* btom;
	Bond** b;
	int g, bg;
	float b_bond_angle;
	LocatedVector lv;
	
	// Get the atom's zero-index bonded atom. Call it btom (because why not overuse a foolish pun?).
    b = atoms[i]->get_bonds();
    if (!b) return error;
	btom = b[0]->btom;
	if (!btom) return error;
	
	bloc = btom->get_location();
	g = atoms[i]->get_geometry();
	bg = btom->get_geometry();
	b_bond_angle = btom->get_geometric_bond_angle();
	
	// Make an imaginary sphere around btom, whose radius equals the optimal bond distance.
	lv.origin = bloc;
	lv.r = InteratomicForce::covalent_bond_radius(atoms[i], btom, b[0]->cardinality);
	float thstep = fiftyseventh*5;
	float besttheta = 0, bestphi = 0, bestscore = -1e9;
	for (lv.theta = -square; lv.theta <= square; lv.theta += thstep)
	{
		float phstep = M_PI/(20.0*(sin(lv.theta) + 1));
		for (lv.phi = 0; lv.phi < (M_PI*2); lv.phi += phstep)
		{
			// At many points along the sphere, evaluate the goodness-of-fit as a function of:
			// Success in conforming to btom's geometry;
			// Success in avoiding clashes with atoms not bonded to self or btom;
			// Success in maintaining optimal binding distances to own bonded atoms.
			// Later, we'll test edge cases where bond strain distorts the usual angles.
			float score = 0;
			
			score -= _SANOM_BOND_ANGLE_WEIGHT*btom->get_bond_angle_anomaly(lv, atoms[i]);
			
			// Avoid clashes with strangers.
			for (j=0; atoms[j]; j++)
			{
				if (j == i) continue;
				if (atoms[j]->is_bonded_to(atoms[i])) continue;
				
				float r = atoms[j]->get_location().get_3d_distance(lv.to_point());
				score -= _SANOM_CLASHES_WEIGHT/fabs(r+0.000000001);
			}
			
			// Seek optimal bond radii.
			for (j=1; b[j]; j++)
			{
				if (!b[j]->btom) continue;
				float optimal = InteratomicForce::covalent_bond_radius(atoms[i], b[j]->btom, b[j]->cardinality);
				float r = b[j]->btom->get_location().get_3d_distance(lv.to_point());
				
				score -= _SANOM_BOND_RAD_WEIGHT * fabs(optimal-r);
			}
			
			if (score > bestscore)
			{
				besttheta = lv.theta;
				bestphi = lv.phi;
				bestscore = score;
			}
		}
	}
	
	if (best_lv)
	{
		best_lv->origin = lv.origin;
		best_lv->r = lv.r;
		best_lv->theta = besttheta;
		best_lv->phi = bestphi;
	}
	
	lv = (SCoord)atoms[i]->get_location().subtract(bloc);
	lv.origin = bloc;
	
	error += _SANOM_BOND_ANGLE_WEIGHT*btom->get_bond_angle_anomaly(lv, atoms[i]);
	
	for (j=0; atoms[j]; j++)
	{
		if (j == i) continue;
		if (atoms[j]->is_bonded_to(atoms[i])) continue;
		
		float r = atoms[j]->get_location().get_3d_distance(lv.to_point());
		error += _SANOM_CLASHES_WEIGHT/fabs(r+0.000000001);
	}
	 
	for (j=1; b[j]; j++)
	{
		if (!b[j]->btom) continue;
		float optimal = InteratomicForce::covalent_bond_radius(atoms[i], b[j]->btom, b[j]->cardinality);
		float r = b[j]->btom->get_location().get_3d_distance(lv.to_point());
		
		error += _SANOM_BOND_RAD_WEIGHT * fabs(optimal-r);
	}
	
	return error+bestscore;
}


#define _DEV_FIX_MSTRUCT 0
float Molecule::correct_structure(int iters)
{
    if (noAtoms(atoms)) return 0;
    int iter, i, j, k;
    Point zero(0,0,0);
	float error = 0;
	Point aloc, bloc;
	Atom* btom;
	Bond** b;
	int g, bg;
	float b_bond_angle;
	LocatedVector lv;
    
	#if _DEV_FIX_MSTRUCT
    // TODO
    if (ringcount)
    {
    	for (i=0; i<ringcount; i++)
    	{
    		if (ring_aromatic[i] || get_ring_num_atoms(i) < 6)
    		{
				make_coplanar_ring(ring_atoms[i], i);
    		}
    	}
    }
    else return 0;			// Non-ring structures work fine. The buggy algorithm is when rings are involved.
	
    for (iter=0; iter<iters; iter++)
    {
    	error = 0;
        for (i=0; atoms[i]; i++)
        {
			// Get the atom's zero-index bonded atom. Call it btom (because why not overuse a foolish pun?).
			b = atoms[i]->get_bonds();
			if (!b) return error;
			btom = b[0]->btom;
			if (!btom) return error;
			
			// TODO
    		if (atoms[i]->num_rings() && atoms[i]->is_pi() && in_same_ring(atoms[i], btom)) continue;
	
			error += get_atom_error(i, &lv);
			if (!iter) continue;
			
			// Once a "best fit" point in space is found, move there.
			if (atoms[i]->num_rings() && atoms[i]->is_pi())
			{
				for (j=0; j<ringcount; j++)
				{
					for (k=0; ring_atoms[j][k]; k++)
					{
						if (ring_atoms[j][k] == atoms[i])
						{
							// Get distance from atom to ring center.
							Point rcen = get_ring_center(j);
							Point aloc = atoms[i]->get_location();
							float rad = rcen.get_3d_distance(aloc);
							
							// Center ring at combined distance from btom.
							LocatedVector lvr = lv;
							lvr.r += rad;
							recenter_ring(j, lvr.to_point());
							
							// Rotate ring, and all assemblies, to align atom to btom.
							Point atarget = lv.to_point();
							aloc = atoms[i]->get_location();
							rcen = get_ring_center(j);
							Rotation rot = align_points_3d(&aloc, &atarget, &rcen);
							rotate_ring(j, rot);
							
							break;
						}
					}
				}
			}
			else
			{
				/*Point pt = lv.to_point();
				atoms[i]->move_assembly(&pt, btom);*/
				atoms[i]->move(lv.to_point());
			}
		}
		// cout << error << " " << atcount << endl;
		// if (error <= 10.0*atcount) break;
	}
	#endif

	error = 0;
    for (i=0; atoms[i]; i++)
    {
    	error += get_atom_error(i, &lv);
    }

	return error;
}












