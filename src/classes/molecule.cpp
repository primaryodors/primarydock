
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

Molecule::Molecule(const char* lname)
{
    name = new char[strlen(lname)+1];
    strcpy(name, lname);

    atoms = nullptr;
    smiles = nullptr;
    ring_atoms = nullptr;
    atcount = ringcount = 0;
    reset_conformer_momenta();
    rotatable_bonds = nullptr;
}

Molecule::Molecule()
{
    atoms = nullptr;
    smiles = nullptr;
    ring_atoms = nullptr;
    atcount = ringcount = 0;
    reset_conformer_momenta();
    rotatable_bonds = 0;
    paren = nullptr; // not sure what a good default is here, but it was not initialized (warning from clang)
}

Molecule::Molecule(const char* lname, Atom** collection)
{
    if (!collection)
    {
        cout << "Temporary molecule creation attempted from nullptr atom pointer array." << endl;
        throw 0xbadca22;
    }

    int i, j;
    for (i=0; collection[i]; i++);
    atoms = new Atom*[i+1];
    for (j=0; j<i; j++)
        atoms[j] = collection[j];
    atoms[i] = 0;
    atcount = i;

    name = new char[strlen(lname)+1];
    strcpy(name, lname);

    smiles = nullptr;
    ring_atoms = nullptr;
    ringcount = 0;
    reset_conformer_momenta();
    rotatable_bonds = 0;
}

Molecule::~Molecule()
{
    /*int i;
    if (atoms)
    {	for (i=0; atoms[i]; i++) delete atoms[i];
    	delete[] atoms;
    }
    if (name) delete[] name;
    if (smiles) delete smiles;
    if (ring_atoms)
    {	for (i=0; ring_atoms[i]; i++) delete[] ring_atoms[i];
    	delete[] ring_atoms;
    }
    if (ring_aromatic) delete[] ring_aromatic;*/
}

void Molecule::delete_atom(Atom* a)
{
	if (!a) return;
	int i, j;
	
	if (atoms)
	{
		for (i=0; atoms[i]; i++)
		{
			if (atoms[i] == a)
			{
				a->unbond_all();
				for (j=i+1; atoms[j]; j++) atoms[j-1] = atoms[j];
				atoms[j-1] = nullptr;
				rotatable_bonds = nullptr;
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
    rotatable_bonds = 0;
}

Atom* Molecule::add_atom(const char* elemsym, const char* aname, const Point* location, Atom* bond_to, const float bcard)
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

Atom* Molecule::add_atom(const char* elemsym, const char* aname, Atom* bondto, const float bcard)
{
    if (!bondto || !bcard)
    {
        Point pt;
        return add_atom(elemsym, aname, &pt, bondto, bcard);
    }

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

void Molecule::hydrogenate()
{
    if (!atoms) return;
    int i, j;

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_metal()) continue;					// No hydrides.
        if (atoms[i]->dnh) continue;

        int valence = atoms[i]->get_valence();
        if (valence > 4) valence = 8 - valence;

        //cout << atoms[i]->name << " has valence " << valence;

        int bcardsum = 0;

        Bond** aib = atoms[i]->get_bonds();

        if (aib)
        {
            for (j=0; aib[j]; j++)
            {
                if (aib[j]->btom) bcardsum += aib[j]->cardinality;
            }
        }

        //cout << " minus existing bonds " << bcardsum ;

        bcardsum -= atoms[i]->get_charge();

        //cout << " given charge makes " << bcardsum << endl;

        for (j=bcardsum; j<valence; j++)
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

Atom* Molecule::get_atom(const char* aname) const
{
    if (!atoms) return 0;

    int i;
    for (i=0; atoms[i]; i++)
        if (!strcmp(atoms[i]->name, aname)) return atoms[i];

    return 0;
}

Point Molecule::get_atom_location(char const * aname)
{	if (!atoms)
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
    if (!atoms) return 0;

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
    if (!atoms) return 0;

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

int Molecule::from_sdf(const char* sdf_dat)
{
    if (!sdf_dat) return 0;
    const char* lines[8192];
    int i,j=0,lncount;

    immobile = false;

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

        if (!fields[0] || !fields[1]) break;
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

            atoms = new Atom*[na+1];
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
        }
        else
        {
            int a1i = atoi(fields[0]);
            int a2i = atoi(fields[1]);

            if (!a1i || !a2i) break;
            atoms[a1i-1]->bond_to(atoms[a2i-1], atoi(fields[2]));
        }

        if (fields) delete[] fields;
    }
    atoms[atcount] = 0;
    if (fields) delete[] fields;

    identify_rings();
    identify_acidbase();
    minclash = get_internal_clashes();
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

        delete[] fields;
    }

    return added;
}

int Molecule::get_bond_count(bool unidirectional) const
{
    int i, j, bc=0;

    if (!atoms) return 0;
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
    if (!atoms) return -1;

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
    fprintf(os, "  podock-%02d%02d%02d%02d%02d%02d3D\n", gmtm->tm_year % 100, gmtm->tm_mon+1, gmtm->tm_mday,
            gmtm->tm_hour, gmtm->tm_min, gmtm->tm_sec
           );

    fprintf(os, "https://github.com/primaryodors/podock\n");

    int ac, bc, chargeds=0;
    ac = get_atom_count();
    bc = get_bond_count(true);

    int i, j, k, l;
    Atom* latoms[65536];
    Bond* lbonds[65536];

    if (atoms)
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
        const char* esym = latoms[i]->get_elem_sym();
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

bool Molecule::ring_is_coplanar(int ringid)
{
    int ringsz;
    if (!ring_atoms) return false;
    if (!ring_atoms[ringid]) return false;
    for (ringsz = 0; ring_atoms[ringid][ringsz]; ringsz++);

    if (ringsz < 4) return true;

    int i;
    float anomaly;
    for (i=3; i<ringsz; i++)
    {
        anomaly = are_points_planar(ring_atoms[ringid][0]->get_location(),
                                    ring_atoms[ringid][1]->get_location(),
                                    ring_atoms[ringid][2]->get_location(),
                                    ring_atoms[ringid][i]->get_location()
                                   );
    }

    return (anomaly < 0.1);
}

bool Molecule::ring_is_aromatic(int ringid)
{
    if (!ring_aromatic) return false;
    if (ringid >= ringcount) return false;
    if (ringid < 0) return false;
    return ring_aromatic[ringid];
}

Point Molecule::get_ring_center(int ringid)
{
    int ringsz;
    Point nothing;
    if (!ring_atoms) return nothing;
    if (!ring_atoms[ringid]) return nothing;
    for (ringsz = 0; ring_atoms[ringid][ringsz]; ringsz++);

    Point* alocs = new Point[ringsz];
    int i;
    for (i=0; i<ringsz; i++)
    {
        alocs[i] = ring_atoms[ringid][i]->get_location();
    }

    Point retval = average_of_points(alocs, ringsz);
    delete[] alocs;
    return retval;
}

SCoord Molecule::get_ring_normal(int ringid)
{
    SCoord nothing;
    if (!ring_atoms) return nothing;
    if (!ring_atoms[ringid]) return nothing;

    if (!ring_aromatic[ringid] && !ring_is_coplanar(ringid))
    {
        // TODO: Some kind of best-fit algorithm using an imaginary circle.
        return nothing;
    }

    Point a, b, c;
    a = ring_atoms[ringid][0]->get_location();
    b = ring_atoms[ringid][1]->get_location();
    c = ring_atoms[ringid][2]->get_location();

    return compute_normal(&a, &b, &c);
}

Atom** Molecule::get_ring_atoms(int ringid)
{
    if (!ring_atoms) return 0;
    if (!ring_atoms[ringid]) return 0;

    int ringsz;
    for (ringsz = 0; ring_atoms[ringid][ringsz]; ringsz++);

    Atom** retval = new Atom*[ringsz+1];

    int i;
    for (i=0; i<=ringsz; i++) retval[i] = ring_atoms[ringid][i];

    return retval;
}

int Molecule::identify_rings()
{
    Atom **ringstmp[256], *a;
    int chainlen[256];
    bool is_ring[256];
    int found_rings=0, chains=0, cnvchain, active, i, j, k, l, m, n, p;
    Atom *cnva, *cnvb;
    Atom *ra, *rb;

    ringcount = 0;
    ring_atoms = new Atom**[atcount];
    ring_aromatic = new bool[atcount];

    // Start at any ato Mark it "used".
    a = atoms[0];
    a->used = true;

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

        for (i=0; i<chains; i++)
            if (ringstmp[i][0])
            {
                //cout << "Begin chain " << i << ": "; Atom::dump_array(ringstmp[i]); cout << endl;
            }

        for (i=0; i<chains; i++)
        {
            if (!ringstmp[i] || !chainlen[i]) continue;

            cnva = cnvb = 0;
            cnvchain = 0;

            if (!(a = ringstmp[i][chainlen[i]-1])) continue;
            //cout << "Chain " << i << " ends in atom " << a->name << endl;

            if (!cnva)
            {
                for (l=0; l<chains; l++)
                {
                    //cout << "Preparing to check chain " << i << " against chain " << l << " for bond convergence...\n";
                    if (l != i && chainlen[l] && ringstmp[l][chainlen[l]-1])
                    {
                        // cout << "Checking chain " << i << " against chain " << l << " for bond convergence...\n";
                        // If two chains converge on a single atom, this is your one converging ato
                        if (ringstmp[l][chainlen[l]-1] == a)
                        {
                            cnva = a;
                            cnvb = 0;
                            //cout << "****** Chains " << i << " and " << l << " converge on atom " << cnva->name << " ******" << endl;
                            cnvchain = l;
                            break;
                        }
                        // If two chains converge on two atoms that are bonded to each other, then there are two converging atoms.
                        else if (ringstmp[l][chainlen[l]-1]->is_bonded_to(a))
                        {
                            cnva = a;
                            cnvb = ringstmp[l][chainlen[l]-1];
                            //cout << "***** Chains " << i << " and " << l << " converge on atoms " << cnva->name << " and " << cnvb->name << " *****" << endl;
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
                // until reach the/a converging ato This is a RING; add it to the list.
                ring_atoms[ringcount] = new Atom*[chainlen[i] + chainlen[cnvchain]];
                n = 0;
                // cout << "Counting backwards from " << chainlen[i]-1 << endl;
                for (m = chainlen[i]-1; m >= 0; m--)
                {
                    ring_atoms[ringcount][n++] = ringstmp[i][m];
                    // cout << "m: " << ringstmp[i][m]->name << " ";
                    l = in_array(reinterpret_cast<void*>(ringstmp[i][m]),
                                 reinterpret_cast<void**>(ringstmp[cnvchain])
                                );
                    if (l >= 0)
                    {
                        l++;
                        // cout << "\nCounting forwards from " << l << endl;
                        for (; ringstmp[cnvchain][l]; l++)
                        {
                            // cout << "l: " << ringstmp[cnvchain][l]->name << " ";
                            ring_atoms[ringcount][n++] = ringstmp[cnvchain][l];
                            // cout << "Building " << n << " membered ring: "; Atom::dump_array(ring_atoms[ringcount]); cout << endl;
                            if (n >= 3 && (ringstmp[cnvchain][l] == cnva || ringstmp[cnvchain][l] == cnvb))
                            {
                                // cout << "Found " << n << " membered ring: "; Atom::dump_array(ring_atoms[ringcount]); cout << endl;
                                ring_atoms[ringcount][n] = 0;
                                bool cp = ring_is_coplanar(ringcount);
                                if (cp)
                                {
                                    //cout << "Ring is coplanar." << endl;

                                    int pi_e = 0, pi_eo = 0;

                                    for (p=0; ra = ring_atoms[ringcount][p]; p++)
                                    {
                                        // Count the number of double bonds in the ring. For each double bond, count 2 pi electrons.
                                        rb = ring_atoms[ringcount][p ? (p-1) : (n-1)];
                                        int card = ra->is_bonded_to(rb);

                                        if (card == 2) pi_e += 2;

                                        // If there is a non-pi tetrel in the ring, without a negative charge, the ring is not aromatic.
                                        int val = ra->get_valence();
                                        if (val == 4 && !ra->is_pi() && (ra->get_charge() >= 0)) goto _not_aromatic;

                                        // If there is a pnictogen in the ring, count it as two optional pi electrons.
                                        int geo = ra->get_geometry();
                                        if (val == 3 && geo >= 4) pi_eo += 2;

                                        // If there is a chalcogen in the ring, count two more optional pi electrons.
                                        if (val == 2 && geo >= 4) pi_eo += 2;
                                    }

                                    // If the total number of pi electrons = 4n+2, then the ring is aromatic.
                                    pi_eo += pi_e;
                                    for (p = pi_e; p <= pi_eo; p += 2)
                                    {
                                        if (!((p-2)&0x3))
                                        {
                                            //cout << "Ring is aromatic." << endl;
                                            ring_aromatic[ringcount] = true;
                                            break;
                                        }
                                    }
                                }
_not_aromatic:

#if _ALLOW_FLEX_RINGS
                                if (cp || ring_aromatic[ringcount])
#else
                                if (1)
#endif
                                {
                                    for (p=0; ra = ring_atoms[ringcount][p]; p++)
                                    {
                                        rb = ring_atoms[ringcount][p ? (p-1) : (n-1)];
                                        int card = ra->is_bonded_to(rb);

                                        Bond* ab = ra->get_bond_between(rb);
                                        ab->can_rotate = false;
                                        ab = rb->get_bond_between(ra);
                                        ab->can_rotate = false;
                                    }
                                }

                                ringcount++;
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
                            //cout << "Another new chain " << chains << ": "; Atom::dump_array(ringstmp[chains]); cout << endl;
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
                            //cout << "Chain " << i << ": "; Atom::dump_array(ringstmp[i]); cout << endl;
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
    for (i=0; ringstmp[i]; i++) delete ringstmp[i];

    for (i=0; atoms[i]; i++) atoms[i]->used = false;

    // Return the number of rings found.
    return found_rings;
}

void Molecule::identify_acidbase()
{
    if (!atoms) return;

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
            if (c = atoms[i]->is_bonded_to("C"))
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
    if (!atoms) return 0;
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

    // cout << bonds << " rotatable bond(s)." << endl;
    rotatable_bonds = new Bond*[bonds+1];
    for (i=0; i<=bonds; i++) rotatable_bonds[i] = btemp[i];
    rotatable_bonds[bonds] = 0;

    return rotatable_bonds;
}

// TODO: There has to be a better way.
Bond** AminoAcid::get_rotatable_bonds()
{
    // cout << name << " AminoAcid::get_rotatable_bonds()" << endl << flush;
    // Return ONLY side chain bonds, from lower to higher Greek. E.g. CA-CB but NOT CB-CA.
    // Exclude CA-N and CA-C as these will be managed by the Protein class.
    if (!atoms) return 0;
    if (rotatable_bonds) return rotatable_bonds;
    if (aadef && aadef->proline_like)
    {
    	// cout << "Proline-like! No rotbonds!" << endl;
    	return nullptr;
	}
    Bond* btemp[65536];
    
    int i,j, bonds=0;
    if (aadef && aadef->aabonds)
    {
    	for (i=0; aadef->aabonds[i]; i++)
    	{
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
    					cout << "Warning: No bond between " << la->residue << ":" << la->name
							 << " and " << aadef->aabonds[i]->bname
							 << endl << flush;
					}
    				else
    				{
						lb->can_rotate =
							aadef->aabonds[i]->can_rotate;
						if (!la->is_backbone
							&&
							la->get_Z() > 1
				            &&
				            greek_from_aname(la->name) == (greek_from_aname(lb->btom->name)+1)
		                   )
		                {
						    btemp[bonds++] = lb->btom->get_bond_between(la);
						    btemp[bonds] = 0;
						}
				    }
				}
    		}
    	}
    	
    	goto _found_aadef;
    }

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
                btemp[bonds++] = lb[j];
                btemp[bonds] = 0;
            }
            else lb[j]->can_rotate = false;
        }
        if (lb) delete[] lb;
    }

	_found_aadef:
    Bond** retval = new Bond*[bonds+8];
    for (i=0; i<=bonds; i++) retval[i] = btemp[i];
    retval[i] = 0;
    rotatable_bonds = retval;

    return retval;
}

Bond** Molecule::get_all_bonds(bool unidirectional)
{
    if (!atoms) return 0;
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

    if (!atoms) return 0;

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
            if (r < avdW + bvdW)
            {
                float lclash = sphere_intersection(avdW, bvdW, r);
                clash += lclash;
                // cout << atoms[i]->name << " clashes with " << atoms[j]->name << " by " << lclash << " cu. A." << endl;
            }
        }
    }

    if (clash < minclash) minclash = clash;

    return clash-minclash;
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

    if (!atoms) return 0;
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
    if (!atoms) return;
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
    if (!atoms) return;
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

Point Molecule::get_barycenter() const
{
    if (!atoms)
    {
        Point pt;
        return pt;
    }

    Point locs[atcount];
    int i;

    for (i=0; i<atcount; i++) locs[i] = atoms[i]->get_location();

    return average_of_points(locs, atcount);
}

void Molecule::recenter(Point nl)
{
    Point loc = get_barycenter();
    Point rel = nl.subtract(&loc);
    SCoord v(&rel);
    move(v);
}

void Molecule::rotate(SCoord* SCoord, float theta)
{
    if (!atoms) return;
    // cout << name << " Molecule::rotate()" << endl;

    Point cen = get_barycenter();

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
    if (!atoms) return;

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
    if (!atoms) return 0;
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

float Molecule::get_intermol_binding(Molecule** ligands)
{
    if (!ligands) return 0;
    if (!ligands[0]) return 0;
    int i, j, l;
    float kJmol = 0;
    kJmol -= get_internal_clashes();

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
                    if (!shielded(atoms[i], ligands[l]->atoms[j])
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
    if (!atoms) return;

    minclash = 0;
    int i, j, iter;
    float clash = get_internal_clashes();

    if (!clash) return;		// Already zero, nothing to decrease to.

    Bond** b = get_rotatable_bonds();
    if (!b) return;
    if (!b[0]) return;		// No bonds to rotate.

    int numrb = 0;
    for (i=0; b[i]; i++) numrb = i;

    float angle[numrb];
    for (i=0; i<numrb; i++) angle[i] = M_PI;

    for (iter=0; iter<500; iter++)
    {
        for (i=0; i<numrb; i++)
        {
            b[i]->rotate(angle[i]);
            minclash = 0;
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

    minclash = clash;
}

void Molecule::intermol_conform(Molecule* ligand, int iters)
{
    Molecule* ligands[4];
    ligands[0] = ligand;
    ligands[1] = nullptr;
    intermol_conform(ligands, iters);
}

void Molecule::intermol_conform(Molecule** ligands, int iters)
{
    Molecule** m = 0;
    intermol_conform(ligands, iters, m);
}

void Molecule::intermol_conform_norecen(Molecule* ligand, int iters, AminoAcid** avcw)
{
    Molecule* ligands[4];
    ligands[0] = ligand;
    ligands[1] = nullptr;
    intermol_conform_norecen(ligands, iters, avcw);
}

void Molecule::intermol_conform_norecen(Molecule** ligands, int iters, AminoAcid** avcw)
{
    Molecule* aa[256];
    int i;
    for (i=0; avcw[i] && i<256; i++) aa[i] = reinterpret_cast<Molecule*>(avcw[i]);
    aa[i] = 0;

    intermol_conform_norecen(ligands, iters, aa);
}

void Molecule::intermol_conform(Molecule* ligand, int iters, Molecule** avcw)
{
    Molecule* ligands[4];
    ligands[0] = ligand;
    ligands[1] = nullptr;
    intermol_conform(ligands, iters, avcw);
}

void Molecule::intermol_conform(Molecule** ligands, int iters, Molecule** avcw)
{
    int iter, i, j;
    Point xpole(1,0,0),
          ypole(0,1,0),
          zpole(0,0,1);
    SCoord xv(&xpole),
           yv(&ypole),
           zv(&zpole);
    // cout << "Conforming " << name << endl << flush;

    // return intermol_conform_norecen(ligand, iters, avcw);

    float bind = get_intermol_binding(ligands);
    float bind1;
    if (avcw) bind -= get_intermol_clashes(avcw);

    for (iter=0; iter<iters; iter++)
    {
        if (echo_iters) cout << "Iteration " << iter << ": " << bind << endl;
        if (!immobile)
        {
            // Entire molecule location.
            Point lmpt(lmx,0,0);
            move(lmpt);
            bind1 = get_intermol_binding(ligands);
            if (avcw) bind -= get_intermol_clashes(avcw);

            /*if (fabs(lmx) > 0.5) lmx *= 0.9;
            if (fabs(lmy) > 0.5) lmy *= 0.9;
            if (fabs(lmz) > 0.5) lmz *= 0.9;*/
            if (fabs(lmx) > _def_lin_momentum) lmx = _def_lin_momentum*sgn(lmx);
            if (fabs(lmy) > _def_lin_momentum) lmy = _def_lin_momentum*sgn(lmy);
            if (fabs(lmz) > _def_lin_momentum) lmz = _def_lin_momentum*sgn(lmz);

            if (bind1 < bind)
            {
                lmpt.x *= -1;
                move(lmpt);
                lmx *= -0.75;
            }
            else
            {
                bind = bind1;
                if (fabs(lmx)<0.25) lmx *= 1.1;
            }

            lmpt.x = 0;
            lmpt.y = lmy;
            float oldy = get_barycenter().y;
            move(lmpt);
            bind1 = get_intermol_binding(ligands);
            if (avcw) bind -= get_intermol_clashes(avcw);

            if (bind1 < bind)
            {
                lmpt.y *= -1;
                move(lmpt);
                lmy *= -0.75;
            }
            else
            {
                float newy = get_barycenter().y;
                /*cout << iter << " Moved " << name << " Y" << (lmy>0 ? "+" : "") << lmy
                	 << "; Y center was " << oldy << " now " << newy
                	 << "; binding was " << bind << " now " << bind1 << endl;*/
                bind = bind1;
                if (fabs(lmy)<0.25) lmy *= 1.1;
            }

            lmpt.y = 0;
            lmpt.z = lmz;
            move(lmpt);
            bind1 = get_intermol_binding(ligands);
            if (avcw) bind -= get_intermol_clashes(avcw);

            if (bind1 < bind)
            {
                lmpt.z *= -1;
                move(lmpt);
                lmz *= -0.75;
            }
            else
            {
                bind = bind1;
                if (fabs(lmz)<0.25) lmz *= 1.1;
            }


            // Entire molecule rotation.
            rotate(&xv, amx);
            bind1 = get_intermol_binding(ligands);
            if (avcw) bind -= get_intermol_clashes(avcw);

            if (bind1 < bind)
            {
                rotate(&xv, -amx);
                amx *= -(0.5 + 0.5 / (bind/10+1));
            }
            else
            {
                bind = bind1;
                if (fabs(amx)<M_PI/4) amx *= 1.5;
            }

            rotate(&yv, amy);
            bind1 = get_intermol_binding(ligands);
            if (avcw) bind -= get_intermol_clashes(avcw);

            if (bind1 < bind)
            {
                rotate(&yv, -amy);
                amy *= -(0.5 + 0.5 / (bind/10+1));
            }
            else
            {
                bind = bind1;
                if (fabs(amy)<M_PI/4) amy *= 1.5;
            }

            rotate(&zv, amz);
            bind1 = get_intermol_binding(ligands);
            if (avcw) bind -= get_intermol_clashes(avcw);

            if (bind1 < bind)
            {
                rotate(&zv, -amz);
                amz *= -(0.5 + 0.5 / (bind/10+1));
                // cout << amz << " ";
            }
            else
            {
                bind = bind1;
                if (fabs(amz)<M_PI/4) amz *= 1.5;
            }
        }

        intermol_conform_flexonly(ligands, iters, avcw, bind);
    }
}

void Molecule::intermol_conform_norecen(Molecule* ligand, int iters, Molecule** avcw)
{
    Molecule* ligands[4];
    ligands[0] = ligand;
    ligands[1] = nullptr;
    intermol_conform_norecen(ligands, iters, avcw);
}

void Molecule::intermol_conform_norecen(Molecule** ligands, int iters, Molecule** avcw)
{
    intermol_conform_norecen(ligands, iters, avcw, nanf(""));
}

void Molecule::intermol_conform_norecen(Molecule* ligand, int iters, Molecule** avcw, float lastbind)
{
    Molecule* ligands[4];
    ligands[0] = ligand;
    ligands[1] = nullptr;
    intermol_conform_norecen(ligands, iters, avcw, lastbind);
}

void Molecule::intermol_conform_norecen(Molecule** ligands, int iters, Molecule** avcw, float lastbind)
{
    int iter, i, j;
    Point xpole(1,0,0),
          ypole(0,1,0),
          zpole(0,0,1);
    SCoord xv(&xpole),
           yv(&ypole),
           zv(&zpole);
    // cout << "Conforming " << name << endl << flush;

    float bind;
    bind = isnan(lastbind) ? get_intermol_binding(ligands) : lastbind;
    float bind1;
    if (avcw) bind -= get_intermol_clashes(avcw);

    for (iter=0; iter<iters; iter++)
    {
        if (!immobile)
        {
            // Entire molecule rotation.
            rotate(&xv, amx);
            bind1 = get_intermol_binding(ligands);
            if (avcw) bind -= get_intermol_clashes(avcw);

            if (bind1 < bind)
            {
                rotate(&xv, -amx);
                amx *= -0.75;
            }
            else
            {
                bind = bind1;
            }

            rotate(&yv, amy);
            bind1 = get_intermol_binding(ligands);
            if (avcw) bind -= get_intermol_clashes(avcw);

            if (bind1 < bind)
            {
                rotate(&yv, -amy);
                amy *= -0.75;
            }
            else
            {
                bind = bind1;
            }

            rotate(&zv, amz);
            bind1 = get_intermol_binding(ligands);
            if (avcw) bind -= get_intermol_clashes(avcw);

            if (bind1 < bind)
            {
                rotate(&zv, -amz);
                amz *= -0.75;
            }
            else
            {
                bind = bind1;
            }
        }

        intermol_conform_flexonly(ligands, iters, avcw, bind);
    }
}

void Molecule::intermol_conform_flexonly(Molecule* ligand, int iters, Molecule** avcw)
{
    intermol_conform_flexonly(ligand, iters, avcw, nanf(""));
}

void Molecule::intermol_conform_flexonly(Molecule** ligands, int iters, Molecule** avcw)
{
    intermol_conform_flexonly(ligands, iters, avcw, nanf(""));
}

void Molecule::intermol_conform_flexonly(Molecule* ligand, int iters, Molecule** avcw, float lastbind)
{
    Molecule* ligands[4];
    ligands[0] = ligand;
    ligands[1] = nullptr;
    intermol_conform_flexonly(ligands, iters, avcw, lastbind);
}

void Molecule::intermol_conform_flexonly(Molecule** ligands, int iters, Molecule** avcw, float lastbind)
{
    int iter, i, j, k, l;
    Point xpole(1,0,0),
          ypole(0,1,0),
          zpole(0,0,1);
    SCoord xv(&xpole),
           yv(&ypole),
           zv(&zpole);
    // cout << "Conforming " << name << endl << flush;

    float bind;
    bind = isnan(lastbind) ? get_intermol_binding(ligands) : lastbind;
    float bind1;
    if (avcw) bind -= get_intermol_clashes(avcw);
    bind -= get_internal_clashes();

    for (iter=0; iter<iters; iter++)
    {
        // Individual bond rotations.
        get_rotatable_bonds();

        if (rotatable_bonds)
        {
            for (j=0; rotatable_bonds[j]; j++)
            {
                rotatable_bonds[j]->rotate(rotatable_bonds[j]->angular_momentum);

                float bener = 0;
                Atom** mwb = rotatable_bonds[j]->get_moves_with_btom();
                if (mwb)
                {
                    for (k=0; mwb[k]; k++)
                    {
                        if (mwb[k]->last_bind_energy > bener) bener = mwb[k]->last_bind_energy;
                    }
                    delete mwb;
                }

                bind1 = get_intermol_binding(ligands);
                if (avcw) bind -= get_intermol_clashes(avcw);
                bind1 -= get_internal_clashes();

                if (bind1 < bind)
                {
                    float reversal = -(0.5 + 0.4 / (bener/10+1));
                    rotatable_bonds[j]->rotate(-rotatable_bonds[j]->angular_momentum);
                    rotatable_bonds[j]->angular_momentum *= reversal;
                    //cout << "Reversal " << rotatable_bonds[j]->angular_momentum << "|" << reversal << endl;
                }
                else
                {
                    bind = bind1;
                }
            }
        }
    }
}



void Molecule::multimol_conform(Molecule** mm, int iters, void (*cb)(int))
{
    if (!mm) return;

    int i, j, k, inplen, iter;
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

            float reversal = -0.666; // TODO: Make this binding-energy dedpendent.
            float accel = 1.1;

            /**** Linear Motion ****/
            if (mm[i]->movability >= MOV_ALL)
            {
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
                    pt.x = -pt.x;
                    mm[i]->move(pt);
                    mm[i]->lmx *= reversal;
                }
                else
                {
                    improvement += (bind1 - bind);
                    if (fabs(mm[i]->lmx) < 0.25) mm[i]->lmx *= accel;
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
                    if (fabs(mm[i]->lmy) < 0.25) mm[i]->lmy *= accel;
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
                    if (fabs(mm[i]->lmz) < 0.25) mm[i]->lmz *= accel;
                    bind = bind1;
                }

                mm[i]->lastbind = bind;
            }
            /**** End Linear Motion ****/

            /**** Axial Tumble ****/
            if (mm[i]->movability >= MOV_NORECEN)
            {
                Point pt(1,0,0);
                SCoord v(pt);

                rad = 0;
                bestfrb = 0;
                bestfrrad = nanf("No good results.");

                if (!(iter % _fullrot_every))
                {
                    while ((M_PI*2-rad) > 1e-3)
                    {
                        mm[i]->rotate(&v, _fullrot_steprad);
                        rad += _fullrot_steprad;

                        bind1 = 0;
                        for (j=0; mm[j]; j++)
                        {
                            if (!nearby[j]) continue;
                            bind1 += mm[i]->get_intermol_binding(mm[j]);
                        }

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
                    mm[i]->rotate(&v, mm[i]->amx);
                    bind1 = 0;
                    for (j=0; mm[j]; j++)
                    {
                        if (!nearby[j]) continue;
                        bind1 += mm[i]->get_intermol_binding(mm[j]);
                    }
                    if (bind1 < bind)
                    {
                        mm[i]->rotate(&v, -mm[i]->amx);
                        mm[i]->amx *= reversal;
                    }
                    else
                    {
                        improvement += (bind1 - bind);
                        bind = bind1;
                    }
                }

                pt.x=0;
                pt.y=1;
                SCoord v1(pt);

                rad = 0;
                bestfrb = 0;
                bestfrrad = nanf("No good results.");

                if (!(iter % _fullrot_every))
                {
                    while ((M_PI*2-rad) > 1e-3)
                    {
                        mm[i]->rotate(&v, _fullrot_steprad);
                        rad += _fullrot_steprad;

                        bind1 = 0;
                        for (j=0; mm[j]; j++)
                        {
                            if (!nearby[j]) continue;
                            bind1 += mm[i]->get_intermol_binding(mm[j]);
                        }

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
                    mm[i]->rotate(&v, mm[i]->amy);
                    bind1 = 0;
                    for (j=0; mm[j]; j++)
                    {
                        if (!nearby[j]) continue;
                        bind1 += mm[i]->get_intermol_binding(mm[j]);
                    }
                    if (bind1 < bind)
                    {
                        mm[i]->rotate(&v, -mm[i]->amy);
                        mm[i]->amy *= reversal;
                    }
                    else
                    {
                        improvement += (bind1 - bind);
                        bind = bind1;
                    }
                }


                pt.y=0;
                pt.z=1;
                SCoord v2(pt);

                rad = 0;
                bestfrb = 0;
                bestfrrad = nanf("No good results.");

                if (!(iter % _fullrot_every))
                {
                    while ((M_PI*2-rad) > 1e-3)
                    {
                        mm[i]->rotate(&v, _fullrot_steprad);
                        rad += _fullrot_steprad;

                        bind1 = 0;
                        for (j=0; mm[j]; j++)
                        {
                            if (!nearby[j]) continue;
                            bind1 += mm[i]->get_intermol_binding(mm[j]);
                        }

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
                    mm[i]->rotate(&v, mm[i]->amz);
                    bind1 = 0;
                    for (j=0; mm[j]; j++)
                    {
                        if (!nearby[j]) continue;
                        bind1 += mm[i]->get_intermol_binding(mm[j]);
                    }
                    if (bind1 < bind)
                    {
                        mm[i]->rotate(&v, -mm[i]->amz);
                        mm[i]->amz *= reversal;
                        //cout << "x";
                    }
                    else
                    {
                        improvement += (bind1 - bind);
                        bind = bind1;
                    }
                }


                mm[i]->lastbind = bind;
            }
            /**** End Axial Tumble ****/

            if ((iter % _fullrot_every)) continue;


            /**** Bond Flexion ****/
            if (mm[i]->movability >= MOV_FLEXONLY)
            {
                mm[i]->get_rotatable_bonds();

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

                // cout << mm[i]->name;
                if (mm[i]->rotatable_bonds)
                {
                    // cout << " has rotbonds ";
                    for (k=0; mm[i]->rotatable_bonds[k]; k++)
                    {
                        // cout << k;

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
                                for (j=0; mm[j]; j++)
                                {
                                    if (!nearby[j]) continue;
                                    bind1 += mm[i]->get_intermol_binding(mm[j]);
                                }

                                if (bind1 > bestfrb)
                                {
                                    bestfrb = bind1;
                                    bestfrrad = rad;
                                }
                            }

                            if (!isnan(bestfrrad))
                                mm[i]->rotatable_bonds[k]->rotate(bestfrrad);
                        }
                        else
                        {
                            float ra = mm[i]->rotatable_bonds[k]->angular_momentum;
                            mm[i]->rotatable_bonds[k]->rotate(ra);

                            bind1 = 0;
                            for (j=0; mm[j]; j++)
                            {
                                if (!nearby[j]) continue;
                                bind1 += mm[i]->get_intermol_binding(mm[j]);
                            }
                            if (bind1 < bind)
                            {
                                mm[i]->rotatable_bonds[k]->rotate(-ra);
                                mm[i]->rotatable_bonds[k]->angular_momentum *= reversal;
                            }
                            else
                            {
                                improvement += (bind1 - bind);
                                bind = bind1;
                            }
                        }
                    }
                }
                // cout << endl;
                mm[i]->lastbind = bind;

                if (!(iter % _fullrot_every)) mm[i]->reset_conformer_momenta();
            }
            /**** End Bond Flexion ****/
        }	// for i
        // cout << "Iteration " << iter << " improvement " << improvement << endl;

        if (cb) cb(iter);
    }	// for iter.
}





bool Molecule::from_smiles(char const * smilesstr)
{
    smlen = strlen(smilesstr);
    paren = new SMILES_Parenthetical[smlen];
    spnum = 0;

    bool retval = from_smiles(smilesstr, nullptr);

    int i;
    for (i=0; i<spnum; i++)
    {
        retval &= from_smiles(paren[i].smilesstr, paren[i].startsfrom);
    }

    delete[] paren;
    correct_structure(50);
    hydrogenate();

    // voxel_computation(5);

    return retval;
}

bool Molecule::from_smiles(char const * smilesstr, Atom* ipreva)
{
    Atom* stack[256];
    Atom* numbered[10];
    Atom* sequence[65536];
    bool seqarom[65536];
    int sqidx[10];
    int sp = 0;
    bool bracket=false, prevarom=false;
    Atom* bracketed=0;

    immobile = false;

    int i, j=1, k=0, l, atno=get_atom_count()+1;
    for (i=0; i<10; i++) numbered[i] = 0;

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
            char* aname = new char[5];
            i++;
            j=0;
            while (smilesstr[i] != '}' && j<5)
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
            j = smilesstr[i] - 48;
            if (!numbered[j])
            {
                numbered[j] = preva;
                sqidx[j] = k-1;
                continue;
            }
            else
            {
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
                
                if (!ring_atoms)
                	ring_atoms = new Atom**[16];
                
                if (!ring_aromatic)
                	ring_aromatic = new bool[16];
                
                if (!ring_atoms[ringcount])
                	ring_atoms[ringcount] = new Atom*[ringsz+2];
                
                for (l=0; l<ringsz; l++)
                {
                	ring_atoms[ringcount][l] = aloop[l];
            	}
                ring_atoms[ringcount][ringsz] = 0;
                
                ring_aromatic[ringcount] = ring_is_aromatic(ringcount);
                
                ringcount++;

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
                if (ringsz<5 || allarom)
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

                    numbered[j] = 0;
                    card = 1;

                    continue;
                }

                float anomaly = close_loop(aloop, card);
                // Ring closure anomaly of more than a tiny amount should always come up as a fail in the tests.
                if (fabs(anomaly) > 0.1) cout << "Ring closure anomaly " << anomaly << endl;

                /*for (l=0; l<dbi; l++)
                {	if (!EZgiven[l] && EZatom0[l] && EZatom1[l])
                	{	Bond* lb = EZatom0[l]->get_bond_between(EZatom1[l]);
                		lb->can_flip = false;
                		lb->can_rotate = false;
                		lb = EZatom1[l]->get_bond_between(EZatom0[l]);
                		if (lb)
                		{	lb->can_flip = false;
                			lb->can_rotate = false;
                		}
                	}
                }*/

                if (card) preva->bond_to(numbered[j], card);
                card = 1;

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
                    bracketed->swap_chirality = !bracketed->swap_chirality;
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

void Molecule::make_coplanar_ring(Atom** ring_members)
{
    if (!ring_members) return;
    int i, j, l, ringsz;

    bool allarom = true;
    for (i=0; ring_members[i]; i++)
    {
        ringsz = i+1;
        Bond** ab = ring_members[i]->get_bonds();
        bool haspi = false;
        for (j=0; ab[j]; j++)
            if (ab[j]->cardinality > 1 && ab[j]->cardinality < 2) haspi = true;
        if (!haspi) allarom = false;
    }

    if (ringsz<3) return;

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
}

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

Atom* Molecule::get_most_bindable()
{
    if (!atoms) return 0;

    int i, j=-1;
    float best=0;
    for (i=0; atoms[i]; i++)
    {
        float score = 0;

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

        if (score > best)
        {
            best = score;
            j = i;
        }
    }

    if (j < 0) return 0;
    else return atoms[j];
}

Point Molecule::get_bounding_box() const
{
    if (!atoms) return 0;

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

void Molecule::voxel_computation(int iters)
{
    if (!atoms) return;
    return;			// THIS FUNCTION DOES NOT WORK.

    // Get the molecule's bounding box, then grow it slightly.
    Point boxvtx = get_bounding_box();
    Point molcen = get_barycenter();

    boxvtx.scale(boxvtx.magnitude() * 1.25);

    // Define a voxel space with 0.1A resolution.
    int xmax = ceil(boxvtx.x * 2 / _voxel_resolution);
    int ymax = ceil(boxvtx.y * 2 / _voxel_resolution);
    int zmax = ceil(boxvtx.z * 2 / _voxel_resolution);
    int numvox = (xmax+1) * (ymax+1) * (zmax+1);
    int voxperx = (ymax+1) * (zmax+1);
    int voxpery = (zmax+1);
#if USE_VOXEL_ARRAY
    float voxelspace[numvox + 8];
#endif

    // For each iteration,
    int iter;
    for (iter=0; iter<iters; iter++)
    {
        int i, j, vx, vy, vz, vxi, vxyi;
        float ax, ay, az;
        float min_energy = 99999;
        Point min_ener_vox;
        Point aloc;

        // For each atom of the molecule,
        for (i=0; atoms[i]; i++)
        {
            cout << iter << "," << i << " / " << iters << endl;

            // Find its energy level for every voxel in the space.
            // For non-bonded atoms, just use an inverse square of the distance.
            // For bonded atoms, use the negative inverse square of 1+distance anomaly.
            // Keep track of the lowest energy level.
            float energy;
            min_energy = 99999;
            min_ener_vox = aloc = atoms[i]->get_location();

            for (vx=0; vx<xmax; vx++)
            {
                vxi = vx * voxperx;
                if (vxi >= numvox) break;
                ax = _voxel_resolution * vx;

                for (vy=0; vy<ymax; vy++)
                {
                    vxyi = vxi + vy * voxpery;
                    if (vxyi >= numvox) break;
                    ay = _voxel_resolution * vy;

                    for (vz=0; vz<zmax; vz++)
                    {
                        az = _voxel_resolution * vz;
                        energy = 0;
                        Point aloctmp(ax, ay, az);
                        aloctmp = aloctmp.add(molcen);

                        for (j=0; atoms[j]; j++)
                        {
                            if (j==i) continue;
                            Bond* b = atoms[i]->get_bond_between(atoms[j]);
                            float r = aloctmp.get_3d_distance(atoms[j]->get_location());
                            if (b)
                            {
                                float optimal_r = InteratomicForce::covalent_bond_radius(atoms[i], atoms[j], b->cardinality);
                                float anomaly = fabs(r - optimal_r);
                                energy -= 1.0/pow(1.0+anomaly, 2);
                            }
                            else
                            {
                                energy += 1.0/pow(0.00001+r, 2);
                            }
                        }

#if USE_VOXEL_ARRAY
                        int vxyzi = vxyi + vz;
                        if (vxyzi >= numvox) break;

                        voxelspace[vxyzi] = energy;
#endif

                        if (energy < min_energy)
                        {
                            min_energy = energy;
                            min_ener_vox = aloctmp;
                        }
                    }	// next vz
                }	// next vy
            }	// next vx

            // Move the atom to the lowest energy voxel.
            atoms[i]->move(min_ener_vox);
        }
    }
}

#define _DEV_FIX_MSTRUCT 0
void Molecule::correct_structure(int iters)
{
    if (!atoms) return;
    int iter, i, j, k;
    Point zero(0,0,0);
    
    if (ringcount)
    {
    	for (i=0; i<ringcount; i++)
    	{
    		if (ring_aromatic[i])
    		{
				// cout << "Ring " << i << endl;
				make_coplanar_ring(ring_atoms[i]);
    		}
    	}
    }

    for (iter=0; iter<iters; iter++)
    {
        for (i=0; atoms[i]; i++)
        {
            Point aloc = atoms[i]->get_location();
            Bond** b = atoms[i]->get_bonds();
            int g = atoms[i]->get_geometry();

            for (j=0; j<g; j++)
            {
                if (b[j]->btom)
                {
                    Point bloc = b[j]->btom->get_location();
                    SCoord v(bloc.subtract(aloc));
                    float optimal = InteratomicForce::covalent_bond_radius(atoms[i], b[j]->btom, b[j]->cardinality);
                    
                    if (g == 4 && b[j]->cardinality > 1 && b[j]->cardinality <= 2) atoms[i]->aromatize();
                    
                    // Failed attempt at bond angles.
                    if (iter < (iters-20))
                    {
		                for (k=0; k<g; k++)
		                {
		                	if (k == j) continue;
		                	if (b[k]->btom)
		                	{
		                		Point cloc = b[k]->btom->get_location();
		                		float f=0, theta = find_3d_angle(bloc, cloc, aloc);
		                		
		                		switch (g)
		                		{
		                			case 2:   f = M_PI          - theta;	break;		                			
		                			case 3:   f = triangular    - theta;	break;		                			
		                			case 4:   f = tetrahedral   - theta;	break;		                			
		                			case 6:   f = square        - theta;	break;		                			
		                			default:  ;
		                		}

	                			#if _DEV_FIX_MSTRUCT
		                		if (fabs(f) > 0.1)
		                		{
				            		cout << atoms[i]->name << "-" << b[j]->btom->name << "(" << b[k]->btom->name << ")"
		                				 << " " << g << " " << theta*fiftyseven << " " << f*fiftyseven << endl;
		                			
				            		Point pt(v);
				            		
				            		SCoord normal = compute_normal(aloc, bloc, cloc);
				            		Point pt0 = rotate3D(pt, zero, normal, -0.1*f);
				            		Point pt1 = rotate3D(pt, zero, normal,  0.1*f);
				            		
				            		float r0 = pt0.get_3d_distance(cloc);
				            		float r1 = pt1.get_3d_distance(cloc);
				            		
				            		v = (r0 > r1) ? pt0 : pt1;
			            		}
			            		#endif
		                	}
		                }
	                }

                    v.r = optimal; // += 0.1 * (optimal-v.r);

                    if (atoms[i]->num_rings())
                    {
                    	// b[j]->btom->move(aloc.add(v));
                    	;
                	}
                    else
                    {
		                Point pt(aloc.add(v));
		                b[j]->btom->move_assembly(&pt, atoms[i]);
	                }
                }
            }
        }
    }

}












