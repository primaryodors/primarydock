
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

Atom* Molecule::nearest_local_atom = nullptr;
Atom* Molecule::nearest_remote_atom = nullptr;
int Molecule::nearest_local_aidx = 0;
float Molecule::nearest_r = 0;

float conformer_momenta_multiplier = 1;
float conformer_tumble_multiplier = 1;

float cavity_stuffing = default_cavity_stuffing;
float clash_fleeing = lmpush;

bool allow_ligand_360_tumble = true;
bool allow_ligand_360_flex = true;

Molecule *worst_clash_1 = nullptr, *worst_clash_2 = nullptr;
float worst_mol_clash = 0;

#if _dbg_improvements_only_rule
Molecule** check_mols = nullptr;
Molecule* check_ligand = nullptr;
bool excuse_deterioration = false;
#endif

float _momentum_rad_ceiling = fiftyseventh * 30;

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

    int j;
    for (j=0; j<10; j++) lastbind_history[j] = 0;
}

Molecule::Molecule()
{
    atoms = nullptr;
    smiles = nullptr;
    rings = nullptr;
    atcount = 0;
    reset_conformer_momenta();
    rotatable_bonds = 0;
    paren = nullptr;

    int j;
    for (j=0; j<10; j++) lastbind_history[j] = 0;
}

Molecule::~Molecule()
{
    if (atoms)
    {
        delete[] atoms;
    }
    if (smiles) delete[] smiles;
    if (rings) delete[] rings;
    if (rotatable_bonds) delete[] rotatable_bonds;
    if (vdw_surface) delete[] vdw_surface;
    if (vdw_vertex_atom) delete[] vdw_vertex_atom;
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
    identify_conjugations();
    rotatable_bonds = 0;
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

Pose::~Pose()
{
    //
}

void Pose::reset()
{
    sz = 0;
    saved_atom_locs.clear();
    saved_atom_Z.clear();
    saved_from = nullptr;
}

void Pose::copy_state(Molecule* m)
{
    int i;
    if (!saved_atom_locs.size() || saved_from != m)
    {
        saved_atom_locs.clear();
        saved_atom_Z.clear();
        saved_from = m;
        if (!m || !m->atoms) return;

        sz = m->get_atom_count();
        for (i=0; i<=sz; i++)
        {
            Point pt;
            saved_atom_locs.push_back(pt);
            saved_atom_Z.push_back(0);
        }
    }

    for (i=0; m->atoms[i] && i<sz; i++)
    {
        saved_atom_locs[i] = m->atoms[i]->get_location();
        saved_atom_Z[i] = m->atoms[i]->get_Z();
    }
    sz = i;
}

void Pose::restore_state(Molecule* m)
{
    if (!m || !m->atoms || !sz) return;
    int i, n;
    if (m != saved_from)
    {
        n = saved_atom_locs.size();
        for (i=0; i<n; i++)
        {
            if (i == n-1 && !m->atoms[i] && (saved_atom_Z[i] == 1)) break;
            if (!m->atoms[i] && !saved_atom_Z[i]) break;
            if (/*n != sz ||*/ !m->atoms[i] || (saved_atom_Z[i] != m->atoms[i]->get_Z()))
            {
                cout << "Attempt to restore pose to incompatible molecule (from " << saved_from->name << " to " << m->name << ")." << endl;
                if (m->is_residue())
                {
                    Star s;
                    s.pmol = m;
                    if (s.paa->conditionally_basic()) return;
                }
                throw -4;
            }
        }

        saved_from = m;
    }

    for (i=0; i<sz && m->atoms[i]; i++)
    {
        m->atoms[i]->move(saved_atom_locs[i]);
    }
}

float Pose::total_atom_motions()
{
    if (!saved_from || !saved_from->atoms || !sz) return 0;
    int i;
    float result = 0;

    for (i=0; i<sz && saved_from->atoms[i]; i++)
    {
        float r = saved_from->atoms[i]->get_location().get_3d_distance(saved_atom_locs[i]);
        result += r;
    }

    return result;
}

void Molecule::delete_atom(Atom* a)
{
    if (!a) return;
    paths = nullptr;

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
                for (atcount=0; atoms[atcount]; atcount++);     // Get count.
                return;
            }
        }
    }

    cout << "Attempt to delete atom " << a->name << " from a molecule it is not part of." << endl << flush;
    throw 0xbada70b;
}

void Molecule::delete_all_atoms()
{
    if (!atoms) return;
    if (paths) delete[] paths;
    paths = nullptr;

    int i;
    for (i=0; atoms[i]; i++)
    {
        delete atoms[i];
    }
    
    delete[] atoms;
    atoms = nullptr;
    atcount = 0;
}

void Molecule::reset_conformer_momenta()
{
    srand(time(nullptr));

    lmx = _def_lin_momentum * sgn(0.5-(rand()&1));
    lmy = _def_lin_momentum * sgn(0.5-(rand()&1));
    lmz = _def_lin_momentum * sgn(0.5-(rand()&1));
    amx = _def_ang_momentum * conformer_momenta_multiplier * conformer_tumble_multiplier * sgn(0.5-(rand()&1));
    amy = _def_ang_momentum * conformer_momenta_multiplier * conformer_tumble_multiplier * sgn(0.5-(rand()&1));
    amz = _def_ang_momentum * conformer_momenta_multiplier * conformer_tumble_multiplier * sgn(0.5-(rand()&1));

    Bond** b = get_rotatable_bonds();
    int i;

    if (b)
    {
        for (i=0; b[i]; i++)
        {
            b[i]->angular_momentum = _def_bnd_momentum * conformer_momenta_multiplier * sgn(0.5-(rand()&1));
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

        if (atcount) cout << "Molecule " << (name?name:"(no name)") << " has " << atcount << " atoms." << endl;

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

void Molecule::add_existing_atom(Atom* a)
{
    if (!atoms) atcount = 0;
    else for (atcount=0; atoms[atcount]; atcount++);     // Get count.

    reallocate();
    atoms[atcount++] = a;
    atoms[atcount] = nullptr;

    if (atcount > 1)
    {
        strcpy(a->aa3let, atoms[0]->aa3let);
        a->residue = atoms[0]->residue;
        a->aaletter = atoms[0]->aaletter;
    }

    clear_all_bond_caches();
}

Atom* Molecule::add_atom(char const* elemsym, char const* aname, const Point* location, Atom* bond_to, const float bcard, const int charge)
{
    paths = nullptr;

    Atom* a = new Atom(elemsym, location, charge);
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

    paths = nullptr;

    if (bondto->is_pi()) bondto->aromatize();
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

        clear_all_bond_caches();

        if ((atcount & 1) && bondto->get_bonded_atoms_count() == 2)
        {
            Bond* b = bondto->get_bond_between(a);
            if (b && b->can_rotate)
            {
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
    if (!atoms) return;
    int i, j;
    for (i=0; atoms[i]; i++)
    {
        Bond* b[16];
        atoms[i]->fetch_bonds(b);
        if (b[0])
        {
            for (j=0; b[j]; j++) b[j]->clear_moves_with_cache();
        }
        atoms[i]->used = 0;
    }
}

int Molecule::is_residue()
{
    if (noAtoms(atoms)) return 0;
    int i;

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->residue) return atoms[i]->residue;
    }
    return 0;
}

int Molecule::get_hydrogen_count()
{
    if (noAtoms(atoms)) return 0;
    int i, retval=0;

    for (i=0; atoms[i]; i++)
        if (atoms[i]->get_Z() == 1)
            retval++;
    
    return retval;
}

int Molecule::count_atoms_by_element(const char* esym)
{
    if (noAtoms(atoms)) return 0;
    int findZ = Atom::Z_from_esym(esym);
    int i, retval=0;

    for (i=0; atoms[i]; i++)
        if (atoms[i]->get_Z() == findZ)
            retval++;
    
    return retval;
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

        #if _dbg_hydrogenate
        cout << atoms[i]->name << " has valence " << valence << endl;
        #endif

        float bcardsum = 0;

        Bond* aib[16];
        atoms[i]->fetch_bonds(aib);
        int db = 0;
        if (aib[0])
        {
            for (j=0; aib[j]; j++)
            {
                if (aib[j]->atom2) bcardsum += aib[j]->cardinality;
                if (aib[j]->cardinality > 1) db++;
            }
        }
        if (!db && steric_only) continue;

        if (bcardsum && atoms[i]->get_Z() == 1) continue;
        bcardsum = ceil(bcardsum);

        #if _dbg_hydrogenate
        cout << " minus existing bonds " << bcardsum;
        #endif

        bcardsum -= atoms[i]->get_charge();
        if (atoms[i]->is_backbone && !strcmp(atoms[i]->name, "N")) bcardsum = 1;

        #if _dbg_hydrogenate
        cout << " given charge makes " << bcardsum << endl;
        #endif

        int fam = atoms[i]->get_family();
        if (fam == PNICTOGEN || fam == CHALCOGEN)
        {
            Atom* C = atoms[i]->is_bonded_to_pi(TETREL, true);
            if (!C) C = atoms[i]->is_bonded_to_pi(PNICTOGEN, true);
            if (C) atoms[i]->aromatize();
        }

        int h_to_add = round(valence - bcardsum);
        if (atoms[i]->aaletter == 'H' && !strcmp(atoms[i]->name, "NE2")) h_to_add++;
        for (j=0; j<h_to_add; j++)
        {
            char hname[15];
            sprintf(hname, "H%d", atcount+1);
            Atom* H = add_atom("H", hname, atoms[i], 1);
            #if _dbg_hydrogenate
            cout << "Adding " << hname << " to " << atoms[i]->name << " whose valence is " << valence << " and has " << bcardsum << " bonds already." << endl;
            #endif

            /*atoms[i]->clear_geometry_cache();
            SCoord v = atoms[i]->get_next_free_geometry(1);
            v.r = InteratomicForce::covalent_bond_radius(atoms[i], H, 1);
            H->move(atoms[i]->get_location().add(v));*/

            if (atoms[i]->get_geometry() == 3)
            {
                Bond* aib = atoms[i]->get_bond_by_idx(1);
                if (1) // aib && aib->total_rotations) // aib->atom2 == H)
                {
                    /*Bond* b0 = atoms[i]->get_bond_by_idx(0);
                    Bond* b1 = atoms[i]->get_bond_by_idx(1);
                    if (!b1 || !b1->atom2) b1 = atoms[i]->get_bond_by_idx(2);*/

                    int k=0;
                    Bond* b0 = atoms[i]->get_bond_by_idx(k++);
                    if (!b0->atom2 || b0->atom2 == H) b0 = atoms[i]->get_bond_by_idx(k++);
                    Bond* b1 = atoms[i]->get_bond_by_idx(k++);
                    if (!b1->atom2 || b1->atom2 == H) b1 = atoms[i]->get_bond_by_idx(k++);

                    if (b0->atom2 && b1->atom2 && (abs((__int64_t)(b1) - (__int64_t)b1->atom2) < memsanity))
                    {
                        Point source = atoms[i]->get_location();
                        /*Point axis = b0->atom2->get_location();
                        Point avoid = b1->atom2->get_location();
                        Point movable = H->get_location();
                        SCoord v(axis.subtract(source));

                        Point visibly_moved = rotate3D(movable, source, v, M_PI);
                        float r_is = movable.get_3d_distance(avoid);
                        float r_wouldbe = visibly_moved.get_3d_distance(avoid);

                        if (r_is < r_wouldbe) H->move(visibly_moved);*/

                        Point movable = b1->atom2->get_location();
                        /*cout << "Getting normal of " << atoms[i]->name << " - "
                        	 << b0->atom2->name << " - " << b1->atom2->name << endl;*/
                        SCoord axis = compute_normal(source, b0->atom2->get_location(), b1->atom2->get_location());
                        Point plus  = rotate3D(movable, source, axis,  triangular);
                        Point minus = rotate3D(movable, source, axis, -triangular);

                        float rp = plus.get_3d_distance(b0->atom2->get_location());
                        float rm = minus.get_3d_distance(b0->atom2->get_location());

                        Point pt = ((rp > rm) ? plus : minus).subtract(source);
                        pt.scale(InteratomicForce::covalent_bond_radius(atoms[i], H, 1));
                        H->move(pt.add(source));
                    }
                }
            }
        }
    }

    for (atcount=0; atoms[atcount]; atcount++);     // Update count.

    clear_all_bond_caches();
}

void Molecule::dehydrogenate()
{
    if (!atoms) return;

    Atom* tmp[atcount];
    int i, j=0;
    for (i=0; i<atcount; i++)
    {
        if (!atoms[i]) continue;
        if (atoms[i]->get_Z() < 2) continue;
        tmp[j++] = atoms[i];
    }

    for (i=0; i<j; i++) atoms[i] = tmp[i];
    atoms[j] = nullptr;
}

char** Molecule::get_atom_names() const
{
    int i;
    char** retval = new char*[atcount+1];

    for (i=0; atoms[i]; i++) retval[i] = atoms[i]->name;
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
{
    if (noAtoms(atoms))
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

std::vector<Atom*> Molecule::longest_dimension()
{
    std::vector<Atom*> retval;
    if (!atoms) return retval;
    int i, j;
    float rmax = 0;

    for (i=0; atoms[i]; i++)
    {
        for (j=i+1; atoms[j]; j++)
        {
            float r = atoms[i]->distance_to(atoms[j]);
            if (r > rmax)
            {
                retval.clear();
                retval.push_back(atoms[i]);
                retval.push_back(atoms[j]);
                rmax = r;
            }
        }
    }

    return retval;
}

float Molecule::total_eclipses()
{
    #if !include_eclipses
    return 0;
    #else
    if (!atoms) return 0;
    float result = 0;
    int i, j, k, l, m, n;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->get_Z() < 2) continue;
        Bond* abt[16];
        atoms[i]->fetch_bonds(abt);
        n = atoms[i]->get_geometry();
        for (j=0; j<n; j++)
        {
            if (!abt[j]) continue;
            if (abt[j]->cardinality > 1) continue;
            if (!abt[j]->can_rotate) continue;
            if (!abt[j]->atom2) continue;
            if (abt[j]->atom2->get_Z() < 2) continue;
            if (atoms[i]->is_pi() && abt[j]->atom2->is_pi()) continue;
            if (atoms[i]->is_backbone && abt[j]->atom2->is_backbone) continue;
            SCoord axis = abt[j]->atom2->get_location().subtract(atoms[i]->get_location());
            Bond* bbt[16];
            abt[j]->atom2->fetch_bonds(bbt);
            m = abt[j]->atom2->get_geometry();
            for (k=0; k<m; k++)
            {
                if (!bbt[k]) continue;
                if (!bbt[k]->atom2) continue;
                if (bbt[k]->atom2->get_Z() > 1) continue;
                // if (bbt[k]->atom2 == atoms[i]) continue;
                for (l=0; l<n; l++)
                {
                    if (l == j) continue;
                    if (!abt[l]) continue;
                    if (!abt[l]->atom2) continue;
                    if (abt[l]->atom2->get_Z() > 1) continue;
                    float theta = find_angle_along_vector(bbt[k]->atom2->get_location(), abt[l]->atom2->get_location(), atoms[i]->get_location(), axis);
                    #if _dbg_eclipses
                    cout << bbt[k]->atom2->name << " is " << (theta*fiftyseven) << "deg from " << abt[l]->atom2->name
                        << " along the " << atoms[i]->name << " - " << abt[j]->atom2->name << " axis."
                        << endl;
                    #endif
                    theta -= M_PI;
                    while (theta < -hexagonal) theta += hexagonal;
                    while (theta > hexagonal) theta -= hexagonal;
                    result += fabs(theta);
                }
            }
        }
    }

    return result*eclipsing_kJmol_per_radian;
    #endif
}

float Molecule::bindability_by_type(intera_type t, bool ib)
{
    if (!atoms) return 0;

    float result = 0;
    int i, heavy = 0;
    float c;
    for (i=0; atoms[i]; i++)
    {
        if (!ib && atoms[i]->is_backbone) continue;
        if (atoms[i]->get_Z() > 1) heavy++;
        switch (t)
        {
            case covalent:
            // TODO
            ;
            break;

            case ionic:
            c = atoms[i]->get_charge();
            if (c) result += c;
            break;

            case hbond:
             if (atoms[i]->get_Z() > 1)
             {
                 result += fabs(atoms[i]->is_polar());
                 result += fabs(atoms[i]->get_charge());
             }
            break;

            case pi:
            result += fabs(atoms[i]->is_pi());
            break;

            case polarpi:
            result += fabs(atoms[i]->is_polar()) + fabs(atoms[i]->is_pi());
            break;

            case mcoord:
            if (atoms[i]->is_metal())
            {
                result += atoms[i]->get_charge();
            }
            else
            {
                int fam = atoms[i]->get_family();
                int Z = atoms[i]->get_Z();

                if (fam == PNICTOGEN && Z <= 15)
                    result -= (atoms[i]->is_pi() ? 0.5 : 1);

                if (fam == CHALCOGEN && Z <= 35)
                    result -= (atoms[i]->is_pi() ? 0.25 : 1);
            }

            break;

            case vdW:
            default:
            result += 1;
            break;
        }
    }

    if (heavy) result /= heavy;
    if (t == hbond)
    {
        if (result < 0.3) result = 0;
    }

    return result;
}

int Molecule::has_hbond_donors()
{
    if (!atoms) return 0;

    int i, result=0;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        if (atoms[i]->get_family() == TETREL) continue;
        if (atoms[i]->is_polar() >= hydrophilicity_cutoff)
        {
            result++;
            // if (get_charge() < 0.5) cout << name << ":" << atoms[i]->name << " is an hbond donor." << endl;
        }
        if (atoms[i]->get_family() == HALOGEN && atoms[i]->is_bonded_to(TETREL)) result++;
    }

    return result;
}

int Molecule::has_hbond_acceptors()
{
    if (!atoms) return 0;

    int i, result=0;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        if (atoms[i]->get_family() == TETREL) continue;
        if (atoms[i]->get_family() == HALOGEN && atoms[i]->is_bonded_to(TETREL)) continue;
        if (atoms[i]->get_bonded_atoms_count() > 3) continue;
        if (atoms[i]->is_polar() <= -hydrophilicity_cutoff) result++;
    }

    return result;
}

int Molecule::has_pi_atoms(bool ib)
{
    if (!atoms) return 0;

    int i, result=0;
    for (i=0; atoms[i]; i++)
    {
        if (!ib && atoms[i]->is_backbone) continue;
        if (atoms[i]->is_pi()) result++;
    }

    return result;
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
    char** words;

    for (j=3; j<lncount; j++)
    {
        char line[1024];
        strncpy(line, lines[j], 1023);
        words = chop_spaced_words(line);

        if (!words || !words[0] || !words[1]) break;
        if (!strcmp(words[1], "END")) break;

        if (!strcmp(words[1], "CHG"))
        {
            for (i=3; words[i] && words[i+1]; i+=2)
            {
                int aidx = atoi(words[i]);
                atoms[aidx-1]->increment_charge(atof(words[i+1]));
            }
        }

        if (j == 3)
        {
            na = atoi(words[0]);
            nb = atoi(words[1]);

            atoms = new Atom*[na+4];
        }
        else if (added < na)
        {
            Point* loc = new Point(atof(words[0]), atof(words[1]), atof(words[2]));
            if (words[3][0] >= 'a' && words[3][0] <= 'z') words[3][0] -= 0x20;
            Atom* a = new Atom(words[3], loc);
            delete loc;
            a->name = new char[16];
            sprintf(a->name, "%s%d", words[3], added+1);
            a->residue = 0;
            atoms[atcount++] = a;
            atoms[atcount] = nullptr;
            added++;
        }
        else
        {
            int a1i = atoi(words[0]);
            int a2i = atoi(words[1]);

            if (!a1i || !a2i) break;
            atoms[a1i-1]->bond_to(atoms[a2i-1], atof(words[2]));
        }

        if (words) delete[] words;
    }
    atoms[atcount] = 0;
    if (words) delete[] words;

    identify_conjugations();
    identify_rings();
    identify_cages();
    identify_acidbase();
    return added;
}

void Molecule::identify_conjugations()
{
    if (!atoms) return;
    int i;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_pi() && !atoms[i]->conjugation)
        {
            Conjugation* conj = new Conjugation(atoms[i]);          // Don't have to save this pointer because the atoms will save it.
        }
    }
}

bool Molecule::check_Greek_continuity()
{
    if (!atoms) return true;
    int i, j, k, l, m;
    for (i=0; atoms[i]; i++)
    {
        if (!atoms[i]->check_Greek_continuity()) return false;

        continue;

        if (!atoms[i]->residue) continue;
        if (atoms[i]->is_backbone) continue;
        if (atoms[i]->get_Z() < 2) continue;

        int n = atoms[i]->get_geometry();
        if (!n) continue;
        Bond* bb[n+8];
        for (j=0; j<n+8; j++) bb[j] = nullptr;
        atoms[i]->fetch_bonds(bb);

        int ig = greek_from_aname(atoms[i]->name);
        if (ig < 0) continue;
        for (j=i+1; atoms[j]; j++)
        {
            if (atoms[j]->get_Z() < 5) continue;
            int jg = greek_from_aname(atoms[j]->name);
            if (jg > ig)
            {
                bool found = false;
                for (k=0; bb[k]; k++)
                {
                    Atom* mwb[256];
                    for (l=0; l<256; l++) mwb[l] = nullptr;
                    bb[k]->fetch_moves_with_atom2(mwb);
                    for (l=0; mwb[l]; l++)
                    {
                        if (mwb[l] == atoms[j])
                        {
                            found = true;
                            goto _exit_mwbsearch;
                        }
                    }
                }
                if (!found) return false;
            }
            _exit_mwbsearch:
            ;
        }
    }

    return true;
}

int Molecule::from_pdb(FILE* is, bool het_only)
{
    /*
    ATOM     55  SG  CYS     4       6.721  -8.103   4.542  1.00001.00           S
    */
    char buffer[1024];
    int added=0;

    while (!feof(is))
    {
        fgets(buffer, 1003, is);
        int charge = 0, offset = (buffer[21] != ' ' && buffer[22] == ' ') ? 1 : 0;
        char** words = chop_spaced_words(buffer);

        if (buffer[78] && buffer[78] > ' ')
        {
            char chgstr[3];
            chgstr[0] = buffer[78];
            chgstr[1] = buffer[79];
            chgstr[2] = 0;

            if 		(!strcmp(chgstr, "+")) charge = 1;
            else if	(!strcmp(chgstr, "++")) charge = 2;
            else if	(!strcmp(chgstr, "-")) charge = -1;
            else if	(!strcmp(chgstr, "--")) charge = -2;
            else if (atoi(chgstr)) charge = atoi(chgstr);
        }

        if (words)
        {
            if (
                  (!strcmp(words[0], "ATOM") && !het_only)
                  ||
                  !strcmp(words[0], "HETATM")
               )
            {
                try
                {
                    char esym[7];
                    if (words[2][0] >= '0' && words[2][0] <= '9')
                        strcpy(esym, &words[2][1]);
                    else
                        strcpy(esym, words[2]);

                    int i;
                    for (i=1; i<6; i++)
                    {
                        if (!esym[i+1]) esym[i] = 0;
                        if (esym[i+1] >= '0' && esym[i+1] <= '9') esym[i]=0;
                        if (i>1) esym[i] = 0;
                        if (!esym[i]) break;
                    }
                    esym[1] &= 0x5f;

                    // cout << buffer[21] << buffer[22] << " " << offset << endl;
                    Point aloc(atof(words[5+offset]), atof(words[6+offset]),atof(words[7+offset]));

                    Atom* a = add_atom(esym, words[2], &aloc, 0, 0, charge);
                    added++;

                    // a->residue = atoi(words[4]);

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

        delete words;
    }

    identify_conjugations();
    return added;
}

int Molecule::get_bond_count(bool unidirectional) const
{
    int i, j, bc=0;

    if (noAtoms(atoms)) return 0;
    for (i=0; i<atcount && atoms[i]; i++)
    {
        Bond* b[16];
        atoms[i]->fetch_bonds(b);
        if (!b[0]) continue;

        for (j=0; b[j]; j++)
        {
            if (b[j]->atom2 > atoms[i] || !unidirectional) bc++;
        }
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
            if (latoms[j] == lbonds[i]->atom1) laidx = j+1;
            if (latoms[j] == lbonds[i]->atom2) lbidx = j+1;
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

void Molecule::save_pdb(FILE* os, int atomno_offset, bool endpdb)
{
    int i;

    for (i=0; atoms[i]; i++)
    {
        atoms[i]->save_pdb_line(os, i+1+atomno_offset);
    }

    if (endpdb) fprintf(os, "\nTER\nEND\n\n");
    else fprintf(os, "\n\n");
}

int Molecule::add_ring(Atom** atoms)
{
    int i, j, l, n, ringcount;

    l = 0;
    Atom* min_atom = nullptr;
    for (i=0; atoms[i]; i++)
    {
        if (!min_atom || atoms[i] < min_atom)
        {
            min_atom = atoms[i];
            l = i;
        }
    }
    n = i;

    i = l-1; if (i<0) i += n;
    j = l+1; if (j>=n) j -= n;
    bool reversed = (atoms[j] > atoms[i]);

    Atom* atoms_ordered[n+4];
    for (i=0; i<n; i++)
    {
        if (reversed) j = l - i;
        else j = i + l;
        if (j < 0) j += n;
        if (j >= n) j -= n;
        atoms_ordered[i] = atoms[j];
    }
    atoms_ordered[n] = nullptr;


    Ring* r = new Ring(atoms_ordered);
    int m = r->get_atom_count();

    if (rings)
    {
        bool already_exists = false;
        for (i=0; rings[i]; i++)
        {
            if (rings[i]->get_overlap_count(r) == m) already_exists = true; 
        }
        ringcount = i;
        if (already_exists)
        {
            // delete r;                // For some reason, uncommenting this line causes segfaults.
            return ringcount;
        }
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

    ringstmp[ringcount++] = r;
    ringstmp[ringcount] = nullptr;
    rings = ringstmp;

    return ringcount-1;
}

int Molecule::identify_rings()
{
    find_paths();
    // return 0;

    Atom** ringstmp[256];
    int ringcount;
    Atom *a;
    int chainlen[256];
    bool is_ring[256];
    int found_rings=0, chains=0, cnvchain, active, i, j, k, l, m, n, p;
    Atom *cnva, *cnvb;
    Atom *ra, *rb;

    if (!rings) return 0;

    for (ringcount = 0; rings[ringcount]; ringcount++)
    {
        #if _dbg_identify_rings
        cout << "Ring number " << ringcount << *rings[ringcount] << ": coplanar? " << rings[ringcount]->is_coplanar() << ", conjugated? " << rings[ringcount]->is_conjugated() << endl;
        #endif

        if (rings[ringcount]->is_coplanar() && rings[ringcount]->is_conjugated())
        {
            rings[ringcount]->aromatize();
            #if _dbg_identify_rings
            cout << "Aromatized." << endl;
            #endif
        }
    }
    return ringcount;
}

int Molecule::path_contains_atom(int path_idx, Atom* a)
{
    if (!paths) return 0;
    if (!paths[path_idx]) return 0;

    int i;
    for (i=0; paths[path_idx][i]; i++)
        if (paths[path_idx][i] == a) return i+1;

    return 0;
}

int Molecule::path_get_length(int path_idx)
{
    if (!paths) return 0;
    if (!paths[path_idx]) return 0;

    int i;
    for (i=0; paths[path_idx][i]; i++);

    return i;
}

Atom* Molecule::path_get_terminal_atom(int path_idx)
{
    if (!paths) return nullptr;
    if (!paths[path_idx]) return nullptr;

    int n = path_get_length(path_idx);
    if (!n) return nullptr;
    return paths[path_idx][n-1];
}

void Molecule::copy_path(int old_idx, int new_idx)
{
    if (!paths) return;
    if (!paths[old_idx]) return;
    if (abs((__int64_t)paths[old_idx][0] - (__int64_t)paths[old_idx]) >= memsanity) return;

    if (!paths[new_idx]) paths[new_idx] = new Atom*[get_atom_count()];
    int i;
    for (i=0; paths[old_idx][i]; i++)
        paths[new_idx][i] = paths[old_idx][i];
    
    paths[new_idx][i] = nullptr;
}

bool Molecule::path_is_subset_of(int short_path, int long_path)
{
    if (!paths) return false;
    if (!paths[long_path]) return false;
    if (!paths[short_path]) return true;

    int i;
    for (i=0; paths[short_path][i] && paths[long_path][i]; i++)
    {
        if (paths[long_path][i] != paths[short_path][i]) return false;
    }

    return true;
}

void Molecule::echo_path(int i)
{
    int j;
    cout << i << ":";
    for (j=0; paths[i][j]; j++)
        cout << " " << paths[i][j]->name;
    
    cout << endl;
}

void Molecule::find_paths()
{
    if (!atoms || !atoms[0]) return;

    if (paths)
    {
        #if _dbg_path_search
        cout << "Paths already set; skipping." << endl;
        #endif
        return;
    }
    #if _dbg_path_search
    cout << "Searching for paths..." << endl;
    #endif

    int h, i, j, k, l, m, n = 0, p, q, limit;
    Bond* b[16];
    for (i=0; atoms[i]; i++)
    {
        n += atoms[i]->get_bonded_heavy_atoms_count();
    }

    // paths = new Atom**[n];
    limit = n*n;
    Atom** lpaths[limit];
    paths = lpaths;
    for (i=0; i<limit; i++) paths[i] = nullptr;

    atcount = get_atom_count();

    Atom* a = atoms[0];
    a->fetch_bonds(b);
    if (!b[0]) return;
    n=0;
    for (i=0; b[i]; i++)
    {
        if (!b[i]->atom2) continue;
        if (b[i]->atom2->get_Z() < 2) continue;
        if (b[i]->atom2->residue && b[i]->atom2->residue != a->residue) continue;
        paths[n] = new Atom*[atcount];
        for (q=0; q<atcount; q++) paths[n][q] = nullptr;
        paths[n][0] = a;
        paths[n][1] = b[i]->atom2;
        paths[n][2] = nullptr;
        n++;
    }

    int num_added, iter;
    for (iter=0; iter<1000; iter++)
    {
        num_added = 0;
        p = n;

        for (i=0; i<p; i++)
        {
            m = path_get_length(i);
            a = path_get_terminal_atom(i);
            if (!a) continue;
            if (!a->name) continue;
            if (abs((__int64_t)(atoms[0]) - (__int64_t)a) >= memsanity) continue;
            a->fetch_bonds(b);
            if (!b[0]) continue;

            k=0;
            for (j=0; b[j]; j++)
            {
                if (abs((__int64_t)(a) - (__int64_t)b[j]) > memsanity) break;
                if (abs((__int64_t)(b[j]) - (__int64_t)b[j]->atom2) > memsanity) break;
                if (!b[j]->atom2) continue;
                if (b[j]->atom2->get_Z() < 2) continue;
                if (b[j]->atom2->get_bonded_heavy_atoms_count() < 2) continue;
                if (b[j]->atom2->residue && b[j]->atom2->residue != a->residue) continue;

                #if _dbg_path_search
                cout << "Trying " << b[j]->atom2->name << "... ";
                #endif

                l = path_contains_atom(i, b[j]->atom2);
                if (l > 0)
                {
                    if ((m-l) > 1)
                    {
                        Atom* ring_atoms[m];
                        for (h=l-1; h<m; h++)
                        {
                            ring_atoms[h-l] = paths[i][h];
                        }
                        ring_atoms[h-l] = b[j]->atom2;
                        h++;
                        ring_atoms[h-l] = nullptr;

                        add_ring(ring_atoms);
                        #if _dbg_path_search
                        cout << "Created ring from ";
                        Atom::dump_array(ring_atoms);
                        #endif
                    }
                }
                else
                {
                    if (paths[n]) delete[] paths[n];
                    paths[n] = new Atom*[atcount];
                    for (q=0; q<atcount; q++) paths[n][q] = nullptr;
                    copy_path(i, n);
                    paths[n][m] = b[j]->atom2;
                    paths[n][m+1] = nullptr;

                    #if _dbg_path_search
                    cout << "Created ";
                    echo_path(n);
                    #endif

                    n++;
                    if (n >= limit) goto _exit_paths;
                    k++;
                    num_added++;
                }
            }

            #if _dbg_path_search
            cout << endl;
            #endif
        }

        for (j=n-2; j>=0; j--)
        {
            if (path_is_subset_of(j, n-1))
            {
                n--;
                int plj = path_get_length(j);
                int pln = path_get_length(n);
                if (plj == pln) num_added--;

                #if _dbg_path_search
                cout << plj << "/" << pln << " ";
                cout << "Replacing ";
                echo_path(j);
                cout << "...with ";
                echo_path(n);
                cout << endl;
                #endif

                copy_path(n, j);
                delete[] paths[n];
                paths[n] = nullptr;
            }
        }

        if (!num_added) break;
    }

    _exit_paths:
    #if _dbg_path_search
    cout << "Paths:" << endl;
    for (i=0; i<limit && paths[i]; i++) echo_path(i);
    #else
    ;
    #endif
}

void Molecule::identify_acidbase()
{
    if (noAtoms(atoms)) return;

    // For every atom in the molecule:
    int i, j, k, l;
    Bond* b[16];

    for (i=0; atoms[i]; i++)
    {
        // If it is a carbon, pi-bonded to a chalcogen, not bonded to a pnictogen,
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
            atoms[i]->fetch_bonds(b);
            if (!b[0]) goto _not_acidic;

            int nb2 = atoms[i]->get_bonded_atoms_count();
            if ((fama == TETREL || fama == PNICTOGEN) && !atoms[i]->is_pi() && nb2 < 4)
            {
                if (nb2 >= 2)
                {
                    l=0;
                    Point planarity_check[4];
                    planarity_check[l++] = atoms[i]->get_location();
                    for (j=0; j<nb2; j++)
                    {
                        if (!b[j]) break;
                        if (!b[j]->atom2) break;
                        if (!b[j]->atom2->get_Z()) break;
                        planarity_check[l++] = b[j]->atom2->get_location();
                    }

                    bool lplanar = false;
                    if (l > 3)
                    {
                        float coplanarity = are_points_planar(planarity_check[0], planarity_check[1], planarity_check[2], planarity_check[3]);
                        if (coplanarity < coplanar_threshold) lplanar = true;
                    }
                    else if (l == 3)
                    {
                        float theta = find_3d_angle(planarity_check[1], planarity_check[2], planarity_check[0]);
                        if (theta > square && theta < M_PI && fabs(theta-triangular) < fabs(theta - tetrahedral)) lplanar = true;
                        #if _dbg_internal_energy
                        if (!atoms[i]->residue) cout << atoms[i]->name << " " << (theta*fiftyseven) << (lplanar ? " pi" : "") << endl;
                        #endif
                    }

                    if (lplanar)
                    {
                        atoms[i]->aromatize();
                        int chalcogens = 0;
                        for (j=0; j<l; j++)
                        {
                            if (!b[j]) break;
                            if (!b[j]->atom2) break;
                            if (!b[j]->atom2->get_Z()) break;
                            int bfam = b[j]->atom2->get_family();
                            if (bfam == PNICTOGEN || bfam == CHALCOGEN)
                            {
                                b[j]->atom2->aromatize();
                                b[j]->cardinality = 1.5;

                                if (bfam == CHALCOGEN && b[j]->atom2->get_bonded_heavy_atoms_count() < 2)
                                {
                                    chalcogens++;
                                    if (chalcogens > 1 && !b[j]->atom2->get_charge())
                                    {
                                        b[j]->atom2->increment_charge(-1);
                                        chalcogens = -65536;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (carbon)
            {
                if (!atoms[i]->is_pi() || !atoms[i]->is_bonded_to_pi(CHALCOGEN, true))
                {
                    goto _not_acidic;
                }
                for (j=0; b[j]; j++)
                {
                    if (!b[j]->atom2) continue;
                    if (b[j]->cardinality == 2)
                    {
                        int fam = b[j]->atom2->get_family();
                        if (fam != CHALCOGEN)
                        {
                            goto _not_acidic;
                        }
                    }
                }
            }
            for (j=0; b[j]; j++)
            {
                if (!b[j]->atom2) continue;
                int fam = b[j]->atom2->get_family();
                if (carbon && fam == PNICTOGEN)
                {
                    goto _not_acidic;
                }
                //cout << "Fam: " << fam << endl;
                if (fam == CHALCOGEN && b[j]->cardinality < 2)
                {
                    if (b[j]->atom2->get_charge() < 0)
                    {
                        sbOH++;
                        break;
                    }
                    else
                    {
                        Bond* b1[16];
                        b[j]->atom2->fetch_bonds(b1);
                        if (!b1[0])
                        {
                            goto _not_acidic;
                        }
                        for (k=0; b1[k]; k++)
                        {
                            if (!b1[k]->atom2) continue;
                            if (b1[k]->atom2->get_Z() == 1)
                            {
                                sbOH++;
                                //cout << "OH" << endl;
                                break;
                            }
                        }
                    }
                }
            }
        }
        if (sbOH)
        {
            atoms[i]->fetch_bonds(b);
            for (j=0; b[j]; j++)
            {
                if (!b[j]->atom2) continue;
                int fam = b[j]->atom2->get_family();
                if (fam == CHALCOGEN)
                {
                    b[j]->atom2->set_acidbase(-1);
                    //cout << "Atom " << b[j]->atom2->name << " is acidic." << endl;
                }
            }
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
                Atom* bto = c->is_bonded_to("O");
                if (bto)
                {
                    if (amide_zwitterionic_amount)
                    {
                        // Amides are weakly zwitterionic.
                        float arity = c->is_bonded_to(bto);
                        if (arity >= 1.5)
                        {
                            atoms[i]->increment_charge(amide_zwitterionic_amount);
                            bto->increment_charge(-amide_zwitterionic_amount);
                        }
                    }

                    goto _not_basic;
                }
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

Bond** Molecule::get_rotatable_bonds(bool icf)
{
    if (noAtoms(atoms)) return 0;
    if (rotatable_bonds) return rotatable_bonds;
    if (mol_typ == MOLTYP_AMINOACID)
    {
        // TODO: There has to be a better way.
        Star s;
        s.pmol = this;
        if (!rotatable_bonds) rotatable_bonds = s.paa->get_rotatable_bonds();
        else if (rotatable_bonds[0] && rotatable_bonds[1] && rotatable_bonds[0]->atom1 == rotatable_bonds[1]->atom1) rotatable_bonds = s.paa->get_rotatable_bonds();
        else if (rotatable_bonds[0] && abs(rotatable_bonds[0]->atom1 - rotatable_bonds[0]->atom2) >= 524288) rotatable_bonds = s.paa->get_rotatable_bonds();
        return rotatable_bonds;
    }
    // cout << name << " Molecule::get_rotatable_bonds()" << endl << flush;

    Bond* btemp[65536];

    int i,j, bonds=0;
    if (!immobile)
        for (i=0; atoms[i]; i++)
        {
            Bond* lb[32];
            atoms[i]->fetch_bonds(lb);
            int g = atoms[i]->get_geometry();
            for (j=0; j<g && lb[j]; j++)
            {
                if (!lb[j]->atom1 || !lb[j]->atom2) continue;

                if (lb[j]->cardinality > 1 && (!lb[j]->atom1->is_pi() || !lb[j]->atom2->is_pi()) && lb[j]->cardinality < 3)
                    lb[j]->cardinality = 1;

                bool pia = lb[j]->atom1->is_pi(),
                     pib = lb[j]->atom2->is_pi();

                int fa = lb[j]->atom1->get_family(),
                    fb = lb[j]->atom2->get_family();

                if (lb[j]->atom1->in_same_ring_as(lb[j]->atom2))
                {
                    #if _ALLOW_FLEX_RINGS
                    lb[j]->can_rotate = false;
                    lb[j]->compute_flip_capability();
                    #else
                    lb[j]->can_rotate = lb[j]->can_flip = false;
                    continue;
                    #endif
                }

                // Generally, a single bond from a pi atom to an amino group cannot rotate.
                if (pia && pib)
                {
                    lb[j]->can_rotate = false;
                    lb[j]->compute_flip_capability();
                }

                // If atoms a and b are pi, and a-b cannot rotate, then a-b can flip.
                if (!lb[j]->can_rotate
                    && lb[j]->atom2->is_bonded_to(CHALCOGEN)
                    && fa != CHALCOGEN
                    && !(lb[j]->atom1->in_same_ring_as(lb[j]->atom2))
                    )
                {
                    lb[j]->compute_flip_capability();
                }

                if (lb[j]->atom2
                        &&
                        lb[j]->atom1 < lb[j]->atom2
                        &&
                        (lb[j]->can_rotate || (icf && lb[j]->can_flip))
                   )
                {
                    btemp[bonds++] = lb[j];
                    btemp[bonds] = 0;
                }
            }
        }
    else
        for (i=0; atoms[i]; i++)
        {
            Bond* lb[16];
            atoms[i]->fetch_bonds(lb);
            int g = atoms[i]->get_geometry();
            for (j=0; j<g; j++)
            {
                // Generally, a single bond between pi atoms cannot rotate.
                // Same if pi atom bonded to a pnictogen or chalcogen without a single bond to other atoms.
                if (lb[j]->atom1 && lb[j]->atom2
                    &&  (
                            (
                                lb[j]->atom1->is_pi() &&
                                (   lb[j]->atom2->is_pi()
                                    ||
                                    (   
                                        !lb[j]->atom2->is_bonded_to_pi(TETREL, false)
                                        &&
                                        (
                                            lb[j]->atom2->get_family() == PNICTOGEN
                                            ||
                                            lb[j]->atom2->get_family() == CHALCOGEN
                                        )
                                    )
                                )
                            )
                            ||
                            (
                                lb[j]->atom2->is_pi() &&
                                (   lb[j]->atom1->is_pi()
                                    ||
                                    (   
                                        !lb[j]->atom1->is_bonded_to_pi(TETREL, false)
                                        &&
                                        (
                                            lb[j]->atom1->get_family() == PNICTOGEN
                                            ||
                                            lb[j]->atom1->get_family() == CHALCOGEN
                                        )
                                    )
                                )
                            )
                        )
                   )
                    lb[j]->can_rotate = false;

                if ((lb[j]->can_rotate || (icf && lb[j]->can_flip))
                        &&
                        lb[j]->atom1 && lb[j]->atom2
                        &&
                        (!lb[j]->atom1->is_backbone || !strcmp(lb[j]->atom1->name, "CA"))
                        &&
                        !lb[j]->atom2->is_backbone
                        &&
                        greek_from_aname(lb[j]->atom1->name) == (greek_from_aname(lb[j]->atom2->name)-1)
                        &&
                        lb[j]->atom2->get_Z() > 1
                   )
                {
                    btemp[bonds++] = lb[j];
                    btemp[bonds] = 0;
                }
            }
        }

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
        if (b[i]->can_rotate)
        {
            float ltheta = frand(-theta, theta)*pow(frand(0,1),2);
            b[i]->rotate(ltheta);
            if (get_internal_clashes() > int_clsh*4) b[i]->rotate(-ltheta);
        }
        else if (b[i]->can_flip)
        {
            float sint = sin(theta);
            if (frand(0,2) <= sint)
            {
                b[i]->rotate(b[i]->flip_angle);
                if (get_internal_clashes() > int_clsh*4) b[i]->rotate(b[i]->flip_angle);
            }
        }
    }
}

void Molecule::clear_cache()
{
    rotatable_bonds = nullptr;
}

// TODO: There has to be a better way.
Bond** AminoAcid::get_rotatable_bonds()
{
    if (rotatable_bonds) return rotatable_bonds;

    // cout << name << " AminoAcid::get_rotatable_bonds()" << endl << flush;
    // Return ONLY side chain bonds, from lower to higher Greek. E.g. CA-CB but NOT CB-CA.
    // Exclude CA-N and CA-C as these will be managed by the Protein class.
    if (noAtoms(atoms)) return 0;

    // TODO: Something is overwriting the cached rotatable_bonds, causing segfaults.
    // So the cache is unusable for amino acids until the problem gets fixed.
    // if (rotatable_bonds) return rotatable_bonds;
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
        for (i=0; aadef->aabonds[i] && aadef->aabonds[i]->Za && aadef->aabonds[i]->Zb; i++)
        {
            // cout << (name ? name : "(no name)") << "." << *(aadef->aabonds[i]) << endl;
            if (    (
                        aadef->aabonds[i]->cardinality == 1
                        &&
                        (aadef->aabonds[i]->can_rotate || aadef->aabonds[i]->can_flip)
                    )
                    ||
                    (
                        (aadef->aabonds[i]->cardinality < 2 && aadef->aabonds[i]->Za == 6 && aadef->aabonds[i]->Zb == 8)
                        ||
                        (aadef->aabonds[i]->cardinality < 2 && aadef->aabonds[i]->Za == 8 && aadef->aabonds[i]->Zb == 6)
                    )
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

                        Bond* lbb[32];
                        la->fetch_bonds(lbb);
                        if (lbb[0])
                        {
                            cout << la->name << " is bonded to:";
                            int o, ag = la->get_geometry();
                            for (o=0; o<ag; o++)
                                if (lbb[o]
                                    && (abs((__int64_t)(this) - (__int64_t)lbb[o]) < memsanity)
                                    && lbb[o]->atom2
                                    && (abs((__int64_t)(this) - (__int64_t)lbb[o]->atom2) < memsanity)
                                    )
                                    cout << " " << lbb[o]->atom2->name;
                            cout << "." << endl;
                        }

                        Atom* lba = get_atom(aadef->aabonds[i]->bname);
                        if (lba)
                        {
                            lba->fetch_bonds(lbb);
                            if (lbb[0])
                            {
                                cout << lba->name << " is bonded to:";
                                int o, ag = lba->get_geometry();
                                for (o=0; o<ag; o++) if (lbb[o]->atom2) cout << " " << lbb[o]->atom2->name;
                                cout << "." << endl;
                            }
                        }
                        else cout << aadef->aabonds[i]->bname << " not found." << endl;
                    }
                    else
                    {
                        // cout << (name ? name : "(no name)") << ":" << *(lb) << endl;
                        // Generally, a single bond between pi atoms cannot rotate.
                        if (lb->atom1->is_pi() && lb->atom2 && lb->atom2->is_pi())
                        {
                            lb->can_rotate = false;
                            lb->compute_flip_capability();
                            lb->flip_angle = M_PI;
                        }

                        lb->can_rotate = aadef->aabonds[i]->can_rotate;

                        if ((!la->is_backbone || !strcmp(la->name, "CA"))
                                &&
                                la->get_Z() > 1
                                &&
                                (	greek_from_aname(la->name) == (greek_from_aname(lb->atom2->name)+1)
                                    ||
                                    greek_from_aname(la->name) == (greek_from_aname(lb->atom2->name)-1)
                                )
                           )
                        {
                            // cout << "Included." << endl;

                            if (greek_from_aname(la->name) < greek_from_aname(lb->atom2->name))
                                btemp[bonds] = la->get_bond_between(lb->atom2);
                            else
                                btemp[bonds] = lb->atom2->get_bond_between(la);

                            btemp[++bonds] = 0;

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
        Bond* lb[16];
        atoms[i]->fetch_bonds(lb);
        int g = atoms[i]->get_geometry();
        for (j=0; j<g; j++)
        {
            if (!lb[j]) break;
            if (lb[j]->can_rotate
                    &&
                    lb[j]->atom1 && lb[j]->atom2
                    &&
                    (!lb[j]->atom1->is_backbone || !strcmp(lb[j]->atom1->name, "CA"))
                    &&
                    !lb[j]->atom2->is_backbone
                    &&
                    greek_from_aname(lb[j]->atom1->name) == (greek_from_aname(lb[j]->atom2->name)-1)
                    &&
                    lb[j]->atom2->get_Z() > 1
               )
            {
                // cout << *lb[j] << " ";
                btemp[bonds++] = lb[j];
                btemp[bonds] = 0;
            }
            else lb[j]->can_rotate = false;
        }
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
    int i, count = 0;
    float total = 0;
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
        Bond* lb[16];
        atoms[i]->fetch_bonds(lb);
        int g = atoms[i]->get_geometry();
        for (j=0; j<g; j++)
        {
            if (!lb[j]) break;
            if (lb[j]->atom1 < lb[j]->atom2
                    ||
                    !unidirectional
               )
            {
                btemp[bonds++] = lb[j];
                btemp[bonds] = 0;
            }
        }
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
            if (atoms[i]->is_bonded_to(atoms[j]) || atoms[j]->is_bonded_to(atoms[i]))
            {
                Bond* ab = atoms[i]->get_bond_between(atoms[j]);
                if (atoms[i] < atoms[j] && atoms[i]->is_pi() && atoms[j]->is_pi() && !atoms[i]->in_same_ring_as(atoms[j]))
                {
                    // https://laney.edu/corlett/wp-content/uploads/sites/234/2012/01/ch17.pdf
                    Atom* c = atoms[i]->get_heaviest_bonded_atom_that_isnt(atoms[j]);
                    Atom* d = atoms[j]->get_heaviest_bonded_atom_that_isnt(atoms[i]);

                    if (c && d)
                    {
                        SCoord axis = atoms[j]->get_location().subtract(atoms[i]->get_location());
                        float theta = find_angle_along_vector(c->get_location(), d->get_location(), atoms[i]->get_location(), axis);
                        float cpartial = 13.5 - 13.5 * cos(theta*2);
                        #if _dbg_internal_energy
                        if (!is_residue())
                            cout << "Conjugated " << atoms[i]->name << "-" << atoms[j]->name
                                << " " << ab->cardinality << " bond theta = " << (theta*fiftyseven) << "deg."
                                << " adding " << cpartial << " kJ/mol."
                                << endl;
                        clash += cpartial;
                        #endif
                    }
                }
                continue;
            }

            Point ptb = atoms[j]->get_location();
            float bvdW = atoms[j]->get_vdW_radius();

            r = pta.get_3d_distance(&ptb);
            if (r >= 0.9 && atoms[i]->shares_bonded_with(atoms[j])) continue;

            if (!r) r += 10e-15;
            if (r < avdW + bvdW)
            {
                float lclash = fmax(InteratomicForce::Lennard_Jones(atoms[i], atoms[j]), 0); // sphere_intersection(avdW, bvdW, r);
                clash += lclash;

                if (false && lclash > 3)
                {
                    cout << atoms[i]->name << " clashes with " << atoms[j]->name << " by " << lclash << " cu. A. resulting in " << clash << endl;
                    int g = atoms[i]->get_geometry();
                    cout << "Geometry: " << g << endl;
                    Bond* b[16];
                    atoms[i]->fetch_bonds(b);
                    int k;
                    for (k=0; k<g; k++)
                        cout << atoms[i]->name << " is bonded to " << hex << b[k]->atom2 << dec << " "
                             << (b[k]->atom2 ? b[k]->atom2->name : "") << "." << endl;
                }
            }
        }
    }

    return clash; // -base_internal_clashes;
}

#if compute_vdw_repulsion
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
                    retval += fabs(sphere_intersection(halflim, halflim, r) * kJmol / asphere);
                }
            }
        }
    }

    return retval;
}
#endif

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

    clash1 = clash2 = nullptr;
    float worst = 0;

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

                if (atoms[i]->is_backbone && ligands[l]->atoms[j]->is_backbone && abs(atoms[i]->residue - ligands[l]->atoms[j]->residue) < 2) continue;

                float f = fmax(InteratomicForce::Lennard_Jones(atoms[i], ligands[l]->atoms[j]), 0);

                if (f > worst)
                {
                    worst = f;
                    clash1 = atoms[i];
                    clash2 = ligands[l]->atoms[j];
                }

                if (f > worst_mol_clash)
                {
                    worst_mol_clash = f;
                    worst_clash_1 = this;
                    worst_clash_2 = ligands[l];
                }

                clash += f;
                continue;
            }
        }
    }

    return clash;
}

float Molecule::total_intermol_clashes(Molecule** ligands)
{
    if (!ligands) return 0;
    int i;
    float clash = 0;
    for (i=0; ligands[i]; i++)
    {
        clash += ligands[i]->get_intermol_clashes(ligands);
    }
    return clash;
}

void Molecule::mutual_closest_atoms(Molecule* mol, Atom** a1, Atom** a2)
{
    if (!a1 || !a2) return;

    *a1 = *a2 = nullptr;

    int i, j, m, n;
    Atom *a, *b;
    float rbest = Avogadro;

    m = get_atom_count();
    n = mol->get_atom_count();
    for (i=0; i<m; i++)
    {
        a = get_atom(i);
        int aZ = a->get_Z();
        for (j=0; j<n; j++)
        {
            b = mol->get_atom(j);
            int bZ = b->get_Z();
            if (aZ == 1 && bZ == 1) continue;

            float r = a->distance_to(b);
            if (r < rbest)
            {
                rbest = r;
                *a1 = a;
                *a2 = b;
            }
        }
    }
}

void Molecule::move(SCoord move_amt, bool override_residue)
{
    if (noAtoms(atoms)) return;
    if (immobile)
    {
        cout << "Warning: Attempt to move \"immobile\" molecule " << name << endl;
        return;
    }
    int i;
    int vvc = vdw_vertex_count;

    #if _dbg_improvements_only_rule
    float before;
    if (check_ligand && check_mols && !excuse_deterioration) before = cfmol_multibind(check_ligand, check_mols);
    #endif

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->residue && !override_residue) return;
        Point loc = atoms[i]->get_location();
        loc = loc.add(&move_amt);
        atoms[i]->move(&loc);
    }

    #if _dbg_improvements_only_rule
    if (check_ligand && check_mols && !excuse_deterioration) 
    {
        float after = cfmol_multibind(check_ligand, check_mols);
        if (after < before) throw 0xb0661;
    }
    excuse_deterioration = false;
    #endif

    vdw_vertex_count = vvc;
}

void Molecule::move(Point move_amt, bool override_residue)
{
    if (noAtoms(atoms)) return;
    if (immobile)
    {
        cout << "Warning: Attempt to move \"immobile\" molecule " << name << endl;
        return;
    }
    int i;
    int vvc = vdw_vertex_count;

    #if _dbg_improvements_only_rule
    float before;
    if (check_ligand && check_mols && !excuse_deterioration) before = cfmol_multibind(check_ligand, check_mols);
    #endif

    for (i=0; atoms[i]; i++)
    {
        // cout << atoms[i]->name << " ";
        if (atoms[i]->residue && !override_residue) return;
        Point loc = atoms[i]->get_location();
        loc = loc.add(&move_amt);
        atoms[i]->move(&loc);
    }

    #if _dbg_improvements_only_rule
    if (check_ligand && check_mols && !excuse_deterioration) 
    {
        float after = cfmol_multibind(check_ligand, check_mols);
        if (after < before) throw 0xb0661;
    }
    excuse_deterioration = false;
    #endif

    vdw_vertex_count = vvc;
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

    for (i=0; atoms[i]; i++)
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

float Molecule::get_charge() const
{
    int i;
    float charge=0;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->get_Z() == 1) continue;
        charge += atoms[i]->get_charge();
    }
    return charge;
}

void Molecule::recenter(Point nl)
{
    if (movability <= MOV_NORECEN) return;
    Point loc = get_barycenter();
    Point rel = nl.subtract(&loc);
    SCoord v(&rel);
    move(v);
}

void Molecule::rotate(SCoord* SCoord, float theta, bool bond_weighted)
{
    if (noAtoms(atoms)) return;
    int vvc = vdw_vertex_count;
    // cout << name << " Molecule::rotate()" << endl;

    if (movability <= MOV_FLEXONLY) return;
    if (movability <= MOV_NORECEN) bond_weighted = false;
    Point cen = get_barycenter(bond_weighted);

    #if _dbg_improvements_only_rule
    float before;
    if (check_ligand && check_mols && !excuse_deterioration) before = cfmol_multibind(check_ligand, check_mols);
    #endif

    int i;
    if (fabs(theta) > hexagonal)
    {
        i = 0;
    }
    for (i=0; atoms[i]; i++)
    {
        // if (atoms[i]->residue) return;
        Point loc = atoms[i]->get_location();
        Point nl  = rotate3D(&loc, &cen, SCoord, theta);
        atoms[i]->move(&nl);
    }

    #if _dbg_improvements_only_rule
    if (check_ligand && check_mols && !excuse_deterioration) 
    {
        float after = cfmol_multibind(check_ligand, check_mols);
        if (after < before) throw 0xb0661;
    }
    excuse_deterioration = false;
    #endif

    vdw_vertex_count = vvc;
}

void Molecule::rotate(LocatedVector lv, float theta)
{
    if (noAtoms(atoms)) return;

    if (movability <= MOV_FLEXONLY) return;
    if (movability <= MOV_NORECEN) lv.origin = get_barycenter();
    int vvc = vdw_vertex_count;

    #if _dbg_improvements_only_rule
    float before;
    if (check_ligand && check_mols && !excuse_deterioration) before = cfmol_multibind(check_ligand, check_mols);
    #endif

    int i;
    if (fabs(theta) > hexagonal)
    {
        i = 0;
    }
    for (i=0; atoms[i]; i++)
    {
        // if (atoms[i]->residue) return;
        Point loc = atoms[i]->get_location();
        Point nl  = rotate3D(&loc, &lv.origin, &lv, theta);
        atoms[i]->move(&nl);
    }

    #if _dbg_improvements_only_rule
    if (check_ligand && check_mols && !excuse_deterioration) 
    {
        float after = cfmol_multibind(check_ligand, check_mols);
        if (after < before) throw 0xb0661;
    }
    excuse_deterioration = false;
    #endif

    vdw_vertex_count = vvc;
}

bool Molecule::shielded(Atom* a, Atom* b) const
{
    int i;
    float r = a->distance_to(b);
    float r6 = r*1.26, r125 = 1.25*r;
    if (r < 2) return false;

    a->shielding_angle = b->shielding_angle = 0;

    Point aloc = a->get_location(), bloc = b->get_location();
    for (i=0; atoms[i]; i++)
    {
        Atom* ai = atoms[i];
        if (!ai) break;
        if (ai == a || ai == b) continue;
        float rai = ai->distance_to(a);
        if (rai > r6) continue;
        float rbi = ai->distance_to(b);
        if (rbi > r6) continue;
        if ((rai+rbi) > r125) continue;
        Point sloc = ai->get_location();
        float f3da = find_3d_angle(&aloc, &bloc, &sloc);
        if (f3da > a->shielding_angle) a->shielding_angle = b->shielding_angle = f3da;
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

float Molecule::pi_stackability(bool ib)
{
    if (noAtoms(atoms)) return 0;
    int i, j=0;
    float result = 0;

    for (i=0; atoms[i]; i++)
    {
        if (!ib && atoms[i]->is_backbone) continue;
        if (atoms[i]->is_pi()) result += 1;
        j++;
    }

    if (j) result /= j;
    return result;
}

float Molecule::get_atom_mol_bind_potential(Atom* a)
{
    if (noAtoms(atoms)) return 0;
    int i, j, n;

    float hydro = 0;
    j = 0;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        if (atoms[i]->get_Z() == 1) continue;

        hydro += fabs(atoms[i]->is_polar());
        j++;
    }

    if (j) hydro /= j;

    float retval=0;
    n = 0;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;

        InteratomicForce* ifs[32];
        InteratomicForce::fetch_applicable(a, atoms[i], ifs);
        if (!ifs) continue;

        for (j=0; ifs[j]; j++)
        {
            float partial;

            if (hydro > hydrophilicity_cutoff && ifs[j]->get_type() == vdW) continue;

            if (ifs[j]->get_type() == ionic)
            {
                if (sgn(a->get_charge()) != -sgn(atoms[i]->get_charge())) continue;
                partial = 60;
            }
            else
            {
                partial = ifs[j]->get_kJmol();
            }

            if (ifs[j]->get_type() == hbond)
            {
                partial *= fmin(fabs(a->is_polar()), fabs(atoms[i]->is_polar()));
            }

            if (ifs[j]->get_type() == polarpi) partial /= 6;            // Config is for benzene rings.

            if (ifs[j]->get_type() == mcoord)
            {
                partial *= (1.0 + 1.0 * cos((a->get_electronegativity() + atoms[i]->get_electronegativity()) / 2 - 2.25));
            }

            retval += partial;
            n++;
        }
    }

    if (n) retval /= n;

    return retval;
}

float Molecule::get_intermol_binding(Molecule* ligand, bool subtract_clashes)
{
    Molecule* ligands[4];
    ligands[0] = ligand;
    ligands[1] = nullptr;
    return get_intermol_binding(ligands, subtract_clashes);
}

void Molecule::clear_atom_binding_energies()
{
    if (!atoms) return;
    int i;
    for (i=0; atoms[i]; i++)
        atoms[i]->last_bind_energy = 0;
}

float Molecule::get_intermol_potential(Molecule* ligand, bool pure)
{
    Molecule* ligands[4];
    ligands[0] = ligand;
    ligands[1] = nullptr;
    return get_intermol_potential(ligands, pure);
}

float Molecule::get_intermol_potential(Molecule** ligands, bool pure)
{
    if (!atoms) return 0;
    if (!ligands) return 0;
    if (!ligands[0]) return 0;
    int i, j, l, n;
    float kJmol = 0;

    for (i=0; atoms[i]; i++)
    {
        if (!atoms[i]) continue;
        Point aloc = atoms[i]->get_location();
        for (l=0; ligands[l]; l++)
        {
            if (ligands[l] == this) continue;
            for (j=0; j<ligands[l]->atcount; j++)
            {
                if (!ligands[l]->atoms[j]) continue;
                float r = ligands[l]->atoms[j]->get_location().get_3d_distance(&aloc);
                float f = 1.0 / r;		// Regular invert rather than inv square so that actual bonding will take over at short range.
                InteratomicForce* iff[32];
                InteratomicForce::fetch_applicable(atoms[i], ligands[l]->atoms[j], iff);

                if (iff) for (n=0; iff[n]; n++)
                {
                    if (iff[n]->get_type() == vdW) continue;

                    if (pure || r < iff[n]->get_distance())
                        kJmol += iff[n]->get_kJmol();
                    else
                        kJmol += iff[n]->get_kJmol()*f;
                }
            }
        }
    }

    return kJmol;
}

float Molecule::get_intermol_binding(Molecule** ligands, bool subtract_clashes)
{
    if (!ligands) return 0;
    if (!ligands[0]) return 0;
    int i, j, l;
    float kJmol = 0;
    if (!atoms) return 0;
    if (subtract_clashes) kJmol -= get_internal_clashes();
    nearest_r = Avogadro;
    nearest_local_aidx = 0;

    lastshielded = 0;
    clash1 = clash2 = nullptr;
    float best_atom_energy = 0;

    #if _dbg_internal_energy
    cout << (name ? name : "") << " base internal clashes: " << base_internal_clashes << "; final internal clashes " << -kJmol << endl;
    #endif

    for (i=0; atoms[i]; i++)
    {
        atoms[i]->last_bind_energy = 0;
        atoms[i]->strongest_bind_energy = 0;
        atoms[i]->strongest_bind_atom = nullptr;
    }

    clash_worst = 0;
    for (i=0; atoms[i]; i++)
    {
        Point aloc = atoms[i]->get_location();
        int Z = atoms[i]->get_Z();
        for (l=0; ligands[l]; l++)
        {
            #if _dbg_51e2_ionic
            if (!is_residue() && atoms[i]->get_family() == CHALCOGEN && ligands[l]->is_residue() == 262)
            {
                j = 0;
            }
            #endif

            if (ligands[l] == this) continue;
            for (j=0; j<ligands[l]->atcount; j++)
            {
                // TODO: Fix this in the hydrogenate function, but for now we'll fix it here and hope for the best. 🤞🏼
                if (!ligands[l]->atoms[j])
                {
                    ligands[l]->atcount = j;
                    break;
                }
                if (atoms[i]->is_backbone && ligands[l]->atoms[j]->is_backbone
                    &&
                    (	(	atoms[i]->residue == ligands[l]->atoms[j]->residue - 1
                            &&
                            !strcmp(atoms[i]->name, "C")
                            &&
                            !strcmp(ligands[l]->atoms[j]->name, "N")
                        )
                        ||
                        (	atoms[i]->residue == ligands[l]->atoms[j]->residue + 1
                            &&
                            !strcmp(atoms[i]->name, "N")
                            &&
                            !strcmp(ligands[l]->atoms[j]->name, "C")
                        )
                    )) continue;			// kludge to prevent adjacent residue false clashes.
                float r = ligands[l]->atoms[j]->get_location().get_3d_distance(&aloc);

                if (Z>0 && r<nearest_r)
                {
                    nearest_r = r;
                    nearest_local_aidx = i;
                    nearest_local_atom = atoms[i];
                    nearest_remote_atom = ligands[l]->atoms[j];
                }

                if (r < _INTERA_R_CUTOFF)
                {
                    if (	!shielded(atoms[i], ligands[l]->atoms[j])
                            &&
                            !ligands[l]->shielded(atoms[i], ligands[l]->atoms[j])
                       )
                    {
                        missed_connection.r = -Avogadro;
                        mc_bpotential = 0;
                        float abind = InteratomicForce::total_binding(atoms[i], ligands[l]->atoms[j]);
                        #if _dbg_internal_energy
                        if (ligands[l] == this)
                        {
                            cout << "Energy between " << atoms[i]->name << "..." << ligands[l]->atoms[j]->name
                                << " = " << abind << " kJ/mol." << endl;
                        }
                        #endif
                        if (abind && !isnan(abind) && !isinf(abind))
                        {
                            if (abind > 0 && minimum_searching_aniso && ligands[l]->priority) abind *= 1.5;
                            kJmol += abind;

                            if (abind > best_atom_energy)
                            {
                                best_atom_energy = abind;
                                best_intera = atoms[i];
                                best_interactor = ligands[l];
                                best_other_intera = ligands[l]->atoms[j];
                            }

                            atoms[i]->last_bind_energy += abind;
                            if (abind > atoms[i]->strongest_bind_energy)
                            {
                                atoms[i]->strongest_bind_energy = abind;
                                atoms[i]->strongest_bind_atom = ligands[l]->atoms[j];
                            }

                            if (abind < 0 && -abind > clash_worst)
                            {
                                clash_worst = -abind;
                                clash1 = atoms[i];
                                clash2 = ligands[l]->atoms[j];
                            }

                            if (abind < 0 && ligands[l]->is_residue() && movability >= MOV_ALL)
                            {
                                Point ptd = aloc.subtract(ligands[l]->atoms[j]->get_location());
                                ptd.multiply(fmin(fabs(-abind) / 1000, 1));
                                if (ptd.magnitude() > speed_limit) ptd.scale(speed_limit);
                                lmx += lmpush * sgn(ptd.x);
                                lmy += lmpush * sgn(ptd.y);
                                lmz += lmpush * sgn(ptd.z);
                            }
                        }

                        if (missed_connection.r > 0)
                        {
                            Point mc = missed_connection;
                            // cout << mc << endl;
                            if (ligands[l]->priority) mc_bpotential *= 3.333;
                            float lc = ligands[l]->atoms[j]->get_charge();
                            if (lc && sgn(lc) == -sgn(atoms[i]->get_charge())) mc_bpotential *= 3.333;
                            float mcrr = missed_connection.r * missed_connection.r;
                            lmx += lmpull * mc.x * mc_bpotential / mcrr;
                            lmy += lmpull * mc.y * mc_bpotential / mcrr;
                            lmz += lmpull * mc.z * mc_bpotential / mcrr;
                            lc = sqrt(lmx*lmx+lmy*lmy+lmz*lmz);
                            if (lc > speed_limit)
                            {
                                lc = speed_limit / lc;
                                lmx *= lc;
                                lmy *= lc;
                                lmz *= lc;
                            }
                        }
                    }
                    else lastshielded += InteratomicForce::total_binding(atoms[i], ligands[l]->atoms[j]);
                }
            }
        }
    }
    // cout << "Total: " << kJmol << endl;
    // cout << endl;

    return kJmol;
}

float Molecule::get_intermol_contact_area(Molecule* ligand, bool hpho)
{
    if (!ligand) return 0;
    if (!atoms) return 0;
    int i, j;
    float result = 0;

    for (i=0; atoms[i]; i++)
    {
        Point aloc = atoms[i]->get_location();
        float avdw = atoms[i]->get_vdW_radius();

        if (hpho && atoms[i]->is_polar()) continue;

        for (j=0; j<ligand->atcount; j++)
        {
            if (hpho && ligand->atoms[j]->is_polar()) continue;

            float r = ligand->atoms[j]->get_location().get_3d_distance(&aloc);
            if (!r) continue;
            float bvdw = ligand->atoms[j]->get_vdW_radius();
            if (r > avdw+bvdw) continue;

            float f = sphere_inter_area(avdw, bvdw, r);
            if (!isnan(f)) result += f;
        }
    }

    return result;
}

float Molecule::distance_to(Molecule* om)
{
    if (!atoms || !om || !om->atoms) return nanf("Bad molecule.");

    int i, j;
    float minr = Avogadro;

    for (i=0; atoms[i]; i++)
    {
        for (j=0; om->atoms[j]; j++)
        {
            float r = atoms[i]->distance_to(om->atoms[j]);

            if (r < minr) minr = r;
        }
    }

    return minr;
}

float Molecule::get_intermol_polar_sat(Molecule* ligand)
{
    if (!ligand) return 0;
    if (!atoms) return 0;
    int i, j, l, n;
    float result = 0;

    for (i=0; atoms[i]; i++)
    {
        Point aloc = atoms[i]->get_location();
        float aapol = fabs(atoms[i]->is_polar());
        for (j=0; j<ligand->atcount; j++)
        {
            float r = ligand->atoms[j]->get_location().get_3d_distance(&aloc);
            if (r > 5) continue;
            float abpol = fabs(ligand->atoms[j]->is_polar());
            float f = (aapol >= 0.5 ? 2 : 1) + (abpol >= 0.5 ? 2 : 1);
            if (!f) continue;
            f = 1.0 / pow(fmax(1,r/f),2);

            if (aapol < 0.2)
            {
                if (abpol < 0.2) result += f;
                else if (abpol >= 0.5) result -= f;
            }
            else if (aapol >= 0.5)
            {
                if (abpol < 0.2) result -= f;
                else if (abpol >= 0.5) result += f;
            }
        }
    }

    return result;
}

void Molecule::minimize_internal_clashes()
{
    // cout << (name ? name : "(no name)");
    if (noAtoms(atoms)) return;
    base_internal_clashes = 0;

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

void Molecule::do_histidine_flip(HistidineFlip* hf)
{
    if (hf->N1->get_charge() >= 0.2 || hf->N2->get_charge() >= 0.2) return;

    Point ptC  = hf->C->get_location();
    Point ptN1 = hf->N1->get_location();
    Point ptN2 = hf->N2->get_location();
    Point ptH  = hf->H->get_location();

    Point arr[2] = {ptN1, ptN2};
    Point Navg = average_of_points(arr, 2);

    Point newloc = rotate3D(ptH, Navg, ptC.subtract(Navg), M_PI);
    hf->H->move(newloc);

    Atom* was_bonded = hf->H->get_bond_by_idx(0)->atom2;
    hf->H->unbond(was_bonded);
    float rN1 = newloc.get_3d_distance(ptN1), rN2 = newloc.get_3d_distance(ptN2);
    Atom* new_bonded = (rN1 > rN2) ? hf->N2 : hf->N1;
    hf->H->bond_to(new_bonded, 1);
    strcpy(hf->H->name+1, new_bonded->name+1);

    #if _DBG_HISFLIP
    cout << hf->H->name << " moved from " << ptH << " to " << newloc << endl;
    #endif
}

float Molecule::get_springy_bond_satisfaction()
{
    if (!springy_bonds) return 0;

    float retval = 0;
    int i;
    for (i=0; i<springy_bondct; i++)
    {
        if (!springy_bonds[i].atom1 || !springy_bonds[i].atom2 || !springy_bonds[i].optimal_radius) continue;
        float r = springy_bonds[i].atom1->get_location().get_3d_distance(springy_bonds[i].atom2->get_location());
        r /= springy_bonds[i].optimal_radius;
        if (r < 1) retval += r;
        else retval += 1.0/(r*r);
    }

    return retval;
}

void Molecule::allocate_mandatory_connections(int mcmax)
{
    delete_mandatory_connections();
    mandatory_connection = new Molecule*[mcmax+4];
    mandatory_connection[0] = nullptr;
    last_mc_binding = new float[mcmax+4];
    int i;
    for (i=0; i<mcmax; i++) last_mc_binding[i] = 0;
}

void Molecule::add_mandatory_connection(Molecule* addmol)
{
    int i;
    for (i=0; mandatory_connection[i]; i++)                 // get count.
    {
        if (mandatory_connection[i] == addmol)
        {
            #if _dbg_mand_conn
            cout << "Already have mandatory connection " << addmol->name;
            if (last_mc_binding) cout << " last binding " << last_mc_binding[i];
            cout << endl;
            #endif
            return;      // already have it.
        }
    }
    mandatory_connection[i] = addmol;
    mandatory_connection[i+1] = nullptr;
    if (last_mc_binding) last_mc_binding[i] = 0;
    #if _dbg_mand_conn
    cout << "Add mandatory connection " << addmol->name << endl;
    #endif
}

void Molecule::zero_mandatory_connection_cache()
{
    if (!last_mc_binding) return;
    if (!mandatory_connection) return;
    int i;
    for (i=0; mandatory_connection[i]; i++) last_mc_binding[i] = 0;
}

void Molecule::remove_mandatory_connection(Molecule* rmvmol)
{
    if (!last_mc_binding) return;
    if (!mandatory_connection) return;
    int i, j;
    for (i=0; mandatory_connection[i]; i++)
    {
        if (mandatory_connection[i] == rmvmol)
        {
            for (j=i; mandatory_connection[j]; j++)
            {
                mandatory_connection[j] = mandatory_connection[j+1];
                last_mc_binding[j] = last_mc_binding[j+1];
            }
            #if _dbg_mand_conn
            cout << "Remove mandatory connection " << rmvmol->name << endl;
            #endif
            return;
        }
    }
}

void Molecule::delete_mandatory_connections()
{
    if (last_mc_binding) delete last_mc_binding;
    last_mc_binding = nullptr;
    if (mandatory_connection) delete mandatory_connection;
    mandatory_connection = nullptr;
    #if _dbg_mand_conn
    cout << "No more mandatory connections." << endl;
    #endif
}

float Molecule::intermol_bind_for_multimol_dock(Molecule* om, bool is_ac)
{
    float lbias = 1.0 + (sgn(is_residue()) == sgn(om->is_residue()) ? 0 : dock_ligand_bias);
    float rawbind = get_intermol_binding(om, !is_ac);
    float lbind = rawbind * lbias;
    // if (!is_residue() && om->is_residue()) lbind += get_intermol_polar_sat(om) * polar_sat_influence_for_dock;
    lbind += get_intermol_contact_area(om, true) * cavity_stuffing;
    float clashes = get_intermol_clashes(om);
    lbind -= clashes * iteration_additional_clash_coefficient;

    int i, j, m, n;
    if (mandatory_connection && rawbind >= 0)                   // Allow pullaway if mols are clashing.
    {
        if (!last_mc_binding)
        {
            for (i=0; mandatory_connection[i]; i++);            // Get count.
            last_mc_binding = new float[i+4];
            for (i=0; mandatory_connection[i]; i++) last_mc_binding[i] = -Avogadro;
        }
        for (i=0; mandatory_connection[i]; i++)
        {
            if (mandatory_connection[i] == om)
            {
                if (lbind < last_mc_binding[i] && lbind < mandatory_coordination_threshold)
                {
                    return -1e9;
                }
                else last_mc_binding[i] = lbind;
            }
        }
    }

    n = atcount;
    m = om->atcount;
    for (i=0; i<n; i++)
    {
        if (!atoms[i])
        {
            atcount = i;
            break;
        }
        if (atoms[i]->get_Z() > 0) continue;
        if (fabs(atoms[i]->is_polar()) < hydrophilicity_cutoff) continue;
        Bond* b = atoms[i]->get_bond_by_idx(0);
        if (!b) continue;
        if (!b->atom2) continue;
        if (!b->can_rotate && !b->can_flip) continue;
        b = b->atom2->get_bond_by_idx(0);
        if (!b) continue;
        b = b->get_reversed();
        if (!b) continue;
        if (!b->atom2) continue;

        for (j=0; j<m; j++)
        {
            if (om->atoms[j]->get_Z() < 2) continue;
            if (fabs(om->atoms[j]->is_polar()) < hydrophilicity_cutoff) continue;
            float r = om->atoms[j]->distance_to(atoms[i]);
            if (r > 4) continue;

            Atom* H = om->atoms[j]->is_bonded_to("H");
            float theta = find_angle_along_vector(atoms[i]->get_location(),
                om->atoms[j]->get_location(), b->atom2->get_location(), b->get_axis());
            if (H && atoms[i]->distance_to(H) < r) theta += M_PI;
            if (b->can_rotate)
            {
                b->rotate(theta, false, true);
            }
            else
            {
                if (fabs(theta) > square) b->rotate(M_PI);
            }

            break;
        }
    }

    return lbind;
}

float Molecule::cfmol_multibind(Molecule* a, Molecule** nearby)
{
    static int nearbya[256];
    static float nearbyr[256];
    int lused = rand();

    a->occlusion = 0;
    if (a->is_residue() && ((AminoAcid*)a)->conditionally_basic()) ((AminoAcid*)a)->set_conditional_basicity(nearby);

    float result = -a->total_eclipses();
    if (a->is_residue()) result += reinterpret_cast<AminoAcid*>(a)->initial_eclipses;

    int j;
    for (j=0; nearby[j]; j++)
    {
        float f = a->intermol_bind_for_multimol_dock(nearby[j], false);
        result += f;
        if (nearest_local_aidx < 256)
        {
            nearbya[nearest_local_aidx] = lused;
            nearbyr[nearest_local_aidx] = nearest_r;
        }
    }
    if (a->mclashables)
    {
        for (j=0; a->mclashables[j]; j++)
        {
            float f = a->intermol_bind_for_multimol_dock(a->mclashables[j], false);
            result += f;
        }
    }

    for (j=0; j<a->atcount; j++) if (nearbya[j] == lused) a->occlusion += 1.0 / fmax(1, nearbyr[j]/4);
    a->occlusion /= max(a->get_heavy_atom_count(), 1);

    return result;
}

void Molecule::conform_molecules(Molecule** mm, Molecule** bkg, int iters, void (*cb)(int, Molecule**),
    void (*group_realign)(Molecule*, std::vector<std::shared_ptr<GroupPair>>),
    void (*progress)(float))
{
    int m, n;

    if (!mm) m=0;
    else for (m=0; mm[m]; m++);         // Get count.

    if (!bkg) n=0;
    else for (n=0; bkg[n]; n++);        // Get count.

    Molecule* all[m+n+8];
    int i, j, l=0;

    for (i=0; i<m; i++)
    {
        bool duplicate = false;
        for (j=0; j<l; j++)
        {
            if (all[j] == mm[i]) duplicate = true;
        }
        if (duplicate) continue;

        all[l++] = mm[i];
        mm[i]->movability = static_cast<MovabilityType>(static_cast<int>(mm[i]->movability & !MOV_BKGRND));
    }

    for (i=0; i<n; i++)
    {
        bool duplicate = false;
        for (j=0; j<l; j++)
        {
            if (all[j] == bkg[i]) duplicate = true;
        }
        if (duplicate) continue;

        all[l++] = bkg[i];
        bkg[i]->movability = static_cast<MovabilityType>(static_cast<int>(bkg[i]->movability | MOV_BKGRND));
    }

    all[l] = nullptr;

    conform_molecules(all, iters, cb, group_realign);

    for (i=0; i<n; i++)
    {
        bkg[i]->movability = static_cast<MovabilityType>(static_cast<int>(bkg[i]->movability & !MOV_BKGRND));
    }
}

void Molecule::conform_atom_to_location(const char* an, Point t, int iters, float od)
{
    if (!atoms) return;
    int i;
    for (i=0; atoms[i]; i++)
    {
        if (!strcmp(atoms[i]->name, an))
        {
            conform_atom_to_location(i, t, iters, od);
            return;
        }
    }
}
    
void Molecule::maintain_contact(Atom* my_atom, Atom* other_atom)
{
    int i = atom_idx_from_ptr(my_atom);
    if (i < 0) throw 0xbadda7a;
    dlt1 = my_atom;
    dlt2 = other_atom;
    dltr = InteratomicForce::optimal_distance(dlt1, dlt2);
}

bool Molecule::contact_maintained()
{
    if (!dlt1 || !dlt2) return true;

    float r = dlt1->distance_to(dlt2);
    if (dlt1->mol && dlt2->mol)
    {
        float c = dlt1->mol->get_intermol_clashes(dlt2->mol);
        if (c > clash_limit_per_aa) dltr *= contact_maintenance_creep;
    }
    else
    {
        float e = InteratomicForce::total_binding(dlt1, dlt2);
        if (e > clash_limit_per_atom) dltr *= contact_maintenance_creep;
    }

    // cout << dlt1->name << " ~ " << dlt2->name << " are " << r << "A apart, limit " << dltr << endl << endl;
    return r <= dltr*contact_maintenance_allowance;
}

const Point* Molecule::obtain_vdW_surface(float d)
{
    if (!atcount || !atoms) return nullptr;
    if (vdw_surface && vdw_vertex_count > 0) return vdw_surface;

    int maxpoints = atcount * d * d / 3 + 256;
    if (!vdw_surface)
    {
        vdw_surface = new Point[maxpoints];
        vdw_vertex_atom = new Atom*[maxpoints];
    }

    float halfstep = M_PI / d;
    float step = halfstep * 2;

    int i, ivdW = 0;
    SCoord v;
    for (i=0; i<atcount && atoms[i]; i++)
    {
        Point aloc = atoms[i]->get_location();
        v.r = atoms[i]->get_vdW_radius();
        float ystep = step / v.r / v.r;
        for (v.theta = -square; v.theta <= square; v.theta += step)
        {
            float xstep = step / v.r / fmax(cos(v.theta), 0.000001);
            float end = M_PI*2-xstep/2;
            for (v.phi = 0; v.phi < end; v.phi += xstep)
            {
                Point pt = aloc.add(v);
                Atom* na = this->get_nearest_atom(pt);
                if (na != atoms[i] && pt.get_3d_distance(na->get_location()) < na->get_vdW_radius()) continue;
                if (!pt.x && !pt.y && !pt.z) pt = Point(-0.001, 0.001, -0.001);
                vdw_vertex_atom[ivdW] = atoms[i];
                vdw_surface[ivdW++] = pt;
                if (ivdW >= maxpoints)
                {
                    cout << "Too many vdW surface vertices. Please increase limit in code." << endl;
                    throw 0xbadc0de;
                }
            }
        }
    }
    vdw_vertex_count = ivdW;

    return vdw_surface;
}

void Molecule::shape_has_changed()
{
    vdw_vertex_count = 0;
}

void mol_shape_has_changed(Molecule* mol)
{
    mol->shape_has_changed();
}

#define _dbg_atom_pointing 0

void Molecule::conform_atom_to_location(int i, Point t, int iters, float od)
{
    if (!(movability & MOV_CAN_FLEX)) return;

    int iter, j, l, circdiv = 18;
    Bond** b = get_rotatable_bonds();
    if (!b) return;
    Atom* a = atoms[i];
    if (!a) return;

    float oc = get_internal_clashes();

    for (iter = 0; iter < iters; iter++)
    {
        float r;
        for (j=0; b[j]; j++)
        {
            float bestr = Avogadro, bestth = 0;
            for (l=0; l<=circdiv; l++)
            {
                b[j]->rotate(M_PI*2.0/circdiv, false, true);
                r = a->get_location().get_3d_distance(t);
                if (od) r = fabs(r-od);
                float c = get_internal_clashes();
                if (r < bestr && c < oc+2.5*clash_limit_per_aa)
                {
                    bestr = r;
                    bestth = M_PI*2.0/circdiv * l;
                }
            }
            b[j]->rotate(bestth);

            #if _dbg_atom_pointing
            cout << name << "." << iter << ": " << circdiv << "|" << bestr << endl;
            #endif
        }
        if (!(iter % 3)) circdiv++;
    }
}

float Molecule::total_intermol_binding(Molecule** l)
{
    int i;
    float f = 0;

    for (i=0; l[i]; i++)
    {
        f += l[i]->get_intermol_binding(l);
    }

    return f;
}

void Molecule::conform_molecules(Molecule** mm, int iters, void (*cb)(int, Molecule**),
    void (*group_realign)(Molecule*, std::vector<std::shared_ptr<GroupPair>>),
    void (*progress)(float))
{
    if (!mm) return;
    int i, j, l, n, iter;

    minimum_searching_aniso = 0.5;

    #if _dbg_improvements_only_rule
    check_mols = mm;
    #endif

    for (iter=0; iter<iters; iter++)
    {
        for (i=0; mm[i]; i++) mm[i]->lastbind = 0;

        for (n=0; mm[n]; n++);      // Get count.
        Molecule* nearby[n+8];
        bool do_full_rotation = _allow_fullrot && ((iter % _fullrot_every) == 0);

        for (i=0; i<n; i++)
        {
            Molecule* a = mm[i];
            bool flipped_rings = false;

            if (!a->iterbegan) a->iterbegan = new Pose(a);
            a->iterbegan->copy_state(a);
            if (!iter) a->iters_without_change = 0;

            #if _dbg_asunder_atoms
            if (!a->check_Greek_continuity()) throw 0xbadc0de;
            #endif

            if (a->movability & MOV_BKGRND) continue;

            Point aloc = a->get_barycenter();

            float benerg = 0;
            l = 0;
            for (j=0; j<n; j++)
            {
                if (j==i) continue;
                Molecule* b = mm[j];
                Point bloc = b->get_barycenter();

                Atom* na = a->get_nearest_atom(bloc);
                if (!na) continue;
                float r = na->distance_to(b->get_nearest_atom(aloc));
                if (r > _INTERA_R_CUTOFF) continue;
                nearby[l++] = b;
            }
            nearby[l] = 0;
            benerg = cfmol_multibind(a, nearby);

            #if _dbg_fitness_plummet
            if (!i) cout << "# mol " << a->name << " iter " << iter << ": initial " << -benerg << " ";
            #endif

            #if _dbg_asunder_atoms
            if (!a->check_Greek_continuity()) throw 0xbadc0de;
            #endif

            float tryenerg;
            Pose pib;
            // pib.copy_state(a);

            /**** Linear Motion ****/
            #if allow_linear_motion
            if ((a->movability & MOV_CAN_RECEN) && !(a->movability & MOV_FORBIDDEN))
            {
                Point motion(a->lmx, a->lmy, a->lmz);
                if (motion.magnitude() > speed_limit) motion.scale(speed_limit);

                benerg = cfmol_multibind(a, nearby);
                if (motion.magnitude() > 0.01*speed_limit)
                {
                    pib.copy_state(a);
                    motion.scale(motion.magnitude()/lmsteps);
                    for (l=0; l<lmsteps; l++)
                    {
                        #if _dbg_improvements_only_rule
                        excuse_deterioration = true;
                        #endif

                        a->move(motion);
                        // if (a->agroups.size() && group_realign) group_realign(a, a->agroups);
                        tryenerg = cfmol_multibind(a, nearby);

                        if (tryenerg > benerg && a->contact_maintained())
                        {
                            benerg = tryenerg;
                            pib.copy_state(a);
                        }
                    }
                    pib.restore_state(a);
                    benerg = cfmol_multibind(a, nearby);
                }
                pib.copy_state(a);

                #if _dbg_linear_motion
                if (!a->is_residue()) cout << iter << "! ";
                #endif

                a->lmx = a->lmy = a->lmz = 0;

                int xyz;
                for (xyz=0; xyz<3; xyz++)
                {
                    motion.scale(0);

                    switch(xyz)
                    {
                        case 0: motion.x = frand(-speed_limit,speed_limit); break;
                        case 1: motion.y = frand(-speed_limit,speed_limit); break;
                        case 2: motion.z = frand(-speed_limit,speed_limit); break;
                        default:
                        ;
                    }

                    #if _dbg_improvements_only_rule
                    excuse_deterioration = true;
                    #endif

                    a->move(motion);

                    tryenerg = cfmol_multibind(a, nearby);

                    #if _dbg_fitness_plummet
                    if (!i) cout << "(linear motion try " << -tryenerg << ") ";
                    #endif

                    if (tryenerg > benerg && a->contact_maintained())
                    {
                        benerg = tryenerg;
                        pib.copy_state(a);

                        #if _dbg_linear_motion
                        if (!a->is_residue()) cout << ">";
                        #endif
                    }
                    else
                    {
                        switch(xyz)
                        {
                            case 0: motion.x *= -2; break;
                            case 1: motion.y *= -2; break;
                            case 2: motion.z *= -2; break;
                            default:
                            ;
                        }

                        #if _dbg_improvements_only_rule
                        excuse_deterioration = true;
                        #endif

                        a->move(motion);

                        tryenerg = cfmol_multibind(a, nearby);

                        if (tryenerg > benerg && a->contact_maintained())
                        {
                            benerg = tryenerg;
                            pib.copy_state(a);

                            #if _dbg_linear_motion
                            if (!a->is_residue()) cout << "<";
                            #endif
                        }
                        else
                        {
                            pib.restore_state(a);

                            #if _dbg_linear_motion
                            if (!a->is_residue()) cout << "-" << flush;
                            #endif
                        }
                    }
                }
            }       // If can recenter.
            #if _dbg_linear_motion
            else
            {
                if (!a->is_residue()) cout << iter << "* ";
            }
            #endif
            #endif

            #if _dbg_fitness_plummet
            if (!i) cout << "linear " << -benerg << " ";
            #endif

            #if _dbg_asunder_atoms
            if (!a->check_Greek_continuity()) throw 0xbadc0de;
            #endif

            /**** Histidine flip ****/
            if (mm[i]->hisflips)
            {
                for (l=0; mm[i]->hisflips[l]; l++)
                {
                    #if _DBG_HISFLIP
                    cout << "Flipping " << mm[i]->name << endl;
                    #endif
                    mm[i]->do_histidine_flip(mm[i]->hisflips[l]);

                    tryenerg = cfmol_multibind(a, nearby);

                    if (tryenerg > benerg && a->contact_maintained())
                    {
                        benerg = tryenerg;
                        pib.copy_state(a);
                    }
                    else
                    {
                        mm[i]->do_histidine_flip(mm[i]->hisflips[l]);
                    }
                }
            }
            /**** End histidine flip ****/
            
            #if _dbg_asunder_atoms
            if (!a->check_Greek_continuity()) throw 0xbadc0de;
            #endif

            #if allow_axial_tumble
            if ((a->movability & MOV_CAN_AXIAL) && !(a->movability & MOV_FORBIDDEN))
            {
                pib.copy_state(a);
                Point ptrnd(frand(-1,1), frand(-1,1), frand(-1,1));
                if (frand(0,1) < 0.4 && a->best_intera && a->best_other_intera)
                {
                    ptrnd = a->best_other_intera->get_location().subtract(a->best_intera->get_location());
                }
                if (ptrnd.magnitude())
                {
                    LocatedVector axis = (SCoord)ptrnd;
                    axis.origin = a->get_barycenter(true);
                    float theta;

                    if (a->movability & MOV_MC_AXIAL && frand(0,1) < 0.2) theta = frand(-M_PI, M_PI);
                    else theta = frand(-0.5, 0.5)*fiftyseventh*min(20, iter);

                    #if _dbg_improvements_only_rule
                    excuse_deterioration = true;
                    #endif

                    a->rotate(&axis, theta);
                    tryenerg = cfmol_multibind(a, nearby);

                    if (tryenerg > benerg && a->contact_maintained())
                    {
                        benerg = tryenerg;
                        pib.copy_state(a);
                    }
                    else
                    {
                        pib.restore_state(a);
                    }
                }
            }       // If can axial rotate.
            #endif

            #if _dbg_fitness_plummet
            if (!i) cout << "axial " << -benerg << " ";
            #endif

            #if _dbg_asunder_atoms
            if (!a->check_Greek_continuity()) throw 0xbadc0de;
            #endif

            #if allow_bond_rots
            pib.copy_state(a);
            #if _dbg_mol_flexion
            bool is_flexion_dbg_mol = (a->is_residue() == 107);
            if (is_flexion_dbg_mol) cout << a->name << " movability " << hex << a->movability << dec << endl << flush;
            #endif
            if (((a->movability & MOV_CAN_FLEX) && !(a->movability & MOV_FORBIDDEN)) || a->movability == MOV_FLXDESEL)
            {
                #if _dbg_asunder_atoms
                if (!a->check_Greek_continuity()) throw 0xbadc0de;
                #endif

                float self_clash = max(1.25*a->base_internal_clashes, clash_limit_per_aa);
                Bond** bb = a->get_rotatable_bonds(true);
                if (bb)
                {
                    int q, rang=0;
                    for (q=0; bb[q]; q++)
                    {
                        if (!bb[q]->atom1 || !bb[q]->atom2) continue;         // Sanity check, otherwise we're sure to get random foolish segfaults.
                        if (bb[q]->atom1->get_Greek() > bb[q]->atom2->get_Greek()) bb[q] = bb[q]->get_reversed();
                        if (!bb[q]->count_moves_with_atom2()) continue;
                        if (bb[q]->atom1->is_backbone && strcmp(bb[q]->atom1->name, "CA")) continue;
                        if (bb[q]->atom2->is_backbone) continue;
                        float theta;
                        int heavy_atoms = bb[q]->count_heavy_moves_with_atom2();
                        if (heavy_atoms && (!(a->movability & MOV_CAN_FLEX) || (a->movability & MOV_FORBIDDEN))) continue;

                        #if _dbg_mol_flexion
                        bool is_flexion_dbg_mol_bond = is_flexion_dbg_mol & !strcmp(bb[q]->atom2->name, "OG");
                        #endif

                        if (do_full_rotation && a->is_residue() /*&& benerg <= 0*/ && bb[q]->can_rotate)
                        {
                            float best_theta = 0;
                            Pose prior_state;
                            prior_state.copy_state(a);
                            for (theta=_fullrot_steprad; theta < M_PI*2; theta += _fullrot_steprad)
                            {
                                bb[q]->rotate(_fullrot_steprad, false);
                                if (a->agroups.size() && group_realign)
                                {
                                    // group_realign(a, a->agroups);
                                }
                                tryenerg = cfmol_multibind(a, nearby);

                                #if _dbg_mol_flexion
                                if (is_flexion_dbg_mol_bond) cout << (theta*fiftyseven) << "deg: " << -tryenerg << endl;
                                #endif

                                if (tryenerg > benerg && a->contact_maintained() && a->get_internal_clashes() <= self_clash)
                                {
                                    benerg = tryenerg;
                                    best_theta = theta;
                                }
                            }
                            if (best_theta)
                            {
                                prior_state.restore_state(a);
                                bb[q]->rotate(best_theta, false);
                                a->been_flexed = true;

                                if (a->agroups.size() && group_realign)
                                {
                                    // group_realign(a, a->agroups);
                                }

                                #if _dbg_mol_flexion
                                if (is_flexion_dbg_mol_bond) cout << "Rotating to " << (best_theta*fiftyseven) << "deg." << endl << endl;
                                #endif
                            }
                        }
                        else
                        {
                            if (a->movability & MOV_MC_FLEX && frand(0,1) < 0.25) theta = frand(-_fullrot_steprad, _fullrot_steprad);
                            else if (!heavy_atoms) theta = frand(-_fullrot_steprad, _fullrot_steprad);
                            else if (bb[q]->count_heavy_moves_with_atom() < heavy_atoms)
                                theta = frand(-0.3, 0.3)*fiftyseventh*min(iter, 20);
                            else theta = frand(-0.3, 0.3)*fiftyseventh*min(iter, 20);

                            if (!bb[q]->can_rotate)
                            {
                                bb[q]->compute_flip_capability();
                                if (!bb[q]->flip_angle) bb[q]->flip_angle = M_PI;
                            }

                            Ring* isra = bb[q]->atom1->in_same_ring_as(bb[q]->atom2);
                            if (isra)
                            {
                                if (rang) continue;
                                isra->flip_atom(bb[q]->atom1);
                                rang++;
                                flipped_rings = true;
                            }
                            else bb[q]->rotate(theta, false);

                            #if bb_realign_flexions
                            if (!a->is_residue() && a->agroups.size() && group_realign) group_realign(a, a->agroups);
                            #endif

                            tryenerg = cfmol_multibind(a, nearby);

                            #if _dbg_mol_flexion
                            if (is_flexion_dbg_mol_bond) cout << "Trying " << (theta*fiftyseven) << "deg rotation...";
                            #endif

                            if (tryenerg > benerg && a->contact_maintained() && a->get_internal_clashes() <= self_clash)
                            {
                                benerg = tryenerg;
                                // pib.copy_state(a);
                                a->been_flexed = true;

                                #if _dbg_mol_flexion
                                if (is_flexion_dbg_mol_bond) cout << " energy now " << -tryenerg << ", keeping." << endl << endl;
                                #endif
                            }
                            else
                            {
                                pib.restore_state(a);

                                #if _dbg_mol_flexion
                                if (is_flexion_dbg_mol_bond) cout << " energy from " << -benerg << " to " << -tryenerg << ", reverting." << endl << endl;
                                #endif
                            }
                        }       // full rotation else
                    }       // each rotatable bond

                    int ii;
                    if (a->is_residue()) for (ii=0; mm[ii]; ii++)
                    {
                        if (!mm[ii]->is_residue()) continue;
                        if (mm[ii] == a) continue;
                        float aaclash = mm[ii]->get_intermol_clashes(a);
                        if (aaclash > clash_limit_per_aa)
                        {
                            pib.restore_state(a);
                        }
                    }
                }       // Rotatable bonds.
                #endif

                #if _dbg_fitness_plummet
                if (!i) cout << "flexion " << -benerg << " ";
                #endif

                mm[i]->lastbind = benerg;
            }   // if MOV_CAN_FLEX

            #if _dbg_fitness_plummet
            if (!i) cout << "final " << -benerg << " " << endl << flush;
            #endif

            #if bb_realign_mol
            if (a->agroups.size() && group_realign)
            {
                Pose pre_realign(a);
                benerg = cfmol_multibind(a, nearby);
                group_realign(a, a->agroups);
                tryenerg = cfmol_multibind(a, nearby);

                if (tryenerg < benerg) pre_realign.restore_state(a);
            }
            #endif

            if (!a->is_residue() && flipped_rings) a->evolve_structure(100);

            if (!i && !a->is_residue())
            {
                float ttl_atom_mtn = a->iterbegan->total_atom_motions() / a->get_heavy_atom_count();
                if (ttl_atom_mtn < iter_lostreturns_threshold) a->iters_without_change++;
                else a->iters_without_change = 0;
                if (a->iters_without_change >= max_iters_without_ligand_change) iter = iters;
            }

            if (!(i%8) && progress)
            {
                float f = (float)i / n;
                float fiter = (f + iter) / iters * 100;
                progress(fiter);
            }
        }       // for i

        #if allow_iter_cb
        if (cb) cb(iter+1, mm);
        #endif

        minimum_searching_aniso *= 0.99;

        #if _dbg_fitness_plummet
        if (!i) cout << endl << flush;
        #endif
    }       // for iter

    #if _dbg_linear_motion
    cout << endl;
    #endif

    minimum_searching_aniso = 0;
}

int Molecule::get_heavy_atom_count() const
{
    if (!atoms) return 0;

    int i, result=0;
    for (i=0; i<atcount; i++)
    {
        if (!atoms[i]) break;
        if (atoms[i]->get_Z() > 1) result++;
    }

    return result;
}

#define dbg_optimal_contact 0
SCoord Molecule::motion_to_optimal_contact(Molecule* l)
{
    Pose p;
    p.copy_state(this);
    int vvc = vdw_vertex_count;

    MovabilityType m = movability, lm = l->movability;
    movability = l->movability = MOV_FLEXONLY;

    SCoord total_motion(0,0,0);
    SCoord incremental_motion;

    float energy = -this->get_intermol_binding(l);

    incremental_motion = this->get_barycenter().subtract(l->get_barycenter());
    incremental_motion.r = sgn(energy);

    int i;

    #if dbg_optimal_contact
    cout << "Energy was " << energy << endl;
    #endif

    for (i=0; i<50; i++)
    {
        this->move(incremental_motion, true);
        float new_energy = -this->get_intermol_binding(l);

        if (new_energy < energy)        // Improvement
        {
            energy = new_energy;
            total_motion = total_motion.add(incremental_motion);

            #if dbg_optimal_contact
            cout << "Energy has improved to " << energy << "; keeping." << endl;
            #endif
        }
        else
        {
            incremental_motion.r *= -1;
            this->move(incremental_motion, true);
            incremental_motion.r *= 0.8;

            #if dbg_optimal_contact
            cout << "Energy tried " << energy << endl;
            #endif
        }

        if (fabs(incremental_motion.r) < 0.01) break;
    }

    p.restore_state(this);
    movability = m;
    l->movability = lm;

    vdw_vertex_count = vvc;
    return total_motion;
}

Atom* numbered[10];
bool ring_warned = false;

bool Molecule::from_smiles(char const * smilesstr, bool use_parser)
{
    if (strchr(smilesstr, '{')) use_parser = false;	    // {AtomName} is a nonstandard feature and cannot be handled by a third party app.
    if (use_parser)
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
                while (buffer[0] != '$' && !feof(pf))
                {
                    fgets(buffer, 1022, pf);
                    lno++;
                    sdfdat += (std::string)buffer;

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

    // hydrogenate(true);
    float anomaly = correct_structure();
    cout << "# Structural anomaly = " << anomaly << endl;
    hydrogenate(false);
    identify_conjugations();
    identify_rings();
    identify_cages();

    return retval;
}

bool Molecule::from_smiles(char const * smilesstr, Atom* ipreva)
{
    if (!smilesstr) return true;

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

    smlen = strlen(smilesstr);
    paren = new SMILES_Parenthetical[smlen];

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
                preva->EZ_flip = EZgiven[dbi-1] = sgn(EZ) == sgn(lastEZ);
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

                add_ring(aloop);
                float anomaly = close_loop(aloop, card);

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

    identify_conjugations();
    return true;
}


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
        Bond* ab[16];
        ring_members[i]->fetch_bonds(ab);
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
    Bond* a0b[16];
    ring_members[0]->fetch_bonds(a0b);
    if (!a0b[0]->atom2 || !a0b[1]->atom2)
    {
        cout << "Attempted to form coplanar ring with starting atom bonded to fewer than two other atoms." << endl;
        throw 0xbad12196;					// If you use your imagination, 12196 spells "ring".
    }

    Point A, B, C;
    A = ring_members[0]->get_location();
    B = a0b[0]->atom2->get_location();
    C = a0b[1]->atom2->get_location();

    normal = compute_normal(&A, &B, &C);
    while (!normal.r)
    {
    	B.x = frand(-0.01, 0.01);
    	B.y = frand(-0.01, 0.01);
    	B.z = frand(-0.01, 0.01);
    	C.x = frand(-0.01, 0.01);
    	C.y = frand(-0.01, 0.01);
    	C.z = frand(-0.01, 0.01);
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
            if (!b2->atom2)
            {
            	// cout << "nullptr atom2." << endl;
                b2=0;
                continue;
            }
            for (j=0; j<ringsz; j++)
            {
            	if (ring_members[j] == b2->atom2)
                {
            		// cout << "atom2 is part of the ring." << endl;
                    b2 = 0;
                    break;
                }
            }
            if (b2) break;
        }
        if (b2 && b2->atom2)
        {
        	// cout << "found " << b2->atom2->name << endl;
            Point ptnew = ring_members[l]->get_location().subtract(ringcen);
            ptnew.scale(InteratomicForce::covalent_bond_radius(ring_members[l], b2->atom2, b2->cardinality));
            ptnew = ring_members[l]->get_location().add(ptnew);
            b2->atom2->move_assembly(&ptnew, ring_members[l]);
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
        path[i]->doing_ring_closure = true;
        Bond* b[16];
        path[i]->fetch_bonds(b);
        if (!b[0]) continue;
        int geo = path[i]->get_geometry();

        for (j=0; j<geo; j++)
        {
            if (b[j] && b[j]->atom2
                    &&
                    (	b[j]->can_rotate
                        ||
                        b[j]->can_flip
                    )
                    &&
                    strcmp(b[j]->atom2->name, "N")
               ) rotables[k++] = b[j];
        }

        last = path[i];
    }
    rotables[k] = 0;
    if (_DBGCLSLOOP) cout << "Close Loop: found " << k << " rotables." << endl;

    if (last == first) return 0;
    last->mirror_geo = -1;
    int ringsize = k;

    float bond_length = InteratomicForce::covalent_bond_radius(first, last, lcard);

    if (ringsize < 5)
    {
        // TODO: Make equilateral ring, except accommodating any differences in bond length.
        // But know that cyclopentane and cyclopropane are not coplanar, but rather puckered because of steric hindrance.
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
            float rr = rotables[i]->atom1->get_location().get_3d_distance(rotables[i]->atom2->get_location());

            if (fabs(rr-bond_length) > 0.01)
            {
                Point aloc = rotables[i]->atom1->get_location();
                Point bloc = rotables[i]->atom1->get_location();

                bloc = bloc.subtract(aloc);
                bloc.scale(bond_length);
                bloc = bloc.add(aloc);

                rotables[i]->atom2->move(bloc);
            }

            // issue_5
            if (rotables[i]->rotate(bondrot[i], true))
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
        path[i]->doing_ring_closure = false;
        path[i]->clear_geometry_cache();
        if (_DBGCLSLOOP) cout << "Reset " << path[i]->name << " geometry." << endl;
    }

    for (i=0; rotables[i]; i++)
        rotables[i]->can_rotate = false;

    return anomaly;
}

Atom** Molecule::get_most_bindable(int max_count)
{
    if (noAtoms(atoms)) return 0;
    if (most_bindable) return most_bindable;

    int i, j=-1, k, l;
    float best[max_count+2];
    most_bindable = new Atom*[max_count+2];

    for (k=0; k<max_count; k++)
    {
        best[k]=0;
        most_bindable[k]=0;
    }

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;

        #if _DBG_MOLBB
        cout << "Probing atom " << atoms[i]->name << endl;
        #endif

        float score = 0;
        atoms[i]->clear_geometry_cache();

        for (k=0; most_bindable[k] && k<max_count; k++)
        {
            if (most_bindable[k] == atoms[i])
            {
                #if _DBG_MOLBB
                cout << "Atom is already in return array." << endl;
                #endif
                goto _resume;
            }
            if (most_bindable[k]->is_bonded_to(atoms[i]) && best[k] >= 10)
            {
                #if _DBG_MOLBB
                cout << "Atom is bonded to " << most_bindable[k]->name << ", already in return array." << endl;
                #endif
                goto _resume;
            }
            if (most_bindable[k]->shares_bonded_with(atoms[i]) && best[k] >= 10)
            {
                #if _DBG_MOLBB
                cout << "Atom shares a bond with " << most_bindable[k]->name << ", already in return array." << endl;
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
                    most_bindable[l] = most_bindable[l-1];
                }

                best[k] = score;
                most_bindable[k] = atoms[i];
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
    most_bindable[max_count] = 0;

    if (j < 0) return 0;
    else return most_bindable;
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
        InteratomicForce* iff[32];
        InteratomicForce::fetch_applicable(atoms[i], for_atom, iff);
        if (!iff[0]) continue;

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
    if (!rings[ringid]) return false;
    return rings[ringid]->is_coplanar();
}

bool Molecule::ring_is_aromatic(int ringid) const
{
    if (!rings) return false;
    if (!rings[ringid]) return false;
    return rings[ringid]->get_type() == AROMATIC;
}

Point Molecule::get_ring_center(int ringid)
{
    if (!rings) return Point(0,0,0);
    if (!rings[ringid]) return Point(0,0,0);
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

int Molecule::get_num_rings() const
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

float Molecule::get_atom_error(int i, LocatedVector* best_lv, bool hemi)
{
    int j;
    float error = 0;
    Point bloc;
    Atom* atom2;
    Bond* b[16];
    int g, bg;
    float b_bond_angle;
    LocatedVector lv;

    // Get the atom's zero-index bonded atom. Call it atom2 (because why not overuse a foolish pun?).
    atoms[i]->fetch_bonds(b);
    if (!b[0]) return error;
    atom2 = b[0]->atom2;
    float card = b[0]->cardinality;
    if (!atom2)
    {
        return error;
    }

    bloc = atom2->get_location();
    g = atoms[i]->get_geometry();
    bg = atom2->get_geometry();
    b_bond_angle = atom2->get_geometric_bond_angle();

    // Make an imaginary sphere around atom2, whose radius equals the optimal bond distance.
    lv.origin = bloc;
    if (atoms[i]->get_Z() == 1 || atom2->get_Z() == 1) card = 1;
    float optimal_radius = InteratomicForce::covalent_bond_radius(atoms[i], atom2, card);

    lv = (SCoord)atoms[i]->get_location().subtract(bloc);
    lv.origin = bloc;

    error += _SANOM_BOND_RADIUS_WEIGHT * (fabs(optimal_radius-lv.r)/optimal_radius);

    error += pow(_SANOM_BOND_ANGLE_WEIGHT*atom2->get_bond_angle_anomaly(lv, atoms[i]), 2);

    float thstep = fiftyseventh*5;
    float besttheta = 0, bestphi = 0, bestscore = 0;
    if (hemi)
    {
        bestscore = -1e9;
        lv.r = InteratomicForce::covalent_bond_radius(atoms[i], atom2, card);
        for (lv.theta = -square; lv.theta <= square; lv.theta += thstep)
        {
            float phstep = M_PI/(20.0*(sin(lv.theta) + 1));
            for (lv.phi = 0; lv.phi < (M_PI*2); lv.phi += phstep)
            {
                // At many points along the sphere, evaluate the goodness-of-fit as a function of:
                // Success in conforming to atom2's geometry;
                // Success in avoiding clashes with atoms not bonded to self or atom2;
                // Success in maintaining optimal binding distances to own bonded atoms.
                // Later, we'll test edge cases where bond strain distorts the usual angles.
                float score = 0;

                score -= _SANOM_BOND_ANGLE_WEIGHT*atom2->get_bond_angle_anomaly(lv, atoms[i]);

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
                    if (!b[j]->atom2) continue;
                    float optimal = InteratomicForce::covalent_bond_radius(atoms[i], b[j]->atom2, b[j]->cardinality);
                    float r = b[j]->atom2->get_location().get_3d_distance(lv.to_point());

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
    }

    for (j=0; atoms[j]; j++)
    {
        if (j == i) continue;
        if (atoms[j]->is_bonded_to(atoms[i])) continue;

        float r = atoms[j]->get_location().get_3d_distance(lv.to_point());
        error += _SANOM_CLASHES_WEIGHT/fabs(r+0.000000001);
    }

    for (j=1; b[j]; j++)
    {
        if (!b[j]->atom2) continue;
        float optimal = InteratomicForce::covalent_bond_radius(atoms[i], b[j]->atom2, b[j]->cardinality);
        float r = b[j]->atom2->get_location().get_3d_distance(lv.to_point());

        error += _SANOM_BOND_RAD_WEIGHT * fabs(optimal-r);
    }

    return fmax(0, error+bestscore);
}


#define _DEV_FIX_MSTRUCT 1
float Molecule::correct_structure(int iters)
{
    if (noAtoms(atoms)) return 0;
    int iter, i, j, k, n;
    Point zero(0,0,0);
    float error = 0;
    Point aloc, bloc;
    Atom* atom2;
    Bond* b[16];
    int g, bg;
    float b_bond_angle;
    LocatedVector lv;

    #if _DEV_FIX_MSTRUCT
    // TODO
    if (n = get_num_rings())            // Assignment, not comparison.
    {
        for (i=0; i<n; i++)
        {
            if (ring_is_aromatic(i) || get_ring_num_atoms(i) < 6)
            {
                make_coplanar_ring(get_ring_atoms(i), i);
            }
        }
    }
    else return 0;			// Non-ring structures work fine. The buggy algorithm is when rings are involved.

    for (iter=0; iter<iters; iter++)
    {
        error = 0;
        for (i=0; atoms[i]; i++)
        {
            // Get the atom's zero-index bonded atom. Call it atom2 (because why not overuse a foolish pun?).
            atoms[i]->fetch_bonds(b);
            if (!b[0]) return error;
            atom2 = b[0]->atom2;
            if (!atom2) return error;

            // TODO
            if (atoms[i]->num_rings() && atoms[i]->is_pi() && in_same_ring(atoms[i], atom2)) continue;

            error += get_atom_error(i, &lv);
            if (!iter) continue;

            // Once a "best fit" point in space is found, move there.
            if (atoms[i]->num_rings() && atoms[i]->is_pi())
            {
                for (j=0; j<n; j++)
                {
                    Atom** ring_atoms_j = get_ring_atoms(j);
                    for (k=0; ring_atoms_j[k]; k++)
                    {
                        if (ring_atoms_j[k] == atoms[i])
                        {
                            // Get distance from atom to ring center.
                            Point rcen = get_ring_center(j);
                            Point aloc = atoms[i]->get_location();
                            float rad = rcen.get_3d_distance(aloc);

                            // Center ring at combined distance from atom2.
                            LocatedVector lvr = lv;
                            lvr.r += rad;
                            recenter_ring(j, lvr.to_point());

                            // Rotate ring, and all assemblies, to align atom to atom2.
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
                atoms[i]->move_assembly(&pt, atom2);*/
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

bool Molecule::is_thiol()
{
    if (!atoms) return false;
    int i, j;

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->get_Z() == 16)
        {
            Atom* H = atoms[i]->is_bonded_to("H");
            if (H) return true;
        }
    }

    return false;
}

void Molecule::identify_cages()
{
    if (!rings) return;

    int i, j, l, n;
    for (i=0; rings[i]; i++)
    {
        Atom** ra = rings[i]->get_atoms();
        if (!ra) continue;
        for (j=0; ra[j]; j++)
        {
            for (l=j+1; ra[l]; l++)
            {
                Ring* other = ra[j]->in_same_ring_as(ra[l], rings[i]);

                #if _dbg_identify_rings
                if (other) cout << ra[j]->name << " in same ring as " << ra[l]->name << ": " << *other << endl;
                else cout << ra[j]->name << " -/- " << ra[l]->name << endl;
                #endif

                if (other)
                {
                    // Any bond between the two atoms cannot flip.
                    Bond* b = ra[j]->get_bond_between(ra[l]);
                    if (b)
                    {
                        b->can_rotate = b->can_flip = false;
                        b = ra[l]->get_bond_between(ra[j]);
                        if (b) b->can_rotate = b->can_flip = false;

                        #if _dbg_identify_rings
                        cout << *b << " cannot flip." << endl;
                        #endif
                    }
                    else
                    {
                        // If no bond between atoms, all bonds in both rings cannot flip.
                        Bond** bb = rings[i]->get_bonds();
                        if (bb) for (n=0; bb[n]; n++)
                        {
                            bb[n]->can_rotate = bb[n]->can_flip = false;
                            bb[n]->caged = true;
                        }
                        delete[] bb;
                        bb = other->get_bonds();
                        if (bb) for (n=0; bb[n]; n++)
                        {
                            bb[n]->can_rotate = bb[n]->can_flip = false;
                            bb[n]->caged = true;
                        }
                        delete[] bb;

                        #if _dbg_identify_rings
                        cout << *rings[i] << " immobilized." << endl;
                        cout << *other << " immobilized." << endl;
                        #endif
                    }
                }
            }
        }
        delete[] ra;
    }
}

float Molecule::get_atom_bond_length_anomaly(Atom* a, Atom* ignore)
{
    if (!atoms) return 0;
    if (!a) return 0;

    Bond* thesebonds[16];
    a->fetch_bonds(thesebonds);
    if (!thesebonds) return 0;
    int i;
    float anomaly = 0;

    int geometry = a->get_geometry();
    for (i=0; i<geometry; i++)
    {
        if (!thesebonds[i]) continue;
        if (thesebonds[i]->atom2)
        {
            if (thesebonds[i]->atom2 == ignore) continue;
            float optimal = InteratomicForce::covalent_bond_radius(a, thesebonds[i]->atom2, thesebonds[i]->cardinality);
            float r = a->distance_to(thesebonds[i]->atom2);
            anomaly += pow(1.0+fabs(r-optimal)/optimal, 4)-1;

            #if _dbg_molstruct_evolution_bond_lengths
            if (fabs(r-optimal) > 0.01) cout << "Atoms " << a->name << " and " << thesebonds[i]->atom2->name
                << ", cardinality = " << thesebonds[i]->cardinality
                << " should be " << optimal << "Å apart, but they're " << r << endl;
            #endif
        }
    }

    return anomaly;
}

float Molecule::evolve_structure(int gens, float mr, int ps)
{
    if (!atoms) return 0;
    int ac = get_atom_count();

    Point parents[2][ac];
    Point population[ps][ac];
    float anomalies[ps];
    int i, j, l, n, gen;
    float r, optimal;

    for (i=0; i<ac; i++)
    {
        parents[0][i] = parents[1][i] = atoms[i]->get_location();
    }

    int ibest, i2best;          // Indices of best and second-best anomalies.
    float fbest, f2best;        // Values of best and second-best anomalies.

    // Main loop
    for (gen=1; gen<=gens; gen++)
    {
        ibest = i2best = -1;

        for (i=0; i<ps; i++)
        {
            for (j=0; j<ac; j++)
            {
                float atom_displacement = _evolution_atom_displacement;
                if (atoms[j]->is_pi()) atom_displacement /= _evolution_aromatic_rigidity;

                switch (i)
                {
                    case 0:
                    case 1:
                    population[i][j] = parents[i][j];
                    break;

                    default:
                    int which_parent = rand() & 0x1;
                    population[i][j].x = parents[which_parent][j].x;
                    if (frand(0, 1) < mr) population[i][j].x += frand(-atom_displacement, atom_displacement);
                    population[i][j].y = parents[which_parent][j].y;
                    if (frand(0, 1) < mr) population[i][j].y += frand(-atom_displacement, atom_displacement);
                    population[i][j].z = parents[which_parent][j].z;
                    if (frand(0, 1) < mr) population[i][j].z += frand(-atom_displacement, atom_displacement);

                    if (frand(0,1) < 0.81)
                    {
                        Bond* b0 = atoms[j]->get_bond_by_idx(0);
                        if (b0)
                        {
                            Atom* aparent = b0->atom2;
                            if (aparent)
                            {
                                optimal = InteratomicForce::covalent_bond_radius(atoms[j], aparent, b0->cardinality);
                                SCoord v = population[i][j].subtract(aparent->get_location());
                                v.r = optimal;
                                population[i][j] = aparent->get_location().add(v);
                            }
                        }
                    }
                }

                atoms[j]->move(population[i][j]);
            }

            float anomaly = 0;
            for (j=0; j<ac; j++)
            {
                // anomaly += get_atom_error(j, nullptr, false);
                Bond* b0 = atoms[j]->get_bond_by_idx(0);
                if (!b0) continue;

                Atom* aparent = b0->atom2;
                if (!aparent) continue;

                anomaly += get_atom_bond_length_anomaly(aparent, atoms[j]);

                SCoord v = atoms[j]->get_location().subtract(aparent->get_location());
                anomaly += aparent->get_bond_angle_anomaly(v, atoms[j]);
            }
            anomalies[i] = anomaly;

            if (ibest < 0 || anomaly < fbest)
            {
                i2best = ibest;
                f2best = fbest;
                ibest = i;
                fbest = anomaly;
            }
            else if (i2best < 0 || anomaly < f2best)
            {
                i2best = i;
                f2best = anomaly;
            }
        }

        for (i=0; i<ac; i++)
        {
            parents[0][i] = population[ibest][i];
            parents[1][i] = population[i2best][i];
        }

        #if _dbg_molstruct_evolutions
        cout << ibest << ':' << fbest << ' ' << i2best << ':' << f2best << endl << flush;
        #endif
    }

    for (i=0; i<ac; i++)
    {
        atoms[i]->move(population[ibest][i]);
    }

    return fbest / get_atom_count();
}









