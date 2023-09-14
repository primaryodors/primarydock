
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <ctime>
#include <stdint.h>
#include <cassert>
#include <strings.h>
#include "molecule.h"
#include "aminoacid.h"

using namespace std;

AADef aa_defs[256];
char* override_aminos_dat=0;
float aa_sim_xref[65536];
AminoAcid* aa_archetypes[256];

void AminoAcid::find_his_flips()
{
    if (!atoms) return;
    int i;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        if (atoms[i]->get_family() == TETREL)
        {
            Atom* a = atoms[i]->is_bonded_to(PNICTOGEN, 2);
            if (a && !a->get_charge())
            {
                Atom* b = atoms[i]->is_bonded_to(a->get_elem_sym(), 1);
                if (b && !b->get_charge())
                {
                    Atom* h = b->is_bonded_to("H");
                    if (!h) h = a->is_bonded_to("H");

                    if (h)
                    {
                        hisflips = new HistidineFlip*[4];
                        hisflips[0] = new HistidineFlip;
                        hisflips[0]->C  = atoms[i];
                        hisflips[0]->N1 = a;
                        hisflips[0]->N2 = b;
                        hisflips[0]->H  = h;
                        hisflips[1] = nullptr;
                        break;
                    }
                }
            }
        }
    }
}

AminoAcid::AminoAcid(FILE* instream, AminoAcid* prevaa, int rno)
{
    if (!aa_defs[0x41]._1let) AminoAcid::load_aa_defs();
    immobile = false; // true;
    movability = MOV_FLEXONLY;
    from_pdb(instream, rno);
    base_internal_clashes = get_internal_clashes();
    mol_typ = MOLTYP_AMINOACID;
    prev_aa = prevaa;
    if (prevaa) prevaa->next_aa = this;
    find_his_flips();
}

AminoAcid::~AminoAcid()
{
    if (hisflips) delete[] hisflips;
    if (atoms) delete[] atoms;
    atoms = nullptr;
}

#define _ALGORITHMIC_GREEK 1
AminoAcid::AminoAcid(const char letter, AminoAcid* prevaa, bool minintc)
{
    if (!aa_defs[0x41]._1let) AminoAcid::load_aa_defs();
    immobile = false; // true;
    movability = MOV_FLEXONLY;
    mol_typ = MOLTYP_AMINOACID;

    if (!prevaa) residue_no = 1;
    else residue_no = prevaa->residue_no + 1;

    int idx = letter;
    aadef = &aa_defs[idx];
    /*if (!aa_defs[idx].aabonds)
    {
        cout << "Cannot load " << letter << " please make sure amino acid exists in aminos.dat and all atoms are defined." << endl;
        throw 0xbadac1d;
    }*/

    name = aa_defs[idx].name;
    prev_aa = prevaa;
    if (prevaa) prevaa->next_aa = this;

    int i, j, k, l, n;

    if (aa_defs[idx].SMILES.length())
    {
        std::string fname = (std::string)"sdf/" + (std::string)name + (std::string)".sdf";
        FILE* pf = fopen(fname.c_str(), "rb");
        if (pf)
        {
            fseek(pf, 0, SEEK_END);
            int fsz = ftell(pf);
            fseek(pf, 0, SEEK_SET);
            char buffer[fsz + 4];
            for (i=0; i<(fsz+4); i++) buffer[i] = 0;

            fread(buffer, 1, fsz, pf);
            fclose(pf);

            from_sdf(buffer);
            ensure_pi_atoms_coplanar();
        }
        else
        {
            from_smiles(aa_defs[idx].SMILES.c_str());
            ensure_pi_atoms_coplanar();

            FILE* pf = fopen(fname.c_str(), "wb");
            if (pf)
            {
                save_sdf(pf);
                fclose(pf);
            }
        }

        #ifdef _ALGORITHMIC_GREEK
        int ac4 = get_atom_count()+4;
        int atom_Greek[ac4];
        int atom_append[ac4];
        int atom_prepend[ac4];
        bool atom_isheavy[ac4];
        int atom_prev[ac4];
        int numgrk[29];			// 24 letters, one based, plus extra space.

        Atom *N = nullptr, *HN = nullptr, *CA = nullptr, *C = nullptr, *O = nullptr;

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
                                CA = atoms[i];
                                CA->name = new char[5];
                                strcpy(CA->name, "CA");
                                CA->is_backbone = true;
                                k = atom_idx_from_ptr(CA);
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
                                    if (ab) delete[] ab;
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

                // I hate doing these nonce kludges. It's nonce hence!
                if (!strcmp(aa_defs[idx]._3let, "Trp") && atom_Greek[i] == 7)
                {
                    atom_append[i] = 2;
                    continue;
                }

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
        // cout << aa_defs[idx]._3let << endl;
        for (i=0; atoms[i]; i++)
        {
            /*cout << i
            	 << "\t" << atom_prev[i]
            	 << "\t" << atoms[i]->get_elem_sym()
            	 << "\t" << atoms[i]->name
            	 << "\t" << atom_prepend[i]
            	 << "\t" << atom_Greek[i]
            	 << "\t" << atom_append[i]
            	 ;*/

            if (atom_Greek[i] >= 1)
            {
                std::string newname = "";
                if (atom_prepend[i]) newname += std::to_string(atom_prepend[i]);
                newname += atoms[i]->get_elem_sym();					// TODO: Convert to upper case.
                newname += strGrk.substr(atom_Greek[i]-1, 1);
                if (atom_append[i]) newname += std::to_string(atom_append[i]);

                if (aa_defs[idx].isoleucine_fix && atom_Greek[i]==4 && !atom_append[i]) newname += "1";

                // cout << atoms[i]->name << " is now " << newname << endl;
                strcpy(atoms[i]->name, newname.c_str());
            }

            // cout << "\t" << atoms[i]->name << endl;
        }

        // Sort all the backbone and sidechain atoms into a new array.
        Atom** new_atoms = new Atom*[atcount+4];
        for (i=0; atoms[i]; i++) new_atoms[i] = nullptr;

        l=0;
        if (N) new_atoms[l++] = N;
        if (HN) new_atoms[l++] = HN;

        if (minintc) minimize_internal_clashes();
        ensure_pi_atoms_coplanar();

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

            float r = atoms[i]->distance_to(CA);
            if (r > aa_defs[idx].reach) aa_defs[idx].reach = r;
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

        new_atoms[l] = 0;
        // delete[] atoms;
        atoms = new_atoms;
        // cout << "Atom count was " << atcount << " now " << l << endl;
        atcount = l;

        #endif

        rotatable_bonds = get_rotatable_bonds();
    }	// if SMILES
    else
    {
        if (letter != '#') cout << "WARNING no atoms for " << letter << " (blank SMILES)." << endl;
        return;
    }

    ensure_pi_atoms_coplanar();

    if (!atoms || !atoms[0]) cout << "WARNING no atoms for " << aa_defs[idx].name << " (" << aa_defs[idx].SMILES << ")" << endl;

    if (!aa_defs[idx].aabonds)
    {
        AABondDef** aabd = new AABondDef*[get_atom_count()+16];
        n=0;
        // cout << *this << endl;
        for (i=0; atoms[i]; i++)
        {
            Bond** bb = atoms[i]->get_bonds();
            int bg = atoms[i]->get_geometry();

            if (!i)
            {
                aabd[n] = new AABondDef();
                strcpy(aabd[n]->aname, atoms[i]->name);
                strcpy(aabd[n]->bname, "<C");
                aabd[n]->Za = atoms[i]->get_Z();
                aabd[n]->Zb = 6;
                aabd[n]->cardinality = 1.5;
                aabd[n]->acharge = 0;
                aabd[n]->can_rotate = false;
                n++;
            }

            for (j=0; j<bg; j++)
            {
                if (bb[j]->btom && bb[j]->btom < bb[j]->atom)
                {
                    aabd[n] = new AABondDef();
                    strcpy(aabd[n]->aname, bb[j]->atom->name);
                    strcpy(aabd[n]->bname, bb[j]->btom->name);
                    aabd[n]->Za = bb[j]->atom->get_Z();
                    aabd[n]->Zb = bb[j]->btom->get_Z();
                    aabd[n]->cardinality = bb[j]->cardinality;
                    aabd[n]->acharge = bb[j]->atom->get_charge();

                    if (!strcmp(bb[j]->atom->name, "OH") && !strcmp(bb[j]->btom->name, "CZ"))
                    {
                        aabd[n]->can_rotate = false;
                        aabd[n]->can_flip = true;
                    }
                    
                    aabd[n]->can_rotate =
                        (	aabd[n]->cardinality <= 1.1
                            &&
                            (	!bb[j]->atom->is_pi() || !bb[j]->btom->is_pi()	)
                            &&
                            (	!bb[j]->atom->is_pi()
                                ||
                                (   bb[j]->btom->get_family() != PNICTOGEN
                                    &&
                                    bb[j]->btom->get_family() != CHALCOGEN
                                )
                                ||
                                bb[j]->btom->is_bonded_to_pi(TETREL, false)
                            )
                            &&
                            (	!bb[j]->btom->is_pi()
                                ||
                                (   bb[j]->atom->get_family() != PNICTOGEN
                                    &&
                                    bb[j]->atom->get_family() != CHALCOGEN
                                )
                                ||
                                bb[j]->atom->is_bonded_to_pi(TETREL, false)
                            )
                            &&
                            (	bb[j]->atom->get_family() != PNICTOGEN || !bb[j]->btom->is_pi()	)
                        );
                    aabd[n]->can_flip =
                        (	aabd[n]->cardinality <= 1.1
                            &&
                            (	!bb[j]->atom->is_pi() || !bb[j]->btom->is_pi()	)
                            &&
                            (	(
                                    bb[j]->atom->is_pi()
                                    &&
                                    (   bb[j]->btom->get_family() == PNICTOGEN || bb[j]->btom->get_family() == CHALCOGEN    )
                                    &&
                                    !bb[j]->btom->is_bonded_to_pi(TETREL, false)
                                )
                                ||
                                (
                                    bb[j]->btom->is_pi()
                                    &&
                                    (   bb[j]->atom->get_family() == PNICTOGEN || bb[j]->atom->get_family() == CHALCOGEN    )
                                    &&
                                    !bb[j]->atom->is_bonded_to_pi(TETREL, false)
                                )
                            )
                        );
                    n++;
                }
            }
        }

        aa_defs[idx].aabonds = new AABondDef*[n+2];
        for (i=0; i<n; i++)
        {
            aa_defs[idx].aabonds[i] = aabd[i];
            // cout << aa_defs[idx].name << ":" << aa_defs[idx].aabonds[i]->aname << " is bonded to " << aa_defs[idx].aabonds[i]->bname << "." << endl;
        }
        aa_defs[idx].aabonds[i] = nullptr;

        delete[] aabd;
    }

    ensure_pi_atoms_coplanar();

    // flatten();
    rotatable_bonds = get_rotatable_bonds();
    if (minintc) minimize_internal_clashes();

    ensure_pi_atoms_coplanar();

    identify_acidbase();
    identify_rings();

    if (rings && !aa_defs[idx].aarings)
    {
        int rc2 = get_num_rings() + 2;
        aa_defs[idx].aarings = new Ring*[rc2];
        for (i=0; rings[i]; i++)
        {
            aa_defs[idx].aarings[i] = rings[i];
        }
        aa_defs[idx].aarings[i] = nullptr;
    }

    aa_defs[idx].hydrophilicity = AminoAcid::hydrophilicity();
    if (capable_of_inter(mcoord)) aa_defs[idx].can_coord_metal = true;
    if (fabs(get_charge()) >= 0.9) aa_defs[idx].charged = sgn(get_charge());

    if (rings)
    {
        for (i=0; rings[i]; i++)
        {
            if (rings[i]->is_conjugated() && rings[i]->is_coplanar())
            {
                aa_defs[idx].aromatic = true;
                break;
            }
        }
    }

    if (!isnan(aa_defs[idx].sidechain_pKa))
    {
        float chg = get_charge();
        for (i=0; atoms[i]; i++)
        {
            if (atoms[i]->is_backbone) continue;
            int fam = atoms[i]->get_family();
            if (fam != PNICTOGEN && fam != CHALCOGEN) continue;
            if (atoms[i]->get_charge() || !chg)
                atoms[i]->pK = aa_defs[idx].sidechain_pKa;
        }
    }

    aa_defs[idx].loaded = true;

    Atom* CA = get_atom("CA");
    Atom* CB = get_atom("CB");

    if (prevaa)
    {
        movability = MOV_ALL;

        Point* pts = prevaa->predict_next_NHCA();
        attach_to_prediction(pts);
        delete[] pts;

        movability = MOV_FLEXONLY;
        immobile = true;
    }

    ensure_pi_atoms_coplanar();

    find_his_flips();
}

void AminoAcid::establish_internal_clash_baseline()
{
    base_internal_clashes = 0;
    base_internal_clashes = get_internal_clashes();
}


Point* AminoAcid::predict_previous_COCA()
{
    Atom* N = get_atom("N");
    if (!N) return nullptr;
    N->aromatize();
    Atom* CA = get_atom("CA");
    Atom* HN  = HN_or_substitute();
    if (!CA || !HN) return nullptr;

    Point prevCloc(1.32,0,0), prevOloc(0,1.2,0);

    prevCloc = N->get_location().add(prevCloc);
    prevOloc = prevOloc.add(prevCloc);

    LocatedVector lv;
    Point* retval = new Point[4];
    Point neighborCA;
    int i;
    for (i=0; i<10; i++)
    {
        // Step 1: H-N-C angle.
        SCoord axis = compute_normal(CA->get_location(), N->get_location(), HN->get_location());
        float theta = find_angle_along_vector(HN->get_location(), prevCloc, N->get_location(), axis);
        float plus  = triangular - theta;
        float minus = triangular*2 - theta;

        // Rotate to get 120 degrees not clashing next.N with curr.CA.
        Point maybeP = rotate3D(prevCloc, N->get_location(), axis, plus);
        Point maybeM = rotate3D(prevCloc, N->get_location(), axis, minus);
        lv.copy(axis);
        lv.origin = N->get_location();
        if (maybeP.get_3d_distance(CA->get_location()) > maybeM.get_3d_distance(CA->get_location()))
        {
            prevCloc = maybeP;
            prevOloc = rotate3D(prevOloc, lv.origin, axis, plus);
        }
        else
        {
            prevCloc = maybeM;
            prevOloc = rotate3D(prevOloc, lv.origin, axis, minus);
        }

        // Step 2: N-C-O angle.
        axis = compute_normal(CA->get_location(), N->get_location(), HN->get_location());
        theta = find_angle_along_vector(N->get_location(), prevOloc, prevCloc, axis);
        plus  = triangular - theta;
        minus = triangular*2 - theta;

        // Rotate to get 120 degrees each way.
        maybeP = rotate3D(prevOloc, prevCloc, axis, plus);
        maybeM = rotate3D(prevOloc, prevCloc, axis, minus);

        // Whichever point is farthest from O is the direction of HN (trans configuration).
        lv.copy(axis);
        lv.origin = N->get_location();
        if (maybeP.get_3d_distance(HN->get_location()) > maybeM.get_3d_distance(HN->get_location()))
        {
            prevOloc = maybeP;
        }
        else
        {
            prevOloc = maybeM;
        }

        // Step 3: Enforce trans configuration 180 degrees.
        axis = prevCloc.subtract(N->get_location());
        theta = find_angle_along_vector(prevOloc, HN->get_location(), prevCloc, axis);

        // Rotate the around next.N to get 180 degrees.
        lv.copy(axis);
        lv.origin = prevCloc;
        prevOloc = rotate3D(prevOloc, lv.origin, axis, M_PI+theta);

        // Rotate around C axis C-O to get 180 degrees about that axis.
        axis = HN->get_location().subtract(N->get_location());
        theta = find_angle_along_vector(CA->get_location(), prevCloc, N->get_location(), axis);
        lv.copy(axis);
        lv.origin = N->get_location();
        prevCloc = rotate3D(prevCloc, lv.origin, axis, M_PI-theta);
        prevOloc = rotate3D(prevOloc, lv.origin, axis, M_PI-theta);

        // Step 4: Obtain other residue's CA location.
        Point pt = CA->get_location();
        /*Rotation rot = align_points_3d(&prevOloc, &pt, &prevCloc);
        neighborCA = rotate3D(prevOloc, prevCloc, rot.v, -rot.a);*/
        SCoord normal = compute_normal(prevCloc, prevOloc, N->get_location());
        neighborCA = rotate3D(prevOloc, prevCloc, normal, triangular);
        if (neighborCA.get_3d_distance(N->get_location()) < 1) neighborCA = rotate3D(prevOloc, prevCloc, normal, -triangular);
        neighborCA = neighborCA.subtract(prevCloc);
        neighborCA.scale(1.54);			// TODO: Make not hard-coded.
        neighborCA = neighborCA.add(prevCloc);
    }

    retval[0] = prevCloc;
    retval[1] = prevOloc;
    retval[2] = neighborCA;

    return retval;
}

Point* AminoAcid::predict_next_NHCA()
{
    Atom* C = get_atom("C");
    if (!C) return nullptr;
    C->aromatize();
    Atom* CA = get_atom("CA");
    Atom* O  = get_atom("O");
    if (!CA || !O) return nullptr;

    Point nextNloc(1.32,0,0), nextHNloc(0,1.0,0);

    nextNloc = C->get_location().add(nextNloc);
    nextHNloc = nextHNloc.add(nextNloc);

    LocatedVector lv;
    Point* retval = new Point[4];
    Point neighborCA;
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

        // Step 4: Obtain other residue's CA location.
        Point pt = CA->get_location();
        /*Rotation rot = align_points_3d(&nextHNloc, &pt, &nextNloc);
        neighborCA = rotate3D(nextHNloc, nextNloc, rot.v, -rot.a);*/
        SCoord normal = compute_normal(nextNloc, nextHNloc, C->get_location());
        neighborCA = rotate3D(nextHNloc, nextNloc, normal, triangular);
        if (neighborCA.get_3d_distance(C->get_location()) < 1) neighborCA = rotate3D(nextHNloc, nextNloc, normal, -triangular);
        neighborCA = neighborCA.subtract(nextNloc);
        neighborCA.scale(1.54);			// TODO: Make not hard-coded.
        neighborCA = neighborCA.add(nextNloc);
    }

    retval[0] = nextNloc;
    retval[1] = nextHNloc;
    retval[2] = neighborCA;

    return retval;
}

#define _dbg_attprdc 0
void AminoAcid::attach_to_prediction(Point* predicted, bool CO)
{
    MovabilityType fmov = movability;
    movability = MOV_ALL;
    float anomaly;

    // Translation: move the entire AA so that the N (or C) corresponds to the predicted origin.
    Point moveby = predicted[0].subtract( CO ? get_atom_location("C") : get_atom_location("N") );
    aamove(moveby);
    anomaly = predicted[0].get_3d_distance( CO ? get_atom_location("C") : get_atom_location("N") );
    #if _dbg_attprdc
    if (anomaly > 0.001) cout << "Error: " << ( CO ? "C" : "N" ) << " anomaly outside tolerance!" << endl << "# Anomaly is " << anomaly << endl;
    #endif

    // Rotation: rotate the entire AA about the predicted origin so that the HN (or O) aligns with the predicted vector.
    Point pt1 = CO ? get_atom_location("O") : HN_or_substitute_location(), pt2 = predicted[1];
    Rotation rot = align_points_3d( &pt1, &pt2, &predicted[0] );
    LocatedVector lv = rot.v;
    lv.origin = predicted[0];
    rotate(lv, rot.a);
    pt1 = CO ? get_atom_location("O") : HN_or_substitute_location();
    anomaly = predicted[1].get_3d_distance(pt1);
    #if _dbg_attprdc
    if (anomaly > 0.1) cout << "Error: " << ( CO ? "O" : "HN" ) << " anomaly outside tolerance!" << endl << "# Anomaly is " << anomaly << endl;
    #endif

    // Rotation: rotate the entire AA about its NH or CO axis to bring the CA in line with the predicted CA.
    lv.origin = CO ? get_atom_location("C") : get_atom_location("N");
    SCoord axis = pt1.subtract(lv.origin);
    float theta = find_angle_along_vector(get_atom_location("CA"), predicted[2], lv.origin, axis);
    lv.copy(axis);
    rotate(lv, -theta);

    // Rotation: rotate the rest of the AA to bring the C-N-CA or N-C-CA angle to 120deg.
    pt2 = get_atom_location("CA");
    rot = align_points_3d( &pt2, &predicted[2], &predicted[0] );
    int i;
    Atom* HN = HN_or_substitute();
    for (i=0; atoms[i]; i++)
    {
        if ( CO && !strcmp(atoms[i]->name, "C")) continue;
        if ( CO && !strcmp(atoms[i]->name, "O")) continue;
        if (!CO && !strcmp(atoms[i]->name, "N")) continue;
        if (!CO && atoms[i] == HN) continue;

        Point pt = rotate3D(atoms[i]->get_location(), lv.origin, rot.v, rot.a);
        atoms[i]->move(pt);
    }

    anomaly = predicted[2].get_3d_distance(get_atom_location("CA"));
    if (anomaly > 0.1)
    {
        rotate(lv, theta*2);
        anomaly = predicted[2].get_3d_distance(get_atom_location("CA"));
    }
    #if _dbg_attprdc
    if (anomaly > 0.1) cout << "Error: CA anomaly outside tolerance!" << endl << "# Anomaly is " << anomaly << "." << endl;
    #endif

    movability = fmov;
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

std::string AminoAcid::printable()
{
    stringstream s;
    s << *this;
    return s.str();
}


void AminoAcid::save_pdb(FILE* os, int atomno_offset)
{
    int i;

    for (i=0; atoms[i]; i++)
    {
        atoms[i]->pdbchain = pdbchain;
        atoms[i]->save_pdb_line(os, i+1+atomno_offset);
    }
}

int AminoAcid::from_pdb(FILE* is, int rno)
{
    /*
              1111111111222222222233333333334444444444555555555566666666667777777777
    01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ATOM     55  SG  CYS     4       6.721  -8.103   4.542  1.00001.00           S
	ATOM   1091  N   LYS A1000     -15.894 -46.862 -74.510  1.00104.21           N  
    */
    char buffer[1024], origbuf[1024], res3let[5];
    int added=0, lasttell=0;
    residue_no=0;
    res3let[0] = 0;

    while (!feof(is))
    {
        lasttell = ftell(is);
        fgets(buffer, 1003, is);
        if (buffer[0] == 'A' &&
            buffer[1] == 'T' &&
            buffer[2] == 'O' &&
            buffer[3] == 'M'
           )
            buffer[16] = ' ';

        int thistell = ftell(is);
        strcpy(origbuf, buffer);

        std::vector<std::string> words;
		int places[20] = {0, 6, 11, 17, 21, 22, 30, 38, 46, 54, 60, 76};
		int i, j, k;
		for (i=11; i>=0; i--)
		{
			// words[i] = new char[35];
			j = places[i];

			// Ltrim.
			while (buffer[j] == ' ' && buffer[j+1] > 0) j++;

			k = j+1;
			while (buffer[k]) k++;			// Find zero.
			// Rtrim.
			for (k--; buffer[k] == ' '; k--) buffer[k] = 0;
			
			// strcpy(words[i], &buffer[j]);
            for (k=words.size(); k<=j; k++) words.push_back("");
            words[i] = &buffer[j];
			buffer[places[i]] = 0;
		}

        try
        {
            if (words.size())
            {
                int offset = 1;
                // cout << words[0] << endl;
				if (!strcmp(words[0].c_str(), "ANISOU")) continue;
                if (!strcmp(words[0].c_str(), "ATOM"))
                {
                    if (!residue_no)
                    {
                        residue_no = atoi(words[4+offset].c_str()) + rno;
                    }

                    if (!res3let[0])
                    {
                        strcpy(res3let, words[3].c_str());
                        res3let[1] |= 0x20;
                        res3let[2] |= 0x20;
                    }

                    if (!atno_offset) atno_offset = atoi(words[1].c_str());

                    if (strcasecmp(res3let, words[3].c_str())
                            ||
                            residue_no != (atoi(words[4+offset].c_str())+rno)
                       )
                    {
                        fseek(is, lasttell, SEEK_SET);
                        goto _return_added;
                    }

                    char esym[7] = {0,0,0,0,0,0,0};
                    if (words[2].c_str()[0] >= '0' && words[2].c_str()[0] <= '9')
                        strcpy(esym, &words[2].c_str()[1]);
                    else
                        strcpy(esym, words[2].c_str());

                    if (!strcmp(words[2].c_str(), "OXT")) strcpy(esym, "O");

                    int i;
                    for (i=1; i<6; i++)
                    {
                        if (!esym[i+1]) esym[i] = 0;
                        if (esym[i+1] >= '0' && esym[i+1] <= '9') esym[i]=0;
                        if (i>1) esym[i] = 0;
                        if (!esym[i]) break;
                    }
                    esym[1] &= 0x5f;

                    Point aloc(atof(words[5+offset].c_str()), atof(words[6+offset].c_str()),atof(words[7+offset].c_str()));

                    Atom* a = add_atom(esym, words[2].c_str(), &aloc, 0, 0);
                    a->pdbidx = atoi(words[1].c_str());
                    added++;

                    if (   !strcmp(a->name, "N")
                            || !strcmp(a->name, "HN")
                            || !strcmp(a->name, "CA")
                            || !strcmp(a->name, "C")
                            || !strcmp(a->name, "O")
                       )
                        a->is_backbone = true;
                    else a->is_backbone = false;

                    a->residue = atoi(words[4+offset].c_str())+rno;
                    strcpy(a->aa3let, words[3].c_str());
                    AADef* aaa=0;

                    name=0;
                    for (i=0; i<256; i++)
                    {
                        if (aa_defs[i]._1let && !strcasecmp(aa_defs[i]._3let, words[3].c_str()))
                        {
                            a->aaletter = aa_defs[i]._1let;
                            aaa = &aa_defs[i];
                            name = new char[10];
                            sprintf(name, "%s%d", aa_defs[i]._3let, atoi(words[4+offset].c_str())+rno);
                            break;
                        }
                    }

                    if (!aaa)
                    {
                        fseek(is, lasttell, SEEK_SET);
                        goto _return_added;
                    }
                    else aadef = aaa;

                    if (aaa && !aaa->aabonds)
                    {
                        AminoAcid* tempaa = new AminoAcid(aaa->_1let, 0, false);
                        delete tempaa;
                    }

                    if (aaa->_1let == 'T')
                    {
                        // Have to hard code this because the side chain's chirality forces the numbering system to differ from both
                        // AlphaFold and ZhangLab, even though they agree.
                        if (!strcmp(a->name, "CG2" )) strcpy(a->name, "CG1" );
                        if (!strcmp(a->name, "OG1" )) strcpy(a->name, "OG2" );
                        if (!strcmp(a->name, "HG1" )) strcpy(a->name, "HG2" );
                        if (!strcmp(a->name, "1HG2")) strcpy(a->name, "1HG1");
                        if (!strcmp(a->name, "2HG2")) strcpy(a->name, "2HG1");
                        if (!strcmp(a->name, "3HG2")) strcpy(a->name, "3HG1");
                    }

                    if (aaa->_1let == 'I')
                    {
                        if (!strcmp(a->name, "1HD1")) strcpy(a->name, "HD1");
                        if (!strcmp(a->name, "2HD1")) strcpy(a->name, "HD2");
                        if (!strcmp(a->name, "3HD1")) strcpy(a->name, "HD3");
                    }

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
                                        // if (b) b->can_flip = aab->can_flip;
                                        b = btom->get_bond_between(a);
                                        if (b) b->can_rotate = false;
                                        // if (b) b->can_flip = aab->can_flip;
                                    }
                                }
                            }
                            else if (!strcmp(aab->bname, a->name))
                            {
                                Atom* btom = get_atom(aab->aname);
                                if (btom && !btom->is_bonded_to(a))
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
                else
                {
                    fseek(is, lasttell, SEEK_SET);
                    goto _return_added;
                }
            }
            else goto _return_added;
        }
        catch (int ex)
        {
            if (ex == ATOM_NOT_OF_AMINO_ACID) throw ex;
            if (ex == NOT_ATOM_RECORD) throw ex;
        }
        buffer[0] = 0;
    }

_return_added:
    if (aadef && aadef->aabonds)
    {
        int i, j;

        for (i=0; aadef->aabonds[i]; i++)
        {
            Atom *a = get_atom(aadef->aabonds[i]->aname),
                  *b = get_atom(aadef->aabonds[i]->bname);

            if (a && b && !a->is_bonded_to(b))
            {
                a->bond_to(b, aadef->aabonds[i]->cardinality);
            }
        }
    }

    if (aadef && aadef->aarings)
    {
        int i, j;
        for (i=0; aadef->aarings[i]; i++);	// Get count.
        rings = new Ring*[i+2];

        for (i=0; aadef->aarings[i]; i++)
        {
            rings[i] = nullptr;
            Atom** ringa = aadef->aarings[i]->get_atoms();
            if (!ringa) continue;

            Atom** lra = new Atom*[aadef->aarings[i]->get_atom_count() + 2];

            for (j=0; ringa[j]; j++)
            {
                Atom* la = get_atom(ringa[j]->name);
                if (la)
                {
                    lra[j] = la;
                }
                else
                {
                    lra[j] = new Atom(ringa[j]->get_elem_sym());
                    lra[j]->name = new char[8];
                    strcpy(lra[j]->name, ringa[j]->name);
                }
            }
            lra[j] = nullptr;

            delete[] ringa;
            rings[i] = new Ring(lra, aadef->aarings[i]->get_type());
        }
        rings[i] = nullptr;
    }
    // identify_rings();

    return added;
}

void AminoAcid::copy_loaded_to_object(char letter, int tbdctr, AABondDef** tmpbdefs, bool proline_like)
{
    int lidx = (letter & 0x5f) - 'A';
    if (lidx<0 || lidx>256) return;
    aa_defs[lidx].aabonds = nullptr;
    int j;
    for (j=0; j<tbdctr; j++)
    {
        aa_defs[lidx].proline_like = proline_like;
    }
}



void AminoAcid::load_aa_defs()
{
    FILE* pf = fopen(override_aminos_dat ?: "data/aminos.dat", "rb");
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
        char* lastwords[16];
        lastwords[0] = 0;
        AABondDef** tmpbdefs=0;
        int tbdctr=0;
        char lastletter = '\0';
        bool isbb = false;
        bool proline_like = false;

        for (i=0; i<65536; i++) aa_sim_xref[i] = -1;
        for (i=0; i<256; i++) aa_archetypes[i] = nullptr;

        while (!feof(pf))
        {
            fgets(buffer, 1011, pf);
            if (buffer[0] != '#' && buffer[0] != '\n')
            {
                char** words = chop_spaced_words(buffer);
                if (!words) continue;

                try
                {
                    for (i=0; words[i]; i++) if (words[i][0] == '%')
                        {
                            words[i] = lastwords[i];
                        }

                    int idx = words[0][0];
                    if (idx == 'X')
                    {
                        cout << "Cannot use X as a letter for amino acids. X means accept any residue in a motif search." << endl;
                        throw 0xbadaadef;
                    }

                    if (!lastletter || words[0][0] != lastletter)
                    {
                        copy_loaded_to_object(lastletter, tbdctr, tmpbdefs, proline_like);
                        tbdctr = 0;
                        proline_like = false;

                        if (tmpbdefs) delete[] tmpbdefs;
                    }

                    if (!aa_defs[idx]._1let)
                    {
                        aa_defs[idx]._1let = idx;
                        // cout << aa_defs[idx]._1let;
                        strcpy(aa_defs[idx]._3let, words[1]);
                        strcpy(aa_defs[idx].name, words[2]);
                        aa_defs[idx].reach = 1.09;

                        aa_defs[idx]._3let[1] |= 0x20;
                        aa_defs[idx]._3let[2] |= 0x20;

                        if (!strcmp(aa_defs[idx].name, "isoleucine")) aa_defs[idx].isoleucine_fix = true;
                    }

                    if (words[3])
                    {
                        aa_defs[idx].SMILES = words[3];
                        if (strstr(aa_defs[idx].SMILES.c_str(), "N1")) aa_defs[idx].proline_like = proline_like = true;

                        if (words[4])
                        {
                            aa_defs[idx].sidechain_pKa = atof(words[4]);
                            if (words[5])
                            {
                                aa_defs[idx].flexion_probability = atof(words[5]);
                                if (aa_defs[idx].flexion_probability < 0 || aa_defs[idx].flexion_probability > 1) throw 0xbadf1ec5;
                            }
                        }
                    }

                    tbdctr++;

                    lastletter = words[0][0];
                }
                catch (int e)
                {
                    cout << "Error while reading aminos.dat, please verify file is in the correct format." << endl;
                    throw 0xbadaadef;
                }

                if (!lastwords[0])
                    for (i=0; i<16; i++)
                        lastwords[i] = new char[256];

                for (i=0; words[i]; i++)
                {
                    if (lastwords[i] != words[i]) strcpy(lastwords[i], words[i]);
                }
            }
            buffer[0] = 0;
        }
        copy_loaded_to_object(lastletter, tbdctr, tmpbdefs, proline_like);
        fclose(pf);
    }
}

bool AminoAcid::can_reach(Atom* other) const
{
    Atom* ca1;
    float r;

    ca1 = get_atom("CA");

    if (!ca1 || !other)
    {
        cout << "Warning: Could not determine reach of " << *this << " - " << *other << " to avoid possibility of clash." << endl;
        return false;
    }

    r = ca1->get_location().get_3d_distance(other->get_location());

    if (r <= 1.15 * (get_reach())) return true;
    else return false;
}

bool AminoAcid::can_reach(Molecule* other) const
{
    Atom* ca1, *ca2;
    float r;

    ca1 = get_atom("CA");
    ca2 = other->get_nearest_atom(ca1->get_location());

    if (!ca1 || !ca2)
    {
        cout << "Warning: Could not determine reach of " << *this << " - " << other->get_name() << " to avoid possibility of clash." << endl;
        return false;
    }

    r = ca1->get_location().get_3d_distance(ca2->get_location());

    if (r <= 1.15 * (get_reach())) return true;
    else return false;
}

bool AminoAcid::can_reach(AminoAcid* other) const
{
    Atom *ca1, *cb1, *ca2;
    float r;

    ca1 = get_atom("CA");
    cb1 = get_atom("CB");
    ca2 = other->get_atom("CA");

    if (!ca1 || !ca2)
    {
        cout << "Warning: Could not determine reach of " << *this << " - " << *other << " to avoid possibility of clash." << endl;
        return false;
    }

    r = ca1->get_location().get_3d_distance(ca2->get_location());

    float f = 1;
    if (cb1)
    {
        f = find_3d_angle(cb1->get_location(), ca2->get_location(), ca1->get_location());
        f = 0.5 + 0.5 * cos(f);
    }

    if (r <= f*(get_reach() + other->get_reach() + 4)) return true;
    else return false;
}

Atom* AminoAcid::get_reach_atom()
{
    if (!atoms) return nullptr;

    int i;
    float maxr = 0;
    Atom* CA = get_atom("CA");
    Atom* retval = nullptr;
    if (!CA) return nullptr;

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        float r = atoms[i]->distance_to(CA);

        if (r > maxr)
        {
            maxr = r;
            retval = atoms[i];
        }
    }

    return retval;
}

float AminoAcid::similarity_to(const char letter)
{
    AminoAcid* a = nullptr;
    if (!aa_defs[letter].SMILES.length()) return 0;

    int i = (int)letter + 256*(int)this->get_letter();
    if (aa_sim_xref[i] >= 0) return aa_sim_xref[i];
    
    a = new AminoAcid(letter);
    float s = similarity_to(a);
    delete a;

    aa_sim_xref[i] = s;
    i = (int)this->get_letter() + 256*(int)letter;
    aa_sim_xref[i] = s;

    return s;
}

float AminoAcid::similarity_to(const AminoAcid* aa)
{
    if (!atoms || !aa->atoms) return 0;

    bool polar1 = (hydrophilicity() > 0.2);
    bool polar2 = (aa->hydrophilicity() > 0.2);

    int i, j;

    i = (int)aa->get_letter() + 256*(int)this->get_letter();
    if (aa_sim_xref[i] >= 0) return aa_sim_xref[i];

    float simil=0, divis=0;
    for (i=0; atoms[i]; i++)
    {
        bool apol = fabs(atoms[i]->is_polar()) > hydrophilicity_cutoff;
        if (apol != polar1) continue;

        for (j=0; aa->atoms[j]; j++)
        {
            bool bpol = fabs(aa->atoms[j]->is_polar()) > hydrophilicity_cutoff;
            if (bpol != polar2) continue;

            float f = (apol && bpol) ? 1 : hydrophilicity_cutoff;
            simil += f*atoms[i]->similarity_to(aa->atoms[j]);
            divis += f;
        }
    }

    if (divis) simil /= divis;

    simil -= 0.5*fabs( fmin(1, fabs(hydrophilicity())) - fmin(1, fabs(aa->hydrophilicity())) );
    simil += 0.5*(sgn(get_charge() * aa->get_charge()));
    simil = fmax(0, fmin(1, simil));

    i = (int)aa->get_letter() + 256*(int)this->get_letter();
    aa_sim_xref[i] = simil;
    i = (int)this->get_letter() + 256*(int)aa->get_letter();
    aa_sim_xref[i] = simil;

    return simil;
}

Ring* AminoAcid::get_most_distal_arom_ring()
{
    if (!rings) return nullptr;

    Point caloc = get_atom_location("CA");
    int i, retidx=-1;
    float r = 0;
    for (i=0; rings[i]; i++)
    {
        if (rings[i]->get_type() == AROMATIC)
        {
            float lr = caloc.get_3d_distance(rings[i]->get_center());
            if ((lr > r) || !i)
            {
                r = lr;
                retidx = i;
            }
        }
    }

    if (retidx < 0) return nullptr;
    return rings[retidx];
}

char AminoAcid::set_pdb_chain(char c)
{
    if (c >= 'A' && c <= 'Z')
    {
        pdbchain = c;
    }

    int i;
    if (atoms) for (i=0; atoms[i]; i++)
    {
        atoms[i]->pdbchain = pdbchain;
    }

    return pdbchain;
}

std::ostream& operator<<(std::ostream& os, const AABondDef& b)
{
    os << b.aname;
    os << cardinality_printable(b.cardinality);
    os << b.bname;

    return os;
}

#define DBG_TYRLIKE 0
bool AminoAcid::is_tyrosine_like()
{
    if (!atoms)
    {
        #if DBG_TYRLIKE
        cout << "No atoms." << endl;
        #endif
        return false;
    }
    if (!rings)
    {
        #if DBG_TYRLIKE
        cout << "No rings." << endl;
        #endif
        return false;
    }

    int i, j;
    bool has_aromatic_ring = false;

    for (i=0; rings[i]; i++)
    {
        if (rings[i]->get_type() == AROMATIC) has_aromatic_ring = true;
        #if DBG_TYRLIKE
        cout << "Ring " << i << " type " << rings[i]->get_type() << endl;
        #endif
    }

    if (!has_aromatic_ring)
    {
        #if DBG_TYRLIKE
        cout << "No aromatic rings." << endl;
        #endif
        return false;
    }

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        if (atoms[i]->get_Z() == 1) continue;
        if (atoms[i]->is_polar() < -hydrophilicity_cutoff)
        {
            bool part_of_arom_ring = false;

            for (j=0; rings[j]; j++)
            {
                bool atom_is_in_ring = atoms[i]->is_in_ring(rings[j]);
                if (!atom_is_in_ring) continue;

                bool ring_is_aromatic = rings[j]->get_type() == AROMATIC;
                if (!ring_is_aromatic) continue;

                #if DBG_TYRLIKE
                if (atom_is_in_ring) cout << atoms[i]->name << " is part of ring " << *rings[j] << endl;
                if (ring_is_aromatic) cout << *rings[j] << " is aromatic." << endl;
                #endif

                /*if (atom_is_in_ring && ring_is_aromatic)*/ part_of_arom_ring = true;
            }
            if (!part_of_arom_ring)
            {
                #if DBG_TYRLIKE
                cout << atoms[i]->name << " seems to be ringless, polarity " << atoms[i]->is_polar()
                     << " yup not missing any important info, derp." << endl;
                #endif

                return true;
            }
        }
    }

    return false;
}

bool AminoAcid::is_glycine()
{
    if (!atoms) return false;
    int i;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        if (atoms[i]->get_Z() > 1) return false;
    }
    return true;
}

float AminoAcid::sc_pKa() const
{
    return aadef->sidechain_pKa;
}

bool AminoAcid::conditionally_basic() const
{
    #if _allow_conditional_basicity
    if (!atoms) return false;
    if (!aadef || !aadef->sidechain_pKa) return false;
    int i;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        if (atoms[i]->get_Z() == 1) continue;
        if (atoms[i]->get_family() == PNICTOGEN)
        {
            if (aadef->sidechain_pKa >= 5 && aadef->sidechain_pKa < 7) return true;
        }
    }
    #endif
    return false;
}

std::ostream& operator<<(std::ostream& os, const AminoAcid& aa)
{
    if (!&aa) return os;

    char c = aa.get_pdb_chain();
    if (c && c != ' ') os << c << ":";

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


Atom* AminoAcid::capable_of_inter(intera_type inter)
{
    if (!atoms) return nullptr;

    // Lysine and arginine cannot coordinate metals. Histidine is not normally protonated so it is free for metal coordination.
    // We see real world instances where a mutation of H to R in a metal binding site renders the entire protein nonfunctional.
    if (inter == mcoord && get_charge() > 0) return nullptr;

    Atom* retval = nullptr;
    float cadist = 0;
    Point caloc = get_atom("CA")->get_location();

    int i, j;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        if (InteratomicForce::atom_is_capable_of(atoms[i], inter)) // return atoms[i];
        {
            float r = atoms[i]->get_location().get_3d_distance(caloc);
            if (r > cadist)
            {
                retval = atoms[i];
                cadist = r;
            }
        }
        /*InteratomicForce** iff = InteratomicForce::get_applicable(atoms[i], atoms[i]);
        if (!iff) continue;
        for (j=0; iff[j]; j++)
            if (iff[j]->get_type() == inter)
            {
                // cout << *this << " is capable of " << inter << " binding because of atom " << atoms[i]->name << endl;
                delete[] iff;
                return true;
            }*/
    }

    return retval;
}

float AminoAcid::get_phi()
{
    if (!prev_aa)
    {
        if (!next_aa) return 0;
        return next_aa->get_phi();
    }

    Atom *C0, *N, *CA, *C1;
    C0 = prev_aa->get_atom("C");
    N  = get_atom("N");
    CA = get_atom("CA");
    C1 = get_atom("C");

    if (!C0 || !N || !CA || !C1) return 0;

    SCoord axis = CA->get_location().subtract(N->get_location());

    return find_angle_along_vector(C0->get_location(), C1->get_location(), CA->get_location(), axis);
}

float AminoAcid::get_psi()
{
    if (!next_aa)
    {
        if (!prev_aa) return 0;
        return prev_aa->get_phi();
    }

    Atom *N0, *CA, *C, *N1;
    N0 = get_atom("N");
    CA = get_atom("CA");
    C  = get_atom("C");
    N1 = next_aa->get_atom("N");

    if (!N0 || !CA || !C || !N1) return 0;

    SCoord axis = C->get_location().subtract(CA->get_location());

    return find_angle_along_vector(N0->get_location(), N1->get_location(), CA->get_location(), axis);
}

float AminoAcid::get_omega()
{
    if (!next_aa)
    {
        if (!prev_aa) return 0;
        return prev_aa->get_phi();
    }

    Atom *CA0, *C, *N, *CA1;
    CA0 = get_atom("CA");
    C   = get_atom("C");
    N   = next_aa->get_atom("N");
    CA1 = next_aa->get_atom("CA");

    if (!CA0 || !C || !N || !CA1) return 0;

    SCoord axis = N->get_location().subtract(C->get_location());

    return find_angle_along_vector(CA0->get_location(), CA1->get_location(), N->get_location(), axis);
}

bond_rotation_fail_reason AminoAcid::rotate_phi(float a)
{
    if (aadef->proline_like) return bf_disallowed_rotation;
    Atom* N  = get_atom("N");
    Atom* CA = get_atom("CA");
    Bond* b = CA->get_bond_between(N);
    if (!b) return bf_bond_not_found;
    b->rotate(a, true, true);
    return b->last_fail;
}

bond_rotation_fail_reason AminoAcid::rotate_psi(float a)
{
    Atom* CA = get_atom("CA");
    Atom* C  = get_atom("C");
    Bond* b = CA->get_bond_between(C);
    if (!b) return bf_bond_not_found;
    b->rotate(a, true, true);
    return b->last_fail;
}

bool AminoAcid::is_alpha_helix()
{
    return is_helix(4);
}

bool AminoAcid::is_helix(int p)
{
    int i, j;
    Atom *a, *b;
    AminoAcid* aa;

    if (!prev_aa && !next_aa) return false;

    float f = (next_aa) ? (get_psi() + next_aa->get_phi()) : (prev_aa->get_psi() + get_psi());
    f *= fiftyseven;
    if (f > 180) f -= 360;

    // cout << residue_no << " has f = " << f << endl;
    float pf = -195.0 + 360.0 / p;
    if (fabs(f - pf) > 18) return false;

    for (j=0; j<2; j++)
    {
        aa = this;
        for (i=0; i<p; i++)
        {
            aa = j ? aa->next_aa : aa->prev_aa;
            if (!aa) break;
        }

        if (aa)
        {
            a = (j ? this : aa)->get_atom("O");
            b = (j ? aa : this)->get_atom("HN");
            if (!b) b = (j ? aa : this)->get_atom("H");

            if (a && b)
            {
                float r = a->distance_to(b);
                if (r < helix_hbond_cutoff) return true;
            }
        }
    }

    return false;
}

float AminoAcid::hydrophilicity() const
{
    int i, count=0;
    float total=0, weight;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        if (!strcmp(atoms[i]->name, "HA")) continue;

        weight = 1;

        int Z = atoms[i]->get_Z();
        int fam = atoms[i]->get_family();
        // if (Z==1) continue;
        if (Z==1) fam = atoms[i]->get_bond_by_idx(0)->btom->get_family();
        if (Z==1 && fam == TETREL) continue;

        if (Z == 7) weight = 2.5;
        else if (Z == 8) weight = 3;
        else if (Z == 15) weight = 1.5;
        else if (Z == 16) weight = 1.25;

        if (fam == PNICTOGEN && conditionally_basic()) total += protonation(sc_pKa())*2;

        float h = atoms[i]->hydrophilicity_rule();
        total += h*weight;
        count++;
    }
    return count ? (total / count) : 0;
}

void AminoAcid::hydrogenate(bool steric_only)
{
    if (!atoms) return;

    Atom* oxt = get_atom("OXT");
    if (oxt)
    {
        Atom* C = get_atom("C");
        if (!oxt->get_bond_between(C)) oxt->bond_to(C, 1);
    }    
    
    int i, j, k, l, n;

    #if hydrogenate_add_missing_heavy_atoms
    if (aadef && get_atom("CB"))
    {
        for (i=0; aadef->aabonds[i]; i++)
        {
            if (aadef->aabonds[i]->Za < 2 || aadef->aabonds[i]->Zb < 2) continue;
            Atom* a = get_atom(aadef->aabonds[i]->aname);
            Atom* b = get_atom(aadef->aabonds[i]->bname);

            if (a && a->is_backbone) continue;
            if (b && b->is_backbone) continue;

            if (!aa_archetypes[aadef->_1let]) aa_archetypes[aadef->_1let] = new AminoAcid(aadef->_1let);
            AminoAcid* at = aa_archetypes[aadef->_1let];
            if (!at->get_atom("CB")) continue;

            at->movability = MOV_ALL;
            at->aamove(get_CA_location().subtract(at->get_CA_location()));
            Rotation rot = align_points_3d(at->get_atom("CB")->get_location(), get_atom("CB")->get_location(), get_CA_location());
            LocatedVector lv;
            lv = rot.v;
            lv.origin = get_CA_location();
            at->rotate(lv, rot.a);
            Atom* c;

            if (!a && !b)
            {
                cout << "Warning: amino acid definition " << aadef->name
                    << " not compatible with hydrogenate_add_missing_heavy_atoms; atoms "
                    << aadef->aabonds[i]->aname << " and " << aadef->aabonds[i]->bname
                    << " are out of sequence." << endl;
            }
            else if (!a)
            {
                a = add_atom( Atom::esym_from_Z(aadef->aabonds[i]->Za), aadef->aabonds[i]->aname, b, aadef->aabonds[i]->cardinality );
                a->residue = b->residue;
                strcpy(a->aa3let, b->aa3let);
                a->increment_charge(aadef->aabonds[i]->acharge);
                c = at->get_atom(aadef->aabonds[i]->aname);
                if (c) a->move(c->get_location());
                added_heavies = true;
            }
            else if (!b)
            {
                b = add_atom( Atom::esym_from_Z(aadef->aabonds[i]->Zb), aadef->aabonds[i]->bname, a, aadef->aabonds[i]->cardinality );
                b->residue = a->residue;
                strcpy(b->aa3let, a->aa3let);
                c = at->get_atom(aadef->aabonds[i]->bname);
                if (c) b->move(c->get_location());
                added_heavies = true;
            }
            else if (!a->is_bonded_to(b))
            {
                a->bond_to(b, aadef->aabonds[i]->cardinality);
            }
        }
    }

    if (added_heavies)
    {
        minimize_internal_clashes();
    }
    #endif


    Molecule::hydrogenate(steric_only);
    int already[128][4];
    Atom* onlyone[128][4];
    const char* alpha = "ABGDEZH";

    for (i=0; i<24; i++) for (j=0; j<4; j++)
    {
        already[i][j] = 0;
        onlyone[i][j] = nullptr;
    }

    Bond** bt;
    Bond* bb;
    Atom* heavy;

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->get_Z() > 1) continue;
        bt = atoms[i]->get_bonds();
        if (!bt) continue;
        bb = bt[0];
        if (!bb) continue;
        heavy = bb->btom;
        if (!heavy) continue;

        atoms[i]->residue = heavy->residue;
        atoms[i]->aaletter = heavy->aaletter;
        strcpy(atoms[i]->aa3let, heavy->aa3let);
        atoms[i]->is_backbone = heavy->is_backbone;

        if (heavy->is_backbone && heavy->get_family() == PNICTOGEN)
            strcpy(atoms[i]->name, "HN");
        else
        {
            char greek = heavy->name[1];
            if (greek >= 'a') greek = heavy->name[2];
            const char* gramma = strchr(alpha, greek);
            if (gramma)
            {
                j = gramma - alpha;
                n = 0;
                for (k=1; heavy->name[k]; k++) if (heavy->name[k] >= '0' && heavy->name[k] <= '9')
                {
                    n = atoi(&heavy->name[k]);
                    break;
                }

                k = ++(already[j][n]);
                if (residue_no == 309)
                {
                    onlyone[j][n]++;
                }
                if (k == 1) onlyone[j][n] = atoms[i];
                else onlyone[j][n] = nullptr;
                char aname[10];
                if (n) sprintf(aname, "%dH%c%d", k, greek, n);
                else   sprintf(aname, "H%c%d", greek, k);
                aname[4] = 0;       // Molecule class imposes a 4-char limit.

                strcpy(atoms[i]->name, aname);
            }
        }
    }

    for (i=0; i<7; i++) for (j=0; j<4; j++)
    {
        if (onlyone[i][j])
        {
            if (onlyone[i][j]->name[0] == '1')
            {
                char tmpbuf[10];
                strcpy(tmpbuf, &onlyone[i][j]->name[1]);
                strcpy(onlyone[i][j]->name, tmpbuf);
            }
            else
            {
                for (k=1; onlyone[i][j]->name[k] && onlyone[i][j]->name[k] != '1'; k++);
                onlyone[i][j]->name[k] = 0;
            }
        }
    }

    Atom* atomtmp[get_atom_count()+8];
    l = 0;
    Atom* cursor = get_atom("N");
    if (cursor)
    {
        atomtmp[l++] = cursor;
        for (i=0; atoms[i]; i++)
            if (atoms[i]->get_Z() == 1)
            {
                bt = atoms[i]->get_bonds();
                if (bt)
                {
                    bb = bt[0];
                    if (bb && bb->btom == cursor)
                        atomtmp[l++] = atoms[i];
                }
            }
    }

    for (j=0; alpha[j]; j++)
    {
        for (i=0; atoms[i]; i++)
        {
            if (atoms[i]->get_Z() > 1)
            {
                if (atoms[i]->is_backbone && strcmp(atoms[i]->name, "CA")) continue;

                char greek = atoms[i]->name[1];
                if (greek >= 'a') greek = atoms[i]->name[2];
                const char* gramma = strchr(alpha, greek);
                if (gramma && *gramma == alpha[j])
                {
                    cursor = atoms[i];
                    if (cursor)
                    {
                        atomtmp[l++] = cursor;
                        for (k=0; atoms[k]; k++)
                            if (atoms[k]->get_Z() == 1)
                            {
                                bt = atoms[k]->get_bonds();
                                if (bt)
                                {
                                    bb = bt[0];
                                    if (bb && bb->btom == cursor)
                                        atomtmp[l++] = atoms[k];
                                }
                            }
                    }
                }
            }
        }
    }

    cursor = get_atom("C");
    if (cursor) atomtmp[l++] = cursor;

    cursor = get_atom("O");
    if (cursor) atomtmp[l++] = cursor;

    cursor = get_atom("OXT");
    if (cursor)
    {
        atomtmp[l++] = cursor;
        for (i=0; atoms[i]; i++)
            if (atoms[i]->get_Z() == 1)
            {
                bt = atoms[i]->get_bonds();
                if (bt)
                {
                    bb = bt[0];
                    if (bb && bb->btom == cursor)
                        atomtmp[l++] = atoms[i];
                }
            }
    }

    atomtmp[l] = 0;

    for (i=0; i<l; i++) atoms[i] = atomtmp[i];
    atoms[i] = nullptr;

    ensure_pi_atoms_coplanar();
}

Point AminoAcid::get_CA_location()
{
    Atom* a = get_atom("CA");
    if (!a) return Point(0,0,0);
    return a->get_location();
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

void AminoAcid::rotate(LocatedVector vec, float theta)
{
    if (!atoms) return;
    // cout << name << " AminoAcid::rotate()" << endl;

    ensure_pi_atoms_coplanar();

    int i;
    for (i=0; atoms[i]; i++)
    {
        Point loc = atoms[i]->get_location();
        Point nl  = rotate3D(&loc, &vec.origin, &vec, theta);
        atoms[i]->move(&nl);
        // cout << *this << ":" << atoms[i]->name << nl << " " << endl;
    }
    // cout << endl;

    ensure_pi_atoms_coplanar();

    // If you have a metal coordination, AND YOU ARE THE FIRST COORDINATING RESIDUE OF THE METAL, move the metal with you.
    if (m_mcoord && m_mcoord->coord_res && m_mcoord->coord_res[0] == this)
    {
        Point loc = m_mcoord->metal->get_location();
        Point nl  = rotate3D(&loc, &vec.origin, &vec, theta);
        m_mcoord->metal->move(&nl);
    }

    ensure_pi_atoms_coplanar();
}

void AminoAcid::renumber(int new_resno)
{
	residue_no = new_resno;
	if (!atoms) return;

	int i;
	for (i=0; atoms[i]; i++)
	{
		atoms[i]->residue = residue_no;
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

Atom* AminoAcid::HN_or_substitute()
{
    Atom* retval = get_atom("HN");
    if (!retval) retval = get_atom("H");
    if (!retval)
    {
        Atom* a = get_atom("N");
        if (!a) return nullptr;
        int i;
        Bond** bb = a->get_bonds();
        if (!bb) return nullptr;
        int g = a->get_geometry();

        for (i=0; i<g; i++)
        {
            if (bb[i]->btom && strcmp(bb[i]->btom->name, "CA")) return bb[i]->btom;
        }
    }
    return retval;
}

Point AminoAcid::HN_or_substitute_location()
{
    Atom* a = HN_or_substitute();
    if (!a) return Point(0,0,0);
    else return a->get_location();
}

LocRotation AminoAcid::enforce_peptide_bond(bool cis)
{
    LocRotation retval;
    if (!atoms) return retval;
    if (!prev_aa) return retval;

    Atom *pC = prev_aa->get_atom("C");
    Atom *pO = prev_aa->get_atom("O");
    Atom* lN = get_atom("N");
    Atom* lH = HN_or_substitute();

    if (!pC || !pO || !lN || !lH) return retval;

    Point ptC = pC->get_location();
    Point ptO = pO->get_location();
    Point ptN = lN->get_location();
    Point ptH = lH->get_location();

    SCoord v = ptN.subtract(ptC);

    // TODO: Refine the algorithm further so that coplanarity is enforced, including CA atoms, and bond angles as well.
    float theta, step=1.0*fiftyseventh, bestr = 0, besttheta=0;
    if (cis) bestr = 999999;
    for (theta=0; theta < M_PI*2; theta+=step)
    {
        Point ptnew = rotate3D(ptH, ptN, v, theta);
        float r = ptnew.get_3d_distance(ptO);
        if ((cis && (r < bestr)) || (!cis && (r > bestr)))
        {
            bestr = r;
            besttheta = theta;
        }
    }

    retval.origin = ptN;
    retval.v = v;
    retval.a = besttheta;

    int i;
    for (i=0; atoms[i]; i++)
    {
        Point ptnew = rotate3D(atoms[i]->get_location(), retval.origin, retval.v, retval.a);
        atoms[i]->move(ptnew);
        // cout << atoms[i]->residue << ":" << atoms[i]->name << "->" << ptnew << " ";
    }
    // cout << endl;

    return retval;
}

LocRotation* AminoAcid::flatten()
{
    Bond* b;
    LocRotation* retval = new LocRotation[5];
    if (m_mcoord) return retval;	// NO.

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
    Atom* localHN = HN_or_substitute();
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

        /*if (aadef) cout << aadef->_3let;
        cout << residue_no;
        cout << ":" << endl;*/
        float planar = 9999, r = 9999;
        for (i=0; i<250; i++)
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

            Point pCA, pC, pO, lN, lCA, lHN, lC, lO;

            pCA = prevCA->get_location();
            pC = prevC->get_location();
            pO = prevO->get_location();

            if (j < 3)
            {
                lN = rotate3D(localN->get_location(), retval[j].origin, retval[j].v, retval[j].a+ad[j]);
                lCA = rotate3D(localCA->get_location(), retval[j].origin, retval[j].v, retval[j].a+ad[j]);
                lHN = rotate3D(localHN->get_location(), retval[j].origin, retval[j].v, retval[j].a+ad[j]);
                lC = rotate3D(localC->get_location(), retval[j].origin, retval[j].v, retval[j].a+ad[j]);
                lO = rotate3D(localO->get_location(), retval[j].origin, retval[j].v, retval[j].a+ad[j]);
            }
            else
            {
                lN = localN->get_location();
                lCA = localCA->get_location();
                lHN = localHN->get_location();
                lC = localC->get_location();
                lO = localO->get_location();
            }

            float lplanar, lr;

            switch (j)
            {
            case 0:
                lplanar = are_points_planar(pC, pO, lN, lHN);
                lr = 4;
                break;

            case 1:
                lplanar = are_points_planar(pC, lN, lHN, lCA);
                lr = 4;
                break;

            case 2:
                lplanar = are_points_planar(pC, pO, lN, pCA);
                lr = 4;
                break;

            case 3:
                lplanar = are_points_planar(lN, lCA, lHN, lC);
                lr = lHN.get_3d_distance(lO);
                break;

            case 4:
                lplanar = are_points_planar(lCA, lHN, lC, lO);
                lr = lHN.get_3d_distance(lO);
                break;

            default:
                ;
            }

            if ( (i < 60 && j >= 3) ? (lr <= r) : (lplanar <= planar) )
            {
                retval[j].a += ad[j];
                if (retval[j].a > M_PI) retval[j].a -= M_PI*2;
                if (fabs(ad[j]) < 0.5) ad[j] *= 1.01;
                planar = lplanar;
                // cout << planar << " ";
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
                ad[j] *= -0.75;
                // cout << "x ";
            }
        }
        // cout << endl << endl;

        LocatedVector lv = retval[j].get_lv();
        switch(j)
        {
        case 0:
        case 1:
        case 2:
            if (retval[j].a > M_PI) retval[j].a -= M_PI*2;
            rotate(lv, retval[j].a);
            break;

        default:
            ;
        }
    }

    return retval;
}

void AminoAcid::ensure_pi_atoms_coplanar()
{
    if (!atoms) return;

    int i, j, n;
    bool dirty[get_atom_count()];
    Point conjugated[get_atom_count()];
    std::string conj_aname[get_atom_count()];
    for (i=0; atoms[i]; i++) dirty[i] = false;

    for (i=0; atoms[i]; i++)
    {
        if (dirty[i]) continue;
        if (!atoms[i]->is_pi()) continue;

        n = 1;
        conj_aname[0] = atoms[i]->name;
        conjugated[0] = atoms[i]->get_location();
        for (j = i+1; atoms[j]; j++)
        {
            if (!atoms[j]->is_pi()) continue;
            if (!atoms[i]->is_conjugated_to(atoms[j])) continue;

            conj_aname[n] = atoms[j]->name;
            conjugated[n++] = atoms[j]->get_location();
            dirty[j] = true;
        }

        if (n < 4) continue;

        for (j = 3; j < n; j++)
        {
            float result = are_points_planar(conjugated[0], conjugated[1], conjugated[2], conjugated[j]);
            if (result >= coplanar_threshold)
            {
                cout << *this << " has non-coplanar pi atoms "
                    << conj_aname[0] << ", " << conj_aname[1] << ", " << conj_aname[2] << ", " << conj_aname[j]
                    << " having anomaly " << result
                    << endl;
            }
        }
    }
}

LocRotation AminoAcid::rotate_backbone_abs(bb_rot_dir dir, float angle)
{
    // For the N-CA bond, of a residue not at the N-terminus, there's one torsion angle that will place
    // the local C atom as far as possible from the previous residue's C atom.
    // For the CA-C bond, of a residue not at the C-terminus, there is a torsion angle that will place
    // the local N atom as far as possible from the next residue's N atom.
    LocRotation retval;
    if (m_mcoord)
    {
        cout << *this << " is a metal-coordinating residue; cannot rotate backbone." << endl;
        return retval;	// NO.
    }
    // return retval;

    if (aadef && aadef->proline_like) return retval;

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
        cout << "Bad direction for " << *this << "; cannot rotate backbone." << endl;
        return retval;
    }
    if (!atom || !btom)
    {
        cout << "Not found axial atom for " << *this << "; cannot rotate backbone." << endl;
        return retval;
    }
    Bond* b = atom->get_bond_between(btom);
    if (!b)
    {
        atom->bond_to(btom, 1);
        b = atom->get_bond_between(btom);
    }
    if (!b)
    {
        cout << "No bond between " << *this << ":" << atom->name << "-" << btom->name << "; cannot rotate backbone." << endl;
        return retval;
    }
    if (!m_mcoord) b->can_rotate = true;

    b->clear_moves_with_cache();

    switch (dir)
    {
    case N_asc:
    case CA_desc:
        atom = HN_or_substitute();
        btom = get_atom("C");
        break;

    case CA_asc:
    case C_desc:
        atom = get_atom("N");
        btom = get_atom("O");
        break;

    default:
        cout << "Bad direction (second attempt) for " << *this << "; cannot rotate backbone." << endl;
        return retval;
    }

    // Use the SCoord as the axis and make an imaginary circle, finding the angle that brings HN and O closest.
    int i;
    float theta, step = 3, oldstep = 360;
    float bestrad=0, bestr;

    if (!atom || !btom)
    {
        cout << "Not found reference atom for " << *this << "; cannot rotate backbone." << endl;
        return retval;
    }

    while (step > 0.01)
    {
        bestr = atom->get_location().get_3d_distance(btom->get_location()) * 1000;
        for (theta=bestrad-oldstep; theta<bestrad+oldstep; theta+=step)
        {
            b->rotate(fiftyseventh*step, true);
            float r = atom->get_location().get_3d_distance(btom->get_location());
            #if 0
            /*if (dir == N_asc)*/ cout << b->atom->name << "-" << b->btom->name << " rotation " << theta << ": "
                << atom->name << "-" << btom->name << " distance = " << r << endl << flush;
            #endif
            if (r < bestr)
            {
                bestrad = fiftyseventh*(theta+step);
                bestr = r;
            }
        }

        if (oldstep < 359.9999) b->rotate(fiftyseventh * -(oldstep+step), true);
        oldstep = step;
        step /= 10;
    }
    // cout << "Best rotation = " << (bestrad*fiftyseven) << " degrees." << endl;

    // To this maximum stretch angle, add the input angle and do the backbone rotation.
    // cout << "Calling " << *this << ".rotate_backbone( " << (bestrad*fiftyseven) << " + " << (angle*fiftyseven) << ")..." << endl;
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
    if (!atoms)
    {
        // cout << "No atoms for " << get_3letter() << residue_no << "; cannot rotate backbone." << endl;
        return retval;
    }
    if (m_mcoord)
    {
        // cout << get_3letter() << residue_no << "has a metal coordination; cannot rotate backbone." << endl;
        return retval;	// NO.
    }
    // return retval;

    if (aadef && aadef->proline_like) return retval;

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
    if (!rotcen)
    {
        // cout << get_3letter() << residue_no << ":" << aname << "not found; cannot rotate backbone." << endl;
        return retval;
    }
    retval.origin = rotcen->get_location();
    btom = get_atom(bname);
    if (!btom)
    {
        // cout << get_3letter() << residue_no << ":" << bname << "not found; cannot rotate backbone." << endl;
        return retval;
    }

    Bond* b = rotcen->get_bond_between(btom);
    if (!b) // || !b->can_rotate)
    {
        // cout << "Non-rotatable bond " << *this << ":" << rotcen->name << "-" << btom->name << "; cannot rotate backbone." << endl;
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
        btom = HN_or_substitute();
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

float AminoAcid::get_intermol_binding(AminoAcid* neighb, bool backbone_atoms_only)
{
    AminoAcid* neighbs[4];
    int i;
    for (i=0; i<4; i++) neighbs[i] = nullptr;
    neighbs[0] = neighb;
    return get_intermol_binding(neighbs, backbone_atoms_only);
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















