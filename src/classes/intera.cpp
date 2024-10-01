
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <ctime>
#include <stdint.h>
#include "intera.h"

using namespace std;

float total_binding_by_type[_INTER_TYPES_LIMIT];
float minimum_searching_aniso = 0;
InteratomicForce* lif = nullptr;
SCoord missed_connection(0,0,0);
float mc_bpotential = 0;

#if _peratom_audit
std::vector<std::string> interaudit;
bool interauditing = false;
#endif

Interaction::Interaction()
{
    attractive = 0;
    repulsive = 0;
}

Interaction::Interaction(double v)
{
    attractive = fmax(v, 0);
    repulsive = fmax(-v, 0);
}

Interaction Interaction::operator+(Interaction const& obj)
{
    Interaction result;
    result.attractive = attractive + obj.attractive;
    result.repulsive = repulsive + obj.repulsive;
    return result;
}
    
Interaction Interaction::operator+=(Interaction const& obj)
{
    attractive += obj.attractive;
    repulsive += obj.repulsive;
    return obj;
}

Interaction Interaction::operator*(double const& f)
{
    Interaction result;
    result.attractive = attractive * f;
    result.repulsive = repulsive * f;
    return result;
}

Interaction Interaction::operator+=(float const& f)
{
    attractive += f;
    return *this;
}

Interaction Interaction::operator-=(float const& f)
{
    repulsive += f;
    return *this;
}

void InteratomicForce::append_by_Z(int Za, int Zb, InteratomicForce* iff)
{
    int i;

    if (!forces_by_Z[Za][Zb])
    {
        forces_by_Z[Za][Zb] = new InteratomicForce*[24];
        for (i=0; i<24; i++)
            forces_by_Z[Za][Zb][i] = 0;
    }

    for (i=0; i<23; i++)
    {
        if (!forces_by_Z[Za][Zb][i])
        {
            forces_by_Z[Za][Zb][i] = iff;
            // if (all_forces[ifcount]->type == mcoord) cout << all_forces[ifcount]->Za << " " << all_forces[ifcount]->Zb << " " << i << " " << all_forces[ifcount]->type << endl;
            break;
        }
    }
}

void InteratomicForce::read_all_forces()
{
    int i, ifcount = 0;
    init_nulls(all_forces, _MAX_NUM_FORCES);

    FILE* pf = fopen("data/bindings.dat", "rb");
    if (!pf)
    {
        cout << "ERROR failed to open bindings.dat, please verify file exists and you have permissions." << endl;
        throw 0xbadf12e;
    }
    else
    {
        reading_forces = true;
        char buffer[1024];
        while (!feof(pf))
        {
            fgets(buffer, 1011, pf);
            if (buffer[0] != '#' && buffer[0] != '\n')
            {
                all_forces[ifcount] = new InteratomicForce();
                all_forces[ifcount]->read_dat_line(buffer);
                /*if (all_forces[ifcount]->type == mcoord)
                	cout << all_forces[ifcount]->Za << "*M*" << all_forces[ifcount]->Zb << endl;*/

                if (all_forces[ifcount]->Za >= 1 && all_forces[ifcount]->Za < 36 
                    &&
                    all_forces[ifcount]->Zb >= 1 && all_forces[ifcount]->Zb < 36
                )
                {
                    append_by_Z(all_forces[ifcount]->Za, all_forces[ifcount]->Zb, all_forces[ifcount]);
                    append_by_Z(all_forces[ifcount]->Zb, all_forces[ifcount]->Za, all_forces[ifcount]);
                }

                if (all_forces[ifcount]->Za == any_element
                    &&
                    all_forces[ifcount]->Zb >= 1 && all_forces[ifcount]->Zb < 36
                )
                {
                    for (i=1; i<36; i++)
                    {
                        append_by_Z(i, all_forces[ifcount]->Zb, all_forces[ifcount]);
                        append_by_Z(all_forces[ifcount]->Zb, i, all_forces[ifcount]);
                    }
                }

                if (all_forces[ifcount]->Zb == any_element
                    &&
                    all_forces[ifcount]->Za >= 1 && all_forces[ifcount]->Za < 36
                )
                {
                    for (i=1; i<36; i++)
                    {
                        append_by_Z(i, all_forces[ifcount]->Za, all_forces[ifcount]);
                        append_by_Z(all_forces[ifcount]->Za, i, all_forces[ifcount]);
                    }
                }

                if (all_forces[ifcount]->Za == any_element && all_forces[ifcount]->Zb == any_element) throw 0x666b75;

                ifcount++;
            }
            buffer[0] = 0;
        }
        fclose(pf);
        reading_forces = false;

        all_forces[ifcount] = 0;
        read_forces_dat = true;
    }
}

InteratomicForce::InteratomicForce()
{
    if (!read_forces_dat && !reading_forces) read_all_forces();
}

void InteratomicForce::read_dat_line(char* line)
{
    char** words = chop_spaced_words(line);
    if (!words) return;
    if (words[0]
            && words[1]
            && words[2]
            && words[3]
            && words[4]
       )
    {
        char ea[3], eab[3], eb[3], ebb[3];

        if (words[0][1] >= 'A' && words[0][1] <= 'Z')
        {
            eab[0] = words[0][0];
            eab[1] = 0;
            strcpy(ea, words[0]+1);
        }
        else if (words[0][1] > 0 && words[0][1] < 'A')
        {
            if (words[0][1] == '-') aritybZa = 1;
            if (words[0][1] == '=') aritybZa = 2;
            if (words[0][1] == '#') aritybZa = 3;
            if (words[0][1] == '$') aritybZa = 4;
            eab[0] = words[0][0];
            eab[1] = 0;
            strcpy(ea, words[0]+2);
        }
        else if (words[0][2] >= 'A' && words[0][2] <= 'Z')
        {
            eab[0] = words[0][0];
            eab[1] = words[0][1];
            eab[2] = 0;
            strcpy(ea, words[0]+2);
        }
        else if (words[0][2] > 0 && words[0][2] < 'A')
        {
            if (words[0][2] == '-') aritybZa = 1;
            if (words[0][2] == '=') aritybZa = 2;
            if (words[0][2] == '#') aritybZa = 3;
            if (words[0][2] == '$') aritybZa = 4;
            eab[0] = words[0][0];
            eab[1] = words[0][1];
            eab[2] = 0;
            strcpy(ea, words[0]+3);
        }
        else
        {
            eab[0] = 0;
            strcpy(ea, words[0]);
        }


        if (words[1][1] >= 'A' && words[1][1] <= 'Z')
        {
            ebb[0] = words[1][0];
            ebb[1] = 0;
            strcpy(eb, words[1]+1);
        }
        else if (words[1][1] > 0 && words[1][1] < 'A')
        {
            if (words[1][1] == '-') aritybZb = 1;
            if (words[1][1] == '=') aritybZb = 2;
            if (words[1][1] == '#') aritybZb = 3;
            if (words[1][1] == '$') aritybZb = 4;
            ebb[0] = words[1][0];
            ebb[1] = 0;
            strcpy(eb, words[1]+2);
        }
        else if (words[1][2] >= 'A' && words[1][2] <= 'Z')
        {
            ebb[0] = words[1][0];
            ebb[1] = words[1][1];
            ebb[2] = 0;
            strcpy(eb, words[1]+2);
        }
        else if (words[1][2] > 0 && words[1][2] < 'A')
        {
            if (words[1][2] == '-') aritybZb = 1;
            if (words[1][2] == '=') aritybZb = 2;
            if (words[1][2] == '#') aritybZb = 3;
            if (words[1][2] == '$') aritybZb = 4;
            ebb[0] = words[1][0];
            ebb[1] = words[1][1];
            ebb[2] = 0;
            strcpy(eb, words[1]+3);
        }
        else
        {
            ebb[0] = 0;
            strcpy(eb, words[1]);
        }


        Za  = Atom::Z_from_esym(ea);
        bZa = Atom::Z_from_esym(eab);
        Zb  = Atom::Z_from_esym(eb);
        bZb = Atom::Z_from_esym(ebb);			// There is only ebb, no flow. Cope.


        if (!strcmp(words[2], "coval"))	type = covalent;
        if (!strcmp(words[2], "ionic"))	type = ionic;
        if (!strcmp(words[2], "hbond"))	type = hbond;
        if (!strcmp(words[2], "pi"   ))	type = pi;
        if (!strcmp(words[2], "plpi" ))	type = polarpi;
        if (!strcmp(words[2], "coord"))	type = mcoord;
        if (!strcmp(words[2], "vdW"  ))	type = vdW;

        if (words[3])
        {
            arity = atof(words[3]);
            if (words[4])
            {
                distance = atof(words[4]);
                if (words[5])
                {
                    kJ_mol = atof(words[5]);
                    if (words[6])
                    {
                        dirprop = atof(words[6]);
                    }
                    else dirprop = 2;
                }
                else kJ_mol = 200;
            }
            else distance=1;
        }
    }
    delete[] words;
}

bool InteratomicForce::atom_is_capable_of(Atom* a, intera_type t)
{
    InteratomicForce** look = all_forces;
    int i, Za = a->get_Z();

    for (i=0; look[i]; i++)
    {
        if (look[i]->type == t)
            if (look[i]->Za == Za || look[i]->Zb == Za || look[i]->Za == any_element || look[i]->Zb == any_element)
            {
                switch (t)
                {
                case covalent:
                    if (a->get_valence()) return true;
                    break;

                case ionic:
                    if (a->get_charge()
                            ||
                            a->get_acidbase()
                       )
                        return true;
                    break;

                case hbond:
                    if (a->is_polar()) return true;
                    break;

                case pi:
                    if (a->is_pi()) return true;
                    break;

                case polarpi:
                    if (a->is_pi() || a->is_metal())
                        return true;
                    break;

                case mcoord:
                    if (a->get_charge() <= 0)		// -NH3+ groups would be unable to coordinate metals because of the steric hindrance of the extra proton.
                        return true;
                    break;

                case vdW:
                    return true;
                    break;

                default:
                    ;
                }
            }
    }
    return false;
}

#define _dbg_applicable 0
void InteratomicForce::fetch_applicable(Atom* a, Atom* b, InteratomicForce** retval)
{
    if (!read_forces_dat && !reading_forces) read_all_forces();
    if (!a || !b)
    {
        retval[0] = nullptr;
        return;
    }

    #if _dbg_applicable
    cout << "Getting forces between " << a->name << " and " << b->name << "..." << endl;
    #endif

    InteratomicForce** look = all_forces;
    int Za = a->get_Z();
    int Zb = b->get_Z();

    if (Za < 36 && Zb < 36)
    {
        if (forces_by_Z[Za][Zb])
        {
            look = forces_by_Z[Za][Zb];
        }
    }

    retval[0] = nullptr;
    int i, j=0;

    // Charged atoms always attract or repel, irrespective of Z. But charges in conjugated systems move to interact at the mutual nearest atoms.
    bool do_ionic = (a->get_charge() && b->get_charge());
    if (do_ionic)
    {
        if (a->conjugation && b->conjugation)
        {
            if (a != a->conjugation->get_nearest_atom(b->conjugation)) do_ionic = false;
            if (b != b->conjugation->get_nearest_atom(a->conjugation)) do_ionic = false;
        }
        else if (a->conjugation && a != a->conjugation->get_nearest_atom(b->get_location())) do_ionic = false;
        else if (b->conjugation && b != b->conjugation->get_nearest_atom(a->get_location())) do_ionic = false;
    }
    if (do_ionic)
    {
        retval[j] = &intertmp;
        retval[j]->Za = a->get_Z();
        retval[j]->Zb = b->get_Z();
        retval[j]->type = ionic;
        retval[j]->kJ_mol = 60; // Do not multiply by sgn charges here or total_binding() will reverse it.
        retval[j]->distance = 0.584 * (a->get_vdW_radius() + b->get_vdW_radius());		// Based on NH...O and the vdW radii of O and H.
        // retval[j]->distance = 0.7 * (a->get_vdW_radius() + b->get_vdW_radius());		// Based on NaCl.
        retval[j]->dirprop = 0;

        j++;
    }
    retval[j] = nullptr;

    #if allow_auto_hydroxy
    Atom *H = nullptr, *O = nullptr;
    Bond *brot = nullptr;
    if (a->get_Z() == 1 && b->get_Z() == 8)
    {
        Atom* HO = a->is_bonded_to("O");
        if (HO)
        {
            Bond* bb[16];
            HO->fetch_bonds(bb);
            if (bb[0])
            {
                for (i=0; bb[i]; i++)
                {
                    if (bb[i]->atom2 != a && bb[i]->can_rotate)
                    {
                        brot = bb[i];
                        H = a;
                        O = b;
                        break;
                    }
                }
            }
        }
    }
    else if (b->get_Z() == 1 && a->get_Z() == 8)
    {
        Atom* HO = b->is_bonded_to("O");
        if (HO)
        {
            Bond* bb[16];
            HO->fetch_bonds(bb);
            if (bb[0])
            {
                for (i=0; bb[i]; i++)
                {
                    if (bb[i]->atom2 != b && bb[i]->can_rotate)
                    {
                        brot = bb[i];
                        H = b;
                        O = a;
                        break;
                    }
                }
            }
        }
    }

    if (H && O && brot)
    {
        float rad, step = 30*fiftyseventh, bestr=999999, bestrad=0;

        for (rad=0; rad < M_PI*2; rad += step)
        {
            float r = H->distance_to(O);
            if (r < bestr && r >= 1.6)			// WARNING: Hard coded value.
            {
                bestrad = rad;
                bestr = r;
            }

            brot->rotate(step);
        }

        if (bestrad) brot->rotate(bestrad);
    }
    #endif

    for (i=0; look[i]; i++)
    {
        if (	(	(look[i]->Za == Za || look[i]->Za == any_element)
                    &&
                    (	!look[i]->bZa
                        ||
                        (!look[i]->aritybZa && a->is_bonded_to(Atom::esym_from_Z(look[i]->bZa)))
                        ||
                        ( look[i]->aritybZa && a->is_bonded_to(Atom::esym_from_Z(look[i]->bZa), look[i]->aritybZa))
                    )
                    &&
                    (look[i]->Zb == Zb || look[i]->Zb == any_element)
                    &&
                    (	!look[i]->bZb
                        ||
                        (!look[i]->aritybZb && b->is_bonded_to(Atom::esym_from_Z(look[i]->bZb)))
                        ||
                        ( look[i]->aritybZb && b->is_bonded_to(Atom::esym_from_Z(look[i]->bZb), look[i]->aritybZb))
                    )
             )
                ||
                (	(look[i]->Zb == Za || look[i]->Zb == any_element)
                    &&
                    (	!look[i]->bZb
                        ||
                        (!look[i]->aritybZb && a->is_bonded_to(Atom::esym_from_Z(look[i]->bZb)))
                        ||
                        ( look[i]->aritybZb && a->is_bonded_to(Atom::esym_from_Z(look[i]->bZb), look[i]->aritybZb))
                    )
                    &&
                    (look[i]->Za == Zb || look[i]->Za == any_element)
                    &&
                    (	!look[i]->bZa
                        ||
                        (!look[i]->aritybZa && b->is_bonded_to(Atom::esym_from_Z(look[i]->bZa)))
                        ||
                        ( look[i]->aritybZa && b->is_bonded_to(Atom::esym_from_Z(look[i]->bZa), look[i]->aritybZa))
                    )
                )
           )
        {
            switch (look[i]->type)
            {
            case covalent:
                if (a->is_bonded_to(b) == look[i]->arity)
                    retval[j++] = look[i];
                break;

            case ionic:
                if (	(a->get_acidbase() && b->get_acidbase())
                        ||
                        (a->get_charge() && b->get_charge())
                        ||
                        ((a == b) && a->get_acidbase())
                   )
                    retval[j++] = look[i];
                break;

            case hbond:
                retval[j++] = look[i];
                break;

            case pi:
                if (a->is_pi() && b->is_pi())
                    retval[j++] = look[i];
                break;

            case polarpi:
                if (  (a->is_pi() && (b->is_polar() || b->is_metal()) )
                        ||
                        (b->is_pi() && (a->is_polar() || a->is_metal()) )
                   )
                    retval[j++] = look[i];
                break;

            case mcoord:
                if (a->is_metal() || b->is_metal())
                    retval[j++] = look[i];
                break;

            case vdW:
                if (!j && !a->get_charge() && !b->get_charge())
                    retval[j++] = look[i];
                break;

            default:
                ;
            }
        }
    }
    retval[j] = nullptr;
}

SCoord* get_geometry_for_pi_stack(SCoord* in_geo)
{
    SCoord* retval = new SCoord[5];
    int i;
    Point pt[5];

    for (i=0; i<3; i++)
    {
        retval[i] = in_geo[i];
        pt[i] = in_geo[i];
    }

    pt[3] = compute_normal(&pt[0], &pt[1], &pt[2]);
    pt[3].scale(1);
    pt[4] = pt[3];
    pt[4].negate();

    retval[3] = pt[3];
    retval[4] = pt[4];

    return retval;
}

float InteratomicForce::metal_compatibility(Atom* a, Atom* b)
{
    // This metal compatibility estimate is based on the electronegativity of the metal and that of the coordinating atom.
    // It is designed to preferentially pair e.g. K/Na/Ca with O, Mg/Zn/Fe with N, Cu/Zn with S/N, Ag with S, but not e.g. Cu with O or Na with S.
    // Real world chelates, including metal binding sites in proteins, show a strong tendency to match electronegativities in this way.
    float ea = a->get_electronegativity();
    float eb = b->get_electronegativity();
    float sum = ea + eb;
    float delta = fabs(4.5 - sum);
    float f = fmax(0, 1.0 - delta);

    #if _dbg_groupsel
    // cout << "Metal compatibility for " << *a << "..." << *b << " = " << f << endl;
    #endif
    return f;
}

float InteratomicForce::potential_binding(Atom* a, Atom* b)
{
    InteratomicForce* forces[32];
    fetch_applicable(a, b, forces);

    int i;
    float potential = 0;

    for (i=0; forces[i]; i++)
    {
        if (!forces[i]->distance) continue;

        float partial = forces[i]->kJ_mol;

        if (forces[i]->type == hbond)
        {
            partial *= fmin(fabs(a->is_polar()), fabs(b->is_polar()));
        }

        if (forces[i]->type == mcoord)
        {
            partial *= metal_compatibility(a, b);
        }

        potential += partial;
    }

    if (!a->is_metal() && !b->is_metal())
    {
        // Oil and water don't mix.
        if ((fabs(a->is_polar()) < hydrophilicity_cutoff && (fabs(b->is_polar()) >= hydrophilicity_cutoff || fabs(b->get_charge())))
            ||
            (fabs(b->is_polar()) < hydrophilicity_cutoff && (fabs(a->is_polar()) >= hydrophilicity_cutoff || fabs(a->get_charge())))
        )
        {
            potential -= ((fabs(a->is_polar()) < hydrophilicity_cutoff && a->is_pi()) || (fabs(b->is_polar()) < hydrophilicity_cutoff && b->is_pi())) ? 66 : 99;
        }
    }

    return potential;
}

#define _num_force_precedences 6
const intera_type force_precedence[_num_force_precedences] = {mcoord, ionic, hbond, pi, polarpi, vdW};
Interaction InteratomicForce::total_binding(Atom* a, Atom* b)
{
    InteratomicForce* forces[32];
    fetch_applicable(a, b, forces);

    int i, j, k;
    Interaction kJmol = 0;
    float partial = 0;

    InteratomicForce* forces_by_type[_num_force_precedences];
    for (i=0; i<_num_force_precedences; i++) forces_by_type[i] = nullptr;
    for (i=0; forces[i]; i++)
    {
        for (j=0; j<_num_force_precedences; j++)
        {
            if (force_precedence[j] == forces[i]->type)
            {
                forces_by_type[j] = forces[i];
                break;
            }
        }
    }
    j=0;
    for (i=0; i<_num_force_precedences; i++)
    {
        if (forces_by_type[i]) forces[j++] = forces_by_type[i];
    }
    forces_by_type[j] = nullptr;


    float r = a->distance_to(b);
    float avdW = a->get_vdW_radius(), bvdW = b->get_vdW_radius();
    float rbind = avdW+bvdW;

    float dp = 2;			// Directional propensity. The anisotropic component from each single vertex
                            // is calculated as a cosine and then raised to this exponent.
    Point center;

    float achg = a->get_charge(), bchg = b->get_charge()
        , apol = a->is_polar(), bpol = b->is_polar();
    int aZ = a->get_Z(), bZ = b->get_Z();
    bool api = a->is_pi(), bpi = b->is_pi();

    #if auto_pK_protonation
    if (a->is_pKa_near_bio_pH() && bchg < 0) achg = 1;
    if (b->is_pKa_near_bio_pH() && achg < 0) bchg = 1;
    #endif

    Atom* aheavy = a;
    if (aheavy->get_Z() == 1)
    {
        Bond* ab = a->get_bond_by_idx(0);
        if (ab->atom2 && ab->atom2->get_Z() > 1) aheavy = ab->atom2;
    }
    Atom* bheavy = b;
    if (bheavy->get_Z() == 1)
    {
        Bond* bb = b->get_bond_by_idx(0);
        if (bb->atom2 && bb->atom2->get_Z() > 1) bheavy = bb->atom2;
    }
    float rheavy = aheavy->distance_to(bheavy);
    float l_heavy_atom_mindist = aheavy->get_vdW_radius() + bheavy->get_vdW_radius();

    float ahcg = aheavy->is_conjugated_to_charge(), bhcg = bheavy->is_conjugated_to_charge();
    if (fabs(ahcg) > fabs(achg)) achg = ahcg;
    if (fabs(bhcg) > fabs(bchg)) bchg = bhcg;
    // if (!strcmp(a->name, "O6") && !strcmp(b->name, "HD1") && b->residue == 180 ) cout << ahcg << " " << bhcg << endl;

    #if _ALLOW_PROTONATE_PNICTOGENS

    // TODO: Increase this value if multiple negative charges are nearby; decrease if positive nearby.
    if (!achg && bchg < 0 && aheavy->get_family() == PNICTOGEN && !aheavy->is_amide() && !isnan(aheavy->pK))
    {
        achg = protonation(aheavy->pK) * fabs(bchg) / pow(fabs(r-1.5)+1, 2);
    }

    if (!bchg && achg < 0 && bheavy->get_family() == PNICTOGEN && !bheavy->is_amide() && !isnan(bheavy->pK))
    {
        bchg = protonation(bheavy->pK) * fabs(achg) / pow(fabs(r-1.5)+1, 2);
    }
    #endif

    if (achg && sgn(achg) == sgn(bchg))
    {
        float repulsion = charge_repulsion * achg*bchg / pow(r, 2);
        kJmol.repulsive += repulsion;
        k = (ionic - covalent) % _INTER_TYPES_LIMIT;
        total_binding_by_type[k] -= repulsion;
    }

    if (achg) apol += achg;
    if (bchg) bpol += bchg;

    if (apol && sgn(apol) == sgn(bpol))
    {
        float pr = polar_repulsion / pow(r, 2) * fabs(apol) * fabs(bpol);

        if (a->get_Z() == 1)
        {
            Bond* prb = a->get_bond_by_idx(0);
            if (prb && prb->atom2)
            {
                float prtheta = find_3d_angle(b->get_location(), prb->atom2->get_location(), a->get_location());
                pr *= 0.5 + 0.5 * cos(prtheta);
            }
        }
        else
        {
            Bond* sb = a->get_bond_closest_to(b->get_location());
            if (sb && sb->atom2 && sb->atom2->distance_to(b) < (r - 0.5 * sb->optimal_radius) ) goto no_polar_repuls;
        }

        if (b->get_Z() == 1)
        {
            Bond* prb = b->get_bond_by_idx(0);
            if (prb && prb->atom2)
            {
                float prtheta = find_3d_angle(a->get_location(), prb->atom2->get_location(), b->get_location());
                pr *= 0.5 + 0.5 * cos(prtheta);
            }
        }
        else
        {
            Bond* sb = b->get_bond_closest_to(a->get_location());
            if (sb && sb->atom2 && sb->atom2->distance_to(a) < (r - 0.5 * sb->optimal_radius) ) goto no_polar_repuls;
        }

        kJmol.repulsive += pr;
        k = (hbond - covalent) % _INTER_TYPES_LIMIT;
        total_binding_by_type[k] -= pr;
    }

    no_polar_repuls:
    ;

    bool atoms_are_bonded = a->is_bonded_to(b);

    if (a->residue && b->residue && a->is_backbone && b->is_backbone)
        if (abs(a->residue - b->residue) == 1)
            atoms_are_bonded = true;

    #if _ALLOW_PROTONATE_PNICTOGENS
    if (forces)
        for (i=0; forces[i]; i++)
        {
            if (forces[i]->type == ionic) goto _has_ionic_already;
        }

    if (!atoms_are_bonded && achg && bchg && sgn(achg) == -sgn(bchg))
    {
        float pcf = 60.0 * fabs(achg)*fabs(bchg) / pow(r/((avdW+bvdW)*0.6), 2);
        kJmol.attractive += pcf;
        k = (ionic - covalent) % _INTER_TYPES_LIMIT;
        total_binding_by_type[k] += pcf;
    }
    #endif

    _has_ionic_already:
    
    if (!forces)
    {
        goto _canstill_clash;
    }

    if (r < 0.5) forces[0] = NULL;

    float atheta, btheta;

    for (i=0; forces[i]; i++)
    {
        if (forces[i]->type == ionic)
        {
            if (a->conjugation && b->conjugation)
            {
                if (a != a->conjugation->get_nearest_atom(b->conjugation)) continue;
                if (b != b->conjugation->get_nearest_atom(a->conjugation)) continue;
            }
            else if (a->conjugation && a != a->conjugation->get_nearest_atom(b->get_location())) continue;
            else if (b->conjugation && b != b->conjugation->get_nearest_atom(a->get_location())) continue;
        }

        // https://chemistry.stackexchange.com/questions/42085/can-an-amide-nitrogen-be-a-hydrogen-bond-acceptor
        if (forces[i]->type == hbond)
        {
            if  (   (aZ == 1 && b->get_family() == PNICTOGEN && (bpi && b->get_bonded_atoms_count() > 2) ) 
                 || (bZ == 1 && a->get_family() == PNICTOGEN && (api && a->get_bonded_atoms_count() > 2) )
                 || (aZ == 1 && bchg > hydrophilicity_cutoff)
                 || (bZ == 1 && achg > hydrophilicity_cutoff)
                )
                continue;
        }

        if (!forces[i]->distance) continue;
        float r1 = r / forces[i]->distance;

        #if _dbg_interatomic_forces
        bool debug_criteria = forces[i]->type != vdW && !a->is_backbone && !b->is_backbone;
        if (debug_criteria) cout << a->name << " ~ " << *forces[i] << " ~ " << b->name << endl;
        #endif

        if (forces[i]->type == covalent) continue;
        
        float rdecayed;
        float asum=0, bsum=0, aniso=1;
        bool stacked_pi_rings = false;

        #if _enhanced_pi_stacking
        if (forces[i]->type == pi && api && bpi)
        {
            Ring** ar = aheavy->get_rings();
            Ring** br = bheavy->get_rings();

            partial = 0;

            if (ar && br)
            {
                for (j=0; ar[j]; j++)
                {
                    if (!ar[j]->is_conjugated() || !ar[j]->is_coplanar()) continue;

                    SCoord anorm = ar[j]->get_normal();
                    for (k=0; br[k]; k++)
                    {
                        if (!br[k]->is_conjugated() || !br[k]->is_coplanar()) continue;
                        float lpart = 0;

                        SCoord bnorm = br[k]->get_normal();
                        dp = forces[i]->get_dp();
                        float norm_theta = find_3d_angle(anorm, bnorm, Point(0,0,0));
                        lpart = forces[i]->kJ_mol * pow(0.5+0.5*cos(norm_theta*2), dp);
                        lpart /= pi_mult_dkytw;

                        // Sandwiched and parallel displaced stacking. https://pubs.acs.org/doi/10.1021/jp912094q
                        if ((a->get_Z() == 1 && b->get_Z() >  1)
                            ||
                            (a->get_Z() >  1 && b->get_Z() == 1)
                           )
                        {
                            Point projnorm = anorm;
                            projnorm.scale(r);
                            Point proj1 = a->get_location().add(projnorm);
                            Point proj2 = a->get_location().subtract(projnorm);
                            float disp = fmin(b->get_location().get_3d_distance(proj1), b->get_location().get_3d_distance(proj2));

                            lpart *= fmax(0, 1.0-disp) * pi_mult_dkytw * pi_CH_dkytw;
                        }

                        // T-shaped stacking, ibid.
                        if (a->get_Z() == 1 && b->get_Z() > 1 && norm_theta > 60*fiftyseventh)
                        {
                            Point bcen = br[k]->get_center();
                            SCoord proj = a->get_location().subtract(aheavy->get_location());
                            proj.r = bcen.get_3d_distance(a->get_location());

                            float disp = a->get_location().add(proj).get_3d_distance(bcen);

                            lpart += forces[i]->kJ_mol * fmax(0, 1.0-disp) * pi_mult_dkytw * pi_HT_dkytw;
                        }
                        else if (b->get_Z() == 1 && a->get_Z() > 1 && norm_theta > 60*fiftyseventh)
                        {
                            Point acen = ar[j]->get_center();
                            SCoord proj = b->get_location().subtract(bheavy->get_location());
                            proj.r = acen.get_3d_distance(b->get_location());

                            float disp = b->get_location().add(proj).get_3d_distance(acen);

                            lpart += forces[i]->kJ_mol * fmax(0, 1.0-disp) * pi_mult_dkytw * pi_HT_dkytw;
                        }

                        stacked_pi_rings = true;
                        rbind = forces[i]->distance;

                        if (r >= rbind) lpart /= pow(r/rbind, 2);
                        else lpart *= (r/rbind);
                        partial += lpart;

                        // cout << "Partial: " << lpart << endl;
                    }
                }
            }

            if (ar) delete[] ar;
            if (br) delete[] br;
        }
        #endif

        if (!stacked_pi_rings)
        {
            if (forces[i]->type == pi && (a->shielding_angle >= _shield_angle_pi || b->shielding_angle >= _shield_angle_pi)) continue;

            dp = 2;
            if (forces[i]->type == vdW) dp = 1;					// van der Waals forces are non-directional, but the C-H bond still shields.
            else if (forces[i]->get_dp()) dp = forces[i]->get_dp();

            #if _dbg_interatomic_forces
            if (debug_criteria) cout << "Anisotropic directional propensity: " << dp << endl;
            #endif

            // Anisotropy.
            SCoord* ageo = a->get_geometry_aligned_to_bonds();
            SCoord* bgeo = b->get_geometry_aligned_to_bonds();
            bool del_ageo=false, del_bgeo=false;
            int ag = a->get_geometry();
            int bg = b->get_geometry();
            int abc = a->get_bonded_atoms_count();
            int bbc = b->get_bonded_atoms_count();

            if (forces[i]->type == pi && ag >= 3 && bg >= 3)
            {
                if (del_ageo) delete[] ageo;
                if (del_bgeo) delete[] bgeo;
                ageo = get_geometry_for_pi_stack(ageo);
                bgeo = get_geometry_for_pi_stack(bgeo);
                ag = bg = 5;
                del_ageo = del_bgeo = true;
            }
            else if ((forces[i]->type == polarpi || forces[i]->type == mcoord) && ag >= 3 && bg >= 3)
            {
                if (!a->is_polar())
                {
                    if (del_ageo) delete[] ageo;
                    ageo = get_geometry_for_pi_stack(ageo);
                    ag = 5;
                    del_ageo = true;
                }
                if (!b->is_polar())
                {
                    if (del_bgeo) delete[] bgeo;
                    bgeo = get_geometry_for_pi_stack(bgeo);
                    bg = 5;
                    del_bgeo = true;
                }
            }

            float dpa, dpb;

            dpa = dpb = dp;

            if (forces[i]->type == ionic || forces[i]->type == hbond || forces[i]->type == mcoord)
            {
                if (a->is_polar() < 0 && b->is_polar() >= 0)
                {
                    dpa = dp;
                    dpb = (forces[i]->type == hbond) ? 3 : 1;
                }
                else if (b->is_polar() < 0 && a->is_polar() >= 0)
                {
                    dpb = dp;
                    dpa = (forces[i]->type == hbond) ? 3 : 1;
                }
                else dpa = dpb = dp;
            }

            // Hydrogens have very low directional propensities.
            // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7325729/bin/NIHMS1601788-supplement-1.pdf
            if (a->get_Z() == 1) dpa = 0.75;
            if (b->get_Z() == 1) dpb = 0.75;

            Point aloc = a->get_location();
            Point bloc = b->get_location();

            int anx = a->get_idx_next_free_geometry();
            int bnx = b->get_idx_next_free_geometry();

            ag = abs(ag);
            bg = abs(bg);
            SCoord avec[ag], bvec[bg];

            Ring *ar = nullptr, *br = nullptr;

            for (j=0; j<ag; j++)
            {
                avec[j] = ageo[j];
                Bond* jb = a->get_bond_by_idx(j);
                if (jb && jb->atom2) avec[j].r = 0;
            }

            for (j=0; j<bg; j++)
            {
                bvec[j] = bgeo[j];
                Bond* jb = b->get_bond_by_idx(j);
                if (jb && jb->atom2) bvec[j].r = 0;
            }

            #if _dbg_259
            /*if (a->get_family() == PNICTOGEN && b->get_Z() == 1)
            {
                for (j=0; j<ag; j++)
                {
                    cout << a->name << " vertex " << j;
                    Bond* b259 = a->get_bond_by_idx(j);
                    if (b259 && b259->atom2) cout << " occupied by " << b259->atom2->name;
                    else cout << " vacant";

                    float th259 = find_3d_angle(aloc.add(ageo[j]), bloc, aloc);
                    cout << ", " << (th259*fiftyseven) << "deg from " << b->name;

                    cout << endl;
                }
            }*/
            #endif

            if (del_ageo) delete[] ageo;
            if (del_bgeo) delete[] bgeo;

            // When pi-bonding to a heavy atom of a conjugated coplanar ring, treat the entire ring as if it were one atom.
            if ((forces[i]->type == pi || forces[i]->type == polarpi)
                    &&
                    a->get_Z() > 1
                    &&
                    a->num_conj_rings()
            )
            {
                // Find the coplanar ring closest to bloc, if any.
                ar = a->closest_arom_ring_to(bloc);
                if (ar)
                {
                    aloc = ar->get_center();
                    avec[0] = ar->get_normal();
                    avec[1] = avec[0];
                    avec[1] = avec[1].negate();
                    for (j=2; j<ag; j++) avec[j] = SCoord(0,0,0);
                    ag = 2;
                }
            }

            if ((forces[i]->type == pi || forces[i]->type == polarpi)
                    &&
                    b->get_Z() > 1
                    &&
                    b->num_conj_rings()
            )
            {
                // Find the coplanar ring closest to aloc, if any.
                br = b->closest_arom_ring_to(aloc);
                if (br)
                {
                    bloc = br->get_center();
                    bvec[0] = br->get_normal();
                    bvec[1] = bvec[0];
                    bvec[1] = bvec[1].negate();
                    for (j=2; j<ag; j++) bvec[j] = SCoord(0,0,0);
                    bg = 2;
                }
            }

            if (forces[i]->type != ionic)
            {
                if (a->get_bonded_atoms_count() >= ag || b->get_bonded_atoms_count() >= bg)
                    continue;

                // Sum up the anisotropic contribution from each geometry vertex of a.
                for (j=0; j<ag; j++)
                {
                    if (!avec[j].r) continue;
                    Point pt(&avec[j]);
                    pt.scale(r);
                    pt = pt.add(aloc);
                    atheta = find_3d_angle(&bloc, &pt, &aloc);
                    if (forces[i]->type != pi && forces[i]->type != polarpi) if (atheta > M_PI/2)
                    {
                        #if _dbg_interatomic_forces
                        if (debug_criteria) cout << a->name << " anisotropic angle " << (atheta*fiftyseven)
                            << "deg, out of range." << endl; 
                        #endif
                        continue;
                    }
                    float contrib = pow(fmax(0,cos(atheta)), dpa);
                    if (!isnan(contrib) && !isinf(contrib)) asum += contrib;

                    #if _dbg_interatomic_forces
                    if (debug_criteria) cout << a->name << " anisotropic angle " << (atheta*fiftyseven)
                        << "deg with dpa = " << dpa
                        << ", contributing " << contrib
                        << endl;
                    #endif
                }

                // Sum up the anisotropic contribution from each geometry vertex of b.
                for (j=0; j<bg; j++)
                {
                    if (!bvec[j].r) continue;
                    Point pt(&bvec[j]);
                    pt.scale(r);
                    pt = pt.add(bloc);
                    btheta = find_3d_angle(&aloc, &pt, &bloc);
                    if (forces[i]->type != pi && forces[i]->type != polarpi) if (btheta > M_PI/2)
                    {
                        #if _dbg_interatomic_forces
                        if (debug_criteria) cout << b->name << " anisotropic angle " << (btheta*fiftyseven)
                            << "deg, out of range." << endl; 
                        #endif
                        continue;
                    }
                    float contrib = pow(fmax(0,cos(btheta)), dpb);
                    if (!isnan(contrib) && !isinf(contrib)) bsum += contrib;

                    #if _dbg_interatomic_forces
                    if (debug_criteria) cout << b->name << " anisotropic angle " << (btheta*fiftyseven)
                        << "deg with dpb = " << dpb
                        << ", contributing " << contrib
                        << endl;
                    #endif
                }

                asum = fmin(1, fmax(0, fabs(asum)));
                bsum = fmin(1, fmax(0, fabs(bsum)));

                // Multiply the two sums.
                aniso = fmax(minimum_searching_aniso, asum * bsum);
            }

            #if _dbg_interatomic_forces
            if (debug_criteria) cout << "Anisotropy: " << aniso << " from " << asum << " * " << bsum << endl;
            #endif

            float force_eff_kJmol = forces[i]->kJ_mol;
            if (forces[i]->type == ionic && sgn(achg) == sgn(bchg)) continue;

            if (r1 >= 0.99)
            {
                rdecayed = (forces[i]->type == vdW
                            ||
                            forces[i]->type == pi
                            ||
                            forces[i]->type == polarpi
                        )
                        ? r1*r1*r1*r1*r1*r1
                        : r1*r1;
                partial = aniso * force_eff_kJmol / rdecayed;

                SCoord mc = bloc.subtract(aloc);
                mc.r = r - forces[i]->distance; // fabs((r1 - 1) * (1.0 - (partial / forces[i]->kJ_mol)) / (r1*r1));

                #if compute_missed_connections
                #if summed_missed_connections
                missed_connection = missed_connection.add(mc);
                mc_bpotential += force_eff_kJmol;
                #else
                if ((mc.r < fabs(missed_connection.r) && force_eff_kJmol >= mc_bpotential)
                    || force_eff_kJmol > mc_bpotential)
                {
                    missed_connection = mc;
                    mc_bpotential = force_eff_kJmol;
                }
                #endif
                #endif
            }
            else
            {
                // rdecayed = r1*r1*r1*r1*r1*r1;
                partial = aniso * force_eff_kJmol - Lennard_Jones(a, b, forces[i]->get_distance());
            }

            if (forces[i]->type == hbond)
            {
                // https://www.sciencedirect.com/science/article/abs/pii/S0009261497011172
                // https://web.archive.org/web/20200305164852id_/https://boris.unibe.ch/134571/1/1NpOH-Hydbond_JCP_resub.pdf
                // TODO: Replace this with better data in bindings.dat. THE FOLLOWING IS A GROSS OVERSIMPLIFICATION.
                if (api && !a->is_amide() && !bpi) partial *= 21.8 / 37.6;
                if (bpi && !b->is_amide() && !api) partial *= 21.8 / 37.6;
            }

            if (forces[i]->type == mcoord)
            {
                partial *= metal_compatibility(a, b);
            }

            // Divide each ring by its number of atoms.
            if (forces[i]->type == pi)
            {
                if (ar) partial /= ar->get_atom_count();
                else partial /= 6;

                if (br) partial /= br->get_atom_count();
                else partial /= 6;
            }
            else if (forces[i]->type == polarpi) partial /= 6;

            if (fabs(partial) > fabs(forces[i]->kJ_mol*2) || partial >= 500)
            {                
                #if ignore_invalid_partial
                return 0;
                #endif
                cout << "Invalid partial! " << partial << " (max " << forces[i]->kJ_mol << ") from "
                    << a->name << "..." << b->name << " r=" << r
                    << " (optimal " << forces[i]->distance << ") rdecayed=" << rdecayed
                    << " aniso=" << aniso << " (" << asum << "*" << bsum << ")" << endl;
                throw 0xbadf0ace;
            }

            if (forces[i]->type == polarpi || forces[i]->type == mcoord || forces[i]->type == ionic)
            {
                partial *= fmax(forces[i]->type == ionic?0:1, fabs(achg)) * fmax(forces[i]->type == ionic?0:1, fabs(bchg));
            }

            if (forces[i]->type == hbond && fabs(apol) && fabs(bpol)) partial *= fabs(apol) * fabs(bpol);

            # if 0
            //if (forces[i]->type == polarpi || forces[i]->type == mcoord)
            if (forces[i]->type == hbond)
            {
                cout << a->name << " ... " << b->name << " " << forces[i]->type << " partial " << partial
                    << " (" << *forces[i] << ") "
                    << " rdecayed " << rdecayed << " aniso " << aniso << " ag " << ag << " bg " << bg
                    << endl;
            }
            #endif
        }
        
        #if _dbg_interatomic_forces
        if (debug_criteria) cout << -partial << endl;
        #endif

        #if _dbg_51e2_ionic
        if (forces[i]->type == ionic && achg < 0 && bchg > 0 && !a->residue && b->residue == 262)
        {
            cout << a->name << " charge " << achg << " ~ "
                << r << " ~ " << b->residue << ":" << b->name
                << ", max energy " << forces[i]->kJ_mol
                << ", optimal distance " << forces[i]->distance
                << "; aniso " << aniso << " = " << partial
                << endl;
        }
        #endif

        if (fabs(partial) >= 500)
        {
            #if ignore_invalid_partial
            return 0;
            #else
            cout << "Invalid partial! " << partial << " (max " << forces[i]->kJ_mol << ") from "
                 << a->name << "...." << b->name << " r=" << r
                 << " (optimal " << forces[i]->distance << ") rdecayed=" << rdecayed
                 << " aniso=" << aniso << " (" << asum << "*" << bsum << ")" << endl;
            throw 0xbadf0ace;
            #endif
        }

        kJmol.attractive += partial;
        if (partial > 0.5 && forces[i]->distance < rbind) 
        {
            float f = forces[i]->distance / rbind;
            avdW *= f;
            bvdW *= f;
            rbind = forces[i]->distance;
        }

        #if _dbg_interatomic_forces
        if (debug_criteria) cout << endl;
        #endif

        #if _peratom_audit
        if (interauditing)
        {
            std::string str;
            if (a->residue) str += (std::string)a->aa3let + to_string(a->residue) + (std::string)":";
            str += (std::string)a->name + (std::string)"-";
            if (b->residue) str += (std::string)b->aa3let + to_string(b->residue) + (std::string)":";
            str += (std::string)b->name + (std::string)" ";
            switch (forces[i]->type)
            {
                case covalent:
                str += (std::string)"cov";
                break;

                case ionic:
                str += (std::string)"ion";
                break;

                case hbond:
                str += (std::string)"hb";
                break;

                case pi:
                str += (std::string)"pi";
                break;

                case polarpi:
                str += (std::string)"ppi";
                break;

                case mcoord:
                str += (std::string)"mtl";
                break;

                case vdW:
                str += (std::string)"vdw";
                break;

                default:
                str += (std::string)"unk";
            }

            str += (std::string)" " + to_string(-partial) + (std::string)" (" + to_string(-kJmol.summed()) + (std::string)")";

            str += (std::string)" theta: " + to_string(atheta*fiftyseven) + (std::string)", " + to_string(btheta*fiftyseven);

            interaudit.push_back(str);
        }
        #endif

        k = (forces[i]->type - covalent) % _INTER_TYPES_LIMIT;
        total_binding_by_type[k] += partial;

        if (forces[i]->type == ionic && achg && bchg)
        {
            break;
        }
    
        if (partial) break;
    }

    if (rbind < 0.7) rbind = 0.7;

_canstill_clash:
    float sigma;
    float local_clash_allowance = global_clash_allowance;

    #if ignore_double_hydrogen_clashes
    if (a->get_Z() == 1 && b->get_Z() == 1 && r > 1.0)
    {
        float rheavy = fmin(a->distance_to(bheavy), b->distance_to(aheavy));
        if (rheavy < r) r = rheavy;
        else goto _finished_clashing;
    }
    #elif ignore_nonpolar_hydrogen_clashes                   // The leading molecular docker seems to.
    if ((a->get_Z() == 1 || b->get_Z() == 1)
        && r > 1.0
        && fabs(a->is_polar()) < hydrophilicity_cutoff
        && fabs(b->is_polar()) < hydrophilicity_cutoff
        )
    {
        float rheavy = fmin(a->distance_to(bheavy), b->distance_to(aheavy));
        if (rheavy < r) r = rheavy;
        else goto _finished_clashing;
    }
    #endif

    if (a->get_Z() == 1 && b->get_Z() == 1) local_clash_allowance *= double_hydrogen_clash_allowance_multiplier;

    sigma = fmin(rbind, avdW+bvdW) - local_clash_allowance;

    float clash = 0;
    if (r < rbind && !atoms_are_bonded) // && (!achg || !bchg || (sgn(achg) != -sgn(bchg))) )
    {
        // if (!strcmp(a->name, "O6") && !strcmp(b->name, "HD1") && b->residue == 180 ) cout << achg << " " << bchg << endl;
        clash = fmax(Lennard_Jones(a, b, sigma), 0);
        kJmol.repulsive += clash;
        k = (vdW - covalent) % _INTER_TYPES_LIMIT;
        total_binding_by_type[k] -= clash;

        #if _peratom_audit
        if (interauditing)
        {
            std::string str;
            if (a->residue) str += (std::string)a->aa3let + to_string(a->residue) + (std::string)":";
            str += (std::string)a->name + (std::string)"-";
            if (b->residue) str += (std::string)b->aa3let + to_string(b->residue) + (std::string)":";
            str += (std::string)b->name + (std::string)" ";
            str += (std::string)"clash " + to_string(clash);

            interaudit.push_back(str);
        }
        #endif
    }

    _finished_clashing:

    return kJmol;
}

float InteratomicForce::Lennard_Jones(Atom* atom1, Atom* atom2, float sigma)
{
    float local_clash_allowance = global_clash_allowance;
    if (atom1->get_Z() == 1 && atom2->get_Z() == 1) local_clash_allowance *= double_hydrogen_clash_allowance_multiplier;

    if (!sigma) sigma = atom1->get_vdW_radius() + atom2->get_vdW_radius() - local_clash_allowance;
    float r = atom1->distance_to(atom2);
    float sigma_r = sigma / r;

    return Lennard_Jones_epsilon_x4 * (pow(sigma_r, 12) - 2.0*pow(sigma_r, 6));
}

bool Interaction::improved(Interaction rel)
{
    if (summed() < -clash_limit_per_aa*10) return rel.repulsive > repulsive;
    return attractive > rel.attractive || repulsive < rel.repulsive;
}

float InteratomicForce::distance_anomaly(Atom* a, Atom* b)
{
    InteratomicForce* forces[32];
    fetch_applicable(a, b, forces);

    int i;
    float anomaly = 0;
    float r = a->distance_to(b);
    for (i=0; forces[i]; i++)
    {
        if (!forces[i]->distance) continue;
        anomaly += fabs(r - forces[i]->distance);
    }

    return anomaly;
}

float InteratomicForce::covalent_bond_radius(Atom* a, Atom* b, float cardinality)
{
    if (!read_forces_dat && !reading_forces) read_all_forces();
    if (!a || !b || !a->get_Z() || !b->get_Z()) return 0;

    int i, j=0;
    for (i=0; all_forces[i]; i++)
    {
        if (all_forces[i]->type == covalent && fabs(all_forces[i]->arity - cardinality) < 0.1)
            if (	(	all_forces[i]->Za == a->get_Z()
                        &&
                        (	!all_forces[i]->bZa || a->is_bonded_to(Atom::esym_from_Z(all_forces[i]->bZa)))
                        &&
                        all_forces[i]->Zb == b->get_Z()
                        &&
                        (	!all_forces[i]->bZb || b->is_bonded_to(Atom::esym_from_Z(all_forces[i]->bZb)))
                 )
                    ||
                    (	all_forces[i]->Za == b->get_Z()
                        &&
                        (	!all_forces[i]->bZa || b->is_bonded_to(Atom::esym_from_Z(all_forces[i]->bZa)))
                        &&
                        all_forces[i]->Zb == a->get_Z()
                        &&
                        (	!all_forces[i]->bZb || a->is_bonded_to(Atom::esym_from_Z(all_forces[i]->bZb)))
                    )
               )
            {
                return all_forces[i]->distance;
            }
    }

    cout << "Error: Not found covalent bond definition " << a->get_elem_sym() << cardinality_printable(cardinality)
         << b->get_elem_sym() << " cardinality " << cardinality << "; please check bindings.dat." << endl;
    throw BOND_DEF_NOT_FOUND;
}

float InteratomicForce::optimal_distance(Atom* a, Atom* b)
{
    InteratomicForce* forces[32];
    fetch_applicable(a, b, forces);

    int i;
    float strongest = 0, optimal = -1;
    for (i=0; forces[i]; i++)
    {
        if (!forces[i]->distance) continue;
        if (forces[i]->kJ_mol > strongest)
        {
            strongest = forces[i]->kJ_mol;
            optimal = forces[i]->distance;
        }
    }

    if (optimal < 0) optimal = a->get_vdW_radius() + b->get_vdW_radius();
    return optimal;
}

float InteratomicForce::coordinate_bond_radius(Atom* a, Atom* b, intera_type btype)
{
    if (!read_forces_dat && !reading_forces) read_all_forces();

    int i, j=0;
    for (i=0; all_forces[i]; i++)
    {
        if (all_forces[i]->type == btype)
            if (	(	all_forces[i]->Za == a->get_Z()
                        &&
                        (	!all_forces[i]->bZa || a->is_bonded_to(Atom::esym_from_Z(all_forces[i]->bZa)))
                        &&
                        all_forces[i]->Zb == b->get_Z()
                        &&
                        (	!all_forces[i]->bZb || b->is_bonded_to(Atom::esym_from_Z(all_forces[i]->bZb)))
                 )
                    ||
                    (	all_forces[i]->Za == b->get_Z()
                        &&
                        (	!all_forces[i]->bZa || b->is_bonded_to(Atom::esym_from_Z(all_forces[i]->bZa)))
                        &&
                        all_forces[i]->Zb == a->get_Z()
                        &&
                        (	!all_forces[i]->bZb || a->is_bonded_to(Atom::esym_from_Z(all_forces[i]->bZb)))
                    )
               )
            {
                return all_forces[i]->distance;
            }
    }

    cout << "WARNING: No bond definition for " << a->get_elem_sym() << " " << btype
        << " " << b->get_elem_sym() << ". Please check data/bindings.dat." << endl << endl;
    // throw BOND_DEF_NOT_FOUND;
    return a->get_vdW_radius() + b->get_vdW_radius();
}

std::string InteratomicForce::get_config_string() const
{
    char buffer[128];
    const char* bondchrs[5] = {"", "-", "=", "#", "$"};
    sprintf(buffer, "%s%s%s %s%s%s %s %f %f %f",
            bZa ? Atom::esym_from_Z(bZa) : "",
            bondchrs[aritybZa],
            Atom::esym_from_Z(Za),
            bZb ? Atom::esym_from_Z(bZb) : "",
            bondchrs[aritybZb],
            Atom::esym_from_Z(Zb),
            (type == covalent) 		? "coval" :
            (type == ionic) 		? "ionic" :
            (type == hbond) 		? "hbond" :
            (type == pi) 		? "pi" :
            (type == polarpi) 	? "plpi" :
            (type == mcoord) 	? "coord" :
            (type == vdW) 	? "vdW" :
            "?",
            arity,
            distance,
            kJ_mol
           );
    std::string returnValue = buffer;
    return returnValue;
}

std::ostream& operator<<(std::ostream& os, const intera_type& it)
{
    switch (it)
    {
    case covalent:
        os << "covalent";
        break;
    case ionic:
        os << "ionic";
        break;
    case hbond:
        os << "hbond";
        break;
    case pi:
        os << "pi";
        break;
    case polarpi:
        os << "polarpi";
        break;
    case mcoord:
        os << "mcoord";
        break;
    case vdW:
        os << "vdW";
        break;
    default:
        os << "unknown";
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const InteratomicForce& f)
{
    os << f.get_config_string();
    return os;
}




