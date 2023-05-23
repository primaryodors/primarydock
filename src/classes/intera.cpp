
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

#if _peratom_audit
std::vector<std::string> interaudit;
bool interauditing = false;
#endif

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

        // cout << *this << endl;
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
InteratomicForce** InteratomicForce::get_applicable(Atom* a, Atom* b)
{
    if (!read_forces_dat && !reading_forces) read_all_forces();
    if (!a || !b)
    {
        return NULL;
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

    InteratomicForce** retval = new InteratomicForce*[16];
    init_nulls(retval, 16);
    int i, j=0;

    // Charged atoms always attract or repel, irrespective of Z.
    if (a->get_charge() && b->get_charge())
    {
        retval[j] = &intertmp;
        retval[j]->Za = a->get_Z();
        retval[j]->Zb = b->get_Z();
        retval[j]->type = ionic;
        retval[j]->kJ_mol = 20; // Do not multiply by sgn charges here or total_binding() will reverse it.
        retval[j]->distance = 0.584 * (a->get_vdW_radius() + b->get_vdW_radius());		// Based on NH...O and the vdW radii of O and H.
        retval[j]->dirprop = 0;

        j++;
    }

    #if allow_auto_hydroxy
    Atom *H = nullptr, *O = nullptr;
    Bond *brot = nullptr;
    if (a->get_Z() == 1 && b->get_Z() == 8)
    {
        Atom* HO = a->is_bonded_to("O");
        if (HO)
        {
            Bond** bb = HO->get_bonds();
            if (bb)
            {
                for (i=0; bb[i]; i++)
                {
                    if (bb[i]->btom != a && bb[i]->can_rotate)
                    {
                        brot = bb[i];
                        H = a;
                        O = b;
                        break;
                    }
                }

                delete[] bb;
            }
        }
    }
    else if (b->get_Z() == 1 && a->get_Z() == 8)
    {
        Atom* HO = b->is_bonded_to("O");
        if (HO)
        {
            Bond** bb = HO->get_bonds();
            if (bb)
            {
                for (i=0; bb[i]; i++)
                {
                    if (bb[i]->btom != b && bb[i]->can_rotate)
                    {
                        brot = bb[i];
                        H = b;
                        O = a;
                        break;
                    }
                }

                delete[] bb;
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
                if (a->get_family() == PNICTOGEN && (a->is_backbone || a->is_amide())) break;
                if (b->get_family() == PNICTOGEN && (b->is_backbone || b->is_amide())) break;
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
                if (!j)
                    retval[j++] = look[i];
                break;

            default:
                ;
            }
        }
    }
    retval[j] = 0;

    return retval;
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
    float f = (1.0 + 1.0 * cos(fmin(fabs(((a->get_electronegativity() + b->get_electronegativity()) / 2 - 2.25)*6), M_PI)));
    #if _dbg_groupsel
    // cout << "Metal compatibility for " << *a << "..." << *b << " = " << f << endl;
    #endif
    return f;
}

float InteratomicForce::potential_binding(Atom* a, Atom* b)
{
    InteratomicForce** forces = get_applicable(a, b);

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
        if ((fabs(a->is_polar()) < 0.333 && (fabs(b->is_polar()) >= 0.333 || fabs(b->get_charge())))
            ||
            (fabs(b->is_polar()) < 0.333 && (fabs(a->is_polar()) >= 0.333 || fabs(a->get_charge())))
        )
        {
            potential -= ((fabs(a->is_polar()) < 0.333 && a->is_pi()) || (fabs(b->is_polar()) < 0.333 && b->is_pi())) ? 66 : 99;
        }
    }

    return potential;
}

float InteratomicForce::total_binding(Atom* a, Atom* b)
{
    InteratomicForce** forces = get_applicable(a, b);

    int i, j, k;
    float kJmol = 0;

    float r = a->distance_to(b);
    float avdW = a->get_vdW_radius()*vdw_clash_allowance, bvdW = b->get_vdW_radius()*vdw_clash_allowance;
    float rbind = avdW+bvdW;

    float dp = 2;			// Directional propensity. The anisotropic component from each single vertex
                            // is calculated as a cosine and then raised to this exponent.
    Point center;

    float achg = a->get_charge(), bchg = b->get_charge()
        , apol = a->is_polar(), bpol = b->is_polar();
    
    if (!achg && a->get_Z() == 1)
    {
        Bond* lb = a->get_bond_by_idx(0);
        if (lb && lb->btom) achg = lb->btom->get_charge();
    }
    if (!bchg && b->get_Z() == 1)
    {
        Bond* lb = b->get_bond_by_idx(0);
        if (lb && lb->btom) bchg = lb->btom->get_charge();
    }
    
    #if _ALLOW_PROTONATE_PNICTOGENS
    Atom* aheavy = a;
    if (aheavy->get_Z() == 1)
    {
        /*Bond** tmpb = a->get_bonds();
        if (tmpb && tmpb[0] && tmpb[0]->btom) aheavy = tmpb[0]->btom;
        if (tmpb) delete[] tmpb;*/
        Bond* ab = a->get_bond_by_idx(0);
        if (ab->btom && ab->btom->get_Z() > 1) aheavy = ab->btom;
    }
    // TODO: Increase this value if multiple negative charges are nearby; decrease if positive nearby.
    if (!achg && bchg < 0 && aheavy->get_family() == PNICTOGEN && !isnan(aheavy->pK) && aheavy->pK > pn_protonation_pKa_min) // || (a->get_Z() == 1 && a->is_bonded_to(PNICTOGEN) )))
    {
        achg = pnictogen_partial_protonation * fabs(bchg) / pow(fabs(r-1.5)+1, 2);
        /* if (a->residue == 243) cout << "Partial protonation for " << a->residue << ":" << a->name
            << " due to proximity of " << b->residue << ":" << b->name << endl; */
    }
    Atom* bheavy = b;
    if (bheavy->get_Z() == 1)
    {
        Bond** tmpb = b->get_bonds();
        if (tmpb && tmpb[0] && tmpb[0]->btom) bheavy = tmpb[0]->btom;
        if (tmpb) delete[] tmpb;
    }
    // if (!bchg && achg < 0 && (b->get_family() == PNICTOGEN || (b->get_Z() == 1 && b->is_bonded_to(PNICTOGEN) )))
    if (!bchg && achg < 0 && bheavy->get_family() == PNICTOGEN && !isnan(bheavy->pK) && bheavy->pK > pn_protonation_pKa_min)
    {
        bchg = pnictogen_partial_protonation * fabs(achg) / pow(fabs(r-1.5)+1, 2);
        /* if (b->residue == 243) cout << "Partial protonation for " << b->residue << ":" << b->name
            << " due to proximity of " << a->residue << ":" << a->name << endl; */
    }
    #endif

    if (achg && sgn(achg) == sgn(bchg)) kJmol -= charge_repulsion * achg*bchg / pow(r, 2);
    if (apol>0 && sgn(apol) == sgn(bpol))
    {
        float pr = polar_repulsion / pow(r, 2) * fabs(apol) * fabs(bpol);

        if (a->get_Z() == 1)
        {
            Bond* prb = a->get_bond_by_idx(0);
            if (prb && prb->btom)
            {
                float prtheta = find_3d_angle(b->get_location(), prb->btom->get_location(), a->get_location());
                pr *= 0.5 + 0.5 * cos(prtheta);
            }
        }

        if (b->get_Z() == 1)
        {
            Bond* prb = b->get_bond_by_idx(0);
            if (prb && prb->btom)
            {
                float prtheta = find_3d_angle(a->get_location(), prb->btom->get_location(), b->get_location());
                pr *= 0.5 + 0.5 * cos(prtheta);
            }
        }

        kJmol -= pr;
    }

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

    if (!atoms_are_bonded && sgn(achg) == -sgn(bchg)) kJmol += 60.0 * fabs(achg)*fabs(bchg) / pow(r/2, 2);
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
        if (!forces[i]->distance) continue;
        float r1 = r / forces[i]->distance;

        if (forces[i]->type == covalent) continue;
        
        float partial, rdecayed;
        float asum=0, bsum=0, aniso=1;
        bool stacked_pi_rings = false;

        #if _enhanced_pi_stacking
        if (forces[i]->type == pi && a->is_pi() && b->is_pi())
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

            Point aloc = a->get_location();
            Point bloc = b->get_location();

            int anx = a->get_idx_next_free_geometry();
            int bnx = b->get_idx_next_free_geometry();

            ag = abs(ag);
            bg = abs(bg);
            SCoord avec[ag], bvec[bg];

            Ring *ar = nullptr, *br = nullptr;

            for (j=0; j<ag; j++)
                avec[j] = SCoord(0,0,0);
            for (j=anx; j<ag; j++)
                avec[j-anx] = ageo[j];

            for (j=0; j<bg; j++)
                bvec[j] = SCoord(0,0,0);
            for (j=bnx; j<bg; j++)
                bvec[j-bnx] = bgeo[j];

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
                if (a->get_bonded_atoms_count() == ag || b->get_bonded_atoms_count() == bg)
                    continue;

                // Sum up the anisotropic contribution from each geometry vertex of a.
                for (j=0; j<ag; j++)
                {
                    if (!avec[j].r) continue;
                    Point pt(&avec[j]);
                    pt.scale(r);
                    pt = pt.add(aloc);
                    atheta = find_3d_angle(&bloc, &pt, &aloc);
                    if (forces[i]->type != pi && forces[i]->type != polarpi) if (atheta > M_PI/2) continue;
                    float contrib = pow(fmax(0,cos(atheta)), dpa);
                    if (!isnan(contrib) && !isinf(contrib)) asum += contrib;
                }

                // Sum up the anisotropic contribution from each geometry vertex of b.
                for (j=0; j<bg; j++)
                {
                    if (!bvec[j].r) continue;
                    Point pt(&bvec[j]);
                    pt.scale(r);
                    pt = pt.add(bloc);
                    btheta = find_3d_angle(&aloc, &pt, &bloc);
                    if (forces[i]->type != pi && forces[i]->type != polarpi) if (btheta > M_PI/2) continue;
                    float contrib = pow(fmax(0,cos(btheta)), dpb);
                    if (!isnan(contrib) && !isinf(contrib)) bsum += contrib;
                }

                asum = fmin(1, fmax(0, fabs(asum)));
                bsum = fmin(1, fmax(0, fabs(bsum)));

                // Multiply the two sums.
                aniso = fmax(minimum_searching_aniso, asum * bsum);
                // cout << aniso << " | " << asum << " | " << bsum << " | " << dpa << " | " << dpb << endl;
            }

            if (r1 >= 1)
            {
                rdecayed = (forces[i]->type == vdW
                            ||
                            forces[i]->type == pi
                            ||
                            forces[i]->type == polarpi
                        )
                        ? r1*r1*r1*r1*r1*r1
                        : r1*r1;
                partial = aniso * forces[i]->kJ_mol / rdecayed;
            }
            else
            {
                rdecayed = r1 * r1 * r1;
                partial = aniso * forces[i]->kJ_mol * rdecayed;
            }

            // TODO: Replace this with a more generalized model of competitive h-bonding as well as ionic, mcoord, etc.
            if (forces[i]->type == hbond && a->is_backbone != b->is_backbone) partial *= 0.5;

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
                cout << "Invalid partial! " << partial << " (max " << forces[i]->kJ_mol << ") from "
                    << a->name << "..." << b->name << " r=" << r
                    << " (optimal " << forces[i]->distance << ") rdecayed=" << rdecayed
                    << " aniso=" << aniso << " (" << asum << "*" << bsum << ")" << endl;
                throw 0xbadf0ace;
            }

            if (forces[i]->type == polarpi || forces[i]->type == mcoord)
            {
                if (a->is_metal()) partial *= a->get_charge();
                if (b->is_metal()) partial *= b->get_charge();
            }

            if (forces[i]->type == ionic && a->get_charge() && b->get_charge())
            {
                partial *= fabs((a->get_charge()) * -(b->get_charge()));

                if (0 && (a->residue == 114 || b->residue == 114))
                    cout << "Ionic interaction between "
                        << a->name << " (" << a->get_charge() << ") and "
                        << b->name << " ("  << b->get_charge() << ")"
                        << " r=" << r
                        << " (optimal " << forces[i]->distance << ") rdecayed=" << rdecayed
                        << " aniso=" << aniso << " (" << asum << "*" << bsum << ") = "
                        << partial
                        << endl;
            }

            /*if (achg && bchg) partial = fabs(partial) * sgn(achg) * -sgn(bchg);
            if (achg && !bchg && bpol) partial = fabs(partial) * sgn(achg) * -sgn(bpol);
            if (!achg && apol && bchg) partial = fabs(partial) * sgn(apol) * -sgn(bchg);
            if (!achg && apol && !bchg && bpol) partial = fabs(partial) * sgn(apol) * -sgn(bpol);*/

            if (forces[i]->type == hbond) partial *= fabs(apol) * fabs(bpol);

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
        
        if (fabs(partial) >= 500)
        {
            cout << "Invalid partial! " << partial << " (max " << forces[i]->kJ_mol << ") from "
                 << a->name << "...." << b->name << " r=" << r
                 << " (optimal " << forces[i]->distance << ") rdecayed=" << rdecayed
                 << " aniso=" << aniso << " (" << asum << "*" << bsum << ")" << endl;
            throw 0xbadf0ace;
        }

        #if active_persistence
        if (partial > 0)
        {
            if (a->residue && !b->residue)
            {
                partial *= residue_binding_multiplier(a->residue);
                #if _DBG_RESBMULT
                if (residue_binding_multiplier(a->residue) > 1) std::cout << *a << "..." << *b << " partial " << partial << " multiplied." << endl;
                #endif
            }
            else if (!a->residue && b->residue)
            {
                partial *= residue_binding_multiplier(b->residue);
                #if _DBG_RESBMULT
                if (residue_binding_multiplier(a->residue) > 1) std::cout << *a << "..." << *b << " partial " << partial << " multiplied." << endl;
                #endif
            }
        }
        #endif

        kJmol += partial;
        if (partial > 0.5 && forces[i]->distance < rbind) rbind = forces[i]->distance;

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

            str += (std::string)" " + to_string(-partial) + (std::string)" (" + to_string(-kJmol) + (std::string)")";

            str += (std::string)" theta: " + to_string(atheta*fiftyseven) + (std::string)", " + to_string(btheta*fiftyseven);

            interaudit.push_back(str);
        }
        #endif

        k = (forces[i]->type - covalent) % _INTER_TYPES_LIMIT;
        total_binding_by_type[k] += partial;

        if (forces[i]->type == ionic && achg && bchg) break;
    }

    if (rbind < 0.7) rbind = 0.7;

_canstill_clash:
    /*float confidence = 2.5;		// TODO: Get this from the PDB.
    float give = 0.5;			// TODO: Compute this from the receptor secondary structure.

    float allowable = give + confidence / sqrt(3);

    r += allowable;*/

    if (r < rbind && !atoms_are_bonded)
    {
        float f = rbind/(avdW+bvdW);
        float clash = pow(fabs(sphere_intersection(avdW*f, bvdW*f, r)*_kJmol_cuA), 4);
        kJmol -= clash;
        
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

    delete[] forces;
    return kJmol;
}


float InteratomicForce::distance_anomaly(Atom* a, Atom* b)
{
    InteratomicForce** forces = get_applicable(a, b);

    int i;
    float anomaly = 0;
    float r = a->distance_to(b);
    for (i=0; forces[i]; i++)
    {
        if (!forces[i]->distance) continue;
        anomaly += fabs(r - forces[i]->distance);
    }

    delete[] forces;
    return anomaly;
}

float InteratomicForce::covalent_bond_radius(Atom* a, Atom* b, float cardinality)
{
    if (!read_forces_dat && !reading_forces) read_all_forces();
    InteratomicForce** retval = new InteratomicForce*[16];
    init_nulls(retval, 16);

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

float InteratomicForce::coordinate_bond_radius(Atom* a, Atom* b, intera_type btype)
{
    if (!read_forces_dat && !reading_forces) read_all_forces();
    InteratomicForce** retval = new InteratomicForce*[16];
    init_nulls(retval, 16);

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

    throw BOND_DEF_NOT_FOUND;
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




