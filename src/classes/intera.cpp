
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

void InteratomicForce::read_all_forces()
{
    int i, ifcount = 0;
    init_nulls(all_forces, _MAX_NUM_FORCES);

    FILE* pf = fopen("bindings.dat", "rb");
    if (!pf)
        cout << "ERROR failed to open bindings.dat, please verify file exists and you have permissions." << endl;
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

                if (all_forces[ifcount]->Za < 36 && all_forces[ifcount]->Zb < 36)
                {
                    if (!forces_by_Z[all_forces[ifcount]->Za][all_forces[ifcount]->Zb])
                    {
                        forces_by_Z[all_forces[ifcount]->Za][all_forces[ifcount]->Zb] = new InteratomicForce*[16];
                        for (i=0; i<16; i++)
                            forces_by_Z[all_forces[ifcount]->Za][all_forces[ifcount]->Zb][i] = 0;
                    }

                    for (i=0; i<15; i++)
                    {
                        if (!forces_by_Z[all_forces[ifcount]->Za][all_forces[ifcount]->Zb][i])
                        {
                            forces_by_Z[all_forces[ifcount]->Za][all_forces[ifcount]->Zb][i] = all_forces[ifcount];
                            // if (all_forces[ifcount]->type == mcoord) cout << all_forces[ifcount]->Za << " " << all_forces[ifcount]->Zb << " " << i << " " << all_forces[ifcount]->type << endl;
                            break;
                        }
                    }

                    if (!forces_by_Z[all_forces[ifcount]->Zb][all_forces[ifcount]->Za])
                    {
                        forces_by_Z[all_forces[ifcount]->Zb][all_forces[ifcount]->Za] = new InteratomicForce*[16];
                        for (i=0; i<16; i++)
                            forces_by_Z[all_forces[ifcount]->Zb][all_forces[ifcount]->Za][i] = 0;
                    }

                    for (i=0; i<15; i++)
                    {
                        if (!forces_by_Z[all_forces[ifcount]->Zb][all_forces[ifcount]->Za][i])
                        {
                            forces_by_Z[all_forces[ifcount]->Zb][all_forces[ifcount]->Za][i] = all_forces[ifcount];
                            break;
                        }
                    }
                }

                ifcount++;
            }
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
    char** fields = chop_spaced_fields(line);
    if (fields[0]
            && fields[1]
            && fields[2]
            && fields[3]
            && fields[4]
       )
    {
        char ea[3], eab[3], eb[3], ebb[3];

        if (fields[0][1] >= 'A' && fields[0][1] <= 'Z')
        {
            eab[0] = fields[0][0];
            eab[1] = 0;
            strcpy(ea, fields[0]+1);
        }
        else if (fields[0][1] > 0 && fields[0][1] < 'A')
        {
            if (fields[0][1] == '-') aritybZa = 1;
            if (fields[0][1] == '=') aritybZa = 2;
            if (fields[0][1] == '#') aritybZa = 3;
            if (fields[0][1] == '$') aritybZa = 4;
            eab[0] = fields[0][0];
            eab[1] = 0;
            strcpy(ea, fields[0]+2);
        }
        else if (fields[0][2] >= 'A' && fields[0][2] <= 'Z')
        {
            eab[0] = fields[0][0];
            eab[1] = fields[0][1];
            eab[2] = 0;
            strcpy(ea, fields[0]+2);
        }
        else if (fields[0][2] > 0 && fields[0][2] < 'A')
        {
            if (fields[0][2] == '-') aritybZa = 1;
            if (fields[0][2] == '=') aritybZa = 2;
            if (fields[0][2] == '#') aritybZa = 3;
            if (fields[0][2] == '$') aritybZa = 4;
            eab[0] = fields[0][0];
            eab[1] = fields[0][1];
            eab[2] = 0;
            strcpy(ea, fields[0]+3);
        }
        else
        {
            eab[0] = 0;
            strcpy(ea, fields[0]);
        }


        if (fields[1][1] >= 'A' && fields[1][1] <= 'Z')
        {
            ebb[0] = fields[1][0];
            ebb[1] = 0;
            strcpy(eb, fields[1]+1);
        }
        else if (fields[1][1] > 0 && fields[1][1] < 'A')
        {
            if (fields[1][1] == '-') aritybZb = 1;
            if (fields[1][1] == '=') aritybZb = 2;
            if (fields[1][1] == '#') aritybZb = 3;
            if (fields[1][1] == '$') aritybZb = 4;
            ebb[0] = fields[1][0];
            ebb[1] = 0;
            strcpy(eb, fields[1]+2);
        }
        else if (fields[1][2] >= 'A' && fields[1][2] <= 'Z')
        {
            ebb[0] = fields[1][0];
            ebb[1] = fields[1][1];
            ebb[2] = 0;
            strcpy(eb, fields[1]+2);
        }
        else if (fields[1][2] > 0 && fields[1][2] < 'A')
        {
            if (fields[1][2] == '-') aritybZb = 1;
            if (fields[1][2] == '=') aritybZb = 2;
            if (fields[1][2] == '#') aritybZb = 3;
            if (fields[1][2] == '$') aritybZb = 4;
            ebb[0] = fields[1][0];
            ebb[1] = fields[1][1];
            ebb[2] = 0;
            strcpy(eb, fields[1]+3);
        }
        else
        {
            ebb[0] = 0;
            strcpy(eb, fields[1]);
        }


        Za  = Atom::Z_from_esym(ea);
        bZa = Atom::Z_from_esym(eab);
        Zb  = Atom::Z_from_esym(eb);
        bZb = Atom::Z_from_esym(ebb);			// There is only ebb, no flow. Cope.


        if (!strcmp(fields[2], "coval"))	type = covalent;
        if (!strcmp(fields[2], "ionic"))	type = ionic;
        if (!strcmp(fields[2], "hbond"))	type = hbond;
        if (!strcmp(fields[2], "pi"   ))	type = pi;
        if (!strcmp(fields[2], "plpi" ))	type = polarpi;
        if (!strcmp(fields[2], "coord"))	type = mcoord;
        if (!strcmp(fields[2], "vdW"  ))	type = vdW;

        if (fields[3])
        {
            arity = atof(fields[3]);
            if (fields[4])
            {
                distance = atof(fields[4]);
                if (fields[5])
                {
                    kJ_mol = atof(fields[5]);
                    if (fields[6])
                    {
                        dirprop = atof(fields[6]);
                    }
                    else dirprop = 2;
                }
                else kJ_mol = 200;
            }
            else distance=1;
        }

        // cout << *this << endl;
    }
    delete[] fields;
}

bool InteratomicForce::atom_is_capable_of(Atom* a, intera_type t)
{
    InteratomicForce** look = all_forces;
    int i, Za = a->get_Z();

    for (i=0; look[i]; i++)
    {
        if (look[i]->type == t)
            if (look[i]->Za == Za || look[i]->Zb == Za)
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

InteratomicForce** InteratomicForce::get_applicable(Atom* a, Atom* b)
{
    if (!read_forces_dat && !reading_forces) read_all_forces();
    if (!a || !b)
    {
        return NULL;
    }

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

    for (i=0; look[i]; i++)
    {
        if (	(	look[i]->Za == Za
                    &&
                    (	!look[i]->bZa
                        ||
                        (!look[i]->aritybZa && a->is_bonded_to(Atom::esym_from_Z(look[i]->bZa)))
                        ||
                        ( look[i]->aritybZa && a->is_bonded_to(Atom::esym_from_Z(look[i]->bZa), look[i]->aritybZa))
                    )
                    &&
                    look[i]->Zb == Zb
                    &&
                    (	!look[i]->bZb
                        ||
                        (!look[i]->aritybZb && b->is_bonded_to(Atom::esym_from_Z(look[i]->bZb)))
                        ||
                        ( look[i]->aritybZb && b->is_bonded_to(Atom::esym_from_Z(look[i]->bZb), look[i]->aritybZb))
                    )
             )
                ||
                (	look[i]->Zb == Za
                    &&
                    (	!look[i]->bZb
                        ||
                        (!look[i]->aritybZb && a->is_bonded_to(Atom::esym_from_Z(look[i]->bZb)))
                        ||
                        ( look[i]->aritybZb && a->is_bonded_to(Atom::esym_from_Z(look[i]->bZb), look[i]->aritybZb))
                    )
                    &&
                    look[i]->Za == Zb
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
                if (sgn(a->get_acidbase()) == -sgn(b->get_acidbase())
                        ||
                        sgn(a->get_charge()) == -sgn(b->get_charge())
                        ||
                        ((a == b) && a->get_acidbase())
                   )
                    retval[j++] = look[i];
                break;

            case hbond:
                // if (sgn(a->is_polar()) == -sgn(b->is_polar()))
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

float InteratomicForce::total_binding(Atom* a, Atom* b)
{
    InteratomicForce** forces = get_applicable(a, b);

    int i, j, k;
    float kJmol = 0;

    float r = a->distance_to(b);
    float avdW = a->get_vdW_radius(), bvdW = b->get_vdW_radius();
    float rbind = avdW+bvdW;

    float dp = 2;			// Directional propensity. The anisotropic component from each single vertex
    // is calculated as a cosine and then raised to this exponent.
    Point center;

    if (!forces) goto _canstill_clash;

    if (r < 0.5) forces[0] = NULL;

    for (i=0; forces[i]; i++)
    {
        if (!forces[i]->distance) continue;
        float r1 = r / forces[i]->distance;

        if (forces[i]->type == covalent) continue;

        if (forces[i]->type == vdW
                ||
                forces[i]->type == pi
                ||
                forces[i]->type == polarpi
           )
        {
            if (r1 > 2.5) continue;
        }
        else if (r1 > 7) continue;

        dp = 2;
        if (forces[i]->type == vdW) dp = 1;					// van der Waals forces are non-directional, but the C-H bond still shields.
        else if (forces[i]->get_dp()) dp = forces[i]->get_dp();

        // Anisotropy.
        SCoord* ageo = a->get_geometry_aligned_to_bonds();
        SCoord* bgeo = b->get_geometry_aligned_to_bonds();
        int ag = a->get_geometry();
        int bg = b->get_geometry();
        int abc = a->get_bonded_atoms_count();
        int bbc = b->get_bonded_atoms_count();
        float asum=0, bsum=0, aniso=1;
        bool del_ageo=false, del_bgeo=false;

        if (forces[i]->type == pi && ag >= 3 && bg >= 3)
        {
            ageo = get_geometry_for_pi_stack(ageo);
            bgeo = get_geometry_for_pi_stack(bgeo);
            ag = bg = 5;
            del_ageo = del_bgeo = true;
        }
        else if ((forces[i]->type == polarpi || forces[i]->type == mcoord) && ag >= 3 && bg >= 3)
        {
            if (!a->is_polar())
            {
                ageo = get_geometry_for_pi_stack(ageo);
                ag = 5;
                del_ageo = true;
            }
            if (!b->is_polar())
            {
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
                dpb = 2;		// Assume same for all donors.
            }
            else if (b->is_polar() < 0 && a->is_polar() >= 0)
            {
                dpb = dp;
                dpa = 2;
            }
            else dpa = dpb = dp;
        }

        Point aloc = a->get_location();
        Point bloc = b->get_location();

        int anx = a->get_idx_next_free_geometry();
        int bnx = b->get_idx_next_free_geometry();

        if (forces[i]->type != ionic)
        {
            if (a->get_bonded_atoms_count() == ag || b->get_bonded_atoms_count() == bg)
                continue;

            // Sum up the anisotropic contribution from each geometry vertex of a.
            for (j=anx; j<ag; j++)
            {
                Point pt(&ageo[j]);
                pt.scale(r);
                pt = pt.add(aloc);
                float theta = find_3d_angle(&bloc, &pt, &aloc);
                if (theta > M_PI/2) continue;
                float contrib = pow(fmax(0,cos(theta)), dpa);
                if (!isnan(contrib) && !isinf(contrib)) asum += contrib;
                // else cout << "Bad contrib! " << cos(theta) << " = cos(" << (theta*180.0/M_PI) << ") exp=" << dpa << endl;
                // if (fabs(contrib) > 10000) cout << "Bad contrib! " << cos(theta) << " = cos(" << (theta*180.0/M_PI) << ") exp=" << dpa << endl;
            }

            // Sum up the anisotropic contribution from each geometry vertex of b.
            for (j=bnx; j<bg; j++)
            {
                Point pt(&bgeo[j]);
                pt.scale(r);
                pt = pt.add(bloc);
                float theta = find_3d_angle(&aloc, &pt, &bloc);
                if (theta > M_PI/2) continue;
                float contrib = pow(fmax(0,cos(theta)), dpb);
                if (!isnan(contrib) && !isinf(contrib)) bsum += contrib;
            }

            if (asum > 1) asum = 1;
            if (bsum > 1) bsum = 1;

            // Multiply the two sums.
            aniso = asum * bsum;
        }

        float partial, rdecayed;
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
            rdecayed = r1 * r1;
            // partial = aniso * forces[i]->kJ_mol * rdecayed;
            partial = aniso * forces[i]->kJ_mol;
        }

        /*if (fabs(partial) >= 10000)
        // if (isnan(partial) || isinf(partial))
        {	cout << "Invalid partial! " << a->name << ":" << b->name << " r=" << r
        		 << " (optimal " << forces[i]->distance << ") rdecayed=" << rdecayed
        		 << " aniso=" << aniso << " (" << asum << "*" << bsum << ")" << endl;
        }*/

        if (fabs(partial) > fabs(forces[i]->kJ_mol) || partial >= 500)
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
            partial *= ((a->get_charge()) * -(b->get_charge()));

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

        if (fabs(partial) >= 500)
        {
            cout << "Invalid partial! " << partial << " (max " << forces[i]->kJ_mol << ") from "
                 << a->name << "...." << b->name << " r=" << r
                 << " (optimal " << forces[i]->distance << ") rdecayed=" << rdecayed
                 << " aniso=" << aniso << " (" << asum << "*" << bsum << ")" << endl;
            throw 0xbadf0ace;
        }

        kJmol += partial;
        if (partial > 0.5 && forces[i]->distance < rbind) rbind = forces[i]->distance;

        k = (forces[i]->type - covalent) % _INTER_TYPES_LIMIT;
        total_binding_by_type[k] += partial;

        if (del_ageo) delete[] ageo;
        if (del_bgeo) delete[] bgeo;

        if (forces[i]->type == ionic) break;
    }

    if (rbind < 0.7) rbind = 0.7;

_canstill_clash:
    ;
    if (r < rbind)
    {
        float f = rbind/(avdW+bvdW);
        kJmol -= pow(fabs(sphere_intersection(avdW*f, bvdW*f, r)*_kJmol_cuA), 4);
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




