
#include "dynamic.h"

using namespace std;

DynamicMotion::DynamicMotion(Protein* ppro)
{
    if (!ppro)
    {
        cout << "Null protein passed into DynamicMotion constructor." << endl;
        throw -1;
    }

    prot = ppro;
    int i;
    for (i=0; i<=MAX_DYN_CONSTRAINTS; i++) constraints[i] = nullptr;
}

bool DynamicMotion::add_constraint(DynamicConstraint* c)
{
    if (!c->depends_on)
    {
        cout << "Null dependency in constraint added to dynamic motion." << endl;
        throw -1;
    }

    int i;
    for (i=0; constraints[i]; i++);                     // Count existing constraints.
    if (i >= MAX_DYN_CONSTRAINTS) return false;         // Array full.
    constraints[i] = c;
    return true;
}

float DynamicMotion::apply_incremental(float amt)
{
    AminoAcid* aa;
    Point fulcrum, ptaxis;
    SCoord axis;
    int i, j, sr, er;
    float lamt, lamt_phi, lamt_psi;
    LocatedVector lv;

    switch (type)
    {
        case dyn_rock:
        
        aa = prot->get_residue(fulcrum_resno);
        if (!aa) throw -1;
        fulcrum = aa->get_CA_location();

        aa = prot->get_residue(axis_resno);
        if (!aa) throw -1;
        ptaxis = aa->get_CA_location();

        axis = ptaxis.subtract(fulcrum);

        i = prot->get_bw50(start_resno.helix_no);
        if (!i) throw -1;
        sr = i - 50 + start_resno.member_no;

        i = prot->get_bw50(end_resno.helix_no);
        if (!i) throw -1;
        er = i - 50 + end_resno.member_no;

        lamt = amt * bias;
        prot->rotate_piece(sr, er, fulcrum, axis, lamt);
        applied += lamt;
        break;


        case dyn_move:
        
        aa = prot->get_residue(fulcrum_resno);
        if (!aa) throw -1;
        fulcrum = aa->get_CA_location();

        aa = prot->get_residue(axis_resno);
        if (!aa) throw -1;
        ptaxis = aa->get_CA_location();

        axis = ptaxis.subtract(fulcrum);

        i = prot->get_bw50(start_resno.helix_no);
        if (!i) throw -1;
        sr = i - 50 + start_resno.member_no;

        i = prot->get_bw50(end_resno.helix_no);
        if (!i) throw -1;
        er = i - 50 + end_resno.member_no;

        lamt = amt * bias;
        axis.r = lamt;
        prot->move_piece(sr, er, axis);
        applied += lamt;
        break;


        case dyn_wind:

        i = prot->get_bw50(start_resno.helix_no);
        if (!i) throw -1;
        sr = i - 50 + start_resno.member_no;

        i = prot->get_bw50(end_resno.helix_no);
        if (!i) throw -1;
        er = i - 50 + end_resno.member_no;

        lamt = amt * bias;
        lamt_phi = lamt/(M_PI*2) * (ALPHA_PHI - -M_PI);
        lamt_psi = lamt/(M_PI*2) * (ALPHA_PSI - -M_PI);

        for (i=sr; i<=er; i++)
        {
            aa = prot->get_residue(i);

            lv = aa->rotate_backbone(N_asc, lamt_phi);
            for (j=i+1; j<=er; j++)
            {
                aa = prot->get_residue(j);
                aa->rotate(lv, lamt_phi);
            }

            lv = aa->rotate_backbone(CA_asc, lamt_psi);
            for (j=i+1; j<=er; j++)
            {
                aa = prot->get_residue(j);
                aa->rotate(lv, lamt_psi);
            }
        }

        break;


        default:
        throw -1;
    }

    return applied;
}

float DynamicMotion::apply_absolute(float amt)
{
    return apply_incremental(amt - applied);
}