
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

void DynamicMotion::read_config_line(const char* ln, DynamicMotion** all)
{
    char buffer[strlen(ln)+16];
    strcpy(buffer, ln);
    char** words = chop_spaced_words(buffer);
    if (!words[0]) throw -1;

    if (!words[1]) throw -1;
    if (!strcmp(words[1], "ROCK")) type = dyn_rock;
    else if (!strcmp(words[1], "BEND")) type = dyn_bend;
    else if (!strcmp(words[1], "MOVE")) type = dyn_move;
    else if (!strcmp(words[1], "WIND")) type = dyn_wind;
    else throw -1;

    if (!words[2]) throw -1;
    name = words[2];

    int i, l = 3;
    if (!words[l]) throw -1;
    start_resno.from_string(words[l]);
    if (!words[++l]) throw -1;
    end_resno.from_string(words[l]);

    if (type == dyn_rock || type == dyn_move)
    {
        if (!words[++l]) throw -1;
        fulcrum_resno.from_string(words[l]);

        if (!words[++l]) throw -1;
        axis_resno.from_string(words[l]);
    }

    if (!words[++l]) throw -1;
    bias = atof(words[l]);

    while (words[++l])
    {
        DynamicConstraint* dc = new DynamicConstraint;
        if (!strcmp(words[l], "MIN")) dc->type = ddep_MIN;
        else if (!strcmp(words[l], "MAX")) dc->type = ddep_MAX;
        else throw -1;

        if (!words[++l]) throw -1;
        if (!all) throw -1;
        bool matched = false;
        for (i=0; all[i]; i++)
        {
            if (!strcmp(all[i]->name.c_str(), words[l]))
            {
                dc->depends_on = all[i];
                matched = true;
                break;
            }
        }

        if (!matched)
        {
            cout << "Motion " << words[l] << " not found for constraint." << endl;
            throw -1;
        }

        add_constraint(dc);
    }
    
    return;
}