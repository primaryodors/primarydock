
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
    for (i=0; i<MAX_DYN_NEARBY; i++) nearby_contacts[i] = nullptr;
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

float DynamicMotion::get_ligand_contact_energy(Molecule* ligand)
{
    if (!prot) throw 0xffff;

    int sr = prot->get_residue(start_resno)->get_residue_no();
    int er = prot->get_residue(end_resno)->get_residue_no();
    if (!sr || !er) throw 0xffff;

    float result = 0;

    int i, j;
    for (i=sr; i<=er; i++)
    {
        AminoAcid* aa = prot->get_residue(i);
        if (!aa) continue;
    
        float e = -((Molecule*)aa)->get_intermol_binding(ligand);
            
        #if _dbg_soft_dynamics
        // cout << *aa << " ~ ligand energy: " << e << endl;
        #endif

        result += e;
    }

    return result;
}

float DynamicMotion::get_nearby_contact_energy()
{
    if (!prot) throw 0xffff;

    int sr = prot->get_residue(start_resno)->get_residue_no();
    int er = prot->get_residue(end_resno)->get_residue_no();
    if (!sr || !er) throw 0xffff;

    if (!nearby_contacts[0]) fill_nearby_contacts();

    float result = 0;

    int i, j;
    for (i=sr; i<=er; i++)
    {
        AminoAcid* aa = prot->get_residue(i);
        if (!aa) continue;
    
        for (j=0; nearby_contacts[j]; j++)
        {
            if (aa->get_CA_location().get_3d_distance(nearby_contacts[j]->get_CA_location()) > (aa->get_reach() + nearby_contacts[j]->get_reach()) ) continue;

            float e = -aa->get_intermol_binding(nearby_contacts[j]);
            
            #if _dbg_internal_clashes
            if (e > 0) cout << *aa << " clashes with " << *(nearby_contacts[j]) << " by " << e << endl;
            #endif

            result += e;
        }
    }

    return result;
}

void DynamicMotion::fill_nearby_contacts()
{
    if (!prot) throw 0xffff;

    int sr = prot->get_residue(start_resno)->get_residue_no();
    int er = prot->get_residue(end_resno)->get_residue_no();
    if (!sr || !er) throw 0xffff;

    int sr5 = sr-5, er5 = er+5;

    int i, j, l, n;
    n = prot->get_end_resno();
    bool resno_used[n+1];

    for (i=1; i<=n; i++) resno_used[i] = false;

    l=0;
    for (i=sr; i<=er; i++)
    {
        AminoAcid** near = prot->get_residues_can_clash(i);
        if (!near) continue;

        for (j=0; near[j]; j++)
        {
            int resno = near[j]->get_residue_no();
            if (resno_used[resno]) continue;
            if (resno >= sr5 && resno <= er5) continue;

            nearby_contacts[l++] = near[j];
            resno_used[resno] = true;
        }
        nearby_contacts[l] = nullptr;
    }
}

float DynamicMotion::apply_incremental(float amt)
{
    int i;

    float new_total = applied + amt;
    for (i=0; constraints[i]; i++)
    {
        if (constraints[i]->type == ddep_MIN)
        {
            if (new_total < constraints[i]->depends_on->applied) amt = constraints[i]->depends_on->applied - applied;
        }
        else if (constraints[i]->type == ddep_MAX)
        {
            if (new_total > constraints[i]->depends_on->applied) amt = constraints[i]->depends_on->applied - applied;
        }
        else if (constraints[i]->type == ddep_SYNC)
        {
            amt = constraints[i]->depends_on->applied - applied;
        }
    }

    new_total = applied + amt;
    if (new_total < -DYN_MAX_OVERAGE) amt = -DYN_MAX_OVERAGE - applied;
    else if (new_total > 1.0+DYN_MAX_OVERAGE) amt = (1.0+DYN_MAX_OVERAGE) - applied;

    new_total = applied + amt;
    if (new_total < minimum) amt = minimum - applied;

    #if debug_dyn_motion
    if (!strcmp(name.c_str(), "bend7")) cout << "Trying " << amt << " for " << name << " for a total of " << (applied + amt) << endl;
    #endif

    return apply_incremental_nochecks(amt);
}

float DynamicMotion::apply_incremental_nochecks(float amt)
{
    AminoAcid* aa;
    Point fulcrum, ptaxis;
    int i, j, sr, er;
    float lamt, lamt_phi, lamt_psi;
    LocatedVector lv;

    switch (type)
    {
        case dyn_rock:
        
        aa = prot->get_residue(fulcrum_resno);
        if (!aa) throw -1;
        fulcrum = aa->get_CA_location();

        if (axis_resno.helix_no && !axis.r)
        {
            aa = prot->get_residue(axis_resno);
            if (!aa) throw -1;
            ptaxis = aa->get_CA_location();

            axis = ptaxis.subtract(fulcrum);
        }

        i = prot->get_bw50(start_resno.helix_no);
        if (!i) throw -1;
        sr = i - 50 + start_resno.member_no;

        i = prot->get_bw50(end_resno.helix_no);
        if (!i) throw -1;
        er = i - 50 + end_resno.member_no;

        lamt = amt * bias;
        last_change = amt;
        prot->rotate_piece(sr, er, fulcrum, axis, lamt*fiftyseventh);
        applied += amt;
        break;


        case dyn_move:
        
        aa = prot->get_residue(fulcrum_resno);
        if (!aa) throw -1;
        fulcrum = aa->get_CA_location();

        if (axis_resno.helix_no && !axis.r)
        {
            aa = prot->get_residue(axis_resno);
            if (!aa) throw -1;
            ptaxis = aa->get_CA_location();

            axis = ptaxis.subtract(fulcrum);
        }

        i = prot->get_bw50(start_resno.helix_no);
        if (!i) throw -1;
        sr = i - 50 + start_resno.member_no;

        i = prot->get_bw50(end_resno.helix_no);
        if (!i) throw -1;
        er = i - 50 + end_resno.member_no;

        lamt = amt * bias;
        last_change = amt;
        axis.r = lamt;
        prot->move_piece(sr, er, axis);
        applied += amt;
        break;


        case dyn_wind:

        i = prot->get_bw50(start_resno.helix_no);
        if (!i) throw -1;
        sr = i - 50 + start_resno.member_no;

        i = prot->get_bw50(end_resno.helix_no);
        if (!i) throw -1;
        er = i - 50 + end_resno.member_no;

        lamt = amt * bias * fiftyseventh;
        last_change = amt;
        lamt_phi = lamt/(M_PI*2) * (ALPHA_PHI - -M_PI);
        lamt_psi = lamt/(M_PI*2) * (ALPHA_PSI - -M_PI);

        for (i=sr; i<=er; i++)
        {
            prot->rotate_backbone_partial(i, er, N_asc, lamt_phi);
            prot->rotate_backbone_partial(i, er, CA_asc, lamt_psi);
        }
        applied += amt;

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
        else if (!strcmp(words[l], "SYNC")) dc->type = ddep_SYNC;
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

void DynamicMotion::make_random_change()
{
    float f = frand(-DYN_RANDOM_RANGE, DYN_RANDOM_RANGE) + (1.0 - applied) * DYN_BIAS_EFFECT;
    apply_incremental(f);
}

void DynamicMotion::undo()
{
    apply_incremental(-last_change);

    #if debug_dyn_motion
    // cout << "Reverting " << name << endl;
    #endif
}
