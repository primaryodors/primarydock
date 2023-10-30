#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <regex>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <algorithm>
#include <unistd.h>
#include "classes/dynamic.h"
#include "classes/group.h"
#include "classes/protein.h"
#include "classes/scoring.h"

using namespace std;

struct AcvBndRot
{
    int resno;
    std::string aname;
    std::string bname;
    Atom* atom = nullptr;
    Atom* btom = nullptr;
    Bond* bond = nullptr;
    float theta;
};



char* get_file_ext(char* filename)
{
    int i = strlen(filename)-1;
    for (; i>=0; i--)
    {
        if (filename[i] == '.') return &filename[i+1];
    }
    return 0;
}

bool output_each_iter = false;
bool progressbar = false;
int movie_offset = 0;
char configfname[256];
char protfname[256];
char protstrand = '\0';
char protafname[256];
char tplfname[256];
char tplrfnam[256];
char ligfname[256];
char smiles[256];
char outfname[256];
Point pocketcen, loneliest, pocketsize, ligbbox;
std::ofstream *output = NULL;

std::vector<int> exclusion;

std::string CEN_buf = "";
std::vector<std::string> pathstrs;
std::vector<std::string> states;

std::vector<std::string> dyn_strings;
std::vector<DynamicMotion> dyn_motions;

bool configset=false, protset=false, tplset=false, tprfset=false, ligset=false, ligcmd=false, smset = false, smcmd = false, pktset=false;

Protein* protein;
Protein* ptemplt;
Protein* ptplref;
int seql = 0;
int mcoord_resno[256];
int addl_resno[256];
Molecule* ligand;
Molecule** waters = nullptr;
Molecule** owaters = nullptr;
Point ligcen_target;
Point size(10,10,10);
SCoord path[256];
int pathnodes = 0;				// The pocketcen is the initial node.
int poses = 10;
int iters = 50;
int maxh2o = 0;
int omaxh2o = 0;
bool flex = true;
float kJmol_cutoff = 0.01;
int pose, nodeno, iter;
bool kcal = false;
float drift = initial_drift;
Molecule** gcfmols = NULL;
int activation_node = -1;		// Default is never.
int found_poses = 0;
int triesleft = 0;				// Default is no retry.
bool echo_progress = false;
bool hydrogenate_pdb = false;
std::string temp_pdb_file = "";
int pid = getpid();
bool append_pdb = false;
bool do_output_colors = false;

AminoAcid*** reaches_spheroid = nullptr;
int sphres = 0;

std::string origbuff = "";
std::string optsecho = "";

// Switch to enable "best-binding" algorithm rather than "tumble spheres" algorithm.
bool use_bestbind_algorithm = default_bestbind;
bool use_prealign = false;
std::string prealign_residues = "";
Bond retain_bindings[4];

#if _use_groups
AtomGroup ligand_groups[3];
ResidueGroup sc_groups[3];
#else
Atom** ligbb = nullptr;
Atom** ligbbh = nullptr;

intera_type lig_inter_typ[5];
Molecule* alignment_aa[5];
#endif

Pose pullaway_undo;
float last_ttl_bb_dist = 0;

std::vector<AcvBndRot> active_bond_rots;
std::vector<int> tripswitch_clashables;
std::vector<ResiduePlaceholder> required_contacts;
std::vector<std::string> bridges;
std::vector<std::string> atomto;
std::string outpdb;
int outpdb_poses = 0;

std::vector<int>flexible_resnos;
std::vector<ResiduePlaceholder>forced_flexible_resnos;
std::vector<ResiduePlaceholder>forced_static_resnos;

#if _dummy_atoms_for_debug
std::vector<Atom> dummies;
#endif


void append_dummy(Point pt)
{
    #if _dummy_atoms_for_debug
    Atom a("He");
    a.move(pt);

    a.name = new char[8];
    int i=dummies.size()+1;
    sprintf(a.name, "HE%i", i);

    strcpy(a.aa3let, "DMY");

    dummies.push_back(a);
    #endif
}

void delete_water(Molecule* mol)
{
    if (!waters) return;
    int i, j;
    for (i=0; waters[i]; i++)
    {
        if (waters[i] == mol)
        {
            for (j=i+1; waters[j]; j++)
            {
                waters[j-1] = waters[j];
            }
            waters[j-1] = nullptr;
            maxh2o--;

            #if _DBG_H2O_TELEPORT
            cout << "Deleted water molecule " << mol << endl;
            #endif
        }
        break;
    }
}

float teleport_water(Molecule* mol)
{
    if (!waters) return -1000;

    int i, j;
    float e;
    for (j=0; waters[j]; j++)
        if (waters[j] == mol) break;
    
    if (!waters[j]) return -1000;

    for (i=0; i<_water_teleport_tries; i++)
    {
        Point teleport(
            ligcen_target.x + frand(-size.x, size.x),
            ligcen_target.y + frand(-size.y, size.y),
            ligcen_target.z + frand(-size.z, size.z)
                      );
        waters[j]->recenter(teleport);
        e = -waters[j]->get_intermol_binding(gcfmols);
        if (e < _water_satisfaction_threshold) break;
    }
    if (e > _water_satisfaction_threshold) delete_water(mol);
    #if _DBG_H2O_TELEPORT
    else cout << "Teleported water molecule " << mol << endl;
    #endif

    return e;
}

MCoord* search_mtlcoords_for_residue(AminoAcid* aa)
{
    int m, n;
    if (!(n = mtlcoords.size())) return nullptr;

    int i, j;
    for (i=0; i<n; i++)
    {
        m = mtlcoords[i].coordres.size();
        for (j=0; j<m; j++)
        {
            if (mtlcoords[i].coordres[j].resno == aa->is_residue()) return &mtlcoords[i];
        }
    }

    return nullptr;
}

int interpret_resno(const char* field)
{
    char buffer[strlen(field)+4];
    strcpy(buffer, field);
    char* dot = strchr(buffer, '.');
    if (dot)
    {
        *(dot++) = 0;
        int b = atoi(buffer);
        int w = atoi(dot);
        int _50 = protein->get_bw50(b);
        if (_50 < 1)
        {
            cout << "Error: unknown BW number " << b << "." << w << ", please ensure PDB file has REMARK 800 SITE BW words." << endl;
            throw 0xbad12e5;
        }
        return _50 + w - 50;
    }
    else return atoi(buffer);
}

void freeze_bridged_residues()
{
    int i, l;

    if (bridges.size())
    {
        for (i=0; i<bridges.size(); i++)
        {
            int resno1 = interpret_resno(bridges[i].c_str());
            const char* r2 = strchr(bridges[i].c_str(), '|');
            if (!r2) throw 0xbadc0de;
            r2++;
            int resno2 = interpret_resno(r2);
            
            AminoAcid *aa1 = protein->get_residue(resno1), *aa2 = protein->get_residue(resno2);
            if (aa1)
            {
                aa1->movability = MOV_PINNED;
                aa1->been_flexed = true;
                Bond** bb = aa1->get_rotatable_bonds();
                if (bb)
                {
                    for (l=0; bb[l]; l++)
                    {
                        bb[l]->can_rotate = false;
                    }
                    // delete bb;
                }
            }
            if (aa2)
            {
                aa2->movability = MOV_PINNED;
                aa2->been_flexed = true;
                Bond** bb = aa2->get_rotatable_bonds();
                if (bb)
                {
                    for (l=0; bb[l]; l++)
                    {
                        bb[l]->can_rotate = false;
                    }
                    // delete bb;
                }
            }
        }
    }

    if (forced_static_resnos.size())
    {
        for (i=0; i<forced_static_resnos.size(); i++)
        {
            forced_static_resnos[i].resolve_resno(protein);
            int resno = forced_static_resnos[i].resno;
            if (!resno) continue;
            
            AminoAcid *aa = protein->get_residue(resno);
            if (!aa) continue;

            aa->movability = MOV_PINNED;
            aa->been_flexed = true;
            Bond** bb = aa->get_rotatable_bonds();
            if (bb)
            {
                for (l=0; bb[l]; l++)
                {
                    bb[l]->can_rotate = false;
                }
                // delete bb;
            }
        }
    }
}

void reconnect_bridges()
{
    int i;
    for (i=0; i<bridges.size(); i++)
    {
        int resno1 = interpret_resno(bridges[i].c_str());
        const char* r2 = strchr(bridges[i].c_str(), '|');
        if (!r2) throw 0xbadc0de;
        r2++;
        int resno2 = interpret_resno(r2);

        #if _dbg_bridges
        cout << "Bridging " << resno1 << " and " << resno2 << "..." << endl;
        #endif

        protein->bridge(resno1, resno2);

        AminoAcid *aa1 = protein->get_residue(resno1), *aa2 = protein->get_residue(resno2);
        if (aa1) aa1->movability = MOV_PINNED;
        if (aa2) aa2->movability = MOV_PINNED; 

        #if _dbg_bridges
        if (!aa1) cout << resno1 << " not found." << endl;
        if (!aa2) cout << resno2 << " not found." << endl;
        if (aa1 && aa2)
        {
            float tb = -aa1->get_intermol_binding(aa2);
            cout << "Bridge energy " << tb << " kJ/mol." << endl;
        }
        #endif
    }

    freeze_bridged_residues();
}

Pose iter_best_pose[1000];
float iter_best_bind;
void iteration_callback(int iter, Molecule** mols)
{
    #if _dbg_linear_motion
    if (ligand->movability & MOV_FORBIDDEN) throw 0xbadbad;
    #endif

    // if (kJmol_cutoff > 0 && ligand->lastbind >= kJmol_cutoff) iter = (iters-1);
    int l;
    float f = 0;

    if (iter == 1)
    {
        for (l=0; mols[l]; l++) iter_best_pose[l].copy_state(mols[l]);
        iter_best_bind = 0;
    }

    for (l=0; mols[l]; l++)
    {
        float lf = mols[l]->get_intermol_binding(mols), ptnl = mols[l]->get_intermol_potential(mols);
        f += lf;

        if (flex
            && iter < 5
            && (lf < 10 || lf < 0.1 * ptnl)
            && mols[l]->movability == MOV_FLXDESEL
            && mols[l]->is_residue()
            && frand(0,1) < 0.05
            )
        {
            mols[l]->movability = MOV_FORCEFLEX;
        }
    }

    #if _dbg_linear_motion
    if (ligand->movability & MOV_FORBIDDEN) throw 0xbadbad;
    #endif

    if (f > iter_best_bind)
    {
        for (l=0; mols[l]; l++) iter_best_pose[l].copy_state(mols[l]);
        iter_best_bind = f;
    }

    if (iter == iters && iter_best_bind > 0)
    {
        for (l=0; mols[l]; l++) iter_best_pose[l].restore_state(mols[l]);
    }
    
    #if _dbg_linear_motion
    if (ligand->movability & MOV_FORBIDDEN) throw 0xbadbad;
    #endif
    
    int i, j, n;

    #if bb_realign_iters
    #if _dbg_bb_realign
    cout << ligand->lastbind << (iter == iters ? " ] " : " ") << flush;
    #endif
    if (use_bestbind_algorithm && global_pairs.size() >= 2 && ligand->lastbind < bb_realign_b_threshold)
    {
        Point scg0 = global_pairs[0]->scg->get_center();
        Point scg1 = global_pairs[1]->scg->get_center();
        Point ag0  = global_pairs[0]->ag->get_center();
        Point ag1  = global_pairs[1]->ag->get_center();

        Rotation rot;

        if (scg0.get_3d_distance(ag0) > bb_realign_threshold_distance)
        {
            rot = align_points_3d(ag0, scg0, ag1);
            if (rot.a < 0)
            {
                rot.a *= -1;
                rot.v.r *= -1;
                rot.v = (SCoord)((Point)rot.v);
            }
            if (rot.a > bb_realign_threshold_angle)
            {
                rot.a = fmin(hexagonal, rot.a*bb_realign_amount) / global_pairs[0]->ag->atoms.size();
                LocatedVector lv = rot.v;
                lv.origin = ag1;
                ligand->rotate(lv, rot.a);
            }
        }

        if (scg1.get_3d_distance(ag1) > bb_realign_threshold_distance)
        {
            rot = align_points_3d(ag1, scg1, ag0);
            if (rot.a < 0)
            {
                rot.a *= -1;
                rot.v.r *= -1;
                rot.v = (SCoord)((Point)rot.v);
            }
            if (rot.a > bb_realign_threshold_angle)
            {
                rot.a = fmin(hexagonal, rot.a*bb_realign_amount) / global_pairs[1]->ag->atoms.size();
                LocatedVector lv = rot.v;
                lv.origin = ag0;
                ligand->rotate(lv, rot.a);
            }
        }

        if (global_pairs.size() >= 3)
        {
            Point scg2 = global_pairs[2]->scg->get_center();
            Point ag2  = global_pairs[2]->ag->get_center();

            if (scg2.get_3d_distance(ag2) > bb_realign_threshold_distance)
            {
                float theta = find_angle_along_vector(ag2, scg2, ag0, ag1.subtract(ag0));

                if (fabs(theta) > bb_realign_threshold_angle)
                {
                    LocatedVector v = (SCoord)ag1.subtract(ag0);
                    v.origin = ag0;

                    ligand->rotate(v, theta * bb_realign_amount / global_pairs[2]->ag->atoms.size());
                }
            }
        }
    }
    #endif
    
    #if _dbg_linear_motion
    if (ligand->movability & MOV_FORBIDDEN) throw 0xbadbad;
    #endif

    Point bary = ligand->get_barycenter();

    int ac = ligand->get_atom_count();
    float bbest = 0;
    Atom *atom, *btom;

    float progress = (float)iter / iters;
    // float lsrca = (1.0 - progress) * soft_rock_clash_allowance;
    float lsrca = (progress < 0.5) ? soft_rock_clash_allowance : 0;

    n = dyn_motions.size();
    if (!(iter % soft_dynamics_every_n_iters) && n)
    {
        #if _dbg_soft_dynamics
        cout << endl;
        #endif

        // Each dynamic motion, try successively smaller increments/decrements, realigning BB pairs each time, until optimal energy.
        // Include contacts between dynamic motion region and nearby residues in the energy calculation.
        for (l=0; l<n; l++)
        {
            float incr = 0.1;
            float before = dyn_motions[l].get_nearby_contact_energy() + dyn_motions[l].get_ligand_contact_energy(ligand);
            float after;
            #if _dbg_soft_dynamics
            cout << dyn_motions[l].fulcrum_resno.helix_no << "." << dyn_motions[l].fulcrum_resno.member_no
                << " starting contact energy " << before << endl;
            #endif

            for (i=0; i<25; i++)
            {
                Pose ligand_was;
                ligand_was.copy_state(ligand);
                #if _dbg_soft_dynamics
                cout << "Copied states." << endl << flush;
                #endif

                dyn_motions[l].apply_incremental(incr);
                #if _dbg_soft_dynamics
                cout << "Applied motion." << endl << flush;
                #endif

                if (use_bestbind_algorithm && global_pairs.size() >= 2)
                    GroupPair::align_groups(ligand, global_pairs, false, 0.3);
                #if _dbg_soft_dynamics
                cout << "Aligned groups." << endl << flush;
                #endif

                after = dyn_motions[l].get_nearby_contact_energy() + dyn_motions[l].get_ligand_contact_energy(ligand);
                #if _dbg_soft_dynamics
                cout << "New contact energy " << after << flush;
                #endif
                

                if (after > before)
                {
                    dyn_motions[l].apply_incremental(-incr);
                    ligand_was.restore_state(ligand);
                    incr *= -0.666;

                    #if _dbg_soft_dynamics
                    cout << " reverting.";
                    #endif
                }
                else
                {
                    before = after;

                    #if _dbg_soft_dynamics
                    cout << " total applied " << dyn_motions[l].get_total_applied();
                    #endif
                }

                #if _dbg_soft_dynamics
                cout << endl << flush;
                #endif

                if (fabs(incr) < 0.01) break;
            }
        }

        #if _dbg_soft_dynamics
        cout << endl;
        #endif
    }
    
    #if _dbg_linear_motion
    if (ligand->movability & MOV_FORBIDDEN) throw 0xbadbad;
    #endif

    if (!iter) goto _oei;
    if (iter == (iters-1)) goto _oei;
    
    #if enforce_no_bb_pullaway
    #if _use_groups
    if (use_bestbind_algorithm && ligand_groups[0].atoms.size())
    #else
    if (ligbb)
    #endif
    {
        float ttl_bb_dist = 0;
        for (l=0; l<3; l++)
        {
            #if _use_groups
            if (global_pairs.size() > l)
            #else
            if (ligbb[l] && alignment_aa[l])
            #endif
            {
                #if _use_groups
                float r = global_pairs[l]->ag->distance_to(global_pairs[l]->scg->get_center()); //  ligand_groups[l].distance_to(sc_groups[l].get_center());
                float r1 = r;
                if (r < 2.5) r = 2.5;
                if (r > _INTERA_R_CUTOFF) r = _INTERA_R_CUTOFF;
                ttl_bb_dist += r * (1.0 + 1.0 / (l+1));
                #if _dbg_bb_pullaway
                cout << pose << ":" << iter << ": Ligand atoms ";
                for (i=0; i<global_pairs[l]->ag->atoms.size(); i++) cout << global_pairs[l]->ag->atoms[i]->name << " ";
                cout << "are " << r1 << " A from residues";
                for (i=0; i<global_pairs[l]->scg->aminos.size(); i++) cout << " " << global_pairs[l]->scg->aminos[i]->get_3letter() << global_pairs[l]->scg->aminos[i]->get_residue_no();
                cout << "." << endl;
                #endif
                #else
                Atom** mbb = alignment_aa[l]->get_most_bindable(1, ligbb[l]);
                Atom* alca = mbb[0];
                delete mbb;             // Delete the pointer array, but not the pointers.

                if (alca)
                {
                    float r = alca->get_location().get_3d_distance(ligbb[l]->get_location());
                    if (r < 2.5) r = 2.5;
                    ttl_bb_dist += r;
                    #if _dbg_bb_pullaway
                    cout << alignment_aa[l]->get_name() << ":" << alca->name << " " << r << "A from " << ligbb[l]->name << "... ";
                    #endif
                }
                #endif
            }
        }
        #if _dbg_bb_pullaway
        cout << endl;
        #endif

        if (!iter || ttl_bb_dist <= last_ttl_bb_dist)
        {
            pullaway_undo.copy_state(ligand);
            last_ttl_bb_dist = ttl_bb_dist;
        }
        else if (iter && ttl_bb_dist > (1.0+bb_pullaway_allowance)*last_ttl_bb_dist)
        {
            pullaway_undo.restore_state(ligand);
        }
    }
    #endif
    
    #if _dbg_linear_motion
    if (ligand->movability & MOV_FORBIDDEN) throw 0xbadbad;
    #endif
    
    bary = ligand->get_barycenter();

    for (i=0; i < ac; i++)
    {
        if (ligand->get_atom(i)->strongest_bind_energy > bbest)
        {
            atom = ligand->get_atom(i);
            bbest = atom->strongest_bind_energy;
            btom = atom->strongest_bind_atom;
        }
    }
    
    #if _dbg_linear_motion
    if (ligand->movability & MOV_FORBIDDEN) throw 0xbadbad;
    #endif

    if (bbest >= 15)
    {
        Point A = btom->get_location();
        Point B = bary;
        Point C = ligcen_target;

        SCoord N = compute_normal(A, B, C);

        float theta = find_3d_angle(B, C, A);
        theta /= 3;
        LocatedVector lv = N;
        lv.origin = A;
        ligand->rotate(lv, theta);

        bary = ligand->get_barycenter();
    }
    else
    {
        #if allow_drift
        #if !pocketcen_is_loneliest
        if (ligand->lastbind <= -100)
        {
            ligcen_target.x += (loneliest.x - ligcen_target.x) * drift;
            ligcen_target.y += (loneliest.y - ligcen_target.y) * drift;
            ligcen_target.z += (loneliest.z - ligcen_target.z) * drift;
        }
        #endif

        if (bary.get_3d_distance(ligcen_target) > size.magnitude())
        {
            //cout << "Wrangle! " << bary << ": " << bary.get_3d_distance(ligcen_target) << " vs. " << size.magnitude() << endl;
            bary = ligcen_target;
            // ligand->reset_conformer_momenta();
        }
        else
        {
            if (ligand->lastbind < 0)
            {
                bary.x += (ligcen_target.x - bary.x) * drift;
                bary.y += (ligcen_target.y - bary.y) * drift;
                bary.z += (ligcen_target.z - bary.z) * drift;
            }
            else drift *= (1.0 - drift_decay_rate/iters);
        }

        ligand->recenter(bary);

        #endif
    }
    
    #if _dbg_linear_motion
    if (ligand->movability & MOV_FORBIDDEN) throw 0xbadbad;
    #endif

    #if _teleport_dissatisfied_waters
    if (waters && (iter % 5) == 4)
    {
        for (i=0; waters[i]; i++)
        {
            float r = waters[i]->get_barycenter().get_3d_distance(ligcen_target);
            if (r > size.magnitude()) teleport_water(waters[i]);
            else
            {
                float e = 0;
                int j;
                for (j=0; j<10; j++)
                {
                    if (!j || waters[i]->lastbind_history[j] > e) e = waters[i]->lastbind_history[j];
                }
                if (e < _water_satisfaction_threshold) teleport_water(waters[i]);
            }
        }
    }
    #endif
    
    #if _dbg_linear_motion
    if (ligand->movability & MOV_FORBIDDEN) throw 0xbadbad;
    #endif

    if (gcfmols && seql)
    {
        Star discrete[SPHREACH_MAX+8];
        /*discrete[0].pmol = gcfmols[0];
        discrete[1].pmol = gcfmols[1];*/

        int offset = 0;

        #if _dbg_anemia
        cout << iter << ": ";
        #endif
        for (i=0; gcfmols[i]; i++)
        {
            if (!gcfmols[i]->is_residue() || !strcmp(gcfmols[i]->get_name(), "MTL") )
            {
                discrete[offset++].pmol = gcfmols[i];
                #if _dbg_anemia
                cout << discrete[offset-1].pmol->get_name() << " ";
                #endif
            }
        }
        #if _dbg_anemia
        cout << endl;
        progressbar = false;
        #endif

        AminoAcid* resphres[SPHREACH_MAX+4];
        for (i=0; i<SPHREACH_MAX+4; i++) resphres[i] = nullptr;
        sphres = protein->get_residues_can_clash_ligand(resphres, ligand, bary, size, addl_resno);
        //cout << "Sphres: " << sphres << endl;
        for (i=0; i<sphres; i++)
        {
            discrete[i+offset].paa = reaches_spheroid[nodeno][i];
        }

        discrete[i+offset].n = 0;

        sphres = i;
        for (i=0; discrete[i].n; i++) gcfmols[i] = discrete[i].pmol;
        gcfmols[i] = nullptr;
    }
    
    #if _dbg_linear_motion
    if (ligand->movability & MOV_FORBIDDEN) throw 0xbadbad;
    #endif

    _oei:
    ;

    #if recapture_ejected_ligand
    Point lig_center = ligand->get_barycenter();
    float r = lig_center.get_3d_distance(ligcen_target);
    float recapture_distance = size.magnitude() / 2;
    if (r >= recapture_distance) ligand->recenter(ligcen_target);
    #endif
    
    #if _dbg_linear_motion
    if (ligand->movability & MOV_FORBIDDEN) throw 0xbadbad;
    #endif

    if (output_each_iter)
    {
        std::string itersfname = (std::string)"tmp/" /*+ (std::string)protein->get_name()*/ + (std::string)"_iters.dock";
        int liter = iter - 1 + movie_offset;
        FILE* fp = fopen(itersfname.c_str(), ((liter == 0 && pose == 1) ? "wb" : "ab") );
        if (fp)
        {
            if (!liter && (pose == 1))
            {
                fprintf(fp, "PDB file: %s\n", protfname);
            }
            fprintf(fp, "Pose: %d\nNode: %d\n\nPDBDAT:\n", pose, liter);
            int foff = 0;

            for (i=0; reaches_spheroid[nodeno][i]; i++)
            {
                reaches_spheroid[nodeno][i]->save_pdb(fp, foff);
                foff += reaches_spheroid[nodeno][i]->get_atom_count();
            }

            ligand->save_pdb(fp, foff);

            fclose(fp);
        }
    }
    
    #if _dbg_linear_motion
    if (ligand->movability & MOV_FORBIDDEN) throw 0xbadbad;
    #endif

    if (progressbar)
    {
        int ni = (pathnodes+1) * iters, pni = poses * ni;
        float percentage = (float)((pose-1) * ni + nodeno * iters + iter) / pni * 100;

        cout << "\033[A|";
        for (i=0; i<80; i++)
        {
            float cmpi = 1.25*i;
            if (cmpi <= percentage) cout << "\u2593";
            else cout << "\u2591";
        }
        i = iter % 4;
        cout << ("|/-\\")[i] << " " << (int)percentage << "%.               " << endl;
    }
    
    #if _dbg_linear_motion
    if (ligand->movability & MOV_FORBIDDEN) throw 0xbadbad;
    #endif
}

Point pocketcen_from_config_words(char** words, Point* old_pocketcen)
{
    int i=1;
    Point local_pocketcen;
    if (!strcmp(words[i], "RES"))
    {
        i++;
        std::vector<int> resnos;
        for (; words[i]; i++)
        {
            int j = interpret_resno(words[i]);
            if (!j) break;
            resnos.push_back(j);
            AminoAcid* aa = protein->get_residue(j);
            aa->priority = true;
        }

        int sz = resnos.size();
        Point foravg[sz + 2];
        for (i=0; i<sz; i++)
        {
            foravg[i] = protein->get_atom_location(resnos[i], "CA");
        }

        return average_of_points(foravg, sz);
    }
    else if (!strcmp(words[i], "REL"))
    {
        if (!old_pocketcen)
        {
            cout << "Error: relative coordinates not supported for CEN." << endl << flush;
            throw 0xbadb19d;
        }
        else
        {
            i++;
            local_pocketcen.x = old_pocketcen->x + atof(words[i++]);
            local_pocketcen.y = old_pocketcen->y + atof(words[i++]);
            local_pocketcen.z = old_pocketcen->z + atof(words[i++]);
            return local_pocketcen;
        }
    }
    else
    {
        if (!strcmp(words[i], "ABS")) i++;
        local_pocketcen.x = atof(words[i++]);
        local_pocketcen.y = atof(words[i++]);
        local_pocketcen.z = atof(words[i++]);
        return local_pocketcen;
    }
}

int interpret_config_line(char** words)
{
    int i;

    optsecho = "";

    if (0) { ; }
    else if (!strcmp(words[0], "ACVBROT"))
    {
        AcvBndRot abr;
        abr.resno = atoi(words[1]);
        abr.aname = words[2];
        abr.bname = words[3];
        abr.theta = atof(words[4]) * fiftyseventh;
        active_bond_rots.push_back(abr);
        optsecho = (std::string)"Active bond rotation " + (std::string)words[2] + (std::string)"-" + (std::string)words[3];
    }
    else if (!strcmp(words[0], "APPENDPROT") || !strcmp(words[0], "OPEND"))
    {
        append_pdb = true;
    }
    else if (!strcmp(words[0], "ATOMTO"))
    {
        atomto.push_back(origbuff);
    }
    else if (!strcmp(words[0], "BRIDGE"))
    {
        std::string str = words[1];
        str += (std::string)"|" + (std::string)words[2];
        bridges.push_back(str);
        return 2;
    }
    else if (!strcmp(words[0], "CEN"))
    {
        CEN_buf = origbuff;
        optsecho = (std::string)"Center " + CEN_buf;
        return 0;
    }
    else if (!strcmp(words[0], "COLORS"))
    {
        do_output_colors = true;
        return 0;
    }
    else if (!strcmp(words[0], "COLORLESS"))
    {
        do_output_colors = false;
        return 0;
    }
    else if (!strcmp(words[0], "DEACVNODE"))
    {
        cout << "Notice: DEACVNODE has been deprecated in favor of the HXR option. Please update your config files." << endl;
        return 0;
    }
    else if (!strcmp(words[0], "DEBUG"))
    {
        if (!words[1])
        {
            cout << "Missing debug file name; check config file." << endl << flush;
            throw 0xbadf12e;
        }
        #if _DBG_STEPBYSTEP
        cout << "Starting a debug outstream." << endl;
        #endif
        debug = new std::ofstream(words[1], std::ofstream::out);
        optsecho = "Debug file: " + (std::string)words[1];
        return 1;
    }
    else if (!strcmp(words[0], "DIFF"))
    {
        differential_dock = true;
        optsecho = "Differential dock.";
    }
    else if (!strcmp(words[0], "ECHO"))
    {
        echo_progress = true;
        optsecho = "Echo on.";
    }
    else if (!strcmp(words[0], "ELIM"))
    {
        kJmol_cutoff = -atof(words[1]);
        optsecho = "Energy limit: " + to_string(-kJmol_cutoff);
        return 1;
    }
    else if (!strcmp(words[0], "EMIN"))
    {
        kJmol_cutoff = atof(words[1]);
        optsecho = "Energy limit: " + to_string(kJmol_cutoff);
        return 1;
    }
    else if (!strcmp(words[0], "EXCL"))
    {
        i=1;
        int excls = atoi(words[i++]);
        int excle = atoi(words[i++]);

        for (i=excls; i<=excle; i++) exclusion.push_back(i);
        optsecho = "Exclude range " + to_string(excls) + (std::string)"-" + to_string(excle);
        return i-1;
    }
    else if (!strcmp(words[0], "FLEX"))
    {
        flex = (atoi(words[1]) != 0);
        optsecho = "Flex: " + (std::string)(flex ? "ON" : "OFF");
        return 1;
    }
    else if (!strcmp(words[0], "FLXR"))
    {
        i = 1;
        while (words[i])
        {
            ResiduePlaceholder rph;
            rph.set(words[i]);
            forced_flexible_resnos.push_back(rph);
            i++;
        }
        return i-1;
    }
    else if (!strcmp(words[0], "STCR"))
    {
        i = 1;
        while (words[i])
        {
            ResiduePlaceholder rph;
            rph.set(words[i]);
            forced_static_resnos.push_back(rph);
            i++;
        }
        return i-1;
    }
    else if (!strcmp(words[0], "H2O"))
    {
        maxh2o = omaxh2o = atoi(words[1]);
        if (maxh2o > 0)
        {
            waters = new Molecule*[maxh2o+2];
            owaters = new Molecule*[maxh2o+2];
            for (i=0; i<maxh2o; i++)
            {
                waters[i] = new Molecule("H2O");
                waters[i]->from_smiles("O");
                owaters[i] = waters[i];

                int j;
                for (j=0; j<3; j++) strcpy(waters[i]->get_atom(j)->aa3let, "H2O");
            }
            waters[i] = nullptr;
            owaters[i] = nullptr;
        }
        optsecho = "Water molecules: " + to_string(maxh2o);
        return 1;
    }
    else if (!strcmp(words[0], "HYDRO"))
    {
        hydrogenate_pdb = true;
    }
    else if (!strcmp(words[0], "ITER") || !strcmp(words[0], "ITERS"))
    {
        iters = atoi(words[1]);
        optsecho = "Iterations: " + to_string(iters);
        return 1;
    }
    else if (!strcmp(words[0], "KCAL"))
    {
        kcal = true;
        optsecho = "Output units switched to kcal/mol.";
        return 0;
    }
    else if (!strcmp(words[0], "LIG"))
    {
        strcpy(ligfname, words[1]);
        // optsecho = "Ligand file is " + (std::string)ligfname;
        ligset = true;
        if (ligcmd) smset = smcmd;
        return 1;
    }
    else if (!strcmp(words[0], "SMILES"))
    {
        strcpy(smiles, words[1]);
        smset = true;
        return 1;
    }
    else if (!strcmp(words[0], "MCOORD"))
    {
        int j=0;
        // optsecho = "Metal coordination on residues ";
        MCoord mcr;
        i=1; if (!words[i]) throw 0xbad372;
        mcr.Z = Atom::Z_from_esym(words[1]);
        if (!mcr.Z) throw 0xbad372;

        i++; if (!words[i]) throw 0xbad372;
        mcr.charge = atoi(words[i]);

        i++; if (!words[i]) throw 0xbad372;
        for (; words[i]; i++)
        {
            if (words[i][0] == '-' && words[i][1] == '-') break;
            ResiduePlaceholder rp;
            rp.set(words[i]);
            mcr.coordres.push_back(rp);
        }

        mtlcoords.push_back(mcr);
        return i-1;
    }
    else if (!strcmp(words[0], "MOVIE"))
    {
        output_each_iter = true;
    }
    else if (!strcmp(words[0], "NODEPDB"))
    {
        activation_node = atoi(words[1]);
        strcpy(protafname, words[2]);
        optsecho = "Active PDB " + (std::string)protafname + " for node " + to_string(activation_node);
    }
    else if (!strcmp(words[0], "OUT"))
    {
        if (!words[1])
        {
            cout << "Missing output file name; check config file." << endl << flush;
            throw 0xbadf12e;
        }
        strcpy(outfname, words[1]);
        optsecho = "Output file is " + (std::string)outfname;
        return 1;
    }
    else if (!strcmp(words[0], "OUTPDB"))
    {
        outpdb_poses = atoi(words[1]);
        outpdb = words[2];
        return 2;
    }
    else if (!strcmp(words[0], "PATH"))
    {
        i=1;
        int nodeno = atoi(words[i]);

        if (!strcmp(words[i], "ABS")) i++;
        if (!strcmp(words[i], "REL")) i++;
        if (!strcmp(words[i], "RES")) i++;

        if (nodeno > 255)
        {
            cout << "Binding path is limited to 255 nodes." << endl << flush;
            throw 0xbad90de;
        }
        if (nodeno)
        {
            if ((nodeno) > pathnodes) pathnodes = nodeno;

            pathstrs.resize(nodeno+1);
            pathstrs[nodeno] = origbuff;
        }
        optsecho = "Path node set #" + to_string(nodeno);
        return i-1;
    }
    else if (!strcmp(words[0], "POSE"))
    {
        poses = atoi(words[1]);
        optsecho = "Number of poses: " + to_string(poses);
        return 1;
    }
    else if (!strcmp(words[0], "PREALIGN"))
    {
        use_prealign = true;
        prealign_residues = origbuff;
        return 1;
    }
    else if (!strcmp(words[0], "PROGRESS"))
    {
        progressbar = true;
        return 0;
    }
    else if (!strcmp(words[0], "CONGRESS"))
    {
        progressbar = false;
        return 0;
    }
    else if (!strcmp(words[0], "PROT"))
    {
        strcpy(protfname, words[1]);
        char* c = strchr(protfname, ':');
        if (c)
        {
            *c = 0;
            protstrand = *(++c);
        }
        else protstrand = 0;
        protset = true;
        // optsecho = "Protein file is " + (std::string)protfname;
        return 1;
    }
    else if (!strcmp(words[0], "REQSR"))
    {
        int reqrnode = atoi(words[1]);
        if (!strcmp(words[1], "all") || !strcmp(words[1], "ALL")) reqrnode = -1;
        for (i=2; words[i]; i++)
        {
            if (words[i][0] == '-' && words[i][1] == '-') break;
            ResiduePlaceholder rp;
            rp.set(words[i]);
            rp.node = reqrnode;
            required_contacts.push_back(rp);
        }
        return i-1;
    }
    else if (!strcmp(words[0], "RETRY"))
    {
        // triesleft = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "RLIM"))
    {
        _INTERA_R_CUTOFF = atof(words[1]);
        optsecho = "Interatomic radius limit: " + to_string(_INTERA_R_CUTOFF);
        return 1;
    }
    else if (!strcmp(words[0], "SEARCH"))
    {
        if (!words[1]) return 0;       // Stay with default.
        if (!strcmp(words[1], "BB"))
        {
            use_bestbind_algorithm = true;
            return 1;
        }
        else if (!strcmp(words[1], "TS"))
        {
            use_bestbind_algorithm = false;
            return 1;
        }
        else
        {
            cout << "Unknown search method " << words[1] << endl;
            throw 0xbad5eec;
        }
    }
    else if (!strcmp(words[0], "SIZE"))
    {
        size.x = atof(words[1]);
        if (words[2])
        {
            size.y = atof(words[2]);
            size.z = atof(words[3]);
        }
        else size.z = size.y = size.x;
        if (!size.x || !size.y || !size.z)
        {
            cout << "Pocket size cannot be zero in any dimension." << endl << flush;
            throw 0xbad512e;
        }
        optsecho = "Interatomic radius limit: " + to_string(size.x) + (std::string)"," + to_string(size.y) + (std::string)"," + to_string(size.z);
        return 3;
    }
    else if (!strcmp(words[0], "SOFT"))
    {
        dyn_strings.push_back(origbuff);
    }
    else if (!strcmp(words[0], "STATE"))
    {
        states.push_back(origbuff);
        optsecho = "Added state " + (std::string)origbuff;
        return 0;
    }
    else if (!strcmp(words[0], "TEMPLATE"))
    {
        tplset = (strcmp(words[1], "off") != 0) && (strcmp(words[1], "OFF") != 0);
        if (tplset)
        {
            strcpy(tplfname, words[1]);
            if (words[2])
            {
                tprfset = true;
                strcpy(tplrfnam, words[2]);
            }
        }
        return 1;
    }
    else if (!strcmp(words[0], "TRIP"))
    {
        optsecho = "Added trip clashables ";
        for (i = 1; words[i]; i++)
        {
            if (words[i][0] == '-' && words[i][1] == '-') break;
            tripswitch_clashables.push_back(atoi(words[i]));
            optsecho += (std::string)words[i] + (std::string)" ";
        }
        return i-1;
    }
    else if (!strcmp(words[0], "WET"))
    {
        wet_environment = true;
    }

    return 0;
}

void read_config_file(FILE* pf)
{
    char buffer[65536];
    int i;

    while (!feof(pf))
    {
        char* wgas = fgets(buffer, 1015, pf);
        origbuff = buffer;
        if (buffer[0] >= ' ' && buffer[0] != '#')
        {
            char** words = chop_spaced_words(buffer);
            if (!words) continue;

            interpret_config_line(words);

            delete words;
        }
        buffer[0] = 0;
    }
}

void prepare_initb()
{
    int i, j;

    if (differential_dock)
    {
        initial_binding = new float[seql+4];
        initial_vdWrepl = new float[seql+4];

        for (i=0; i<seql+4; i++) initial_binding[i] = initial_vdWrepl[i] = 0;

        std::vector<AminoAcid*> preres = protein->get_residues_near(pocketcen, pre_ligand_multimol_radius);
        int qpr = preres.size();
        AminoAcid* preaa[seql+4];
        MovabilityType aamov[seql+4];

        for (i=0; i<seql+4; i++) preaa[i] = nullptr;
        for (i=0; i<qpr; i++)
        {
            preaa[i] = preres[i];
            Atom* CA = preaa[i]->get_atom("CA");
            if (!CA)
            {
                cout << "Residue " << preaa[i]->get_residue_no() << " is missing its CA." << endl << flush;
                throw 0xbad12e5;
            }
            float r = CA->get_location().get_3d_distance(pocketcen);

            aamov[i] = preaa[i]->movability;

            if (r > pre_ligand_flex_radius) preaa[i]->movability = MOV_NONE;
        }

        Molecule* prem[seql+maxh2o+4];
        for (i=0; i<qpr; i++) prem[i] = reinterpret_cast<Molecule*>(preaa[i]);
        if (waters)
        {
            for (i=0; i<maxh2o; i++) prem[i+qpr] = waters[i];
        }
        int qpm = qpr + maxh2o;
        prem[qpm] = nullptr;

        #if !flexion_selection
        bool preconform;

        #if preconform_protein
        preconform = flex;
        #else
        preconform = flex && (waters != nullptr) && differential_dock;
        #endif

        if (preconform && pre_ligand_iteration_ratio)
        {
            Molecule** delete_me;
            Molecule::conform_molecules(
                prem /*reinterpret_cast<Molecule**>(preaa)*/,
                delete_me = protein->all_residues_as_molecules(),
                iters*pre_ligand_iteration_ratio
            );
            delete delete_me;
        }
        #endif

        for (i=0; i<_INTER_TYPES_LIMIT; i++) total_binding_by_type[i] = 0;

        for (i=0; i<qpr; i++)
        {
            int resno = preaa[i]->get_residue_no();
            #if _DBG_TOOLARGE_DIFFNUMS
            std::string ibdbg = to_string(resno) + (std::string)" ibdbg:\n";
            #endif

            for (j=0; j<qpm; j++)
            {
                if (j == i) continue;
                float f = reinterpret_cast<Molecule*>(preaa[i])->get_intermol_binding(reinterpret_cast<Molecule*>(prem[j]), j==0);

                #if _DBG_TOOLARGE_DIFFNUMS
                if (f) ibdbg += to_string(preaa[j]->get_residue_no()) + (std::string)" " + to_string(f) + (std::string)"\n";
                #endif

                initial_binding[resno] += f;
                initial_vdWrepl[resno] += preaa[i]->get_vdW_repulsion(prem[j]);
            }

            #if _DBG_TOOLARGE_DIFFNUMS
            if (fabs(initial_binding[resno]) >= 200) cout << ibdbg << endl;
            #endif
        }

        for (i=0; i<_INTER_TYPES_LIMIT; i++) init_total_binding_by_type[i] = total_binding_by_type[i];

        for (i=0; i<qpr; i++) preaa[i]->movability = aamov[i];
    }
}

void prepare_acv_bond_rots()
{
    int i;
    if (active_bond_rots.size())
    {
        for (i=0; i<active_bond_rots.size(); i++)
        {
            active_bond_rots[i].atom = protein->get_atom(active_bond_rots[i].resno, active_bond_rots[i].aname.c_str());
            active_bond_rots[i].btom = protein->get_atom(active_bond_rots[i].resno, active_bond_rots[i].bname.c_str());

            if (!active_bond_rots[i].atom) cout << "WARNING: " << active_bond_rots[i].resno << ":" << active_bond_rots[i].aname
                << " not found in protein!" << endl;
            if (!active_bond_rots[i].btom) cout << "WARNING: " << active_bond_rots[i].resno << ":" << active_bond_rots[i].bname
                << " not found in protein!" << endl;

            if (active_bond_rots[i].atom && active_bond_rots[i].btom)
            {
                active_bond_rots[i].bond = active_bond_rots[i].atom->get_bond_between(active_bond_rots[i].btom);
                #if _debug_active_bond_rot
                active_bond_rots[i].bond->echo_on_rotate = true;
                #endif
            }
            else active_bond_rots[i].bond = nullptr;
        }
    }
}

void do_tumble_spheres(Point l_pocket_cen)
{
    int i, j, l, n;
    float lig_min_int_clsh = ligand->get_internal_clashes();

    // Begin tumble sphere behavior.
    std::vector<AminoAcid*> tsphres = protein->get_residues_near(l_pocket_cen, size.magnitude()+6);
    int tsphsz = tsphres.size();
    float outer_sphere[tsphsz+4], inner_sphere[tsphsz+4];

    for (i=0; i<tsphsz+4; i++) outer_sphere[i] = inner_sphere[i] = 0;

    pocketsize = protein->estimate_pocket_size(tsphres);
    ligbbox = ligand->get_bounding_box();

    for (i=0; !ligbbox.fits_inside(pocketsize) && i<100; i++)
    {
        ligand->crumple(fiftyseventh*30);
    }

    for (i=0; i<tsphsz; i++)
    {
        #if use_exclusions
        if (exclusion.size()
                &&
                std::find(exclusion.begin(), exclusion.end(), tsphres[i]->get_residue_no())!=exclusion.end()
        )
        {
            tsphres.erase(tsphres.begin()+i);
            tsphsz--;
            continue;
        }
        #endif

        // TODO: Algorithmically determine more accurate values based on interaction type, etc.
        outer_sphere[i] = tsphres[i]->get_reach() + 2.5;
        inner_sphere[i] = tsphres[i]->get_reach() / 3 + 1;
    }

    const SCoord xaxis = Point(1,0,0), yaxis = Point(0,1,0), zaxis = Point(0,0,1);
    float loneliness=0, blone=0, xrad, yrad, zrad, lrad, step, bestxr, bestyr, bestzr, score, worth, weight, bestscore;
    const int ac = ligand->get_atom_count();
    Pose besp(ligand);
    #if _DBG_TUMBLE_SPHERES
    std::string tsdbg = "", tsdbgb = "";
    #endif

    if (ligbbox.x > ligbbox.y && pocketsize.x < pocketsize.y) ligand->rotate(zaxis, square);
    if (ligbbox.x > ligbbox.z && pocketsize.x < pocketsize.z) ligand->rotate(yaxis, square);
    if (ligbbox.y > ligbbox.x && pocketsize.y < pocketsize.x) ligand->rotate(zaxis, square);
    if (ligbbox.y > ligbbox.z && pocketsize.y < pocketsize.z) ligand->rotate(xaxis, square);
    if (ligbbox.z > ligbbox.x && pocketsize.z < pocketsize.x) ligand->rotate(yaxis, square);
    if (ligbbox.z > ligbbox.y && pocketsize.z < pocketsize.y) ligand->rotate(xaxis, square);

    step = fiftyseventh*30;
    bestscore = -Avogadro;
    float lonely_step = 1.0 / loneliest.get_3d_distance(l_pocket_cen);
    #if _DBG_LONELINESS
    cout << "Loneliest point " << loneliest << " is " << loneliest.get_3d_distance(l_pocket_cen) << "A from pocketcen " << l_pocket_cen << "." << endl;
    cout << "Pocket size is " << pocketsize << " vs. ligand bounding box " << ligbbox << endl;
    #endif
    if (isnan(lonely_step) || lonely_step < 0.1) lonely_step = 0.1;

    #if pocketcen_is_loneliest
    if (1)
    {
        ligand->recenter(l_pocket_cen);
    #else
    for (loneliness=0; loneliness <= 1; loneliness += lonely_step)
    {
        float centeredness = 1.0 - loneliness;
        Point tmpcen(loneliest.x * loneliness + l_pocket_cen.x * centeredness,
                     loneliest.y * loneliness + l_pocket_cen.y * centeredness,
                     loneliest.z * loneliness + l_pocket_cen.z * centeredness
                    );
        ligand->recenter(tmpcen);
    #endif

        #if _DBG_LONELINESS && !pocketcen_is_loneliest
        cout << "Ligand is " << loneliness << " lonely centered at " << tmpcen << "." << endl;
        #endif

        for (xrad=0; xrad <= M_PI*2; xrad += step)
        {
            for (yrad=0; yrad <= M_PI*2; yrad += step)
            {
                for (zrad=0; zrad <= M_PI*2; zrad += step)
                {
                    ligbbox = ligand->get_bounding_box();

                    if (ligbbox.x > ligbbox.y && pocketsize.x < pocketsize.y) continue;
                    if (ligbbox.x > ligbbox.z && pocketsize.x < pocketsize.z) continue;
                    if (ligbbox.y > ligbbox.x && pocketsize.y < pocketsize.x) continue;
                    if (ligbbox.y > ligbbox.z && pocketsize.y < pocketsize.z) continue;
                    if (ligbbox.z > ligbbox.x && pocketsize.z < pocketsize.x) continue;
                    if (ligbbox.z > ligbbox.y && pocketsize.z < pocketsize.y) continue;

                    Bond** rb = ligand->get_rotatable_bonds();

                    if (!rb) n = 0;
                    else for (n=0; rb[n]; n++);		// Get count.

                    l = 0;
                    lrad = 0;
                _xyzl_loop:
                    if (ligand->get_internal_clashes() >= lig_min_int_clsh*5+5) goto _xyzl_skip_loop;

                    score = 0;
                    #if _DBG_TUMBLE_SPHERES
                    tsdbg = "";
                    // cout << ligand->get_internal_clashes() << " vs. " << lig_min_int_clsh << endl;
                    #endif
                    for (i=0; i<ac; i++)
                    {
                        Atom* a = ligand->get_atom(i);
                        intera_type it = vdW;

                        for (j=0; j<tsphsz; j++)
                        {
                            worth = 0.4;
                            if (a->get_charge() && tsphres[j]->get_charge()
                                    &&
                                    sgn(a->get_charge()) == -sgn(tsphres[j]->get_charge())
                            )
                            {
                                it = ionic;
                                worth = 100;
                            }
                            else if (a->get_charge() || a->is_polar())
                            {
                                it = hbond;
                                worth = 40;
                            }
                            else if (a->is_pi())
                            {
                                it = pi;
                                worth = 7;
                            }

                            #if active_persistence
                            worth *= residue_binding_multiplier(tsphres[j]->get_residue_no());
                            #endif

                            if (tsphres[j]->capable_of_inter(it))
                            {
                                float r = a->get_location().get_3d_distance(tsphres[j]->get_atom_location("CA"));
                                if (r <= outer_sphere[j])
                                {
                                    if (r > inner_sphere[j])
                                    {
                                        weight = 1;

                                        if (tsphres[j]->priority)
                                        {
                                            weight = 1.25;		// Extra weight for residues mentioned in a CEN RES or PATH RES parameter.
                                        }

                                        #if !tumble_spheres_include_vdW
                                        if ((worth*weight) < 1) continue;
                                        #endif

                                        score += worth * weight;
                                        #if _DBG_TUMBLE_SPHERES
                                        tsdbg += std::string("+ ")
                                                +  std::string(a->name) + std::string(" reaches ") + std::string(tsphres[j]->get_3letter())
                                                +  std::to_string(tsphres[j]->get_residue_no()) + std::string(".  ");
                                        #endif
                                    }
                                    else
                                    {
                                        score -= 200;
                                        #if _DBG_TUMBLE_SPHERES
                                        tsdbg += std::string("- ")
                                                +  std::string(a->name) + std::string(" clashes ") + std::string(tsphres[j]->get_3letter())
                                                +  std::to_string(tsphres[j]->get_residue_no()) + std::string(".  ");
                                        #endif
                                    }
                                }
                            }
                        }
                    }

                    #if !pocketcen_is_loneliest
                    if (score > 0) score *= 1.0 + 0.1 * centeredness;
                    #endif

                    if (score > bestscore)
                    {
                        besp.copy_state(ligand);
                        blone = loneliness;
                        bestxr = xrad;
                        bestyr = yrad;
                        bestzr = zrad;
                        bestscore = score;

                        #if _DBG_TUMBLE_SPHERES
                        tsdbgb = tsdbg;

                        cout << "Tumble score " << score << " for ligand box " << ligand->get_bounding_box() << endl;


                        #if output_tumble_debug_docs
                        int u, v, w;
                        char protfttl[1000];
                        strcpy(protfttl, protfname);

                        char** lwords = chop_spaced_words(protfttl, '/');

                        for (u=0; lwords[u]; u++);
                        u--;

                        char fname[1000];
                        sprintf(fname, "output/tumble_%s_%d_%d_%d_%f.dock",
                                lwords[u],
                                (int)(xrad*fiftyseven),
                                (int)(yrad*fiftyseven),
                                (int)(zrad*fiftyseven),
                                score);
                        cout << fname << endl;
                        std::ofstream tspdbdat(fname, std::ofstream::out);

                        tspdbdat << "PDB file: " << protfname << endl;
                        tspdbdat << "Pose: 1\nNode: 0\nPDBDAT:\n";

                        int lac = ligand->get_atom_count();
                        for (u=0; u<lac; u++) ligand->get_atom(u)->stream_pdb_line(tspdbdat, 9000+u);

                        int pseql = protein->get_seq_length();
                        v = 1;
                        for (u = 1; u < pseql; u++)
                        {
                            AminoAcid* dbgaa = protein->get_residue(u);
                            if (dbgaa)
                            {
                                int aaac = dbgaa->get_atom_count();
                                for (w=0; w<aaac; w++)
                                {
                                    Atom* dbga = dbgaa->get_atom(w);
                                    if (!strcmp(dbga->name, "CA") || !strcmp(dbga->name, "CB")) dbga->stream_pdb_line(tspdbdat, v++);
                                }
                            }
                        }
                        tspdbdat << "END" << endl;
                        tspdbdat.close();
                        #endif
                        #endif
                    }

                _xyzl_skip_loop:

                    if (rb && rb[l])
                    {
                        rb[l]->rotate(step);

                        lrad += step;
                        if (lrad >= M_PI*2)
                        {
                            l++;
                            if (l < n) goto _xyzl_loop;
                        }
                        else goto _xyzl_loop;
                    }

                    ligand->rotate(zaxis, step);
                }		// zrad.
                ligand->rotate(yaxis, step);
            }			// yrad.
            ligand->rotate(xaxis, step);
        }				// xrad.

        #if !pocketcen_is_loneliest
        if (bestscore >= (ligand->get_atom_count()*13)) break;
        #endif

        #if _DBG_LONELINESS
        cout << "Best score: " << bestscore << endl;
        #endif
    }					// loneliness.

    #if _DBG_TUMBLE_SPHERES
    cout << "Tumble sphere best score " << bestscore << " for "
        << "x" << bestxr*fiftyseven << "deg, "
        << "y" << bestyr*fiftyseven << "deg, "
        << "z" << bestzr*fiftyseven << "deg."
        << " (" << blone << " lonely)."
        << endl;
    cout << tsdbgb << endl;
    #endif

    // Load the best ligand conformer.
    besp.restore_state(ligand);

    // Minimize ligand clashes.
    #if prerot_sidechains_from_ligand
    for (i=0; i<tsphsz; i++)
    {
        Bond** tsphb = tsphres[i]->get_rotatable_bonds();
        if (tsphb)
        {
            for (j=0; tsphb[j]; j++)
            {
                float rad=0, bestrad=0, clash, bestclash=6.25e24;
                for (rad=0; rad < M_PI*2; rad += step)
                {
                    clash = tsphres[i]->get_intermol_clashes(ligand);

                    if (clash < bestclash)
                    {
                        bestrad = rad;
                        bestclash = clash;
                    }

                    tsphb[j]->rotate(step);
                }

                tsphb[j]->rotate(bestrad);
            }
            // delete[] tsphb;
        }
    }
    #endif
    // End tumble sphere behavior.
}

void apply_protein_specific_settings(Protein* p)
{
    int i, n;

    n = atomto.size();
    for (i=0; i<n; i++)
    {
        char buffer[1024];
        strcpy(buffer, atomto[i].c_str());
        char** words = chop_spaced_words(buffer);

        if (!words[1]) throw -1;
        AminoAcid* aa = protein->get_residue_bw(words[1]);
        if (!words[2]) throw -1;
        char* aname = words[2];
        if (!words[3]) throw -1;
        AminoAcid* target = protein->get_residue_bw(words[3]);
        if (words[4]) throw -1;

        if (!aa)
        {
            cout << "Warning: residue " << words[1] << " not found." << endl;
            continue;
        }

        if (!target)
        {
            cout << "Warning: residue " << words[3] << " not found." << endl;
            continue;
        }

        Atom* a = aa->get_atom(aname);
        if (!strcmp("EXTENT", aname)) a = aa->get_reach_atom();
        if (!a)
        {
            cout << "Warning: atom not found " << *aa << ":" << aname << endl;
            continue;
        }

        MovabilityType aamov = aa->movability;
        aa->movability = MOV_FLEXONLY;
        aa->conform_atom_to_location(a->name, target->get_CA_location());
        aa->movability = aamov;
    }

    dyn_motions.clear();
    n = dyn_strings.size();
    for (i=0; i<n; i++)
    {
        DynamicMotion dyn(protein);
        char buffer[1024];
        strcpy(buffer, dyn_strings[i].c_str());
        char** words = chop_spaced_words(buffer);
        AminoAcid* aa1 = protein->get_residue_bw(words[1]);
        AminoAcid* aa2 = protein->get_residue_bw(words[2]);
        int resno1 = aa1->get_residue_no();
        int resno2 = aa2->get_residue_no();
        dyn.name = (std::string)words[1] + (std::string)"-" + (std::string)words[2];
        dyn.type = dyn_bend;
        dyn.start_resno.from_string(resno1<resno2 ? words[1] : words[2]);
        dyn.end_resno.from_string(resno1>resno2 ? words[1] : words[2]);
        dyn.fulcrum_resno.from_string(words[1]);
        dyn.axis = compute_normal(aa1->get_CA_location(), pocketcen, aa2->get_CA_location());
        dyn.bias = 30*fiftyseventh;
        dyn_motions.push_back(dyn);
    }
}

int main(int argc, char** argv)
{
    char buffer[65536];
    int i, j;

    _momentum_rad_ceiling = fiftyseventh * 5;

    for (i=0; i<65536; i++) buffer[i] = 0;
    #if active_persistence
    for (i=0; i<active_persistence_limit; i++) active_persistence_resno[i] = 0;

    #if _DBG_RESBMULT
    cout << "Cleared active persistence resnos." << endl;
    #endif
    #endif

    for (i=0; i<256; i++)
        configfname[i] = protfname[i] = protafname[i] = ligfname[i] = 0;

    time_t began = time(NULL);

    strcpy(configfname, "primarydock.config");

    smcmd = false;
    for (i=1; i<argc; i++)
    {
        if (argv[i][0] == '-' && argv[i][1] == '-')
        {
            argv[i] += 2;
            for (j=0; argv[i][j]; j++) if (argv[i][j] >= 'a' && argv[i][j] <= 'z') argv[i][j] &= 0x5f;
            j = interpret_config_line(&argv[i]);
            if (ligset) ligcmd = true;
            if (smset) smcmd = true;
            // if (optsecho.size()) cout << optsecho << endl;
            argv[i] -= 2;
            i += j;
        }
        else
        {
            strcpy(configfname, argv[i]);
            configset = true;
        }
    }

    FILE* pf = fopen(configfname, "r");
    if (!pf)
    {
        cout << "Config file not found: " << configfname << ", exiting." << endl;
        return 0xbadf12e;
    }

    read_config_file(pf);
    fclose(pf);

    for (i=1; i<argc; i++)
    {
        if (argv[i][0] == '-' && argv[i][1] == '-')
        {
            argv[i] += 2;
            for (j=0; argv[i][j]; j++) if (argv[i][j] >= 'a' && argv[i][j] <= 'z') argv[i][j] &= 0x5f;
            j = interpret_config_line(&argv[i]);
            if (optsecho.size()) cout << optsecho << endl;
            argv[i] -= 2;
            i += j;
        }
    }

    char* pcntp = strstr(outfname, "%p");
    if (pcntp)
    {
        char tmp[512], protn[64];
        *(pcntp++) = 0;
        *(pcntp++) = 0;
        strcpy(protn, strrchr(protfname, '/')+1);
        char* dot = strchr(protn, '.');
        if (dot) *dot = 0;
        sprintf(tmp, "%s%s%s", outfname, protn, pcntp);
        strcpy(outfname, tmp);
    }

    char* pcntl = strstr(outfname, "%l");
    if (pcntl)
    {
        char tmp[512], lign[64];
        *(pcntl++) = 0;
        *(pcntl++) = 0;
        strcpy(lign, strrchr(ligfname, '/')+1);
        char* dot = strchr(lign, '.');
        if (dot) *dot = 0;
        sprintf(tmp, "%s%s%s", outfname, lign, pcntl);
        strcpy(outfname, tmp);
    }

    #if _DBG_SPACEDOUT
    cout << "Starting a file outstream: " << outfname << endl;
    #endif
    output = new std::ofstream(outfname, std::ofstream::out);
    if (!output) return -1;

    pre_ligand_flex_radius = size.magnitude();
    pre_ligand_multimol_radius = pre_ligand_flex_radius + (default_pre_ligand_multimol_radius - default_pre_ligand_flex_radius);

    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Loaded config file." << endl;
    #endif

    if (kcal) kJmol_cutoff /= _kcal_per_kJ;
    drift = 1.0 / (iters/25+1);

    char protid[255];
    char* slash = strrchr(protfname, '/');
    if (!slash) slash = strrchr(protfname, '\\');
    strcpy(protid, slash ? slash+1 : protfname );
    char* dot = strchr(protid, '.');
    if (dot) *dot = 0;

    Protein pose_proteins[poses];
    Molecule pose_ligands[poses+1];
    protein = &pose_proteins[0]; // new Protein(protid);
    pf = fopen(protfname, "r");
    if (!pf)
    {
        cout << "Error trying to read " << protfname << endl;
        return 0xbadf12e;
    }
    protein->load_pdb(pf, 0, protstrand ?: 'A');
    apply_protein_specific_settings(protein);
    fclose(pf);
    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Loaded protein." << endl;
    #endif

    Pose tmp_pdb_residue[poses+1][protein->get_end_resno()+1];
    Pose tmp_pdb_waters[poses+1][omaxh2o+1];
    Pose tmp_pdb_ligand[poses+1];
    Point tmp_pdb_metal_locs[poses+1][mtlcoords.size()+1];

    if (tplset)
    {
        ptemplt = new Protein("template");
        pf = fopen(tplfname, "r");
        if (!pf)
        {
            cout << "Error trying to read " << tplfname << endl;
            return 0xbadf12e;
        }
        ptemplt->load_pdb(pf);
        fclose(pf);

        #if _dbg_homology
        cout << "Homology template is " << tplfname << endl;
        #endif

        if (tprfset)
        {
            ptplref = new Protein("reference");
            pf = fopen(tplrfnam, "r");
            if (!pf)
            {
                cout << "Error trying to read " << tplrfnam << endl;
                return 0xbadf12e;
            }
            ptplref->load_pdb(pf);
            fclose(pf);
            protein->homology_conform(ptemplt, ptplref);
        }
        else protein->homology_conform(ptemplt, protein);

        temp_pdb_file = (std::string)"tmp/" + std::to_string(pid) + (std::string)"_homolog.pdb";

        pf = fopen(temp_pdb_file.c_str(), "wb");
        protein->save_pdb(pf);
        fclose(pf);
    }

    if (hydrogenate_pdb)
    {
        int resno, endres = protein->get_end_resno();
        for (resno=1; resno<=endres; resno++)
        {
            AminoAcid* res = protein->get_residue(resno);
            if (res) res->hydrogenate();
        }

        temp_pdb_file = (std::string)"tmp/" + std::to_string(pid) + (std::string)"_hydro.pdb";

        pf = fopen(temp_pdb_file.c_str(), "wb");
        protein->save_pdb(pf);
        fclose(pf);
    }

    int l, j1, i2, miter;

    if (mtlcoords.size())
    {
        for (i=0; i<mtlcoords.size(); i++)
        {
            int charge_left = mtlcoords[i].charge;
            Point lpt;
            Molecule** lmc = new Molecule*[mtlcoords[i].coordres.size()+4];
            lmc[0] = new Molecule("lcm");
            Atom* lmtl = lmc[0]->add_atom(Atom::esym_from_Z(mtlcoords[i].Z), Atom::esym_from_Z(mtlcoords[i].Z), nullptr, 0);
            lmtl->increment_charge(mtlcoords[i].charge);
            mtlcoords[i].mtl = lmtl;
            lmc[0]->movability = MOV_ALL;

            l = 1;
            for (j=0; j<mtlcoords[i].coordres.size(); j++)
            {
                mtlcoords[i].coordres[j].resolve_resno(protein);
                AminoAcid* aa = protein->get_residue(mtlcoords[i].coordres[j].resno);
                if (aa)
                {
                    aa->movability = MOV_FLEXONLY;
                    lmc[l++] = (Molecule*)aa;
                    Atom** Ss = aa->get_most_bindable(1, lmtl);
                    lpt = lpt.add((Ss && Ss[0]) ? Ss[0]->get_location() : aa->get_barycenter());

                    // If cysteine, make thiolate form.
                    if (aa->is_thiol() && charge_left)
                    {
                        if (Ss)
                        {
                            Atom* S = Ss[0];
                            Atom* H = S->is_bonded_to("H");
                            if (H)
                            {
                                aa->delete_atom(H);
                                S->increment_charge(-1);
                                charge_left--;
                            }
                        }
                    }

                    if (l <= 2) aa->add_existing_atom(lmtl);
                    aa->coordmtl = lmtl;
                }
            }
            lmc[l] = nullptr;
            if (l > 1)
            {
                l--;
                lpt.x /= l; lpt.y /= l; lpt.z /= l;
                l++;

                lmtl->move(lpt);
            }

            lmtl->aaletter = '\0';
            strcpy(lmtl->aa3let, "MTL");
            lmtl->residue = 0;

            Molecule::conform_molecules(lmc, 50);

            AminoAcid* can_reach_metal[256];
            int num_can_reach = protein->get_residues_can_clash_ligand(can_reach_metal, lmc[0], lmtl->get_location(), Point(2.5,2.5,2.5), nullptr);
            bool cr_eq_mc[num_can_reach];

            for (miter=0; miter<20; miter++)
            {
                for (j1=0; j1<num_can_reach; j1++)
                {
                    bool found = false;
                    if (!miter)
                    {
                        for (i2=0; i2<mtlcoords[i].coordres.size(); i2++)
                        {
                            if (mtlcoords[i].coordres[i2].resno == can_reach_metal[j1]->get_residue_no()) found = true;
                        }
                        cr_eq_mc[j1] = found;
                    }
                    else found = cr_eq_mc[j1];

                    if (found) continue;

                    Atom* a = can_reach_metal[j1]->get_nearest_atom(lmtl->get_location());
                    float r = a->distance_to(lmtl);
                    float vdW = lmtl->get_vdW_radius() + a->get_vdW_radius();
                    if (r < vdW)
                    {
                        SCoord to_move = lmtl->get_location().subtract(a->get_location());
                        to_move.r = vdW - r;
                        lmtl->move_rel(&to_move);
                    }
                }
            }
            Molecule::conform_molecules(lmc, 50);

            for (j=0; j<mtlcoords[i].coordres.size(); j++)
            {
                AminoAcid* aa = protein->get_residue(mtlcoords[i].coordres[j].resno);
                if (aa) aa->movability = MOV_PINNED;
            }
        }

        temp_pdb_file = (std::string)"tmp/" + std::to_string(pid) + (std::string)"_metal.pdb";

        pf = fopen(temp_pdb_file.c_str(), "wb");
        protein->save_pdb(pf);
        fclose(pf);
    }

    if (bridges.size())
    {
        reconnect_bridges();

        temp_pdb_file = (std::string)"tmp/" + std::to_string(pid) + (std::string)"_bridged.pdb";

        pf = fopen(temp_pdb_file.c_str(), "wb");
        protein->save_pdb(pf);
        fclose(pf);
    }

    prepare_acv_bond_rots();

    for (i=0; i<required_contacts.size(); i++)
    {
        required_contacts[i].resolve_resno(protein);
    }

    if (!CEN_buf.length())
    {
        cout << "Error: no binding pocket center defined." << endl;
        return 0xbadb19d;
    }

    strcpy(buffer, CEN_buf.c_str());
    char** words = chop_spaced_words(buffer);
    pocketcen = pocketcen_from_config_words(words, nullptr);
    loneliest = protein->find_loneliest_point(pocketcen, size);

    #if pocketcen_is_loneliest
    pocketcen = loneliest;
    #endif


    #if !_use_groups
    for (i=0; i<3; i++)
    {
        alignment_aa[i] = nullptr;
    }
    #endif

    if (waters)
    {
        float szscale = 0.5;
        for (i=0; i<maxh2o; i++)
        {
            waters[i]->recenter(Point(
                pocketcen.x + frand(-size.x * szscale, size.x * szscale),
                pocketcen.y + frand(-size.y * szscale, size.y * szscale),
                pocketcen.z + frand(-size.z * szscale, size.z * szscale)
            ));
        }
    }

    pktset = true;

    protein->mcoord_resnos = mcoord_resno;

    for (i=0; mcoord_resno[i]; i++) addl_resno[i] = mcoord_resno[i];
    for (l=0; l < tripswitch_clashables.size(); l++) addl_resno[i+l] = tripswitch_clashables[l];
    addl_resno[i+l] = 0;

    // Load the ligand or return an error.
    // Molecule m(ligfname);
    // ligand = &m;
    for (l=0; l<=poses; l++)
    {
        ligand = &pose_ligands[l];
        ligand->set_name(ligfname);
        char* ext = get_file_ext(ligfname);
        if (!ext)
        {
            cout << "Ligand file is missing its extension! " << ligfname << endl;
            return 0xbadf12e;
        }

        for (i=0; i<65536; i++) buffer[i] = 0;

        size_t wgaf;
        if (smset) ligand->from_smiles(smiles);
        else switch (ext[0])
        {
        case 's':
        case 'S':
            // SDF
            pf = fopen(ligfname, "r");
            if (!pf)
            {
                cout << "Error trying to read " << ligfname << endl;
                return 0xbadf12e;
            }
            wgaf = fread(buffer, 1, 65535, pf);
            fclose(pf);
            ligand->from_sdf(buffer);
            break;

        case 'p':
        case 'P':
            pf = fopen(ligfname, "r");
            if (!pf)
            {
                cout << "Error trying to read " << ligfname << endl;
                return 0xbadf12e;
            }
            ligand->from_pdb(pf);
            fclose(pf);
            break;

        default:
            cout << "Unrecognized ligand file extension: " << ext << endl;
            return 0xbadf12e;
        }

        ligand->minimize_internal_clashes();
    }
    ligand = &pose_ligands[0];


    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Loaded ligand." << endl;
    #endif

    Point box = ligand->get_bounding_box();

    if (debug) *debug << "Ligand bounding box corner (centered at zero): " << box.printable() << endl;
    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Ligand bounding box." << endl;
    #endif

    // Identify the ligand atom with the greatest potential binding.
    int k, n;

    #if !flexion_selection

    if (use_prealign) use_bestbind_algorithm = false;

    if (use_prealign)
    {        
        strcpy(buffer, prealign_residues.c_str());
        words = chop_spaced_words(buffer);

        for (n=0; words[n]; n++);      // Get length.

        Molecule** prealign_res = new Molecule*[n+4];
        Molecule** lig_grp = new Molecule*[n+4];

        for (i=0; i<n; i++)
        {
            int resno = atoi(words[i]);
            prealign_res[i] = protein->get_residue(resno);
        }

        ligand->recenter(pocketcen);

        // First line up ligand to stationary residues.
        lig_grp[0] = ligand;
        lig_grp[1] = nullptr;
        ligand->movability = MOV_NORECEN;
        Molecule::conform_molecules(lig_grp, prealign_res, prealign_iters, nullptr);

        // Then line up residues to ligand.
        if (flex) Molecule::conform_molecules(prealign_res, lig_grp, prealign_iters, nullptr);
        ligand->movability = MOV_ALL;

        delete words;
        delete prealign_res;
        delete lig_grp;
    }

    #endif

    // Best-Binding Code

    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Identified best binding ligand atoms." << endl;
    #endif

    Point nodecen = pocketcen;
    seql = protein->get_seq_length();
    int rstart = protein->get_start_resno();

    // Filter residues according to which ones are close enough to the spheroid to "reach" it.
    nodecen = pocketcen;

    // When docking with a metalloprotein, use this temporary Molecule for interactions the same as
    // we use AminoAcid objects, except don't attempt to flex the metals object.
    Molecule* met = protein->metals_as_molecule();
    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Created metals molecule." << endl;
    #endif

    float bclash = 0;

    for (l=0; l<seql; l++)
    {
        AminoAcid* laa = protein->get_residue(l+rstart);
        if (!laa) continue;
        for (n=l+1; n<seql; n++)
        {
            AminoAcid* naa = protein->get_residue(n+rstart);
            if (!naa) continue;
            bclash += laa->get_intermol_clashes(naa);
        }

        bclash += ligand->get_intermol_clashes(laa);

        if (met) bclash += laa->get_intermol_clashes(met);
    }
    if (met) bclash += ligand->get_intermol_clashes(met);
    if (debug) *debug << "Initial clashes: " << bclash << endl;

    // TODO: Output some basic stats: receptor, ligand, etc.
    cout << "PDB file: " << protfname << endl;
    if (output) *output << "PDB file: " << protfname << endl;
    cout << "Ligand: " << ligfname << endl;
    if (output) *output << "Ligand: " << ligfname << endl;
    if (differential_dock)
    {
        cout << "Differential dock." << endl;
        if (output) *output << "Differential dock." << endl;
    }
    cout << endl;
    if (output) *output << endl;

    i = poses*(triesleft+1)+8;
    j = pathnodes+2;
    DockResult dr[i][j];
    #if _dbg_find_blasted_segfault
    cout << "dr[" << i << "][" << j << "] allocated. " << dr << endl;
    #endif

    float rgnxform_r[i][pathnodes+2][PROT_MAX_RGN], rgnxform_theta[i][pathnodes+2][PROT_MAX_RGN], rgnxform_y[i][pathnodes+2][PROT_MAX_RGN];
    float rgnrot_alpha[i][pathnodes+2][PROT_MAX_RGN], rgnrot_w[i][pathnodes+2][PROT_MAX_RGN], rgnrot_u[i][pathnodes+2][PROT_MAX_RGN];

    int drcount = 0, qpr;

    if (pathnodes)
    {
        cout << pathnodes << " path node" << (pathnodes == 1 ? "" : "s") << "." << endl;
        if (output) *output << pathnodes << " path node" << (pathnodes == 1 ? "" : "s") << "." << endl;
    }
    else
    {
        cout << "Static dock - no path nodes." << endl;
        if (output) *output << "Static dock - no path nodes." << endl;
    }

    reaches_spheroid = new AminoAcid**[pathnodes+2];
    for (i=0; i<=pathnodes; i++) reaches_spheroid[i] = new AminoAcid*[SPHREACH_MAX];

    found_poses = 0;
    int wrote_acvmx = -1, wrote_acvmr = -1;
    float l_individual_clash_limit = individual_clash_limit - kJmol_cutoff;

_try_again:
    // srand(0xb00d1cca);
    srand(time(NULL));
    Point nodecens[pathnodes+1];
    for (pose = 1; pose <= poses; pose++)
    {
        ligand = &pose_ligands[pose];
        ligand->movability = MOV_ALL;

        last_ttl_bb_dist = 0;
        ligand->minimize_internal_clashes();
        float lig_min_int_clsh = ligand->get_internal_clashes();
        ligand->crumple(triangular);

        for (i=0; i<dyn_motions.size(); i++) dyn_motions[i].apply_absolute(0);

        int rcn = required_contacts.size();
        if (rcn)
        {
            ligand->delete_mandatory_connections();
            ligand->allocate_mandatory_connections(rcn);

            for (i=0; i<rcn; i++)
            {
                if (required_contacts[i].node >= 0 && required_contacts[i].node != nodeno) continue;
                Star s;
                s.paa = protein->get_residue(required_contacts[i].resno);
                if (s.n) ligand->add_mandatory_connection(s.pmol);
            }
        }


        // delete protein;
        // protein = new Protein(protfname);
        protein = &pose_proteins[pose-1];

        if (temp_pdb_file.length())
        {
            pf = fopen(temp_pdb_file.c_str(), "r");
            protein->load_pdb(pf);
            fclose(pf);
            apply_protein_specific_settings(protein);

            if (mtlcoords.size())
            {
                for (i=0; i<mtlcoords.size(); i++)
                {
                    for (j=0; j<mtlcoords[i].coordres.size(); j++)
                    {
                        AminoAcid* aa = protein->get_residue(mtlcoords[i].coordres[j].resno);
                        if (aa)
                        {
                            aa->coordmtl = mtlcoords[i].mtl;
                        }
                    }
                }
            }
        }
        else
        {
            pf = fopen(protfname, "r");
            protein->load_pdb(pf, 0, protstrand ?: 'A');
            fclose(pf);
            apply_protein_specific_settings(protein);
        }

        strcpy(buffer, CEN_buf.c_str());
        words = chop_spaced_words(buffer);
        pocketcen = pocketcen_from_config_words(words, nullptr);

        freeze_bridged_residues();
        prepare_initb();

        ligand->recenter(pocketcen);
        // cout << "Centered ligand at " << pocketcen << endl;
        std::vector<std::shared_ptr<AtomGroup>> agc = AtomGroup::get_potential_ligand_groups(ligand, mtlcoords.size() > 0);

        if (!use_bestbind_algorithm && !use_prealign)
        {
            do_tumble_spheres(pocketcen);

            #if debug_stop_after_tumble_sphere
            return 0;
            #endif
        }

        #if _DBG_STEPBYSTEP
        if (debug) *debug << "Pose " << pose << endl;
        #endif
        nodecen = pocketcen;
        nodecen.weight = 1;

        #if _dummy_atoms_for_debug
        dummies.clear();
        #endif

        for (nodeno=0; nodeno<=pathnodes; nodeno++)
        {
            movie_offset = iters * (nodeno /*+ (pose-1)*(pathnodes+1)*/);

            if (waters)
            {
                for (i = 0; i <= omaxh2o; i++)
                {
                    waters[i] = owaters[i];
                }
                maxh2o = omaxh2o;
            }

            // if (pathstrs.size() < nodeno) break;
            drift = initial_drift;

            if (echo_progress) cout << (time(NULL) - began) << " seconds: starting pose " << pose << " node " << nodeno << "..." << endl;

            #if internode_momentum_only_on_activation 
            conformer_momenta_multiplier = 1;
            #else
            conformer_momenta_multiplier = nodeno ? internode_momentum_mult : 1;
            #endif
            conformer_tumble_multiplier = 1;

            allow_ligand_360_tumble = (nodes_no_ligand_360_tumble ? (nodeno == 0) : true) && !use_prealign && !use_bestbind_algorithm;
            allow_ligand_360_flex   = (nodes_no_ligand_360_flex   ? (nodeno == 0) : true) /*&& !use_bestbind_algorithm*/;

            if (use_prealign) conformer_tumble_multiplier *= prealign_momenta_mult;
            if (use_bestbind_algorithm) conformer_tumble_multiplier *= prealign_momenta_mult;

            if (strlen(protafname) && nodeno == activation_node)
            {
                #if internode_momentum_only_on_activation 
                conformer_momenta_multiplier = nodeno ? internode_momentum_mult : 1;
                #endif

                #if prevent_ligand_360_on_activate
                allow_ligand_360_tumble = allow_ligand_360_flex = false;
                #endif

                // Persist the flexions of the side chains. 
                // TODO: Do not persist those residues whose positions are important to activation.
                float* sidechain_bondrots[seql+4];
                int sidechain_bondrotq[seql+4];
                for (i=0; i<seql+4; i++)
                {
                    sidechain_bondrots[i] = nullptr;
                    sidechain_bondrotq[i] = 0;
                }
                for (i=1; i<=seql; i++)
                {
                    Bond** b = protein->get_residue(i)->get_rotatable_bonds();
                    if (b)
                    {
                        int bq;
                        for (bq=0; b[bq]; bq++);                // Get count.
                        sidechain_bondrots[i] = new float[bq+2];
                        for (j=0; j<bq; j++)
                        {
                            sidechain_bondrots[i][j] = b[j]->total_rotations;
                        }
                        sidechain_bondrotq[i] = j;

                        // delete[] b;
                    }
                }

                // delete protein;
                // protein = new Protein(protafname);
                protein = &pose_proteins[pose-1];

                pf = fopen(protafname, "r");
                protein->load_pdb(pf);
                fclose(pf);
                apply_protein_specific_settings(protein);

                freeze_bridged_residues();
                prepare_initb();

                strcpy(buffer, CEN_buf.c_str());
                words = chop_spaced_words(buffer);
                pocketcen = pocketcen_from_config_words(words, nullptr);

                for (i=1; i<=seql; i++)
                {
                    Bond** b = protein->get_residue(i)->get_rotatable_bonds();
                    if (b)
                    {
                        int bq;

                        for (j=0; j<sidechain_bondrotq[i]; j++)
                        {
                            if (!b[j]) break;
                            b[j]->clear_moves_with_cache();
                            b[j]->rotate(sidechain_bondrots[i][j]);
                        }

                        // delete[] b;
                    }

                    delete sidechain_bondrots[i];
                }
            }

            #if active_persistence
            for (j=0; j<active_persistence_limit; j++) active_persistence_resno[j] = 0;

            #if _DBG_RESBMULT
            cout << "Cleared active persistence resnos." << endl;
            #endif

            #if active_persistence_noflex
            allow_ligand_flex = true;
            #endif
            #endif

            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Pose " << pose << endl << "Node " << nodeno << endl;
            #endif
            if (nodeno)
            {
                for (i=0; i<states.size(); i++)
                {
                    strcpy(buffer, states[i].c_str());
                    words = chop_spaced_words(buffer);
                    if (atoi(words[1]) == nodeno)
                    {
                        int sr = atoi(words[2]), er = atoi(words[3]);
                        float theta = atof(words[4]) * fiftyseventh;

                        Point sloc = protein->get_atom_location(sr, "CA"),
                              eloc = protein->get_atom_location(er, "CA");

                        LocatedVector lv = (SCoord)(sloc.subtract(eloc));
                        lv.origin = sloc;

                        int resno;
                        for (resno = sr; resno <= er; resno++)
                        {
                            AminoAcid* aa = protein->get_residue(resno);
                            if (aa)
                            {
                                MovabilityType mt = aa->movability;
                                aa->movability = MOV_ALL;

                                aa->rotate(lv, theta);

                                aa->movability = mt;
                            }
                        }
                    }
                }

                // nodecen = nodecen.add(&path[nodeno]);
                strcpy(buffer, pathstrs[nodeno].c_str());
                if (!strlen(buffer))
                {
                    cout << "Error in config file: path node " << nodeno << " is missing." << endl;
                    return 0xbadc09f;
                }
                words = chop_spaced_words(buffer);
                nodecen = pocketcen_from_config_words(&words[1], &nodecen);

                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Added whatever points together." << endl;
                #endif
            }

            loneliest = protein->find_loneliest_point(nodecen, size);
            // cout << "Loneliest is " << loneliest << endl;

            #if pocketcen_is_loneliest
            nodecen = loneliest;
            #endif

            Point lastnodecen = nodecen;
            ligcen_target = nodecen;
            nodecens[nodeno] = ligcen_target;

            #if redo_tumble_spheres_every_node
            if (!use_bestbind_algorithm && !use_prealign && (!prevent_ligand_360_on_activate)) do_tumble_spheres(ligcen_target);
            #endif

            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Saved last nodecen." << endl;
            #endif

            #if recenter_ligand_each_node
            // Move the ligand to the new node center.
            ligand->recenter(nodecen);
            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Molecule recenter (or not)." << endl;
            #endif
            ligand->reset_conformer_momenta();
            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Conformer momenta reset." << endl;
            #endif
            #endif

            sphres = protein->get_residues_can_clash_ligand(reaches_spheroid[nodeno], ligand, nodecen, size, addl_resno);
            for (i=sphres; i<SPHREACH_MAX; i++) reaches_spheroid[nodeno][i] = NULL;

            // Flexion Selection
            if (flex && !nodeno)
            {
                if (forced_static_resnos.size())
                {
                    for (i=0; i<forced_static_resnos.size(); i++)
                    {
                        forced_static_resnos[i].resolve_resno(protein);
                        AminoAcid* mvaa = protein->get_residue(forced_static_resnos[i].resno);
                        if (mvaa)
                        {
                            mvaa->movability = MOV_PINNED;
                            #if _dbg_flexion_selection
                            cout << mvaa->get_name() << " forced static." << endl;
                            #endif
                        }
                    }
                }
                #if flexion_selection

                flexible_resnos.clear();
                j = protein->get_end_resno();
                for (i=protein->get_start_resno(); i<=j; i++)
                {
                    AminoAcid* mvaa = protein->get_residue(i);
                    if (mvaa) mvaa->movability = min(MOV_FLXDESEL, mvaa->movability);
                }

                #if _dbg_null_flexions
                bool another_flex = false;
                #elif no_zero_flexions
                bool another_flex = true;
                #else
                bool another_flex = (frand(0,1) < 0.6);
                #endif

                while (another_flex)
                {
                    float bestwt = 0;
                    int besti = -1;
                    for (j=0; j<100; j++)
                        for (i=0; i<sphres; i++)
                        {
                            if (reaches_spheroid[nodeno][i]->movability != MOV_FLXDESEL) continue;
                            float weight = reaches_spheroid[nodeno][i]->get_aa_definition()->flexion_probability;
                            if (!weight) continue;

                            // Multiply weight by unrealized ligand binding potential.
                            float potential = reaches_spheroid[nodeno][i]->get_intermol_potential(ligand, true);
                            float adjusted_potential = fmin(1, potential / 1000);
                            Atom *nearest1, *nearest2;
                            nearest1 = reaches_spheroid[nodeno][i]->get_nearest_atom(ligand->get_barycenter());
                            if (!nearest1) throw 0xbadd157;
                            nearest2 = ligand->get_nearest_atom(nearest1->get_location());
                            if (!nearest2) throw 0xbadd157;
                            float nearr = fmax(1, nearest1->distance_to(nearest2) / 2);
                            adjusted_potential *= nearr;

                            // weight = (1.0 - ((1.0 - weight) / adjusted_potential)) / 2;
                            weight *= sqrt(adjusted_potential);

                            #if _dbg_flexion_selection
                            if (reaches_spheroid[nodeno][i]->get_residue_no() == 9262)
                                cout << reaches_spheroid[nodeno][i]->get_name() << " has weight " << weight << endl;
                            #endif

                            if ( /*weight >= bestwt &&*/ frand(0,100) < weight )
                            {
                                besti = i;
                                bestwt = weight;
                            }
                        }
                    if (besti >= 0)
                    {
                        if (!(reaches_spheroid[nodeno][besti]->movability & MOV_PINNED))
                            reaches_spheroid[nodeno][besti]->movability = MOV_FLEXONLY;
                        flexible_resnos.push_back(reaches_spheroid[nodeno][besti]->get_residue_no());
                        #if _dbg_flexion_selection
                        cout << "Selected " << reaches_spheroid[nodeno][besti]->get_name()
                                << " for flexion with a weight of " << bestwt << endl;
                        #endif
                    }
                    another_flex = (frand(0,1) < 0.6);
                }
                
                #if _dbg_flexion_selection
                cout << flexible_resnos.size() << " residues selected for flexion." << endl;
                #endif

                freeze_bridged_residues();

                if (forced_flexible_resnos.size())
                {
                    for (i=0; i<forced_flexible_resnos.size(); i++)
                    {
                        forced_flexible_resnos[i].resolve_resno(protein);
                        AminoAcid* mvaa = protein->get_residue(forced_flexible_resnos[i].resno);
                        if (mvaa)
                        {
                            mvaa->movability = MOV_FORCEFLEX;
                            flexible_resnos.push_back(mvaa->get_residue_no());
                            #if _dbg_flexion_selection
                            cout << mvaa->get_name() << " forced flexible." << endl;
                            #endif
                        }
                    }
                }
                #endif
            }

            #if _dbg_flexion_selection
            cout << endl;
            #endif

            for (i=0; i<_INTER_TYPES_LIMIT; i++) total_binding_by_type[i] = 0;
            for (i=0; i<sphres; i++)
            {
                #if use_exclusions
                if (exclusion.size()
                        &&
                        std::find(exclusion.begin(), exclusion.end(), reaches_spheroid[nodeno][i]->get_residue_no())!=exclusion.end()
                   )
                {
                    for (j=i; j<sphres; j++) reaches_spheroid[nodeno][j] = reaches_spheroid[nodeno][j+1];
                    sphres--;
                    reaches_spheroid[nodeno][sphres] = nullptr;
                }
                #endif
            }


            if (!nodeno)
            {
                float alignment_distance[5];
                for (l=0; l<3; l++)
                {
                    #if !_use_groups
                    alignment_aa[l]=0;
                    #endif
                    alignment_distance[l]=0;
                }

                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Initialize null AA pointer." << endl;
                #endif

                protein->set_conditional_basicities();
                if (use_bestbind_algorithm)
                {
                    std::vector<std::shared_ptr<ResidueGroup>> scg = ResidueGroup::get_potential_side_chain_groups(reaches_spheroid[nodeno], ligcen_target);
                    global_pairs = GroupPair::pair_groups(agc, scg, ligcen_target);

                    if (global_pairs.size() > 2)
                    {
                        // TODO: If the 2nd group is closer to the 1st group than the 3rd group is, swap the 2nd and 3rd groups.
                        Point grpcen1 = global_pairs[0]->ag->get_center(), grpcen2 = global_pairs[1]->ag->get_center(), grpcen3 = global_pairs[2]->ag->get_center();

                        if (grpcen1.get_3d_distance(grpcen2) > grpcen1.get_3d_distance(grpcen3))
                        {
                            std::shared_ptr<GroupPair> tmpg = global_pairs[2];
                            global_pairs[2] = global_pairs[1];
                            global_pairs[1] = tmpg;
                        }
                    }
                }
                ligand->recenter(ligcen_target);
                if (use_bestbind_algorithm)
                {
                    GroupPair::align_groups(ligand, global_pairs, false, 1);

                    #if _dbg_groupsel
                    cout << endl;
                    #endif

                    int gpn = global_pairs.size();
                    for (l=0; l<3 && l<gpn; l++)
                    {
                        ligand_groups[l] = *(global_pairs[l]->ag);
                        sc_groups[l] = *(global_pairs[l]->scg);
                    }
                }

                // Best-Binding Algorithm
                // Find a binding pocket feature with a strong potential binding to the ligand.
                std::string alignment_name = "";
                
                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Selected an alignment AA." << endl;
                #endif

                if (use_bestbind_algorithm && met)
                {
                    #if !_use_groups
                    alignment_aa[2] = alignment_aa[1];
                    alignment_aa[1] = alignment_aa[0];
                    alignment_aa[0] = met;
                    #endif
                    // alignment_name = "metal";
                }
                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Alignment AA." << endl;
                #endif

                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Aligned ligand to AA." << endl;
                cout << endl;
                #endif

                freeze_bridged_residues();
            }

            // float driftamt = 1.0 / (iters/25+1);
            // cout << pose << ":" << nodeno << " drift " << driftamt << endl;
            int iters_div = iters*0.259;

            Molecule* cfmols[SPHREACH_MAX+4];
            for (i=0; i<SPHREACH_MAX+4; i++) cfmols[i] = nullptr;
            gcfmols = cfmols;
            i=0;
            ligand->movability = MOV_ALL;
            cfmols[i++] = ligand;
            if (met)
            {
                met->movability = MOV_NONE;
                cfmols[i++] = met;
            }
            if (waters)
            {
                for (j=0; j<maxh2o; j++)
                {
                    waters[j]->movability = MOV_ALL;
                    waters[j]->reset_conformer_momenta();
                    cfmols[i++] = waters[j];
                }
            }


            for (j=0; j<sphres; j++)
            {
                #if ! flexion_selection
                if (reaches_spheroid[nodeno][j]->movability >= MOV_FLEXONLY) reaches_spheroid[nodeno][j]->movability = MOV_FLEXONLY;
                #endif
                if (!flex) reaches_spheroid[nodeno][j]->movability = MOV_FLXDESEL;
                cfmols[i++] = reaches_spheroid[nodeno][j];
                protein->get_residues_can_clash(reaches_spheroid[nodeno][j]->get_residue_no());
            }

            int cfmolqty = i;
            for (; i<SPHREACH_MAX; i++) cfmols[i] = NULL;

            ligand->reset_conformer_momenta();

            // Molecule** delete_me;
            int trsz = tripswitch_clashables.size();
            Molecule* trip[j = trsz+4];

            for (; j; j--) trip[j-1] = nullptr;

            for (j=0; j<trsz; j++)
                trip[j] = (Molecule*)protein->get_residue(tripswitch_clashables[j]);
            
            if (rcn)
            {
                for (i=0; i<rcn; i++)
                {
                    Star s;
                    s.paa = protein->get_residue(required_contacts[i].resno);
                    if (required_contacts[i].node >= 0 && required_contacts[i].node != nodeno) ligand->remove_mandatory_connection(s.pmol);
                    else ligand->add_mandatory_connection(s.pmol);
                }
            }

            int mcn;
            Molecule lm("MTL");
            if (mcn = mtlcoords.size())         // Assignment, not comparison.
            {
                for (i=0; i<mcn; i++)
                {
                    if (!mtlcoords[i].mtl) continue;                    
                    lm.add_existing_atom(mtlcoords[i].mtl);
                }

                lm.movability = MOV_NONE;
                cfmols[cfmolqty++] = &lm;
                cfmols[cfmolqty] = nullptr;
            }

            protein->find_residue_initial_bindings();
            freeze_bridged_residues();

            n = mtlcoords.size();
            if (n)
            {
                int j1;
                for (j=0; j<n; j++)
                {
                    int n1 = mtlcoords[j].coordres.size();
                    for (j1=0; j1<n1; j1++)
                    {
                        AminoAcid* aa = protein->get_residue(mtlcoords[j].coordres[j1].resno);
                        if (aa) aa->movability = MOV_PINNED;
                    }
                }
            }

            ligand->movability = (MovabilityType)(MOV_ALL - MOV_MC_AXIAL);
            if (!flex) for (j=0; j<sphres; j++)
            {
                reaches_spheroid[nodeno][j]->movability = MOV_FLXDESEL;
            }
            ligand->agroups = global_pairs;
            Molecule::conform_molecules(cfmols, iters, &iteration_callback, &GroupPair::align_groups_noconform);

            if (!nodeno) // && outpdb.length())
            {
                protein->get_internal_clashes(1, protein->get_end_resno(), true);

                n = protein->get_end_resno();
                for (j=1; j <= n; j++)
                {
                    AminoAcid* aa = protein->get_residue(j);
                    if (aa)
                    {
                        tmp_pdb_residue[pose][j].copy_state(aa);
                        #if _dbg_residue_poses
                        cout << "tmp_pdb_residue[" << pose << "][" << j << "].copy_state(" << aa->get_name() << ")" << endl;
                        #endif
                    }
                }
                if (waters)
                {
                    for (j=0; waters[j]; j++)
                    {
                        tmp_pdb_waters[pose][j].copy_state(waters[j]);
                    }
                }
                n = mtlcoords.size();
                for (j=0; j < n; j++)
                {
                    tmp_pdb_metal_locs[pose][j] = mtlcoords[j].mtl->get_location();
                }
                tmp_pdb_ligand[pose].copy_state(ligand);
            }

            #if active_persistence
            for (j=0; j<active_persistence_limit; j++) active_persistence_resno[j] = 0;

            #if _DBG_RESBMULT
            cout << "Cleared active persistence resnos." << endl;
            #endif
            #endif
            
            #if active_persistence_noflex
            allow_ligand_flex = true;
            #endif

            // Add the current pose/path sequentially to the dr[][] array.
            // If the path node # is zero:
            // If it is the first (zeroth) entry, set the pose number to 1.
            // Otherwise, go through all of the preceding entries and:
            // Any entry with a smaller kJ/mol, increment its pose# but remember the smallest pre-increment pose #
            // from the lot of them;
            // Claim that new smallest pose# (which might be 1) as your own.
            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Preparing output." << endl;
            #endif

            // Set the dock result properties and allocate the arrays.
            protein->set_conditional_basicities();
            #if _dbg_cond_basic
            AminoAcid* aadbg = protein->get_residue(155);
            cout << aadbg->get_name() << " charge = " << aadbg->get_charge() << endl;
            #endif
            dr[drcount][nodeno] = DockResult(protein, ligand, size, addl_resno, drcount, differential_dock);
            float btot = dr[drcount][nodeno].kJmol;
            float pstot = dr[drcount][nodeno].polsat;

            dr[drcount][nodeno].proximity = ligand->get_barycenter().get_3d_distance(nodecen);

            if (use_bestbind_algorithm)
            {
                dr[drcount][nodeno].miscdata += (std::string)"Best-Binding Pairs:\n";
                for (i=0; i<3 && i<global_pairs.size(); i++)
                {
                    n = global_pairs[i]->ag->atoms.size();
                    int j2;
                    for (j2=0; j2<n; j2++)
                        dr[drcount][nodeno].miscdata += (std::string)global_pairs[i]->ag->atoms[j2]->name + (std::string)" ";
                    
                    dr[drcount][nodeno].miscdata += (std::string)"- ";

                    n = global_pairs[i]->scg->aminos.size();
                    for (j2=0; j2<n; j2++)
                        dr[drcount][nodeno].miscdata += (std::string)global_pairs[i]->scg->aminos[j2]->get_name() + (std::string)" ";
                    
                    dr[drcount][nodeno].miscdata += (std::string)"\n";
                }
                dr[drcount][nodeno].miscdata += (std::string)"\n";
            }

            if (dyn_motions.size())
            {
                dr[drcount][nodeno].miscdata += (std::string)"Dynamic Motions:\n";
                for (i=0; i<dyn_motions.size(); i++)
                {
                    dr[drcount][nodeno].miscdata += (std::string)dyn_motions[i].name
                        + (std::string)" "
                        + std::to_string(dyn_motions[i].get_total_applied())
                        + (std::string)"\n";
                }
                dr[drcount][nodeno].miscdata += (std::string)"\n";
            }

            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Allocated memory." << endl;
            #endif
            
            std::ostringstream pdbdat;

            // Prepare a partial PDB of the ligand atoms and all involved residue sidechains.
            n = ligand->get_atom_count();
            int offset = n;
            for (l=0; l<n; l++) ligand->get_atom(l)->stream_pdb_line(pdbdat, 9000+l);
            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Prepared ligand PDB." << endl;
            #endif

            if (waters)
            {
                for (k=0; k<maxh2o; k++)
                {
                    for (l=0; l<3; l++) waters[k]->get_atom(l)->stream_pdb_line(pdbdat, 9000+offset+l+3*k);
                }
            }

            #if _dummy_atoms_for_debug
            if (dummies.size())
            {
                for (k=0; k<dummies.size(); k++)
                {
                    dummies[k].stream_pdb_line(pdbdat, 9900+offset+l+3*k);
                }
            }
            #endif

            sphres = protein->get_residues_can_clash_ligand(reaches_spheroid[nodeno], ligand, pocketcen, size, addl_resno);

            if (flex)
            {
                int en = protein->get_end_resno();
                int resno;
                for (resno = protein->get_start_resno(); resno <= en; resno++)
                {
                    AminoAcid* laa = protein->get_residue(resno);
                    if (!laa) continue;
                    if (!laa->been_flexed)
                    {
                        if (laa->distance_to(ligand) > 5) continue;
                        for (k=0; reaches_spheroid[nodeno][k]; k++)
                        {
                            if (!protein->aa_ptr_in_range(reaches_spheroid[nodeno][k])) continue;
                            if (reaches_spheroid[nodeno][k] == laa) goto _afterall;
                        }
                        continue;
                    }
                    _afterall:
                    n = laa->get_atom_count();
                    for (l=0; l<n; l++)
                    {
                        laa->get_atom(l)->stream_pdb_line(
                            pdbdat,
                            laa->atno_offset+l
                        );
                    }
                }
                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Prepared flex PDBs." << endl;
                #endif
            }

            if (mtlcoords.size())
            {
                for (l=0; l<mtlcoords.size(); l++)
                {
                    mtlcoords[l].mtl->stream_pdb_line(
                        pdbdat,
                        9900+l
                    );
                }
            }

            dr[drcount][nodeno].pdbdat = pdbdat.str();
            if (debug) *debug << "Prepared the PDB strings." << endl;

            if (!nodeno)
            {
                if (dr[drcount][nodeno].ligand_self < -individual_clash_limit)
                {
                    // cout << "Internal ligand energy " << -dr[drcount][nodeno].ligand_self << " out of range." << endl << endl;
                    break;          // Exit nodeno loop.
                }
                // else cout << "Internal ligand energy " << -dr[drcount][nodeno].ligand_self << " satisfactory." << endl << endl;

                if (dr[drcount][nodeno].worst_energy > l_individual_clash_limit)
                {
                    // cout << "Least favorable binding energy " << dr[drcount][nodeno].worst_energy << " out of range." << endl << endl;
                    break;          // Exit nodeno loop.
                }
                // else cout << "Least favorable binding energy " << dr[drcount][nodeno].worst_energy << " satisfactory." << endl << endl;

                if (pose==1) dr[drcount][nodeno].pose = pose;
                else
                {
                    int bestpose = pose;
                    for (i=0; i<drcount; i++)
                    {
                        if ((	differential_dock
                                &&
                                (dr[i][0].kJmol - dr[i][0].ikJmol + dr[i][0].polsat * polar_sat_influence_for_scoring)
                                <
                                (dr[drcount][nodeno].kJmol - dr[drcount][nodeno].ikJmol + dr[drcount][nodeno].polsat * polar_sat_influence_for_scoring)
                            )
                            ||
                            (	!differential_dock
                                &&
                                (dr[i][0].kJmol + dr[i][0].polsat * polar_sat_influence_for_scoring)
                                <
                                (btot + pstot * polar_sat_influence_for_scoring)
                            ))
                        {
                            if (dr[i][0].pose < bestpose || bestpose < 0) bestpose = dr[i][0].pose;
                            dr[i][0].pose++;
                        }
                    }
                    dr[drcount][nodeno].pose = bestpose;
                    // cout << "Around the posie: "; for (i=0; i<=drcount; i++) cout << dr[i][nodeno].pose << " "; cout << endl;
                }
                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Added pose to output array." << endl;
                #endif

                #if _DBG_MAX_CLASHES
                cout << "Pose " << dr[drcount][nodeno].pose << " maxclash " << maxclash << " kJmol " << dr[drcount][nodeno].kJmol << endl;
                #endif
            }

            // For performance reasons, once a path node (including #0) fails to meet the binding energy threshold, discontinue further
            // calculations for this pose.
            if (btot < kJmol_cutoff && !differential_dock)
            {
                drcount++;
                break;
            }
            else if (nodeno == pathnodes) drcount++;
        }	// nodeno loop.
    } // pose loop.
    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Finished poses." << endl;
    #endif

    if (progressbar)
    {
        cout << "\033[A\033[K";
    }

    // Output the dr[][] array in order of increasing pose number.
    cout << endl;
    if (output) *output << endl;

    const float energy_mult = kcal ? _kcal_per_kJ : 1;
    pose = 1;
    for (i=1; i<=poses; i++)
    {
        for (j=0; j<poses; j++)
        {
            protein = &pose_proteins[j];
            ligand = &pose_ligands[j+1];
            if (dr[j][0].pose == i && dr[j][0].pdbdat.length())
            {
                if (differential_dock || dr[j][0].kJmol >= kJmol_cutoff)
                {
                    for (k=0; k<=pathnodes; k++)
                    {
                        // If pathnode is not within kJ/mol cutoff, abandon it and all subsequent pathnodes of the same pose.
                        if (dr[j][k].kJmol < kJmol_cutoff)
                        {
                            cout << "Pose " << pose << " node " << k
                                 << " energy " << -dr[j][k].kJmol*energy_mult
                                 << " is outside of limit; aborting path nodes." << endl;
                            if (output) *output << "Pose " << pose << " node " << k
                                                << " energy " << -dr[j][k].kJmol*energy_mult
                                                << " is outside of limit; aborting path nodes." << endl;
                            break;
                        }

                        /*if (flex && !dr[j][k].pdbdat.length())
                        {
                            cout << "Pose " << j << " node " << k << " is missing." << endl;
                            if (output) *output << "Pose " << j << " node " << k << " is missing." << endl;
                            continue;
                        }*/

                        cout << "Pose: " << pose << endl << "Node: " << k << endl;
                        if (output) *output << "Pose: " << pose << endl << "Node: " << k << endl;

                        dr[j][k].energy_mult = energy_mult;
                        dr[j][k].do_output_colors = do_output_colors;
                        dr[j][k].include_pdb_data = (output == nullptr);
                        cout << dr[j][k];
                        dr[j][k].do_output_colors = false;
                        dr[j][k].include_pdb_data = true;
                        if (output) *output << dr[j][k];


                        if (!k && outpdb.length() && pose <= outpdb_poses)
                        {
                            char protn[64];
                            strcpy(protn, strrchr(protfname, '/')+1);
                            char* dot = strchr(protn, '.');
                            if (dot) *dot = 0;

                            char lign[64];
                            strcpy(lign, strrchr(ligfname, '/')+1);
                            dot = strchr(lign, '.');
                            if (dot) *dot = 0;

                            // std::string temp_pdb_fn = (std::string)"tmp/pose" + std::to_string(j+1) + (std::string)".pdb";
                            std::string out_pdb_fn = std::regex_replace(outpdb, std::regex("[%][p]"), protn);
                            out_pdb_fn = std::regex_replace(out_pdb_fn, std::regex("[%][l]"), lign);
                            out_pdb_fn = std::regex_replace(out_pdb_fn, std::regex("[%][o]"), to_string(pose));

                            // FILE* pftmp = fopen(temp_pdb_fn.c_str(), "rb");
                            FILE* pfout = fopen(out_pdb_fn.c_str(), "wb");
                            if (/*!pftmp ||*/ !pfout)
                            {
                                // if (!pftmp) cout << "Failed to open " << temp_pdb_fn << " for reading." << endl;
                                if (!pfout) cout << "Failed to open " << out_pdb_fn << " for writing." << endl;
                                return -1;
                            }


                            int j1;
                            int n1 = protein->get_end_resno();
                            for (j1=1; j1 <= n1; j1++)
                            {
                                AminoAcid* aa = protein->get_residue(j1);
                                if (aa /* && aa->been_flexed */)
                                {
                                    // tmp_pdb_residue[j+1][j1].restore_state(aa);
                                    #if _dbg_residue_poses
                                    cout << "tmp_pdb_residue[" << (j+1) << "][" << j1 << "].restore_state(" << aa->get_name() << ")" << endl;
                                    #endif
                                }
                            }
                            // tmp_pdb_ligand[j+1].restore_state(ligand);

                            #if recapture_ejected_ligand
                            float r = ligand->get_barycenter().get_3d_distance(nodecens[k]);
                            if (r > size.magnitude()/2) goto _next_pose;
                            #endif

                            protein->save_pdb(pfout, ligand);

                            int atno_offset = protein->last_saved_atom_number;
                            if (waters)
                            {
                                for (j1=0; waters[j1]; j1++)
                                {
                                    tmp_pdb_waters[pose][j1].restore_state(waters[j1]);
                                    waters[j1]->save_pdb(pfout, atno_offset);
                                    atno_offset += waters[j1]->get_atom_count();
                                }
                            }

                            n1 = mtlcoords.size();
                            for (j1=0; j1 < n1; j1++)
                            {
                                mtlcoords[j1].mtl->move(tmp_pdb_metal_locs[pose][j1]);
                                mtlcoords[j1].mtl->save_pdb_line(pfout, ++atno_offset);
                            }

                            protein->end_pdb(pfout);

                            fclose(pfout);
                        }


                        if (!k) found_poses++;
                    }
                    pose++;
                }
                else
                {
                    if (i == 1)
                    {
                        if (kcal)
                        {
                            cout << "No poses found within kcal/mol limit." << endl;
                            if (output) *output << "No poses found within kcal/mol limit." << endl;
                        }
                        else
                        {
                            cout << "No poses found within kJ/mol limit." << endl;
                            if (output) *output << "No poses found within kJ/mol limit." << endl;
                        }
                    }
                    cout << "Exiting." << endl;
                    goto _exitposes;
                }

                break;
            }

            _next_pose:
            ;
        }
    }

_exitposes:
    if (found_poses < poses && triesleft)
    {
        triesleft--;
        goto _try_again;
    }

    cout << found_poses << " pose(s) found." << endl;
    if (output) *output << found_poses << " pose(s) found." << endl;
    if (debug) *debug << found_poses << " pose(s) found." << endl;

    if (met) delete met;

    time_t finished = time(NULL);
    cout << "\nCalculation took: " << (finished-began) << " seconds." << endl;
    if (output) *output << "\nCalculation took: " << (finished-began) << " seconds." << endl;
    if (debug) *debug << "\nCalculation took: " << (finished-began) << " seconds." << endl;

    if (output) output->close();
    if (append_pdb)
    {
        if (output)
        {
            pf = fopen(temp_pdb_file.length() ? temp_pdb_file.c_str() : protfname, "r");
            if (!pf)
            {
                cout << "Error trying to read " << protfname << endl;
                return 0xbadf12e;
            }
            protein->load_pdb(pf);
            fclose(pf);
            FILE* pf = fopen(outfname, "ab");
            fprintf(pf, "\nOriginal PDB:\n");
            protein->save_pdb(pf);
            fclose(pf);
            cout << "PDB appended to output file." << endl;
        }
        else cout << "ERROR: Append PDB can only be used when specifying an output file." << endl;
    }

    if (temp_pdb_file.length()) std::remove(temp_pdb_file.c_str());

    if (debug) debug->close();

    return 0;
}


















