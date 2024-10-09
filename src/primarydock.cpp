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
#include "classes/search.h"
#include "classes/scoring.h"
#include "classes/cavity.h"

using namespace std;

struct AcvBndRot
{
    int resno;
    std::string aname;
    std::string bname;
    Atom* atom1 = nullptr;
    Atom* atom2 = nullptr;
    Bond* bond = nullptr;
    float theta;
};

std::string get_fttl_no_path_or_ext(const char* filename)
{
    char buffer[256];
    strcpy(buffer, filename);
    int i, n = strlen(buffer);
    for (i=n; i>0; i--)
    {
        if (buffer[i] == '.')
        {
            buffer[i] = 0;
            break;
        }
    }
    for (; i>0; i--)
    {
        if (buffer[i] == '/')
        {
            std::string result = &buffer[i+1];
            return result;
        }
    }

    std::string result = filename;
    return result;
}

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
char cvtyfname[256];
Point pocketcen;
std::ofstream *output = NULL;

std::vector<std::string> CEN_buf;
int cenbuf_idx = 0;
std::vector<std::string> pathstrs;
std::vector<std::string> states;

Cavity cvtys[256];
int ncvtys = 0;
std::vector<std::string> dyn_strings;
std::vector<DynamicMotion> dyn_motions;

bool configset=false, protset=false, tplset=false, tprfset=false, ligset=false, ligcmd=false, smset = false, smcmd = false, pktset=false;

bool out_per_res_e = true;
bool out_per_btyp_e = true;
float out_itemized_e_cutoff = 0.01;
bool out_lig_int_e = true;
bool out_bb_pairs = false;
bool out_lig_pol_sat = false;
bool out_prox = false;
bool out_pro_clash = false;
bool out_mc = false;
bool out_vdw_repuls = false;
bool out_pdbdat_lig = true;
bool out_pdbdat_res = true;

Protein* protein;
Protein* ptemplt;
Protein* ptplref;
Protein* metald_prot = nullptr;
int seql = 0;
int addl_resno[256];
const Region* regions;
SCoord region_clashes[85][3];
Molecule* ligand;
std::vector<std::string> isomers;
Molecule** waters = nullptr;
Molecule** owaters = nullptr;
Point ligcen_target;
SCoord path[256];
int pathnodes = 0;				// The pocketcen is the initial node.
int poses = 10;
int iters = 50;
int maxh2o = 0;
int omaxh2o = 0;
bool flex = true;
float kJmol_cutoff = 0.01;
int pose=1, nodeno=0, iter=0;
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

bool softdock = false;
float softness = 0;
std::vector<Region> softrgns;
std::vector<Atom*> softrgpiva;
std::vector<int> softrgn_allowed;
std::vector<float> softrgn_initclash;
std::vector<ResiduePlaceholder> soft_nodelete_start, soft_nodelete_end;
char splash[16384];

AminoAcid*** reaches_spheroid = nullptr;
int sphres = 0;

std::string origbuff = "";
std::string optsecho = "";

// Switch to enable "best-binding" algorithm rather than "tumble spheres" algorithm.
PoseSearchType pdpst = default_search_algorithm;
std::string copyfrom_filename;
char copyfrom_ligname[5] = {0,0,0,0,0};
int copyfrom_resno = -1;
Bond retain_bindings[4];
std::vector<int> center_resnos;
std::vector<int> priority_resnos;

Atom* pivotal_hbond_aaa = nullptr;
Atom* pivotal_hbond_la = nullptr;
float pivotal_hbond_r = 0;

Pose pullaway_undo;
float last_ttl_bb_dist = 0;

std::vector<AcvBndRot> active_bond_rots;
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
        e = -waters[j]->get_intermol_binding(gcfmols).summed();
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

    char* offset = buffer;
    while (*offset >= 'A' && *offset <= 'Z') offset++;

    int retval = 0;
    char* dot = strchr(buffer, '.');
    char* bang = strchr(buffer, '!');
    if (dot)
    {
        *(dot++) = 0;
        int b = atoi(offset);
        int w = atoi(dot);
        int _50 = protein->get_bw50(b);
        if (_50 < 1)
        {
            cout << "Error: unknown BW number " << b << "." << w << ", please ensure PDB file has REMARK 800 SITE BW words." << endl;
            throw 0xbad12e5;
        }
        if (_50 < 1) return 0;
        else retval = _50 + w - 50;
    }
    else retval = atoi(offset);

    AminoAcid* aa = protein->get_residue(retval);
    if (!aa) return 0;

    if (offset == buffer)
    {
        if (bang)
        {
            aa->priority = true;
            priority_resnos.push_back(aa->get_residue_no());
        }
        return retval;
    }

    int i;
    for (i=0; buffer[i] >= 'A' && buffer[i] <= 'Z'; i++)
    {
        if (buffer[i] == aa->get_letter())
        {
            if (bang)
            {
                aa->priority = true;
                priority_resnos.push_back(aa->get_residue_no());
            }
            return retval;
        }
    }

    return bang ? retval : 0;
}

void freeze_bridged_residues()
{
    int i, l;

    if (bridges.size())
    {
        for (i=0; i<bridges.size(); i++)
        {
            int resno1 = interpret_resno(bridges[i].c_str());
            if (!resno1) continue;
            const char* r2 = strchr(bridges[i].c_str(), '|');
            if (!r2) throw 0xbadc0de;
            r2++;
            int resno2 = interpret_resno(r2);
            if (!resno2) continue;
            
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
        if (!resno1) continue;
        const char* r2 = strchr(bridges[i].c_str(), '|');
        if (!r2) throw 0xbadc0de;
        r2++;
        int resno2 = interpret_resno(r2);
        if (!resno2) continue;

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
            float tb = -aa1->get_intermol_binding(aa2).summed();
            cout << "Bridge energy " << tb << " kJ/mol." << endl;
        }
        #endif
    }

    freeze_bridged_residues();
}

void do_pivotal_hbond_rot_and_scoot()
{
    // Rotate ligand so atom faces side chain atom.
    ligand->movability = MOV_ALL;
    float r = pivotal_hbond_aaa->distance_to(pivotal_hbond_la);
    Point cen = ligand->get_barycenter();
    Rotation rot = align_points_3d(pivotal_hbond_la->get_location(), pivotal_hbond_aaa->get_location(), cen);
    LocatedVector lv = rot.v;
    lv.origin = cen;
    ligand->rotate(lv, rot.a);

    #if _dbg_priority_hbond
    cout << "Rotated ligand " << (rot.a*fiftyseven) << "deg." << endl;
    cout << "Atoms were " << r << "Å apart, now " << pivotal_hbond_aaa->distance_to(pivotal_hbond_la) << "Å" << endl;
    #endif

    // Measure the distance between atoms and move the ligand to optimize that distance.
    SCoord scooch = pivotal_hbond_aaa->get_location().subtract(pivotal_hbond_la->get_location());
    if (scooch.r > 2.0)
    {
        scooch.r -= 2.0;
        ligand->move(scooch);

        #if _dbg_priority_hbond
        cout << "Scooched ligand " << scooch.r << "Å." << endl;
        cout << "Ligand was centered at " << cen << "; now " << ligand->get_barycenter() << endl;
        cout << "Atoms are now " << pivotal_hbond_aaa->distance_to(pivotal_hbond_la) << "Å apart" << endl;
        #endif
    }

    // Rotate the ligand about the hbond in order to minimize the intermolecular energy.
    SCoord v = pivotal_hbond_la->get_location().subtract(pivotal_hbond_aaa->get_location());
    lv = v;
    lv.origin = pivotal_hbond_aaa->get_location();
    float theta = 0, th, step = M_PI/50, clash;
    AminoAcid* reaches[SPHREACH_MAX+4];
    protein->get_residues_can_clash_ligand(reaches, ligand, ligand->get_barycenter(), Point(5,5,5), nullptr);
    for (th=0; th<circle; th += step)
    {
        float c = ligand->get_intermol_clashes((Molecule**)reaches);

        if (!th || c < clash)
        {
            clash = c;
            theta = th;
        }

        ligand->rotate(lv, step);
    }

    ligand->rotate(lv, theta);
}

void output_iter(int iter, Molecule** mols)
{
    std::string itersfname = (std::string)"tmp/" + (std::string)"_iters.dock";
    int i, liter = iter + movie_offset;
    FILE* fp = fopen(itersfname.c_str(), ((liter == 0 && pose == 1) ? "wb" : "ab") );
    if (fp)
    {
        if (!liter && (pose == 1))
        {
            fprintf(fp, "PDB file: %s\n", protfname);
        }
        fprintf(fp, "Pose: %d\nNode: %d\n\n", pose, liter);
        int foff = 0;

        DockResult ldr(protein, ligand, size, nullptr, pose, waters);
        ldr.include_pdb_data = false;
        ldr.display_clash_atoms = true;
        std::stringstream stst;
        stst << ldr;
        fprintf(fp, "%s\n", stst.str().c_str());

        fprintf(fp, "\nPDBDAT:\n");

        for (i=0; reaches_spheroid[nodeno][i]; i++)
        {
            reaches_spheroid[nodeno][i]->save_pdb(fp, foff);
            foff += reaches_spheroid[nodeno][i]->get_atom_count();
        }

        for (i=0; mols[i]; i++)
        {
            if (mols[i]->is_residue()) continue;
            mols[i]->save_pdb(fp, foff, false);
            foff += mols[i]->get_atom_count();
        }

        protein->end_pdb(fp);

        fclose(fp);
    }
}

Pose iter_best_pose[1000];
float iter_best_bind;
void iteration_callback(int iter, Molecule** mols)
{
    // if (kJmol_cutoff > 0 && ligand->lastbind >= kJmol_cutoff) iter = (iters-1);
    int i, j, l, n;
    float f = 0;

    if (pdpst == pst_constrained)
    {
        Point resca = cs_res[cs_idx]->get_CA_location();
        Point agcen = cs_lag[cs_idx]->get_center();
        Point barcn = ligand->get_barycenter();
        Rotation csrot = align_points_3d(agcen, resca, barcn);
        csrot.a *= cs_ligand_rotation;
        Pose cswas(ligand);
        cswas.copy_state(ligand);
        Interaction e1 = ligand->get_intermol_binding(reinterpret_cast<Molecule**>(reaches_spheroid[nodeno]));
        ligand->rotate(&csrot.v, csrot.a, true);
        Interaction e2 = ligand->get_intermol_binding(reinterpret_cast<Molecule**>(reaches_spheroid[nodeno]));
        if (!e2.improved(e1)) cswas.restore_state(ligand);

        Atom *ra, *la;
        cs_res[cs_idx]->mutual_closest_atoms(ligand, &ra, &la);
        SCoord r = ra->get_location().subtract(la->get_location());
        r.r -= 2.5;
        r.r *= frand(0.333,0.666);
        cswas.copy_state(ligand);
        e1 = ligand->get_intermol_binding(reinterpret_cast<Molecule**>(reaches_spheroid[nodeno]));
        ligand->move(r);
        e2 = ligand->get_intermol_binding(reinterpret_cast<Molecule**>(reaches_spheroid[nodeno]));
        if (!e2.improved(e1)) cswas.restore_state(ligand);
    }

    // Initialization for best-iteration saving for pose output.
    if (iter == 1)
    {
        for (l=0; mols[l]; l++) iter_best_pose[l].copy_state(mols[l]);
        iter_best_bind = 0;
    }

    // Stochastically force flexion on some side chains that get clashes.
    for (l=0; mols[l]; l++)
    {
        Interaction li = mols[l]->get_intermol_binding(mols);
        float lf = li.summed(), lc = li.repulsive, ptnl = mols[l]->get_intermol_potential(mols);
        f += lf;

        int lres = mols[l]->is_residue();
        float prob = 1;
        if (lres)
        {
            AminoAcid* aa = protein->get_residue(lres);
            if (aa) prob = aa->get_aa_definition()->flexion_probability;
        }

        if (flex
            // && iter < 5
            && (lf < -10 || lf < 0.1 * ptnl || (lc > 0 && lf < 5))
            && mols[l]->movability == MOV_FLXDESEL
            && lres
            && frand(0,1) < flexion_probability_multiplier * prob
            )
        {
            mols[l]->movability = MOV_FORCEFLEX;
        }
    }

    if (f > iter_best_bind)
    {
        for (l=0; mols[l]; l++) iter_best_pose[l].copy_state(mols[l]);
        iter_best_bind = f;
    }

    // Attempt to connect hydrogen bonds to ligand.
    if (flex)
    {
        n = ligand->get_atom_count();
        for (l=0; mols[l]; l++)
        {
            if (!mols[l]->is_residue()) continue;
            if (fabs(mols[l]->hydrophilicity()) < hydrophilicity_cutoff) continue;
            if (mols[l]->movability & MOV_PINNED) continue;

            AminoAcid* hbaa = reinterpret_cast<AminoAcid*>(mols[l]);
            Atom* reach = hbaa->get_reach_atom();
            if (!reach) continue;
            if (fabs(reach->is_polar()) < hydrophilicity_cutoff) continue;
            Bond* rapb = reach->get_bond_by_idx(0);
            if (!rapb) continue;
            Atom* prev = rapb->atom2;
            if (!prev) continue;
            Atom* CB = hbaa->get_atom("CB");
            if (!CB) continue;
            int reachz = reach->get_Z();
            Atom* target = nullptr;
            float nearest = Avogadro;

            for (i=0; i<n; i++)
            {
                Atom* la = ligand->get_atom(i);
                if (!la) continue;
                if (fabs(la->is_polar()) < hydrophilicity_cutoff) continue;
                if (reachz > 1 && la->get_Z() > 1) continue;
                float rCB = la->distance_to(CB);
                if (rCB > _DEFAULT_INTERA_R_CUTOFF + hbaa->get_reach()) continue;
                float rCA = la->get_location().get_3d_distance(hbaa->get_CA_location());
                if (rCB > rCA) continue;
                float rp = la->distance_to(prev);
                if (rp > _DEFAULT_INTERA_R_CUTOFF) continue;
                float r = fmin(rp, rCB/1.25);
                if (r < nearest)
                {
                    nearest = r;
                    target = la;
                }
            }

            if (target)
            {
                if (reachz == 1 && target->get_Z() == 1)
                {
                    if (rand() & 1) reach = reach->get_bond_by_idx(0)->atom2;
                    else target = target->get_bond_by_idx(0)->atom2;
                }

                if (target->distance_to(reach) < 3) continue;

                Point pttgt = target->get_location();
                Pose hbwas(mols[l]);
                hbwas.copy_state(mols[l]);
                Interaction before = ligand->get_intermol_binding(mols);
                hbaa->conform_atom_to_location(reach->name, pttgt, 10, 2);
                Interaction after = ligand->get_intermol_binding(mols);
                if (!after.improved(before)) hbwas.restore_state(mols[l]);
            }
        }
    }

    #if use_best_binding_iteration
    if (iter == iters && iter_best_bind > 0)
    {
        for (l=0; mols[l]; l++) iter_best_pose[l].restore_state(mols[l]);
    }
    #endif

    if (pivotal_hbond_aaa && pivotal_hbond_la) do_pivotal_hbond_rot_and_scoot();

    Point bary = ligand->get_barycenter();

    int ac = ligand->get_atom_count();
    float bbest = 0;
    Atom *atom1, *atom2;

    float progress = (float)iter / iters;

    n = softrgns.size();
    if (n && iter>8) for (i=0; i<n; i++)
    {
        Point softpush(0,0,0);
        Point foravg[softrgns[i].end-softrgns[i].start+16];
        float pushmax = 0;
        l = 0;

        for (j=softrgns[i].start; j<=softrgns[i].end; j++)
        {
            AminoAcid* aa = protein->get_residue(j);
            if (!aa) continue;
            float c = aa->get_intermol_clashes(ligand);
            if (c > clash_limit_per_aa*10)
            {
                Point A = aa->get_CA_location();
                foravg[l++] = A;
                SCoord AB = A.subtract(ligand->get_nearest_atom(A)->get_location());
                AB.r = pow(c, 1.0/3);
                softpush = softpush.add(AB);
                if (AB.r > pushmax) pushmax = AB.r;
            }
        }

        if (l)
        {
            SCoord AB = softpush;
            AB.r = pushmax;
            Point A = average_of_points(foravg, l);
            Point B = A.add(AB);
            Point C = softrgpiva[i]->get_location();
            Rotation rot = align_points_3d(A, B, C);
            rot.a = fmin(rot.a, 0.1*softness*fiftyseventh);

            float c1;
            if (i >= softrgn_initclash.size())
            {
                c1 = protein->get_internal_clashes(softrgns[i].start, softrgns[i].end);
                softrgn_initclash.push_back(c1);
            }
            else c1 = softrgn_initclash[i];

            protein->rotate_piece(softrgns[i].start, softrgns[i].end, C, rot.v, rot.a);
            float c2 = protein->get_internal_clashes(softrgns[i].start, softrgns[i].end, true, 20);
            if (c2 > c1 + clash_limit_per_aa*2) protein->rotate_piece(softrgns[i].start, softrgns[i].end, C, rot.v, -rot.a);
        }
    }

    if (!iter) goto _oei;
    if (iter == (iters-1)) goto _oei;

    if (pdpst == pst_best_binding && ligand_groups[0].atct)
    {
        l = 1;
        Point agcen = global_pairs[l]->ag->get_center();
        Atom* scgna = global_pairs[l]->scg->get_nearest_atom(agcen);
        float r = global_pairs[l]->ag->distance_to(scgna->get_location());
        if (r > 2.5)
        {
            Pose was(ligand);
            was.copy_state(ligand);
            ligand->conform_atom_to_location(global_pairs[l]->ag->atoms[0]->name, scgna->get_location(), 10, frand(2, r));
            if (was.total_atom_motions() > 3.5*ligand->get_heavy_atom_count()) was.restore_state(ligand);
        }
    }

    #if enforce_no_bb_pullaway
    if (pdpst == pst_best_binding && ligand_groups[0].atct)
    {
        float ttl_bb_dist = 0;
        for (l=0; l<3; l++)
        {
            if (global_pairs.size() > l)
            {
                float r = global_pairs[l]->ag->distance_to(global_pairs[l]->scg->get_center()); //  ligand_groups[l].distance_to(sc_groups[l].get_center());
                float r1 = r;
                if (r < 2.5) r = 2.5;
                if (r > _INTERA_R_CUTOFF) r = _INTERA_R_CUTOFF;
                ttl_bb_dist += r * (1.0 + 1.0 / (l+1));
                #if _dbg_bb_pullaway
                cout << pose << ":" << iter << ": Ligand atoms ";
                for (i=0; i<global_pairs[l]->ag->atct; i++) cout << global_pairs[l]->ag->atoms[i]->name << " ";
                cout << "are " << r1 << " A from residues";
                for (i=0; i<global_pairs[l]->scg->aminos.size(); i++) cout << " " << global_pairs[l]->scg->aminos[i]->get_3letter() << global_pairs[l]->scg->aminos[i]->get_residue_no();
                cout << "." << endl;
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
    
    bary = ligand->get_barycenter();

    for (i=0; i < ac; i++)
    {
        if (ligand->get_atom(i)->strongest_bind_energy > bbest)
        {
            atom1 = ligand->get_atom(i);
            bbest = atom1->strongest_bind_energy;
            atom2 = atom1->strongest_bind_atom;
        }
    }

    if (bbest >= 15)
    {
        bary = ligand->get_barycenter();
    }
    else
    {
        #if allow_drift

        Pose predrift(ligand);
        Interaction predrift_binding = ligand->get_intermol_binding(mols);

        if (predrift_binding.summed() < -drift_energy_threshold)
        {
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
        }

        #endif
    }

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

    _oei:
    ;

    #if recapture_ejected_ligand
    Point lig_center = ligand->get_barycenter();
    float r = lig_center.get_3d_distance(ligcen_target);
    float recapture_distance = size.magnitude() / 2;
    if (r >= recapture_distance) ligand->recenter(ligcen_target);
    #endif

    if (output_each_iter) output_iter(iter, mols);
}

int spinchr = 0;
float hueoffset = 0;
void update_progressbar(float percentage)
{
    percentage = percentage/poses + (float)(pose-1)*100.0/poses;
    if (percentage > 100) percentage = 100;
    cout << "\033[A|";
    int i;
    for (i=0; i<100; i++)
    {
        float cmpi = i;
        if (cmpi <= percentage)
        {
            float h = M_PI*2 * cmpi / 46 + hueoffset;
            int r, g, b;
            r =  96 +  24 * sin(h-0.333);
            g = 128 +  26 * sin(h+0.333);
            b = 224 +  31 * sin(h);
            colorrgb(r, g, b);
            cout << "\u2593";
            colorless();
        }
        else cout << "\u2591";
    }
    cout << ("|/-\\")[spinchr] << " " << (int)percentage << "%.               " << endl;
    spinchr++;
    if (spinchr >= 4) spinchr = 0;
    hueoffset += 0.3;
}

void erase_progressbar()
{
    cout << "\033[A\033[K";
}

Point pocketcen_from_config_words(char** words, Point* old_pocketcen)
{
    int i=1;
    Point local_pocketcen;
    center_resnos.clear();
    priority_resnos.clear();
    if (!strcmp(words[i], "RES"))
    {
        i++;
        for (; words[i]; i++)
        {
            int j = interpret_resno(words[i]);
            if (!j) continue;

            center_resnos.push_back(j);
        }

        int sz = center_resnos.size(), div=0;
        Point foravg[sz + 2];
        for (i=0; i<sz; i++)
        {
            #if pocketcen_from_reach_atoms
            AminoAcid* aa = protein->get_residue(center_resnos[i]);
            if (aa)
            {
                foravg[i] = aa->get_reach_atom_location();
                div++;
            }
            #else
            foravg[div++] = protein->get_atom_location(center_resnos[i], "CA");
            #endif
        }

        return average_of_points(foravg, div?:1);
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
        CEN_buf.push_back(origbuff);
        optsecho = (std::string)"Center " + (std::string)origbuff;
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
        ligset = true;
        if (ligcmd) smset = smcmd;
        return 1;
    }
    else if (!strcmp(words[0], "ISO"))
    {
        isomers.push_back(words[1]);
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
            if (words[i][0] == '#') break;
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
    else if (!strcmp(words[0], "PERRES"))
    {
        out_per_res_e = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "PERBTYP"))
    {
        out_per_btyp_e = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "ELIMITEM"))
    {
        out_itemized_e_cutoff = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "LIGINTE"))
    {
        out_lig_int_e = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "OUTBBP"))
    {
        out_bb_pairs = words[1] ? atoi(words[1]) : true;
        return 1;
    }
    else if (!strcmp(words[0], "OUTLPS"))
    {
        out_lig_pol_sat = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "OUTPROX"))
    {
        out_prox = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "OUTPCLSH"))
    {
        out_pro_clash = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "OUTMC"))
    {
        out_mc = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "OUTVDWR"))
    {
        out_vdw_repuls = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "OUTPDBL"))
    {
        out_pdbdat_lig = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "OUTPDBR"))
    {
        out_pdbdat_res = atoi(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "POSE"))
    {
        poses = atoi(words[1]);
        optsecho = "Number of poses: " + to_string(poses);
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
            pdpst = pst_best_binding;
            return 1;
        }
        else if (!strcmp(words[1], "TS"))
        {
            pdpst = pst_tumble_spheres;
            return 1;
        }
        else if (!strcmp(words[1], "CS"))
        {
            pdpst = pst_constrained;
            return 1;
        }
        else if (!strcmp(words[1], "CF"))
        {
            pdpst = pst_cavity_fit;
            return 1;
        }
        else if (!strcmp(words[1], "CP"))
        {
            int lf = 1;
            pdpst = pst_copyfrom;
            if (!words[2])
            {
                cout << "ERROR: Search mode CP without source file." << endl;
                throw 0xbad19b07;
            }
            copyfrom_filename = words[2];
            if (words[3])
            {
                if (strlen(words[3]) > 3) words[3][3] = 0;
                strcpy(copyfrom_ligname, words[3]);
                lf++;
                if (words[4])
                {
                    copyfrom_resno = atoi(words[4]);
                    lf++;
                }
            }
            return lf;
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
        softdock = true;
        softness = 1;

        for (i=1; words[i]; i++)
        {
            if (strchr(words[i], '.'))
            {
                softness = atof(words[i]);
                #if _dbg_soft
                cout << "Softness: " << softness << endl;
                #endif
            }
            else
            {
                int j = atoi(words[i]);
                if (j) softrgn_allowed.push_back(j);
                #if _dbg_soft
                cout << "Allowed soft region: " << j << endl;
                #endif
            }
        }
    }
    else if (!strcmp(words[0], "NODEL"))
    {
        ResiduePlaceholder rpsr, rper;
        rpsr.set(words[1]);
        rper.set(words[2]);
        soft_nodelete_start.push_back(rpsr);
        soft_nodelete_end.push_back(rper);
        return 2;
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
    else if (!strcmp(words[0], "CAVS"))
    {
        cavity_stuffing = atof(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "CFLEE"))
    {
        clash_fleeing = atof(words[1]);
        return 1;
    }
    else if (!strcmp(words[0], "VCVTY"))
    {
        strcpy(cvtyfname, words[1]);
        return 1;
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

            delete[] words;
        }
        buffer[0] = 0;
    }
}

void prepare_acv_bond_rots()
{
    int i;
    if (active_bond_rots.size())
    {
        for (i=0; i<active_bond_rots.size(); i++)
        {
            active_bond_rots[i].atom1 = protein->get_atom(active_bond_rots[i].resno, active_bond_rots[i].aname.c_str());
            active_bond_rots[i].atom2 = protein->get_atom(active_bond_rots[i].resno, active_bond_rots[i].bname.c_str());

            if (!active_bond_rots[i].atom1) cout << "WARNING: " << active_bond_rots[i].resno << ":" << active_bond_rots[i].aname
                << " not found in protein!" << endl;
            if (!active_bond_rots[i].atom2) cout << "WARNING: " << active_bond_rots[i].resno << ":" << active_bond_rots[i].bname
                << " not found in protein!" << endl;

            if (active_bond_rots[i].atom1 && active_bond_rots[i].atom2)
            {
                active_bond_rots[i].bond = active_bond_rots[i].atom1->get_bond_between(active_bond_rots[i].atom2);
                #if _debug_active_bond_rot
                active_bond_rots[i].bond->echo_on_rotate = true;
                #endif
            }
            else active_bond_rots[i].bond = nullptr;
        }
    }
}

void attempt_priority_hbond()
{
    int i, j, m, n;

    n = priority_resnos.size();
    for (i=0; i<n; i++)
    {
        AminoAcid* aa = protein->get_residue(priority_resnos[i]);
        if (fabs(aa->hydrophilicity()) >= hydrophilicity_cutoff)
        {
            // Find atom of side chain capable of hydrogen bond.
            Atom* aaa = aa->capable_of_inter(hbond);
            if (!aaa) continue;

            #if _dbg_priority_hbond
            cout << aa->get_name() << ":" << aaa->name << " found." << endl;
            #endif

            // Find nearest atom of ligand capable of hbond.
            m = ligand->get_atom_count();
            float r = Avogadro;
            Atom* la = nullptr;
            for (j=0; j<m; j++)
            {
                if (frand(0, 1) < 0.29) continue;
                Atom* a1 = ligand->get_atom(j);
                if (fabs(a1->is_polar()) >= hydrophilicity_cutoff)
                {
                    float r1 = a1->distance_to(aaa);
                    if (r1 < r)
                    {
                        r = r1;
                        la = a1;
                    }
                }
            }
            if (!la) continue;

            #if _dbg_priority_hbond
            cout << "ligand:" << la->name << " found." << endl;
            #endif

            // Ensure one atom is a donor and the other an acceptor. If one donor and one acceptor cannot be found, skip.
            if (sgn(aaa->is_polar()) == sgn(la->is_polar()))
            {
                if (aaa->get_Z() == 1)
                {
                    Atom* heavy = aaa->get_bond_by_idx(0)->atom1;
                    if (heavy->get_family() == CHALCOGEN) aaa = heavy;
                }
                else
                {
                    Atom* hyd = aaa->is_bonded_to("H");
                    if (hyd) aaa = hyd;
                }
            }
            if (sgn(aaa->is_polar()) == sgn(la->is_polar()))
            {
                if (la->get_Z() == 1)
                {
                    Atom* heavy = la->get_bond_by_idx(0)->atom1;
                    if (heavy->get_family() == CHALCOGEN) la = heavy;
                    else if (heavy->get_family() == PNICTOGEN && !heavy->get_charge()) la = heavy;
                }
                else
                {
                    Atom* hyd = la->is_bonded_to("H");
                    if (hyd) la = hyd;
                }
            }
            if (sgn(aaa->is_polar()) == sgn(la->is_polar())) continue;

            #if _dbg_priority_hbond
            cout << aa->get_name() << ":" << aaa->name << " and ligand:" << la->name << endl;
            #endif

            // Pin the ligand so it can only rotate about the hbond atom.
            pivotal_hbond_aaa = aaa;
            pivotal_hbond_la = la;
            pivotal_hbond_r = 2.0;

            do_pivotal_hbond_rot_and_scoot();

            #if _dbg_priority_hbond
            cout << endl;
            #endif
            return;
        }
    }
}

void choose_cen_buf()
{
    int i, n;

    n = protein->get_end_resno();
    for (i=1; i<=n; i++)
    {
        AminoAcid* aa = protein->get_residue(i);
        if (aa) aa->priority = false;
    }

    priority_resnos.clear();

    n = CEN_buf.size();
    if (pose <= n) cenbuf_idx = pose-1;
    else
    {
        for (i=0; i<10; i++)
        {
            cenbuf_idx = rand() % n;
            if (strchr(CEN_buf[cenbuf_idx].c_str(), '!')) return;
            else if (frand(0,1) < 0.1) return;
        }
    }
}

void apply_protein_specific_settings(Protein* p)
{
    int i, j, n;
    choose_cen_buf();

    char buffer[1024];
    strcpy(buffer, CEN_buf[cenbuf_idx].c_str());
    char** words = chop_spaced_words(buffer);
    pocketcen = pocketcen_from_config_words(words, nullptr);
    loneliest = protein->find_loneliest_point(pocketcen, size);
    delete[] words;

    if (n = priority_resnos.size()) for (i=0; i<n; i++)
    {
        AminoAcid* aa = protein->get_residue(priority_resnos[i]);
        if (aa) aa->priority = true;
    }

    if (!p->get_metals_count() && metald_prot) p->copy_mcoords(metald_prot);

    j=0;

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

        delete[] words;
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

    if (softdock)
    {
        int sndn = soft_nodelete_start.size();
        for (i=0; i<sndn; i++)
        {
            if (!soft_nodelete_start[i].resno) soft_nodelete_start[i].resolve_resno(protein);
            if (!soft_nodelete_end[i].resno) soft_nodelete_end[i].resolve_resno(protein);
        }

        n = protein->get_end_resno();
        for (i=1; i<=n; i++)
        {
            bool helixed = false;
            AminoAcid* aa = protein->get_residue(i);
            if (!aa) continue;
            if (aa->priority) continue;

            // TODO: This should not be hard coded to 6.48-6.59, and it should not depend on the protein having BW numbers.
            BallesterosWeinstein bw = protein->get_bw_from_resno(aa->is_residue());
            if (bw.helix_no == 6 && bw.member_no >= 48 && bw.member_no <= 59) continue;

            std::string rgname = (bw.helix_no < 8 ? "TMR" : (bw.helix_no == 8 ? "HXR" : 
                ((bw.helix_no == 23 || bw.helix_no == 45 || bw.helix_no == 67) ? "EXR" : "CYT")));
            if (bw.helix_no < 10) rgname += std::to_string(bw.helix_no);
            else rgname += std::to_string((int)(bw.helix_no+10)/20);
            if (i >= protein->get_region_start(rgname.c_str()) && i <= protein->get_region_end(rgname.c_str())) continue;

            bool snfound = false;
            for (j=0; j<sndn; j++)
            {
                if (i >= soft_nodelete_start[j].resno && i <= soft_nodelete_end[j].resno)
                {
                    snfound = true;
                    break;
                }
            }
            if (snfound) continue;

            if (aa->is_alpha_helix()) helixed = true;
            else for (j=i-2; j<=i+2; j++)
            {
                if (j < 1) continue;
                if (j == i) continue;
                if (j > n) break;
                aa = protein->get_residue(j);
                if (!aa) continue;
                if (aa->is_alpha_helix() || aa->priority) helixed = true;
                if (helixed) break;
            }

            if (!helixed) protein->delete_residue(i);
        }

        if (!softrgns.size()) for (i=1; i<=n; i++)
        {
            softrgn_initclash.clear();
            if (protein->get_residue(i))
            {
                Region r;
                r.start = i;
                BallesterosWeinstein bwi = protein->get_bw_from_resno(i);
                for (j=i+1; j<=n; j++)
                {
                    BallesterosWeinstein bwj = protein->get_bw_from_resno(j);
                    if (j==n || !protein->get_residue(j) || bwi.helix_no != bwj.helix_no)
                    {
                        r.end = j-1;
                        if (r.end > r.start+2)
                        {
                            bool allowed = true;
                            int l = (i+j-1)/2;
                            BallesterosWeinstein bw = protein->get_bw_from_resno(l);

                            int m;
                            if (m = softrgn_allowed.size())         // assignment not comparison
                            {
                                allowed = false;
                                for (l=0; l<m; l++) if (bw.helix_no == softrgn_allowed[l]) allowed = true;
                            }

                            if (allowed)
                            {
                                Atom* a = protein->region_pivot_atom(r);
                                if (a)
                                {
                                    softrgns.push_back(r);
                                    softrgpiva.push_back(a);
                                    #if _dbg_soft
                                    cout << "Region " << r.start << "-" << r.end << " bw " << bw.helix_no << " pivots about " << a->residue << ":" << a->name << endl;
                                    #endif
                                }
                                #if _dbg_soft
                                else
                                {
                                    cout << "Region " << r.start << "-" << r.end << " bw " << bw.helix_no << " cannot soft pivot." << endl;
                                }
                                #endif
                            }
                            #if _dbg_soft
                            else
                            {
                                cout << "Region " << r.start << "-" << r.end << " bw " << bw.helix_no << " is not allowed to soft pivot." << endl;
                            }
                            #endif
                        }
                        i=j;
                        break;
                    }
                }
            }
        }
    }
}

int main(int argc, char** argv)
{
    strcpy(splash, "\n                                                                                      __       ____  \npppp                                            ddd                               ,-_/  `-_--_/    \\  \np   p         i                                 d  d                 k            )             /  (__   \np   p                                           d   d                k           )   ()   \\__   \\__   )   \npppp  r rrr  iii  mmm mm   aaaa   r rrr  y   y  d   d   ooo    ccc   k   k      /      \\__/  \\__/    /  \np     rr      i   m  m  m      a  rr     y   y  d   d  o   o  c   c  k  k      (       /  \\__/  \\   (  \np     r       i   m  m  m   aaaa  r      y   y  d   d  o   o  c      blm        \\    ()        _     )  \np     r       i   m  m  m  a   a  r      y   y  d  d   o   o  c   c  k  k        )     __     / \\   /  \np     r      iii  m  m  m   aaaa  r       yyyy  ddd     ooo    ccc   k   k       \\____/  `---'   \\__)  \n                                             y\nprimaryodors.org molecular docker      yyyyyy\n\n");
    char buffer[65536];
    int i, j;

    _momentum_rad_ceiling = fiftyseventh * 5;

    for (i=0; i<65536; i++) buffer[i] = 0;

    for (i=0; i<256; i++)
        configfname[i] = protfname[i] = protafname[i] = ligfname[i] = cvtyfname[i] = 0;

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

    cout << splash << endl << endl;

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

    if (out_bb_pairs && pdpst != pst_best_binding)
    {
        cout << "Warning: OUTBBP is enabled but the search mode is not the best-binding algorithm." << endl;
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

    if (strlen(cvtyfname))
    {
        if (!file_exists(cvtyfname))
        {
            cout << "ERROR file not found: " << cvtyfname << endl;
            return -7;
        }
        FILE* fp = fopen(cvtyfname, "rb");
        if (!fp)
        {
            cout << "FAILED to open " << cvtyfname << " for reading." << endl;
            return -7;
        }

        cout << "Reading " << cvtyfname << "..." << endl;
        char buffer[1024];
        while (!feof(fp))
        {
            fgets(buffer, 1022, fp);
            CPartial cp;
            int cno = cp.from_cvty_line(buffer);
            cvtys[cno].add_partial(cp);
            if (cno+1 > ncvtys) ncvtys = cno+1;
        }
        fclose(fp);
        cout << "Read " << ncvtys << " cavities." << endl;
    }

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

    #if _dbg_A100
    float init_A100 = protein->A100();
    #endif

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

        if (temp_pdb_file.length()) std::remove(temp_pdb_file.c_str());
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

        if (temp_pdb_file.length()) std::remove(temp_pdb_file.c_str());
        temp_pdb_file = (std::string)"tmp/" + std::to_string(pid) + (std::string)"_hydro.pdb";

        pf = fopen(temp_pdb_file.c_str(), "wb");
        protein->save_pdb(pf);
        fclose(pf);
    }

    int l, j1, i2, miter;

    if (mtlcoords.size())
    {
        protein->pocketcen = pocketcen;
        mtlcoords = protein->coordinate_metal(mtlcoords);
        metald_prot = protein;

        if (temp_pdb_file.length()) std::remove(temp_pdb_file.c_str());
        temp_pdb_file = (std::string)"tmp/" + std::to_string(pid) + (std::string)"_metal.pdb";

        pf = fopen(temp_pdb_file.c_str(), "wb");
        protein->save_pdb(pf);
        fclose(pf);
    }

    if (bridges.size())
    {
        reconnect_bridges();

        if (temp_pdb_file.length()) std::remove(temp_pdb_file.c_str());
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

    if (!CEN_buf.size())
    {
        cout << "Error: no binding pocket centers defined." << endl;
        return 0xbadb19d;
    }

    #if pocketcen_is_loneliest
    pocketcen = loneliest;
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

    l=0;
    addl_resno[l] = 0;

    // Load the ligand or return an error.
    // Molecule m(ligfname);
    // ligand = &m;
    for (l=0; l<=poses; l++)
    {
        ligand = &pose_ligands[l];
        std::string ligname = get_fttl_no_path_or_ext(ligfname);
        std::string lligfname = ligfname;
        ligand->set_name(ligname.c_str());
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
            if (isomers.size())
            {
                i = rand() % isomers.size();
                ligname = get_fttl_no_path_or_ext(isomers[i].c_str());
                ligand->set_name(ligname.c_str());
                lligfname = isomers[i].c_str();
            }
            pf = fopen(lligfname.c_str(), "rb");
            if (!pf)
            {
                cout << "Error trying to read " << lligfname << endl;
                return 0xbadf12e;
            }
            wgaf = fread(buffer, 1, 65535, pf);
            fclose(pf);
            ligand->from_sdf(buffer);
            break;

        case 'p':
        case 'P':
            if (isomers.size())
            {
                i = rand() % isomers.size();
                ligname = get_fttl_no_path_or_ext(isomers[i].c_str());
                ligand->set_name(ligname.c_str());
                lligfname = isomers[i].c_str();
            }
            pf = fopen(lligfname.c_str(), "rb");
            if (!pf)
            {
                cout << "Error trying to read " << lligfname << endl;
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
    cout << endl;
    if (output) *output << endl;

    i = poses*(triesleft+1)+8;
    j = pathnodes+2;
    DockResult dr[i][j];

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
    if (progressbar) cout << endl;

    reaches_spheroid = new AminoAcid**[pathnodes+2];
    for (i=0; i<=pathnodes; i++) reaches_spheroid[i] = new AminoAcid*[SPHREACH_MAX+4];

    found_poses = 0;
    int wrote_acvmx = -1, wrote_acvmr = -1;
    float l_atom_clash_limit = clash_limit_per_atom; // - kJmol_cutoff;

    std::vector<std::shared_ptr<AtomGroup>> lagc;
    if (pdpst == pst_constrained || pdpst == pst_cavity_fit)
    {
        ligand = &pose_ligands[1];
        lagc = AtomGroup::get_potential_ligand_groups(ligand, mtlcoords.size() > 0);
        agqty = lagc.size();
        if (agqty > MAX_CS_RES-2) agqty = MAX_CS_RES-2;
        for (i=0; i<agqty; i++)
            agc[i] = lagc.at(i).get();

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
                        aa->priority = true;
                        priority_resnos.push_back(aa->get_residue_no());
                    }
                }
            }
        }

        Search::prepare_constrained_search(protein, ligand, pocketcen);
    }

_try_again:
    srand(time(NULL));
    Point nodecens[pathnodes+1];

    /////////////////////////////////////////////////////////////////////////////////
    // Main loop.
    /////////////////////////////////////////////////////////////////////////////////

    regions = protein->get_regions();
    for (i=0; regions[i].start; i++)
    {
        region_clashes[i][0] = region_clashes[i][1] = region_clashes[i][2] = SCoord(0,0,0);
    }

    n = mtlcoords.size();
    Point metal_initlocs[n+4];
    for (i=0; i<n; i++)
    {
        metal_initlocs[i] = mtlcoords[i].mtl->get_location();
    }

    float best_energy = 0, best_acc_energy = 0, best_worst_clash = 0;
    for (pose = 1; pose <= poses; pose++)
    {
        ligand = &pose_ligands[pose];
        ligand->movability = MOV_ALL;

        if (agqty)
        {
            for (i=0; i<agqty; i++) agc[i]->update_atom_pointers(ligand);
        }

        last_ttl_bb_dist = 0;
        ligand->minimize_internal_clashes();
        float lig_min_int_clsh = ligand->get_internal_clashes();
        // if (frand(0,1) < 0.666) ligand->crumple(frand(0, hexagonal));

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

        n = mtlcoords.size();
        for (i=0; i<n; i++)
        {
            mtlcoords[i].mtl->move(metal_initlocs[i]);
        }

        freeze_bridged_residues();

        ligand->recenter(pocketcen);
        // cout << "Centered ligand at " << pocketcen << endl << endl << flush;

        if (pdpst == pst_tumble_spheres)
        {
            Search::do_tumble_spheres(protein, ligand, pocketcen);
            attempt_priority_hbond();
            pocketcen = ligand->get_barycenter();

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

            allow_ligand_360_tumble = nodes_no_ligand_360_tumble && pdpst != pst_best_binding;
            allow_ligand_360_flex   = nodes_no_ligand_360_flex;

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

            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Pose " << pose << endl << "Node " << nodeno << endl;
            #endif
            if (nodeno)
            {
                for (i=0; i<states.size(); i++)
                {
                    strcpy(buffer, states[i].c_str());
                    char** words = chop_spaced_words(buffer);
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
                char** words = chop_spaced_words(buffer);
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
            
            if (pdpst == pst_tumble_spheres && (!prevent_ligand_360_on_activate))
            {
                Search::do_tumble_spheres(protein, ligand, ligcen_target);
                attempt_priority_hbond();
            }
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

            #if _dbg_groupsel
            cout << "Priority resnos: ";
            #endif

            if (n = priority_resnos.size()) for (i=0; i<n; i++)
            {
                AminoAcid* aa = protein->get_residue(priority_resnos[i]);
                if (aa) aa->priority = true;
                #if _dbg_groupsel
                cout << *aa << " ";
                #endif
            }

            #if _dbg_groupsel
            cout << endl << endl;
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

            #if _dbg_groupsel
            cout << "Candidate binding residues: ";
            for (i=0; i<sphres; i++)
            {
                cout << *reaches_spheroid[nodeno][i];
                if (reaches_spheroid[nodeno][i]->priority) cout << "!";
                cout << " ";
            }
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
                    alignment_distance[l]=0;
                }

                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Initialize null AA pointer." << endl;
                #endif

                protein->set_conditional_basicities();
                if (pdpst == pst_best_binding)
                {
                    Search::do_best_binding(protein, ligand, ligcen_target, reaches_spheroid[nodeno]);

                    int gpn = global_pairs.size();
                    for (l=0; l<3 && l<gpn; l++)
                    {
                        ligand_groups[l] = *(global_pairs[l]->ag);
                        sc_groups[l] = *(global_pairs[l]->scg);

                        if (out_bb_pairs)
                        {
                            n = global_pairs[l]->ag->atct;
                            int j2;
                            for (j2=0; j2<n; j2++)
                                cout << global_pairs[l]->ag->atoms[j2]->name << " ";

                            cout << "- ";

                            n = global_pairs[l]->scg->aminos.size();
                            for (j2=0; j2<n; j2++)
                                cout << global_pairs[l]->scg->aminos[j2]->get_name() << " ";

                            cout << endl;
                        }
                    }

                    if (out_bb_pairs) cout << endl;
                }
                else if (pdpst == pst_constrained)
                {
                    int csiter, ultimate_csidx=0;
                    float lbest_energy;
                    Pose best_cslig;
                    best_cslig.copy_state(ligand);
                    for (csiter=0; csiter<5; csiter++)
                    {
                        Search::do_constrained_search(protein, ligand);
                        float cse = ligand->get_intermol_binding(reinterpret_cast<Molecule**>(reaches_spheroid[nodeno])).summed();
                        if (!csiter || cse > lbest_energy)
                        {
                            lbest_energy = cse;
                            best_cslig.copy_state(ligand);
                            ultimate_csidx = cs_idx;
                        }
                        if (cs_bt[cs_idx] == ionic) break;
                    }
                    best_cslig.restore_state(ligand);
                    cs_idx = ultimate_csidx;

                    #if _dbg_groupsel
                    cout << "Binding constraint:\n" << cs_res[ultimate_csidx]->get_name()
                        << " ~ " << cs_bt[ultimate_csidx]
                        << " ~ " << *cs_lag[ultimate_csidx] << endl << endl << flush;
                    #endif
                }
                else if (pdpst == pst_cavity_fit)
                {
                    float bestc = 0;
                    int bestl = 0;
                    Pose bestp(ligand);
                    bestp.copy_state(ligand);
                    for (l=0; l<ncvtys; l++)
                    {
                        float ctainmt = cvtys[l].find_best_containment(ligand, true);
                        if (!l || ctainmt > bestc)
                        {
                            bestp.copy_state(ligand);
                            bestc = ctainmt;
                            bestl = l;
                        }
                    }
                    bestp.restore_state(ligand);
                    // erase_progressbar(); cout << "Using cavity " << cvtys[bestl].resnos_as_string(protein) << endl << endl;

                    int csidx = Search::choose_cs_pair(protein, ligand);

                    Atom* mtl = (cs_bt[csidx] == mcoord) ? cs_res[csidx]->coordmtl : nullptr;
                    ligand->find_mutual_max_bind_potential(cs_res[csidx]);
                    if (mtl) ligand->stay_close_other = mtl;

                    ligand->movability = MOV_ALL;
                    ligand->enforce_stays();
                }
                else if (pdpst == pst_copyfrom)
                {
                    Search::copy_ligand_position_from_file(protein, ligand, copyfrom_filename.c_str(), copyfrom_ligname, copyfrom_resno);
                }

                // else ligand->recenter(ligcen_target);

                // Best-Binding Algorithm
                // Find a binding pocket feature with a strong potential binding to the ligand.
                std::string alignment_name = "";
                
                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Selected an alignment AA." << endl;
                #endif

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
            for (i=0; i<=SPHREACH_MAX; i++) cfmols[i] = nullptr;
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

            if (n = priority_resnos.size())
            {
                for (j=0; j<n; j++)
                {
                    AminoAcid* aa = protein->get_residue(priority_resnos[j]);
                    if (aa) cfmols[i++] = (Molecule*)aa;
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
            for (; i<=SPHREACH_MAX; i++) cfmols[i] = NULL;

            if (pdpst == pst_cavity_fit && ligand->stay_close_mine && ligand->stay_close_other)
            {
                LocatedVector axis = (SCoord)ligand->stay_close_other->get_location().subtract(ligand->stay_close_mine->get_location());
                axis.origin = ligand->stay_close_mine->get_location();
                float theta, step = fiftyseventh*2;
                Pose bestp(ligand);
                Interaction bestb = Molecule::total_intermol_binding(cfmols);

                for (theta=0; theta<M_PI*2; theta += step)
                {
                    ligand->rotate(axis, step);
                    Interaction linter = Molecule::total_intermol_binding(cfmols);
                    if (linter.attractive > bestb.attractive)
                    {
                        bestb = linter;
                        bestp.copy_state(ligand);
                    }
                }

                bestp.restore_state(ligand);
            }

            ligand->reset_conformer_momenta();
            
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

            #if _dbg_improvements_only_rule
            check_ligand = ligand;
            #endif

            ligand->movability = MOV_ALL;
            if (!flex) for (j=0; j<sphres; j++)
            {
                reaches_spheroid[nodeno][j]->movability = MOV_FLXDESEL;
            }
            ligand->agroups = global_pairs;
            if (output_each_iter) output_iter(0, cfmols);
            if (pdpst == pst_best_binding) ligand->movability = (MovabilityType)(MOV_CAN_AXIAL | MOV_CAN_RECEN | MOV_CAN_FLEX);
            Molecule::conform_molecules(cfmols, iters, &iteration_callback, &GroupPair::align_groups_noconform, progressbar ? &update_progressbar : nullptr);

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
            dr[drcount][nodeno] = DockResult(protein, ligand, size, addl_resno, drcount, waters);
            dr[drcount][nodeno].out_per_res_e = out_per_res_e;
            dr[drcount][nodeno].out_per_btyp_e = out_per_btyp_e;
            dr[drcount][nodeno].out_itemized_e_cutoff = out_itemized_e_cutoff;
            dr[drcount][nodeno].out_lig_int_e = out_lig_int_e;
            dr[drcount][nodeno].out_lig_pol_sat = out_lig_pol_sat;
            dr[drcount][nodeno].out_prox = out_prox;
            dr[drcount][nodeno].out_pro_clash = out_pro_clash;
            dr[drcount][nodeno].out_mc = out_mc;
            dr[drcount][nodeno].out_vdw_repuls = out_vdw_repuls;
            float btot = dr[drcount][nodeno].kJmol;
            float pstot = dr[drcount][nodeno].polsat;
            if (isomers.size()) dr[drcount][nodeno].isomer = ligand->get_name();

            if ((pose==1 && !nodeno) || best_energy > -btot) best_energy = -btot;
            if ((pose==1 && !nodeno) || best_worst_clash > dr[drcount][nodeno].worst_energy) best_worst_clash = dr[drcount][nodeno].worst_energy;

            #if compute_clashdirs
            n = protein->get_end_resno();
            for (i=1; i<=n; i++)
            {
                if (dr[drcount][nodeno].residue_clash[i])
                {
                    for (j=0; regions[j].start; j++)
                    {
                        if (i >= regions[j].start && i <= regions[j].end)
                        {
                            int rgcen, hxno = atoi(&regions[j].name[3]);
                            if (hxno >= 1 && hxno <= 7) rgcen = protein->get_bw50(hxno);
                            else rgcen = (regions[j].start + regions[j].end)/2;
                            k = abs(i - rgcen);
                            if (k < 5) k = 1;
                            else k = (i < rgcen) ? 0 : 2;
                            region_clashes[j][k] = region_clashes[j][k].add(dr[drcount][nodeno].res_clash_dir[i]);
                        }
                    }
                }
            }
            #endif

            dr[drcount][nodeno].proximity = ligand->get_barycenter().get_3d_distance(nodecen);

            if (pdpst == pst_best_binding && out_bb_pairs)
            {
                dr[drcount][nodeno].miscdata += (std::string)"Best-Binding Pairs:\n";
                for (i=0; i<3 && i<global_pairs.size(); i++)
                {
                    n = global_pairs[i]->ag->atct;
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

            if (pdpst == pst_constrained)
            {
                dr[drcount][nodeno].miscdata += (std::string)"Binding constraint:\n";
                dr[drcount][nodeno].miscdata += (std::string)cs_res[cs_idx]->get_name() + (std::string)" ~ ";
                std:stringstream stst;
                stst << cs_bt[cs_idx] << " ~ " << *cs_lag[cs_idx] << endl;
                dr[drcount][nodeno].miscdata += stst.str();
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

            dr[drcount][nodeno].auth = pose;

            if (!nodeno)
            {
                if ((dr[drcount][nodeno].ligand_self + ligand->total_eclipses()) < -clash_limit_per_aa*2)
                {
                    #if _dbg_worst_energy
                    cout << "Internal ligand energy " << -dr[drcount][nodeno].ligand_self << " out of range." << endl << endl;
                    #endif

                    break;          // Exit nodeno loop.
                }
                // else cout << "Internal ligand energy " << -dr[drcount][nodeno].ligand_self << " satisfactory." << endl << endl;

                #if !_dbg_allow_excessive_aa_clashes
                // cout << dr[drcount][nodeno].kJmol << " / " << dr[drcount][nodeno].worst_energy << endl << endl;
                if (dr[drcount][nodeno].worst_energy > l_atom_clash_limit || dr[drcount][nodeno].worst_nrg_aa > clash_limit_per_aa)
                {
                    #if _dbg_worst_energy
                    cout << "Total binding energy " << dr[drcount][nodeno].kJmol
                        << " and worst energy " << dr[drcount][nodeno].worst_energy;
                    if (dr[drcount][nodeno].worst_clash_1 && dr[drcount][nodeno].worst_clash_2)
                    {
                        cout << " (" << dr[drcount][nodeno].worst_clash_1->residue << ":" << dr[drcount][nodeno].worst_clash_1->name
                            << "-" << dr[drcount][nodeno].worst_clash_2->residue << ":" << dr[drcount][nodeno].worst_clash_2->name
                            << ") ";
                    }
                    cout << "; skipping." << endl << endl;
                    #endif

                    // cout << "Least favorable binding energy " << dr[drcount][nodeno].worst_energy << " out of range." << endl << endl;
                    break;          // Exit nodeno loop.
                }
                // else cout << "Least favorable binding energy " << dr[drcount][nodeno].worst_energy << " satisfactory." << endl << endl;
                #endif

                if (pose==1) dr[drcount][nodeno].pose = pose;
                else
                {
                    int bestpose = pose;
                    for (i=0; i<drcount; i++)
                    {
                        if ((	(dr[i][0].kJmol - dr[i][0].ikJmol + dr[i][0].polsat * polar_sat_influence_for_scoring)
                                <
                                (dr[drcount][nodeno].kJmol - dr[drcount][nodeno].ikJmol + dr[drcount][nodeno].polsat * polar_sat_influence_for_scoring)
                            )
                            ||
                            (	(dr[i][0].kJmol + dr[i][0].polsat * polar_sat_influence_for_scoring)
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
            if (btot < kJmol_cutoff)
            {
                #if _dbg_worst_energy
                cout << "Total binding energy " << -btot << "; skipping." << endl << endl;
                #endif
                drcount++;
                break;
            }
            else if (nodeno == pathnodes) drcount++;
        }	// nodeno loop.
    } // pose loop.

    /////////////////////////////////////////////////////////////////////////////////
    // End main loop.
    /////////////////////////////////////////////////////////////////////////////////

    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Finished poses." << endl;
    #endif

    if (progressbar)
    {
        erase_progressbar();
    }

    // Output the dr[][] array in order of increasing pose number.
    cout << endl;
    if (output) *output << endl;

    const float energy_mult = kcal ? _kcal_per_kJ : 1;
    pose = 1;
    std::string auths;
    for (i=1; i<=poses; i++)
    {
        for (j=0; j<poses; j++)
        {
            protein = &pose_proteins[j];
            ligand = &pose_ligands[j+1];

            if (dr[j][0].disqualified) continue;
            dr[j][0].ligpos.restore_state(ligand);

            if (ncvtys)
            {
                dr[j][0].disqualified = true;
                int cno;
                for (cno = 0; cno < ncvtys; cno++)
                {
                    if (!cvtys[cno].count_partials()) continue;
                    CPartial* cp;
                    if (cp = cvtys[cno].point_inside_pocket(ligand->get_barycenter()))
                    {
                        /*cout << -dr[j][0].kJmol << " " << ligand->get_barycenter() << " is inside " << cp->s.center
                            << " of " << cvtys[cno].resnos_as_string(protein) << endl;*/
                        dr[j][0].disqualified = false;
                        break;
                    }
                }
            }
            if (dr[j][0].disqualified) continue;

            if (dr[j][0].pose == i && dr[j][0].pdbdat.length())
            {
                if (dr[j][0].kJmol >= kJmol_cutoff)
                {
                    if (dr[j][0].proximity > size.magnitude()) continue;
                    if (dr[j][0].worst_nrg_aa > clash_limit_per_aa) continue;

                    auths += (std::string)" " + std::to_string(dr[j][0].auth);

                    if (!best_acc_energy) best_acc_energy = -dr[j][0].kJmol;

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
                                if (aa) aa->set_pdb_chain('A');
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

    cout << "Best candidate pose energy: " << (kcal ? best_energy/_kcal_per_kJ : best_energy) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;
    if (output) *output << "Best candidate pose energy: " << (kcal ? best_energy/_kcal_per_kJ : best_energy) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;
    if (debug) *debug << "Best candidate pose energy: " << (kcal ? best_energy/_kcal_per_kJ : best_energy) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;

    cout << "Best accepted pose energy: " << (kcal ? best_acc_energy/_kcal_per_kJ : best_acc_energy) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;
    if (output) *output << "Best accepted pose energy: " << (kcal ? best_acc_energy/_kcal_per_kJ : best_acc_energy) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;
    if (debug) *debug << "Best accepted pose energy: " << (kcal ? best_acc_energy/_kcal_per_kJ : best_acc_energy) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;

    cout << "Best worst clash: " << (kcal ? best_worst_clash/_kcal_per_kJ : best_worst_clash) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;
    if (output) *output << "Best worst clash: " << (kcal ? best_worst_clash/_kcal_per_kJ : best_worst_clash) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;
    if (debug) *debug << "Best worst clash: " << (kcal ? best_worst_clash/_kcal_per_kJ : best_worst_clash) << (kcal ? " kcal/mol." : " kJ/mol.") << endl;

    #if compute_clashdirs
    if (regions)
    {
        cout << endl;
        if (output) *output << endl;
        for (i=0; regions[i].start; i++)
        {
            for (j=0; j<3; j++) if (region_clashes[i][j].r)
            {
                cout << regions[i].name;
                if (output) *output << regions[i].name;
                if (!j)
                {
                    cout << ".nseg";
                    if (output) *output << ".nseg";
                }
                else if (j==1)
                {
                    cout << ".center";
                    if (output) *output << ".center";
                }
                else if (j==2)
                {
                    cout << ".cseg";
                    if (output) *output << ".cseg";
                }
                cout << ".clashdir = " << (Point)region_clashes[i][j] << endl;
                if (output) *output << ".clashdir = " << (Point)region_clashes[i][j] << endl;
            }
        }
        cout << endl;
        if (output) *output << endl;
    }
    #endif

    if (met) delete met;

    time_t finished = time(NULL);
    int seconds = finished-began;
    int minutes = seconds/60;
    seconds -= 60*minutes;
    int hours = minutes/60;
    minutes -= 60*hours;

    std::string elapsed;
    if (hours)
    {
        elapsed += std::to_string(hours);
        elapsed += (std::string)":";
        if (minutes < 10) elapsed += (std::string)"0";
    }
    elapsed += std::to_string(minutes);
    elapsed += (std::string)":";
    if (seconds < 10) elapsed += (std::string)"0";
    elapsed += std::to_string(seconds);

    cout << "\nCalculation time: " << elapsed << "." << endl;
    if (output) *output << "\nCalculation time: " << elapsed << "." << endl;
    if (debug) *debug << "\nCalculation time: " << elapsed << "." << endl;

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
            apply_protein_specific_settings(protein);
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


















