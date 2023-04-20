#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <algorithm>
#include "classes/protein.h"

using namespace std;

struct ResiduePlaceholder
{
    int node = 0;
    int resno = 0;
    std::string bw;

    void set(const char* str)
    {
        if (strchr(str, '.')) bw = str;
        else resno = atoi(str);
    }

    void resolve_resno(Protein* prot)
    {
        int hxno = atoi(bw.c_str());
        const char* dot = strchr(bw.c_str(), '.');
        int bwpos = atoi(dot+1);
        resno = prot->get_bw50(hxno) + bwpos - 50;
    }
};

struct DockResult
{
    int pose;
    float kJmol;
    float ikJmol;
    char** metric;
    float* mkJmol;
    float* imkJmol;
    float* mvdWrepl;
    float* imvdWrepl;
    std::string softrock;
    std::string pdbdat;
    float bytype[_INTER_TYPES_LIMIT];
    float ibytype[_INTER_TYPES_LIMIT];
    float proximity;                    // How far the ligand center is from the node's center.
    float tripswitch;                   // Effect of the ligand on the receptor's trip switch.
    float polsat;
    float protclash;
};

enum HxRotAxisType
{
    hxrat_Cartesian,
    hxrat_region,
    hxrat_BW,
    hxrat_atom,
    hxrat_helical
};

struct AcvHxRot
{
    std::string regname;
    int nodeno;
    int start_resno;
    int end_resno;
    Point transform;
    int origin_resno;
    SCoord axis;
    float theta;
    bool soft;
    float dtheta;
    HxRotAxisType axis_type = hxrat_Cartesian;
    std::string start_resno_str;
    std::string end_resno_str;
    std::string origin_resno_str;
    std::string axis_str;
};

struct AcvBndRot
{
    int resno;
    std::string aname;
    std::string bname;
    Atom* atom;
    Atom* btom;
    Bond* bond;
    float theta;
};

struct SoftBias
{
    std::string region_name;
    float radial_transform = 0;                 // Motion away from or towards the pocket center.
    float angular_transform = 0;                // Motion towards and away from neighboring helices.
    float vertical_transform = 0;               // Motion in the extracelllar or cytoplasmic direction.
    float helical_rotation = 0;                 // Rotation about the helical axis.
    float radial_rotation = 0;                  // Rotation about the imaginary line to the pocket center.
    float transverse_rotation = 0;              // Rotation about the imaginary line perpendicular to the pocket center.
};

std::vector<int> extra_wt;

#if _use_gloms
struct AtomGlom
{
    std::vector<Atom*> atoms;

    Point get_center()
    {
        int i;
        float mass = 0;
        Point result(0,0,0);
        int atct = atoms.size();
        if (!atct) return Point(0,0,0);
        for (i=0; i<atct; i++)
        {
            Point pt = atoms[i]->get_location();
            float m = atoms[i]->get_atomic_weight();
            pt.scale(pt.magnitude() * m);
            result = result.add(pt);
            mass += m;
        }
        if (mass) result.scale(result.magnitude() / mass);
        return result;
    }

    float get_pi()
    {
        int atct = atoms.size();
        if (!atct) return 0;
        int i;
        float result = 0;
        for (i=0; i<atct; i++)
        {
            if (atoms[i]->is_pi()) result += 1;
        }
        return result;
    }

    float get_polarity()
    {
        int atct = atoms.size();
        if (!atct) return 0;
        int i;
        float result = 0;
        for (i=0; i<atct; i++)
        {
            result += fabs(atoms[i]->is_polar());
        }
        return result;
    }

    float get_ionic()
    {
        int atct = atoms.size();
        if (!atct) return 0;
        int i;
        float result = 0;
        for (i=0; i<atct; i++)
        {
            float c = atoms[i]->get_charge();
            if (c) result += c;
            else
            {
                if (atoms[i]->get_family() == PNICTOGEN && !atoms[i]->is_amide())
                {
                    result += 0.5;
                }
            }
        }
        return result;
    }

    float get_sum()
    {
        return get_ionic()*60 + get_polarity()*25 + get_pi()*2;
    }

    float distance_to(Point pt)
    {
        int atct = atoms.size();
        if (!atct) return -1;
        int i;
        float result = 0;
        for (i=0; i<atct; i++)
        {
            float f = atoms[i]->get_location().get_3d_distance(pt);
            if (!i || !f || f < result) result = f;
        }
        return result;
    }

    float bounds()
    {
        int atct = atoms.size();
        if (!atct) return -1;
        int i;
        Point ptmin(0,0,0), ptmax(0,0,0);
        for (i=0; i<atct; i++)
        {
            Point aloc = atoms[i]->get_location();
            if (!i)
            {
                ptmin = aloc;
                ptmax = aloc;
            }
            else
            {
                if (aloc.x < ptmin.x) ptmin.x = aloc.x;
                if (aloc.y < ptmin.y) ptmin.y = aloc.y;
                if (aloc.z < ptmin.z) ptmin.z = aloc.z;
                if (aloc.x > ptmax.x) ptmax.x = aloc.x;
                if (aloc.y > ptmax.y) ptmax.y = aloc.y;
                if (aloc.z > ptmax.z) ptmax.z = aloc.z;
            }
        }

        return ptmax.get_3d_distance(ptmin);
    }

    float compatibility(AminoAcid* aa)
    {
        int atct = atoms.size();
        if (!atct) return 0;

        float result = 0;
        float lgi = get_ionic(), lgh = get_polarity(), lgp = get_pi();

        float aachg = aa->get_charge();
        if (aa->conditionally_basic()) aachg += 0.5;
        if (lgi && aachg && sgn(lgi) != -sgn(aachg)) return 0;

        if (aa->hydrophilicity() > 0.25)
        {
            if ((lgh / atct) < 0.19) return 0;
        }
        else
        {
            if ((lgh / atct) > 0.33333) return 0;
        }

        if (lgh)
        {
            int i;
            int atct = atoms.size();
            if (!aa->has_hbond_acceptors())
            {
                lgh = 0;
                if (atct)
                for (i=0; i<atct; i++)
                {
                    if (atoms[i]->is_polar() < 0) lgh++;
                }
            }
            else if (!aa->has_hbond_donors())
            {
                lgh = 0;
                if (atct)
                for (i=0; i<atct; i++)
                {
                    if (atoms[i]->is_polar() > 0) lgh++;
                }
            }
        }
 
        if (lgi) result += lgi * -aa->bindability_by_type(ionic) * 1000;
        if (lgh) result += fabs(lgh) * fabs(aa->bindability_by_type(hbond)) * 30;
        if (lgp) result += aa->bindability_by_type(pi);

        return result;
    }
};

struct ResidueGlom
{
    std::vector<AminoAcid*> aminos;

    Point get_center()
    {
        int amsz = aminos.size();
        if (!amsz) return Point(0,0,0);
        int i, j;
        j = 0;
        Point result(0,0,0);
        for (i=0; i<amsz; i++)
        {
            Atom** aa = aminos[i]->get_most_bindable(1);
            Atom* a = aa[0]; // = aminos[i]->get_atom("CB");
            delete aa;
            // if (!a) a = aminos[i]->get_atom("CA");      // even though glycine probably shouldn't be part of a glom.
            if (a)
            {
                Point pt = a->get_location();
                result = result.add(pt);
                j++;
            }
        }
        if (j) result.scale(result.magnitude() / j);
        return result;
    }

    float distance_to(Point pt)
    {
        return pt.get_3d_distance(get_center());
    }

    float compatibility(AtomGlom* ag)
    {
        int amsz = aminos.size();
        if (!amsz) return 0;
        float result = 0;
        int i;
        bool has_acids = false, has_his = false;
        for (i=0; i<amsz; i++)
        {
            if (aminos[i]->get_charge() < 0) has_acids = true;
            if (aminos[i]->conditionally_basic()) has_his = true;

            float f = ag->compatibility(aminos[i]);

            if (extra_wt.size()
                    &&
                    std::find(extra_wt.begin(), extra_wt.end(), aminos[i]->get_residue_no())!=extra_wt.end()
            )
            {
                f *= 1.25;		// Extra weight for residues mentioned in a CEN RES or PATH RES parameter.
            }

            result += f;
        }
        if (has_acids && has_his && ag->get_ionic() > 0) result -= ag->get_ionic()*30;

        return result;
    }
};

#endif

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
int movie_offset = 0;
char configfname[256];
char protfname[256];
char protafname[256];
char tplfname[256];
char ligfname[256];
char outfname[256];
Point pocketcen, loneliest, pocketsize, ligbbox;
std::ofstream *output = NULL;

std::vector<int> exclusion;

std::string CEN_buf = "";
std::vector<std::string> pathstrs;
std::vector<std::string> states;

bool configset=false, protset=false, tplset=false, ligset=false, pktset=false;

Protein* protein;
Protein* ptemplt;
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
std::string temp_pdb_file;
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

#if _use_gloms
AtomGlom ligand_gloms[3];
ResidueGlom sc_gloms[3];
#else
Atom** ligbb = nullptr;
Atom** ligbbh = nullptr;

intera_type lig_inter_typ[5];
Molecule* alignment_aa[5];
#endif

Pose pullaway_undo;
float last_ttl_bb_dist = 0;

float* initial_binding;
float* initial_vdWrepl;
float init_total_binding_by_type[_INTER_TYPES_LIMIT];

Point active_matrix_n[16], active_matrix_c[16], active_matrix_m[16];
int active_matrix_count = 0, active_matrix_node = -1, active_matrix_type = 0, deactivate_node = -1;
std::vector<AcvHxRot> active_helix_rots;
std::vector<AcvHxRot> orig_active_helix_rots;
std::vector<AcvBndRot> active_bond_rots;
std::vector<int> tripswitch_clashables;
std::vector<ResiduePlaceholder> required_contacts;

bool soft_pocket = false;
std::string soft_names;
std::vector<Region> soft_rgns;
std::vector<SoftBias> soft_biases;
std::vector<int>flexible_resnos;
std::vector<ResiduePlaceholder>forced_flexible_resnos;
std::vector<ResiduePlaceholder>forced_static_resnos;

float *g_rgnxform_r = nullptr, *g_rgnxform_theta = nullptr, *g_rgnxform_y = nullptr;
float *g_rgnrot_alpha = nullptr, *g_rgnrot_w = nullptr, *g_rgnrot_u = nullptr;

#if _dummy_atoms_for_debug
std::vector<Atom> dummies;
#endif

void colorize(float f)
{
    float red, green, blue;

    if (f >= 0)
    {
        f = sqrt(f/5);
        blue = 128 + 128 * f;
        green = fmax(48, (f-1) * 255);
        red = fmax(64, (f-2) * 255);
    }
    else
    {
        f = sqrt(-f)*3;
        f = fmax(0,fmin(128,f*16));
        red = 128+f;
        blue = 128-f;
        green = 0.333 * red + 0.666 * blue;
    }

    int r, g, b;

    r = max(0, min(255, (int)red));
    g = max(0, min(255, (int)green));
    b = max(0, min(255, (int)blue));

    cout << "\x1b[38;2;" << r << ";" << g << ";" << b << "m";
}

void colorless()
{
    cout << "\x1b[0m";
}

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

SoftBias* get_soft_bias_from_region(const char* region)
{
    int sz = soft_biases.size();
    if (!sz) return nullptr;
    int i;
    for (i=0; i<sz; i++)
    {
        if (soft_biases[i].region_name == (std::string)region) return &soft_biases[i];
    }
    return nullptr;
}

void iteration_callback(int iter)
{
    // if (kJmol_cutoff > 0 && ligand->lastbind >= kJmol_cutoff) iter = (iters-1);
    int l;
    float prebind;

    if (soft_pocket && iter >= 10 && g_rgnrot_alpha && g_rgnrot_u && g_rgnrot_w && g_rgnxform_r && g_rgnxform_theta && g_rgnxform_y)
    {
        int sz = soft_rgns.size();
        if (sz)
        {
            for (l=0; l<sz; l++)
            {
                SoftBias* sb = get_soft_bias_from_region(soft_rgns[l].name.c_str());
                if (!l) prebind = protein->get_intermol_binding(ligand)*soft_ligand_importance + protein->get_internal_binding()*_kJmol_cuA;         // /'kʒmɑɫ.kju.ə/
                
                #if _dbg_soft
                cout << iter << ": from " << prebind;
                #endif

                int tweak = rand() % 6;
                float amount = nanf("unbiased");

                if (sb)
                {
                    switch (tweak)
                    {
                        case 0:
                        amount = sb->radial_transform;
                        break;

                        case 1:
                        amount = sb->angular_transform;
                        break;

                        case 2:
                        amount = sb->vertical_transform;
                        break;

                        case 3:
                        amount = sb->helical_rotation;
                        break;

                        case 4:
                        amount = sb->radial_rotation;
                        break;

                        case 5:
                        amount = sb->transverse_rotation;
                        break;

                        default:
                        ;
                    }
                }

                if (isnan(amount) || !amount) amount = frand(-1, 1);
                else
                {
                    if (amount > 0) amount = frand(-amount*soft_bias_overlap, amount);
                    else amount = frand(amount, -amount*soft_bias_overlap);
                }

                Point ptrgn = protein->get_region_center(soft_rgns[l].start, soft_rgns[l].end);
                SCoord r = ptrgn.subtract(loneliest);
                r.theta = 0;
                r.r = amount;
                Point pr1 = loneliest.add(r);
                Point pr2 = pr1;
                pr2.y += 20;
                SCoord normal = compute_normal(loneliest, pr1, pr2);
                normal.r = amount;
                SCoord alpha = protein->get_region_axis(soft_rgns[l].start, soft_rgns[l].end);
                alpha.r = amount;
                Point rgncen = protein->get_region_center(soft_rgns[l].start, soft_rgns[l].end);

                switch (tweak)
                {
                    case 0:
                    protein->move_piece(soft_rgns[l].start, soft_rgns[l].end, r);
                    break;

                    case 1:
                    protein->move_piece(soft_rgns[l].start, soft_rgns[l].end, normal);
                    break;

                    case 2:
                    protein->move_piece(soft_rgns[l].start, soft_rgns[l].end, alpha);
                    break;

                    case 3:
                    protein->rotate_piece(soft_rgns[l].start, soft_rgns[l].end, rgncen, alpha, amount/10);
                    break;

                    case 4:
                    protein->rotate_piece(soft_rgns[l].start, soft_rgns[l].end, rgncen, r, amount/50);
                    break;

                    case 5:
                    protein->rotate_piece(soft_rgns[l].start, soft_rgns[l].end, rgncen, normal, amount/30);
                    break;

                    default:
                    ;
                }

                float postbind = protein->get_intermol_binding(ligand)*soft_ligand_importance + protein->get_internal_binding()*_kJmol_cuA;
                #if _dbg_soft
                cout << " to " << postbind;
                #endif

                switch (tweak)
                {
                    case 0:
                    if (postbind > prebind) g_rgnxform_r[l] += amount;
                    else
                    {
                        r.r = -r.r;
                        protein->move_piece(soft_rgns[l].start, soft_rgns[l].end, r);
                        #if _dbg_soft
                        cout << " reverting r.";
                        #endif
                    }
                    break;

                    case 1:
                    if (postbind > prebind) g_rgnxform_theta[l] += amount;
                    else
                    {
                        normal.r = -normal.r;
                        protein->move_piece(soft_rgns[l].start, soft_rgns[l].end, normal);
                        #if _dbg_soft
                        cout << " reverting normal.";
                        #endif
                    }
                    break;

                    case 2:
                    if (postbind > prebind) g_rgnxform_y[l] += amount;
                    else
                    {
                        alpha.r = -alpha.r;
                        protein->move_piece(soft_rgns[l].start, soft_rgns[l].end, alpha);
                        #if _dbg_soft
                        cout << " reverting alpha.";
                        #endif
                    }
                    break;

                    case 3:
                    if (postbind > prebind) g_rgnrot_alpha[l] += amount/10;
                    else
                    {
                        protein->rotate_piece(soft_rgns[l].start, soft_rgns[l].end, rgncen, alpha, -amount/10);
                        #if _dbg_soft
                        cout << " reverting alpha rot.";
                        #endif
                    }
                    break;

                    case 4:
                    if (postbind > prebind) g_rgnrot_w[l] += amount/50;
                    else
                    {
                        protein->rotate_piece(soft_rgns[l].start, soft_rgns[l].end, rgncen, r, -amount/50);
                        #if _dbg_soft
                        cout << " reverting r rot.";
                        #endif
                    }
                    break;

                    case 5:
                    if (postbind > prebind) g_rgnrot_u[l] += amount/30;
                    else
                    {
                        protein->rotate_piece(soft_rgns[l].start, soft_rgns[l].end, rgncen, normal, -amount/30);
                        #if _dbg_soft
                        cout << " reverting normal rot.";
                        #endif
                    }
                    break;

                    default:
                    ;
                }

                if (postbind > prebind) prebind = postbind;

                #if _dbg_soft
                cout << endl;
                #endif
            }
        }
    }


    Point bary = ligand->get_barycenter();

    int ac = ligand->get_atom_count();
    float bbest = 0;
    Atom *atom, *btom;
    int i, j;

    float progress = (float)iter / iters;
    // float lsrca = (1.0 - progress) * soft_rock_clash_allowance;
    float lsrca = (progress < 0.5) ? soft_rock_clash_allowance : 0;

    l = active_helix_rots.size();
    if (l && nodeno >= active_matrix_node)
    {
        for (i=0; i<l; i++)
        {
            if (active_helix_rots[i].soft)
            {
                /*AminoAcid* resphres[SPHREACH_MAX+4];
                for (i=0; i<SPHREACH_MAX+4; i++) resphres[i] = nullptr;
                int sphres = protein->get_residues_can_clash_ligand(resphres, ligand, bary, size, addl_resno);*/

                // TODO: Check that any residues of reaches_spheroid[nodeno] are in range for rocking motion
                // and output a warning if not.

                float before = ligand->get_intermol_binding(reinterpret_cast<Molecule**>(reaches_spheroid[nodeno]));
                float pic = protein->get_internal_clashes(active_helix_rots[i].start_resno, active_helix_rots[i].end_resno, repack_on_hxr);
                before -= pic * _kJmol_cuA * soft_rock_clash_penalty;
                
                for (j=0; j<sphres; j++)
                    if (reaches_spheroid[nodeno][j])
                    {
                        float lclash = reaches_spheroid[nodeno][j]->get_intermol_clashes(
                            reinterpret_cast<Molecule**>(reaches_spheroid[nodeno])) * _kJmol_cuA;
                        if (lclash > soft_rock_clash_allowance) before -= (lclash - soft_rock_clash_allowance);
                    }
                int sr = active_helix_rots[i].start_resno;
                int er = active_helix_rots[i].end_resno;
                int mr = active_helix_rots[i].origin_resno;

                #if _dbg_softrock
                Point wasat = protein->get_atom_location(sr, "CA");
                #endif

                protein->rotate_piece(sr, er, protein->get_atom_location(mr, "CA"),
                    active_helix_rots[i].axis, active_helix_rots[i].dtheta);

                #if _dbg_softrock
                cout << "Rotating residues " << sr << "-" << er << " by "
                        << (active_helix_rots[i].dtheta * fiftyseven) << "deg." << endl;
                cout << sr << ":CA moved from " << wasat << " to " << protein->get_atom_location(sr, "CA") << endl;
                #endif

                float after = ligand->get_intermol_binding(reinterpret_cast<Molecule**>(reaches_spheroid[nodeno]));
                float pic1 = protein->get_internal_clashes(active_helix_rots[i].start_resno, active_helix_rots[i].end_resno, repack_on_hxr);
                pic1 -= soft_rock_clash_allowance;
                after -= pic1 * _kJmol_cuA * soft_rock_clash_penalty;
                #if _dbg_rock_pic
                cout << "Iteration " << iter << " angle " << ((active_helix_rots[i].theta + active_helix_rots[i].dtheta) * fiftyseven)
                     << " internal clashes was " << pic << " now " << pic1 << endl;
                #endif
                
                for (j=0; j<sphres; j++)
                    if (reaches_spheroid[nodeno][j])
                    {
                        float lclash = reaches_spheroid[nodeno][j]->get_intermol_clashes(
                            reinterpret_cast<Molecule**>(reaches_spheroid[nodeno])) * _kJmol_cuA;
                        if (lclash > soft_rock_clash_allowance) after -= (lclash - soft_rock_clash_allowance);
                    }

                if (after > before)
                {
                    active_helix_rots[i].theta += active_helix_rots[i].dtheta;
                    active_helix_rots[i].dtheta *= 1.1;
                    #if _dbg_softrock
                    cout << "Ligand binding improved from " << before << " to " << after << ". "
                         << active_helix_rots[i].regname << " rock updated to " << (active_helix_rots[i].theta*fiftyseven)
                         << "deg." << endl;
                    #endif
                }
                else
                {
                    protein->rotate_piece(sr, er, protein->get_atom_location(mr, "CA"),
                        active_helix_rots[i].axis, -active_helix_rots[i].dtheta);
                    active_helix_rots[i].dtheta *= -0.8;
                    #if _dbg_softrock
                    cout << "Ligand binding was " << before << " now " << after << ". "
                         << active_helix_rots[i].regname << " rock maintained at " << (active_helix_rots[i].theta*fiftyseven)
                         << "deg." << endl;
                    #endif
                }
            }
        }
    }

    if (!iter) goto _oei;
    if (iter == (iters-1)) goto _oei;
    
    #if enforce_no_bb_pullaway
    #if _use_gloms
    if (ligand_gloms[0].atoms.size())
    #else
    if (ligbb)
    #endif
    {
        float ttl_bb_dist = 0;
        for (l=0; l<3; l++)
        {
            #if _use_gloms
            if (ligand_gloms[l].atoms.size() && sc_gloms[l].aminos.size())
            #else
            if (ligbb[l] && alignment_aa[l])
            #endif
            {
                #if _use_gloms
                float r = ligand_gloms[l].distance_to(sc_gloms[l].get_center());
                if (r < 2.5) r = 2.5;
                if (r > _INTERA_R_CUTOFF) r = _INTERA_R_CUTOFF;
                ttl_bb_dist += r;
                #if _dbg_bb_pullaway
                cout << "Ligand atoms ";
                for (i=0; i<ligand_gloms[l].atoms.size(); i++) cout << ligand_gloms[l].atoms[i]->name << " ";
                cout << "are " << r << " A from residues";
                for (i=0; i<sc_gloms[l].aminos.size(); i++) cout << " " << sc_gloms[l].aminos[i]->get_3letter() << sc_gloms[l].aminos[i]->get_residue_no();
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

        if (!iter || ttl_bb_dist < (1.0+bb_pullaway_allowance)*last_ttl_bb_dist)
        {
            pullaway_undo.copy_state(ligand);
            last_ttl_bb_dist = ttl_bb_dist;
        }
        else
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
            atom = ligand->get_atom(i);
            bbest = atom->strongest_bind_energy;
            btom = atom->strongest_bind_atom;
        }
    }

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
        Star discrete[SPHREACH_MAX+4];
        /*discrete[0].pmol = gcfmols[0];
        discrete[1].pmol = gcfmols[1];*/

        int offset = 2;             // For some strange reason, if this is set to 1 the TAAR8 test fails.

        if (waters) offset += maxh2o;

        for (i=0; i<offset; i++) discrete[i].pmol = gcfmols[i];

        /*AminoAcid* resphres[SPHREACH_MAX+4];
        for (i=0; i<SPHREACH_MAX+4; i++) resphres[i] = nullptr;
        int sphres = protein->get_residues_can_clash_ligand(resphres, ligand, bary, size, addl_resno);*/
        //cout << "Sphres: " << sphres << endl;
        for (i=0; i<sphres; i++)
        {
            discrete[i+offset].paa = reaches_spheroid[nodeno][i];
        }
        discrete[sphres+offset].n = 0;

        sphres += offset;
        for (i=0; i<sphres; i++) gcfmols[i] = discrete[i].pmol;
        gcfmols[sphres] = nullptr;
    }

    _oei:
    #if _dbg_glomsel
    cout << "." << flush;
    if (iter == iters-1) cout << endl << endl;
    #endif

    if (output_each_iter)
    {
        std::string itersfname = (std::string)"tmp/" + (std::string)protein->get_name() + (std::string)"_iters.dock";
        int liter = iter + movie_offset;
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
    else if (!strcmp(words[0], "ACVHXR"))
    {
        AcvHxRot ahr;
        int n = 1;
        ahr.regname              = words[n++];
        ahr.nodeno               = -1;
        ahr.start_resno_str      = words[n];
        ahr.start_resno          = atoi(words[n++]);
        ahr.end_resno_str        = words[n];
        ahr.end_resno            = atoi(words[n++]);
        ahr.origin_resno_str     = words[n];
        ahr.origin_resno         = atoi(words[n++]);
        ahr.transform            = Point( atof(words[n]), atof(words[n+1]), atof(words[n+2]) ); n += 3;
        ahr.axis                 = Point( atof(words[n]), atof(words[n+1]), atof(words[n+2]) ); n += 3;
        ahr.theta                = atof(words[n]) * fiftyseventh;
        ahr.soft                 = strchr(words[n++], '?') ? true : false;
        if (ahr.soft) ahr.dtheta = 0.01;
        active_helix_rots.push_back(ahr);
        orig_active_helix_rots.push_back(ahr);
        optsecho = (std::string)"Active helix rotation " + to_string(ahr.start_resno) + (std::string)"-" + to_string(ahr.end_resno);
    }
    else if (!strcmp(words[0], "ACVMX"))
    {
        if (!words[1] || !words[1][3])
        {
            cout << "Missing region identifier for active matrix." << endl << flush;
            throw 0xbad512e;
        }
        if (    words[1][0] != 'T'
            ||  words[1][1] != 'M'
            ||  words[1][2] != 'R'
            )
        {
            cout << "Unknown region identifier for active matrix: " << words[1] << endl << flush;
            throw 0xbad512e;
        }

        if (!active_matrix_count)
        {
            for (i=0; i<16; i++) active_matrix_c[i] = active_matrix_n[i] = active_matrix_m[i] = Point(0,0,0);
        }

        int n = atoi(&words[1][3]);
        if (n > active_matrix_count) active_matrix_count = n;

        active_matrix_n[n].x = atof(words[2]);
        active_matrix_n[n].y = atof(words[3]);
        active_matrix_n[n].z = atof(words[4]);
        if (words[8] && words[9] && words[10])
        {
            active_matrix_type = 9;
            active_matrix_m[n].x = atof(words[5]);
            active_matrix_m[n].y = atof(words[6]);
            active_matrix_m[n].z = atof(words[7]);
            active_matrix_c[n].x = atof(words[8]);
            active_matrix_c[n].y = atof(words[9]);
            active_matrix_c[n].z = atof(words[10]);
        }
        else
        {
            active_matrix_type = 6;
            active_matrix_c[n].x = atof(words[5]);
            active_matrix_c[n].y = atof(words[6]);
            active_matrix_c[n].z = atof(words[7]);
        }
        optsecho = (std::string)"Active matrix for TMR" + to_string(n);
    }
    else if (!strcmp(words[0], "ACVNODE"))
    {
        active_matrix_node = atoi(words[1]);
        optsecho = (std::string)"Active node is " + to_string(active_matrix_node);
    }
    else if (!strcmp(words[0], "APPENDPROT") || !strcmp(words[0], "OPEND"))
    {
        append_pdb = true;
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
    else if (!strcmp(words[0], "HXR"))
    {
        AcvHxRot ahr;
        int n = 1;
        ahr.nodeno               = atoi(words[n++]);
        ahr.regname              = words[n++];
        ahr.start_resno_str      = words[n];
        ahr.start_resno          = atoi(words[n++]);
        ahr.end_resno_str        = words[n];
        ahr.end_resno            = atoi(words[n++]);
        ahr.origin_resno_str     = words[n];
        ahr.origin_resno         = atoi(words[n++]);
        ahr.transform            = Point( atof(words[n]), atof(words[n+1]), atof(words[n+2]) ); n += 3;
        if (!strcmp(words[n], "region"))
        {
            ahr.axis_type = hxrat_region;
            ahr.axis_str = words[n+1];
            n += 2;
        }
        else if (!strcmp(words[n], "bw"))
        {
            ahr.axis_type = hxrat_BW;
            ahr.axis_str = words[n+1];
            n += 2;
        }
        else if (!strcmp(words[n], "atom"))
        {
            ahr.axis_type = hxrat_atom;
            ahr.axis_str = words[n+1];
            n += 2;
        }
        else if (!strcmp(words[n], "helical"))
        {
            ahr.axis_type = hxrat_helical;
            n += 1;
        }
        else
        {
            ahr.axis             = Point( atof(words[n]), atof(words[n+1]), atof(words[n+2]) );
            n += 3;
        }
        ahr.theta                = atof(words[n]) * fiftyseventh;
        ahr.soft                 = strchr(words[n++], '?') ? true : false;
        if (ahr.soft) ahr.dtheta = 0.01;
        active_helix_rots.push_back(ahr);
        orig_active_helix_rots.push_back(ahr);
        optsecho = (std::string)"Helix rotation " + to_string(ahr.start_resno) + (std::string)"-" + to_string(ahr.end_resno);
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
        return 1;
    }
    else if (!strcmp(words[0], "MCOORD"))
    {
        int j=0;
        optsecho = "Metal coordination on residues ";
        for (i=1; words[i]; i++)
        {
            if (words[i][0] == '-' && words[i][1] == '-') break;
            mcoord_resno[j++] = atoi(words[i]);
            optsecho += to_string(atoi(words[i])) + (std::string)" ";
        }
        mcoord_resno[j] = 0;
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
    else if (!strcmp(words[0], "PROT"))
    {
        strcpy(protfname, words[1]);
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
        soft_pocket = true;
        optsecho = "Soft regions ";
        soft_names = "";
        for (i=1; words[i]; i++)
        {
            if (words[i][0] == '-' && words[i][1] == '-') break;
            soft_names += (std::string)words[i] + (std::string)" ";
        }
        optsecho += soft_names + (std::string)" ";
        return i-1;
    }
    else if (!strcmp(words[0], "SOFTBIAS"))
    {
        SoftBias lbias;
        i=1;
        if (!words[i]) throw 0xbad50f7e;
        lbias.region_name = words[i++];

        if (!words[i]) throw 0xbad50f7e;
        lbias.radial_transform = atof(words[i++]);
        if (!words[i]) throw 0xbad50f7e;
        lbias.angular_transform = atof(words[i++]);
        if (!words[i]) throw 0xbad50f7e;
        lbias.vertical_transform = atof(words[i++]);

        if (!words[i]) throw 0xbad50f7e;
        lbias.helical_rotation = atof(words[i++]);
        if (!words[i]) throw 0xbad50f7e;
        lbias.radial_rotation = atof(words[i++]);
        if (!words[i]) throw 0xbad50f7e;
        lbias.transverse_rotation = atof(words[i++]);

        soft_biases.push_back(lbias);
        return i-1;
    }
    else if (!strcmp(words[0], "STATE"))
    {
        states.push_back(origbuff);
        optsecho = "Added state " + (std::string)origbuff;
        return 0;
    }
    else if (!strcmp(words[0], "TEMPLATE"))
    {
        strcpy(tplfname, words[1]);
        tplset = true;
        // optsecho = "Protein file is " + (std::string)protfname;
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

            delete[] words;
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
            Molecule::multimol_conform(
                prem /*reinterpret_cast<Molecule**>(preaa)*/,
                delete_me = protein->all_residues_as_molecules(),
                iters*pre_ligand_iteration_ratio
            );
            delete[] delete_me;
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
    float loneliness, blone=0, xrad, yrad, zrad, lrad, step, bestxr, bestyr, bestzr, score, worth, weight, bestscore;
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

                                        if (extra_wt.size()
                                                &&
                                                std::find(extra_wt.begin(), extra_wt.end(), tsphres[j]->get_residue_no())!=extra_wt.end()
                                        )
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

                        cout << "Tumble score " << score << " for ligand box " << m.get_bounding_box() << endl;


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
                    clash = tsphres[i]->get_intermol_clashes(&m);

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

    for (i=1; i<argc; i++)
    {
        if (argv[i][0] == '-' && argv[i][1] == '-')
        {
            argv[i] += 2;
            for (j=0; argv[i][j]; j++) if (argv[i][j] >= 'a' && argv[i][j] <= 'z') argv[i][j] &= 0x5f;
            j = interpret_config_line(&argv[i]);
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

    // Load the protein or return an error.
    /* Protein p(protfname);
    protein = &p; */

    char protid[255];
    char* slash = strrchr(protfname, '/');
    if (!slash) slash = strrchr(protfname, '\\');
    strcpy(protid, slash ? slash+1 : protfname );
    char* dot = strchr(protid, '.');
    if (dot) *dot = 0;

    protein = new Protein(protid);
    pf = fopen(protfname, "r");
    if (!pf)
    {
        cout << "Error trying to read " << protfname << endl;
        return 0xbadf12e;
    }
    protein->load_pdb(pf);
    fclose(pf);
    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Loaded protein." << endl;
    #endif

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

        protein->homology_conform(ptemplt);

        temp_pdb_file = "tmp/homolog.pdb";

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

        temp_pdb_file = "tmp/hydro.pdb";

        pf = fopen(temp_pdb_file.c_str(), "wb");
        protein->save_pdb(pf);
        fclose(pf);
    }

    prepare_acv_bond_rots();

    for (i=0; i<required_contacts.size(); i++)
    {
        required_contacts[i].resolve_resno(protein);
    }

    int l;

    if (soft_pocket)
    {
        if (soft_names.size())
        {
            char lbuff[soft_names.size()+4];
            strcpy(lbuff, soft_names.c_str());
            char** words = chop_spaced_words(lbuff);
            const Region* prgn = protein->get_regions();
            for (i=0; i<PROT_MAX_RGN; i++)
            {
                if (!prgn[i].start || !prgn[i].end) break;
                for (l=0; words[l]; l++)
                {
                    if (!strcmp(words[l], prgn[i].name.c_str())) soft_rgns.push_back(prgn[i]);
                }
            }
        }
        else
        {
            const Region* prgn = protein->get_regions();
            for (i=0; i<PROT_MAX_RGN; i++)
            {
                if (!prgn[i].start || !prgn[i].end) break;
                soft_rgns.push_back(prgn[i]);
            }
        }

        if (!soft_rgns.size())
        {
            cout << "Error: no regions in PDB or specified regions not found." << endl;
            throw 0xbad5e697;
        }
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

    #if !_use_gloms
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

    if (!strcmp(words[1], "RES"))
    {
        for (i=2; words[i]; i++)
        {
            extra_wt.push_back(interpret_resno(words[i]));
        }
    }
    pktset = true;

    protein->mcoord_resnos = mcoord_resno;

    for (i=0; mcoord_resno[i]; i++) addl_resno[i] = mcoord_resno[i];
    for (l=0; l < tripswitch_clashables.size(); l++) addl_resno[i+l] = tripswitch_clashables[l];
    addl_resno[i+l] = 0;

    // Load the ligand or return an error.
    Molecule m(ligfname);
    ligand = &m;
    char* ext = get_file_ext(ligfname);
    if (!ext)
    {
        cout << "Ligand file is missing its extension! " << ligfname << endl;
        return 0xbadf12e;
    }

    for (i=0; i<65536; i++) buffer[i] = 0;

    size_t wgaf;
    switch (ext[0])
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
        m.from_sdf(buffer);
        break;

    case 'p':
    case 'P':
        pf = fopen(ligfname, "r");
        if (!pf)
        {
            cout << "Error trying to read " << ligfname << endl;
            return 0xbadf12e;
        }
        m.from_pdb(pf);
        fclose(pf);
        break;

    default:
        cout << "Unrecognized ligand file extension: " << ext << endl;
        return 0xbadf12e;
    }

    m.minimize_internal_clashes();

    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Loaded ligand." << endl;
    #endif

    Point box = m.get_bounding_box();

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
        Molecule::multimol_conform(lig_grp, prealign_res, prealign_iters, nullptr);

        // Then line up residues to ligand.
        if (flex) Molecule::multimol_conform(prealign_res, lig_grp, prealign_iters, nullptr);
        ligand->movability = MOV_ALL;

        delete[] words;
        delete prealign_res;
        delete lig_grp;
    }

    #endif

    // Best-Binding Code
    if (use_bestbind_algorithm)
    {
        #if _use_gloms
        AtomGlom glomtmp;
        int types[3] = { pi, ionic, hbond };

        int lac = ligand->get_atom_count();
        bool dirty[lac+4], dirttmp[lac+4];

        for (n=0; n<3; n++)
        {
            for (i=0; i<lac; i++) dirty[i] = false;
            for (i=0; i<lac; i++)
            {
                if (dirty[i]) continue;
                for (j=i+1; j<lac; j++) dirttmp[j] = false;

                Atom* a = ligand->get_atom(i);
                Atom* b;
                float ac = a->get_charge();
                float apol = a->is_polar();
                float api = a->is_pi() && (a->get_Z() > 1);
                float rab;

                glomtmp.atoms.clear();

                switch (types[n])
                {
                    case ionic:
                    if (ac)
                    {
                        glomtmp.atoms.push_back(a);
                        dirttmp[i] = true;

                        for (j=i+1; j<lac; j++)
                        {
                            if (dirty[j] || dirttmp[j]) continue;

                            b = ligand->get_atom(j);
                            if (b->is_pi() != api) continue;
                            rab = glomtmp.distance_to(b->get_location());
                            if ((sgn(b->get_charge()) == sgn(ac))
                                &&
                                rab < 2.6
                               )
                            {
                                glomtmp.atoms.push_back(b);
                                dirttmp[j] = true;
                            }
                        }
                    }
                    break;

                    case hbond:
                    if (apol && !ac)
                    {
                        glomtmp.atoms.push_back(a);
                        dirttmp[i] = true;

                        for (j=i+1; j<lac; j++)
                        {
                            if (dirty[j] || dirttmp[j]) continue;

                            b = ligand->get_atom(j);
                            if (b->is_pi() != api) continue;
                            if (b->get_charge()) continue;
                            rab = glomtmp.distance_to(b->get_location());
                            if (b->is_polar()
                                &&
                                (
                                    rab < 2.6
                                    ||
                                    (
                                        rab < 3.7
                                        &&
                                        !b->is_pi()
                                    )
                                )
                               )
                            {
                                glomtmp.atoms.push_back(b);
                                dirttmp[j] = true;
                            }
                        }
                    }
                    break;

                    case pi:
                    if (api)
                    {
                        glomtmp.atoms.push_back(a);
                        dirttmp[i] = true;

                        for (j=i+1; j<lac; j++)
                        {
                            if (dirty[j] || dirttmp[j]) continue;

                            b = ligand->get_atom(j);
                            if (b->get_Z() == 1) continue;
                            rab = glomtmp.distance_to(b->get_location());
                            if (b->is_pi()
                                &&
                                b->is_conjugated_to(a)
                               )
                            {
                                glomtmp.atoms.push_back(b);
                                dirttmp[j] = true;
                            }
                        }
                    }
                    break;

                    default:
                    ;
                }

                bool taken = false;
                for (j=0; j<3; j++) if (glomtmp.atoms == ligand_gloms[j].atoms) taken = true;
                if (taken) continue;

                bool too_similar = false;
                for (j=0; j<n; j++)
                {
                    if (!ligand_gloms[j].atoms.size()) continue;
                    // if (glomtmp.get_pi() >= 3.0*ligand_gloms[j].get_pi()) continue;
                    float r = ligand_gloms[j].get_center().get_3d_distance(glomtmp.get_center());
                    if (r < 2) too_similar = true;
                }
                if (too_similar) continue;

                float tsum = glomtmp.get_sum();
                if (fabs(tsum) > fabs(ligand_gloms[0].get_sum()))
                {
                    ligand_gloms[2] = ligand_gloms[1];
                    ligand_gloms[1] = ligand_gloms[0];
                    ligand_gloms[0] = glomtmp;
                    dirty[i] = true;
                    for (j=0; j<lac; j++) if (dirttmp[j]) dirty[j] = true;
                }
                else if (fabs(tsum) > fabs(ligand_gloms[1].get_sum()))
                {
                    ligand_gloms[2] = ligand_gloms[1];
                    ligand_gloms[1] = glomtmp;
                    dirty[i] = true;
                    for (j=0; j<lac; j++) if (dirttmp[j]) dirty[j] = true;
                }
                else if (fabs(tsum) > fabs(ligand_gloms[2].get_sum()))
                {
                    ligand_gloms[2] = glomtmp;
                    dirty[i] = true;
                    for (j=0; j<lac; j++) if (dirttmp[j]) dirty[j] = true;
                }
            }            
        }
        
        if (ligand_gloms[2].atoms.size())
        {
            Point a = ligand_gloms[0].get_center();
            Point b = ligand_gloms[1].get_center();
            Point c = ligand_gloms[2].get_center();

            float rab = a.get_3d_distance(b);
            float rac = a.get_3d_distance(c);

            if (rab < 0.75*rac)
            {
                glomtmp = ligand_gloms[1];
                ligand_gloms[1] = ligand_gloms[2];
                ligand_gloms[2] = glomtmp;
            }
        }

        #else
        ligbb = m.get_most_bindable(3);
        ligbbh = new Atom*[5];

        for (i=0; i<5; i++)
        {
            ligbbh[i] = nullptr;
            lig_inter_typ[i] = vdW;
        }

        for (i=0; i<3; i++)
        {
            if (!ligbb[i]) continue;
            if (ligbb[i]->get_Z() == 1)
            {
                ligbbh[i] = ligbb[i];
                ligbb[i] = ligbbh[i]->get_bond_by_idx(0)->btom;
            }
            else
            {
                ligbbh[i] = ligbb[i]->is_bonded_to("H");
            }

            if (fabs(ligbb[i]->get_charge()) >= 1
                    ||
                    ligbb[i]->get_acidbase()
                    ||
                    (ligbbh[i] && fabs(ligbbh[i]->get_charge()) >= 1)
                    ||
                    (ligbbh[i] && ligbbh[i]->get_acidbase())
               )
                lig_inter_typ[i] = ionic;
            else if (fabs(ligbb[i]->is_polar()) >= 0.5
                     ||
                     (ligbbh[i] && fabs(ligbbh[i]->is_polar()) >= 0.5)
                    )
                lig_inter_typ[i] = hbond;
            else if (ligbb[i]->is_pi())
                lig_inter_typ[i] = pi;
        }
        #endif
    }

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

        bclash += m.get_intermol_clashes(laa);

        if (met) bclash += laa->get_intermol_clashes(met);
    }
    if (met) bclash += m.get_intermol_clashes(met);
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

    if (use_bestbind_algorithm) for (i=0; i<3; i++)
        {
            #if _use_gloms
            #if _dbg_glomsel
            cout << "Ligand's ";
            if (i == 0) cout << "primary";
            else if (i == 1) cout << "secondary";
            else if (i == 2) cout << "tertiary";
            else cout << "accessory";
            #endif

            int lgsz = ligand_gloms[i].atoms.size();

            #if _dbg_glomsel
            cout << " atom group has " << lgsz << " atoms." << endl;

            if (lgsz)
            {
                for (j = 0; j < lgsz; j++)
                {
                    cout << ligand_gloms[i].atoms[j]->name << " ";
                }
                cout << endl;
            }
            #endif

            #if _dbg_glomsel
            cout << "Ionic: " << ligand_gloms[i].get_ionic()
                 << " H-bond: " << ligand_gloms[i].get_polarity()
                 << " pi-stack: " << ligand_gloms[i].get_pi()
                ;
            cout << endl << endl;
            #endif
            #else
            if (ligbb[i])
            {
                cout << "# Best binding heavy atom " << i << " of ligand" << endl << "# LBBA: " << ligbb[i]->name
                    << " type: " << lig_inter_typ[i] << endl;
                if (output) *output << "# Best binding heavy atom " << i << " of ligand" << endl << "LBBA: " << ligbb[i]->name
                    << " type: " << lig_inter_typ[i] << endl;
            }
            if (ligbbh[i])
            {
                cout << "# Best binding hydrogen " << i << " of ligand" << endl << "# LBBH: " << ligbbh[i]->name << endl;
                if (output) *output << "# Best binding hydrogen " << i << " of ligand" << endl << "LBBH: " << ligbbh[i]->name << endl;
            }
            cout << endl;
            #endif
        }

    i = poses*(triesleft+1)+8;
    DockResult dr[i][pathnodes+2];
    // cout << "dr[" << i << "] allocated." << endl;
    for (i=0; i<poses; i++) dr[i][0].kJmol = 0;

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

    #if active_persistence
    float res_kJmol[seql+8];
    for (i=0; i<seql+8; i++) res_kJmol[i] = 0;
    #endif

    found_poses = 0;
    int wrote_acvmx = -1, wrote_acvmr = -1;
_try_again:
    // srand(0xb00d1cca);
    srand(time(NULL));
    for (pose = 1; pose <= poses; pose++)
    {
        last_ttl_bb_dist = 0;
        ligand->minimize_internal_clashes();
        float lig_min_int_clsh = ligand->get_internal_clashes();
        ligand->crumple(fiftyseventh*44);

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


        delete protein;
        protein = new Protein(protfname);

        if (hydrogenate_pdb || tplset)
        {
            pf = fopen(temp_pdb_file.c_str(), "r");
            protein->load_pdb(pf);
            fclose(pf);
        }
        else
        {
            pf = fopen(protfname, "r");
            protein->load_pdb(pf);
            fclose(pf);
        }
        prepare_initb();

        ligand->recenter(pocketcen);
        // cout << "Centered ligand at " << pocketcen << endl;

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

            if (soft_pocket)
            {
                for (i=0; i<PROT_MAX_RGN; i++)
                {
                    rgnxform_r[pose][nodeno][i] = rgnxform_theta[pose][nodeno][i] = rgnxform_y[pose][nodeno][i]
                        = rgnrot_alpha[pose][nodeno][i] = rgnrot_u[pose][nodeno][i] = rgnrot_w[pose][nodeno][i]
                        = 0;
                }
    
                g_rgnxform_r = rgnxform_r[pose][nodeno];
                g_rgnxform_theta = rgnxform_theta[pose][nodeno];
                g_rgnxform_y = rgnxform_y[pose][nodeno];
                g_rgnrot_alpha = rgnrot_alpha[pose][nodeno];
                g_rgnrot_u = rgnrot_u[pose][nodeno];
                g_rgnrot_w = rgnrot_w[pose][nodeno];
            }

            if (waters)
            {
                for (i = 0; i <= omaxh2o; i++)
                {
                    waters[i] = owaters[i];
                }
                maxh2o = omaxh2o;
            }

            if (pathstrs.size() < nodeno) break;
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

                        delete[] b;
                    }
                }

                delete protein;
                protein = new Protein(protafname);
                
                pf = fopen(protafname, "r");
                protein->load_pdb(pf);
                fclose(pf);

                prepare_initb();

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

                        delete[] b;
                    }

                    delete[] sidechain_bondrots[i];
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


            if (active_helix_rots.size())
            {
                for (j=0; j<active_helix_rots.size(); j++)
                {
                    if (active_helix_rots[j].nodeno == -1) active_helix_rots[j].nodeno = active_matrix_node;
                    if (active_helix_rots[j].nodeno == nodeno)
                    {
                        float acvdirection = 1;
                        int wroteoff = 0;

                        active_helix_rots[j].start_resno  = interpret_resno(active_helix_rots[j].start_resno_str.c_str());
                        active_helix_rots[j].end_resno    = interpret_resno(active_helix_rots[j].end_resno_str.c_str());
                        active_helix_rots[j].origin_resno = interpret_resno(active_helix_rots[j].origin_resno_str.c_str());

                        if (active_helix_rots[j].axis_type == hxrat_region)
                        {
                            int sr = protein->get_region_start( active_helix_rots[j].axis_str );
                            int er = protein->get_region_end( active_helix_rots[j].axis_str );

                            Point pt = protein->get_region_center(sr, er);
                            Point rel = pt.subtract(protein->get_region_center(active_helix_rots[j].start_resno, active_helix_rots[j].end_resno));

                            rel.y = 0;
                            active_helix_rots[j].axis = rotate3D(rel, Point(0,0,0), Point(0,1,0), square);
                            #if _dbg_hxrax
                            cout << "Helix rotation set for axis " << (Point)active_helix_rots[j].axis << endl;
                            #endif
                        }
                        else if (active_helix_rots[j].axis_type == hxrat_BW)
                        {
                            const char* bwno = active_helix_rots[j].axis_str.c_str();
                            const char* bwdot = strchr(bwno, '.');
                            int bw50 = protein->get_bw50(atoi(bwno));
                            int pettia = atoi(bwdot+1);
                            int resno = bw50 + pettia - 50;

                            AminoAcid* aa = protein->get_residue(resno);
                            if (!aa)
                            {
                                cout << "Bad helix rotation axis " << bwno << endl;
                                throw 0xbad12e5;
                            }

                            Point pt = aa->get_atom_location("CA");
                            Point rel = pt.subtract(protein->get_region_center(active_helix_rots[j].start_resno, active_helix_rots[j].end_resno));

                            rel.y = 0;
                            active_helix_rots[j].axis = rotate3D(rel, Point(0,0,0), Point(0,1,0), square);
                            #if _dbg_hxrax
                            cout << "Helix rotation set for axis " << (Point)active_helix_rots[j].axis << endl;
                            #endif
                        }
                        else if (active_helix_rots[j].axis_type == hxrat_atom)
                        {
                            const char* bwno = active_helix_rots[j].axis_str.c_str();
                            const char* bwdot = strchr(bwno, '.');
                            const char* colon = strchr(bwno, ':');
                            int bw50 = protein->get_bw50(atoi(bwno));
                            int pettia = atoi(bwdot+1);
                            int resno = bw50 + pettia - 50;

                            AminoAcid* aa = protein->get_residue(resno);
                            if (!aa)
                            {
                                cout << "Bad helix rotation axis " << bwno << endl;
                                throw 0xbad12e5;
                            }

                            Point pt = aa->get_atom_location(colon+1);
                            Point rel = pt.subtract(protein->get_region_center(active_helix_rots[j].start_resno, active_helix_rots[j].end_resno));

                            rel.y = 0;
                            active_helix_rots[j].axis = rotate3D(rel, Point(0,0,0), Point(0,1,0), -square);
                            #if _dbg_hxrax
                            cout << "Helix rotation set for axis " << (Point)active_helix_rots[j].axis << endl;
                            #endif
                        }
                        else if (active_helix_rots[j].axis_type == hxrat_helical)
                        {
                            Point Navg, Cavg;
                            int Ndiv=0, Cdiv=0;
                            int lres;

                            for (lres = active_helix_rots[j].start_resno+3; lres >= active_helix_rots[j].start_resno; lres--)
                            {
                                AminoAcid* aa = protein->get_residue(lres);
                                if (aa)
                                {
                                    Navg = Navg.add(aa->get_CA_location());
                                    Ndiv++;
                                }
                            }
                            if (!Ndiv)
                            {
                                cout << "Not enough N terminus residues for helical rotation." << endl;
                                throw 0xbad1207;
                            }
                            else Navg.scale(Navg.magnitude() / Ndiv);

                            for (lres = active_helix_rots[j].end_resno-3; lres <= active_helix_rots[j].end_resno; lres++)
                            {
                                AminoAcid* aa = protein->get_residue(lres);
                                if (aa)
                                {
                                    Cavg = Cavg.add(aa->get_CA_location());
                                    Cdiv++;
                                }
                            }
                            if (!Cdiv)
                            {
                                cout << "Not enough C terminus residues for helical rotation." << endl;
                                throw 0xbad1207;
                            }
                            else Cavg.scale(Cavg.magnitude() / Cdiv);

                            active_helix_rots[j].axis = Cavg.subtract(Navg);
                        }

                        int sr = active_helix_rots[j].start_resno;
                        int er = active_helix_rots[j].end_resno;
                        int mr = active_helix_rots[j].origin_resno;

                        protein->move_piece(sr, er,
                                (acvdirection > 0)
                                ? protein->get_region_center(sr, er).add(active_helix_rots[j].transform)
                                : protein->get_region_center(sr, er).subtract(active_helix_rots[j].transform)
                            );

                        protein->rotate_piece(sr, er, protein->get_atom_location(mr, "CA"),
                            active_helix_rots[j].axis, active_helix_rots[j].theta*acvdirection);

                        protein->get_internal_clashes(active_helix_rots[j].start_resno, active_helix_rots[j].end_resno, repack_on_hxr);

                        if (wrote_acvmr < (j+wroteoff) )
                        {
                            Point ptaxis = active_helix_rots[j].axis;
                            Point ptorigin = protein->get_atom_location(mr, "CA");
                            // Write an active matrix to the dock.
                            cout << "ACR " << nodeno << " " << active_helix_rots[j].regname << " " << sr << " " << er << " "
                                << active_helix_rots[j].transform.x*acvdirection << " "
                                << active_helix_rots[j].transform.y*acvdirection << " "
                                << active_helix_rots[j].transform.z*acvdirection << " "
                                << ptorigin.x << " " << ptorigin.y << " " << ptorigin.z << " "
                                << ptaxis.x << " " << ptaxis.y << " " << ptaxis.z << " "
                                << active_helix_rots[j].theta*acvdirection << endl;
                            if (output) *output << "ACR " << nodeno << " " << active_helix_rots[j].regname << " " << sr << " " << er << " "
                                << active_helix_rots[j].transform.x*acvdirection << " "
                                << active_helix_rots[j].transform.y*acvdirection << " "
                                << active_helix_rots[j].transform.z*acvdirection << " "
                                << ptorigin.x << " " << ptorigin.y << " " << ptorigin.z << " "
                                << ptaxis.x << " " << ptaxis.y << " " << ptaxis.z << " "
                                << active_helix_rots[j].theta*acvdirection << endl;

                            wrote_acvmr = j+wroteoff;
                        }
                    }
                }
            }

            if (nodeno == active_matrix_node
                ||
                (nodeno == deactivate_node && active_matrix_node >= 0 && deactivate_node > active_matrix_node)
               )
            {
                float acvdirection = (nodeno == active_matrix_node) ? 1 : -1;
                int wroteoff = (nodeno == deactivate_node) ? 100 : 0;

                #if internode_momentum_only_on_activation 
                conformer_momenta_multiplier = nodeno ? internode_momentum_mult : 1;
                #endif

                #if prevent_ligand_360_on_activate
                allow_ligand_360_tumble = allow_ligand_360_flex = false;
                #endif

                #if active_persistence
                j=0;
                for (i=1; i<=seql; i++)
                {
                    #if _DBG_RESBMULT
                    cout << "res_kJmol[" << i << "] = " << res_kJmol[i] << endl << flush;
                    #endif
                    if (res_kJmol[i] >= active_persistence_threshold)
                    {
                        active_persistence_resno[j] = i;
                        #if _DBG_RESBMULT
                        cout << "Added " << i << " to active persistence resnos." << endl << flush;
                        #endif
                        j++;
                        if (j >= active_persistence_limit) break;
                    }
                }
                #endif

                #if active_persistence && active_persistence_noflex
                allow_ligand_flex = false;
                #endif

                #if active_persistence_follow
                Point residue_follow[active_persistence_limit];
                for (i=0; active_persistence_resno[i]; i++)
                {
                    Atom* CA = protein->get_atom(active_persistence_resno[i], "CA");
                    if (!CA) residue_follow[i] = Point(0,0,0);
                    else residue_follow[i] = CA->get_location();
                }
                #endif

                // Each TMR:
                for (i=1; i<=active_matrix_count; i++)
                {
                    std::string regname = (std::string)"TMR" + std::to_string(i);
                    int sr = protein->get_region_start(regname), er = protein->get_region_end(regname);

                    if (!sr || !er)
                    {
                        cout << "Cannot activate " << regname << "; region not found in protein." << endl << flush;
                        throw 0xbadf12e;
                    }

                    LocRotation lrot, nlrot, clrot;

                    if (active_matrix_type == 6)
                    {
                        // Get the CA location of the C-terminus of the helix;
                        // Add the active_matrix_c to the result;
                        Point calign;
                        calign = (acvdirection > 0)
                            ? protein->get_atom_location(er, "CA").add(active_matrix_c[i])
                            : protein->get_atom_location(er, "CA").subtract(active_matrix_c[i])
                            ;

                        // Move the entire helix by the values of active_matrix_n.
                        protein->move_piece(sr, er,
                            (acvdirection > 0)
                            ? active_matrix_n[i].add(protein->get_region_center(sr, er))
                            : active_matrix_n[i].subtract(protein->get_region_center(sr, er))
                            );

                        // Call protein->rotate_piece() to align the C-terminus residue with the result, using the N-terminus residue as the pivot res.
                        lrot = protein->rotate_piece(sr, er, er, calign, sr);
                        lrot.v.r = 1;
                    }
                    else if (active_matrix_type == 9)
                    {
                        protein->move_piece(sr, er,
                            (acvdirection > 0)
                            ? protein->get_region_center(sr, er).add(active_matrix_m[i])
                            : protein->get_region_center(sr, er).subtract(active_matrix_m[i])
                            );
                        int half = (sr + er) / 2;

                        Point calign, nalign;
                        if (acvdirection > 0)
                        {
                            calign = protein->get_atom_location(er, "CA").add(active_matrix_c[i]);
                            nalign = protein->get_atom_location(sr, "CA").add(active_matrix_n[i]);
                        }
                        else
                        {
                            calign = protein->get_atom_location(er, "CA").subtract(active_matrix_c[i]);
                            nalign = protein->get_atom_location(sr, "CA").subtract(active_matrix_n[i]);
                        }

                        nlrot = protein->rotate_piece(sr, half, sr, nalign, half);
                        clrot = protein->rotate_piece(half, er, er, calign, half);
                        nlrot.v.r = 1;
                        clrot.v.r = 1;
                    }

                    if (wrote_acvmx < (i+wroteoff))
                    {
                        #if write_activation_matrix
                        // Output the activation for the viewer to recognize.
                        cout << "ACM " << nodeno << " " << regname << " " << sr << " " << er << " "
                            << active_matrix_n[i].x*acvdirection << " " << active_matrix_n[i].y*acvdirection << " " << active_matrix_n[i].z*acvdirection << " "
                            << active_matrix_c[i].x*acvdirection << " " << active_matrix_c[i].y*acvdirection << " " << active_matrix_c[i].z*acvdirection
                            << endl;
                        if (output) *output << "ACM " << nodeno << " " << regname << " " << sr << " " << er << " "
                            << active_matrix_n[i].x*acvdirection << " " << active_matrix_n[i].y*acvdirection << " " << active_matrix_n[i].z*acvdirection << " "
                            << active_matrix_c[i].x*acvdirection << " " << active_matrix_c[i].y*acvdirection << " " << active_matrix_c[i].z*acvdirection
                            << endl;
                        #endif

                        #if write_active_rotation
                        if (active_matrix_type == 9)
                        {
                            int half = (sr + er) / 2;

                            Point nlrv(nlrot.v);
                            cout << "ACR " << nodeno << " " << regname << "n " << sr << " " << half << " "
                                << active_matrix_m[i].x << " " << active_matrix_m[i].y << " " << active_matrix_m[i].z << " "
                                << nlrot.origin.x << " " << nlrot.origin.y << " " << nlrot.origin.z << " "
                                << nlrv.x << " " << nlrv.y << " " << nlrv.z << " "
                                << nlrot.a*acvdirection << endl;
                            if (output) *output << "ACR " << nodeno << " " << regname << "n " << sr << " " << half << " "
                                << active_matrix_m[i].x << " " << active_matrix_m[i].y << " " << active_matrix_m[i].z << " "
                                << nlrot.origin.x << " " << nlrot.origin.y << " " << nlrot.origin.z << " "
                                << nlrv.x << " " << nlrv.y << " " << nlrv.z << " "
                                << nlrot.a*acvdirection << endl;

                            Point clrv(clrot.v);
                            cout << "ACR " << nodeno << " " << regname << "c " << half << " " << er << " "
                                << active_matrix_m[i].x << " " << active_matrix_m[i].y << " " << active_matrix_m[i].z << " "
                                << clrot.origin.x << " " << clrot.origin.y << " " << clrot.origin.z << " "
                                << clrv.x << " " << clrv.y << " " << clrv.z << " "
                                << clrot.a*acvdirection << endl;
                            if (output) *output << "ACR " << nodeno << " " << regname << "c " << half << " " << er << " "
                                << active_matrix_m[i].x << " " << active_matrix_m[i].y << " " << active_matrix_m[i].z << " "
                                << clrot.origin.x << " " << clrot.origin.y << " " << clrot.origin.z << " "
                                << clrv.x << " " << clrv.y << " " << clrv.z << " "
                                << clrot.a*acvdirection << endl;
                        }
                        else
                        {
                            Point lrv(lrot.v);
                            cout << "ACR " << nodeno << " " << regname << " " << sr << " " << er << " "
                                << active_matrix_n[i].x << " " << active_matrix_n[i].y << " " << active_matrix_n[i].z << " "
                                << lrot.origin.x << " " << lrot.origin.y << " " << lrot.origin.z << " "
                                << lrv.x << " " << lrv.y << " " << lrv.z << " "
                                << lrot.a*acvdirection << endl;
                            if (output) *output << "ACR " << nodeno << " " << regname << " " << sr << " " << er << " "
                                << active_matrix_n[i].x << " " << active_matrix_n[i].y << " " << active_matrix_n[i].z << " "
                                << lrot.origin.x << " " << lrot.origin.y << " " << lrot.origin.z << " "
                                << lrv.x << " " << lrv.y << " " << lrv.z << " "
                                << lrot.a*acvdirection << endl;
                        }
                        #endif

                        wrote_acvmx = i+wroteoff;
                    }
                }

                // If there are any active bond rotations, perform them but ensure the angle is relative to the *original* position
                // from the PDB file.
                if (active_bond_rots.size())
                {
                    for (j=0; j<active_bond_rots.size(); j++)
                    {
                        if (active_bond_rots[j].bond)
                        {
                            active_bond_rots[j].bond->rotate(active_bond_rots[j].theta*acvdirection - active_bond_rots[j].bond->total_rotations);
                        }
                    }
                }


                #if active_persistence_follow
                for (i=0; active_persistence_resno[i]; i++)
                {
                    Atom* CA = protein->get_atom(active_persistence_resno[i], "CA");
                    if (CA) residue_follow[i] = CA->get_location().subtract(residue_follow[i]);
                }

                if (i)
                {
                    Point catch_up = average_of_points(residue_follow, i);
                    ligand->move(catch_up);
                }
                #endif

                #if save_active_protein
                if (pose == 1 && acvdirection>0)
                {
                    FILE* f = fopen("tmp/active.pdb", "wb");
                    if (f)
                    {
                        protein->save_pdb(f);
                        protein->end_pdb(f);
                        fclose(f);
                    }
                }
                #endif

                // TODO: #if !recenter_ligand_each_node, recenter the ligand here to keep up with the residues that it was coordinated to.
                // Perhaps also multimol it for 10 or so iterations with all flexions (ligand and residue) globally disabled.
            }

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
                if (!strcmp(words[2], "RES"))
                {
                    extra_wt.clear();
                    for (i=2; words[i]; i++)
                    {
                        extra_wt.push_back(interpret_resno(words[i]));
                    }
                }

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

            #if redo_tumble_spheres_on_activation
            if (nodeno == active_matrix_node)
            {
                if (!use_bestbind_algorithm && !use_prealign) do_tumble_spheres(ligcen_target);
            }
            #endif

            #if redo_tumble_spheres_every_node
            if (!use_bestbind_algorithm && !use_prealign && (!prevent_ligand_360_on_activate || (nodeno != active_matrix_node))) do_tumble_spheres(ligcen_target);
            #endif

            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Saved last nodecen." << endl;
            #endif

            #if recenter_ligand_each_node
            // Move the ligand to the new node center.
            m.recenter(nodecen);
            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Molecule recenter (or not)." << endl;
            #endif
            m.reset_conformer_momenta();
            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Conformer momenta reset." << endl;
            #endif
            #endif

            sphres = protein->get_residues_can_clash_ligand(reaches_spheroid[nodeno], &m, nodecen, size, addl_resno);
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
                            mvaa->movability = MOV_NONE;
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

                            // If residue is within any active_helix_rots region, increase the odds.
                            if (active_helix_rots.size())
                            {
                                int l, n = active_helix_rots.size(), resno = reaches_spheroid[nodeno][i]->get_residue_no();
                                for (l=0; l<n; l++)
                                {
                                    if (resno >= active_helix_rots[l].start_resno && resno <= active_helix_rots[l].end_resno)
                                    {
                                        weight *= 1.5;
                                    }
                                }
                            }

                            if (extra_wt.size())
                            {
                                int l, n = extra_wt.size(), resno = reaches_spheroid[nodeno][i]->get_residue_no();
                                for (l=0; l<n; l++)
                                {
                                    if (resno == extra_wt[l])
                                    {
                                        weight *= 20;
                                        #if _dbg_flexion_selection
                                        // cout << resno << " boosted." << endl;
                                        #endif
                                    }
                                }
                            }

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

                if (forced_flexible_resnos.size())
                {
                    for (i=0; i<forced_flexible_resnos.size(); i++)
                    {
                        forced_flexible_resnos[i].resolve_resno(protein);
                        AminoAcid* mvaa = protein->get_residue(forced_flexible_resnos[i].resno);
                        if (mvaa)
                        {
                            mvaa->movability = MOV_FORCEFLEX;
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

            /*cout << "Dock residues for node " << nodeno << ": " << endl;
            if (output) *output << "Dock residues for node " << nodeno << ": " << endl;
            for (i=0; i<sphres; i++)
            {	cout << *reaches_spheroid[nodeno][i] << " ";
            	if (output) *output << *reaches_spheroid[nodeno][i] << " ";
            }
            cout << endl << endl;
            if (output) *output << endl << endl;*/


            if (!nodeno)
            {
                float alignment_distance[5];
                for (l=0; l<3; l++)
                {
                    #if !_use_gloms
                    alignment_aa[l]=0;
                    #endif
                    alignment_distance[l]=0;
                }

                for (l=0; l<=pathnodes; l++)
                {
                    dr[drcount][l].metric = 0;
                }

                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Initialize null AA pointer." << endl;
                #endif

                // Best-Binding Algorithm
                // Find a binding pocket feature with a strong potential binding to the ligand.
                std::string alignment_name = "";
                if (use_bestbind_algorithm)
                {
                    for (l=0; l<3; l++)
                    {
                        sc_gloms[l].aminos.clear();
                        #if _use_gloms
                        // Find the strongest loneliest suitable side chain glom for ligand_gloms[l] and set sc_gloms[l] to equal it.
                        ResidueGlom glomtmp;
                        for (i=0; reaches_spheroid[nodeno][i]; i++)
                        {
                            glomtmp.aminos.clear();
                            if (ligand_gloms[l].compatibility(reaches_spheroid[nodeno][i]))
                            {
                                #if _dbg_glomsel
                                cout << "Considering " << reaches_spheroid[nodeno][i]->get_name() << " for glom " << l << "..." << endl;
                                #endif
                                glomtmp.aminos.push_back(reaches_spheroid[nodeno][i]);
                                for (j=i+1; reaches_spheroid[nodeno][j]; j++)
                                {
                                    if (!ligand_gloms[l].compatibility(reaches_spheroid[nodeno][j])) continue;
                                    Atom* cb = reaches_spheroid[nodeno][j]->get_atom("CB");
                                    if (!cb) continue;
                                    if (frand(0,1) < bb_stochastic) continue;                       // stochastic component.
                                    float r = glomtmp.distance_to(cb->get_location());
                                    if (reaches_spheroid[nodeno][j]->hydrophilicity() >= 0.333)
                                        r -= (reaches_spheroid[nodeno][j]->get_reach() + 2.5);
                                    else
                                        r -= 3.5;
                                    if (r < ligand_gloms[l].bounds())
                                    {
                                        glomtmp.aminos.push_back(reaches_spheroid[nodeno][j]);
                                        #if _dbg_glomsel
                                        cout << "Adding " << reaches_spheroid[nodeno][j]->get_name() << " to candidate." << endl;
                                        #endif
                                    }
                                }

                                bool too_similar = false;
                                for (j=0; j<l; j++)
                                {
                                    if (sc_gloms[j].aminos == glomtmp.aminos) too_similar = true;
                                    if (sc_gloms[j].aminos.size() == 1
                                        &&
                                        std::find(glomtmp.aminos.begin(), glomtmp.aminos.end(), sc_gloms[j].aminos[0]) != glomtmp.aminos.end()
                                    )
                                        too_similar = true;
                                    if (glomtmp.aminos.size() == 1
                                        &&
                                        std::find(sc_gloms[j].aminos.begin(), sc_gloms[j].aminos.end(), glomtmp.aminos[0]) != sc_gloms[j].aminos.end()
                                    )
                                        too_similar = true;
                                    // if (sc_gloms[j].get_center().get_3d_distance(glomtmp.get_center()) < 0.01) too_similar = true;
                                    float rlg = ligand_gloms[j].get_center().get_3d_distance(ligand_gloms[l].get_center());
                                    float rsg = ligand_gloms[j].get_center().get_3d_distance(glomtmp.get_center());
                                    if (rsg < 0.9 * rlg) too_similar = true;
                                }

                                #if _dbg_glomsel
                                if (too_similar) cout << "Too similar to an existing glom." << endl << endl;
                                #endif

                                if (!too_similar)
                                {
                                    float r = glomtmp.distance_to(loneliest);
                                    float rr = sc_gloms[l].distance_to(loneliest);
                                    float rg = ligand_gloms[l].distance_to(ligand->get_barycenter());

                                    r -= rg; if (r < 1) r = 1;
                                    rr -= rg; if (rr < 1) rr = 1;

                                    float tcptbl = glomtmp.compatibility(&ligand_gloms[l]) / (r*r);
                                    float ptcptbl = sc_gloms[l].compatibility(&ligand_gloms[l]) / (rr*rr);
                                    if (!i || tcptbl > ptcptbl)
                                    {
                                        sc_gloms[l] = glomtmp;
                                        #if _dbg_glomsel
                                        cout << "Accepted with compatibility " << glomtmp.compatibility(&ligand_gloms[l]) << " over distance " << r
                                            << " squared = " << tcptbl << "." << endl << endl;
                                        #endif
                                    }
                                    #if _dbg_glomsel
                                    else
                                    {
                                        cout << "Rejected because compatibility " << glomtmp.compatibility(&ligand_gloms[l]) << " over distance " << r
                                            << " squared is " << tcptbl << " not greater than previous value of " << ptcptbl << "." << endl << endl;
                                    }
                                    #endif
                                }
                            }
                        }

                        #else
                        retain_bindings[l].cardinality = 0;
                        if (!ligbb[l]) continue;
                        retain_bindings[l].atom = ligbb[l];
                        ligand->springy_bonds = retain_bindings;
                        ligand->springy_bondct = l+1;
                        float alignment_potential = 0;
                        for (i=0; reaches_spheroid[nodeno][i]; i++)
                        {
                            if (!protein->aa_ptr_in_range(reaches_spheroid[nodeno][i]))
                            {
                                reaches_spheroid[nodeno][i] = NULL;
                                continue;
                            }

                            if (l && reaches_spheroid[nodeno][i] == alignment_aa[l-1]) continue;
                            if (l>1 && reaches_spheroid[nodeno][i] == alignment_aa[l-2]) continue;

                            if (lig_inter_typ[l] == vdW && reaches_spheroid[nodeno][i]->hydrophilicity() >= 0.3) continue;
                            if (lig_inter_typ[l] == hbond && reaches_spheroid[nodeno][i]->hydrophilicity() < 0.2) continue;
                            if (lig_inter_typ[l] == ionic && reaches_spheroid[nodeno][i]->hydrophilicity() < 0.3) continue;

                            if (reaches_spheroid[nodeno][i]->is_glycine()) continue;

                            #if _DBG_STEPBYSTEP
                            if (debug)
                            {
                                *debug << "Check capable of inter (" << i << ") ";
                                *debug << lig_inter_typ[l];
                                *debug << flush;
                                Star s;
                                s.paa = reaches_spheroid[nodeno][i];
                                *debug << *s.paa << " " << flush;
                                *debug << *reaches_spheroid[nodeno][i];
                                *debug << endl;
                            }
                            #endif

                            float pottmp = reaches_spheroid[nodeno][i]->get_atom_mol_bind_potential(ligbb[l]);
                            if (ligbbh[l]) pottmp += reaches_spheroid[nodeno][i]->get_atom_mol_bind_potential(ligbbh[l]);

                            if (lig_inter_typ[l] == hbond)
                            {
                                pottmp *= (1.0 + reaches_spheroid[nodeno][i]->hydrophilicity() * hydrophilicity_boost);
                            }
                            else if (lig_inter_typ[l] == pi || lig_inter_typ[l] == vdW)
                            {
                                pottmp *= fmax(0, 1.0 - reaches_spheroid[nodeno][i]->hydrophilicity() * hydrophilicity_boost);
                            }

                            pottmp *= (1.0 + frand(-best_binding_stochastic, best_binding_stochastic) );

                            if (extra_wt.size()
                                    &&
                                    std::find(extra_wt.begin(), extra_wt.end(), reaches_spheroid[nodeno][i]->get_residue_no())!=extra_wt.end()
                               )
                            {
                                pottmp *= 1.25;		// Extra weight for residues mentioned in a CEN RES or PATH RES parameter.
                            }
                            else
                            {
                                pottmp /= pocketcen.get_3d_distance(reaches_spheroid[nodeno][i]->get_barycenter());
                            }
                            // cout << reaches_spheroid[nodeno][i]->get_3letter() << reaches_spheroid[nodeno][i]->get_residue_no() << " " << pottmp << endl;
                            Atom* coi = reaches_spheroid[nodeno][i]->capable_of_inter(lig_inter_typ[l]);
                            if (coi
                                    &&
                                    (	!alignment_aa[l]
                                        ||
                                        pottmp > alignment_potential
                                        ||
                                        ((pose>1) && pottmp > (0.9 * alignment_potential) && !(rand() % sphres))
                                    )
                               )
                            {
                                alignment_aa[l] = reaches_spheroid[nodeno][i];
                                // alignment_name += std::to_string("|") + std::to_string(reaches_spheroid[nodeno][i]->get_residue_no());
                                alignment_potential = pottmp;
                                alignment_distance[l] = potential_distance;
                                retain_bindings[l].btom = coi;
                                retain_bindings[l].cardinality = 0.25;
                                retain_bindings[l].type = lig_inter_typ[l];
                                
                                bool opml_known = false;
                                try
                                {
                                    retain_bindings[l].optimal_radius = InteratomicForce::coordinate_bond_radius(ligbb[l], coi, lig_inter_typ[l]);
                                    retain_bindings[l].optimal_radius = InteratomicForce::coordinate_bond_radius(
                                        retain_bindings[l].atom,
                                        retain_bindings[l].btom,
                                        retain_bindings[l].type
                                        );
                                    opml_known = true;
                                }
                                catch (int ex)
                                {
                                    ;
                                }
                                if (!opml_known && ligbbh[l] && (retain_bindings[l].btom->get_Z() > 1 || retain_bindings[l].type == vdW))
                                {
                                    try
                                    {
                                        retain_bindings[l].optimal_radius = InteratomicForce::coordinate_bond_radius(
                                            ligbbh[l],
                                            retain_bindings[l].btom,
                                            retain_bindings[l].type
                                            );
                                        retain_bindings[l].atom = ligbbh[l];
                                        opml_known = true;
                                    }
                                    catch (int ex)
                                    {
                                        retain_bindings[l].optimal_radius = 2.5;
                                    }
                                }

                                alignment_aa[l]->springy_bonds = retain_bindings;
                                alignment_aa[l]->springy_bondct = l+1;
                            }
                            #if _DBG_STEPBYSTEP
                            if (debug) *debug << "Candidate alignment AA." << endl;
                            #endif
                        }
                        #endif
                        _found_alignaa:
                        ;
                    }
                }
                
                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Selected an alignment AA." << endl;
                #endif

                if (use_bestbind_algorithm && met)
                {
                    #if !_use_gloms
                    alignment_aa[2] = alignment_aa[1];
                    alignment_aa[1] = alignment_aa[0];
                    alignment_aa[0] = met;
                    #endif
                    // alignment_name = "metal";
                }
                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Alignment AA." << endl;
                #endif

                if (use_bestbind_algorithm)
                {
                    ligand->recenter(loneliest);
                    for (l=0; l<_bb_maxglom; l++)
                    {
                        #if _use_gloms
                        Point xform;
                        Point zcen;
                        Point axis;
                        LocatedVector lv;
                        Rotation rot;
                        float theta, clash;

                        Atom* alca;
                        Atom** alcaa;
                        #if !flexion_selection
                        if (flex && sc_gloms[l].aminos.size())
                        {
                            for (i=0; i<sc_gloms[l].aminos.size(); i++)
                            {
                                if (sc_gloms[l].aminos[i]->movability >= MOV_FLEXONLY)
                                {
                                    Bond** rbb = sc_gloms[l].aminos[i]->get_rotatable_bonds();
                                    if (rbb)
                                    {
                                        alcaa = sc_gloms[l].aminos[i]->get_most_bindable(1, ligand_gloms[l].atoms[0]);
                                        if (!alcaa) continue;
                                        alca = alcaa[0];
                                        delete alcaa;

                                        if (alca)
                                        {
                                        	float br = alca->get_location().get_3d_distance(loneliest);
		                                    float brad = 0;
		                                    for(j=0; rbb[j]; j++)
		                                    {
		                                        float rad;
		                                        for (rad=0; rad<M_PI*2; rad += square/4)
		                                        {
		                                            rbb[j]->rotate(square);
		                                            float lr = alca->get_location().get_3d_distance(loneliest);
		                                            if (lr < br)
		                                            {
		                                                brad = rad;
		                                                br = lr;
		                                            }
		                                        }
		                                        if (brad) rbb[j]->rotate(brad);
		                                    }
		                                }
                                    }
                                }
                            }
                        }
                        #endif

                        switch (l)
                        {
                            case 0:
                            // Move ligand to center ligand_gloms[0] at the center of sc_gloms[0].
                            // If there is only one residue in the sc glom, then move the ligand
                            // 2A towards loneliest.
                            n = sc_gloms[l].aminos.size();
                            if (!n) goto _deadglob;
                            #if _dbg_glomsel
                            cout << "Moving primary atom group to vicinity of";
                            for (i=0; i<n; i++) cout << " " << sc_gloms[l].aminos[i]->get_3letter() << sc_gloms[l].aminos[i]->get_residue_no();
                            cout << "." << endl;
                            #endif

                            xform = sc_gloms[l].get_center();
                            if (n < 4)
                            {
                                Point ptmp = loneliest.subtract(xform);
                                ptmp.scale(n+1);
                                xform = xform.add(ptmp);
                            }
                            xform = xform.subtract(ligand_gloms[l].get_center());

                            ligand->movability = MOV_ALL;
                            ligand->move(xform);
                            
                            // Slowly back the ligand away from whatever it may be clashing into.
                            reaches_spheroid[nodeno][sphres] = nullptr;
                            clash = ligand->get_intermol_clashes(reinterpret_cast<Molecule**>(reaches_spheroid[nodeno]));
                            n = 15;
                            while (clash > 50)
                            {
                                xform = loneliest.subtract(ligand->get_barycenter());
                                xform.scale(0.1);
                                ligand->move(xform);
                                n--;
                                if (!n) break;          // Prevent infinite loops.
                            }
                            break;

                            case 1:
                            // Rotate ligand about ligand_gloms[0] center to get ligand_gloms[1] center
                            // as close as possible to sc_gloms[1] center.
                            n = sc_gloms[l].aminos.size();
                            if (!n)
                            {
                                #if _dbg_glomsel
                                cout << "Aligning center of ligand towards center of pocket.";
                                #endif
                                zcen = ligand_gloms[0].get_center();
                                rot = align_points_3d(ligand->get_barycenter(), loneliest, zcen);
                                lv = rot.v;
                                lv.origin = zcen;
                                ligand->rotate(lv, rot.a);
                                goto _deadglob;
                            }
                            #if _dbg_glomsel
                            cout << "Aligning secondary atom group towards";
                            for (i=0; i<n; i++) cout << " " << sc_gloms[l].aminos[i]->get_3letter() << sc_gloms[l].aminos[i]->get_residue_no();
                            cout << "." << endl;
                            #endif

                            zcen = ligand_gloms[0].get_center();
                            rot = align_points_3d(ligand_gloms[l].get_center(), sc_gloms[l].get_center(), zcen);
                            lv = rot.v;
                            lv.origin = zcen;
                            ligand->rotate(lv, rot.a);
                            break;

                            case 2:
                            // "Rotisserie" rotate ligand about the imaginary line between
                            // ligand_gloms[0] center and ligand_gloms[1] center, to bring
                            // ligand_gloms[2] center as close as possible to sc_gloms[2] center.
                            n = sc_gloms[l].aminos.size();
                            if (!n) goto _deadglob;
                            #if _dbg_glomsel
                            cout << "\"Rotisserie\"-aligning tertiary atom group towards";
                            for (i=0; i<n; i++) cout << " " << sc_gloms[l].aminos[i]->get_3letter() << sc_gloms[l].aminos[i]->get_residue_no();
                            cout << "." << endl;
                            #endif

                            zcen = ligand_gloms[0].get_center();
                            axis = ligand_gloms[1].get_center().subtract(zcen);
                            lv = (SCoord)axis;
                            lv.origin = zcen;
                            theta = find_angle_along_vector(ligand_gloms[l].get_center(), sc_gloms[l].get_center(), zcen, axis);
                            ligand->rotate(lv, theta);
                            break;

                            default:
                            ;
                        }

                        _deadglob:
                        ;
                        #else
                        if (alignment_aa[l])
                        {
                            cout << "# Aligning " << ligbb[l]->name << " to " << alignment_aa[l]->get_name() << "..." << endl;
                            Atom* alca;

                            if (retain_bindings[l].btom) alca = retain_bindings[l].btom;
                            else
                            {
                                if (alignment_aa[l] == met)
                                    alca = alignment_aa[l]->get_nearest_atom(ligbb[l]->get_location());
                                else
                                {
                                    Atom** mbb = alignment_aa[l]->get_most_bindable(1);
                                    alca = mbb[0];
                                    delete mbb;             // Delete the pointer array, but not the pointers.
                                }
                            }
                            #if _DBG_STEPBYSTEP
                            if (debug) *debug << "Got alignment atom." << endl;
                            #endif

                            if (alca)
                            {
                                #if _preflex_alignment_res
                                // If alignment aa is not pinned, and flexion is enabled,
                                // flex alignment aa so that alca is nearest to loneliest.
                                if (flex && alignment_aa[l]->movability >= MOV_FLEXONLY)
                                {
                                    Bond** rbb = alignment_aa[l]->get_rotatable_bonds();
                                    if (rbb)
                                    {
                                        float br = alca->get_location().get_3d_distance(loneliest);
                                        float brad = 0;
                                        for(j=0; rbb[j]; j++)
                                        {
                                            float rad;
                                            for (rad=0; rad<M_PI*2; rad += square)
                                            {
                                                rbb[j]->rotate(square);
                                                float lr = alca->get_location().get_3d_distance(loneliest);
                                                if (lr < br)
                                                {
                                                    brad = rad;
                                                    br = lr;
                                                }
                                            }
                                            if (brad) rbb[j]->rotate(brad);
                                        }
                                    }
                                }
                                #endif

                                Point pt, al, cen;
                                al	= alca->get_location();

                                cen	= (l==1) ? ligbb[0]->get_location() : m.get_barycenter();

                                pt	= ligbb[l]->get_location();
                                if (ligbbh[l])
                                {
                                    Point pth = ligbbh[l]->get_location();
                                    pt.x += 0.5*(pth.x-pt.x);
                                    pt.y += 0.5*(pth.y-pt.y);
                                    pt.z += 0.5*(pth.z-pt.z);
                                }
                                

                                Rotation rot;
                                Point origin = ligbb[0]->get_location();
                                SCoord axis;
                                LocatedVector lv;
                                float theta;
                                float besttheta = 0, bestr = 100000;
                                float rstep = fiftyseventh*30;
                                switch (l)
                                {
                                case 0:
                                    // Pivot about molcen.
                                    /*append_dummy(pt);
                                    append_dummy(al);
                                    append_dummy(cen);*/
                                    rot = align_points_3d(&pt, &al, &cen);
                                    m.rotate(&rot.v, rot.a);
                                    ligand->recenter(cen);
                                    #if _dbg_bb_rots
                                    cout << "# Pivoted ligand " << (rot.a*fiftyseven) << "deg about ligand molcen." << endl;
                                    #endif

                                    if (false && !l)
                                    {
                                        Point ptr = alca->get_location().subtract(pt);
                                        SCoord v(ptr);
                                        v.r -= alignment_distance[l];
                                        v.r *= 0.5;
                                        m.move(v);
                                        cout << "# Moved ligand " << v.r << "A towards " << alignment_aa[l]->get_name()
                                             << ":" << alca->name << "." << endl;
                                    }
                                    break;

                                case 1:
                                    // Pivot about bb0.
                                    origin = ligbb[0]->get_location();
                                    lv = compute_normal(pt, al, origin);
                                    lv.origin = origin;
                                    rot.a = find_angle_along_vector(pt, al, origin, (SCoord)lv);
                                    m.rotate(lv, rot.a);
                                    #if _dbg_bb_rots
                                    cout << "# Pivoted ligand " << (rot.a*fiftyseven) << "deg about ligand " << ligbb[0]->name << "." << endl;
                                    #endif
                                    break;

                                case 2:
                                    // Rotisserie.
                                    axis = ligbb[1]->get_location().subtract(origin);

                                    for (theta=0; theta < M_PI*2; theta += rstep)
                                    {
                                        lv.copy(axis);
                                        lv.origin = origin;
                                        m.rotate(lv, rstep);

                                        float r2 = alca->get_location().get_3d_distance(ligbb[l]->get_location());
                                        if (r2 < bestr)
                                        {
                                            bestr = r2;
                                            besttheta = theta;
                                        }
                                    }

                                    lv.copy(axis);
                                    lv.origin = origin;
                                    m.rotate(lv, besttheta);
                                    ligand->recenter(cen);
                                    #if _dbg_bb_rots
                                    cout << "# Pivoted ligand " << (besttheta*fiftyseven) << "deg about ligand " << ligbb[0]->name
                                         << "-" << ligbb[1]->name << " axis." << endl;
                                    #endif
                                    break;

                                default:
                                    ;
                                }

                                #if preemptively_minimize_intermol_clashes
                                Molecule* mtmp[3];
                                mtmp[0] = &m;
                                mtmp[1] = flex ? alignment_aa[l] : nullptr;
                                mtmp[2] = nullptr;
                                m.movability = MOV_FLEXONLY;
                                alignment_aa[l]->movability = MOV_FLEXONLY;
                                Molecule::multimol_conform(mtmp);
                                m.movability = MOV_ALL;
                                // m.intermol_conform_norecen(alignment_aa[l], iters, reaches_spheroid[nodeno]);
                                // alignment_aa[l]->intermol_conform_norecen(&m, iters, reaches_spheroid[nodeno]);
                                if (debug) *debug << "Alignment atom " << l << " is "
                                                      << alignment_aa[l]->get_name() << ":" << alca->name
                                                      << " Z " << alca->get_Z() << endl;
                                #endif

                            }
                        }
                        #endif
                    }

                    #if enforce_no_bb_pullaway && _use_gloms
                    last_ttl_bb_dist = 0;
                    for (l=0; l<_bb_maxglom; l++)
                    {
                        if (ligand_gloms[l].atoms.size() && sc_gloms[l].aminos.size())
                        {
                            float r = ligand_gloms[l].distance_to(sc_gloms[l].get_center());
                            if (r < 2.5) r = 2.5;
                            last_ttl_bb_dist += r;
                        }
                    }
                    pullaway_undo.copy_state(ligand);
                    #endif

                    #if !_use_gloms
                    Molecule* mtmp[4], *mbkg[2];
                    mbkg[0] = ligand;
                    mbkg[1] = nullptr;
                    mtmp[0] = alignment_aa[0];
                    mtmp[1] = alignment_aa[1];
                    mtmp[2] = alignment_aa[2];
                    mtmp[3] = nullptr;
                    m.movability = MOV_FLEXONLY;
                    if (flex) Molecule::multimol_conform(mtmp);
                    m.movability = MOV_ALL;
                    #endif

                    #if _DBG_STEPBYSTEP
                    cout << endl;
                    #endif
                }
                
                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Aligned ligand to AA." << endl;
                cout << endl;
                #endif
            }

            // float driftamt = 1.0 / (iters/25+1);
            // cout << pose << ":" << nodeno << " drift " << driftamt << endl;
            int iters_div = iters*0.259;

            Molecule* cfmols[SPHREACH_MAX+4];
            for (i=0; i<SPHREACH_MAX+4; i++) cfmols[i] = nullptr;
            gcfmols = cfmols;
            i=0;
            m.movability = MOV_ALL;
            cfmols[i++] = &m;
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

            if (flex)
            {
                #if flexion_selection
                for (j=0; j<flexible_resnos.size(); j++)
                {
                    cfmols[i++] = protein->get_residue(flexible_resnos[j]);
                }
                #else
                for (j=0; j<sphres; j++)
                {
                    if (reaches_spheroid[nodeno][j]->movability >= MOV_FLEXONLY) reaches_spheroid[nodeno][j]->movability = MOV_FLEXONLY;
                    cfmols[i++] = reaches_spheroid[nodeno][j];
                }
                #endif
            }
            for (; i<SPHREACH_MAX; i++) cfmols[i] = NULL;

            ligand->reset_conformer_momenta();

            Molecule** delete_me;
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

            protein->find_residue_initial_bindings();
            Molecule::multimol_conform(
                cfmols,
                delete_me = protein->all_residues_as_molecules_except(cfmols),
                trip,
                iters,
                &iteration_callback
            );
            delete[] delete_me;

            /*time_t jlgsux = time(NULL);
            cout << "\nIterations took: " << (jlgsux-preiter) << " seconds." << endl;*/

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

            char metrics[protein->get_seq_length()+8][10];
            float mkJmol[protein->get_seq_length()+8];
            float imkJmol[protein->get_seq_length()+8];
            float mvdWrepl[protein->get_seq_length()+8];
            float imvdWrepl[protein->get_seq_length()+8];
            int metcount = 0;
            float btot = 0;
            float pstot = 0;

            for (i=0; i<_INTER_TYPES_LIMIT; i++) total_binding_by_type[i] = 0;

            if (debug) *debug << "Pose " << pose << " pathnode " << nodeno /*<< " clashes " << clash*/ << endl;

            m.clear_atom_binding_energies();

            if (met)
            {
                met->clear_atom_binding_energies();
                float lb = m.get_intermol_binding(met);
                strcpy(metrics[metcount], "Metals");
                mkJmol[metcount] = lb;
                imkJmol[metcount] = 0;								// TODO

                mvdWrepl[metcount] = 0;
                mvdWrepl[metcount] += m.get_vdW_repulsion(met);		// TODO: Include repulsions with non-mcoord side chains.

                imvdWrepl[metcount] = 0;							// TODO

                metcount++;

                btot += lb;
                // cout << "Metal adds " << lb << " to btot, making " << btot << endl;
            }

            float final_binding[seql+4];
            float final_vdWrepl[seql+4];
            for (i=0; i<seql+4; i++) final_binding[i] = final_vdWrepl[i] = 0;

            std::vector<AminoAcid*> allres = protein->get_residues_near(pocketcen, 10000);
            qpr = allres.size();
            Molecule* postaa[seql+8];
            postaa[0] = ligand;
            for (i=0; i<qpr; i++)
            {
                postaa[i+1] = reinterpret_cast<Molecule*>(allres[i]);
            }


            float tripclash = 0;


            if (differential_dock)
            {
                for (i=0; i<_INTER_TYPES_LIMIT; i++) total_binding_by_type[i] = 0;

                for (i=0; i<qpr+1; i++)
                {
                    int resno = i ? (allres[i-1]->get_residue_no()) : 0;
                    #if _DBG_TOOLARGE_DIFFNUMS
                    std::string ibdbg = to_string(resno) + (std::string)" ibdbg:\n";
                    #endif

                    bool is_trip_i = false;
                    for (n=0; n<tripswitch_clashables.size(); n++)
                        if (tripswitch_clashables[n] == resno)
                        {
                            is_trip_i = true;
                            break;
                        }

                    for (j=0; j<qpr+1; j++)
                    {
                        if (j == i) continue;
                        int jres = j ? allres[j-1]->get_residue_no() : 0;

                        bool is_trip_j = false;
                        if (is_trip_i && j)
                            for (n=0; n<tripswitch_clashables.size(); n++)
                                if (tripswitch_clashables[n] == jres)
                                {
                                    is_trip_j = true;
                                    break;
                                }

                        float f = postaa[i]->get_intermol_binding(postaa[j], j==0);
                        if (f < 0 && is_trip_j)
                        {
                            tripclash -= f;
                            f = 0;
                        }
                        final_binding[resno] += f;

                        #if _DBG_TOOLARGE_DIFFNUMS
                        if (f) ibdbg += to_string(postaa[j]->get_atom(0)->residue) + (std::string)" " + to_string(f) + (std::string)"\n";
                        #endif

                        final_vdWrepl[resno] += postaa[i]->get_vdW_repulsion(postaa[j]);
                    }

                    #if _DBG_TOOLARGE_DIFFNUMS
                    if (fabs(final_binding[resno]) >= 200) cout << ibdbg << endl;
                    #endif
                }
            }
            else
            {                
                for (i=0; i<qpr; i++)
                {
                    int resno = allres[i]->get_residue_no();

                    bool is_trip_i = false;
                    for (n=0; n<tripswitch_clashables.size(); n++)
                        if (tripswitch_clashables[n] == resno)
                        {
                            is_trip_i = true;
                            break;
                        }

                    for (j=0; j<qpr; j++)
                    {
                        if (j == i) continue;
                        int jres = allres[j]->get_residue_no();

                        bool is_trip_j = false;
                        if (is_trip_i)
                            for (n=0; n<tripswitch_clashables.size(); n++)
                                if (tripswitch_clashables[n] == jres)
                                {
                                    is_trip_j = true;
                                    break;
                                }
                        if (!is_trip_j) continue;

                        float f = postaa[i]->get_intermol_binding(postaa[j], j==0);
                        if (f < 0)
                        {
                            tripclash -= f;
                            f = 0;
                        }
                    }
                }
            }

            float fin_total_binding_by_type[_INTER_TYPES_LIMIT];
            for (i=0; i<_INTER_TYPES_LIMIT; i++) fin_total_binding_by_type[i] = total_binding_by_type[i];

            #if active_persistence
            for (i=0; i<=seql; i++) res_kJmol[i] = 0;
            #endif

            sphres = protein->get_residues_can_clash_ligand(reaches_spheroid[nodeno], &m, m.get_barycenter(), size, addl_resno);
            // cout << "sphres " << sphres << endl;
            float maxclash = 0;
            for (i=0; i<sphres; i++)
            {
                if (!reaches_spheroid[nodeno][i]) continue;
                if (!protein->aa_ptr_in_range(reaches_spheroid[nodeno][i])) continue;
                reaches_spheroid[nodeno][i]->clear_atom_binding_energies();
                int resno = reaches_spheroid[nodeno][i]->get_residue_no();

                float lb = m.get_intermol_binding(reaches_spheroid[nodeno][i], false);
                if (lb < -maxclash) maxclash -= lb;

                if (differential_dock)
                {
                    mkJmol[metcount] = final_binding[resno] + lb;
                }
                else
                {
                    if (lb > 90) lb = 0;
                    mkJmol[metcount] = lb;
                }

                #if active_persistence
                res_kJmol[resno] = lb;
                #endif

                sprintf(metrics[metcount], "%s%d", reaches_spheroid[nodeno][i]->get_3letter(), resno);
                // cout << metrics[metcount] << ": " << lb << " . ";

                if (differential_dock)
                {
                    imkJmol[metcount] = initial_binding[resno];
                    mvdWrepl[metcount] = final_vdWrepl[resno];
                    imvdWrepl[metcount] = initial_vdWrepl[resno];
                }
                else
                {
                    mvdWrepl[metcount] = 0;
                    mvdWrepl[metcount] += m.get_vdW_repulsion(reaches_spheroid[nodeno][i]);
                    /*for (j=0; j<sphres; j++)
                    {
                    	if (j == i) continue;
                    	mvdWrepl[metcount] += reaches_spheroid[nodeno][i]->get_vdW_repulsion(reaches_spheroid[nodeno][j]);
                    }*/
                    imkJmol[metcount] = 0;
                    imvdWrepl[metcount] = 0;
                }
                metcount++;
                btot += lb;
                // cout << *(reaches_spheroid[nodeno][i]) << " adds " << lb << " to btot, making " << btot << endl;

                float lf = m.get_intermol_polar_sat(reaches_spheroid[nodeno][i]);
                pstot += lf;

                #if _dbg_polsat
                cout << *(reaches_spheroid[nodeno][i]) << " adds " << lf << " to pstot, making " << pstot << endl;
                #endif
            }
            // cout << btot << endl;

            if (btot > 60*m.get_atom_count()) btot = 0;
            if (differential_dock && (maxclash > individual_clash_limit)) btot = -Avogadro;

            // drcount = pose-1+found_poses;

            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Prepared metrics." << endl;
            #endif

            // Set the dock result properties and allocate the arrays.
            dr[drcount][nodeno].kJmol = (differential_dock && (maxclash > individual_clash_limit)) ? -Avogadro : btot;
            dr[drcount][nodeno].ikJmol = 0;
            dr[drcount][nodeno].polsat  = pstot;
            dr[drcount][nodeno].metric   = new char*[metcount+4];
            dr[drcount][nodeno].mkJmol    = new float[metcount];
            dr[drcount][nodeno].imkJmol    = new float[metcount];
            dr[drcount][nodeno].mvdWrepl    = new float[metcount];
            dr[drcount][nodeno].imvdWrepl    = new float[metcount];
            dr[drcount][nodeno].tripswitch  = tripclash;
            dr[drcount][nodeno].proximity  = ligand->get_barycenter().get_3d_distance(nodecen);
            dr[drcount][nodeno].protclash = protein->get_rel_int_clashes();
            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Allocated memory." << endl;
            #endif

            int ahsz = active_helix_rots.size();
            if (l)
            {
                for (i=0; i<ahsz; i++)
                {
                    if (active_helix_rots[i].soft)
                    {
                        dr[drcount][nodeno].softrock += active_helix_rots[i].regname;
                        dr[drcount][nodeno].softrock += " active theta: ";
                        float deg = active_helix_rots[i].theta * fiftyseven * 100;
                        deg = round(deg);
                        deg /= 100;
                        dr[drcount][nodeno].softrock += std::to_string(deg);
                        dr[drcount][nodeno].softrock += " degrees\n";
                    }
                }
            }

            int itn;
            for (itn=0; itn<_INTER_TYPES_LIMIT; itn++)
            {
                i = itn;
                dr[drcount][nodeno].bytype[i] = differential_dock ? fin_total_binding_by_type[i] : total_binding_by_type[i];
                dr[drcount][nodeno].ibytype[i] = init_total_binding_by_type[i];
                dr[drcount][nodeno].ikJmol += init_total_binding_by_type[i];
                dr[drcount][nodeno].kJmol += fin_total_binding_by_type[i];
                // cout << drcount << "|" << i << " ";
            }
            // cout << endl;
            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Filled btypes." << endl;
            #endif

            // Populate the array.
            for (i=0; i<metcount; i++)
            {
                dr[drcount][nodeno].metric[i] = new char[max(8,(int)strlen(metrics[i])+4)];
                strcpy(dr[drcount][nodeno].metric[i], metrics[i]);
                dr[drcount][nodeno].mkJmol[i] = mkJmol[i];
                dr[drcount][nodeno].imkJmol[i] = imkJmol[i];
                dr[drcount][nodeno].mvdWrepl[i] = mvdWrepl[i];
                dr[drcount][nodeno].imvdWrepl[i] = imvdWrepl[i];
                // cout << "*" << dr[drcount][nodeno].metric[i] << ": " << dr[drcount][nodeno].mkJmol[i] << endl;
            }

            // Terminate with an empty string and a null pointer.
            dr[drcount][nodeno].metric[i] = new char[1];
            dr[drcount][nodeno].metric[i][0] = 0;
            dr[drcount][nodeno].metric[i+1] = 0;
            #if _DBG_STEPBYSTEP
            if (debug) *debug << "More metrics or something idfk." << endl;
            #endif

            std::ostringstream pdbdat;

            // Prepare a partial PDB of the ligand atoms and all involved residue sidechains.
            n = m.get_atom_count();
            int offset = n;
            for (l=0; l<n; l++) m.get_atom(l)->stream_pdb_line(pdbdat, 9000+l);
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

            dr[drcount][nodeno].pdbdat = pdbdat.str();
            // cout << "Attempt " << drcount << " node " << nodeno << " pdbdat is " << dr[drcount][nodeno].pdbdat.length() << " chars." << endl;
            if (debug) *debug << "Prepared the PDB strings." << endl;

            if (!nodeno)
            {
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

    // Output the dr[][] array in order of increasing pose number.
    cout << endl;
    if (output) *output << endl;

    const float energy_mult = kcal ? _kcal_per_kJ : 1;
    pose = 1;
    for (i=1; i<=poses; i++)
    {
        for (j=0; j<poses; j++)
        {
            if (dr[j][0].pose == i)
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

                        if (differential_dock)
                        {
                            cout << "# Binding energies: delta = with ligand minus without ligand." << endl;
                        }
                        else
                        {
                            cout << "# Binding energies:" << endl;
                        }
                        cout << "BENERG:" << endl;
                        if (output) *output << "# Binding energies" << endl << "BENERG:" << endl;
                        for (	l=0;

                                dr[j][k].metric
                                && dr[j][k].metric[l]
                                && dr[j][k].metric[l][0];

                                l++
                            )
                        {
                            if (differential_dock)
                            {
                                cout << dr[j][k].metric[l]
                                     << ": " << -(dr[j][k].mkJmol[l] - dr[j][k].imkJmol[l])*energy_mult
                                     << " = " << -dr[j][k].mkJmol[l]*energy_mult
                                     << " minus " << -dr[j][k].imkJmol[l]*energy_mult
                                     << endl;
                                if (output && dr[j][k].metric[l]) *output << dr[j][k].metric[l]
                                            << ": " << -(dr[j][k].mkJmol[l] - dr[j][k].imkJmol[l])*energy_mult
                                            << " = " << -dr[j][k].mkJmol[l]*energy_mult
                                            << " minus " << -dr[j][k].imkJmol[l]*energy_mult
                                            << endl;
                            }
                            else
                            {
                                if (output && do_output_colors) colorize(dr[j][k].mkJmol[l]);
                                cout << dr[j][k].metric[l] << ": " << -dr[j][k].mkJmol[l]*energy_mult << endl;
                                if (do_output_colors) colorless();
                                if (output && dr[j][k].metric[l]) *output << dr[j][k].metric[l] << ": " << -dr[j][k].mkJmol[l]*energy_mult << endl;
                            }
                        }
                        cout << endl;
                        if (output) *output << endl;

                        for (l=0; l<_INTER_TYPES_LIMIT; l++)
                        {
                            char lbtyp[64];
                            switch (l+covalent)
                            {
                            case covalent:
                                continue; /*strcpy(lbtyp, "Total covalent: ");		break;*/
                            case ionic:
                                strcpy(lbtyp, "Total ionic: ");
                                break;
                            case hbond:
                                strcpy(lbtyp, "Total H-bond: ");
                                break;
                            case pi:
                                strcpy(lbtyp, "Total pi stack: ");
                                break;
                            case polarpi:
                                strcpy(lbtyp, "Total polar-pi and cation-pi: ");
                                break;
                            case mcoord:
                                strcpy(lbtyp, "Total metal coordination: ");
                                break;
                            case vdW:
                                strcpy(lbtyp, "Total van der Waals: ");
                                break;
                            default:
                                goto _btyp_unassigned;
                            }

                            if (differential_dock)
                            {
                                cout << lbtyp << -(dr[j][k].bytype[l] - dr[j][k].ibytype[l])*energy_mult
                                     << " = " << -dr[j][k].bytype[l]*energy_mult
                                     << " minus " << -dr[j][k].ibytype[l]*energy_mult
                                     << endl;
                                if (output) *output << lbtyp << -(dr[j][k].bytype[l] - dr[j][k].ibytype[l])*energy_mult
                                                        << " = " << -dr[j][k].bytype[l]*energy_mult
                                                        << " minus " << -dr[j][k].ibytype[l]*energy_mult
                                                        << endl;
                            }
                            else
                            {
                                if (output)
                                {
                                    *output << lbtyp << -dr[j][k].bytype[l]*energy_mult << endl;
                                    if (do_output_colors) colorize(dr[j][k].bytype[l]);
                                    cout << lbtyp << -dr[j][k].bytype[l]*energy_mult << endl;
                                    if (do_output_colors) colorless();
                                }
                                else
                                {
                                    cout << lbtyp << -dr[j][k].bytype[l]*energy_mult << endl;
                                }
                            }
                        }
                        cout << endl;
                        if (output) *output << endl;

                    _btyp_unassigned:

                        if (differential_dock)
                        {
                            if (output) *output << "Total: " << -(dr[j][k].kJmol - dr[j][k].ikJmol)*energy_mult
                                                    << " = " << -dr[j][k].kJmol*energy_mult
                                                    << " minus " << -dr[j][k].ikJmol*energy_mult
                                                    << endl << endl;
                            cout << "Total: " << -(dr[j][k].kJmol - dr[j][k].ikJmol)*energy_mult
                                 << " = " << -dr[j][k].kJmol*energy_mult
                                 << " minus " << -dr[j][k].ikJmol*energy_mult
                                 << endl << endl;
                        }
                        else
                        {
                            if (output) *output << "Total: " << -dr[j][k].kJmol*energy_mult << endl << endl;
                            if (output && do_output_colors) colorize(dr[j][k].kJmol);
                            cout << "Total: " << -dr[j][k].kJmol*energy_mult << endl << endl;
                            if (output && do_output_colors) colorless();
                        }

                        if (dr[j][k].softrock.size())
                        {
                            cout << "Active Helix Soft Rotations:" << endl
                                 << dr[j][k].softrock << endl;
                            if (output ) *output << "Active Helix Soft Rotations:" << endl
                                 << dr[j][k].softrock << endl;
                        }

                        cout << "Ligand polar satisfaction: " << dr[j][k].polsat << endl;
                        if (output) *output << "Ligand polar satisfaction: " << dr[j][k].polsat << endl;
                        cout << endl;
                        if (output) *output << endl;

                        if (output) *output << "Proximity: " << dr[j][k].proximity << endl << endl;
                        cout << "Proximity: " << dr[j][k].proximity << endl << endl;

                        if (output) *output << "Protein clashes: " << dr[j][k].protclash << endl << endl;
                        cout << "Protein clashes: " << dr[j][k].protclash << endl << endl;

                        if (tripswitch_clashables.size())
                        {
                            if (output) *output << "Trip switch: " << dr[j][k].tripswitch << endl << endl;
                            cout << "Trip switch: " << dr[j][k].tripswitch << endl << endl;
                        }

                        if (differential_dock)
                        {
                            cout << "# van der Waals repulsion: delta = with ligand minus without ligand." << endl;
                        }
                        else
                        {
                            cout << "# van der Waals repulsion:" << endl;
                        }
                        cout << "vdWRPL:" << endl;
                        if (output) *output << "# van der Waals repulsion" << endl << "vdWRPL:" << endl;
                        for (	l=0;

                                dr[j][k].metric
                                && dr[j][k].metric[l]
                                && dr[j][k].metric[l][0];

                                l++
                            )
                        {
                            if (fabs(dr[j][k].mvdWrepl[l]) < 0.001) continue;

                            if (differential_dock)
                            {
                                cout << dr[j][k].metric[l]
                                     << ": " << (dr[j][k].mvdWrepl[l] - dr[j][k].imvdWrepl[l])*energy_mult
                                     << " = " << dr[j][k].mvdWrepl[l]*energy_mult
                                     << " minus " << dr[j][k].imvdWrepl[l]*energy_mult
                                     << endl;
                                if (output && dr[j][k].metric[l]) *output << dr[j][k].metric[l]
                                            << ": " << (dr[j][k].mvdWrepl[l] - dr[j][k].imvdWrepl[l])*energy_mult
                                            << " = " << dr[j][k].mvdWrepl[l]*energy_mult
                                            << " minus " << dr[j][k].imvdWrepl[l]*energy_mult
                                            << endl;
                            }
                            else
                            {
                                cout << dr[j][k].metric[l] << ": " << dr[j][k].mvdWrepl[l]*energy_mult << endl;
                                if (output && dr[j][k].metric[l]) *output << dr[j][k].metric[l] << ": " << dr[j][k].mvdWrepl[l]*energy_mult << endl;
                            }
                        }
                        cout << endl;
                        if (output) *output << endl;

                        if (soft_pocket)
                        {
                            cout << "Soft transformations:" << endl;
                            if (output) *output << "Soft transformations:" << endl;
                            for (l=0; l<soft_rgns.size(); l++)
                            {
                                cout << soft_rgns[l].name << ".Δr: " << rgnxform_r[j][k][l] << endl;
                                cout << soft_rgns[l].name << ".Δθ: " << rgnxform_theta[j][k][l] << endl;
                                cout << soft_rgns[l].name << ".Δy: " << rgnxform_y[j][k][l] << endl;
                                cout << soft_rgns[l].name << ".Δα: " << rgnrot_alpha[j][k][l]*fiftyseven << endl;
                                cout << soft_rgns[l].name << ".Δφw: " << rgnrot_w[j][k][l]*fiftyseven << endl;
                                cout << soft_rgns[l].name << ".Δφu: " << rgnrot_u[j][k][l]*fiftyseven << endl;
                                if (output) *output << soft_rgns[l].name << ".Δr: " << rgnxform_r[j][k][l] << endl;
                                if (output) *output << soft_rgns[l].name << ".Δθ: " << rgnxform_theta[j][k][l] << endl;
                                if (output) *output << soft_rgns[l].name << ".Δy: " << rgnxform_y[j][k][l] << endl;
                                if (output) *output << soft_rgns[l].name << ".Δα: " << rgnrot_alpha[j][k][l]*fiftyseven << endl;
                                if (output) *output << soft_rgns[l].name << ".Δφw: " << rgnrot_w[j][k][l]*fiftyseven << endl;
                                if (output) *output << soft_rgns[l].name << ".Δφu: " << rgnrot_u[j][k][l]*fiftyseven << endl;
                            }
                            cout << endl;
                            if (output) *output << endl;
                        }

                        /*if (flex)
                        {*/
                            if (!dr[j][k].pdbdat.length())
                            {
                                cout << "WARNING: Failed to generate PDB data." << endl;
                                if (output) *output << "(Missing PDB data.)" << endl;
                            }
                            else
                            {
                                if (!output || echo_pdb_data) cout << "# PDB Data" << endl << "PDBDAT:" << endl;
                                if (output) *output << "# PDB Data" << endl << "PDBDAT:" << endl;

                                if (output) *output << dr[j][k].pdbdat << endl;
                                if (!output || echo_pdb_data) cout << dr[j][k].pdbdat << endl;

                                if (!output || echo_pdb_data) cout << "TER" << endl << "END" << endl << endl << endl;
                                if (output) *output << "TER" << endl << "END" << endl << endl << endl;
                            }
                        // }

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
            hydrogenate_pdb = false;
            pf = fopen(tplset ? temp_pdb_file.c_str() : protfname, "r");
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
            cout << (tplset ? "Homologized " : "Original ") << "PDB appended to output file." << endl;
        }
        else cout << "ERROR: Append PDB can only be used when specifying an output file." << endl;
    }

    if (debug) debug->close();

    return 0;
}


















