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
    std::string pdbdat;
    float bytype[_INTER_TYPES_LIMIT];
    float ibytype[_INTER_TYPES_LIMIT];
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

char configfname[256];
char protfname[256];
char ligfname[256];
Point pocketcen, loneliest, pocketsize, ligbbox;
std::ofstream *output = NULL;

std::vector<int> exclusion;

std::string CEN_buf = "";
std::vector<std::string> pathstrs;
std::vector<std::string> states;
std::vector<int> extra_wt;

bool configset=false, protset=false, ligset=false, pktset=false;

Protein* protein;
int seql = 0;
int mcoord_resno[256];
Molecule* ligand;
Point ligcen_target;
Point size(10,10,10);
SCoord path[256];
int pathnodes=1;		// The pocketcen is the initial node.
int poses = 10;
int iters = 50;
bool flex=true;
float kJmol_cutoff = 0.01;
bool kcal = false;
float drift = 0.333;
Molecule** gcfmols = NULL;

std::string origbuff = "";

// Switch to enable older "best binding" algorithm instead of newer "tumble spheres".
// Generally, tumble spheres give better results but let's leave best bind code
// in place because we'll be reviving it later.
bool use_bestbind_algorithm = default_bestbind;

void iteration_callback(int iter)
{
    #if !allow_iter_cb
    return;
    #endif

    // if (kJmol_cutoff > 0 && ligand->lastbind >= kJmol_cutoff) iter = (iters-1);
    if (iter == (iters-1)) return;

    Point bary = ligand->get_barycenter();

    #if allow_drift
    #if !pocketcen_is_loneliest
    if (ligand->lastbind > 100)
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
        bary.x += (ligcen_target.x - bary.x) * drift;
        bary.y += (ligcen_target.y - bary.y) * drift;
        bary.z += (ligcen_target.z - bary.z) * drift;
    }

    ligand->recenter(bary);

    drift *= (1.0 - 0.5/iters);
    #endif

    if (gcfmols && seql)
    {
        Star discrete[SPHREACH_MAX+4];
        discrete[0].pmol = gcfmols[0];
        discrete[1].pmol = gcfmols[1];

        int i;
        AminoAcid* resphres[SPHREACH_MAX+4];
        for (i=0; i<SPHREACH_MAX+4; i++) resphres[i] = nullptr;
        int sphres = protein->get_residues_can_clash_ligand(resphres, ligand, bary, size, mcoord_resno);
        //cout << "Sphres: " << sphres << endl;
        for (i=0; i<sphres; i++)
        {
            discrete[i+2].paa = resphres[i];
        }
        discrete[sphres+2].n = 0;

        sphres += 2;
        for (i=0; i<sphres; i++) gcfmols[i] = discrete[i].pmol;
        gcfmols[sphres] = nullptr;
    }
}

Point pocketcen_from_config_fields(char** fields, Point* old_pocketcen)
{
    int i=1;
    Point local_pocketcen;
    if (!strcmp(fields[i], "RES"))
    {
        i++;
        std::vector<int> resnos;
        for (; fields[i]; i++)
        {
            int j = atoi(fields[i]);
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
    else if (!strcmp(fields[i], "REL"))
    {
        if (!old_pocketcen)
        {
            cout << "Error: relative coordinates not supported for CEN." << endl;
            throw 0xbadb19d;
        }
        else
        {
            i++;
            local_pocketcen.x = old_pocketcen->x + atof(fields[i++]);
            local_pocketcen.y = old_pocketcen->y + atof(fields[i++]);
            local_pocketcen.z = old_pocketcen->z + atof(fields[i++]);
            return local_pocketcen;
        }
    }
    else
    {
        if (!strcmp(fields[i], "ABS")) i++;
        local_pocketcen.x = atof(fields[i++]);
        local_pocketcen.y = atof(fields[i++]);
        local_pocketcen.z = atof(fields[i++]);
        return local_pocketcen;
    }
}

int interpret_config_line(char** fields)
{
	int i;
	
	if (!strcmp(fields[0], "PROT"))
    {
        strcpy(protfname, fields[1]);
        protset = true;
        return 1;
    }
    else if (!strcmp(fields[0], "LIG"))
    {
        strcpy(ligfname, fields[1]);
        ligset = true;
        return 1;
    }
    else if (!strcmp(fields[0], "CEN"))
    {
        CEN_buf = origbuff;
        return 0;
    }
    else if (!strcmp(fields[0], "EXCL"))
    {
        i=1;
        int excls = atoi(fields[i++]);
        int excle = atoi(fields[i++]);

        for (i=excls; i<=excle; i++) exclusion.push_back(i);
        return i-1;
    }
    else if (!strcmp(fields[0], "PATH"))
    {
        i=1;
        if (!strcmp(fields[i], "ABS")) i++;
        if (!strcmp(fields[i], "REL")) i++;
        if (!strcmp(fields[i], "RES")) i++;

        int nodeno = atoi(fields[i]);
        if (nodeno > 255)
        {
            cout << "Binding path is limited to 255 nodes." << endl;
            throw 0xbad90de;
        }
        if (nodeno)
        {
            if ((nodeno) > pathnodes) pathnodes = nodeno;

            pathstrs.resize(nodeno+1);
            pathstrs[nodeno] = origbuff;
        }
        return i-1;
    }
    else if (!strcmp(fields[0], "STATE"))
    {
        states.push_back(origbuff);
        return 0;
    }
    else if (!strcmp(fields[0], "SIZE"))
    {
        size.x = atof(fields[1]);
        if (fields[2])
        {
            size.y = atof(fields[2]);
            size.z = atof(fields[3]);
        }
        else size.z = size.y = size.x;
        if (!size.x || !size.y || !size.z)
        {
            cout << "Pocket size cannot be zero in any dimension." << endl;
            throw 0xbad512e;
        }
        return 3;
    }
    else if (!strcmp(fields[0], "POSE"))
    {
        poses = atoi(fields[1]);
        return 1;
    }
    else if (!strcmp(fields[0], "EMIN"))
    {
        kJmol_cutoff = atof(fields[1]);
        return 1;
    }
    else if (!strcmp(fields[0], "ELIM"))
    {
        kJmol_cutoff = -atof(fields[1]);
        return 1;
    }
    else if (!strcmp(fields[0], "DIFF"))
    {
        differential_dock = true;
    }
    else if (!strcmp(fields[0], "FLEX"))
    {
        flex = (atoi(fields[1]) != 0);
        return 1;
    }
    else if (!strcmp(fields[0], "KCAL"))
    {
        kcal = true;
        return 0;
    }
    else if (!strcmp(fields[0], "RLIM"))
    {
        _INTERA_R_CUTOFF = atof(fields[1]);
        return 1;
    }
    else if (!strcmp(fields[0], "ITER") || !strcmp(fields[0], "ITERS"))
    {
        iters = atoi(fields[1]);
        return 1;
    }
    else if (!strcmp(fields[0], "MCOORD"))
    {
        int j=0;
        for (i=1; fields[i]; i++)
        {
        	if (fields[i][0] == '-' && fields[i][1] == '-') break;
            mcoord_resno[j++] = atoi(fields[i]);
        }
        mcoord_resno[j] = 0;
        return i-1;
    }
    else if (!strcmp(fields[0], "DEBUG"))
    {
        if (!fields[1])
        {
            cout << "Missing debug file name; check config file." << endl;
            throw 0xbadf12e;
        }
        #if _DBG_STEPBYSTEP
        cout << "Starting a debug outstream." << endl;
        #endif
        debug = new std::ofstream(fields[1], std::ofstream::out);
        return 1;
    }
    else if (!strcmp(fields[0], "OUT"))
    {
        if (!fields[1])
        {
            cout << "Missing output file name; check config file." << endl;
            throw 0xbadf12e;
        }
        #if _DBG_STEPBYSTEP
        cout << "Starting a file outstream." << endl;
        #endif
        output = new std::ofstream(fields[1], std::ofstream::out);
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
        fgets(buffer, 1015, pf);
        origbuff = buffer;
        if (buffer[0] >= ' ' && buffer[0] != '#')
        {
            char** fields = chop_spaced_fields(buffer);
            if (!fields) continue;

            interpret_config_line(fields);
        }
        buffer[0] = 0;
    }
}

int main(int argc, char** argv)
{
    char buffer[65536];
    int i, j;

    for (i=0; i<65536; i++) buffer[i] = 0;

    for (i=0; i<256; i++)
        configfname[i] = protfname[i] = ligfname[i] = 0;

    time_t began = time(NULL);

    strcpy(configfname, "primarydock.config");
    
    for (i=1; i<argc; i++)
    {
    	if (argv[i][0] == '-' && argv[i][1] == '-')
    	{
    		argv[i] += 2;
    		for (j=0; argv[i][j]; j++) if (argv[i][j] >= 'a' && argv[i][j] <= 'z') argv[i][j] &= 0x5f;
    		j = interpret_config_line(&argv[i]);
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
    		argv[i] -= 2;
    		i += j;
    	}
    }

    pre_ligand_flex_radius = size.magnitude();
    pre_ligand_multimol_radius = pre_ligand_flex_radius + (default_pre_ligand_multimol_radius - default_pre_ligand_flex_radius);

    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Loaded config file." << endl;
    #endif

    if (kcal) kJmol_cutoff /= _kcal_per_kJ;
    drift = 1.0 / (iters/25+1);

    // Load the protein or return an error.
    Protein p(protfname);
    protein = &p;
    pf = fopen(protfname, "r");
    if (!pf)
    {
        cout << "Error trying to read " << protfname << endl;
        return 0xbadf12e;
    }
    p.load_pdb(pf);
    fclose(pf);
    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Loaded protein." << endl;
    #endif

    if (!CEN_buf.length())
    {
        cout << "Error: no binding pocket center defined." << endl;
        return 0xbadb19d;
    }

    strcpy(buffer, CEN_buf.c_str());

    char** fields = chop_spaced_fields(buffer);
    pocketcen = pocketcen_from_config_fields(fields, nullptr);
    loneliest = p.find_loneliest_point(pocketcen, size);

    #if pocketcen_is_loneliest
    pocketcen = loneliest;
    #endif

    if (!strcmp(fields[1], "RES"))
    {
        for (i=2; fields[i]; i++)
        {
            extra_wt.push_back(atoi(fields[i]));
        }
    }
    pktset = true;

    p.mcoord_resnos = mcoord_resno;

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
        fread(buffer, 1, 65535, pf);
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
    int k, l, n;
    Atom** ligbb = m.get_most_bindable(3);
    Atom** ligbbh = new Atom*[5];
    intera_type lig_inter_typ[5];

    if (use_bestbind_algorithm)
    {
        for (i=0; i<5; i++)
        {
            ligbbh[i] = NULL;
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
            else if (fabs(ligbb[i]->is_polar()) >= 1
                     ||
                     (ligbbh[i] && fabs(ligbbh[i]->is_polar()) >= 1)
                    )
                lig_inter_typ[i] = hbond;
            else if (ligbb[i]->is_pi())
                lig_inter_typ[i] = pi;
        }
    }

    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Identified best binding ligand atoms." << endl;
    #endif

    int pose, nodeno, iter;
    Point nodecen = pocketcen;
    seql = p.get_seq_length();
    int rstart = p.get_start_resno();
    AminoAcid* reaches_spheroid[pathnodes+2][SPHREACH_MAX];
    int sphres = 0;

    // Filter residues according to which ones are close enough to the spheroid to "reach" it.
    nodecen = pocketcen;

    // When docking with a metalloprotein, use this temporary Molecule for interactions the same as
    // we use AminoAcid objects, except don't attempt to flex the metals object.
    Molecule* met = p.metals_as_molecule();
    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Created metals molecule." << endl;
    #endif

    float bclash = 0;

    for (l=0; l<seql; l++)
    {
        AminoAcid* laa = p.get_residue(l+rstart);
        if (!laa) continue;
        for (n=l+1; n<seql; n++)
        {
            AminoAcid* naa = p.get_residue(n+rstart);
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
            if (ligbb[i])
            {
                cout << "# Best binding heavy atom " << i << " of ligand" << endl << "# LBBA: " << ligbb[i]->name << endl;
                if (output) *output << "# Best binding heavy atom " << i << " of ligand" << endl << "LBBA: " << ligbb[i]->name << endl;
            }
            if (ligbbh[i])
            {
                cout << "# Best binding hydrogen " << i << " of ligand" << endl << "# LBBH: " << ligbbh[i]->name << endl;
                if (output) *output << "# Best binding hydrogen " << i << " of ligand" << endl << "LBBH: " << ligbbh[i]->name << endl;
            }
            cout << endl;
        }

    DockResult dr[poses+2][pathnodes+2];
    for (i=0; i<poses; i++) dr[i][0].kJmol = 0;
    int drcount = 0, qpr;

    srand(0xb00d1cca);
    // srand(time(NULL));
    float initial_binding[seql+4];
    float initial_vdWrepl[seql+4];
    float init_total_binding_by_type[_INTER_TYPES_LIMIT];
    for (pose=1; pose<=poses; pose++)
    {
        ligand->minimize_internal_clashes();
        float lig_min_int_clsh = m.get_internal_clashes();
        ligand->crumple(fiftyseventh*44);

        if (pose > 1)
        {
            // TODO: Revert to saved original locations for the side chain atoms instead of reloading the protein.
            pf = fopen(protfname, "r");
            p.load_pdb(pf);
            fclose(pf);
        }

        if (differential_dock)
        {
            std::vector<AminoAcid*> preres = p.get_residues_near(pocketcen, pre_ligand_multimol_radius);
            qpr = preres.size();
            AminoAcid* preaa[seql+4];
            MovabilityType aamov[seql+4];

            for (i=0; i<seql+4; i++) preaa[i] = nullptr;
            for (i=0; i<qpr; i++)
            {
                preaa[i] = preres[i];
                Atom* CA = preaa[i]->get_atom("CA");
                if (!CA)
                {
                    cout << "Residue " << preaa[i]->get_residue_no() << " is missing its CA." << endl;
                    throw 0xbad12e5;
                }
                float r = CA->get_location().get_3d_distance(pocketcen);
                
                aamov[i] = preaa[i]->movability;
                
                if (r > pre_ligand_flex_radius) preaa[i]->movability = MOV_NONE;
            }

            for (i=0; i<seql+4; i++) initial_binding[i] = initial_vdWrepl[i] = 0;

            if (pre_ligand_iteration_ratio)
            {
                Molecule** delete_me;
                Molecule::multimol_conform(reinterpret_cast<Molecule**>(preaa), delete_me = p.all_residues_as_molecules(), iters*pre_ligand_iteration_ratio);
                delete[] delete_me;
            }

            preres = p.get_residues_near(pocketcen, 10000);
            qpr = preres.size();
            for (i=0; i<qpr; i++)
            {
                preaa[i] = preres[i];
            }

            for (i=0; i<_INTER_TYPES_LIMIT; i++) total_binding_by_type[i] = 0;

            for (i=0; i<qpr; i++)
            {
                int resno = preaa[i]->get_residue_no();
                // std::string ibdbg = to_string(resno) + (std::string)" ibdbg:\n";
                for (j=0; j<qpr; j++)
                {
                    if (j == i) continue;
                    float f = reinterpret_cast<Molecule*>(preaa[i])->get_intermol_binding(reinterpret_cast<Molecule*>(preaa[j]), j==0);
                    // if (f) ibdbg += to_string(preaa[j]->get_residue_no()) + (std::string)" " + to_string(f) + (std::string)"\n";
                    initial_binding[resno] += f;
                    initial_vdWrepl[resno] += preaa[i]->get_vdW_repulsion(preaa[j]);
                }
                // if (fabs(initial_binding[resno]) >= 50) cout << ibdbg << endl;
            }

            for (i=0; i<_INTER_TYPES_LIMIT; i++) init_total_binding_by_type[i] = total_binding_by_type[i];

            for (i=0; i<qpr; i++) preaa[i]->movability = aamov[i];
        }

        if (pose <= 1)
        {
            if (!use_bestbind_algorithm)
            {
                // Begin tumble sphere behavior.
                std::vector<AminoAcid*> tsphres = p.get_residues_near(pocketcen, size.magnitude()+4);
                int tsphsz = tsphres.size();
                float outer_sphere[tsphsz+4], inner_sphere[tsphsz+4];

                pocketsize = p.estimate_pocket_size(tsphres);
                ligbbox = m.get_bounding_box();

                for (i=0; !ligbbox.fits_inside(pocketsize) && i<100; i++)
                {
                    m.crumple(fiftyseventh*30);
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
                const int ac = m.get_atom_count();
                Pose besp(&m);
                #if _DBG_TUMBLE_SPHERES
                std::string tsdbg = "", tsdbgb = "";
                #endif

                if (ligbbox.x > ligbbox.y && pocketsize.x < pocketsize.y) m.rotate(zaxis, square);
                if (ligbbox.x > ligbbox.z && pocketsize.x < pocketsize.z) m.rotate(yaxis, square);
                if (ligbbox.y > ligbbox.x && pocketsize.y < pocketsize.x) m.rotate(zaxis, square);
                if (ligbbox.y > ligbbox.z && pocketsize.y < pocketsize.z) m.rotate(xaxis, square);
                if (ligbbox.z > ligbbox.x && pocketsize.z < pocketsize.x) m.rotate(yaxis, square);
                if (ligbbox.z > ligbbox.y && pocketsize.z < pocketsize.y) m.rotate(xaxis, square);

                step = fiftyseventh*30;
                bestscore = -Avogadro;
                float lonely_step = 1.0 / loneliest.get_3d_distance(pocketcen);
                #if _DBG_LONELINESS
                cout << "Loneliest point " << loneliest << " is " << loneliest.get_3d_distance(pocketcen) << "A from pocketcen " << pocketcen << "." << endl;
                cout << "Pocket size is " << pocketsize << " vs. ligand bounding box " << ligbbox << endl;
                #endif
                if (isnan(lonely_step) || lonely_step < 0.1) lonely_step = 0.1;

                #if pocketcen_is_loneliest
                if (1)
                {
                    m.recenter(pocketcen);
                #else
                for (loneliness=0; loneliness <= 1; loneliness += lonely_step)
                {
                    float centeredness = 1.0 - loneliness;
                    Point tmpcen(loneliest.x * loneliness + pocketcen.x * centeredness,
                                 loneliest.y * loneliness + pocketcen.y * centeredness,
                                 loneliest.z * loneliness + pocketcen.z * centeredness
                                );
                    m.recenter(tmpcen);
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
                                ligbbox = m.get_bounding_box();

                                if (ligbbox.x > ligbbox.y && pocketsize.x < pocketsize.y) continue;
                                if (ligbbox.x > ligbbox.z && pocketsize.x < pocketsize.z) continue;
                                if (ligbbox.y > ligbbox.x && pocketsize.y < pocketsize.x) continue;
                                if (ligbbox.y > ligbbox.z && pocketsize.y < pocketsize.z) continue;
                                if (ligbbox.z > ligbbox.x && pocketsize.z < pocketsize.x) continue;
                                if (ligbbox.z > ligbbox.y && pocketsize.z < pocketsize.y) continue;

                                Bond** rb = m.get_rotatable_bonds();

                                if (!rb) n = 0;
                                else for (n=0; rb[n]; n++);		// Get count.

                                l = 0;
                                lrad = 0;
                            _xyzl_loop:
                                if (m.get_internal_clashes() >= lig_min_int_clsh*5+5) goto _xyzl_skip_loop;

                                score = 0;
                                #if _DBG_TUMBLE_SPHERES
                                tsdbg = "";
                                // cout << m.get_internal_clashes() << " vs. " << lig_min_int_clsh << endl;
                                #endif
                                for (i=0; i<ac; i++)
                                {
                                    Atom* a = m.get_atom(i);
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
                                    besp.copy_state(&m);
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

                                    char** lfields = chop_spaced_fields(protfttl, '/');

                                    for (u=0; lfields[u]; u++);
                                    u--;

                                    char fname[1000];
                                    sprintf(fname, "output/tumble_%s_%d_%d_%d_%f.dock",
                                            lfields[u],
                                            (int)(xrad*fiftyseven),
                                            (int)(yrad*fiftyseven),
                                            (int)(zrad*fiftyseven),
                                            score);
                                    cout << fname << endl;
                                    std::ofstream tspdbdat(fname, std::ofstream::out);

                                    tspdbdat << "PDB file: " << protfname << endl;
                                    tspdbdat << "Pose: 1\nNode: 0\nPDBDAT:\n";

                                    int lac = m.get_atom_count();
                                    for (u=0; u<lac; u++) m.get_atom(u)->stream_pdb_line(tspdbdat, 9000+u);

                                    int pseql = p.get_seq_length();
                                    v = 1;
                                    for (u = 1; u < pseql; u++)
                                    {
                                        AminoAcid* dbgaa = p.get_residue(u);
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

                                m.rotate(zaxis, step);
                            }		// zrad.
                            m.rotate(yaxis, step);
                        }			// yrad.
                        m.rotate(xaxis, step);
                    }				// xrad.

                    #if !pocketcen_is_loneliest
                    if (bestscore >= (m.get_atom_count()*13)) break;
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
                besp.restore_state(&m);

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

                #if debug_stop_after_tumble_sphere
                return 0;
                #endif
                // End tumble sphere behavior.
            }
        }

        #if _DBG_STEPBYSTEP
        if (debug) *debug << "Pose " << pose << endl;
        #endif
        nodecen = pocketcen;
        nodecen.weight = 1;

        for (nodeno=0; nodeno<=pathnodes; nodeno++)
        {
            if (pathstrs.size() < nodeno) break;

            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Pose " << pose << endl << "Node " << nodeno << endl;
            #endif
            if (nodeno)
            {
                for (i=0; i<states.size(); i++)
                {
                    strcpy(buffer, states[i].c_str());
                    fields = chop_spaced_fields(buffer);
                    if (atoi(fields[1]) == nodeno)
                    {
                        int sr = atoi(fields[2]), er = atoi(fields[3]);
                        float theta = atof(fields[4]) * fiftyseventh;

                        Point sloc = p.get_atom_location(sr, "CA"),
                              eloc = p.get_atom_location(er, "CA");

                        LocatedVector lv = (SCoord)(sloc.subtract(eloc));
                        lv.origin = sloc;

                        int resno;
                        for (resno = sr; resno <= er; resno++)
                        {
                            AminoAcid* aa = p.get_residue(resno);
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
                fields = chop_spaced_fields(buffer);
                nodecen = pocketcen_from_config_fields(fields, &nodecen);
                if (!strcmp(fields[1], "RES"))
                {
                    extra_wt.clear();
                    for (i=2; fields[i]; i++)
                    {
                        extra_wt.push_back(atoi(fields[i]));
                    }
                }

                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Added whatever points together." << endl;
                #endif
            }
            Point lastnodecen = nodecen;
            ligcen_target = nodecen;

            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Saved last nodecen." << endl;
            #endif

            // Move the ligand to the new node center.
            // m.recenter(nodecen);
            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Molecule recenter (or not)." << endl;
            #endif
            m.reset_conformer_momenta();
            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Conformer momenta reset." << endl;
            #endif


            sphres = p.get_residues_can_clash_ligand(reaches_spheroid[nodeno], &m, nodecen, size, mcoord_resno);
            for (i=sphres; i<SPHREACH_MAX; i++) reaches_spheroid[nodeno][i] = NULL;

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
                Molecule* alignment_aa[5];
                float alignment_distance[5];
                for (l=0; l<3; l++)
                {
                    alignment_aa[l]=0;
                    alignment_distance[l]=0;
                }

                for (l=0; l<=pathnodes; l++)
                {
                    dr[drcount][l].metric = 0;
                }

                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Initialize null AA pointer." << endl;
                #endif

                // Find a binding pocket feature with a strong potential binding to the ligand.
                std::string alignment_name = "";
                if (use_bestbind_algorithm) for (l=0; l<3; l++)
                {
                    if (!ligbb[l]) continue;
                    float alignment_potential = 0;
                    for (i=0; reaches_spheroid[nodeno][i]; i++)
                    {
                        if (!p.aa_ptr_in_range(reaches_spheroid[nodeno][i]))
                        {
                            reaches_spheroid[nodeno][i] = NULL;
                            continue;
                        }

                        if (l && reaches_spheroid[nodeno][i] == alignment_aa[l-1]) continue;
                        if (l>1 && reaches_spheroid[nodeno][i] == alignment_aa[l-2]) continue;

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
                        if (reaches_spheroid[nodeno][i]->capable_of_inter(lig_inter_typ[l])
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
                        }
                        #if _DBG_STEPBYSTEP
                        if (debug) *debug << "Candidate alignment AA." << endl;
                        #endif
                    }
                }
                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Selected an alignment AA." << endl;
                #endif

                if (use_bestbind_algorithm && met)
                {
                    alignment_aa[0] = met;
                    // alignment_name = "metal";
                }
                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Alignment AA." << endl;
                #endif

                if (use_bestbind_algorithm)	for (l=0; l<3; l++)
                {
                    if (alignment_aa[l])
                    {
                        cout << "Aligning " << ligbb[l]->name << " to " << alignment_aa[l]->get_name() << endl;
                        Atom* alca;
                        if (alignment_aa[l] == met)
                            alca = alignment_aa[l]->get_nearest_atom(ligbb[l]->get_location());
                        else
                        {
                            // alca = alignment_aa->get_atom("CA");
                            alca = alignment_aa[l]->get_most_bindable(1)[0];
                        }
                        #if _DBG_STEPBYSTEP
                        if (debug) *debug << "Got alignment atom." << endl;
                        #endif

                        if (alca)
                        {
                            Point pt, al, cen;
                            pt	= ligbb[l]->get_location();
                            if (ligbbh[l])
                            {
                                Point pth = ligbbh[l]->get_location();
                                pt.x += 0.5*(pth.x-pt.x);
                                pt.y += 0.5*(pth.y-pt.y);
                                pt.z += 0.5*(pth.z-pt.z);
                            }
                            al	= alca->get_location();
                            cen	= (l==1) ? ligbb[0]->get_location() : m.get_barycenter();

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
                                rot = align_points_3d(&pt, &al, &cen);
                                m.rotate(&rot.v, rot.a);
                                cout << "# Pivoted ligand " << (rot.a*fiftyseven) << "deg about ligand molcen." << endl;

                                if (!l)
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
                                rot = align_points_3d(&pt, &al, &origin);
                                lv.copy(rot.v);
                                lv.origin = origin;
                                m.rotate(lv, rot.a);
                                cout << "# Pivoted ligand " << (rot.a*fiftyseven) << "deg about ligand " << ligbb[0]->name << "." << endl;
                                break;

                            case 2:
                                // Rotisserie.
                                axis = ligbb[1]->get_location().subtract(origin);

                                for (theta=0; theta < M_PI*2; theta += rstep)
                                {
                                    lv.copy(axis);
                                    lv.origin = origin;
                                    m.rotate(lv, rstep);

                                    float r2 = alca->get_location().get_3d_distance(pt);
                                    if (r2 < bestr)
                                    {
                                        bestr = r2;
                                        besttheta = theta;
                                    }
                                }

                                lv.copy(axis);
                                lv.origin = origin;
                                m.rotate(lv, besttheta);
                                cout << "# Pivoted ligand " << (besttheta*fiftyseven) << "deg about ligand " << ligbb[0]->name
                                     << "-" << ligbb[1]->name << " axis." << endl;
                                break;

                            default:
                                ;
                            }

                            // Preemptively minimize intermol clashes.
                            Molecule* mtmp[3];
                            mtmp[0] = &m;
                            mtmp[1] = alignment_aa[l];
                            mtmp[2] = NULL;
                            m.movability = MOV_ALL;
                            alignment_aa[l]->movability = MOV_FLEXONLY;
                            Molecule::multimol_conform(mtmp);
                            // m.intermol_conform_norecen(alignment_aa[l], iters, reaches_spheroid[nodeno]);
                            // alignment_aa[l]->intermol_conform_norecen(&m, iters, reaches_spheroid[nodeno]);
                            if (debug) *debug << "Alignment atom " << l << " is "
                                                  << alignment_aa[l]->get_name() << ":" << alca->name
                                                  << " Z " << alca->get_Z() << endl;
                        }
                    }
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
            for (j=0; j<sphres; j++)
            {
                if (reaches_spheroid[nodeno][j]->movability >= MOV_FLEXONLY) reaches_spheroid[nodeno][j]->movability = MOV_FLEXONLY;
                cfmols[i++] = reaches_spheroid[nodeno][j];
            }
            for (; i<SPHREACH_MAX; i++)
                cfmols[i] = NULL;

            // time_t preiter = time(NULL);
            if (differential_dock)
            {
                Molecule** delete_me;
                Molecule::multimol_conform(
                    cfmols,
                    delete_me = p.all_residues_as_molecules(),
                    iters,
                    &iteration_callback
                );
                delete[] delete_me;
            }
            else
            {
                Molecule::multimol_conform(
                    cfmols,
                    iters,
                    &iteration_callback
                );
            }
            /*time_t jlgsux = time(NULL);
            cout << "\nIterations took: " << (jlgsux-preiter) << " seconds." << endl;*/

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

            char metrics[p.get_seq_length()+8][10];
            float mkJmol[p.get_seq_length()+8];
            float imkJmol[p.get_seq_length()+8];
            float mvdWrepl[p.get_seq_length()+8];
            float imvdWrepl[p.get_seq_length()+8];
            int metcount=0;
            float btot=0;

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

            std::vector<AminoAcid*> allres = p.get_residues_near(pocketcen, 10000);
            qpr = allres.size();
            Molecule* postaa[seql+8];
            postaa[0] = ligand;
            for (i=0; i<qpr; i++)
            {
                postaa[i+1] = reinterpret_cast<Molecule*>(allres[i]);
            }

            for (i=0; i<_INTER_TYPES_LIMIT; i++) total_binding_by_type[i] = 0;

            if (differential_dock)
            {
                for (i=0; i<qpr+1; i++)
                {
                    int resno = i ? (allres[i-1]->get_residue_no()) : 0;
                    for (j=0; j<qpr+1; j++)
                    {
                        if (j == i) continue;
                        final_binding[resno] += postaa[i]->get_intermol_binding(postaa[j], j==0);
                        final_vdWrepl[resno] += postaa[i]->get_vdW_repulsion(postaa[j]);
                    }
                }
            }

            float fin_total_binding_by_type[_INTER_TYPES_LIMIT];
            for (i=0; i<_INTER_TYPES_LIMIT; i++) fin_total_binding_by_type[i] = total_binding_by_type[i];

            sphres = p.get_residues_can_clash_ligand(reaches_spheroid[nodeno], &m, m.get_barycenter(), size, mcoord_resno);
            // cout << "sphres " << sphres << endl;
            for (i=0; i<sphres; i++)
            {
                if (!reaches_spheroid[nodeno][i]) continue;
                if (!p.aa_ptr_in_range(reaches_spheroid[nodeno][i])) continue;
                reaches_spheroid[nodeno][i]->clear_atom_binding_energies();
                int resno = reaches_spheroid[nodeno][i]->get_residue_no();

                float lb;
                if (differential_dock)
                {
                    lb = final_binding[resno];
                }
                else
                {
                    lb = m.get_intermol_binding(reaches_spheroid[nodeno][i], false /* i==0 */);
                }
                if (lb > 90) lb = 0;

                sprintf(metrics[metcount], "%s%d", reaches_spheroid[nodeno][i]->get_3letter(), resno);
                // cout << metrics[metcount] << ": " << lb << " . ";
                mkJmol[metcount] = lb;
                imkJmol[metcount] = initial_binding[resno];

                if (differential_dock)
                {
                    mvdWrepl[metcount] = final_vdWrepl[resno];
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
                }

                imvdWrepl[metcount] = initial_vdWrepl[resno];
                metcount++;
                btot += lb;
                // cout << *(reaches_spheroid[nodeno][i]) << " adds " << lb << " to btot, making " << btot << endl;
            }
            // cout << btot << endl;

            if (btot > 15*m.get_atom_count()) btot = 0;

            drcount = pose-1;

            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Prepared metrics." << endl;
            #endif

            // Allocate the array.
            dr[drcount][nodeno].kJmol = 0; //btot;
            dr[drcount][nodeno].ikJmol = 0;
            dr[drcount][nodeno].metric  = new char*[metcount+4];
            dr[drcount][nodeno].mkJmol   = new float[metcount];
            dr[drcount][nodeno].imkJmol   = new float[metcount];
            dr[drcount][nodeno].mvdWrepl   = new float[metcount];
            dr[drcount][nodeno].imvdWrepl   = new float[metcount];
            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Allocated memory." << endl;
            #endif

            for (i=0; i<_INTER_TYPES_LIMIT; i++)
            {
                dr[drcount][nodeno].bytype[i] = fin_total_binding_by_type[i];
                dr[drcount][nodeno].ibytype[i] = init_total_binding_by_type[i];
                dr[drcount][nodeno].ikJmol += init_total_binding_by_type[i];
                dr[drcount][nodeno].kJmol += fin_total_binding_by_type[i];
            }
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
            for (l=0; l<n; l++) m.get_atom(l)->stream_pdb_line(pdbdat, 9000+l);
            #if _DBG_STEPBYSTEP
            if (debug) *debug << "Prepared ligand PDB." << endl;
            #endif

            if (flex)
            {
                for (k=0; reaches_spheroid[nodeno][k]; k++)
                {
                    if (!p.aa_ptr_in_range(reaches_spheroid[nodeno][k])) continue;
                    n = reaches_spheroid[nodeno][k]->get_atom_count();
                    for (l=0; l<n; l++)
                    {
                        reaches_spheroid[nodeno][k]->get_atom(l)->stream_pdb_line(
                            pdbdat,
                            reaches_spheroid[nodeno][k]->atno_offset+l
                        );
                    }
                }
                #if _DBG_STEPBYSTEP
                if (debug) *debug << "Prepared flex PDBs." << endl;
                #endif
            }

            dr[drcount][nodeno].pdbdat = pdbdat.str();
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
                        		(dr[i][0].kJmol - dr[i][0].ikJmol) < (dr[drcount][nodeno].kJmol - dr[drcount][nodeno].ikJmol)
                           )
                           ||
                           (	!differential_dock
                        		&&
                        		dr[i][0].kJmol < btot
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
            }

            drcount = pose;

            // For performance reasons, once a path node (including #0) fails to meet the binding energy threshold, discontinue further
            // calculations for this pose.
            if (btot < kJmol_cutoff) break;
        }	// nodeno loop.
    } // pose loop.
    #if _DBG_STEPBYSTEP
    if (debug) *debug << "Finished poses." << endl;
    #endif

    // Output the dr[][] array in order of increasing pose number.
    cout << endl;
    if (output) *output << endl;

    const float energy_mult = kcal ? _kcal_per_kJ : 1;
    for (i=1; i<=poses; i++)
    {
        for (j=0; j<poses; j++)
        {
            if (dr[j][0].pose == i)
            {
                if (dr[j][0].kJmol >= kJmol_cutoff)
                {
                    for (k=0; k<=pathnodes; k++)
                    {
                        // If pathnode is not within kJ/mol cutoff, abandon it and all subsequent pathnodes of the same pose.
                        if (dr[j][k].kJmol < kJmol_cutoff) break;

                        cout << "Pose: " << i << endl << "Node: " << k << endl;
                        if (output) *output << "Pose: " << i << endl << "Node: " << k << endl;

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
                                cout << dr[j][k].metric[l] << ": " << -dr[j][k].mkJmol[l]*energy_mult << endl;
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
                                cout << lbtyp << -dr[j][k].bytype[l]*energy_mult << endl;
                                if (output) *output << lbtyp << -dr[j][k].bytype[l]*energy_mult << endl;
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
                            cout << "Total: " << -dr[j][k].kJmol*energy_mult << endl << endl;
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

                        if (!dr[j][k].pdbdat.length())
                        {
                            cout << "WARNING: Failed to generate PDB data." << endl;
                            if (output) *output << "(Missing PDB data.)" << endl;
                        }
                        else
                        {
                            cout << "# PDB Data" << endl << "PDBDAT:" << endl;
                            if (output) *output << "# PDB Data" << endl << "PDBDAT:" << endl;

                            if (output) *output << dr[j][k].pdbdat << endl;
                            cout << dr[j][k].pdbdat << endl;

                            cout << "TER" << endl << "END" << endl << endl << endl;
                            if (output) *output << "TER" << endl << "END" << endl << endl << endl;
                        }
                    }
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
    cout << (i-1) << " pose(s) found." << endl;
    if (output) *output << (i-1) << " pose(s) found." << endl;
    if (debug) *debug << (i-1) << " pose(s) found." << endl;

    if (met) delete met;

    time_t finished = time(NULL);
    cout << "\nCalculation took: " << (finished-began) << " seconds." << endl;
    if (output) *output << "\nCalculation took: " << (finished-began) << " seconds." << endl;
    if (debug) *debug << "\nCalculation took: " << (finished-began) << " seconds." << endl;

    if (debug) debug->close();
    if (output) output->close();

    return 0;
}


















