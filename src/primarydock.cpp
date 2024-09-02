
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <regex>
#include "classes/cavity.h"
#include "classes/group.h"
#include "classes/scoring.h"

using namespace std;

Protein* protein;
Molecule* ligand;
Molecule* waters;
int nwater = 0;
bool progressbar = true;
bool output_each_iter = false;
bool flex = true;
int movie_offset = 0;
int poses = 10;
int pose = 1;
int nodeno = 0;
float elim = -0.01;
char* protfname;
char* ligfname;
AminoAcid* reaches_spheroid[10][SPHREACH_MAX+16];
AminoAcid* priorities[64];
int spinchr = 0;
float hueoffset = 0;
bool kcal = false;
char outfile[1024];
Atom *ligbba, *resbba;
bool output_missed_connections = false;
bool output_vdW_repulsions = false;

int interpret_resno(const char* field)
{
    char buffer[strlen(field)+4];
    strcpy(buffer, field);

    char* offset = buffer;
    while (*offset >= 'A' && *offset <= 'Z') offset++;

    int retval = 0;
    char* dot = strchr(buffer, '.');
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

    if (offset == buffer) return retval;

    AminoAcid* aa = protein->get_residue(retval);
    if (!aa) return 0;

    int i;
    for (i=0; buffer[i] >= 'A' && buffer[i] <= 'Z'; i++)
    {
        if (buffer[i] == aa->get_letter()) aa->priority = true;
    }
    return retval;
}

void output_iter(int iter, Molecule** mols)
{
    Point size(_INTERA_R_CUTOFF, _INTERA_R_CUTOFF, _INTERA_R_CUTOFF);
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

        DockResult ldr(protein, ligand, size, nullptr, pose);
        ldr.include_pdb_data = false;
        ldr.display_clash_atom1 = true;
        ldr.display_clash_atom2 = true;
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

void iteration_callback(int iter, Molecule** mols)
{
    int i;
    for (i=0; priorities[i]; i++)
    {
        Atom *a, *b;
        priorities[i]->mutual_closest_atoms(ligand, &a, &b);
        float r = a->distance_to(b);
        float opt = InteratomicForce::optimal_distance(a, b);
        if (r > opt)
        {
            Point pt = a->get_location().subtract(b->get_location());
            pt.multiply(0.333);
            ligand->increment_lm(pt);
        }
    }

    if (output_each_iter) output_iter(iter, mols);
}

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

int main(int argc, char** argv)
{
    // Splash
    cout << "\n                                                                                      __       ____  \npppp                                            ddd                               ,-_/  `-_--_/    \\  \np   p         i                                 d  d                 k            )                (__   \np   p                                           d   d                k           )   ()    __/        )   \npppp  r rrr  iii  mmm mm   aaaa   r rrr  y   y  d   d   ooo    ccc   k   k      /      \\__/  \\__/    /  \np     rr      i   m  m  m      a  rr     y   y  d   d  o   o  c   c  k  k      (       /  \\__/  \\   (  \np     r       i   m  m  m   aaaa  r      y   y  d   d  o   o  c      blm        \\    ()        _     )  \np     r       i   m  m  m  a   a  r      y   y  d  d   o   o  c   c  k  k        )     __     / \\   /  \np     r      iii  m  m  m   aaaa  r       yyyy  ddd     ooo    ccc   k   k       \\____/  `---'   \\__)  \n                                             y\n                                       yyyyyy\n\n";

    Point size(_INTERA_R_CUTOFF, _INTERA_R_CUTOFF, _INTERA_R_CUTOFF);
    int iters = 50;
    Protein p("TheProtein");
    Molecule m("ligand");
    char buffer[65536];
    protein = &p;
    ligand = &m;
    outfile[0] = 0;
    char outpdb[1024];
    int outpdb_poses = 0;
    bool test = false;

    protfname = ligfname = nullptr;
    waters = nullptr;

    FILE* fp;

    int i, j, k, l, n;
    for (i=0; i<64; i++) priorities[i] = nullptr;
    buffer[0] = 0;
    for (i=1; i<argc; i++)
    {
        if (buffer[0] != '-' || !buffer[1]) strcpy(buffer, argv[i]);
        if (!strcmp(buffer, "-p") || !strcmp(buffer, "--prot") || !strcmp(buffer, "--protein"))
        {
            i++;
            if (file_exists(argv[i]))
            {
                protfname = argv[i];
                fp = fopen(argv[i], "rb");
                p.load_pdb(fp);
            }
            else
            {
                cout << "File not found: " << argv[i] << endl;
                return -1;
            }
        }
        else if ((buffer[0] == '-' && buffer[1] == 'l') || !strcmp(buffer, "--lig") || !strcmp(buffer, "--ligand"))
        {
            i++;
            if (file_exists(argv[i]))
            {
                ligfname = argv[i];
                fp = fopen(argv[i], "rb");
                if (fp)
                {
                    fread(buffer, 1, 65536, fp);
                    fclose(fp);
                    m.from_sdf(buffer);
                }
                else
                {
                    cout << "Failed to open " << argv[i] << " for reading." << endl;
                    return -1;
                }
            }
            else
            {
                cout << "File not found: " << argv[i] << endl;
                return -1;
            }
        }
        else if ((buffer[0] == '-' && buffer[1] == 'd') || !strcmp(buffer, "--model"))
        {
            i++;
            outpdb_poses = atoi(argv[i]);

            i++;
            strcpy(outpdb, argv[i]);
        }
        else if ((buffer[0] == '-' && buffer[1] == 'b') || !strcmp(buffer, "--bsr"))           // Binding Site Residues.
        {
            if (!p.get_seq_length())
            {
                cout << "Protein must be specified before binding site residues." << endl;
                return -1;
            }

            i++;
            l = 0;
            while (i < argc)
            {
                if (argv[i][0] == '-')
                {
                    i--;
                    break;
                }
                int resno = interpret_resno(argv[i]);
                if (!resno) continue;
                AminoAcid* aa = p.get_residue(resno);
                if (aa)
                {
                    aa->priority = true;
                    if (l<60) priorities[l++] = aa;
                    priorities[l] = nullptr;
                }
                i++;
            }
        }
        else if ((buffer[0] == '-' && buffer[1] == 'm') || !strcmp(buffer, "--metal"))
        {
            if (!p.get_seq_length())
            {
                cout << "Protein must be specified before metal coordination." << endl;
                return -1;
            }

            i++;
            char* esym = argv[i];
            i++;
            int charge = atoi(argv[i]);
            i++;
            int mcr[16], mcrq=0;
            Point mcraloc[16];

            while (i < argc)
            {
                if (argv[i][0] == '-')
                {
                    i--;
                    break;
                }
                int resno = interpret_resno(argv[i]);
                if (!resno) continue;
                AminoAcid* aa = p.get_residue(resno);
                if (!aa) continue;
                Atom* a = aa->get_one_most_bindable(mcoord);
                if (!a) continue;
                mcraloc[mcrq] = a->get_location();
                mcr[mcrq++] = resno;
                i++;
                if (mcrq > 16)
                {
                    cout << "Maximum number of metal coordination residues exceeded." << endl;
                    return -1;
                }
            }

            if (mcrq)
            {
                Point mcen = average_of_points(mcraloc, mcrq);
                cout << "Performing metal coordination (3 residues)..." << endl << flush;
                std::vector<MCoord> mtlcoords;
                MCoord mc;
                mc.Z = Atom::Z_from_esym(esym);
                mc.charge = charge;
                mc.mtl = new Atom(esym, &mcen, charge);
                mc.mtl->name = new char[16];
                strcpy(mc.mtl->name, "MTL");

                for (l=0; l<mcrq; l++)
                {
                    ResiduePlaceholder rp;
                    rp.resno = mcr[l];
                    mc.coordres.push_back(rp);
                    AminoAcid* aa = p.get_residue(rp.resno);
                    Atom* a = aa->get_one_most_bindable(mcoord);
                    aa->conform_atom_to_location(a->name, mcen);
                }

                mtlcoords.push_back(mc);

                p.coordinate_metal(mtlcoords);
            }
        }
        else if ((buffer[0] == '-' && buffer[1] == 'i') || !strcmp(buffer, "--iter") || !strcmp(buffer, "--iters"))
        {
            i++;
            iters = atoi(argv[i]);
        }
        else if ((buffer[0] == '-' && buffer[1] == 'n') || !strcmp(buffer, "--pose") || !strcmp(buffer, "--poses"))
        {
            i++;
            poses = atoi(argv[i]);
        }
        else if ((buffer[0] == '-' && buffer[1] == 'o') || !strcmp(buffer, "--out"))
        {
            i++;
            strcpy(outfile, argv[i]);
        }
        else if ((buffer[0] == '-' && buffer[1] == 'r') || !strcmp(buffer, "--rigid"))
        {
            flex = false;
        }
        else if ((buffer[0] == '-' && buffer[1] == 'v') || !strcmp(buffer, "--movie"))
        {
            output_each_iter = true;
        }
        else if ((buffer[0] == '-' && buffer[1] == 't') || !strcmp(buffer, "--test"))
        {
            test = true;
        }
        else if ((buffer[0] == '-' && buffer[1] == 'w') || !strcmp(buffer, "--water"))
        {
            i++;
            nwater = atoi(argv[i]);
            waters = new Molecule[n+4];
            for (j=0; j<nwater; j++)
            {
                waters[j].from_smiles("O");
            }
        }
        else if ((buffer[0] == '-' && buffer[1] == 's') || !strcmp(buffer, "--soft"))
        {
            // TODO: Soft dock.
            ;
        }
        else if ((buffer[0] == '-' && buffer[1] == 'x') || !strcmp(buffer, "--iso"))
        {
            // TODO: Ligand enantiomers.
            ;
        }
        else if ((buffer[0] == '-' && buffer[1] == 'a') || !strcmp(buffer, "--atomto"))
        {
            if (!p.get_seq_length())
            {
                cout << "Protein must be specified before pointing side chains." << endl;
                return -1;
            }

            // Conform side chain atoms to location.
            i++;
            int resno = interpret_resno(argv[i]);
            AminoAcid* aa = resno ? p.get_residue(resno) : nullptr;

            i++;
            Atom* pointing;
            if (aa && !strcasecmp(argv[i], "ext") || !strcasecmp(argv[i], "extent")) pointing = aa->get_reach_atom();
            else pointing = aa ? aa->get_atom(argv[i]) : nullptr;

            i++;
            int tres = interpret_resno(argv[i]);
            AminoAcid* target = tres ? p.get_residue(tres) : nullptr;
            Atom* tatom = target ? target->get_nearest_atom(pointing->get_location()) : nullptr;

            if (aa && pointing && target) aa->conform_atom_to_location(pointing->name, tatom->get_location(), 20, InteratomicForce::optimal_distance(pointing, tatom));
        }
        else if (!strcmp(buffer, "--outmc"))
        {
            output_missed_connections = true;
        }
        else if (!strcmp(buffer, "--outvdwr"))
        {
            output_vdW_repulsions = true;
        }
        else if ((buffer[0] == '-' && buffer[1] == 'q') || !strcmp(buffer, "--congress"))
        {
            progressbar = false;
        }
        else if ((buffer[0] == '-' && buffer[1] == 'e') || !strcmp(buffer, "--elim") || !strcmp(buffer, "--energy"))
        {
            i++;
            elim = atof(argv[i]);
        }
        else if ((buffer[0] == '-' && buffer[1] == 'c') || !strcmp(buffer, "--kcal"))
        {
            kcal = true;
        }
        else if ((buffer[0] == '-' && buffer[1] == 'h') || !strcmp(buffer, "--help"))
        {
            cout << "Usage:" << endl << endl;
            cout << "bin/primarydock -p path/to/protein.pdb -l path/to/ligand.sdf [-b binding site residues [-i iters [other options]]]" << endl << endl;
            cout << "Options can occur in any sequence except -b and -m can only occur after -p." << endl << endl;
            cout << "OPTIONS:" << endl;

            cout << "-a, --atomto\tSpecify side chain atoms to move toward other features." << endl;
            cout << "\t\tThe syntax is -a residue atom target, where atom is the name of a" << endl;
            cout << "\t\tside chain atom and target is a residue to point toward. In place" << endl;
            cout << "\t\tof a named atom, the word extent can be used to select the most" << endl;
            cout << "\t\tdistant atom from the residue's α carbon. The residue and target" << endl;
            cout << "\t\tcan be either a residue number or a Ballesteros-Weinstein number." << endl;
            cout << "\t\tThe target residue's nearest atom will be the effective target." << endl << endl;

            cout << "-b, --bsr\tSpecify binding site residues." << endl;
            cout << "\t\tCan be specified as residue numbers, e.g. -b 104 255, or as " << endl;
            cout << "\t\tBallesteros-Weinstein numbers, e.g. -b 3.33 6.54. Can be the " << endl;
            cout << "\t\tpreceded by one or more amino acid letters, meaning residue " << endl;
            cout << "\t\tmust match one of the specified aminos or it will not be " << endl;
            cout << "\t\tcounted as a BSR. Residues listed in this option are given " << endl;
            cout << "\t\tspecial priority for ligand binding." << endl << endl;

            cout << "-c, --kcal\tOutput energies in kcal/mol. Default is kJ/mol." << endl << endl;

            cout << "-d, --model\tSpecifies the filename of output PDB models for dock " << endl;
            cout << "\t\tposes. The syntax is -d num_poses path/to/file.pdb, where " << endl;
            cout << "\t\tnum_poses is a positive integer. The special codes %p and %l " << endl;
            cout << "\t\tcan be used as placeholders for the protein name and ligand " << endl;
            cout << "\t\tname, respectively, and PrimaryDock will automatically fill " << endl;
            cout << "\t\tthem in. Similarly, the %o placeholder will be auto-replaced " << endl;
            cout << "\t\twith the pose number in decreasing order of energetic " << endl;
            cout << "\t\tfavorability. This filename normally ends with a .pdb " << endl;
            cout << "\t\textension." << endl << endl;

            cout << "-e, --energy\tSet maximum energy for output poses." << endl;
            cout << "\t\tCandidate poses not meeting this requirement will be omitted." << endl << endl;

            cout << "-h, --help\tShows this help screen." << endl << endl;

            cout << "-i, --iter\tSets the number of dock iterations." << endl;
            cout << "\t\tHigher numbers will produce better results, but at diminishing " << endl;
            cout << "\t\tgains and longer processing time. Too few iterations will " << endl;
            cout << "\t\tproduce inaccurate results and/or failed docks." << endl << endl;

            cout << "-l, --ligand\tSpecifies a file in SDF format for the ligand to be docked." << endl << endl;

            cout << "-m, --metal\tSpecifies a metal coordination site. Can occur multiple " << endl;
            cout << "\t\ttimes. The format is: -m element charge coordinating_residues. " << endl;
            cout << "\t\tExample: -m Cu +1 M5.36 C5.42 C5.43." << endl;
            cout << "\t\tEither sequence numbers or BW numbers can be used, and the " << endl;
            cout << "\t\tamino acid letters are optional." << endl << endl;

            cout << "-n, --poses\tSets the maximum number of output poses." << endl << endl;

            cout << "-o, --out\tSpecifies an output filename for dock results. The special " << endl;
            cout << "\t\tcodes %p and %l will be replaced with the name of the protein " << endl;
            cout << "\t\tand ligand, respectively. This filename normally ends with a " << endl;
            cout << "\t\t.dock extension." << endl << endl;

            cout << "--outmc\tIncludes \"missed connections\", residues that could not optimally " << endl;
            cout << "\t\tbind to the ligand, in the output." << endl << endl;

            cout << "--outvdwr\tIncludes van der Waals repulsions between the ligand and nearby " << endl;
            cout << "\t\tresidues in the output." << endl << endl;

            cout << "-p, --protein\tSpecifies a file in PDB format for the protein to be docked." << endl << endl;

            cout << "-q\t\tSuppresses the progress bar." << endl << endl;

            cout << "-v, --movie\tCreates a dock file in the tmp/ folder showing each iteration " << endl;
            cout << "\t\tof the dock algorithm. Used for debugging." << endl << endl;

            cout << "-w, --water\tAdds a specified number of water molecules to the binding pocket, " << endl;
            cout << "\t\te.g. -w 5 to add 5 molecules." << endl << endl;

            return 0;
        }

        if (buffer[0] == '-' && buffer[1] != '-')
        {
            strcpy(buffer+1, buffer+2);
            i--;
        }
    }

    if (!p.get_seq_length())
    {
        cout << "No protein specified." << endl;
        return -1;
    }
    int seqlen = p.get_end_resno();

    cout << "PDB file: " << protfname << endl;
    if (ligfname) cout << "Ligand: " << ligfname << endl;

    char* pcntp = strstr(outfile, "%p");
    if (pcntp)
    {
        char tmp[4096], protn[64];
        *(pcntp++) = 0;
        *(pcntp++) = 0;
        strcpy(protn, strrchr(protfname, '/')+1);
        char* dot = strchr(protn, '.');
        if (dot) *dot = 0;
        sprintf(tmp, "%s%s%s", outfile, protn, pcntp);
        strcpy(outfile, tmp);
    }

    char* pcntl = strstr(outfile, "%l");
    if (pcntl)
    {
        char tmp[4096], lign[64];
        *(pcntl++) = 0;
        *(pcntl++) = 0;
        strcpy(lign, strrchr(ligfname, '/')+1);
        char* dot = strchr(lign, '.');
        if (dot) *dot = 0;
        sprintf(tmp, "%s%s%s", outfile, lign, pcntl);
        strcpy(outfile, tmp);
    }

    if (strlen(outfile)) cout << "Output file: " << outfile << endl;
    cout << endl;


    ////////////////////////////////////////////////////////////////////////////
    // Phase I: Cavity Search.                                                //
    ////////////////////////////////////////////////////////////////////////////

    time_t began = time(NULL);

    Cavity cavities[1029];
    int qfound = Cavity::scan_in_protein(&p, cavities, 1024);

    if (test) cout << "Found " << qfound << " cavit" << (qfound == 1 ? "y." : "ies.") << endl;

    fp = fopen("tmp/cavities.js", "w");
    if (fp)
    {
        for (i=0; i<qfound; i++)
        {
            cavities[i].output_ngl_js(fp);
            fprintf(fp, "\n\n\n\n\n");
        }
        fclose(fp);
    }


    ////////////////////////////////////////////////////////////////////////////
    // Phase II: Rotamer Search.                                             //
    ////////////////////////////////////////////////////////////////////////////

    Pose best(&m), secondbest(&m), thirdbest(&m);
    if (m.get_atom_count())
    {
        Pose rotamers[ligand_rotamer_count];
        ligand->movability = MOV_ALL;
        ligand->minimize_internal_clashes();
        for (i=0; i<ligand_rotamer_count; i++)
        {
            if (i >= 10) m.crumple(0.11*fiftyseventh);
            rotamers[i].copy_state(&m);
        }

        float bestviol = Avogadro, bestviol2 = Avogadro, bestviol3 = Avogadro;
        cout << "Trying ligand rotamers in pockets..." << flush;
        for (j=0; j<ligand_rotamer_count; j++)
        {
            ligand->movability = MOV_ALL;
            rotamers[j].restore_state(&m);
            for (i=0; i<qfound; i++)
            {
                Point cavcen = cavities[i].get_center();
                if (cavcen.y < 5 || cavcen.y > 18) continue;
                float viol = cavities[i].find_best_containment(&m);
                if (test) cout << endl << "Cavity centered at " << cavcen << " has a score of " << viol;
                if (viol < bestviol)
                {
                    thirdbest = secondbest;
                    secondbest = best;
                    best.copy_state(&m);
                    bestviol3 = bestviol2;
                    bestviol2 = bestviol;
                    bestviol = viol;
                }
                else if (viol < bestviol2)
                {
                    thirdbest = secondbest;
                    secondbest.copy_state(&m);
                    bestviol3 = bestviol2;
                    bestviol2 = viol;
                }
                else if (viol < bestviol3)
                {
                    thirdbest.copy_state(&m);
                    bestviol3 = viol;
                }
                cout << "." << flush;
            }
        }
        cout << endl;

        if (bestviol > 5e22)
        {
            cout << "No suitable pockets found." << endl;
            return 0;
        }

        best.restore_state(&m);
        if (test)
        {
            cout << "Best pre-iteration pose has a score of " << bestviol << " (negative = favorable)." << endl;
            cout << "Ligand centered at " << m.get_barycenter() << endl;
            fp = fopen("tmp/predock.pdb", "wb");
            if (!fp) return -1;
            p.save_pdb(fp, &m);
            cout << "Saved output file in tmp/." << endl;
        }
    }
    else return 0;


    ////////////////////////////////////////////////////////////////////////////
    // Phase III: Refinement.                                                 //
    ////////////////////////////////////////////////////////////////////////////

    Pose residue_init[seqlen+1];
    Pose tmp_pdb_waters[poses+1][nwater+1];
    Point tmp_pdb_metal_locs[poses+1][mtlcoords.size()+1];
    for (i=1; i<=seqlen; i++)
    {
        AminoAcid* aa = protein->get_residue(i);
        if (aa) residue_init[i].copy_state(aa);
    }

    cout << "Generating poses and refining..." << endl << endl << flush;
    Pose candidates[poses];
    Pose residue_candidates[poses][seqlen+1];
    for (pose=1; pose<=poses; pose++)
    {
        j = (rand() % 3);
        ligand->movability = MOV_ALL;

        if (j==0) best.restore_state(ligand);
        else if (j==1) secondbest.restore_state(ligand);
        else thirdbest.restore_state(ligand);
        for (i=1; i<=seqlen; i++)
        {
            AminoAcid* aa = protein->get_residue(i);
            if (aa) residue_init[i].restore_state(aa);
        }

        Point ligcen = ligand->get_barycenter();
        int sphres = p.get_residues_can_clash_ligand(reaches_spheroid[nodeno], ligand, ligcen, size);
        
        Molecule* cfmols[sphres+8];
        j=0;
        cfmols[j++] = ligand;
        Molecule* met = p.metals_as_molecule();
        if (met) cfmols[j++] = met;
        if (nwater && waters) for (i=0; i<nwater; i++) cfmols[j++] = &waters[i];
        for (i=0; i<sphres; i++) cfmols[j++] = reinterpret_cast<Molecule*>(reaches_spheroid[nodeno][i]);
        cfmols[j] = nullptr;

        if (nwater && waters) for (i=0; i<nwater; i++)
        {
            for (l=nwater+2; l<j; l++)
            {
                if (fabs(cfmols[l]->hydrophilicity()) > hydrophilicity_cutoff && frand(0,1) < 0.333)
                {
                    Point pt = cfmols[l]->get_single_most_bindable(hbond)->get_location();
                    SCoord voffset(2.5, frand(0, M_PI), frand(0, M_PI*2));
                    waters[i].recenter(pt.add(voffset));
                }
            }
        }

        // Identify the mutual best atom pair for binding and preemptively move the ligand in place.
        intera_type bbt = vdW;
        if (met && met->get_atom_count())
        {
            bbt = mcoord;
            ligbba = ligand->get_single_most_bindable(bbt);
            resbba = met->get_nearest_atom(ligbba->get_location());
            ligand->maintain_contact(ligbba, resbba);
        }
        else
        {
            float nr = Avogadro;
            for (i=0; i<sphres; i++)
            {
                Atom* ra = reaches_spheroid[nodeno][i]->get_nearest_atom(ligcen);
                if (!ra) continue;
                Atom* la = ligand->get_nearest_atom(ra->get_location());
                float r = ra->distance_to(la);
                if (reaches_spheroid[nodeno][i]->priority) r /= 10;

                intera_type lt = Molecule::best_binding_type(reinterpret_cast<Molecule*>(reaches_spheroid[nodeno][i]), ligand);
                if (lt == mcoord) r /= 200;
                else if (lt == ionic) r /= 60;
                else if (lt == hbond) r /= 25;
                else if (lt == pi) r /= 12;
                else if (lt == polarpi) r /= 8;

                if (r < nr)
                {
                    nr = r;
                    resbba = ra;
                    ligbba = la;
                    bbt = lt;
                }
            }

            if (resbba && ligbba)
            {
                AminoAcid* aa = p.get_residue(resbba->residue);
                resbba = aa->get_one_most_bindable(bbt);
                ligbba = ligand->get_single_most_bindable(bbt);
                float opt = InteratomicForce::optimal_distance(ligbba, resbba);
                // aa->conform_atom_to_location(resbba->name, ligbba->get_location(), opt);
                ligand->conform_atom_to_location(ligbba->name, resbba->get_location(), opt);
            }
        }
        if (ligbba && resbba)
        {
            ligcen = ligand->get_barycenter();
            Rotation rot = align_points_3d(ligbba->get_location(), resbba->get_location(), ligcen);
            LocatedVector lv = rot.v;
            lv.origin = ligcen;
            ligand->rotate(lv, rot.a);

            float r = ligbba->distance_to(resbba);
            float rcen = ligcen.get_3d_distance(resbba->get_location());
            float ll = ligcen.get_3d_distance(ligbba->get_location());
            float opt = InteratomicForce::optimal_distance(ligbba, resbba);
            SCoord motion = resbba->get_location().subtract(ligcen);
            motion.r = (rcen > ll) ? (r - opt*2) : (-opt*2);
            ligand->move(motion);

            ligand->maintain_contact(ligbba, resbba);
            if (resbba->mol) resbba->mol->maintain_contact(resbba, ligbba);

            if (test) cout << "Ligand " << ligbba->name << " can " << bbt << " with "
                << (resbba->mol ? resbba->mol->get_name() : "UNKNOWN")
                << ":" << resbba->name << "." << endl << endl;
        }
        else if (test) cout << "Mutual nearest atom adjustment failed." << endl << endl;

        ligand->movability = MOV_ALL;
        if (output_each_iter) output_iter(0, cfmols);
        Molecule::conform_molecules(cfmols, iters, &iteration_callback, &GroupPair::align_groups_noconform, progressbar ? &update_progressbar : nullptr);

        n = ligand->get_atom_count();
        for (i=0; i<n= i++)
        {
            Atom* a = ligand->get_atom(i);
            if (a) a->residue = 0;
        }
        candidates[pose-1].copy_state(ligand);
        for (i=1; i<=seqlen; i++)
        {
            AminoAcid* aa = protein->get_residue(i);
            if (aa) residue_candidates[pose-1][i].copy_state(aa);
        }
        if (nwater)
        {
            for (j=0; j<nwater; j++)
            {
                tmp_pdb_waters[pose-1][j].copy_state(&waters[j]);
            }
        }
        n = mtlcoords.size();
        for (j=0; j < n; j++)
        {
            tmp_pdb_metal_locs[pose][j] = mtlcoords[j].mtl->get_location();
        }
    }
    cout << "                                      \033[A                                                                                                                                      " << endl << "                                                                                                                                             " << endl;
 

    ////////////////////////////////////////////////////////////////////////////
    // Phase IV: Scoring.                                                     //
    ////////////////////////////////////////////////////////////////////////////

    DockResult* result[poses];
    for (pose=1; pose<=poses; pose++)
    {
        candidates[pose-1].restore_state(ligand);
        for (i=1; i<=seqlen; i++)
        {
            AminoAcid* aa = protein->get_residue(i);
            if (aa) residue_candidates[pose-1][i].restore_state(aa);
        }

        result[pose-1] = new DockResult(protein, ligand, size, nullptr, pose);
        result[pose-1]->energy_mult = kcal ? _kcal_per_kJ : 1;
        result[pose-1]->out_vdw_repuls = output_vdW_repulsions;
        result[pose-1]->out_mc = output_missed_connections;
        
        std::ostringstream pdbdat;
        n = ligand->get_atom_count();
        int offset = n;

        for (l=0; l<n; l++)
        {
            Atom* a = ligand->get_atom(l);
            if (!a) continue;
            a->residue = pose;
            a->stream_pdb_line(pdbdat, 9000+l, true);
        }

        int sphres = protein->get_residues_can_clash_ligand(reaches_spheroid[nodeno], ligand, ligand->get_barycenter(), size);

        int resno;
        for (resno = protein->get_start_resno(); resno <= seqlen; resno++)
        {
            AminoAcid* laa = protein->get_residue(resno);
            if (!laa) continue;
            if (!flex || !laa->been_flexed)
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

        if (mtlcoords.size())
        {
            for (l=0; l<mtlcoords.size(); l++)
            {
                mtlcoords[l].mtl->stream_pdb_line(
                    pdbdat,
                    9000+offset+l
                );
            }
            offset += l;
        }

        result[pose-1]->pdbdat = pdbdat.str();
    }


    ////////////////////////////////////////////////////////////////////////////
    // Phase V: Filtering.                                                    //
    ////////////////////////////////////////////////////////////////////////////

    bool meets_criteria[poses];
    for (pose=1; pose<=poses; pose++)
    {
        meets_criteria[pose-1] = (-result[pose-1]->kJmol*result[pose-1]->energy_mult) < elim;
        // if (-result[pose-1]->worst_nrg_aa > clash_limit_per_aa) meets_criteria[pose-1] = false;
    }


    ////////////////////////////////////////////////////////////////////////////
    // Phase VI: Sorting.                                                     //
    ////////////////////////////////////////////////////////////////////////////

    int sortidx[poses];
    float sorted[poses];
    for (i=0; i<poses; i++)
    {
        sortidx[i] = -1;
        sorted[i] = Avogadro;
        for (j=0; j<poses; j++)
        {
            if (!meets_criteria[j]) continue;
            float e = -result[j]->kJmol;
            if (   ((!j) || (e < sorted[i]))
                && ((!i) || (e >= sorted[i-1]))
                && ((!i) || (j != sortidx[i-1]))
               )
            {
                sorted[i] = e;
                sortidx[i] = j;
            }
        }
    }

    // for (i=0; i<poses; i++) cout << i << ": index " << sortidx[i] << " energy " << sorted[i] << endl;


    ////////////////////////////////////////////////////////////////////////////
    // Phase VII: Output generation.                                          //
    ////////////////////////////////////////////////////////////////////////////

    std::ofstream *output = nullptr;
    if (outfile[0]) output = new std::ofstream(outfile, std::ofstream::out);
    if (output) *output << "PDB file: " << protfname << endl;
    if (output) *output << "Ligand: " << ligfname << endl;
    if (output) *output << endl;

    int num_poses = 0;
    for (i=0; i<poses; i++)
    {
        pose = i+1;
        j = sortidx[i];
        if (j<0) break;

        cout << "Pose: " << pose << endl << "Node: " << 0 << endl << endl;
        if (output) *output << "Pose: " << pose << endl << "Node: " << 0 << endl << endl;

        result[j]->include_pdb_data = (output == nullptr);
        cout << *result[j] << endl << endl;
        result[j]->include_pdb_data = true;
        if (output) *output << *result[j];

        num_poses++;

        if (pose <= outpdb_poses)
        {
            char protn[64];
            strcpy(protn, strrchr(protfname, '/')+1);
            char* dot = strchr(protn, '.');
            if (dot) *dot = 0;

            char lign[64];
            strcpy(lign, strrchr(ligfname, '/')+1);
            dot = strchr(lign, '.');
            if (dot) *dot = 0;

            std::string out_pdb_fn = std::regex_replace(outpdb, std::regex("[%][p]"), protn);
            out_pdb_fn = std::regex_replace(out_pdb_fn, std::regex("[%][l]"), lign);
            out_pdb_fn = std::regex_replace(out_pdb_fn, std::regex("[%][o]"), to_string(pose));

            FILE* pfout = fopen(out_pdb_fn.c_str(), "wb");
            if (!pfout) cout << "Failed to open " << out_pdb_fn << " for writing." << endl;
            else
            {
                protein->save_pdb(pfout, ligand);

                int atno_offset = protein->last_saved_atom_number;
                if (nwater)
                {
                    for (j=0; j<nwater; j++)
                    {
                        tmp_pdb_waters[pose-1][j].restore_state(&waters[j]);
                        waters[j].save_pdb(pfout, atno_offset);
                        atno_offset += waters[j].get_atom_count();
                    }
                }

                n = mtlcoords.size();
                for (j=0; j < n; j++)
                {
                    mtlcoords[j].mtl->move(tmp_pdb_metal_locs[pose][j]);
                    mtlcoords[j].mtl->save_pdb_line(pfout, ++atno_offset);
                }

                protein->end_pdb(pfout);

                fclose(pfout);
            }
        }
    }

    cout << num_poses << " pose" << (i==1?"":"s") << " found." << endl;
    if (output) *output << num_poses << " pose" << (i==1?"":"s") << " found." << endl;
    if (num_poses)
    {
        cout << "Best energy: " << (-result[sortidx[0]]->kJmol*result[sortidx[0]]->energy_mult) << " "
            << (kcal?"kcal/mol":"kJ/mol") << "." << endl;
        if (output) *output << "Best energy: " << (-result[sortidx[0]]->kJmol*result[sortidx[0]]->energy_mult) << " "
            << (kcal?"kcal/mol":"kJ/mol") << "." << endl;
    }

    time_t finished = time(NULL);
    cout << "\nCalculation took: " << (finished-began) << " seconds." << endl;
    if (output) *output << "\nCalculation took: " << (finished-began) << " seconds." << endl;

    if (output)
    {
        output->close();
        for (i=1; i<=seqlen; i++)
        {
            AminoAcid* aa = protein->get_residue(i);
            if (aa) residue_candidates[max(0,sortidx[0])][i].restore_state(aa);
        }

        FILE* pf = fopen(outfile, "ab");
        fprintf(pf, "\nOriginal PDB:\n");
        protein->save_pdb(pf);
        fclose(pf);
        cout << "PDB appended to output file." << endl;
    }

    return 0;
}


