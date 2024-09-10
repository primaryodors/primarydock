
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

int main(int argc, char** argv)
{
    Point size(_INTERA_R_CUTOFF*0.5, _INTERA_R_CUTOFF*0.5, _INTERA_R_CUTOFF*0.5);
    int iters = 50;
    Protein p("TheProtein");
    Molecule m("ligand");
    char buffer[65536];
    protein = &p;
    ligand = &m;
    outfile[0] = 0;
    char outpdb[1024];
    int outpdb_poses = 0;

    protfname = ligfname = nullptr;
    waters = nullptr;

    FILE* fp;

    int i, j, k, l, n;
    for (i=0; i<64; i++) priorities[i] = nullptr;
    buffer[0] = 0;
    for (i=1; i<argc; i++)
    {
        if (buffer[0] != '-' || !buffer[1]) strcpy(buffer, argv[i]);
        if ((buffer[0] == '-' && buffer[1] == 'p') || !strcmp(buffer, "--prot") || !strcmp(buffer, "--protein"))
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
        else if ((buffer[0] == '-' && buffer[1] == 'h') || !strcmp(buffer, "--help"))
        {
            cout << "Usage:" << endl << endl;
            cout << "bin/cavity_search -p path/to/protein.pdb -l path/to/ligand.sdf [-b binding site residues [-i iters [other options]]]" << endl << endl;
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
            if (buffer[1]) i--;
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

    cout << "Found " << qfound << " cavit" << (qfound == 1 ? "y." : "ies.") << endl;

    fp = fopen("cavities.js", "w");
    if (fp)
    {
        for (i=0; i<qfound; i++)
        {
            cavities[i].output_ngl_js(fp);
            fprintf(fp, "\n\n\n\n\n");
        }
        fclose(fp);
    }

    return 0;
}


