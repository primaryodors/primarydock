
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "classes/cavity.h"
#include "classes/group.h"
#include "classes/scoring.h"

using namespace std;

Protein* protein;
Molecule* ligand;
bool progressbar = true;
bool output_each_iter = false;
bool flex = true;
int movie_offset = 0;
int poses = 10;
int pose = 1;
int nodeno = 0;
float elim = -0.01;
char* protfname;
AminoAcid* reaches_spheroid[10][SPHREACH_MAX+16];
AminoAcid* priorities[64];
int spinchr = 0;
float hueoffset = 0;
bool kcal = false;
char outfile[1024];

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
    for (i=0; i<80; i++)
    {
        float cmpi = 1.25*i;
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
    cout << "\n                                                                                      __       ____  \npppp                                            ddd                  k            ,-_/  `-_--_/    \\  \np   p         i                                 d  d                 k            )                (__   \np   p                                           d   d                k           )   ()    __/        )   \npppp  r rrr  iii  mmm mm   aaaa   r rrr  y   y  d   d   ooo    ccc   k   k      /      \\__/  \\__/    /  \np     rr      i   m  m  m      a  rr     y   y  d   d  o   o  c   c  k  k      (       /  \\__/  \\   (  \np     r       i   m  m  m   aaaa  r      y   y  d   d  o   o  c      blm        \\    ()        _     )  \np     r       i   m  m  m  a   a  r      y   y  d  d   o   o  c   c  k  k        )     __     / \\   /  \np     r      iii  m  m  m   aaaa  r       yyyy  ddd     ooo    ccc   k   k       \\____/  `---'   \\__)  \n                                             y\n                                       yyyyyy\n\n";

    Point size(_INTERA_R_CUTOFF, _INTERA_R_CUTOFF, _INTERA_R_CUTOFF);
    bool verbose = false;
    int iters = 50;
    Protein p("TheProtein");
    Molecule m("ligand");
    char buffer[65536];
    protein = &p;
    ligand = &m;
    outfile[0] = 0;

    FILE* fp;

    int i, j, k, l, n;
    for (i=0; i<64; i++) priorities[i] = nullptr;
    for (i=1; i<argc; i++)
    {
        if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--prot") || !strcmp(argv[i], "--protein"))
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
        else if (!strcmp(argv[i], "-l") || !strcmp(argv[i], "--lig") || !strcmp(argv[i], "--ligand"))
        {
            i++;
            if (file_exists(argv[i]))
            {
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
        else if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "--bsr"))           // Binding Site Residues.
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
        else if (!strcmp(argv[i], "-m") || !strcmp(argv[i], "--metal"))
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
        else if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--iter") || !strcmp(argv[i], "--iters"))
        {
            i++;
            iters = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "-n") || !strcmp(argv[i], "--pose") || !strcmp(argv[i], "--poses"))
        {
            i++;
            poses = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--out"))
        {
            i++;
            strcpy(outfile, argv[i]);
        }
        else if (!strcmp(argv[i], "-r") || !strcmp(argv[i], "--rigid"))
        {
            flex = false;
        }
        else if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--movie"))
        {
            output_each_iter = true;
        }
        else if (!strcmp(argv[i], "-q") || !strcmp(argv[i], "--congress"))
        {
            progressbar = false;
        }
        else if (!strcmp(argv[i], "-e") || !strcmp(argv[i], "--elim") || !strcmp(argv[i], "--energy"))
        {
            i++;
            elim = atof(argv[i]);
        }
        else if (!strcmp(argv[i], "-c") || !strcmp(argv[i], "--kcal"))
        {
            kcal = true;
        }
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help"))
        {
            cout << "Usage:" << endl << endl;
            cout << "bin/primarydock -p path/to/protein.pdb -l path/to/ligand.sdf [-b binding site residues [-i iters [other options]]]" << endl << endl;
            cout << "Options can occur in any sequence except -p must occur before -b and -m." << endl << endl;
            cout << "OPTIONS:" << endl;
            cout << "-b, --bsr\tSpecify binding site residues." << endl;
            cout << "\t\tCan be specified as residue numbers, e.g. -b 104 255," << endl;
            cout << "\t\tor as Ballesteros-Weinstein numbers, e.g. -b 3.33 6.54." << endl;
            cout << "\t\tCan be preceded by one or more amino acid letters, meaning" << endl;
            cout << "\t\tresidue must match one of the specified letters or it will" << endl;
            cout << "\t\tnot be counted as a BSR." << endl << endl;

            cout << "-c, --kcal\tOutput energies in kcal/mol. Default is kJ/mol." << endl << endl;

            cout << "-e, --energy\tSet maximum energy for output poses." << endl;
            cout << "\t\tCandidate poses not meeting this requirement will be omitted." << endl << endl;

            cout << "-h, --help\tShows this help screen." << endl << endl;

            cout << "-i, --iter\tSets the number of dock iterations." << endl;
            cout << "\t\tHigher numbers will produce better results, but at diminishing" << endl;
            cout << "\t\tgains and longer processing time. Too few iterations will produce" << endl;
            cout << "\t\tinaccurate results and/or failed docks." << endl << endl;

            cout << "-l, --ligand\tSpecifies a file in SDF format for the ligand to be docked." << endl << endl;

            cout << "-m, --metal\tSpecifies a metal coordination site. Can occur multiple times." << endl;
            cout << "\t\tThe format is: -m element_symbol charge coordinating_residues." << endl;
            cout << "\t\tExample: -m Cu +1 M5.36 C5.42 C5.43." << endl;
            cout << "\t\tEither sequence numbers or BW numbers can be used, and the amino acid" << endl;
            cout << "\t\tletters are optional." << endl << endl;

            cout << "-n, --poses\tSets the maximum number of output poses." << endl << endl;

            cout << "-o, --out\tSpecifies an output filename for dock results." << endl << endl;

            cout << "-p, --protein\tSpecifies a file in PDB format for the protein to be docked." << endl << endl;

            cout << "-q\t\tSuppresses the progress bar." << endl << endl;

            cout << "-v, --movie\tCreates a dock file in the tmp/ folder showing each" << endl;
            cout << "\t\titeration of the dock algorithm. Used for debugging." << endl << endl;

            return 0;
        }
    }

    if (!p.get_seq_length())
    {
        cout << "No protein specified." << endl;
        return -1;
    }
    int seqlen = p.get_end_resno();

    ////////////////////////////////////////////////////////////////////////////
    // Phase I: Cavity Search.                                                //
    ////////////////////////////////////////////////////////////////////////////

    time_t began = time(NULL);

    Cavity cavities[1029];
    int qfound = Cavity::scan_in_protein(&p, cavities, 1024);

    cout << "Found " << qfound << " cavit" << (qfound == 1 ? "y." : "ies.") << endl;

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
    // Phase II: Fleximer Search.                                             //
    ////////////////////////////////////////////////////////////////////////////

    Pose best(&m), secondbest(&m), thirdbest(&m);
    if (m.get_atom_count())
    {
        Pose fleximers[10];
        for (i=0; i<10; i++)
        {
            m.crumple(7.5*fiftyseventh);
            fleximers[i].copy_state(&m);
        }

        float bestviol = Avogadro, bestviol2 = Avogadro, bestviol3 = Avogadro;
        cout << "Trying ligand in pockets..." << flush;
        for (j=0; j<10; j++)
        {
            fleximers[j].restore_state(&m);
            for (i=0; i<qfound; i++)
            {
                Point cavcen = cavities[i].get_center();
                if (cavcen.y < 5 || cavcen.y > 18) continue;
                float viol = cavities[i].find_best_containment(&m);
                if (verbose) cout << endl << "Cavity centered at " << cavcen << " has a score of " << viol;
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
        cout << "Best pre-iteration pose has a score of " << bestviol << "." << endl;
        cout << "Ligand centered at " << m.get_barycenter() << endl;
        fp = fopen("tmp/quickdock.pdb", "wb");
        if (!fp) return -1;
        p.save_pdb(fp, &m);
        cout << "Saved output file in tmp/." << endl;
    }
    else return 0;


    ////////////////////////////////////////////////////////////////////////////
    // Phase III: Refinement.                                                 //
    ////////////////////////////////////////////////////////////////////////////

    Pose residue_init[seqlen+1];
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

        if (j==0) best.restore_state(ligand);
        else if (j==1) secondbest.restore_state(ligand);
        else thirdbest.restore_state(ligand);
        for (i=1; i<=seqlen; i++)
        {
            AminoAcid* aa = protein->get_residue(i);
            if (aa) residue_init[i].restore_state(aa);
        }

        int sphres = p.get_residues_can_clash_ligand(reaches_spheroid[nodeno], ligand, ligand->get_barycenter(),
            Point(_INTERA_R_CUTOFF, _INTERA_R_CUTOFF, _INTERA_R_CUTOFF));
        
        Molecule* cfmols[sphres+8];
        j=0;
        cfmols[j++] = ligand;
        Molecule* met = p.metals_as_molecule();
        if (met) cfmols[j++] = met;
        for (i=0; i<sphres; i++) cfmols[j++] = reinterpret_cast<Molecule*>(reaches_spheroid[nodeno][i]);
        cfmols[j] = nullptr;

        if (output_each_iter) output_iter(0, cfmols);
        Molecule::conform_molecules(cfmols, iters, &iteration_callback, &GroupPair::align_groups_noconform, progressbar ? &update_progressbar : nullptr);

        candidates[pose-1].copy_state(ligand);
        for (i=1; i<=seqlen; i++)
        {
            AminoAcid* aa = protein->get_residue(i);
            if (aa) residue_candidates[pose-1][i].copy_state(aa);
        }
    }

    /*
    fp = fopen("tmp/refined.pdb", "wb");
    if (!fp) return -1;
    p.save_pdb(fp, &m);
    cout << "Saved output file in tmp/." << endl;
    */


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
        if (-result[pose-1]->worst_nrg_aa > clash_limit_per_aa) meets_criteria[pose-1] = false;
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
    for (i=0; i<poses; i++)
    {
        pose = i+1;
        j = sortidx[i];
        if (j<0) break;

        cout << "Pose: " << pose << endl << "Node: " << 0 << endl;
        if (output) *output << "Pose: " << pose << endl << "Node: " << 0 << endl;

        result[j]->include_pdb_data = (output == nullptr);
        cout << *result[j] << endl << endl;
        result[j]->include_pdb_data = true;
        if (output) *output << *result[j] << endl << endl;
    }

    cout << i << " pose" << (i==1?"":"s") << " found." << endl;
    if (i) cout << "Best energy: " << (-result[sortidx[0]]->kJmol*result[sortidx[0]]->energy_mult) << " "
        << (kcal?"kcal/mol":"kJ/mol") << "." << endl;

    time_t finished = time(NULL);
    cout << "\nCalculation took: " << (finished-began) << " seconds." << endl;

    return 0;
}


