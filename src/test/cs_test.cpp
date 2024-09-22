
#include "../classes/search.h"

int main(int argc, char** argv)
{
    Protein p("Test Receptor");
    Molecule m("Test Ligand");
    BallesterosWeinstein pocketcen_res[256];
    int nbsr=0;
    std::string atoms_of_interest[256];
    int naoi=0;
    bool priorities[256];

    bool save_tmp_pdbs = false;

    int i;
    for (i=0; i<256; i++) priorities[i] = false;

    FILE* fp;
    char buffer[49152];
    for (i=1; i<argc; i++)
    {
        char* dot;
        dot = strrchr(argv[i], '.');
        if (dot && !strcasecmp(dot, ".pdb"))
        {
            fp = fopen(argv[i], "rb");
            if (!fp)
            {
                cout << "Failed to open " << argv[i] << " for reading." << endl;
                return -1;
            }
            p.load_pdb(fp);
            fclose(fp);
        }
        else if (dot && !strcasecmp(dot, ".sdf"))
        {
            fp = fopen(argv[i], "rb");
            if (!fp)
            {
                cout << "Failed to open " << argv[i] << " for reading." << endl;
                return -1;
            }
            fread(buffer, 1, 49150, fp);
            fclose(fp);
            m.from_sdf(buffer);
        }
        else if (!strcasecmp(argv[i], "aoi"))
        {
            for (i++; i<argc; i++)
            {
                if (!strcasecmp(argv[i], "ioa")) break;
                else if (!strcasecmp(argv[i], "*")) break;
                else atoms_of_interest[naoi++] = argv[i];
            }
        }
        else if (!strcasecmp(argv[i], "cen"))
        {
            for (i++; i<argc; i++)
            {
                if (!strcasecmp(argv[i], "nec")) break;
                else if (!strcasecmp(argv[i], "*")) break;
                else
                {
                    if (argv[i][strlen(argv[i])-1] == '!') priorities[nbsr] = true;
                    pocketcen_res[nbsr++] = argv[i];
                }
            }
        }
        else if (!strcasecmp(argv[i], "savtmp"))
        {
            save_tmp_pdbs = true;
        }
        else cout << "Warning: unknown command argument " << argv[i] << endl;
    }

    if (!p.get_end_resno())
    {
        cout << "ERROR: No protein." << endl;
        return -2;
    }

    if (!m.get_atom_count())
    {
        cout << "ERROR: No ligand." << endl;
        return -2;
    }

    Point pocketcen = p.get_region_center(1, p.get_end_resno());
    size = Point(100,100,100);

    if (p.get_bw50(6))
    {
        if (!nbsr)
        {
            pocketcen_res[nbsr++] = "3.33";
            pocketcen_res[nbsr++] = "3.37";
            pocketcen_res[nbsr++] = "5.43";
            pocketcen_res[nbsr++] = "6.48";
            pocketcen_res[nbsr++] = "7.38";
        }

        AminoAcid* aa;
        Point pt4avg[nbsr+2];
        
        int j=0;
        for (i=0; i<nbsr; i++)
        {
            aa = p.get_residue(pocketcen_res[i]);
            if (aa)
            {
                aa->priority = priorities[i];
                pt4avg[j++] = aa->get_CA_location();
            }
        }

        pocketcen = average_of_points(pt4avg, j);
        size = size_of_point_space(pt4avg, j);
    }

    loneliest = p.find_loneliest_point(pocketcen, size);
    cout << "Loneliest point = " << loneliest << endl;

    AtomGroup** lagc = AtomGroup::get_potential_ligand_groups(&m, mtlcoords[0].Z > 0);
    agqty = ptrarray_length((void**)lagc);
    if (agqty > MAX_CS_RES-2) agqty = MAX_CS_RES-2;
    for (i=0; i<agqty; i++)
        agc[i] = lagc[i];
    int n = p.get_end_resno();
    Search::prepare_constrained_search(&p, &m, loneliest);
    n = cs_res_qty;
    for (i=0; i<n; i++)
    {
        cout << "Candidate: " << cs_res[i]->get_name() << (cs_res[i]->priority ? "!" : "") << " ~ " << cs_bt[i] << " ~ " << *cs_lag[i] << endl;
    }

    int iter;

    for (iter=1; iter<=10; iter++)
    {
        cout << "Iteration: " << iter << endl << flush;
        Search::do_constrained_search(&p, &m);
        float r = cs_res[cs_idx]->get_CA_location().get_3d_distance(loneliest) - cs_res[cs_idx]->get_reach();
        cout << "Chose: " << cs_res[cs_idx]->get_name() << " ~ " << cs_bt[cs_idx] << " ~ " << *cs_lag[cs_idx] << " r = " << r << endl;

        n = naoi;
        for (i=0; i<n; i++)
        {
            Atom* a = m.get_atom(atoms_of_interest[i].c_str());
            if (!a) cout << "ERROR: No atom named " << atoms_of_interest[i] << " found." << endl;
            else
            {
                AminoAcid** nessami = p.get_residues_near(a->get_location(), _INTERA_R_CUTOFF, false);
                int j, m = ptrarray_length((void**)nessami);
                if (!m) cout << a->name << " void | ";
                else
                {
                    cout << a->name << " near";
                    for (j=0; j<m; j++)
                    {
                        cout << " " << nessami[j]->get_letter() << nessami[j]->get_residue_no();
                    }
                    cout << " | ";
                }
                delete nessami;
            }
        }
        if (n) cout << endl;

        if (save_tmp_pdbs)
        {
            std::string ofn = (std::string)"tmp/cstest" + std::to_string(iter) + (std::string)".pdb";
            fp = fopen(ofn.c_str(), "wb");
            p.save_pdb(fp, &m);
            p.end_pdb(fp);
            fclose(fp);
        }
    }

    return 0;
}
