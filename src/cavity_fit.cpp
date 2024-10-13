
#include "classes/search.h"
#include "classes/cavity.h"

int main(int argc, char** argv)
{
    Protein p("Test Receptor");
    Protein* protein = &p;
    Molecule m("Test Ligand");
    Molecule* ligand = &m;
    Cavity cvtys[256];
    int ncvtys = 0;
    bool priorities[256];
    std::string outfname = "cavfit.pdb";

    bool save_tmp_pdbs = false;

    int i, j, l, n;
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
        else if (dot && !strcasecmp(dot, ".cvty"))
        {
            if (!file_exists(argv[i]))
            {
                cout << "ERROR file not found: " << argv[i] << endl;
                return -7;
            }
            FILE* fp = fopen(argv[i], "rb");
            if (!fp)
            {
                cout << "FAILED to open " << argv[i] << " for reading." << endl;
                return -7;
            }

            cout << "Reading " << argv[i] << "..." << endl;
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
        else if (!strcmp(argv[i], "--bsr"))
        {
            for (i++; argv[i] && argv[i][0] >= '0' && argv[i][0] <= '9'; i++)
            {
                dot = strrchr(argv[i], '.');
                AminoAcid* aa;
                if (dot) aa = p.get_residue_bw(argv[i]);
                else aa = p.get_residue(atoi(argv[i]));
                if (aa) aa->priority = true;
            }
            if (argv[i]) i--;
        }
        else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--out"))
        {
            i++;
            outfname = argv[i];
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

    if (!ncvtys)
    {
        cout << "ERROR: No cavities." << endl;
        return -2;
    }


    std::vector<std::shared_ptr<AtomGroup>> lagc;
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
                }
            }
        }
    }

    size = Point(999,999,999);
    Search::prepare_constrained_search(protein, ligand, Point(0,0,0));

    float bestc = 0;
    int bestl = 0;
    Pose bestp(ligand);
    bestp.copy_state(ligand);
    for (l=0; l<ncvtys; l++)
    {
        cvtys[l].prot = protein;
        float ctainmt = cvtys[l].find_best_containment(ligand, true) * frand(0.5, 1);
        if (!l || ctainmt > bestc)
        {
            bestp.copy_state(ligand);
            bestc = ctainmt;
            bestl = l;
        }
    }
    bestp.restore_state(ligand);

    int csidx = Search::choose_cs_pair(protein, ligand);

    Atom* mtl = (cs_bt[csidx] == mcoord) ? cs_res[csidx]->coordmtl : nullptr;
    ligand->find_mutual_max_bind_potential(cs_res[csidx]);
    if (mtl) ligand->stay_close_other = mtl;

    ligand->movability = MOV_ALL;
    // ligand->enforce_stays();

    for (i=0; i<ncvtys; i++)
    {
        n = cvtys[i].count_partials();
        for (j=0; j<n; j++)
        {
            CPartial* part = cvtys[i].get_partial_by_idx(j);
            sprintf(buffer, "REMARK 821 %4d %8.3f %8.3f %8.3f %7.3f %c%c%c%c%c%c%c %s\n", i, 
                part->s.center.x,
                part->s.center.y,
                part->s.center.z,
                part->s.radius,
                part->metallic ? 'M' : '_',
                part->chargedn ? '-' : '_',
                part->chargedp ? '+' : '_',
                part->polar    ? 'H' : '_',
                part->thio     ? 'S' : '_',
                part->pi       ? 'P' : '_',
                part->priority ? '!' : '_',
                part->resnos_as_string(&p).c_str()
                );
            p.add_remark(buffer);
        }
    }

    fp = fopen(outfname.c_str(), "w");
    if (!fp)
    {
        cout << "FAILED to open " << outfname << " for writing." << endl;
        return -3;
    }
    p.save_pdb(fp, &m);
    fclose(fp);
    cout << "Saved " << outfname << endl;

    return 0;
}
