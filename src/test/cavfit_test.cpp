
#include "../classes/search.h"
#include "../classes/cavity.h"

int main(int argc, char** argv)
{
    Protein p("Test Receptor");
    Molecule m("Test Ligand");
    Cavity cvtys[256];
    int ncvtys = 0;
    bool priorities[256];

    bool save_tmp_pdbs = false;

    int i, j, n;
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

    float bestc = 0;
    Pose bestp(&m);
    bestp.copy_state(&m);
    for (i=0; i<ncvtys; i++)
    {
        float ctainmt = cvtys[i].find_best_containment(&m);
        if (!i || ctainmt > bestc)
        {
            bestp.copy_state(&m);
            bestc = ctainmt;
        }
    }

    bestp.restore_state(&m);

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

    std::string outf = "cavfit_test.pdb";
    fp = fopen(outf.c_str(), "w");
    if (!fp)
    {
        cout << "FAILED to open output file for writing." << endl;
        return -3;
    }
    p.save_pdb(fp, &m);
    fclose(fp);
    cout << "Saved " << outf << endl;

    return 0;
}
