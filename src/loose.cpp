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
#include <sys/stat.h>
#include "classes/protein.h"

using namespace std;

int iters = 50;

// https://stackoverflow.com/a/6039648/18753403
long GetFileSize(std::string filename)
{
    struct stat stat_buf;
    int rc = stat(filename.c_str(), &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}

void iteration_callback(int iter, Molecule** mols)
{
    float progress = (float)iter / iters;
    float percentage = progress * 100;

    int i;
    Point lbc = mols[0]->get_barycenter();
    for (i=1; mols[i]; i++)
    {
        Point pt = mols[i]->get_barycenter();
        SCoord v = lbc.subtract(pt);
        v.r = 0.333;
        mols[i]->move(v);
    }

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

int main(int argc, char** argv)
{
    std::string dock_fname, pdb_fname, sdf_fname, out_fname;
    std::vector<int> docked_resnos;
    std::vector<Point> orig_CA_locs, new_CA_locs;
    Protein p("TheReceptor");
    Molecule lig("TheLigand");

    int i, j, l, m, n, pose = 0;
    out_fname = "tmp/loose.pdb";

    for (i=1; i<argc; i++)
    {
        if (argv[i][0] == '-')
        {
            if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--iters")) iters = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--out")) out_fname = argv[++i];
            else
            {
                cout << "Unrecognized option " << argv[i] << endl;
                return -1;
            }
        }
        else
        {
            dock_fname = argv[i];
        }
    }

    FILE* fp = fopen(dock_fname.c_str(), "r");
    char buffer[1024];
    while (fgets(buffer, 1020, fp))
    {
        for (j=0; buffer[j]; j++) if (buffer[j] < ' ') buffer[j] = 0;
        std::string ln = buffer;
        if (!strcmp(ln.substr(0, 10).c_str(), "PDB file: "))
        {
            pdb_fname = ln.substr(10);
            FILE* fpdb = fopen(pdb_fname.c_str(), "r");
            p.load_pdb(fpdb);
            fclose(fpdb);
            cout << "Loaded " << pdb_fname << endl;

            for (j=1; j<=7; j++)
            {
                std::string regno = (std::string)"TMR" + std::to_string(j);
                int sr = p.get_region_start(regno), er = p.get_region_end(regno);
                if (sr && er)
                {
                    cout << regno << " begin " << sr << " @ " << p.get_atom_location(sr, "CA");
                    int bw50 = p.get_bw50(j);
                    cout << " n.50 " << bw50 << " @ " << p.get_atom_location(bw50, "CA");
                    cout << " end " << er << " @ " << p.get_atom_location(er, "CA");
                    cout << endl;
                }
            }

            cout << endl << flush;
            continue;
        }

        if (!strcmp(ln.substr(0, 8).c_str(), "Ligand: "))
        {
            sdf_fname = ln.substr(8);
            n = GetFileSize(sdf_fname);
            char sdfdat[n + 4];
            FILE* fsdf = fopen(sdf_fname.c_str(), "r");
            fread(sdfdat, 1, n, fsdf);
            lig.from_sdf(sdfdat);
            fclose(fsdf);
            cout << "Loaded " << sdf_fname << endl << flush;
            continue;
        }

        if (!strcmp(ln.substr(0, 6).c_str(), "Pose: "))
        {
            pose = atoi(ln.substr(6).c_str());
            if (pose > 1) break;                    // for now.
            continue;
        }

        if (!strcmp(ln.substr(0, 7).c_str(), "HETATM "))
        {
            char** words = chop_spaced_words(buffer);

            Atom* a = lig.get_atom(words[2]);
            if (!a) continue;
            float x = atof(words[5]), y = atof(words[6]), z = atof(words[7]);
            a->move(Point(x,y,z));

            delete[] words;
        }

        if (!strcmp(ln.substr(0, 7).c_str(), "ATOM   "))
        {
            char** words = chop_spaced_words(buffer);

            int resno = atoi(words[4]);
            AminoAcid* aa = p.get_residue(resno);
            if (!aa) continue;

            if (!strcmp(words[2], "CA"))
            {
                if (std::find(docked_resnos.begin(), docked_resnos.end(), resno) == docked_resnos.end())
                {
                    docked_resnos.push_back(resno);
                    orig_CA_locs.push_back(aa->get_CA_location());
                    new_CA_locs.push_back(aa->get_CA_location());
                }
            }

            Atom* a = aa->get_atom(words[2]);
            if (!a) continue;
            float x = atof(words[5]), y = atof(words[6]), z = atof(words[7]);
            a->move(Point(x,y,z));

            delete[] words;
        }
    }
    fclose(fp);
    cout << "Dock file loaded." << endl;

    n = docked_resnos.size();
    Molecule* mols[n + 4];
    j = 0;
    lig.movability = MOV_PINNED;
    mols[j++] = &lig;
    for (i=0; i<n; i++)
    {
        AminoAcid* aa = p.get_residue(docked_resnos[i]);
        if (!aa) continue;
        aa->movability = MOV_ALL;
        mols[j++] = aa;
    }
    n = j;
    mols[n] = nullptr;

    Molecule::conform_molecules(mols, iters, &iteration_callback);

    m = docked_resnos.size();
    for (i=0; i<n; i++)
    {
        l = mols[i]->is_residue();
        if (!l) continue;
        AminoAcid* aa = p.get_residue(l);
        if (!aa) continue;

        for (j=0; j<m; j++)
        {
            if (docked_resnos[j] == l)
            {
                Point newloc = aa->get_CA_location();
                new_CA_locs[j] = newloc;

                Point delta = newloc.subtract(orig_CA_locs[j]);
                if (delta.magnitude() < 0.1) continue;
                cout << l << ": " << delta;

                Point lbc = lig.get_barycenter();
                lbc.y = newloc.y = orig_CA_locs[j].y;

                float theta = find_3d_angle(newloc, lbc, orig_CA_locs[j]);
                cout << " ligand angle " << (theta*fiftyseven) << "deg.";

                cout << endl;
                break;
            }
        }
    }

    const Region* reg = p.get_regions();
    m = docked_resnos.size();
    for (i=0; reg[i].start && reg[i].end; i++)
    {
        l = 0;
        Point pt, rgcen;
        for (j=0; j<m; j++)
        {
            if (docked_resnos[j] >= reg[i].start && docked_resnos[j] <= reg[i].end)
            {
                rgcen = rgcen.add(new_CA_locs[j]);
                Point delta = new_CA_locs[j].subtract(orig_CA_locs[j]);
                pt = pt.add(delta);
                l++;
            }
        }

        if (!l) continue;
        pt.multiply(1.0/l);
        rgcen.multiply(1.0/l);

        float rcm = p.region_can_move(reg[i].start, reg[i].end, pt);
        if (rcm > pt.magnitude()) rcm = pt.magnitude();
        else pt.scale(rcm);
        if (rcm > 0)
        {
            p.move_piece(reg[i].start, reg[i].end, (SCoord)pt);
            for (j=0; j<m; j++)
            {
                if (docked_resnos[j] >= reg[i].start && docked_resnos[j] <= reg[i].end)
                {
                    orig_CA_locs[j] = orig_CA_locs[j].add(pt);
                }
            }
        }

        int hxno = atoi(((std::string)reg[i].name).substr(3).c_str());
        if (hxno < 1 || hxno > 7) continue;

        int bw50 = p.get_bw50(hxno);

        SCoord axis;
        float theta = 0;
        Point topdir, oldtop;
        l = 0;
        for (j=0; j<m; j++)
        {
            if (abs(docked_resnos[j] - bw50) < 3) continue;
            if (docked_resnos[j] >= reg[i].start && docked_resnos[j] <= reg[i].end)
            {
                if (!l || new_CA_locs[j].y > oldtop.y) oldtop = new_CA_locs[j];

                Point delta = new_CA_locs[j].subtract(rgcen);
                delta.multiply(1.0 / delta.y);
                delta.y = 0;
                topdir = topdir.add(delta);
                l++;
            }
        }

        if (!l) continue;
        topdir.multiply(1.0/l);
        topdir = topdir.add(oldtop);

        axis = compute_normal(oldtop, topdir, rgcen);
        theta = fmin(find_angle_along_vector(oldtop, topdir, rgcen, axis)
            , p.region_can_rotate(reg[i].start, reg[i].end, axis)
            );

        if (theta) p.rotate_piece(reg[i].start, reg[i].end, rgcen, axis, theta);
    }

    fp = fopen(out_fname.c_str(), "w");
    p.save_pdb(fp, &lig);
    p.end_pdb(fp);
    fclose(fp);

    return 0;
}
