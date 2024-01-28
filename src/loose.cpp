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
    Pose putitback;

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

    putitback.copy_state(&lig);

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


    Protein p1;
    FILE* fpdb = fopen(pdb_fname.c_str(), "r");
    p1.load_pdb(fpdb);
    fclose(fpdb);

    putitback.restore_state(&lig);


    const Region* reg = p1.get_regions();
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

        float rcm = p1.region_can_move(reg[i].start, reg[i].end, pt);
        if (rcm > pt.magnitude()) rcm = pt.magnitude();
        else pt.scale(rcm);
        if (rcm > 0)
        {
            p1.move_piece(reg[i].start, reg[i].end, (SCoord)pt);
            for (j=0; j<m; j++)
            {
                if (docked_resnos[j] >= reg[i].start && docked_resnos[j] <= reg[i].end)
                {
                    orig_CA_locs[j] = orig_CA_locs[j].add(pt);
                }
            }
        }

        // continue;               // Debug moving regions without rotating.

        int hxno = atoi(((std::string)reg[i].name).substr(3).c_str());
        if (hxno < 1 || hxno > 7) continue;

        int bw50 = p1.get_bw50(hxno);

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
            , p1.region_can_rotate(reg[i].start, reg[i].end, axis)
            );

        if (theta) p1.rotate_piece(reg[i].start, reg[i].end, rgcen, axis, theta);
    }


    AminoAcid *aa4x60 = p1.get_residue_bw("4.60");
    AminoAcid *aa45x51 = p1.get_residue_bw("45.51");
    AminoAcid *aa45x53 = p1.get_residue_bw("45.53");
    AminoAcid *aa5x58 = p1.get_residue_bw("5.58");
    AminoAcid *aa6x48 = p1.get_residue_bw("6.48");
    AminoAcid *aa6x49 = p1.get_residue_bw("6.49");
    AminoAcid *aa6x55 = p1.get_residue_bw("6.55");
    AminoAcid *aa6x59 = p1.get_residue_bw("6.59");
    AminoAcid *aa7x53 = p1.get_residue_bw("7.53");

    int n4x60 = aa4x60->get_residue_no();
    int n45x51 = aa45x51->get_residue_no();
    int n45x53 = aa45x53->get_residue_no();
    int n5x58 = aa5x58->get_residue_no();
    int n6x48 = aa6x48->get_residue_no();
    int n6x49 = aa6x49->get_residue_no();
    int n6x55 = aa6x55->get_residue_no();
    int n6x59 = aa6x59->get_residue_no();
    int n7x53 = aa7x53->get_residue_no();

    char l45x51 = aa45x51->get_letter();
    char l45x53 = aa45x53->get_letter();
    char l5x58 = aa5x58->get_letter();
    char l6x48 = aa6x48->get_letter();
    char l6x49 = aa6x49->get_letter();
    char l6x55 = aa6x55->get_letter();
    char l6x59 = aa6x59->get_letter();
    char l7x53 = aa7x53->get_letter();

    int npivot = n6x49;
    if (l6x59 == 'R')
    {
        Region rg = p1.get_region("TMR6");
        Point pt = aa4x60->get_CA_location();
        Point origin = aa6x48->get_CA_location();
        SCoord axis = compute_normal(aa6x59->get_CA_location(), pt, origin);
        float theta = fmin(find_angle_along_vector(aa6x59->get_CA_location(), pt, origin, axis)/3.5
            , p1.region_can_rotate(rg.start, rg.end, axis)
            );
        if (theta)
        {
            cout << "R6.55 Rock6 " << (theta*fiftyseven) << "deg." << endl;
            p1.rotate_piece(rg.start, rg.end, origin, axis, theta);
        }
        aa6x59->conform_atom_to_location("NE", pt);
    }
    else if ((l45x51 == 'D' || l45x51 == 'E') && l6x55 == 'Y')
    {
        Region rg = p1.get_region("TMR6");
        Point pt = aa45x51->get_CA_location().subtract(aa6x55->get_CA_location());
        Atom* a51 = aa45x51->get_nearest_atom(aa6x55->get_CA_location());
        if (a51)
        {
            Atom* a55 = aa6x55->get_nearest_atom(a51->get_location());
            if (a55)
            {
                float r = a51->distance_to(a55);
                if (r > 2.0)
                {
                    pt.scale(r-2.0);
                    pt = pt.add(aa6x55->get_CA_location());
                    Point origin = aa6x48->get_CA_location();
                    SCoord axis = compute_normal(aa6x55->get_CA_location(), pt, origin);
                    float theta = fmin(find_angle_along_vector(aa6x55->get_CA_location(), pt, origin, axis)
                        , p1.region_can_rotate(rg.start, rg.end, axis)
                        );
                    if (theta)
                    {
                        cout << "Y~D/E Rock6 " << (theta*fiftyseven) << "deg." << endl;
                        p1.rotate_piece(rg.start, rg.end, origin, axis, theta);
                    }
                    p1.bridge(n45x51, n6x55);
                }
            }
        }
    }
    else if ((l45x53 == 'N' || l45x53 == 'Q') && (l6x55 == 'D' || l6x55 == 'E'))
    {
        Region rg = p1.get_region("TMR6");
        Point pt = aa45x53->get_CA_location().subtract(aa6x55->get_CA_location());
        Atom* a53 = aa45x53->get_nearest_atom(aa6x55->get_CA_location());
        if (a53)
        {
            Atom* a55 = aa6x55->get_nearest_atom(a53->get_location());
            if (a55)
            {
                float r = a53->distance_to(a55);
                if (r > 2.0)
                {
                    pt.scale(r-2.0);
                    pt = pt.add(aa6x55->get_CA_location());
                    Point origin = aa6x48->get_CA_location();
                    SCoord axis = compute_normal(aa6x55->get_CA_location(), pt, origin);
                    float theta = fmin(find_angle_along_vector(aa6x55->get_CA_location(), pt, origin, axis)
                        , p1.region_can_rotate(rg.start, rg.end, axis)
                        );
                    if (theta)
                    {
                        cout << "D/E~N/Q Rock6 " << (theta*fiftyseven) << "deg." << endl;
                        p1.rotate_piece(rg.start, rg.end, origin, axis, theta);
                    }
                    p1.bridge(n45x53, n6x55);
                }
            }
        }
    }


    // TMR6 motion.
    Region rg5 = p1.get_region("TMR5");
    Region rg6 = p1.get_region("TMR6");
    float theta6_phi_initial = 10;
    LocatedVector lv6 = aa6x49->get_phi_vector();
    p1.rotate_piece(rg6.start, n6x49, lv6.origin, lv6, -fiftyseventh*theta6_phi_initial);
    float theta6_psi_initial = 40;
    lv6 = aa6x49->get_psi_vector();
    p1.rotate_piece(rg6.start, n6x49, lv6.origin, lv6, -fiftyseventh*theta6_psi_initial);

    // 5.58 ~ 7.53 contact.
    p1.bridge(n5x58, n7x53);

    int n5x50 = p1.get_bw50(5);
    AminoAcid* aa5x33 = p1.get_residue(rg5.start), *aa5x68 = p1.get_residue(rg5.end), *aa5x50 = p1.get_residue(n5x50);
    AminoAcid* aa6x28 = p1.get_residue(rg6.start);
    LocatedVector lv5 = compute_normal(aa5x68->get_CA_location(), aa6x28->get_CA_location(), aa5x50->get_CA_location());
    lv5.origin = aa5x50->get_CA_location();

    float theta5 = fmin(15*fiftyseven, p1.region_can_rotate(n5x50, rg5.end, lv5));
    p1.rotate_piece(n5x50, rg5.end, lv5.origin, lv5, theta5);

    lv5 = compute_normal(aa5x58->get_CA_location(), aa7x53->get_CA_location(), aa5x50->get_CA_location());
    theta5 = p1.region_can_rotate(n5x50, rg5.end, lv5, false, 0, rg6.start, n6x59);
    Point pt;
    float r = p1.contact_dist(n5x58, n7x53);
    r = fmax(0, r - contact_r_5x58_7x53);
    if (r)
    {
        pt = aa7x53->get_reach_atom()->get_location().subtract(aa5x58->get_reach_atom()->get_location());
        Point pt5x58 = aa5x58->get_reach_atom()->get_location();
        float f = find_3d_angle(pt5x58, pt5x58.add(pt), lv5.origin);
        if (f < theta5) theta5 = f;

        p1.rotate_piece(n5x50, rg5.end, lv5.origin, lv5, theta5);
        cout << "TMR5 bends " << theta5*fiftyseven << "deg, limited by " << *(p1.stop1) << "->" << *(p1.stop2) << endl;
    }
    else cout << "5.58~7.53 bridge already met." << endl;

    // Adjust TMR6.
    AminoAcid* aapivot = p1.get_residue(npivot);
    float theta6 = p1.region_can_rotate(rg6.start, npivot, lv6);
    Point pt6x28 = rotate3D(aa6x28->get_CA_location(), aapivot->get_CA_location(), lv6, theta6);
    float r56 = pt6x28.get_3d_distance(aa5x68->get_CA_location());
    cout << "r56 = " << r56 << endl;
    if (r56 > 4)
    {
        theta6 = fmin(theta6, find_3d_angle(aa5x68->get_CA_location(), aa6x28->get_CA_location(), aapivot->get_CA_location()));
        p1.stop1 = aa5x68;
        p1.stop2 = aa6x28;
    }
    cout << "TMR6 bends " << (theta6_phi_initial + theta6_psi_initial - theta6 * fiftyseven) << "deg limited by " << *(p1.stop1) << "->" << *(p1.stop2) << endl;
    p1.rotate_piece(rg6.start, npivot, lv6.origin, lv6, theta6);

    p1.bridge(n5x58, n7x53);
    r = p1.contact_dist(n5x58, n7x53);
    cout << "5~7 contact distance: " << r << "A." << endl;


    fp = fopen(out_fname.c_str(), "w");
    p1.save_pdb(fp);
    p1.end_pdb(fp);
    fclose(fp);
    cout << "Saved " << out_fname << endl;

    return 0;
}
