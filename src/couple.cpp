#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include "classes/protein.h"
#include "classes/group.h"

using namespace std;

Protein *ggpcr, *ggnax;
AminoAcid *pivot1, *pivot2, *pivot3, *pivot4, *pivot5, *pivot6, *pivot7;
Molecule** gcm;

const float montecarlo = fiftyseventh * 5;

void show_usage()
{
    cout << "Usage:" << endl;
    cout << "couple path/to/GPCR.pdb path/to/G-protein.pdb";
    cout << endl;
}

float residue_energy()
{
    int i;
    float e = 0;

    for (i=0; gcm[i]; i++)
    {
        e += gcm[i]->get_intermol_binding(gcm);
    }

    return e;
}

void iteration_callback(int iter)
{
    cout << iter << " " << flush;

    int tmrno = (iter % 7) + 1;
    int sr, er;
    Point pivot, axis(0,0,0);

    char buffer[8];
    sprintf(buffer, "TMR%d", tmrno);
    sr = ggpcr->get_region_start(buffer);
    er = ggpcr->get_region_end(buffer);

    switch (tmrno)
    {
        case 1:        pivot = pivot1->get_CA_location();        break;
        case 2:        pivot = pivot2->get_CA_location();        break;
        case 3:        pivot = pivot3->get_CA_location();        break;
        case 4:        pivot = pivot4->get_CA_location();        break;
        case 5:        pivot = pivot5->get_CA_location();        break;
        case 6:        pivot = pivot6->get_CA_location();        break;
        case 7:        pivot = pivot7->get_CA_location();        break;
    
        default:
        break;
    }

    float e = residue_energy(), e1 = 0, theta;

    int i;

    for (i=0; i<3; i++)
    {
        axis = Point( i==0 ? 1000 : 0, i==1 ? 1000 : 0, i==2 ? 1000 : 0 );

        theta = frand(-montecarlo, montecarlo);
        ggpcr->rotate_piece(sr, er, pivot, axis, theta);
        e1 = residue_energy();
        if (e1 >= e)
        {
            e = e1;
        }
        else
        {
            ggpcr->rotate_piece(sr, er, pivot, axis, -theta);
        }
    }

}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        show_usage();
        return -1;
    }

    Protein gpcr("GPCR");
    Protein gnax("GNAX");
    FILE* fp;

    ggpcr = &gpcr;
    ggnax = &gnax;

    fp = fopen(argv[1], "rb");
    if (!fp)
    {
        cout << "Please ensure " << argv[1] << " exists and is readable." << endl;
        return -1;
    }
    gpcr.load_pdb(fp);
    fclose(fp);
    cout << "GPCR: loaded " << gpcr.get_seq_length() << " residues." << endl;

    fp = fopen(argv[2], "rb");
    if (!fp)
    {
        cout << "Please ensure " << argv[2] << " exists and is readable." << endl;
        return -1;
    }
    gnax.load_pdb(fp);
    fclose(fp);
    cout << "G-protein: loaded " << gnax.get_seq_length() << " residues." << endl;

    int bw1_50 = gpcr.get_bw50(1),
        bw2_50 = gpcr.get_bw50(2),
        bw3_50 = gpcr.get_bw50(3),
        bw4_50 = gpcr.get_bw50(4),
        bw45_50 = gpcr.get_bw50(45),
        bw5_50 = gpcr.get_bw50(5),
        bw6_50 = gpcr.get_bw50(6),
        bw7_50 = gpcr.get_bw50(7);

    int c438, c359, c756;

    if (bw3_50 > 0) c359 = bw3_50 + 9;
    else
    {
        int q = gpcr.search_sequence(100, 200, "DRYXAICXPLXY");
        if (q < 1)
        {
            cout << "Cannot find 3.59. Coupling fail." << endl;
            return -1;
        }
        else c359 = q + 10;
    }

    if (bw4_50 > 0) c438 = bw4_50 - 12;
    else c438 = c359 + 6;

    if (bw7_50 > 0) c756 = bw7_50 + 6;
    else
    {
        int q = gpcr.search_sequence(250, 400, "PMLNPLIYSLRNKD");
        if (q < 1)
        {
            cout << "Cannot find 7.56. Coupling fail." << endl;
            return -1;
        }
        else c756 = q + 10;
    }

    cout << "GPCR contact residues: " << c359 << " " << c438 << " " << c756 << endl;

    int q35 = gnax.search_sequence(1, 1000, "EKQLQKD", 6);
    int d215 = gnax.search_sequence(120, 1000, "GIFETKFQVD");
    int e392 = gnax.search_sequence(gnax.get_end_resno()-10, 1000, "YELL");

    if (q35 < 1 || d215 < 1 || e392 < 1)
    {
        cout << "Cannot find G-protein motifs. Coupling fail." << endl;                 // My whole adult life has been coupling fails.
        return -1;
    }

    q35 += 8;
    d215 += 9;
    e392 += 1;

    cout << "G-protein contact residues: " << q35 << " " << d215 << " " << e392 << endl;


    // Upright both proteins.
    cout << "Turning GPCR upright..." << endl;
    gpcr.upright();

    cout << "Turning G-protein upright..." << endl;
    Point pt35x4(0,0,0);
    int i, j=0;

    for (i=q35; i>1 && j<4; i--)
    {
        AminoAcid* aa = gnax.get_residue(i);
        if (aa)
        {
            pt35x4 = pt35x4.add(aa->get_CA_location());
            j++;
        }
    }

    if (j)
    {
        pt35x4.x /= j;
        pt35x4.y /= j;
        pt35x4.z /= j;
    }

    Point pt1x4(0,0,0);

    j=0;
    for (i=1; i<q35 && j<4; i++)
    {
        AminoAcid* aa = gnax.get_residue(i);
        if (aa)
        {
            pt1x4 = pt1x4.add(aa->get_CA_location());
            j++;
        }
    }

    if (j)
    {
        pt1x4.x /= j;
        pt1x4.y /= j;
        pt1x4.z /= j;
    }

    Point axis = pt1x4.subtract(pt35x4);
    Rotation rot = align_points_3d(axis, Point(-1000,0,0), Point(0,0,0));
    gnax.rotate_piece(1, 99999, rot, e392);

    axis = gnax.get_residue(e392)->get_CA_location().subtract(pt35x4);
    axis.x = 0;
    rot = align_points_3d(axis, Point(0,1000,0), Point(0,0,0));
    gnax.rotate_piece(1, 99999, rot, e392);

    rot.v = (SCoord)Point(-1000,0,0);
    rot.a = fiftyseventh*54;
    gnax.rotate_piece(1, 99999, rot, e392);

    gnax.move_piece(1, 99999, Point(0, -50, 0));
    gpcr.move_piece(1, 99999, Point(0, 0, 0));



    // Move TMR6 out of the way sooner than later.
    cout << "Moving TMR6..." << endl;

    AminoAcid* aa;

    AminoAcid *sbb = nullptr, *sba = nullptr;
    for (i=6; i<=11; i++)
    {
        aa = gpcr.get_residue(bw6_50+i);
        if (aa)
        {
            char c = aa->get_letter();
            if (c == 'R' || c == 'K')
            {
                AminoAcid* aa1 = gpcr.get_residue(bw45_50 + 1);
                if (aa1)
                {
                    c = aa1->get_letter();
                    if (c == 'D' || c == 'E')
                    {
                        sbb = aa;
                        sba = aa1;
                        cout << "Found potential salt bridge between " << *aa << " and " << *aa1 << endl;
                        break;
                    }
                }
            }
        }
    }

    Point ptdest, ptsrc;
    Point ptcen(0,0,0);

    pivot6 = nullptr;

    float salt_reach = 0;
    float salt_angle = M_PI;
    if (sba && sbb)
    {
        salt_angle = find_3d_angle(sbb->get_atom("CB")->get_location(), sba->get_CA_location(), sbb->get_CA_location());
        salt_reach = sba->get_reach() + sbb->get_reach() * (0.5 + 0.5 * cos(salt_angle)) + 2;
    }

    if (sba && sbb)
    {
        ptsrc = sbb->get_CA_location();
        float salt_gap = sbb->get_CA_location().get_3d_distance(sba->get_CA_location());
        float salt_delta = fmax(salt_gap - salt_reach, 0);
        // cout << "Moving " << *sbb << " " << salt_delta << "A." << endl;
        ptdest = sba->get_CA_location();
        ptdest = ptdest.multiply_3d_distance(&ptsrc, salt_delta/salt_gap);
        // cout << "Acid at " << sba->get_CA_location() << " so moving base from " << ptsrc << " to " << ptdest << endl;
        gpcr.rotate_piece(gpcr.get_region_start("TMR6"), gpcr.get_region_end("TMR6"), sbb->get_residue_no(), ptdest, bw6_50-2);

        Molecule* tobridge[4];
        tobridge[0] = (Molecule*)sba;
        tobridge[1] = (Molecule*)sbb;
        tobridge[2] = nullptr;

        cout << "Forming salt bridge..." << endl;
        Molecule::conform_molecules(tobridge, 50);

        pivot6 = sbb;
    }
    else
    {
        aa = gpcr.get_residue(bw5_50);
        if (!aa) cout << "No BW5.50 residue." << endl;
        axis = aa->get_CA_location();

        aa = gpcr.get_residue(bw7_50-5);
        if (!aa) cout << "No BW7.45 residue." << endl;
        axis = axis.subtract(aa->get_CA_location());

        rot.v = axis;
        rot.a = fiftyseventh * 15;
        gpcr.rotate_piece(gpcr.get_region_start("TMR6"), gpcr.get_region_end("TMR6"), rot, bw6_50-2);

        pivot6 = gpcr.get_residue(bw6_50-2);
    }


    // Line up CYT3 to at least point to TMR6.
    cout << "Moving CYT3..." << endl;
    gpcr.rotate_piece(gpcr.get_region_end("TMR5")+1, gpcr.get_region_start("TMR6")-1, gpcr.get_region_start("TMR6")-1,
        gpcr.get_residue(gpcr.get_region_start("TMR6"))->get_CA_location(), gpcr.get_region_end("TMR5"));


    // Superimpose YELL with 7.56.
    cout << "Aligning C-terminus of G-protein..." << endl;

    ptdest = gpcr.get_residue(c756)->get_CA_location();
    ptsrc  = gnax.get_residue(e392)->get_CA_location();

    // Move YELL partway towards the midpoint of 5.64 and 6.33.
    j = 4;
    if (bw5_50 > 0 && bw6_50 > 0)
    {
        ptdest = ptdest.multiply_3d_distance(&ptcen, 4);

        aa = gpcr.get_residue(bw5_50 + 14);
        if (aa)
        {
            ptdest = ptdest.add(aa->get_CA_location());
            j++;
        }

        aa = gpcr.get_residue(bw6_50 - 17);
        if (aa)
        {
            ptdest = ptdest.add(aa->get_CA_location());
            j++;
        }

        ptdest = ptdest.multiply_3d_distance(&ptcen, 1.0/j);
    }

    gnax.move_piece(1, 99999, (SCoord)ptdest.subtract(ptsrc));


    // Rotate to align Q35 with *4.38.
    cout << "Aligning rotation of G-protein..." << endl;
    ptdest = gpcr.get_residue(c438)->get_CA_location();
    ptsrc  = gnax.get_residue(q35)->get_CA_location();
    ptcen = gnax.get_residue(e392)->get_CA_location();
    ptsrc.y = ptdest.y = ptcen.y;
    rot = align_points_3d(ptsrc, ptdest, ptcen);
    gnax.rotate_piece(1, 99999, rot, e392);



    // Do the positional fine tuning.
    gpcr.set_pdb_chain('A');
    gnax.set_pdb_chain('B');

    std::vector<AminoAcid*> cr = gpcr.get_contact_residues(&gnax);

    int n = cr.size();
    cout << "Contact residues: ";
    for (i=0; i<n; i++) cout << *cr[i] << " ";
    cout << endl;

    pivot1 = gpcr.get_residue(bw1_50 - 9);
    pivot2 = gpcr.get_residue(bw2_50 + 7);
    pivot3 = gpcr.get_residue(bw3_50 - 14);
    pivot4 = gpcr.get_residue(bw4_50);
    pivot5 = gpcr.get_residue(bw5_50);
    pivot7 = gpcr.get_residue(bw7_50 - 8);

    Molecule* contact_mols[n+4];

    for (i=0; i<n; i++)
    {
        contact_mols[i] = (Molecule*)cr[i];
    }
    contact_mols[i] = nullptr;
    gcm = contact_mols;

    cout << "Fine tuning positions..." << endl;
    Molecule::conform_molecules(contact_mols, 200, &iteration_callback);
    cout << endl;



    // Go ahead and write the output file now, to see where we're at with development.
    fp = fopen("output/coupled.pdb", "wb");
    if (!fp)
    {
        cout << "Failed to open output file." << endl;
        return -1;
    }

    gpcr.set_pdb_chain('A');
    gpcr.save_pdb(fp);
    gnax.set_pdb_chain('B');
    gnax.save_pdb(fp);
    gnax.end_pdb(fp);
    fclose(fp);
    cout << "Wrote output file." << endl;
}