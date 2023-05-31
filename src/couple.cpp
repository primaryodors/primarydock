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

void show_usage()
{
    cout << "Usage:" << endl;
    cout << "couple path/to/GPCR.pdb path/to/G-protein.pdb";
    cout << endl;
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

    int bw3_50 = gpcr.get_bw50(3),
        bw4_50 = gpcr.get_bw50(4),
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


    // Superimpose YELL with 7.56.
    cout << "Aligning C-terminus of G-protein..." << endl;

    Point ptdest = gpcr.get_residue(c756)->get_CA_location();
    Point ptsrc  = gnax.get_residue(e392)->get_CA_location();
    Point ptcen(0,0,0);

    // Move YELL halfway to the midpoint of 5.64 and 6.33.
    AminoAcid* aa;
    j = 2;
    if (bw5_50 > 0 && bw6_50 > 0)
    {
        ptdest = ptdest.multiply_3d_distance(&ptcen, 2);

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



    // Go ahead and write the output file now, to see where we're at with development.
    fp = fopen("output/coupled.pdb", "wb");
    if (!fp)
    {
        cout << "Failed to open output file." << endl;
        return -1;
    }
    gpcr.save_pdb(fp);
    gnax.save_pdb(fp);
    gnax.end_pdb(fp);
    fclose(fp);
    cout << "Wrote output file." << endl;
}