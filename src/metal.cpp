
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "classes/protein.h"

using namespace std;

int xition_no = 1;

// This is used only for debugging. Once the utility is perfected, it can be removed.
void save_transitional_pdb(Protein* p)
{
    char filename[64];
    sprintf(filename, "temp.%d.pdb", xition_no++);
    FILE* pf = fopen(filename, "wb");
    p->save_pdb(pf);
    p->end_pdb(pf);
    fclose(pf);
    cout << "Wrote " << filename << endl;
}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        cout << "No input file." << endl;
        return -1;
    }

    char* inpfile = NULL;
    char* outfile = NULL;
    int mcoordres[8];
    char* mcoorda[8];
    char* esym = NULL;
    int i, j, charge=2;
    Point pocketcen;
    bool pocketset = false;
    bool dohelix = true;
    char* rcpfile = NULL;
    char* rcpinfo = NULL;

    j=0;
    for (i=1; i<argc; i++)
    {
        if 		(!strcmp(argv[i], "-p")) inpfile = argv[++i];
        else if (!strcmp(argv[i], "-r"))
        {
            mcoordres[j] = atoi(argv[++i]);
            mcoorda[j++] = strstr(argv[i], ":")+1;
        }
        else if (!strcmp(argv[i], "-o")) outfile = argv[++i];
        /*else if (!strcmp(argv[i], "-b"))
        {	pocketcen.x = atof(argv[++i]);
        	pocketcen.y = atof(argv[++i]);
        	pocketcen.z = atof(argv[++i]);
        	pocketset = true;
        }*/
        else if (!strcmp(argv[i], "-c")) charge = atoi(argv[++i]);
        // else if (!strcmp(argv[i], "--nh")) dohelix = false;
        else if (!strcmp(argv[i], "-e")) esym = argv[++i];
        else if (!strcmp(argv[i], "-i"))
        {
            if (argc < i+2 || !argv[i+1]) throw 0xbadf12e;
            rcpfile = argv[i+1];
        }
    }

    int numres = j;
    if (!inpfile
            ||	(	(!numres || !mcoordres[0] || (mcoorda[0] < argv[1]) || !esym )
                    && !rcpfile
               )
            || !outfile)
    {
        cout << "Usage:" << endl;
        cout << "metal -p {path/to/input.pdb} -r {resno}:{aname} -r {resno}:{aname} -r {resno}:{aname} -e {elem_sym} [-b {pocket.x} {pocket.y} {pocket.z}] -o {path/to/output.pdb}" << endl;
        cout << endl << "Example:" << endl;
        cout << "metal -p ";
        cout << (inpfile ? inpfile : "../pdbs/OR1A1/OR1A1.upright.pdb");
        cout << " -r 176:OD1 -r 179:SG -r 180:OD1 -o output/OR1A1.metal.pdb -b -4.93883 5.08067 -4.91533 -e Zn" << endl;
        cout << endl << "Optionally, you may indicate the metal's charge, e.g. -c 1 for a charge of +1. The default charge, if -c is omitted, is +2." << endl;
        cout << endl << "By default, a helix will be created and moved above the binding pocket. If this is not a helix-bound metal ion,"
             << " or if its helix already exists in the protein, use the --nh option to leave the protein structure intact." << endl;
        cout << endl;
        return -1;
    }

    Protein p(inpfile);
    FILE* pf = fopen(inpfile, "r");
    if (!pf)
    {
        cout << "Error trying to read " << inpfile << endl;
        return 0xbadf12e;
    }
    p.load_pdb(pf);
    fclose(pf);

    if (rcpfile)
    {
        pf = fopen(rcpfile, "rb");
        if (!pf) throw 0xbadf12e;
        while (!feof(pf))
        {
            char buffer[16384];
            fgets(buffer, 16361, pf);
            if (buffer[0] >= ' ')
            {
                char** fields = chop_spaced_fields(buffer, '|');

                if (!strcmp(fields[0], "REGION"))
                {
                    p.set_region(fields[1], atoi(fields[2]), atoi(fields[3]));
                }
                else if (!strcmp(fields[0], "METAL") && !j)
                {
                    esym = new char[strlen(fields[1])+4];
                    strcpy(esym, fields[1]);

                    for (i=2; fields[i]; i++)
                    {
                        if (strstr(fields[i], ":"))
                        {
                            mcoordres[j] = atoi(&fields[i][3]);
                            mcoorda[j] = new char[8];
                            strcpy(mcoorda[j], strstr(fields[i], ":")+1);
                            j++;
                            mcoordres[j] = 0;
                            mcoorda[j] = NULL;
                        }
                        else if (strstr(fields[i], "HELIX="))
                        {
                            dohelix = (fields[i][6] - '0');
                        }
                    }
                }
                else if (!strcmp(fields[0], "POCKET") && fields[1])
                {
                    if (!strcmp(fields[1], "type=A"))
                    {
                        char* pktf = strstr(fields[2], "=")+1;
                        pocketcen.x = atof(pktf);
                        pktf = strstr(fields[3], "=")+1;
                        pocketcen.y = atof(pktf);
                        pktf = strstr(fields[4], "=")+1;
                        pocketcen.z = atof(pktf);
                    }
                }
            }
        }
        fclose(pf);
        numres = j;

        for (i=1; i<=7; i++)
        {
            char rgname[10];
            sprintf(rgname, "TMR%d", i);
            Region rgn = p.get_region(rgname);
            float turn = (i & 1) ? -360 : 0;
            if (rgn.start) cout << rgname << " angle: " << (p.get_helix_orientation(rgn.start, rgn.end)*fiftyseven+turn) << endl;
        }
    }




    int startres, hxstart, hxend, midpoint, endres, res;
    float stepht = 6;

    Point mtgt = pocketcen, wayuphigh = pocketcen;
    mtgt.y = 5000;
    wayuphigh.y += 4000;
    AminoAcid* aa;
    Atom *Cend, *Oend;
    Point Cpt, Opt;
    if (dohelix)
    {
        hxstart = mcoordres[0]-1;
        hxend = mcoordres[numres-1]+3;
        Region rgn = p.get_region("TMR4");
        startres = rgn.end+1;
        rgn = p.get_region("TMR5");
        endres = rgn.start-1;

        aa = p.get_residue(endres);
        Cend = aa->get_atom("C");
        Oend = aa->get_atom("O");
        Cpt = Cend->get_location();
        Opt = Oend->get_location();


        //cout << "Extending strand." << endl;
        //p.conform_backbone(startres, endres, Cend, wayuphigh, 50);
        //save_transitional_pdb(&p);


        cout << "Generating an alpha helix." << endl;
        p.make_helix(hxstart, hxend, endres, ALPHA_PHI, ALPHA_PSI);
        save_transitional_pdb(&p);
    }

    cout << "Coordinating metal ion." << endl;
    Atom* m = new Atom(esym, &mtgt);
    m->name = esym;
    strcpy(m->aa3let, "MTL");
    m->increment_charge(charge);
    p.coordinate_metal(m, numres, mcoordres, mcoorda);
    save_transitional_pdb(&p);
    mtgt.y = pocketcen.y + 15;

    if (dohelix)
    {
        cout << "Moving the start of the helix to above TMR5." << endl;
        Point foravg[8];
        int avgi = 0;
        foravg[avgi++] = p.get_atom_location(endres+1, "CA");
        foravg[avgi++] = p.get_atom_location(endres+2, "CA");
        foravg[avgi++] = p.get_atom_location(endres+3, "CA");
        foravg[avgi++] = p.get_atom_location(endres+4, "CA");
        Point pt5 = average_of_points(foravg, avgi);
        pt5.y += stepht*2;

        int end3 = p.get_region_start("TMR3");
        avgi = 0;
        foravg[avgi++] = p.get_atom_location(end3, "CA");
        foravg[avgi++] = p.get_atom_location(end3+1, "CA");
        foravg[avgi++] = p.get_atom_location(end3+2, "CA");
        foravg[avgi++] = p.get_atom_location(end3+3, "CA");
        Point pt3 = average_of_points(foravg, avgi);
        pt3.y += stepht*2;

        avgi = 0;
        foravg[avgi++] = pt5;
        foravg[avgi++] = pt3;
        Point pt35 = average_of_points(foravg, avgi);

        // p.conform_backbone(startres, endres, p.get_atom(hxstart, "CA"), pt5, 100);
        p.move_piece(hxstart, hxend, pt35);
        save_transitional_pdb(&p);

        cout << "Extending strand..." << endl;
        p.conform_backbone(hxend+1, endres, p.get_atom(endres, "C"), wayuphigh, 100);
        save_transitional_pdb(&p);

        cout << "Extending shelf..." << endl;
        Point extreme5 = pt5.subtract(pt3);
        extreme5.scale(100000);
        extreme5 = extreme5.add(pt5);
        p.conform_backbone(hxend+1, endres, p.get_atom(endres, "C"), extreme5, 100);
        save_transitional_pdb(&p);

        // Just from the point found in the previous step, reconnect the loose end to TMR5.
        midpoint = (hxend + endres)/2;
        cout << "Reconnecting broken pieces." << endl;
        p.conform_backbone(midpoint, endres, Cend, Cpt, Oend, Opt, 100);
        save_transitional_pdb(&p);

    }

    pf = fopen(outfile, "wb");

    p.save_pdb(pf);
    p.end_pdb(pf);
    fclose(pf);

    return 0;
}
