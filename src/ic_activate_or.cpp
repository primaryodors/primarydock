#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <unistd.h>
#include <fstream>
#include <sys/stat.h>
#include <sys/wait.h>
#include "classes/protein.h"
#include "classes/dynamic.h"

using namespace std;

bool file_exists(std::string fname)
{
    struct stat s;
    if (stat(fname.c_str(), &s) == 0) return true;
    else return false;
}

int main(int argc, char** argv)
{
    if (argc < 2) return -1;
    int i, j, l, m, n;

    // One argument, in the form of a protein ID, e.g. OR51E2.
    // If the ID does not conform to the format of OR{family}{subfamily}{member}, return an error code.
    if (argv[1][0] != 'O' || argv[1][1] != 'R') return -2;
    l = 2;
    std::string tmp = argv[1];
    int fam = atoi(tmp.substr(l, 2).c_str());
    if (!fam) return -2;
    l++;
    if (fam >= 10) l++;
    std::string sub = tmp.substr(l, 1);
    l++;
    if (argv[1][l] >= 'A' && argv[1][l] <= 'Z')
    {
        sub = tmp.substr(l-1, 2);
        l++;
    }
    int mem = atoi(&argv[1][l]);
    if (!mem) return -2;


    // Locate and load the .upright.pdb in the folder structure, or throw an error if not found.
    std::string path = (std::string)"pdbs/OR" + std::to_string(fam) + (std::string)"/";
    std::string orid = (std::string)"OR" + std::to_string(fam) + sub + std::to_string(mem);
    std::string in_filename = path + orid + (std::string)".upright.pdb";
    if (!file_exists(in_filename)) return -3;

    Protein p(orid.c_str());
    FILE* fp = fopen(in_filename.c_str(), "rb");
    if (!fp) return -3;
    p.load_pdb(fp);
    fclose(fp);


    // If there is an R at position 6.59, perform a slight helix unwind and bring 6.59's side chain inward.
    DynamicMotion unwind6(&p);
    DynamicMotion bend45(&p);
    AminoAcid* aa6x59 = p.get_residue_bw("6.59");
    AminoAcid* aa6x58 = p.get_residue_bw("6.58");
    AminoAcid* aa4x60 = p.get_residue_bw("4.60");
    AminoAcid* aa45x53 = p.get_residue_bw("45.53");
    if (aa6x59 && aa6x58 && aa45x53 && (aa6x59->get_letter() == 'R'))
    {
        // Have to move 45.53 and 45.54 out of the way, but leave 45.52 where it is.
        bend45.start_resno.from_string("45.53");
        bend45.end_resno.from_string("45.54");
        bend45.type = dyn_bend;
        bend45.fulcrum_resno.from_string("45.53");
        bend45.axis_resno.from_string("45.54");
        bend45.bias = 30;
        bend45.apply_absolute(1);
        aa45x53->conform_atom_to_location(aa45x53->get_reach_atom()->name, Point(0,-1000,0));

        cout << "R6.59" << endl;
        unwind6.start_resno.from_string("6.58");
        unwind6.end_resno.from_string("6.59");
        unwind6.type = dyn_wind;
        unwind6.bias = -13;
        unwind6.apply_absolute(1);

        aa6x59->conform_atom_to_location("CZ", aa4x60->get_CA_location());
    }


    float TMR6ex = 0, TMR6ey;
    AminoAcid *aa6x48 = p.get_residue_bw("6.48");
    AminoAcid *aa6x55 = p.get_residue_bw("6.55");
    AminoAcid *aa45x51 = p.get_residue_bw("45.51");
    AminoAcid *aa5x39 = p.get_residue_bw("5.39");
    char l45x51 = aa45x51->get_letter();
    char l6x58 = aa6x58->get_letter();
    char l45x53 = aa45x53->get_letter();
    SCoord axis6;

    // If S6.58 and P45.53 and D/E45.51, measure how far S6.58 would have to move to make contact with 45.51.
    if ((l6x58 == 'S' || l6x58 == 'T')
        &&
        (l45x53 == 'P' || l45x53 == 'A' || l45x53 == 'G' || l45x53 == 'C' || l45x53 == 'S')
        &&
        (l45x51 == 'D' || l45x51 == 'E')
        )
    {
        axis6 = compute_normal(aa6x48->get_CA_location(), aa6x58->get_CA_location(), aa45x51->get_CA_location());
        p.bridge(aa6x58->get_residue_no(), aa45x51->get_residue_no());
        Atom* a = aa6x58->get_reach_atom();
        if (a->get_Z() < 3) a = a->get_bond_by_idx(0)->btom;
        float cr = a->distance_to(aa45x51->get_nearest_atom(a->get_location()));
        TMR6ex = fmax(0, cr - 3);
        TMR6ey = aa6x48->distance_to(aa6x58);
    }

    // If there is no Y6.55 or no D/E45.51, measure how far 6.58-6.59 can move toward 5.39 without clashing. Call it TMR6ex.
    else if (aa6x55->get_letter() != 'Y' || (l45x51 != 'D' && l45x51 != 'E'))
    {
        SCoord TMR6edir = aa5x39->get_CA_location().subtract(p.get_region_center(aa6x58->get_residue_no(), aa6x59->get_residue_no()));
        TMR6ex = p.region_can_move(aa6x55->get_residue_no(), aa6x58->get_residue_no(), TMR6edir, true);
        TMR6ey = aa6x48->distance_to(aa6x59);
        cout << "Rock6 can move " << TMR6ex << "A in the EXR domain." << endl;

        axis6 = compute_normal(aa6x48->get_CA_location(), aa6x59->get_CA_location(), aa5x39->get_CA_location());
    }

    // If there is Y6.55 and D/E45.51 but they are not in contact, measure how far 6.55 must move toward 45.51 to make contact. Call it TMR6ex.
    else
    {
        axis6 = compute_normal(aa6x48->get_CA_location(), aa6x55->get_CA_location(), aa45x51->get_CA_location());
        p.bridge(aa6x55->get_residue_no(), aa45x51->get_residue_no());
        float cr = aa6x55->get_atom("OH")->distance_to(aa45x51->get_nearest_atom(aa6x55->get_atom("OH")->get_location()));
        if (cr > 3)
        {
            TMR6ex = cr - 3;
            /*SCoord TMR6edir = aa5x39->get_CA_location().subtract(p.get_region_center(aa6x58->get_residue_no(), aa6x59->get_residue_no()));
            TMR6ex = min(TMR6ex, p.region_can_move(aa6x55->get_residue_no(), aa6x58->get_residue_no(), TMR6edir, true));*/
            TMR6ey = aa6x55->distance_to(aa6x48);
            cout << "Hybrid6 can move " << TMR6ex << "A in the EXR domain." << endl;
        }

        // Otherwise set TMR6ex to zero.
        else
        {
            TMR6ex = 0;
            cout << "Bend6 is stationary in the EXR domain." << endl;
        }
    }


    // Measure how far 5.43-5.54 can move toward 6.44-6.55 without clashing. Call it TMR5ez.
    // TODO: There's also a slight movement of TMR3 that nudges 3.40 out of 5.50's way.
    AminoAcid *aa5x33 = p.get_residue_bw("5.33");
    AminoAcid *aa5x43 = p.get_residue_bw("5.43");
    AminoAcid *aa5x54 = p.get_residue_bw("5.54");
    AminoAcid *aa5x68 = p.get_residue_bw("5.68");
    AminoAcid *aa6x44 = p.get_residue_bw("6.44");
    SCoord TMR5edir = p.get_region_center(aa6x44->get_residue_no(), aa6x55->get_residue_no()).subtract(p.get_region_center(aa5x43->get_residue_no(), aa5x54->get_residue_no()));
    float TMR5ez = p.region_can_move(aa5x43->get_residue_no(), aa5x54->get_residue_no(), TMR5edir, true);
    cout << "TMR5 moves " << TMR5ez << " toward TMR6." << endl;

    // Perform the TMR5ez slide.
    TMR5edir.r = TMR5ez;
    p.move_piece(aa5x33->get_residue_no(), aa5x68->get_residue_no(), TMR5edir);

    // If TMR6ex is nonzero:
    float theta;
    Point was;
    AminoAcid *aa56x50 = p.get_residue_bw("56.50");
    AminoAcid *aa6x28 = p.get_residue_bw("6.28");
    AminoAcid *aa5x50 = p.get_residue_bw("5.50");
    SCoord TMR6c, axis5;
    if (TMR6ex)
    {
        // Compute the angle to rotate about 6.48 and perform the rotation. Measure how far 56.50 moved; call it TMR6c.
        theta = atan2(TMR6ex, TMR6ey);
        was = aa56x50->get_CA_location();
        p.rotate_piece(aa56x50->get_residue_no(), aa6x59->get_residue_no(), aa6x48->get_CA_location(), axis6, theta);
        TMR6c = aa56x50->get_CA_location().subtract(was);

        // Compute the axis and angle to rotate TMR5 about 5.33 to match 56.49 to TMR6c, and perform the rotation.
        axis5 = compute_normal(aa5x50->get_CA_location(), aa5x68->get_CA_location(), aa5x68->get_CA_location().add(TMR6c));
        theta = find_3d_angle(aa5x68->get_CA_location(), aa5x68->get_CA_location().add(TMR6c), aa5x50->get_CA_location());
        p.rotate_piece(aa5x50->get_residue_no(), aa56x50->get_residue_no()-1, aa5x50->get_CA_location(), axis5, theta);
        cout << "TMR5 rotation about " << axis5 << " by " << (theta*fiftyseven) << "deg." << endl;
    }

    // If there is Y5.58 and Y7.53:
    AminoAcid *aa5x58 = p.get_residue_bw("5.58");
    AminoAcid *aa6x40 = p.get_residue_bw("6.40");
    AminoAcid *aa7x48 = p.get_residue_bw("7.48");
    AminoAcid *aa7x52 = p.get_residue_bw("7.52");
    AminoAcid *aa7x53 = p.get_residue_bw("7.53");
    AminoAcid *aa7x56 = p.get_residue_bw("7.56");
    char l5x58 = aa5x58->get_letter();
    char l7x53 = aa7x53->get_letter();
    float bridge57, scooch6x40, TMR7cz;
    SCoord TMR5cdir, TMR6cdir, TMR7cdir, axis7;
    Point pt_tmp;
    
    if (l5x58 == 'Y' && l7x53 == 'Y')
    {
        // Attempt to bridge 5.58~7.53 and measure the distance necessary to complete the contact. Call it Bridge57.
        p.bridge(aa5x58->get_residue_no(), aa7x53->get_residue_no());
        bridge57 = fmax(0, aa5x58->get_atom("OH")->distance_to(aa7x53->get_atom("OH")) - 3);
        cout << "5.58~7.53 bridge must move " << bridge57 << "A together." << endl;

        // Move the side chain of 6.40 to face 7.53:CA and ensure that 6.40 is not clashing with 7.53.
        aa6x40->movability = MOV_ALL;
        aa6x40->conform_atom_to_location(aa6x40->get_reach_atom()->name, aa7x53->get_CA_location());
        
        // If 6.40 is still clashing with 7.53, compute the distance to rotate 6.48 thru 56.50 about 6.48 to eliminate the clash.
        // Then perform the rotation, then compute the rotation of TMR5 about 5.33 to keep up, and perform that rotation.
        // Then re-form the 5.58~7.53 bridge.
        float clash = aa6x40->get_intermol_clashes(aa7x53);
        if (clash > homology_clash_peraa/4)
        {
            Pose pose6x40;
            pose6x40.copy_state(aa6x40);
            TMR6cdir = aa6x40->get_barycenter().subtract(aa7x53->get_atom_location("HH"));
            TMR6cdir.r = 0.1;

            aa6x40->movability = MOV_ALL;
            for (l=0; l<50; l++)
            {
                aa6x40->aamove(TMR6cdir);
                scooch6x40 += TMR6cdir.r;
                pt_tmp = aa6x40->get_CA_location();
                clash = aa6x40->get_intermol_clashes(aa7x53);
                if (clash < homology_clash_peraa/4) break;
            }

            pose6x40.restore_state(aa6x40);
            cout << "Moving 6.40 " << scooch6x40 << " out of the way." << endl;

            axis6 = compute_normal(aa6x48->get_CA_location(), aa6x40->get_CA_location(), pt_tmp);
            theta = find_3d_angle(aa6x40->get_CA_location(), pt_tmp, aa6x48->get_CA_location());

            was = aa56x50->get_CA_location();
            p.rotate_piece(aa56x50->get_residue_no(), aa6x48->get_residue_no(), aa6x48->get_CA_location(), axis6, theta);
            TMR6c = aa56x50->get_CA_location().subtract(was);
            cout << "Bending TMR6 " << (theta*fiftyseven) << "deg in the CYT domain." << endl;

            // Compute the axis and angle to rotate TMR5 about 5.33 to match 56.49 to TMR6c, and perform the rotation.
            axis5 = compute_normal(aa5x50->get_CA_location(), aa5x68->get_CA_location(), aa5x68->get_CA_location().add(TMR6c));
            theta = find_3d_angle(aa5x68->get_CA_location(), aa5x68->get_CA_location().add(TMR6c), aa5x50->get_CA_location());
            p.rotate_piece(aa5x50->get_residue_no(), aa56x50->get_residue_no()-1, aa5x50->get_CA_location(), axis5, theta);
        
            p.bridge(aa5x58->get_residue_no(), aa7x53->get_residue_no());
        }
        
        // Re-measure Bridge57. If it is nonzero, determine how far 7.53 can move toward 5.58 without clashing with 3.43. Call it TMR7cz.
        bridge57 = fmax(0, aa5x58->get_atom("OH")->distance_to(aa7x53->get_atom("OH")) - 3);
        cout << "5.58~7.53 bridge must move " << bridge57 << "A together." << endl;
        TMR7cdir = aa5x58->get_CA_location().subtract(aa7x53->get_CA_location());
        TMR7cz = p.region_can_move(aa7x53->get_residue_no(), aa7x53->get_residue_no(), TMR7cdir, false);
        
        // Compute and execute a bend of TMR7 at 7.48 to move 7.53 the minimum distance of TMR7cz and Bridge57 toward 5.58.
        TMR7cdir.r = fmin(TMR7cz, bridge57);
        axis7 = compute_normal(aa7x48->get_CA_location(), aa7x53->get_CA_location(), aa7x53->get_CA_location().add(TMR7cdir));
        theta = find_3d_angle(aa7x53->get_CA_location(), aa7x53->get_CA_location().add(TMR7cdir), aa7x48->get_CA_location());
        p.rotate_piece(aa7x48->get_residue_no(), aa7x56->get_residue_no(), aa7x48->get_CA_location(), axis7, theta);
        
        // Re-measure Bridge57. If it is nonzero, compute and execute a pivot of TMR5 from 5.33 to move 5.58 the rest of the way to make contact with 7.53.
        // Then compute and execute a y-axis rotation of TMR6 to bring 6.28 as far along horizontally as 5.68 moved.
        // Translate the CYT3 region (BW numbers 56.x) to stay with 5.68 and 6.28 as smoothly as possible.
        bridge57 = fmax(0, aa5x58->get_atom("OH")->distance_to(aa7x53->get_atom("OH")) - 3);
        cout << "5.58~7.53 bridge must move " << bridge57 << "A together." << endl;
        TMR5cdir = aa7x53->get_CA_location().subtract(aa5x58->get_CA_location());
        TMR5cdir.r = bridge57;
        axis5 = compute_normal(aa5x33->get_CA_location(), aa5x58->get_CA_location(), aa5x58->get_CA_location().add(TMR5cdir));
        theta = find_3d_angle(aa5x58->get_CA_location(), aa5x58->get_CA_location().add(TMR5cdir), aa5x33->get_CA_location());
        was = aa56x50->get_CA_location();
        p.rotate_piece(aa5x33->get_residue_no(), aa56x50->get_residue_no(), aa5x33->get_CA_location(), axis5, theta);

        TMR6cdir = aa56x50->get_CA_location().subtract(was);
        axis6 = compute_normal(aa6x59->get_CA_location(), was, aa56x50->get_CA_location());
        theta = find_3d_angle(was, aa56x50->get_CA_location(), aa6x59->get_CA_location());
        p.rotate_piece(aa56x50->get_residue_no()+1, aa6x59->get_residue_no(), aa6x59->get_CA_location(), axis6, theta);

        p.bridge(aa5x58->get_residue_no(), aa7x53->get_residue_no());
        
        // If 6.40 is once again clashing with 7.53, compute the distance to rotate 6.48 thru 56.50 about 6.48 to eliminate the clash.
        // Then perform the rotation, then compute the rotation of TMR5 about 5.33 to keep up, and perform that rotation.
        // Then re-form the 5.58~7.53 bridge.
        scooch6x40 = 0;
        clash = aa6x40->get_intermol_clashes(aa7x52) + aa6x40->get_intermol_clashes(aa7x53);
        if (clash > homology_clash_peraa/4)
        {
            Pose pose6x40;
            pose6x40.copy_state(aa6x40);
            TMR6cdir = aa6x40->get_CA_location().subtract(aa7x53->get_CA_location());
            TMR6cdir.r = 0.1;

            aa6x40->movability = MOV_ALL;
            for (l=0; l<50; l++)
            {
                aa6x40->aamove(TMR6cdir);
                scooch6x40 += TMR6cdir.r;
                pt_tmp = aa6x40->get_CA_location();
                clash = aa6x40->get_intermol_clashes(aa7x52) + aa6x40->get_intermol_clashes(aa7x53);
                if (clash < homology_clash_peraa/4) break;
            }

            pose6x40.restore_state(aa6x40);
            cout << "Moving 6.40 again " << scooch6x40 << " out of the way." << endl;

            axis6 = compute_normal(aa6x48->get_CA_location(), aa6x40->get_CA_location(), pt_tmp);
            theta = find_3d_angle(aa6x40->get_CA_location(), pt_tmp, aa6x48->get_CA_location());

            was = aa56x50->get_CA_location();
            p.rotate_piece(aa56x50->get_residue_no(), aa6x48->get_residue_no(), aa6x48->get_CA_location(), axis6, theta);
            TMR6c = aa56x50->get_CA_location().subtract(was);
            cout << "Bending TMR6 a further " << (theta*fiftyseven) << "deg in the CYT domain." << endl;

            // Compute the axis and angle to rotate TMR5 about 5.33 to match 56.49 to TMR6c, and perform the rotation.
            axis5 = compute_normal(aa5x58->get_CA_location(), aa5x68->get_CA_location(), aa5x68->get_CA_location().add(TMR6c));
            theta = find_3d_angle(aa5x68->get_CA_location(), aa5x68->get_CA_location().add(TMR6c), aa5x58->get_CA_location());
            p.rotate_piece(aa5x58->get_residue_no(), aa56x50->get_residue_no()-1, aa5x58->get_CA_location(), axis5, theta);
        
            p.bridge(aa5x58->get_residue_no(), aa7x53->get_residue_no());
        }
    }

    // If there is R6.59, adjust its side chain to keep pointing inward while avoiding clashes with other nearby side chains.

    // Now save the output file. It will be the same name as the input file except it will end in .icactive.pdb instead of .upright.pdb.
    std::string out_filename = path + orid + (std::string)".icactive.pdb";
    fp = fopen(out_filename.c_str(), "wb");
    if (!fp) return -3;
    p.save_pdb(fp);
    p.end_pdb(fp);
    fclose(fp);
    cout << "Saved " << out_filename << endl;

    // Finally, write all BRIDGE, STCR, and FLXR parameters to a constraints file ending in .params instead of .upright.pdb.
}

