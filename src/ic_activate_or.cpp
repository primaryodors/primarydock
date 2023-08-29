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

enum TMR6ActivationType
{
    Bend6,
    Swing6,
    Hybrid6,
    Rock6x59,
    Rock6Other
};

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

    std::vector<std::string> constraints;

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


    AminoAcid *aa3x21 = p.get_residue_bw("3.21");
    AminoAcid *aa3x40 = p.get_residue_bw("3.40");
    AminoAcid *aa3x50 = p.get_residue_bw("3.50");
    AminoAcid *aa3x56 = p.get_residue_bw("3.56");
    AminoAcid *aa4x53 = p.get_residue_bw("4.53");
    AminoAcid *aa4x60 = p.get_residue_bw("4.60");
    AminoAcid *aa4x61 = p.get_residue_bw("4.61");
    AminoAcid *aa4x64 = p.get_residue_bw("4.64");
    AminoAcid *aa45x51 = p.get_residue_bw("45.51");
    AminoAcid *aa45x52 = p.get_residue_bw("45.52");
    AminoAcid *aa45x53 = p.get_residue_bw("45.53");
    AminoAcid *aa5x33 = p.get_residue_bw("5.33");
    AminoAcid *aa5x36 = p.get_residue_bw("5.36");
    AminoAcid *aa5x39 = p.get_residue_bw("5.39");
    AminoAcid *aa5x43 = p.get_residue_bw("5.43");
    AminoAcid *aa5x50 = p.get_residue_bw("5.50");
    AminoAcid *aa5x54 = p.get_residue_bw("5.54");
    AminoAcid *aa5x58 = p.get_residue_bw("5.58");
    AminoAcid *aa5x68 = p.get_residue_bw("5.68");
    AminoAcid *aa56x50 = p.get_residue_bw("56.50");
    AminoAcid *aa6x28 = p.get_residue_bw("6.28");
    AminoAcid *aa6x40 = p.get_residue_bw("6.40");
    AminoAcid *aa6x44 = p.get_residue_bw("6.44");
    AminoAcid *aa6x48 = p.get_residue_bw("6.48");
    AminoAcid *aa6x49 = p.get_residue_bw("6.49");
    AminoAcid *aa6x55 = p.get_residue_bw("6.55");
    AminoAcid *aa6x58 = p.get_residue_bw("6.58");
    AminoAcid *aa6x59 = p.get_residue_bw("6.59");
    AminoAcid *aa7x48 = p.get_residue_bw("7.48");
    AminoAcid *aa7x52 = p.get_residue_bw("7.52");
    AminoAcid *aa7x53 = p.get_residue_bw("7.53");
    AminoAcid *aa7x56 = p.get_residue_bw("7.56");
    
    int n3x21 = aa3x21->get_residue_no();
    int n3x56 = aa3x56->get_residue_no();
    int n4x53 = aa4x53->get_residue_no();
    int n4x64 = aa4x64->get_residue_no();
    int n45x51 = aa45x51->get_residue_no();
    int n6x28 = aa6x28->get_residue_no();
    int n6x48 = aa6x48->get_residue_no();
    int n6x55 = aa6x55->get_residue_no();
    int n6x58 = aa6x58->get_residue_no();
    int n6x59 = aa6x59->get_residue_no();

    char l3x40 = aa3x40->get_letter();
    char l45x51 = aa45x51->get_letter();
    char l45x52 = aa45x52->get_letter();
    char l45x53 = aa45x53->get_letter();
    char l5x58 = aa5x58->get_letter();
    char l6x48 = aa6x48->get_letter();
    char l6x55 = aa6x55->get_letter();
    char l6x58 = aa6x58->get_letter();
    char l6x59 = aa6x59->get_letter();
    char l7x53 = aa7x53->get_letter();


    // If there is an R at position 6.59, perform a slight helix unwind and bring 6.59's side chain inward.
    DynamicMotion unwind6(&p);
    DynamicMotion bend45(&p);
    if (aa6x59 && aa6x58 && aa45x53 && (aa6x59->get_letter() == 'R'))
    {
        if (l45x53 == 'Q' || l45x53 == 'N' || l45x53 == 'R' || l45x53 == 'K' || l45x53 == 'E' || l45x53 == 'D')
        {
            // Have to move 45.53 and 45.54 out of the way, but leave 45.52 where it is.
            bend45.start_resno.from_string("45.53");
            bend45.end_resno.from_string("45.54");
            bend45.type = dyn_bend;
            bend45.fulcrum_resno.from_string("45.53");
            bend45.axis_resno.from_string("45.54");
            bend45.bias = 60;
            bend45.apply_absolute(1.5);
            aa45x53->conform_atom_to_location(aa45x53->get_reach_atom()->name, Point(0,-1000,0));
        }
        else if (l45x53 == 'M' || l45x53 == 'Y')
        {
            // OR52D1 fix, hopefully.
            constraints.push_back((std::string)"ATOMTO 45.53 EXTENT 5.30");
            constraints.push_back((std::string)"STCR 45.53");
        }

        cout << "R6.59" << endl;
        unwind6.start_resno.from_string("6.54");
        unwind6.end_resno.from_string("6.59");
        unwind6.type = dyn_wind;
        unwind6.bias = -13;
        unwind6.apply_absolute(1);

        aa6x59->conform_atom_to_location("CZ", aa4x60->get_CA_location());
        constraints.push_back((std::string)"STCR 6.59");
    }

    float theta6 = 0;
    SCoord axis6;
    TMR6ActivationType type6;

    if (l6x55 == 'Y' && (l45x51 == 'D' || l45x51 == 'E'))
    {
        // TODO: In OR8H1, F3.32 impinges on where this bridge would form, so it would be more sterically favorable to first point
        // 6.55's EXTENT toward 45.53's CA and then make the bridge. Since this depends on the exact positioning of atoms around the
        // neighborhood of 45.51 as well as the side chains of 45.52 and 45.53, it will be different for each receptor.
        p.bridge(n6x55, n45x51);
        constraints.push_back((std::string)"BRIDGE 6.55 45.51");
        float r = aa6x55->get_atom("OH")->distance_to(aa45x51->get_nearest_atom(aa6x55->get_atom_location("HH")));

        if (r <= 3)
        {
            if (l6x48 == 'Y' && (l3x40 == 'S' || l3x40 == 'T' && l3x40 == 'N' || l3x40 == 'Q' || l3x40 == 'E' || l3x40 == 'D'))
            {
                type6 = Bend6; // No motion of TMR6 in the EXR domain.
                theta6 = 0;
                cout << "Bend6 type activation." << endl;
            }
            else
            {
                type6 = Swing6;
                axis6 = compute_normal(aa6x55->get_CA_location(), aa3x40->get_CA_location(), aa6x48->get_CA_location());
                theta6 = fiftyseventh * 15;      // TODO:
                cout << "Swing6 type activation." << endl;
            }
        }
        else
        {
            type6 = Hybrid6;
            axis6 = compute_normal(aa6x48->get_CA_location(), aa6x55->get_CA_location(), aa45x51->get_CA_location());
            SCoord span = aa45x51->get_CA_location().subtract(aa6x55->get_CA_location());
            span.r = r - 3;
            Point pt = aa6x55->get_CA_location().add(span);
            LocatedVector lv = axis6;
            lv.origin = aa6x48->get_CA_location();
            theta6 = fmin(p.region_can_rotate(n6x48, n6x55, lv, true),
                find_3d_angle(pt, aa6x55->get_CA_location(), lv.origin));
            cout << "Hybrid6 type activation with a " << (fiftyseven * theta6) << "deg EXR component." << endl;
        }
    }
    else
    {
        if (l6x59 == 'R' || l6x59 == 'K')
        {
            type6 = Rock6x59;
            axis6 = compute_normal(aa6x48->get_CA_location(), aa6x59->get_CA_location(), aa5x36->get_CA_location());
            LocatedVector lv = axis6;
            lv.origin = aa6x48->get_CA_location();
            theta6 = p.region_can_rotate(n6x48, n6x59, lv, true, 5);
            cout << "Rock6 type activation with a " << (fiftyseven * theta6) << "deg basic 6.59." << endl;
        }
        else
        {
            type6 = Rock6Other;
            axis6 = compute_normal(aa6x48->get_CA_location(), aa6x58->get_CA_location(), aa45x53->get_CA_location());
            LocatedVector lv = axis6;
            lv.origin = aa6x48->get_CA_location();
            theta6 = p.region_can_rotate(n6x48, n6x58, lv, true);
            cout << "Rock6 type activation with a " << (fiftyseven * theta6) << "deg pivot." << endl;
        }
    }


    // Perform a slight rotation of TMR3 to give TMR5 more room to move.
    
    LocatedVector axis3 = (SCoord)aa3x50->get_CA_location().subtract(aa5x58->get_CA_location());
    axis3.origin = aa3x50->get_CA_location();
    float rock3 = p.region_can_rotate(n3x21, n3x56, axis3, true);
    cout << "TMR3 can rotate " << (rock3 * fiftyseven) << " degrees about 3.50 to make room for TMR5." << endl;

    // Perform the rotation.
    // TODO: This should not be hard coded.
    p.rotate_piece(n3x21, n3x56, aa3x50->get_CA_location(), axis3, rock3*2);

    LocatedVector axis4 = (SCoord)aa4x53->get_CA_location().subtract(aa6x49->get_CA_location());
    axis4.origin = aa4x53->get_CA_location();
    float bend4 = p.region_can_rotate(n4x53, n4x64, axis4, true);
    cout << "TMR4 can bend " << (bend4 * fiftyseven) << " degrees about 4.64 to meet TMR3." << endl;

    // Perform the rotation.
    p.rotate_piece(n4x53, n4x64, aa4x53->get_CA_location(), axis4, bend4);


    // Measure how far 5.43-5.54 can move toward 6.44-6.55 without clashing. Call it TMR5ez.
    
    SCoord TMR5edir = p.get_region_center(aa6x44->get_residue_no(), aa6x55->get_residue_no()).subtract(p.get_region_center(aa5x43->get_residue_no(), aa5x54->get_residue_no()));
    float TMR5ez = p.region_can_move(aa5x43->get_residue_no(), aa5x54->get_residue_no(), TMR5edir, true);
    cout << "TMR5 moves " << TMR5ez << " toward TMR6." << endl;

    // Perform the TMR5ez slide.
    TMR5edir.r = TMR5ez;
    p.move_piece(aa5x33->get_residue_no(), aa5x68->get_residue_no(), TMR5edir);

    // If TMR6ex is nonzero:
    float theta;
    Point was;
    SCoord TMR6c, axis5;
    if (theta6)
    {
        // Compute the angle to rotate about 6.48 and perform the rotation. Measure how far 56.50 moved; call it TMR6c.
        was = aa56x50->get_CA_location();
        p.rotate_piece(aa56x50->get_residue_no(), aa6x59->get_residue_no(), aa6x48->get_CA_location(), axis6, theta6);
        TMR6c = aa56x50->get_CA_location().subtract(was);

        // Compute the axis and angle to rotate TMR5 about 5.50 to match 56.49 to TMR6c, and perform the rotation.
        axis5 = compute_normal(aa5x50->get_CA_location(), aa5x68->get_CA_location(), aa5x68->get_CA_location().add(TMR6c));
        theta = find_3d_angle(aa5x68->get_CA_location(), aa5x68->get_CA_location().add(TMR6c), aa5x50->get_CA_location());
        p.rotate_piece(aa5x50->get_residue_no(), aa56x50->get_residue_no()-1, aa5x50->get_CA_location(), axis5, theta);
        cout << "TMR5 rotation about " << axis5 << " by " << (theta*fiftyseven) << "deg." << endl;
    }

    // If there is Y5.58 and Y7.53:
    float bridge57, scooch6x40 = 0, TMR7cz;
    SCoord TMR5cdir, TMR6cdir, TMR7cdir, axis7;
    Point pt_tmp;
    
    if (l5x58 == 'Y' && l7x53 == 'Y')
    {
        constraints.push_back((std::string)"BRIDGE 5.58 7.53");

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

        #if 0

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
        bridge57 = fmax(0, aa5x58->get_atom("OH")->distance_to(aa7x53->get_atom("OH")) - 3);
        cout << "5.58~7.53 bridge must move " << bridge57 << "A together." << endl;

        #endif
        
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

    cout << "Minimizing internal clashes..." << endl;
    p.get_internal_clashes(1, p.get_end_resno(), true);

    // If there is R6.59, adjust its side chain to keep pointing inward while avoiding clashes with other nearby side chains.
    if (l6x59 == 'R' || l6x59 == 'K')
    {
        cout << "Optimizing 6.59 side chain..." << endl;

        Molecule m("Acid");
        m.from_smiles("[Cl-]");
        Point pt = aa4x61->get_CA_location().add(aa6x59->get_CA_location()); //.add(aa4x60->get_CA_location());
        pt.scale(pt.magnitude()/2);
        m.move(pt);
        m.movability = MOV_PINNED;

        Molecule* mols[256];
        AminoAcid** can_clash = p.get_residues_can_clash(aa6x59->get_residue_no());
        l = 0;

        mols[l++] = (Molecule*)aa6x59;
        aa6x59->movability = MOV_FLEXONLY;
        mols[l++] = &m;

        for (i=0; can_clash[i]; i++)
        {
            can_clash[i]->movability = MOV_FLEXONLY;
            mols[l++] = (Molecule*)can_clash[i];
            if (l >= 255) break;
        }
        mols[l] = nullptr;
        Molecule::conform_molecules(mols, 50);

        mols[3] = nullptr;
        Molecule::conform_molecules(mols, 20);
        
        /*fp = fopen("tmp/dbg6x59.pdb", "wb");
        if (!fp) return -3;
        p.save_pdb(fp, &m);
        p.end_pdb(fp);
        fclose(fp);*/
    }

    // Now save the output file. It will be the same name as the input file except it will end in .active.pdb instead of .upright.pdb.
    std::string out_filename = path + orid + (std::string)".active.pdb";
    fp = fopen(out_filename.c_str(), "wb");
    if (!fp) return -3;
    p.save_pdb(fp);
    p.end_pdb(fp);
    fclose(fp);
    cout << "Saved " << out_filename << endl;

    // Finally, write all BRIDGE, STCR, and FLXR parameters to a constraints file ending in .params instead of .upright.pdb.
    std::string cns_filename = path + orid + (std::string)".params";
    fp = fopen(cns_filename.c_str(), "wb");
    if (!fp) return -3;
    n = constraints.size();
    for (l=0; l<n; l++) fprintf(fp, "%s\n", constraints[l].c_str());
    fclose(fp);
}

