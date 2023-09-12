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

    ////////////////////////////////////////////////////////////////////////////////
    // Read program arguments
    ////////////////////////////////////////////////////////////////////////////////

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


    ////////////////////////////////////////////////////////////////////////////////
    // Load protein
    ////////////////////////////////////////////////////////////////////////////////

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


    ////////////////////////////////////////////////////////////////////////////////
    // Set up residue vars
    ////////////////////////////////////////////////////////////////////////////////

    AminoAcid *aa1x32 = p.get_residue_bw("1.32");
    AminoAcid *aa1x50 = p.get_residue_bw("1.50");
    AminoAcid *aa2x49 = p.get_residue_bw("2.49");
    AminoAcid *aa2x50 = p.get_residue_bw("2.50");
    AminoAcid *aa2x66 = p.get_residue_bw("2.66");
    AminoAcid *aa3x21 = p.get_residue_bw("3.21");
    AminoAcid *aa3x32 = p.get_residue_bw("3.32");
    AminoAcid *aa3x33 = p.get_residue_bw("3.33");
    AminoAcid *aa3x36 = p.get_residue_bw("3.36");
    AminoAcid *aa3x39 = p.get_residue_bw("3.39");
    AminoAcid *aa3x40 = p.get_residue_bw("3.40");
    AminoAcid *aa3x50 = p.get_residue_bw("3.50");
    AminoAcid *aa3x56 = p.get_residue_bw("3.56");
    AminoAcid *aa4x53 = p.get_residue_bw("4.53");
    AminoAcid *aa4x55 = p.get_residue_bw("4.55");
    AminoAcid *aa4x56 = p.get_residue_bw("4.56");
    AminoAcid *aa4x60 = p.get_residue_bw("4.60");
    AminoAcid *aa4x61 = p.get_residue_bw("4.61");
    AminoAcid *aa4x64 = p.get_residue_bw("4.64");
    AminoAcid *aa45x51 = p.get_residue_bw("45.51");
    AminoAcid *aa45x52 = p.get_residue_bw("45.52");
    AminoAcid *aa45x53 = p.get_residue_bw("45.53");
    AminoAcid *aa45x54 = p.get_residue_bw("45.54");
    AminoAcid *aa5x33 = p.get_residue_bw("5.33");
    AminoAcid *aa5x35 = p.get_residue_bw("5.35");
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
    AminoAcid *aa6x53 = p.get_residue_bw("6.53");
    AminoAcid *aa6x54 = p.get_residue_bw("6.54");
    AminoAcid *aa6x55 = p.get_residue_bw("6.55");
    AminoAcid *aa6x58 = p.get_residue_bw("6.58");
    AminoAcid *aa6x59 = p.get_residue_bw("6.59");
    AminoAcid *aa7x31 = p.get_residue_bw("7.31");
    AminoAcid *aa7x37 = p.get_residue_bw("7.37");
    AminoAcid *aa7x41 = p.get_residue_bw("7.41");
    AminoAcid *aa7x42 = p.get_residue_bw("7.42");
    AminoAcid *aa7x43 = p.get_residue_bw("7.43");
    AminoAcid *aa7x48 = p.get_residue_bw("7.48");
    AminoAcid *aa7x49 = p.get_residue_bw("7.49");
    AminoAcid *aa7x52 = p.get_residue_bw("7.52");
    AminoAcid *aa7x53 = p.get_residue_bw("7.53");
    AminoAcid *aa7x56 = p.get_residue_bw("7.56");
    
    int n1x32 = aa1x32->get_residue_no();
    int n1x50 = aa1x50->get_residue_no();
    int n2x50 = aa2x50->get_residue_no();
    int n2x66 = aa2x66->get_residue_no();
    int n3x21 = aa3x21->get_residue_no();
    int n3x39 = aa3x39->get_residue_no();
    int n3x40 = aa3x40->get_residue_no();
    int n3x56 = aa3x56->get_residue_no();
    int n4x53 = aa4x53->get_residue_no();
    int n4x64 = aa4x64->get_residue_no();
    int n45x51 = aa45x51->get_residue_no();
    int n56x50 = aa56x50->get_residue_no();
    int n6x28 = aa6x28->get_residue_no();
    int n6x48 = aa6x48->get_residue_no();
    int n6x55 = aa6x55->get_residue_no();
    int n6x58 = aa6x58->get_residue_no();
    int n6x59 = aa6x59->get_residue_no();
    int n7x31 = aa7x31->get_residue_no();
    int n7x41 = aa7x41->get_residue_no();
    int n7x43 = aa7x43->get_residue_no();
    int n7x49 = aa7x49->get_residue_no();

    char l3x40 = aa3x40->get_letter();
    char l4x56 = aa4x56->get_letter();
    char l45x51 = aa45x51->get_letter();
    char l45x52 = aa45x52->get_letter();
    char l45x53 = aa45x53->get_letter();
    char l5x58 = aa5x58->get_letter();
    char l6x48 = aa6x48->get_letter();
    char l6x55 = aa6x55->get_letter();
    char l6x58 = aa6x58->get_letter();
    char l6x59 = aa6x59->get_letter();
    char l7x41 = aa7x41->get_letter();
    char l7x53 = aa7x53->get_letter();

    float dist67 = aa6x59->get_CA_location().get_3d_distance(p.get_residue(n6x59+1)->get_CA_location());


    ////////////////////////////////////////////////////////////////////////////////
    // Set up vars for processing.
    ////////////////////////////////////////////////////////////////////////////////

    Pose pose6x55;
    bool preserve6x55 = false;

    Bond* b;
    bool bcr;
    bool stiff6x55 = false;

    LocatedVector axis1, axis2, axis3, axis4, axis7;
    float theta;
    Point was;
    SCoord TMR6c, axis5;
    float theta6 = 0;
    SCoord axis6;
    float bridge57, scooch6x40 = 0, TMR7cz;
    SCoord TMR5cdir, TMR6cdir, TMR7cdir;
    Point pt_tmp;
    TMR6ActivationType type6;

    // If there is an R at position 6.59, perform a slight helix unwind and bring 6.59's side chain inward.
    DynamicMotion unwind6(&p);
    DynamicMotion bend45(&p);


    ////////////////////////////////////////////////////////////////////////////////
    // Receptor-specific fixes.
    ////////////////////////////////////////////////////////////////////////////////

    if (!strcmp(orid.c_str(), "OR2M3"))
    {
        // aa7x41->conform_atom_to_location("OH", aa6x48->get_CA_location());
        p.bridge(n3x39, n7x41);
        aa6x53->get_atom("CA")->get_bond_between(aa6x53->get_atom("CB"))->rotate(60*fiftyseventh);
        // aa6x54->conform_atom_to_location(aa6x54->get_reach_atom()->name, aa45x51->get_CA_location());

        unwind6.start_resno.from_string("6.54");
        unwind6.end_resno.from_string("6.55");
        unwind6.type = dyn_wind;
        unwind6.bias = -60;
        unwind6.apply_absolute(1);

        pt_tmp = aa1x32->get_CA_location();
        SCoord delta = aa7x37->get_CA_location().subtract(aa45x51->get_CA_location());
        delta.r = 2;
        pt_tmp = pt_tmp.add(delta);
        axis7 = compute_normal(aa7x41->get_CA_location(), aa7x31->get_CA_location(), pt_tmp);
        LocatedVector lv = axis7;
        lv.origin = aa7x41->get_CA_location();
        theta = fmin(p.region_can_rotate(n7x31, n7x41, lv, true), fiftyseventh*10);
        p.rotate_piece(n7x31, n7x41, lv.origin, axis7, theta);

        axis6 = compute_normal(aa6x48->get_CA_location(), aa6x55->get_CA_location(), aa7x31->get_CA_location());
        lv = axis6;
        lv.origin = aa6x48->get_CA_location();
        theta = fmin(p.region_can_rotate(n6x48, n6x59, lv, true), fiftyseventh*2.5);
        p.rotate_piece(n6x48, n6x59, lv.origin, axis7, theta);

        // aa3x33->conform_atom_to_location("OG2", aa5x43->get_CA_location());
    }
    else if (!strcmp(orid.c_str(), "OR8H1"))
    {
        stiff6x55 = true;
    }
    else if (l7x41 == 'Y')
    {
        aa7x41->conform_atom_to_location(aa7x41->get_reach_atom()->name, aa45x51->get_CA_location());
    }

    if (l4x56 == 'D' || l4x56 == 'E')
    {
        const char* atom_name = (l4x56 == 'D') ? "CG" : "CD";
        SCoord v = aa4x56->get_CA_location().subtract(aa3x33->get_CA_location());
        v.r = 10000;
        aa4x56->conform_atom_to_location(atom_name, v);
        constraints.push_back((std::string)"STCR 4.56");
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Arg6.59 unwind.
    ////////////////////////////////////////////////////////////////////////////////

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
        // constraints.push_back((std::string)"STCR 6.59");
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Determine TMR6 activation type.
    ////////////////////////////////////////////////////////////////////////////////

    if (l6x55 == 'Y' && (l45x51 == 'D' || l45x51 == 'E'))
    {
        aa6x55->conform_atom_to_location(aa6x55->get_reach_atom()->name, aa45x51->get_reach_atom()->get_location());
        aa45x51->conform_atom_to_location(aa45x51->get_reach_atom()->name, aa6x55->get_reach_atom()->get_location());
        aa6x55->conform_atom_to_location(aa6x55->get_reach_atom()->name, aa45x51->get_reach_atom()->get_location());
        aa45x51->conform_atom_to_location(aa45x51->get_reach_atom()->name, aa6x55->get_reach_atom()->get_location());

        // constraints.push_back((std::string)"BRIDGE 6.55 45.51");
        float r = aa6x55->get_atom("OH")->distance_to(aa45x51->get_nearest_atom(aa6x55->get_atom_location("HH")));
        aa6x55->conform_atom_to_location(aa6x55->get_reach_atom()->name, Point(0,10000,0));

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
            // SCoord span = aa45x51->get_CA_location().subtract(aa6x55->get_CA_location());
            SCoord span = aa45x51->get_reach_atom()->get_location().subtract(aa6x55->get_reach_atom()->get_location());
            span.r = r - 3;
            Point pt = aa6x55->get_CA_location().add(span);
            LocatedVector lv = axis6;
            lv.origin = aa6x48->get_CA_location();
            theta6 = fmin(p.region_can_rotate(n6x48, n6x55, lv, true), find_3d_angle(pt, aa6x55->get_CA_location(), lv.origin));
            
            cout << "Hybrid6 type activation with a " << (fiftyseven * theta6) << "deg EXR component limited by "
                << *(p.stop1) << "->" << *(p.stop2) << "." << endl;
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



    ////////////////////////////////////////////////////////////////////////////////
    // Slight bend of TMR3 to allow TMR5 motion.
    ////////////////////////////////////////////////////////////////////////////////
    
    axis3 = (SCoord)aa3x50->get_CA_location().subtract(aa5x58->get_CA_location());
    axis3.origin = aa3x40->get_CA_location();
    float rock3 = p.region_can_rotate(n3x21, n3x40, axis3, true);
    cout << "TMR3 can rotate " << (rock3 * fiftyseven) << " degrees about 3.40 limited by " << *(p.stop1) << "->" << *(p.stop2)
        << " to make room for TMR5." << endl;

    // Perform the rotation.
    p.rotate_piece(n3x21, n3x56, aa3x50->get_CA_location(), axis3, rock3*0.3);

    aa3x33->conform_atom_to_location(aa3x33->get_reach_atom()->name, aa5x35->get_CA_location());


    ////////////////////////////////////////////////////////////////////////////////
    // TMR4 bend to meet TMR3.
    ////////////////////////////////////////////////////////////////////////////////

    axis4 = (SCoord)aa4x53->get_CA_location().subtract(aa6x49->get_CA_location());
    axis4.origin = aa4x53->get_CA_location();
    float bend4 = p.region_can_rotate(n4x53, n4x64, axis4, true);
    cout << "TMR4 can bend " << (bend4 * fiftyseven) << " degrees about 4.64 to meet TMR3." << endl;

    // Perform the rotation.
    p.rotate_piece(n4x53, n4x64, aa4x53->get_CA_location(), axis4, bend4);


    ////////////////////////////////////////////////////////////////////////////////
    // TMR2 and TMR1 bend also.
    ////////////////////////////////////////////////////////////////////////////////

    axis2 = compute_normal(aa2x50->get_CA_location(), aa2x66->get_CA_location(), aa3x21->get_CA_location());
    axis2.origin = aa2x50->get_CA_location();
    float bend2 = bend4 * 1.2;
    p.rotate_piece(n2x50, n2x66, axis2.origin, axis2, -bend2);

    axis1 = compute_normal(aa1x50->get_CA_location(), aa1x32->get_CA_location(), aa2x66->get_CA_location());
    axis1.origin = aa1x50->get_CA_location();
    float bend1 = bend4 * 1.8;
    p.rotate_piece(n1x32, n1x50, axis1.origin, axis1, bend1);

    axis1 = (SCoord)(aa2x49->get_CA_location().subtract(aa1x50->get_CA_location()));
    bend1 = bend4 * 2;
    p.rotate_piece(n1x32, n1x50, axis1.origin, axis1, -bend1);


    ////////////////////////////////////////////////////////////////////////////////
    // Tyr6.55 verticality.
    ////////////////////////////////////////////////////////////////////////////////

    if (l6x55 == 'Y')
    {
        Point he1 = aa6x55->get_atom_location("HE1");
        Point he2 = aa6x55->get_atom_location("HE2");

        const char* up_atom = (he1.y > he2.y) ? "HE1" : "HE2";

        b = aa6x55->get_atom("CA")->get_bond_between(aa6x55->get_atom("CB"));
        b->can_rotate = false;
        b = aa6x55->get_atom("CB")->get_bond_between(aa6x55->get_atom("CG"));
        b->can_rotate = true;
        aa6x55->conform_atom_to_location(up_atom, Point(0,10000,0));
        b->can_rotate = false;
        b = aa6x55->get_atom("CA")->get_bond_between(aa6x55->get_atom("CB"));
        b->can_rotate = true;
        cout << "Vertical 6.55 ring." << endl;

        // Will be using this later.
        /* fp = fopen("tmp/vertical.pdb", "wb");
        p.save_pdb(fp);
        p.end_pdb(fp);
        fclose(fp); */
    }


    ////////////////////////////////////////////////////////////////////////////////
    // TMR5 motion toward TMR6.
    ////////////////////////////////////////////////////////////////////////////////

    // Measure how far 5.43-5.54 can move toward 6.44-6.55 without clashing. Call it TMR5ez.
    
    SCoord TMR5edir = p.get_region_center(aa6x44->get_residue_no(), aa6x55->get_residue_no()).subtract(p.get_region_center(aa5x43->get_residue_no(), aa5x54->get_residue_no()));
    float TMR5ez = p.region_can_move(aa5x43->get_residue_no(), aa5x54->get_residue_no(), TMR5edir, true);
    cout << "TMR5 moves " << TMR5ez << " toward TMR6 limited by " << *(p.stop1) << "->" << *(p.stop2) << "." << endl;

    // Perform the TMR5ez slide.
    TMR5edir.r = TMR5ez;
    p.move_piece(aa5x33->get_residue_no(), aa5x68->get_residue_no(), TMR5edir);


    ////////////////////////////////////////////////////////////////////////////////
    // ********** TMR6 motion. **********
    ////////////////////////////////////////////////////////////////////////////////

    // If TMR6ex is nonzero:
    if (theta6)
    {
        // Compute the angle to rotate about 6.48 and perform the rotation. Measure how far 56.50 moved; call it TMR6c.
        was = aa56x50->get_CA_location();
        p.rotate_piece(n56x50,
            (type6 == Swing6) ? n6x55 : n6x59,
            ((type6 == Swing6) ? aa6x55 : aa6x48)->get_CA_location(),
            axis6, theta6);
        cout << "TMR6 rotation about " << *((type6 == Swing6) ? aa6x55 : aa6x48) << " by " << (theta6*fiftyseven) << "deg." << endl;
        TMR6c = aa56x50->get_CA_location().subtract(was);

        // Compute the axis and angle to rotate TMR5 about 5.50 to match 56.49 to TMR6c, and perform the rotation.
        axis5 = compute_normal(aa5x50->get_CA_location(), aa5x68->get_CA_location(), aa5x68->get_CA_location().add(TMR6c));
        theta = find_3d_angle(aa5x68->get_CA_location(), aa5x68->get_CA_location().add(TMR6c), aa5x50->get_CA_location());
        p.rotate_piece(aa5x50->get_residue_no(), aa56x50->get_residue_no()-1, aa5x50->get_CA_location(), axis5, theta);
        cout << "TMR5 rotation about " << axis5 << " by " << (theta*fiftyseven) << "deg." << endl;

        // TMR6 takes TMR7 with it.
        float dist67_new = aa6x59->get_CA_location().get_3d_distance(p.get_residue(n6x59+1)->get_CA_location());
        int mid67 = (n6x59 + n7x31) / 2;
        AminoAcid* aamid67 = p.get_residue(mid67);
        pt_tmp = aamid67->get_CA_location();
        SCoord v = aa45x51->get_CA_location().subtract(pt_tmp);
        v.r *= 0.3;
        pt_tmp = pt_tmp.add(v);
        axis7 = (SCoord)(aa7x49->get_CA_location().subtract(pt_tmp));
        axis7.origin = aa7x49->get_CA_location();
        theta = find_angle_along_vector(p.get_residue(n6x59+1)->get_CA_location(), aa6x59->get_CA_location(), aa7x31->get_CA_location(), axis7);
        theta *= (dist67_new - dist67) / dist67_new;
        p.rotate_piece(n6x59+1, n7x49, axis7.origin, axis7, theta);
        cout << "TMR7 and EXR3 rotation about 7.43 by " << (theta*fiftyseven) << "deg to follow TMR6." << endl;
    }


    ////////////////////////////////////////////////////////////////////////////////
    // 6.55~45.51 bridge.
    ////////////////////////////////////////////////////////////////////////////////

    if (l6x55 == 'Y' && (l45x51 == 'D' || l45x51 == 'E'))
    {
        if (stiff6x55)
        {
            b = aa6x55->get_atom("CB")->get_bond_between(aa6x55->get_atom("CG"));
            b->can_rotate = false;
            constraints.push_back((std::string)"STCR 6.55 45.51");
        }
        else
        {
            constraints.push_back((std::string)"BRIDGE 6.55 45.51");
        }

        p.bridge(n6x55, n45x51);
        aa6x55->movability = aa45x51->movability = MOV_PINNED;
        cout << "Bridged 6.55 and 45.51." << endl;

        // In OR8H1, F3.32 impinges on where this bridge would form, so it would be more sterically favorable to first point
        // 6.55's EXTENT toward 45.53's CA and then make the bridge. Since this depends on the exact positioning of atoms around the
        // neighborhood of 45.51 as well as the side chains of 45.52 and 45.53, it will be different for each receptor.
        if (aa6x55->get_intermol_clashes(aa3x32) > 0)
        {
            cout << "3.32 fix for Y6.55 ~ D/E45.51 bridge." << endl;
            aa6x55->movability = MOV_FLEXONLY;

            if (l45x52 == 'S' || l45x52 == 'T' || l45x52 == 'A' || l45x52 == 'G')
            {
                cout << "Small 52." << endl;
                aa6x55->conform_atom_to_location("OH", aa45x54->get_CA_location());
                p.bridge(n6x55, n45x51);
                aa6x55->movability = MOV_PINNED;
                aa45x51->movability = MOV_FLEXONLY;
                aa45x51->conform_atom_to_location(aa45x51->get_reach_atom()->name, aa6x55->get_atom_location("HH"));
                aa45x51->movability = MOV_PINNED;
                preserve6x55 = true;
                pose6x55.copy_state(aa6x55);
            }
        }
    }


    ////////////////////////////////////////////////////////////////////////////////
    // TMR5-TMR7 cytoplasmic side contact.
    ////////////////////////////////////////////////////////////////////////////////

    // TODO: There should be a water molecule between 5.58 and 7.53, and R3.50 should coordinate to Y5.58.
    // Source: Manglik and Kruse, 2017 doi:10.1021/acs.biochem.7b00747.
    
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
    p.get_internal_clashes(n3x21, n3x56, true);


    ////////////////////////////////////////////////////////////////////////////////
    // Arg6.59 positioning.
    ////////////////////////////////////////////////////////////////////////////////

    if (type6 == Rock6x59 && (l6x59 == 'R' || l6x59 == 'K'))
    {
        cout << "Optimizing 6.59 side chain..." << endl;

        Molecule m("Acid");
        m.from_smiles("[Cl-]");
        Point pt = aa4x61->get_CA_location().add(aa6x59->get_CA_location()); //.add(aa4x60->get_CA_location());
        pt.scale(pt.magnitude()/2);
        m.move(pt);
        m.movability = MOV_PINNED;

        Molecule* mols[512];
        AminoAcid** can_clash = p.get_residues_can_clash(aa6x59->get_residue_no());
        l = 0;

        mols[l++] = (Molecule*)aa6x59;
        aa6x59->movability = MOV_FLEXONLY;
        mols[l++] = &m;

        for (i=0; can_clash[i]; i++)
        {
            can_clash[i]->movability = MOV_FLEXONLY;
            mols[l++] = (Molecule*)can_clash[i];
            if (l >= 32) break;
        }
        mols[l] = nullptr;
        Molecule::conform_molecules(mols, 20);

        mols[3] = nullptr;
        Molecule::conform_molecules(mols, 15);
        
        /*fp = fopen("tmp/dbg6x59.pdb", "wb");
        if (!fp) return -3;
        p.save_pdb(fp, &m);
        p.end_pdb(fp);
        fclose(fp);*/
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Save output file.
    ////////////////////////////////////////////////////////////////////////////////

    save_and_exit:
    if (preserve6x55) pose6x55.restore_state(aa6x55);

    std::string out_filename = path + orid + (std::string)".active.pdb";
    fp = fopen(out_filename.c_str(), "wb");
    if (!fp) return -3;
    p.save_pdb(fp);
    p.end_pdb(fp);
    fclose(fp);
    cout << "Saved " << out_filename << endl;


    ////////////////////////////////////////////////////////////////////////////////
    // Save parameters file.
    ////////////////////////////////////////////////////////////////////////////////

    std::string cns_filename = path + orid + (std::string)".params";
    fp = fopen(cns_filename.c_str(), "wb");
    if (!fp) return -3;
    n = constraints.size();
    for (l=0; l<n; l++) fprintf(fp, "%s\n", constraints[l].c_str());
    fclose(fp);
}

