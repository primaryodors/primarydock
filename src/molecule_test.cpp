
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "classes/molecule.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule m("nothing");
    cout << "Created empty molecule named " << m.get_name() << ".\n";

    char buffer[65536];
    char tstname[1024] = "BZN.sdf";

    if (argc > 1) strcpy(tstname, argv[1]);

    Molecule m1(tstname);
    cout << "Created molecule named " << m1.get_name() << ".\n";
    int nloaded;

    FILE* pf = fopen(tstname, "rb");
    if (!pf)
    {
        m1.from_smiles(argv[1]);
    }
    else
    {
        fread(buffer, 1, 65535, pf);
        fclose(pf);
        cout << "Read data.\n";
        nloaded = m1.from_sdf(buffer);
    }


    cout << "Loaded " << nloaded << " atoms of " << tstname << " into molecule.\n";

    int rc = m1.get_num_rings();
    if (rc) cout << "Found " << rc << " ring(s)." << endl;

    int i;
    for (i=0; i<rc; i++)
    {
        cout << "Ring atoms: ";
        Atom** ra = m1.get_ring_atoms(i);
        Atom::dump_array(ra);
        cout << endl;
        delete[] ra;

        bool cp = false;
        if (m1.ring_is_aromatic(i))
        {
            cout << "Ring " << i << " is aromatic." << endl;
            cp = true;
        }
        else if (m1.ring_is_coplanar(i))
        {
            cout << "Ring " << i << " is coplanar." << endl;
            cp = true;
        }

        Point rc = m1.get_ring_center(i);
        cout << "Ring centered at [" << rc.x << "," << rc.y << "," << rc.z << "]." << endl;

        if (cp)
        {
            Vector normal = m1.get_ring_normal(i);
            cout << "Ring normal: φ=" << (normal.phi * 180.0/M_PI) << "° θ=" << (normal.theta * 180.0/M_PI) << "°." << endl;
        }
    }

    cout << "Internal clashes: " << m1.get_internal_clashes() << " cu. A." << endl;

    char** anames = m1.get_atom_names();
    if (anames)
    {
        //cout << "Name of first atom " << anames[i] << endl << endl;
        for (i=0; anames[i]; i++)
        {
            Atom* a = m1.get_atom(anames[i]);
            if (a)
            {
                float ab = a->get_acidbase();
                if (ab) cout << "Atom " << anames[i] << " has basicity " << ab << endl;
            }
        }
        delete[] anames;
    }
    else cout << "No names.\n";

    Bond** b = m1.get_rotatable_bonds();
    if (b)
    {
        for (i=0; b[i]; i++)
        {
            cout << "Bond " << b[i]->atom->name << " - " << b[i]->btom->name << " can rotate." << endl;
            b[i]->rotate(30.0 * M_PI / 180);
        }
    }

    Point pt1(0.2,-0.5,-0.3);
    m1.recenter(pt1);
    Point loc = m1.get_barycenter();
    cout << "Molecule moved to [" << loc.x << "," << loc.y << "," << loc.z << "]." << endl;

    // const char* tstlig = "BZN.sdf";
    Molecule m2("test");
    /*
    pf = fopen(tstlig, "rb");
    fread(buffer, 1, 65535, pf);
    fclose(pf);
    m2.from_sdf(buffer);
    */

    if (argc > 2)
    {
        m2.from_smiles(argv[2]);
    }
    else
    {
        Atom* C1 = m2.add_atom("C", "C1", 0, 0);
        strcpy(buffer, C1->get_location().printable());
        cout << "Added a carbon atom. Its location is " << buffer  << "." << endl;

        Atom* C2 = m2.add_atom("C", "C2", C1, 1);
        strcpy(buffer, C2->get_location().printable());
        cout << "Added another carbon atom. Its location is " << buffer << "." << endl;

        Atom* O3 = m2.add_atom("O", "O3", C2, 1);
        strcpy(buffer, O3->get_location().printable());
        cout << "Added an oxygen atom. Its location is " << buffer << "." << endl;

        m2.hydrogenate();
    }

    cout << "Loaded test ligand. Intermol clashes: " << m1.get_intermol_clashes(&m2) << " cu. A." << endl;



    Vector v1(&pt1);
    float rotdeg = -30;
    m2.rotate(&v1, rotdeg * M_PI/180);
    cout << "Rotated molecule 2 by " << rotdeg << " degrees. Intermol clashes: " << m1.get_intermol_clashes(&m2) << " cu. A." << endl;

    Point pt(0,0,1.0);
    Vector v(&pt);
    float ttlmv = 0;
    while (m1.get_intermol_clashes(&m2))
    {
        m2.move(v);
        ttlmv += v.r;
        cout << v.r << " ";
    }
    cout << "Moved molecule 2 by " << ttlmv << " A." << endl;

    cout << "Intermol clashes: " << m1.get_intermol_clashes(&m2) << " cu. A." << endl;
    cout << "Intermol energy level: " << m1.get_intermol_binding(&m2) << " kJ/mol." << endl;

    m1.reset_conformer_momenta();
    m2.reset_conformer_momenta();
    Molecule* mols[3];
    mols[0] = &m1;
    mols[1] = &m2;
    mols[3] = NULL;
    Molecule::multimol_conform(mols, 50);
    float energyLevel = m1.get_intermol_binding(&m2);
    cout << "\n# Post-conformation intermol energy level: " << energyLevel << " kJ/mol." << endl;
    const float energyLevelThreshold = -1.0;
    if(energyLevel > energyLevelThreshold)
        cout << "Energy level above threshold, SUCCESS.\n";
    else
        cout << "Energy level below threshold, FAIL.\n";
 
    const char* tstoutf = "output.sdf";
    pf = fopen(tstoutf, "wb");
    Molecule* ligands[2];
    ligands[0] = &m2;
    ligands[1] = 0;
    if (m1.save_sdf(pf, ligands)) cout << "Saved " << tstoutf << endl;
    else cout << "Failed to save " << tstoutf << endl;
    fclose(pf);

    const char* tstoutf1 = "output.pdb";
    pf = fopen(tstoutf1, "wb");
    m1.save_pdb(pf);
    fclose(pf);



    return 0;
}













