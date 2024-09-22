#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../classes/molecule.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule m1("benzene"), m2("benzene too");

    char buffer[65536];
    FILE* pf = fopen("sdf/benzene.sdf", "rb");
    fread(buffer, 1, 65535, pf);
    fclose(pf);
    m1.from_sdf(buffer);
    m2.from_sdf(buffer);

    Molecule* ligands[2];
    ligands[0] = &m2;
    ligands[1] = nullptr;
    float e;

    m1.recenter(Point(0,0,0));
    m2.recenter(Point(0,0,0));


    // Sandwich configuration.
    m2.move(Point(0,0,3.87));
    e = -m1.get_intermol_binding(&m2).summed();
    cout << "Sandwich stack energy: " << e << " kJ/mol." << endl;
    pf = fopen("output/sandwich.sdf", "wb");
    if (!pf) { cout << "FAILED to open output file, please check output dir for write access." << endl; return -1; }
    m1.save_sdf(pf, ligands);
    fclose(pf);


    // Parallel displaced.
    m2.move(Point(0,1.09,0));
    e = -m1.get_intermol_binding(&m2).summed();
    cout << "Parallel displaced stack energy: " << e << " kJ/mol." << endl;
    pf = fopen("output/pdisp.sdf", "wb");
    if (!pf) { cout << "FAILED to open output file, please check output dir for write access." << endl; return -1; }
    m1.save_sdf(pf, ligands);
    fclose(pf);


    // T-stacked.
    m2.recenter(Point(0,0,4.96));
    Vector axis = Point(1,0,0);
    m2.rotate(&axis, M_PI/2, false);
    e = -m1.get_intermol_binding(&m2).summed();
    cout << "T-stack energy: " << e << " kJ/mol." << endl;
    pf = fopen("output/tstack.sdf", "wb");
    if (!pf) { cout << "FAILED to open output file, please check output dir for write access." << endl; return -1; }
    m1.save_sdf(pf, ligands);
    fclose(pf);


    return 0;
}
