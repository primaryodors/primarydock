#include "../classes/intera.h"

using namespace std;

int main(int argc, char** argv)
{
    Atom F("F"), Cl("Cl"), Br("Br"), I("I");
    Atom O("O"), S("S"), Se("Se");
    Atom N("N"), P("P");
    Atom K("K"), Na("Na"), Li("Li"), Ca("Ca"), Mg("Mg"), Zn("Zn"), Cd("Cd"), Fe("Fe"), Cu("Cu"), Ag("Ag"), Pb("Pb");

    cout << "Metal compatibility test:" << endl;

    cout << "K ... F: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&K, &F)) << endl;
    cout << "K ... Cl: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&K, &Cl)) << endl;
    cout << "K ... Br: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&K, &Br)) << endl;
    cout << "K ... I: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&K, &I)) << endl;
    cout << "K ... O: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&K, &O)) << endl;
    cout << "K ... N: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&K, &N)) << endl;
    cout << "K ... S: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&K, &S)) << endl;

    cout << "Na ... F: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Na, &F)) << endl;
    cout << "Na ... Cl: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Na, &Cl)) << endl;
    cout << "Na ... Br: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Na, &Br)) << endl;
    cout << "Na ... I: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Na, &I)) << endl;
    cout << "Na ... O: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Na, &O)) << endl;
    cout << "Na ... N: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Na, &N)) << endl;
    cout << "Na ... S: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Na, &S)) << endl;

    cout << "Li ... F: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Li, &F)) << endl;
    cout << "Li ... Cl: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Li, &Cl)) << endl;
    cout << "Li ... Br: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Li, &Br)) << endl;
    cout << "Li ... I: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Li, &I)) << endl;
    cout << "Li ... O: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Li, &O)) << endl;
    cout << "Li ... N: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Li, &N)) << endl;
    cout << "Li ... S: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Li, &S)) << endl;

    cout << "Ca ... Cl: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Ca, &Cl)) << endl;
    cout << "Ca ... Br: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Ca, &Br)) << endl;
    cout << "Ca ... I: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Ca, &I)) << endl;
    cout << "Ca ... O: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Ca, &O)) << endl;
    cout << "Ca ... N: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Ca, &N)) << endl;
    cout << "Ca ... S: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Ca, &S)) << endl;

    cout << "Mg ... Cl: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Mg, &Cl)) << endl;
    cout << "Mg ... Br: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Mg, &Br)) << endl;
    cout << "Mg ... I: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Mg, &I)) << endl;
    cout << "Mg ... O: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Mg, &O)) << endl;
    cout << "Mg ... N: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Mg, &N)) << endl;
    cout << "Mg ... P: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Mg, &P)) << endl;
    cout << "Mg ... S: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Mg, &S)) << endl;

    cout << "Zn ... Cl: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Zn, &Cl)) << endl;
    cout << "Zn ... Br: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Zn, &Br)) << endl;
    cout << "Zn ... I: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Zn, &I)) << endl;
    cout << "Zn ... O: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Zn, &O)) << endl;
    cout << "Zn ... N: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Zn, &N)) << endl;
    cout << "Zn ... P: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Zn, &P)) << endl;
    cout << "Zn ... S: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Zn, &S)) << endl;
    cout << "Zn ... Se: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Zn, &Se)) << endl;

    cout << "Cd ... Br: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Cd, &Br)) << endl;
    cout << "Cd ... I: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Cd, &I)) << endl;
    cout << "Cd ... O: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Cd, &O)) << endl;
    cout << "Cd ... N: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Cd, &N)) << endl;
    cout << "Cd ... P: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Cd, &P)) << endl;
    cout << "Cd ... S: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Cd, &S)) << endl;
    cout << "Cd ... Se: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Cd, &Se)) << endl;
    cout << "Cd ... Br: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Cd, &Br)) << endl;

    cout << "Fe ... O: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Fe, &O)) << endl;
    cout << "Fe ... N: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Fe, &N)) << endl;
    cout << "Fe ... P: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Fe, &P)) << endl;
    cout << "Fe ... S: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Fe, &S)) << endl;
    cout << "Fe ... Se: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Fe, &Se)) << endl;

    cout << "Cu ... O: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Cu, &O)) << endl;
    cout << "Cu ... N: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Cu, &N)) << endl;
    cout << "Cu ... P: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Cu, &P)) << endl;
    cout << "Cu ... S: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Cu, &S)) << endl;
    cout << "Cu ... Se: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Cu, &Se)) << endl;

    cout << "Ag ... O: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Ag, &O)) << endl;
    cout << "Ag ... N: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Ag, &N)) << endl;
    cout << "Ag ... P: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Ag, &P)) << endl;
    cout << "Ag ... S: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Ag, &S)) << endl;
    cout << "Ag ... Se: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Ag, &Se)) << endl;

    cout << "Pb ... O: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Pb, &O)) << endl;
    cout << "Pb ... N: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Pb, &N)) << endl;
    cout << "Pb ... P: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Pb, &P)) << endl;
    cout << "Pb ... S: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Pb, &S)) << endl;
    cout << "Pb ... Se: " << 0.01*round(100.0*InteratomicForce::metal_compatibility(&Pb, &Se)) << endl;

    return 0;
}
