
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/stat.h>
#include "classes/group.h"

using namespace std;

int main(int argc, char** argv)
{
    Protein p("TheProtein");
    FILE* fp;

    fp = fopen("pdbs/OR1/OR1E1.upright.pdb", "rb");
    if (!fp) throw 0xffff;
    p.load_pdb(fp);
    fclose(fp);

    AminoAcid* aa109 = p.get_residue(109);
    AminoAcid* aa155 = p.get_residue(155);
    AminoAcid* aa159 = p.get_residue(159);

    float hbond = aa155->get_intermol_binding(aa109);
    if (hbond < 2 && aa155->hisflips)
    {
        aa155->do_histidine_flip(aa155->hisflips[0]);
        hbond = aa155->get_intermol_binding(aa109);
    }
    cout << aa109->get_name() << "-" << aa155->get_name() << " energy: " << -hbond << endl;
}
