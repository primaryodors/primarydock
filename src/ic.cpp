
#include "classes/protein.h"

int main(int argc, char** argv)
{
    if (argc < 1) return -1;

    Protein p("testing");
    FILE* fp = fopen(argv[1], "rb");

    if (argc > 2) p.load_pdb(fp, 0, argv[2][0]);
    else p.load_pdb(fp);

    p.set_name_from_pdb_name(argv[1]);
    fclose(fp);

    float threshold = 1.5;
    if (argc > 3) threshold = atof(argv[3]);

    int i, j, n = p.get_end_resno();
    for (i=1; i<=n; i++)
    {
        AminoAcid* aa = p.get_residue(i);
        if (!aa) continue;
        AminoAcid** cl = p.get_residues_can_clash(i);
        if (!cl) continue;

        for (j=0; cl[j]; j++)
        {
            if (cl[j] < aa) continue;
            if (fabs(cl[j]->get_residue_no() - i) < 3) continue;
            // float r = aa->get_CA_location().get_3d_distance(cl[j]->get_CA_location());
            Atom *a, *b;
            aa->Molecule::mutual_closest_atoms(cl[j], &a, &b);
            if (!a || !b) continue;
            float r = a->distance_to(b);
            if (r > _INTERA_R_CUTOFF) continue;

            float f = aa->get_intermol_binding(cl[j], false).summed();
            if (fabs(f) >= threshold) cout << *aa << "-" << *cl[j] << ": r=" << r << "Å" << endl;
        }
    }
}
