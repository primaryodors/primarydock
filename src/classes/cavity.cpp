
#include "cavity.h"

int Cavity::scan_in_protein(Protein* p, Cavity* cavs, int cmax)
{
    if (!p || !cavs) return 0;
    if (cmax < 1) return 0;

    int i, j, l, sr = p->get_start_resno(), er = p->get_end_resno();
    float x, y, z, step;
    Point pcen = p->get_region_center(sr, er), pbox = p->get_region_bounds(sr, er);
    step = pbox.magnitude() / 100;
    AminoAcid* can_clash[SPHREACH_MAX+4];
    Molecule dummy("DUMMY");
    dummy.add_atom("H", "H1", nullptr, 0);
    Point size(_INTERA_R_CUTOFF, _INTERA_R_CUTOFF, _INTERA_R_CUTOFF);
    CPartial parts[65536];
    j=0;

    for (x = pcen.x - pbox.x; x <= pcen.x + pbox.x; x += step)
    {
        for (y = pcen.y - pbox.y; y <= pcen.y + pbox.y; y += step)
        {
            for (z = pcen.z - pbox.z; z <= pcen.z + pbox.z; z += step)
            {
                Point pt(x,y,z);
                dummy.recenter(pt);
                int sphres = p->get_residues_can_clash_ligand(can_clash, &dummy, pt, size, nullptr);
                float rmin;
                CPartial working;
                for (i=0; i<sphres; i++)
                {
                    Atom* a = can_clash[i]->get_nearest_atom(pt);
                    float r = a->get_location().get_3d_distance(pt);
                    if (r < rmin || !i) rmin = r;
                    if (r < min_partial_radius) break;

                    if (can_clash[i]->coordmtl) working.metallic = true;
                    if (can_clash[i]->get_charge() < -hydrophilicity_cutoff) working.chargedn = true;
                    if (can_clash[i]->get_charge() >  hydrophilicity_cutoff) working.chargedp = true;
                    if (can_clash[i]->pi_stackability() >= 0.2) working.pi = true;
                    if (fabs(can_clash[i]->hydrophilicity()) > hydrophilicity_cutoff) working.polar = true;
                    if (can_clash[i]->count_atoms_by_element("S")) working.thio = true;
                }

                if (rmin >= min_partial_radius)
                {
                    working.s.center = pt;
                    working.s.radius = rmin;
                    parts[j++] = working;
                }
            }
        }
    }

    // Now consolidate all partials into glommed cavities.
    int pmax = j;
    j=0;
    for (i=0; i<pmax; i++)
    {
        bool glommed = false;
        for (l=0; l<j; l++)
        {
            float inter = cavs[l].partial_intersects_cavity(parts[i]);
            if (inter > 1)
            {
                cavs[l].add_partial(parts[i]);
                glommed = true;
            }
        }
        if (!glommed)
        {
            cavs[j++].add_partial(parts[i]);
        }
    }

    return j;
}

float Cavity::partial_intersects_cavity(CPartial p)
{
    if (!pallocd) return 0;
    int i;
    float result = 0;
    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        float inter = sphere_intersection(partials[i].s.radius, p.s.radius, p.s.center.get_3d_distance(partials[i].s.center));
        if (inter > result) result = inter;
    }

    return result;
}

void Cavity::add_partial(CPartial p)
{
    if (!pallocd)
    {
        pallocd = 256;
        partials = new CPartial[pallocd+4];
    }

    int i, j;
    for (i=0; i<pallocd; i++) if (partials[i].s.radius < min_partial_radius) break;             // Get count.

    if (i >= pallocd-4)
    {
        CPartial* old = partials;
        partials = new CPartial[pallocd+260];
        for (j=0; j<i; j++) partials[i] = old[i];
        delete[] old;
        pallocd += 256;
    }

    partials[i] = p;
    partials[i+1].s.radius = 0;
}

void Cavity::output_ngl_js(FILE* fp)
{
    if (!fp) return;
    if (!pallocd) return;

    int i;
    fprintf(fp, "var shape = new NGL.Shape( \"shape\" );\n");
    fprintf(fp, "var sphereBuffer = new NGL.SphereBuffer(\n{\n");
    fprintf(fp, "    position: new Float32Array( [ ");
    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        if (i) fprintf(fp, ", ");
        fprintf(fp, "%f, %f, %f", partials[i].s.center.x, partials[i].s.center.y, partials[i].s.center.z);
    }
    fprintf(fp, " ] ),\n");
    fprintf(fp, "    color: new Float32Array( [ ");
    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        if (i) fprintf(fp, ", ");
        if (partials[i].metallic) fprintf(fp, "0.7, 0.5, 0.3");
        else if (partials[i].chargedp && partials[i].chargedn) fprintf(fp, "1, 0, 1");
        else if (partials[i].chargedp) fprintf(fp, "0.1, 0.1, 1");
        else if (partials[i].chargedn) fprintf(fp, "1, 0.1, 0.1");
        else if (partials[i].polar) fprintf(fp, "0.1, 1, 1");
        else if (partials[i].thio) fprintf(fp, "1, 0.8, 0.1");
        else if (partials[i].pi) fprintf(fp, "0.8, 0.6, 0.8");
        else fprintf(fp, "0.6, 0.6, 0.6");
    }
    fprintf(fp, " ] ),\n");
    fprintf(fp, "    radius: new Float32Array( [ ");
    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        if (i) fprintf(fp, ", ");
        fprintf(fp, "%f", partials[i].s.radius);
    }
    fprintf(fp, " ] ),\n");
    fprintf(fp, "} );\n");
    fprintf(fp, "shape.addBuffer( sphereBuffer );\n");
    fprintf(fp, "var shapeComp = stage.addComponentFromObject( shape );\n");
    fprintf(fp, "shapeComp.addRepresentation( \"buffer\" );\n");
    fprintf(fp, "shapeComp.autoView();\n");
    fprintf(fp, "\n");
}
