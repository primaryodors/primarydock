
#include "cavity.h"

int Cavity::scan_in_protein(Protein* p, Cavity* cavs, int cmax)
{
    if (!p || !cavs) return 0;
    if (cmax < 1) return 0;

    int i, j, l, sr = p->get_start_resno(), er = p->get_end_resno();
    float x, y, z, step, yoff = 0, zoff = 0;
    Point pcen = p->get_region_center(sr, er), pbox = p->get_region_bounds(sr, er);
    step = cav_xyz_step;
    AminoAcid* can_clash[SPHREACH_MAX+4];
    Molecule dummy("DUMMY");
    dummy.add_atom("H", "H1", nullptr, 0);
    Point size(_INTERA_R_CUTOFF/1.5, _INTERA_R_CUTOFF/1.5, _INTERA_R_CUTOFF/1.5);
    CPartial parts[65536];
    j=0;
    bool any_priority = false;

    cout << "Cavity search in progress..." << flush;
    for (x = pcen.x - pbox.x + min_dist_bounding_box; x <= pcen.x + pbox.x - min_dist_bounding_box; x += step)
    {
        cout << "." << flush;
        yoff = yoff ? 0 : step/2;
        for (y = pcen.y - pbox.y + min_dist_bounding_box - yoff; y <= pcen.y + pbox.y - min_dist_bounding_box; y += step)
        {
            zoff = zoff ? 0 : step/2;
            for (z = pcen.z - pbox.z + min_dist_bounding_box - zoff; z <= pcen.z + pbox.z - min_dist_bounding_box; z += step)
            {
                Point pt(x,y,z);
                dummy.recenter(pt);
                int sphres = p->get_residues_can_clash_ligand(can_clash, &dummy, pt, size, nullptr, true);
                if (sphres < 8) continue;          // Too isolated.
                float rmin;
                CPartial working;
                for (i=0; i<sphres; i++)
                {
                    Atom* a = can_clash[i]->get_nearest_atom(pt);
                    float r = a->get_location().get_3d_distance(pt);
                    if (r < rmin || !i)
                    {
                        rmin = r;
                        working.resno = a->residue;
                    }
                    if (r < min_partial_radius) break;

                    Atom* CB = can_clash[i]->get_atom("CB");
                    if (CB && a->get_location().get_3d_distance(CB->get_location()) 
                        < a->get_location().get_3d_distance(can_clash[i]->get_CA_location()))
                    {
                        if (can_clash[i]->priority)
                        {
                            // cout << can_clash[i]->get_name() << " has priority." << endl;
                            any_priority = working.priority = true;
                        }
                        if (can_clash[i]->coordmtl) working.metallic = true;
                        if (can_clash[i]->get_charge() < -hydrophilicity_cutoff) working.chargedn = true;
                        if (can_clash[i]->get_charge() >  hydrophilicity_cutoff) working.chargedp = true;
                        if (can_clash[i]->pi_stackability() >= 0.2) working.pi = true;
                        if (fabs(can_clash[i]->hydrophilicity()) > hydrophilicity_cutoff) working.polar = true;
                        if (can_clash[i]->count_atoms_by_element("S")) working.thio = true;
                    }
                }

                if (rmin >= min_partial_radius)
                {
                    working.s.center = pt;
                    working.s.radius = rmin;
                    // cout << "Found partial at " << pt << " radius " << rmin << endl << flush;
                    parts[j++] = working;
                }
            }
        }
    }
    cout << endl;

    // Now consolidate all partials into glommed cavities.
    Cavity tmpcav[4096];
    int pmax = j;
    j=0;
    for (i=0; i<pmax; i++)
    {
        if (parts[i].s.center.magnitude() == 0) break;
        bool glommed = false;
        for (l=0; l<j; l++)
        {
            float inter = tmpcav[l].partial_intersects_cavity(parts[i]);
            if (/*(parts[i].priority && tmpcav[l].priority) ||*/ inter >= cav_linking_threshold)
            {
                tmpcav[l].add_partial(parts[i]);
                glommed = true;
                /*cout << "Partial at " << parts[i].s.center << " radius " << parts[i].s.radius << " belongs to cavity #" << l
                    << " with intersection " << inter << endl << flush;*/
                break;
            }
        }
        if (!glommed)
        {
            tmpcav[j++].add_partial(parts[i]);
        }
        if (j >= 4090) break;
    }

    // if (any_priority) cout << "Priority residues found." << endl;
    l=j;
    j=0;
    for (i=0; i<l; i++)
    {
        if (tmpcav[i].count_partials() >= cav_min_partials
            && (!any_priority || tmpcav[i].priority)
            ) cavs[j++] = tmpcav[i];
        if (j >= cmax-1) break;
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
        // float inter = sphere_intersection(partials[i].s.radius, p.s.radius, p.s.center.get_3d_distance(partials[i].s.center));
        float inter = partials[i].s.radius + p.s.radius - p.s.center.get_3d_distance(partials[i].s.center);
        // inter /= p.s.volume();
        if (inter > result) result = inter;
        // result += inter;
    }

    return result;
}

void Cavity::add_partial(CPartial p)
{
    if (!pallocd)
    {
        pallocd = 256;
        partials = new CPartial[pallocd+4];
        priority = false;
    }

    if (!p.s.center.x && !p.s.center.y && !p.s.center.z) return;

    int i, j;
    for (i=0; i<pallocd; i++) if (partials[i].s.radius < min_partial_radius) break;             // Get count.

    if (i >= pallocd-4)
    {
        CPartial* old = partials;
        partials = new CPartial[pallocd+260];
        for (j=0; j<i; j++) partials[j] = old[j];
        delete[] old;
        pallocd += 256;
    }

    partials[i] = p;
    partials[i+1].s.radius = 0;
    if (p.priority) this->priority = true;
    // cout << "Cavity " << this << " added partial at " << partials[i].s.center << " radius " << partials[i].s.radius << endl << flush;
}

void Cavity::output_ngl_js(FILE* fp)
{
    if (!fp) return;
    if (!pallocd) return;

    int i;
    if (priority) fprintf(fp, "// PRIORITY:\n");
    fprintf(fp, "var shape = new NGL.Shape( \"shape\" );\n");
    fprintf(fp, "var sphereBuffer = new NGL.SphereBuffer(\n{\n");
    fprintf(fp, "    position: new Float32Array( [ ");
    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        if (i) fprintf(fp, ", ");
        fprintf(fp, "%f, %f, %f", partials[i].s.center.x, partials[i].s.center.y, partials[i].s.center.z);
        if (!(i%5)) fprintf(fp, "\n\t\t");
    }
    fprintf(fp, " ] ),\n");
    fprintf(fp, "    color: new Float32Array( [ ");
    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        if (i) fprintf(fp, ", ");
        if (partials[i].priority) fprintf(fp, "0.0, 1.0, 0.0");
        else if (partials[i].metallic) fprintf(fp, "0.7, 0.5, 0.3");
        else if (partials[i].chargedp && partials[i].chargedn) fprintf(fp, "1, 0, 1");
        else if (partials[i].chargedp) fprintf(fp, "0.1, 0.1, 1");
        else if (partials[i].chargedn) fprintf(fp, "1, 0.1, 0.1");
        else if (partials[i].polar) fprintf(fp, "0.1, 1, 1");
        else if (partials[i].thio) fprintf(fp, "1, 0.8, 0.1");
        else if (partials[i].pi) fprintf(fp, "0.8, 0.6, 0.8");
        else fprintf(fp, "0.6, 0.6, 0.6");
        if (!(i%5)) fprintf(fp, "\n\t\t");
    }
    fprintf(fp, " ] ),\n");
    fprintf(fp, "    radius: new Float32Array( [ ");
    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        if (i) fprintf(fp, ", ");
        fprintf(fp, "%f", partials[i].s.radius);
        if (!(i%5)) fprintf(fp, "\n\t\t");
    }
    fprintf(fp, " ] ),\n");
    fprintf(fp, "} );\n");
    fprintf(fp, "shape.addBuffer( sphereBuffer );\n");
    fprintf(fp, "var shapeComp = stage.addComponentFromObject( shape );\n");
    fprintf(fp, "shapeComp.addRepresentation( \"buffer\", { opacity: 0.3 } );\n");
    // fprintf(fp, "shapeComp.autoView();\n");
    fprintf(fp, "\n");
}

int Cavity::count_partials()
{
    if (!pallocd) return 0;
    int i;
    for (i=0; i<pallocd; i++) if (partials[i].s.radius < min_partial_radius) return i;
    return 0;
}

Point Cavity::get_center()
{
    if (!pallocd) return Point(0,0,0);
    int i, j=0;
    Point foravg[pallocd+4];
    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        foravg[j] = partials[i].s.center;
        foravg[j].weight = partials[i].s.radius;
        j++;
    }

    return average_of_points(foravg, j);
}

Point Cavity::nearest_surface_vertex(Point pt)
{
    if (!vdw_vertex_count) return Point(0,0,0);
    int i;
    float bestr;
    Point result;
    for (i=0; i<vdw_vertex_count; i++)
    {
        float r = pt.get_3d_distance(vdw_surface[i]);
        if (!i || r < bestr)
        {
            bestr = r;
            result = vdw_surface[i];
        }
    }

    return result;
}

CPartial* Cavity::get_nearest_partial(Point pt)
{
    if (!pallocd || !partials) return nullptr;

    int i;
    float bestr;
    CPartial* result;
    for (i=0; i<pallocd && partials[i].s.radius; i++)
    {
        float r = partials[i].s.center.get_3d_distance(pt);
        if (!i || r < bestr)
        {
            bestr = r;
            result = &partials[i];
        }
    }

    return result;
}

void Cavity::compute_vdW_surface(float d)
{
    if (!pallocd || !partials) return;

    int maxpoints = pallocd * d * d / 3 + 256;
    if (!vdw_surface)
    {
        vdw_surface = new Point[maxpoints];
        vdw_vertex_partial = new CPartial*[maxpoints];
    }

    float halfstep = M_PI / d;
    float step = halfstep * 2;

    int i, ivdW = 0;
    Vector v;
    for (i=0; i<pallocd && partials[i].s.radius; i++)
    {
        Point ploc = partials[i].s.center;
        v.r = partials[i].s.radius;
        float ystep = step / v.r / v.r;
        for (v.theta = -square; v.theta <= square; v.theta += step)
        {
            float xstep = step / v.r / fmax(cos(v.theta), 0.000001);
            float end = M_PI*2-xstep/2;
            for (v.phi = 0; v.phi < end; v.phi += xstep)
            {
                Point pt = ploc.add(v);
                CPartial* np = this->get_nearest_partial(pt);
                if (np != &partials[i] && pt.get_3d_distance(np->s.center) < np->s.radius) continue;
                if (!pt.x && !pt.y && !pt.z) pt = Point(-0.001, 0.001, -0.001);
                vdw_vertex_partial[ivdW] = &partials[i];
                vdw_surface[ivdW++] = pt;
                if (ivdW >= maxpoints)
                {
                    cout << "Too many cavity surface ligand_vertices. Please increase limit in code." << endl;
                    throw 0xbadc0de;
                }
            }
        }
    }
    vdw_vertex_count = ivdW;
}

CPartial* Cavity::point_inside_pocket(Point pt)
{
    if (!pallocd) return nullptr;
    int i;
    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        float r = partials[i].s.center.get_3d_distance(pt);
        if (r < partials[i].s.radius) return &partials[i];
    }

    return nullptr;
}

const Point* ligand_vertices;
float Cavity::containment_violations(Molecule* m, float simt)
{
    if (!ligand_vertices) return 0;
    Atom** va = m->get_vdW_vertex_atoms();

    int i;
    float viol = m->total_eclipses()*33;
    for (i=0; ligand_vertices[i].x || ligand_vertices[i].y || ligand_vertices[i].z; i++)
    {
        CPartial* cp = point_inside_pocket(ligand_vertices[i]);
        int Z = va[i]->get_Z();
        if (!cp)
        {
            viol += (Z > 1) ? 1 : 0.5;

            if (viol > simt) return viol+simt;
        }
        else
        {
            float coefficient = cp->priority ? 0.05 : 0.001;
            if (cp->metallic && Z==16) viol -= 1.00*coefficient;
            if (cp->chargedn && va[i]->get_charge() >  hydrophilicity_cutoff) viol -= .50*coefficient;
            if (cp->chargedp && va[i]->get_charge() < -hydrophilicity_cutoff) viol -= .50*coefficient;
            if (cp->polar && fabs(va[i]->is_polar()) > hydrophilicity_cutoff) viol -= .25*coefficient;
        }
    }

    return viol;
}

float Cavity::find_best_containment(Molecule* m)
{
    ligand_vertices = m->obtain_vdW_surface(10);
    Point cen = get_center();
    m->recenter(cen);
    Pose best(m);
    float bestviol = containment_violations(m);

    Point axes[3];
    axes[0] = Point(1,0,0);
    axes[1] = Point(0,1,0);
    axes[2] = Point(0,0,1);
    int i, j, k, l, n = m->get_atom_count();
    float theta, besttheta;
    for (l=0; l<5; l++)
    {
        for (j=0; j<3; j++)
        {
            LocatedVector lv = (Vector)axes[j];
            lv.origin = cen;
            besttheta=0;
            for (theta=0; theta < M_PI*2; theta += cav_360_step)
            {
                int atoms_outside_cavity = 0;
                std::string ldbg = "";
                for (i=0; i<n; i++)
                {
                    Atom* a = m->get_atom(i);
                    if (a->get_Z() < 2) continue;
                    Point aloc = a->get_location();
                    CPartial* inside = point_inside_pocket(aloc);
                    if (!inside) atoms_outside_cavity++;
                    else ldbg += (std::string)"Atom " + (std::string)a->name + (std::string)" is inside partial centered at "
                        + std::to_string(inside->s.center.x) + std::to_string(inside->s.center.y) + std::to_string(inside->s.center.z)
                        + (std::string)" radius " + std::to_string(inside->s.radius) + (std::string)"\n";
                }

                if (!atoms_outside_cavity)
                {
                    float viol = containment_violations(m, fmax(0, bestviol));
                    if (viol < bestviol)
                    {
                        best.copy_state(m);
                        bestviol = viol;
                        besttheta = theta;
                    }
                }

                m->rotate(lv, cav_360_step);
            }
            best.restore_state(m);
        }

        for (j=0; j<=26; j++)
        {
            Point maybe = cen;
            k=0;
            _retry_linear_motion:
            maybe.x += 0.5 * (j%3-1);
            maybe.y += 0.5 * ((j/3)%3-1);
            maybe.z += 0.5 * j/9;

            m->recenter(maybe);
            float viol = containment_violations(m);
            if (viol < bestviol)
            {
                best.copy_state(m);
                bestviol = viol;
                k++;
                if (k < 5) goto _retry_linear_motion;
            }
        }
        best.restore_state(m);
    }

    return bestviol;
}
