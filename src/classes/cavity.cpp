
#include "cavity.h"

float cav_xmax = Avogadro, cav_xmin = -Avogadro, cav_ymax = Avogadro, cav_ymin = -Avogadro, cav_zmax = Avogadro, cav_zmin = -Avogadro;
float cav_xyrlim = Avogadro, cav_xzrlim = Avogadro, cav_yzrlim = Avogadro;
int cav_resmin = -99999, cav_resmax = 99999;

int Cavity::scan_in_protein(Protein* p, Cavity* cavs, int cmax)
{
    if (!p || !cavs) return 0;
    if (cmax < 1) return 0;

    int i, j, l, sr = max(p->get_start_resno(), cav_resmin), er = min(p->get_end_resno(), cav_resmax);
    float x, y, z, step, yoff = 0, zoff = 0;
    Point pcen = p->get_region_center(sr, er), pbox = p->get_region_bounds(sr, er);

    int priorities[1024], pqty;
    priorities[0] = pqty = 0;
    l = p->get_end_resno();
    for (i=1; i<=l; i++)
    {
        AminoAcid* a = p->get_residue(i);
        if (a && a->priority) priorities[pqty++] = i;
    }
    priorities[pqty] = 0;

    step = cav_xyz_step;
    AminoAcid* can_clash[SPHREACH_MAX+4];
    Molecule dummy("DUMMY");
    dummy.add_atom("H", "H1", nullptr, 0);
    Point size(_INTERA_R_CUTOFF, _INTERA_R_CUTOFF, _INTERA_R_CUTOFF);
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
                int sphres = p->get_residues_can_clash_ligand(can_clash, &dummy, pt, size, priorities, true);
                if (sphres < 8+pqty) continue;          // Too isolated.
                float rmin;
                CPartial working;
                for (i=0; i<sphres; i++)
                {
                    Atom* a = can_clash[i]->get_nearest_atom(pt);
                    float r = a->get_location().get_3d_distance(pt);
                    if (r < rmin || !i)
                    {
                        rmin = r;
                    }
                    if (r < min_partial_radius) break;
                }

                for (i=0; i<sphres; i++)
                {
                    Atom* a = can_clash[i]->get_nearest_atom(pt);
                    float r = a->get_location().get_3d_distance(pt);
                    if (r > rmin+1.5 && !can_clash[i]->priority) continue;
                    if (r > rmin+6) continue;

                    Atom* CB = can_clash[i]->get_atom("CB");
                    if (CB
                        && a->get_location().get_3d_distance(CB->get_location()) < a->get_location().get_3d_distance(can_clash[i]->get_CA_location())
                        && (!pqty || can_clash[i]->priority)
                        )
                    {
                        if (can_clash[i]->priority)
                        {
                            // cout << pt << can_clash[i]->get_name() << " distance " << r << " has priority." << endl;
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
        // if (parts[i].resno && (parts[i].resno < sr || parts[i].resno > er)) continue;
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
        Point cen = tmpcav[i].get_center();
        if (cen.x < cav_xmin || cen.x > cav_xmax) continue;
        if (cen.y < cav_ymin || cen.y > cav_ymax) continue;
        if (cen.z < cav_zmin || cen.z > cav_zmax) continue;

        float r = sqrt(pow(cen.x - pcen.x, 2) + pow(cen.y - pcen.y, 2));
        if (r > cav_xyrlim) continue;
        r = sqrt(pow(cen.x - pcen.x, 2) + pow(cen.z - pcen.z, 2));
        if (r > cav_xzrlim) continue;
        r = sqrt(pow(cen.y - pcen.y, 2) + pow(cen.z - pcen.z, 2));
        if (r > cav_yzrlim) continue;

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
    SCoord v;
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

float Cavity::sphere_inside_pocket(Sphere s, CPartial** p)
{
    if (p) *p = nullptr;
    if (!pallocd) return 0;

    int i;
    float result = 0;
    for (i=0; i<pallocd; i++)
    {
        if (partials[i].s.radius < min_partial_radius) break;
        if (partials[i].s.radius < s.radius) break;
        float r = partials[i].s.center.get_3d_distance(s.center);
        if (r < (partials[i].s.radius - s.radius))
        {
            if (p) *p = &partials[i];
            return 1;
        }
        else if (r < partials[i].s.radius)
        {
            // float f = sphere_intersection(partials[i].s.radius, s.radius, r) / sphere_intersection(partials[i].s.radius, s.radius, 0);
            float f = 1.0 - (r - (partials[i].s.radius - s.radius)) / s.radius;
            if (f > result)
            {
                result = f;
                if (p) *p = &partials[i];   
            }
        }
    }

    return result;
}

const Point* ligand_vertices;
float Cavity::containment_violations(Molecule* m, float simt)
{
    int i, n = m->get_atom_count();
    float viol = 0;
    for (i=0; i<n; i++)
    {
        Atom* a = m->get_atom(i);
        float f = sphere_inside_pocket(a->get_sphere());
        int Z = a->get_Z();

        viol += ((Z > 1) ? 1 : 0.5) * (1.0 - f);
        if ((simt >= 0) && (viol > simt)) return viol;
    }

    return viol;
}

float Cavity::find_best_containment(Molecule* m, bool mbt)
{
    ligand_vertices = m->obtain_vdW_surface(10);
    Point cen = get_center();
    m->recenter(cen);
    Pose best(m);
    float bestviol = Avogadro;

    Point axes[3];
    axes[0] = Point(1,0,0);
    axes[1] = Point(0,1,0);
    axes[2] = Point(0,0,1);
    int i, j, k, l, n = m->get_atom_count();
    float theta, besttheta;
    for (l=0; l<15; l++)
    {
        for (j=0; j<3; j++)
        {
            LocatedVector lv = (SCoord)axes[j];
            lv.origin = cen;
            besttheta=0;
            for (theta=0; theta < M_PI*2; theta += cav_360_step)
            {
                float atoms_outside_cavity = 0;
                std::string ldbg = "";
                for (i=0; i<n; i++)
                {
                    Atom* a = m->get_atom(i);
                    if (a->get_Z() < 2) continue;
                    Point aloc = a->get_location();

                    CPartial* inside;
                    float f = sphere_inside_pocket(a->get_sphere(), &inside);
                    if (mbt && inside)
                    {
                        float e = 0;
                        if (inside->chargedp && a->is_conjugated_to_charge() < -hydrophilicity_cutoff) e = 0.6;
                        else if (inside->chargedn && a->is_conjugated_to_charge() > hydrophilicity_cutoff) e = 0.6;
                        else if (inside->metallic && (a->get_family() == CHALCOGEN || a->get_family() == PNICTOGEN) && a->get_Z() != 8) e = 0.8;
                        else if (inside->polar && fabs(a->is_polar()) > hydrophilicity_cutoff) e = 0.3;
                        else if (inside->thio && a->get_family() == CHALCOGEN && a->get_Z() != 8) e = 0.15;
                        else if (inside->pi && a->is_pi()) e = 0.12;
                        if (inside->priority) e *= 2;
                        f += e;
                    }
                    atoms_outside_cavity += (1.0-f);
                }

                if (atoms_outside_cavity < bestviol)
                {
                    best.copy_state(m);
                    bestviol = atoms_outside_cavity;
                    besttheta = theta;
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
            float viol = containment_violations(m) + m->total_eclipses()*33;
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

std::string Cavity::resnos_as_string(Protein* p)
{
    if (!partials || !pallocd) return (std::string)"";
    int i, j, n=p->get_end_resno();
    bool resincluded[n+4];
    for (i=0; i<=n; i++) resincluded[i] = false;
    for (i=0; i<pallocd; i++)
    {
        if (!partials[i].s.radius) break;
        std::string pres = partials[i].resnos_as_string(p);
        char buffer[1024];
        strcpy(buffer, pres.c_str());
        char** words = chop_spaced_words(buffer);
        if (!words) continue;
        for (j=0; words[j]; j++) resincluded[atoi(words[j])] = true;
    }

    std::string result = "";
    for (i=1; i<=n; i++)
    {
        if (resincluded[i])
        {
            if (result.length()) result += " ";
            result += std::to_string(i);
        }
    }

    return result;
}

std::string CPartial::resnos_as_string(Protein* p)
{
    int i, j, n = p->get_end_resno();
    bool intersect[n+4];
    for (i=1; i<=n; i++)
    {
        AminoAcid* aa = p->get_residue(i);
        if (!aa)
        {
            intersect[i] = false;
            continue;
        }
        j = aa->atoms_inside_sphere(s, nullptr, aa->priority ? 4.0 : 1.1);
        intersect[i] = (j>0);
    }

    std::string result;
    for (i=1; i<=n; i++) if (intersect[i]) result += (std::string)(result.length() ? " " : "") + std::to_string(i);
    return result;
}

int CPartial::from_cvty_line(char* lndata)
{
    int cno;

    //           1111111111222222222233333333334444444444
    // 01234567890123456789012345678901234567890123456789
    //    2   -4.228   22.449    7.041   2.569  -+HSP! 96 99 157 158 161 182

    lndata[4] = 0;
    cno = atoi(lndata);
    lndata[13] = 0;
    s.center.x = atof(lndata+5);
    lndata[22] = 0;
    s.center.y = atof(lndata+14);
    lndata[31] = 0;
    s.center.z = atof(lndata+23);
    lndata[39] = 0;
    s.radius = atof(lndata+32);
    chargedn = (lndata[41] == '-');
    chargedp = (lndata[42] == '+');
    polar    = (lndata[43] == 'H');
    thio     = (lndata[44] == 'S');
    pi       = (lndata[45] == 'P');
    priority = (lndata[46] == '!');

    return cno;
}
