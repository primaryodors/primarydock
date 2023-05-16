
#include "glom.h"

std::vector<int> extra_wt;
std::vector<MCoord> mtlcoords;

void ResiduePlaceholder::set(const char* str)
{
    if (strchr(str, '.')) bw = str;
    else resno = atoi(str);
}

void ResiduePlaceholder::resolve_resno(Protein* prot)
{
    int hxno = atoi(bw.c_str());
    const char* dot = strchr(bw.c_str(), '.');
    int bwpos = atoi(dot+1);
    resno = prot->get_bw50(hxno) + bwpos - 50;
}

Point AtomGlom::get_center()
{
    int i;
    float mass = 0;
    Point result(0,0,0);
    int atct = atoms.size();
    if (!atct) return Point(0,0,0);
    for (i=0; i<atct; i++)
    {
        Point pt = atoms[i]->get_location();
        float m = atoms[i]->get_atomic_weight();
        pt.scale(pt.magnitude() * m);
        result = result.add(pt);
        mass += m;
    }
    if (mass) result.scale(result.magnitude() / mass);
    return result;
}

float AtomGlom::get_pi()
{
    int atct = atoms.size();
    if (!atct) return 0;
    int i;
    float result = 0;
    for (i=0; i<atct; i++)
    {
        if (atoms[i]->is_pi()) result += 1;
    }
    return result;
}

float AtomGlom::get_polarity()
{
    int atct = atoms.size();
    if (!atct) return 0;
    int i;
    float result = 0;
    for (i=0; i<atct; i++)
    {
        result += fabs(atoms[i]->is_polar());
    }
    return result;
}

float AtomGlom::get_ionic()
{
    int atct = atoms.size();
    if (!atct) return 0;
    int i;
    float result = 0;
    for (i=0; i<atct; i++)
    {
        float c = atoms[i]->get_charge();
        if (c) result += c;
        else
        {
            if (atoms[i]->get_family() == PNICTOGEN && !atoms[i]->is_amide())
            {
                result += 0.5;
            }
        }
    }
    return result;
}

float AtomGlom::get_mcoord()
{
    int atct = atoms.size();
    if (!atct) return 0;
    int i;
    float result = 0;
    for (i=0; i<atct; i++)
    {
        if (atoms[i]->get_family() == PNICTOGEN && atoms[i]->get_charge() <= 0) result += 1;
        if (atoms[i]->get_family() == CHALCOGEN) result += 1;
    }
    return result;
}

float AtomGlom::get_sum()
{
    float retval = fabs(get_ionic()*60) + get_polarity()*25 + get_pi()*2;
    if (mtlcoords.size()) retval += get_mcoord()*60 - get_polarity()*40;
    return retval;
}

float AtomGlom::get_avg_elecn()
{
    int atct = atoms.size();
    if (!atct) return 0;
    int i, j=0;
    float result = 0;
    for (i=0; i<atct; i++)
    {
        int aif = atoms[i]->get_family();
        if (aif == PNICTOGEN || aif == CHALCOGEN || (aif == HALOGEN && atoms[i]->is_polar()))
        {
            result += atoms[i]->get_electronegativity();
            j++;
            #if _dbg_glomsel
            cout << atoms[i]->name << " has electronegativity " << atoms[i]->get_electronegativity() << endl;
            #endif
        }
    }
    if (j) result /= j;
    return result;
}

int AtomGlom::contains_element(const char* esym)
{
    int i;
    int findZ = Atom::Z_from_esym(esym);
    int atct = atoms.size();
    int result = 0;
    for (i=0; i<atct; i++)
    {
        if (atoms[i]->get_Z() == findZ) result++;
    }
    return result;
}

float AtomGlom::distance_to(Point pt)
{
    int atct = atoms.size();
    if (!atct) return -1;
    int i;
    float result = 0;
    for (i=0; i<atct; i++)
    {
        float f = atoms[i]->get_location().get_3d_distance(pt);
        if (!i || !f || f < result) result = f;
    }
    return result;
}

float AtomGlom::bounds()
{
    int atct = atoms.size();
    if (!atct) return -1;
    int i;
    Point ptmin(0,0,0), ptmax(0,0,0);
    for (i=0; i<atct; i++)
    {
        Point aloc = atoms[i]->get_location();
        if (!i)
        {
            ptmin = aloc;
            ptmax = aloc;
        }
        else
        {
            if (aloc.x < ptmin.x) ptmin.x = aloc.x;
            if (aloc.y < ptmin.y) ptmin.y = aloc.y;
            if (aloc.z < ptmin.z) ptmin.z = aloc.z;
            if (aloc.x > ptmax.x) ptmax.x = aloc.x;
            if (aloc.y > ptmax.y) ptmax.y = aloc.y;
            if (aloc.z > ptmax.z) ptmax.z = aloc.z;
        }
    }

    return ptmax.get_3d_distance(ptmin);
}

float AtomGlom::compatibility(AminoAcid* aa)
{
    int atct = atoms.size();
    if (!atct) return 0;

    float result = 0;
    float lgi = get_ionic(), lgh = get_polarity(), lgp = get_pi();

    float aachg = aa->get_charge();
    if (aa->conditionally_basic()) aachg += 0.5;
    if (lgi && aachg && sgn(lgi) != -sgn(aachg)) return 0;

    if (aa->hydrophilicity() > 0.25)
    {
        if ((lgh / atct) < 0.19) return 0;
    }
    else
    {
        if ((lgh / atct) > 0.33333) return 0;
    }

    if (lgh)
    {
        int i;
        int atct = atoms.size();
        if (!aa->has_hbond_acceptors())
        {
            lgh = 0;
            if (atct)
            for (i=0; i<atct; i++)
            {
                if (atoms[i]->is_polar() < 0.333) lgh++;
            }
        }
        else if (!aa->has_hbond_donors())
        {
            lgh = 0;
            if (atct)
            for (i=0; i<atct; i++)
            {
                if (atoms[i]->is_polar() > 0.333) lgh++;
            }
        }
    }

    if (lgi) result += lgi * -aa->bindability_by_type(ionic) * 100;
    if (lgh) result += fabs(lgh) * fabs(aa->bindability_by_type(hbond)) * 30;
    if (lgp) result += aa->bindability_by_type(pi);

    return result;
}



Point ResidueGlom::get_center()
{
    int amsz = aminos.size();
    if (!amsz) return Point(0,0,0);

    if (metallic) return metal->get_location();

    int i, j;
    j = 0;
    Point result(0,0,0);
    for (i=0; i<amsz; i++)
    {
        Atom** aa = aminos[i]->get_most_bindable(1);
        Atom* a = aa[0]; // = aminos[i]->get_atom("CB");
        delete aa;
        // if (!a) a = aminos[i]->get_atom("CA");      // even though glycine probably shouldn't be part of a glom.
        if (a)
        {
            Point pt = a->get_location();
            result = result.add(pt);
            j++;
        }
    }
    if (j) result.scale(result.magnitude() / j);
    return result;
}

float ResidueGlom::distance_to(Point pt)
{
    return pt.get_3d_distance(get_center());
}

float ResidueGlom::compatibility(AtomGlom* ag)
{
    int amsz = aminos.size();
    if (!amsz) return 0;
    float result = 0;
    int i;
    bool has_acids = false, has_his = false;
    for (i=0; i<amsz; i++)
    {
        if (aminos[i]->get_charge() < 0) has_acids = true;
        if (aminos[i]->conditionally_basic()) has_his = true;

        float f = ag->compatibility(aminos[i]);

        if (extra_wt.size()
                &&
                std::find(extra_wt.begin(), extra_wt.end(), aminos[i]->get_residue_no())!=extra_wt.end()
        )
        {
            f *= 2;		// Extra weight for residues mentioned in a CEN RES or PATH RES parameter.
        }

        result += f;
    }
    if (has_acids && has_his && ag->get_ionic() > 0) result -= ag->get_ionic()*30;

    if (metallic)
    {
        result /= 10;
        float lmc = 0;
        lmc += 20.0 * fabs(ag->get_mcoord());
        // lmc += 1.0 * fmax(0, -ag->get_ionic());
        // lmc -= 1.0 * fabs(ag->get_polarity());

        float kmims = (ag->get_avg_elecn() + metal->get_electronegativity()) / 2 - 2.25;
        float lmm = cos(kmims*2);
        lmc *= lmm;

        #if _dbg_glomsel
        cout << "Metal multiplier: " << lmm << " from kmims " << kmims << endl;
        #endif

        lmc += 1.0 * fabs(ag->get_pi());
        result += lmc;
    }

    return result;
}

float ResidueGlom::glom_reach()
{
    int n = aminos.size();
    if (!n) return 0;
    if (metallic) return 0;
    int i;
    float retval = 0;
    for (i=0; i<n; i++)
    {
        float r = aminos[i]->get_reach();
        if (r > retval) retval = r;
    }

    return retval;
}

std::vector<std::shared_ptr<AtomGlom>> AtomGlom::get_potential_ligand_gloms(Molecule* mol)
{
    std::vector<std::shared_ptr<AtomGlom>> retval;
    if (!mol) return retval;
    int n = mol->get_atom_count();
    if (!n) return retval;

    int i, j, l;
    bool dirty[n+4];
    for (i=0; i<n; i++) dirty[i] = false;

    for (i=0; i<n; i++)
    {
        if (dirty[i]) continue;
        std::shared_ptr<AtomGlom> g(new AtomGlom());
        Atom* a = mol->get_atom(i);
        if (!a) continue;

        if (!a->is_polar() && !a->get_charge() && !a->is_pi()) continue;

        Atom* a_ = a;
        g->atoms.push_back(a);
        #if _dbg_glomsel
        cout << "Creating glom from " << a->name << "..." << endl;
        #endif

        dirty[i] = true;
        for (j=0; j<n; j++)
        {
            if (j==i) continue;

            Atom* b = mol->get_atom(j);
            if (!b) continue;

            float r = fmax(0, g->get_center().get_3d_distance(b->get_location()) - 1.5);
            if (r > 2.5)
            {
                continue;
                #if _dbg_glomsel
                //cout << "Rejected " << b->name << " too far away." << endl;
                #endif
            }
            float simil = fmax(a->similarity_to(b), a_->similarity_to(b));

            if (simil >= 20)
            {
                g->atoms.push_back(b);
                #if _dbg_glomsel
                cout << "Adding " << b->name << " with distance " << r << " and similarity " << simil << endl;
                #endif
                a_ = b;
                dirty[j] = true;
            }
            else
            {
                ;
                #if _dbg_glomsel
                //cout << "Rejected " << b->name << " similarity " << simil << endl;
                #endif
            }
        }

        bool added = false;
        for (l=0; l < retval.size(); l++)
        {
            if (g->get_sum() > retval[l]->get_sum())
            {
                std::vector<std::shared_ptr<AtomGlom>>::iterator it;
                it = retval.begin();
                retval.insert(it+l, g);
                added = true;
                break;
            }
        }
        if (!added) retval.push_back(g);
        #if _dbg_glomsel
        cout << "Glom complete." << endl << endl;
        #endif
    }

    return retval;
}

std::vector<std::shared_ptr<ResidueGlom>> ResidueGlom::get_potential_side_chain_gloms(AminoAcid** aalist, Point pcen)
{
    std::vector<std::shared_ptr<ResidueGlom>> retval;
    if (!aalist) return retval;
    int i, j, n;
    for (n=0; aalist[n]; n++);          // Get count.
    bool dirty[n+4];
    for (i=0; i<n; i++) dirty[i] = false;

    for (i=0; i<n; i++)
    {
        if (dirty[i]) continue;

        std::shared_ptr<ResidueGlom> g(new ResidueGlom());
        AminoAcid* aa = aalist[i];

        Atom* CB = aa->get_atom("CB");
        if (CB)
        {
            float a3d = find_3d_angle(CB->get_location(), pcen, aa->get_CA_location());
            if (a3d > fiftyseventh*120)
            {
                dirty[i] = true;
                continue;
            }
        }

        g->aminos.push_back(aa);
        dirty[i] = true;
        #if _dbg_glomsel
        cout << "Building side chain glom from " << aa->get_name() << endl;
        #endif

        for (j=i+1; j<n; j++)
        {
            if (dirty[j]) continue;
            AminoAcid* bb = aalist[j];

            // Debug trap.
            if (bb->get_residue_no() == 106)
            {
                dirty[j] = false;
            }

            CB = bb->get_atom("CB");
            if (CB)
            {
                float a3d = find_3d_angle(CB->get_location(), pcen, bb->get_CA_location());
                if (a3d > fiftyseventh*120)
                {
                    dirty[j] = true;
                    continue;
                }
            }

            if (aa->coordmtl && aa->coordmtl == bb->coordmtl)
            {
                g->aminos.push_back(bb);
                dirty[j] = true;
                g->metallic = true;
                g->metal = aa->coordmtl;
                #if _dbg_glomsel
                cout << "Adding " << bb->get_name() << " because metal coord." << endl;
                #endif
                continue;
            }
            else if (aa->coordmtl || bb->coordmtl)
            {
                #if _dbg_glomsel
                cout << "Rejected " << bb->get_name() << " metal coord mismatch." << endl;
                #endif
                continue;
            }

            float r = fmax(0, aa->get_CA_location().get_3d_distance(bb->get_CA_location()) - fmax(aa->get_reach(), bb->get_reach()));
            if (r > 2.5)
            {
                #if _dbg_glomsel
                // cout << "Rejected " << bb->get_name() << " distance " << r << endl;
                #endif
                continue;
            }

            float simil = aa->similarity_to(bb);
            if (simil >= 4)
            {
                g->aminos.push_back(bb);
                dirty[j] = true;
                #if _dbg_glomsel
                cout << "Adding " << bb->get_name() << " distance " << r << " similarity " << simil << endl;
                #endif
            }
            else
            {
                ;
                #if _dbg_glomsel
                cout << "Rejected " << bb->get_name() << " similarity " << simil << endl;
                #endif
            }
        }

        #if _dbg_glomsel
        cout << "Completed glom." << endl << endl;
        #endif

        retval.push_back(g);
    }

    return retval;
}

std::ostream& operator<<(std::ostream& os, const AtomGlom& ag)
{
    if (!&ag) return os;
    try
    {
        os << "atom_glom[ ";
        int i;
        for (i=0; i<ag.atoms.size(); i++) os << *ag.atoms[i] << " ";
        os << "]";
    }
    catch (int ex)
    {
        ;
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const ResidueGlom& scg)
{
    if (!&scg) return os;
    try
    {
        os << "residue_glom[ ";
        int i;
        for (i=0; i<scg.aminos.size(); i++) os << *scg.aminos[i] << " ";
        os << "]";
    }
    catch (int ex)
    {
        ;
    }
    return os;
}

float GlomPair::get_potential()
{
    if (potential) return potential;
    else
    {
        int m = ag->atoms.size(), n = scg->aminos.size();
        if (!m || !n) return 0;

        int i, j, q=0;
        for (i=0; i<m; i++)
        {
            Atom* a = ag->atoms[i];
            if (a->get_Z() == 6 && !a->is_polar()) continue;
            for (j=0; j<n; j++)
            {
                AminoAcid* aa = scg->aminos[j];
                float partial;
                if (aa->coordmtl) partial = InteratomicForce::potential_binding(a, aa->coordmtl);
                else
                {
                    partial = aa->get_atom_mol_bind_potential(a);
                    if (fabs(a->is_polar()) > 0.333 && aa->is_tyrosine_like()) partial /= 3;

                    if (extra_wt.size()
                            &&
                            std::find(extra_wt.begin(), extra_wt.end(), aa->get_residue_no())!=extra_wt.end()
                    )
                    {
                        partial *= 2.5;		// Extra weight for residues mentioned in a CEN RES or PATH RES parameter.
                    }

                    #if _dbg_glomsel
                    cout << "Potential for " << *a << "..." << *aa << " = " << partial << endl;
                    #endif
                }
                potential += partial;
                q++;
            }
        }

        if (q) potential /= q;

        float r = pocketcen.get_3d_distance(scg->get_center());
        potential /= fmax(1, r-1.5);

        return potential;
    }
}

std::vector<std::shared_ptr<GlomPair>> GlomPair::pair_gloms(std::vector<std::shared_ptr<AtomGlom>> ag, std::vector<std::shared_ptr<ResidueGlom>> scg, Point pcen)
{
    std::vector<std::shared_ptr<GlomPair>> retval;

    int m = ag.size(), n = scg.size();
    if (!m || !n) return retval;

    int i, j, l;
    bool adirty[m+4], sdirty[n+4];

    for (i=0; i<m; i++) adirty[i] = false;
    for (j=0; j<n; j++) sdirty[j] = false;

    for (i=0; i<m; i++)
    {
        if (adirty[i]) continue;
        int j1 = -1;
        float p = 0;
        for (j=0; j<n; j++)
        {
            if (sdirty[j]) continue;
            GlomPair gp;
            gp.ag = ag[i];
            gp.scg = scg[j];
            gp.pocketcen = pcen;

            float p1 = gp.get_potential() * frand(1.0-best_binding_stochastic, 1.0+best_binding_stochastic);

            int r = retval.size();
            for (l=0; l<r; l++)
            {
                float ra = retval[l]->ag->get_center().get_3d_distance(ag[i]->get_center());
                float rs = retval[l]->scg->get_center().get_3d_distance(scg[j]->get_center());
                float dr = fabs(ra - rs);
                float deff = fmax(1, dr-scg[j]->glom_reach());
                #if _dbg_glomsel
                cout << "Glom spacing for " << *ag[i] << "-" << *scg[j] << " has delta " << dr << " effectively " << deff << endl;
                #endif
                p1 /= deff;
            }

            #if _dbg_glomsel
            cout << "Potential for " << *ag[i] << "-" << *scg[j] << " = " << p1 << endl << endl;
            #endif

            if (p1 > p)
            {
                j1 = j;
                p = gp.potential;
            }
        }

        if (j1 < 0) continue;

        std::shared_ptr<GlomPair> gp(new GlomPair());
        gp->ag = ag[i];
        gp->scg = scg[j1];

        adirty[i] = true;
        sdirty[j1] = true;

        #if _dbg_glomsel
        cout << "Strongest match for " << *ag[i] << " is " << *scg[j1] << endl;
        #endif
        
        bool added = false;
        int r = retval.size();
        if (!r)
        {
            retval.push_back(gp);
            added = true;
            #if _dbg_glomsel
            cout << "Beginning result with " << *gp->ag << "-" << *gp->scg << endl;
            #endif
        }
        else for (l=0; l<r; l++)
        {
            if (gp->get_potential() > retval[l]->get_potential())
            {
                std::vector<std::shared_ptr<GlomPair>>::iterator it;
                it = retval.begin();
                retval.insert(it+l, gp);
                added = true;
                
                #if _dbg_glomsel
                cout << "Inserting " << *gp->ag << "-" << *gp->scg << " before " << *retval[l+1]->ag << "-" << *retval[l+1]->scg << endl;
                #endif

                break;
            }
        }
        if (!added)
        {
            retval.push_back(gp);
            added = true;
            #if _dbg_glomsel
            cout << "Appending to result " << *gp->ag << "-" << *gp->scg << endl;
            #endif
        }
    }

    #if _dbg_glomsel
    cout << endl << endl << "Final pair assignments:" << endl;
    for (i=0; i<retval.size(); i++)
    {
        cout << *retval[i]->ag << "-" << *retval[i]->scg << endl;
    }
    #endif

    return retval;
}


void GlomPair::align_gloms(Molecule* lig, std::vector<std::shared_ptr<GlomPair>> gp)
{
    int n = gp.size();

    if (n < 1) return;
    Rotation rot = align_points_3d(gp[0]->ag->get_center(), gp[0]->scg->get_center(), lig->get_barycenter());
    LocatedVector lv = rot.v;
    lv.origin = lig->get_barycenter();
    #if _dbg_glomsel
    cout << "Rotating " << *gp[0]->ag << " in the direction of " << *gp[0]->scg << endl;
    #endif
    lig->rotate(lv, rot.a);


    // Scooch.
    float r = gp[0]->ag->get_center().get_3d_distance(gp[0]->scg->get_center()) - gp[0]->scg->glom_reach();
    if (r > 0)
    {
        Point rel = gp[0]->scg->get_center().subtract(gp[0]->ag->get_center());
        rel.scale(r);
        #if _dbg_glomsel
        cout << "Scooching " << *gp[0]->ag << " " << r << "Å into the reach of " << *gp[0]->scg << endl;
        #endif
        lig->move(rel);
    }

    if (n < 2) return;
    rot = align_points_3d(gp[1]->ag->get_center(), gp[1]->scg->get_center(), gp[0]->ag->get_center());
    lv = rot.v;
    lv.origin = gp[0]->ag->get_center();
    #if _dbg_glomsel
    cout << "Rotating " << *gp[1]->ag << " in the direction of " << *gp[1]->scg << endl;
    #endif
    lig->rotate(lv, rot.a);

    if (n < 3) return;
    Point zcen = gp[0]->ag->get_center();
    SCoord axis = gp[1]->ag->get_center().subtract(zcen);
    lv = (SCoord)axis;
    lv.origin = zcen;
    float theta = find_angle_along_vector(gp[2]->ag->get_center(), gp[2]->scg->get_center(), zcen, axis);
    #if _dbg_glomsel
    cout << "\"Rotisserie\" aligning " << *gp[2]->ag << " in the direction of " << *gp[2]->scg << endl;
    #endif
    lig->rotate(lv, theta);
}
