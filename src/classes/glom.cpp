
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

std::vector<AtomGlom> AtomGlom::get_potential_ligand_gloms(Molecule* mol)
{
    std::vector<AtomGlom> retval;
    if (!mol) return retval;
    int n = mol->get_atom_count();
    if (!n) return retval;

    int i, j, l;
    bool dirty[n+4];
    for (i=0; i<n; i++) dirty[i] = false;

    for (i=0; i<n; i++)
    {
        if (dirty[i]) continue;
        AtomGlom g;
        Atom* a = mol->get_atom(i);
        if (!a) continue;
        Atom* a_ = a;
        g.atoms.push_back(a);
        #if _dbg_glomsel
        cout << "Creating glom from " << a->name << "..." << endl;
        #endif

        dirty[i] = true;
        for (j=i+1; j<n; j++)
        {
            Atom* b = mol->get_atom(j);
            if (!b) continue;

            float r = fmax(0, g.get_center().get_3d_distance(b->get_location()) - 1.5);
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
                g.atoms.push_back(b);
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
            if (g.get_sum() > retval[l].get_sum())
            {
                std::vector<AtomGlom>::iterator it;
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

std::vector<ResidueGlom> ResidueGlom::get_potential_side_chain_gloms(AminoAcid** aalist)
{
    std::vector<ResidueGlom> retval;
    if (!aalist) return retval;
    int i, j, n;
    for (n=0; aalist[n]; n++);          // Get count.
    bool dirty[n+4];
    for (i=0; i<n; i++) dirty[i] = false;

    for (i=0; i<n; i++)
    {
        if (dirty[i]) continue;

        ResidueGlom g;
        AminoAcid* aa = aalist[i];
        g.aminos.push_back(aa);
        dirty[i] = true;
        #if _dbg_glomsel
        cout << "Building side chain glom from " << aa->get_name() << endl;
        #endif

        for (j=i+1; j<n; j++)
        {
            AminoAcid* bb = aalist[j];

            if (aa->coordmtl && aa->coordmtl == bb->coordmtl)
            {
                g.aminos.push_back(bb);
                dirty[j] = true;
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
                g.aminos.push_back(bb);
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
            for (j=0; j<n; j++)
            {
                AminoAcid* aa = scg->aminos[j];
                float partial;
                if (aa->coordmtl) partial += InteratomicForce::potential_binding(a, aa->coordmtl);
                else partial = aa->get_atom_mol_bind_potential(a);
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

std::vector<GlomPair> GlomPair::pair_gloms(std::vector<AtomGlom> ag, std::vector<ResidueGlom> scg, Point pcen)
{
    std::vector<GlomPair> retval;

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
            gp.ag = &ag[i];
            gp.scg = &scg[j];
            gp.pocketcen = pcen;

            float p1 = gp.get_potential() * frand(1.0-best_binding_stochastic, 1.0+best_binding_stochastic);

            if (p1 > p)
            {
                j1 = j;
                p = gp.potential;
            }
        }

        if (j1 < 0) continue;

        GlomPair gp;
        gp.ag = &ag[i];
        gp.scg = &scg[j1];

        #if _dbg_glomsel
        cout << "Strongest match for " << ag[i] << " is " << scg[j1] << endl;
        #endif
        
        bool added = false;
        int r = retval.size();
        if (!r)
        {
            retval.push_back(gp);
            added = true;
            #if _dbg_glomsel
            cout << "Beginning result with " << *gp.ag << "-" << *gp.scg << endl;
            #endif
        }
        else for (l=0; l<r; l++)
        {
            if (gp.get_potential() > retval[l].get_potential())
            {
                std::vector<GlomPair>::iterator it;
                it = retval.begin();
                retval.insert(it+l, gp);
                added = true;
                
                #if _dbg_glomsel
                cout << "Inserting " << *gp.ag << "-" << *gp.scg << " before " << *retval[l+1].ag << "-" << *retval[l+1].scg << endl;
                #endif

                break;
            }
        }
        if (!added)
        {
            retval.push_back(gp);
            added = true;
            #if _dbg_glomsel
            cout << "Appending to result " << *gp.ag << "-" << *gp.scg << endl;
            #endif
        }
    }

    #if _dbg_glomsel
    cout << endl << endl << "Final pair assignments:" << endl;
    for (i=0; i<retval.size(); i++)
    {
        cout << *retval[i].ag << "-" << *retval[i].scg << endl;
    }
    #endif

    return retval;
}



