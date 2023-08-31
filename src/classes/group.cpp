
#include "group.h"

std::vector<MCoord> mtlcoords;
std::vector<std::shared_ptr<GroupPair>> global_pairs;

void ResiduePlaceholder::set(const char* str)
{
    if (strchr(str, '.')) bw = str;
    else resno = atoi(str);
}

void ResiduePlaceholder::resolve_resno(Protein* prot)
{
    int hxno = atoi(bw.c_str());
    const char* dot = strchr(bw.c_str(), '.');
    if (!dot)
    {
        resno = atoi(bw.c_str());
        return;
    }
    int bwpos = atoi(dot+1);
    resno = prot->get_bw50(hxno) + bwpos - 50;
}

Point AtomGroup::get_center()
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

float AtomGroup::get_pi()
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

float AtomGroup::get_polarity()
{
    int atct = atoms.size();
    if (!atct) return 0;
    int i;
    float result = 0;
    for (i=0; i<atct; i++)
    {
        if (atoms[i]->get_family() == CHALCOGEN) result += 1;
        else if (atoms[i]->get_family() == PNICTOGEN) result += 1;
        else result += fabs(atoms[i]->is_polar());
    }
    return result;
}

float AtomGroup::get_ionic()
{
    int atct = atoms.size();
    if (!atct) return 0;
    int i;
    float result = 0;
    for (i=0; i<atct; i++)
    {
        float c = atoms[i]->get_charge();
        if (c) result += c;
    }
    return result;
}

float AtomGroup::get_mcoord()
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

float AtomGroup::get_sum()
{
    float retval = fabs(get_ionic()*60) + get_polarity()*25 + get_pi()*2;
    if (mtlcoords.size()) retval += get_mcoord()*60 - get_polarity()*40;
    return retval;
}

float AtomGroup::get_avg_elecn()
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
            #if _dbg_groupsel
            cout << atoms[i]->name << " has electronegativity " << atoms[i]->get_electronegativity() << endl;
            #endif
        }
    }
    if (j) result /= j;
    return result;
}

float AtomGroup::hydrophilicity()
{
    int atct = atoms.size();
    if (!atct) return 0;
    int i, j=0, num_heavy_atoms;
    float result = 0, divisor = 0;

    num_heavy_atoms = 0;
    for (i=0; i<atct; i++) if (atoms[i]->get_Z() > 1) num_heavy_atoms++;

    for (i=0; i<atct; i++)
    {
        if (atoms[i]->get_Z() == 1) continue;
        float h = fabs(atoms[i]->is_polar());
        if (num_heavy_atoms == 1 && atoms[i]->is_pi())
        {
            int f = atoms[i]->get_family();
            if (f == PNICTOGEN || f == CHALCOGEN) h = 1;
        }
        if (h > hydrophilicity_cutoff)
        {
            result += h;
            divisor += 1;
        }
        else
        {
            divisor += 1.0/3;
        }
    }

    if (divisor) result /= divisor;
    return result;
}

int AtomGroup::contains_element(const char* esym)
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

float AtomGroup::distance_to(Point pt)
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

float AtomGroup::bounds()
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

float AtomGroup::compatibility(AminoAcid* aa)
{
    int atct = atoms.size();
    if (!atct) return 0;

    float result = 0;
    float lgi = get_ionic(), lgh = get_polarity(), lgp = get_pi();

    float aachg = aa->get_charge();
    if (aa->conditionally_basic()) aachg += protonation(aa->sc_pKa());
    if (lgi && aachg && sgn(lgi) != -sgn(aachg)) return 0;

    if (aa->hydrophilicity() > 0.25)
    {
        if ((lgh / atct) < 0.19) return 0;
    }
    else
    {
        if ((lgh / atct) > hydrophilicity_cutoff) return 0;
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
                if (atoms[i]->is_polar() < hydrophilicity_cutoff) lgh++;
            }
        }
        else if (!aa->has_hbond_donors())
        {
            lgh = 0;
            if (atct)
            for (i=0; i<atct; i++)
            {
                if (atoms[i]->is_polar() > hydrophilicity_cutoff) lgh++;
            }
        }
    }

    if (lgi) result += lgi * -aa->bindability_by_type(ionic) * 100;
    if (lgh) result += fabs(lgh) * fabs(aa->bindability_by_type(hbond)) * 30;
    if (lgp) result += aa->bindability_by_type(pi);

    return result;
}



Point ResidueGroup::get_center()
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
        Atom* a = aa[0];
        delete aa;

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

float ResidueGroup::distance_to(Point pt)
{
    return pt.get_3d_distance(get_center());
}

float ResidueGroup::hydrophilicity()
{
    int amsz = aminos.size();
    if (!amsz) return 0;
    float result = 0, divisor = 0;
    int i;
    bool has_acids = false, has_his = false;
    for (i=0; i<amsz; i++)
    {
        float h = fabs(aminos[i]->hydrophilicity());

        if (h >= 0.2)
        {
            result += h;
            divisor += 1;
        }
        else
        {
            divisor += 1.5;
        }
    }

    if (divisor) result /= divisor;
    return result;
}

float ResidueGroup::compatibility(AtomGroup* ag)
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

        if (aminos[i]->priority)
        {
            f *= priority_weight_group;		// Extra weight for residues mentioned in a CEN RES or PATH RES parameter.
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

        #if _dbg_groupsel
        cout << "Metal multiplier: " << lmm << " from kmims " << kmims << endl;
        #endif

        lmc += 1.0 * fabs(ag->get_pi());
        result += lmc;
    }

    return result;
}

float ResidueGroup::group_reach()
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

void ResidueGroup::conform_to(Molecule* mol)
{
    if (metallic) return;

    int i, n;

    n = aminos.size();
    if (!n) return;

    for (i=0; i<n; i++)
    {
        MovabilityType mt = aminos[i]->movability;
        if (!(mt & MOV_FORBIDDEN))
        {
            aminos[i]->movability = MOV_FORCEFLEX;
            Star s;
            int j, n1;

            Molecule* ll[3];
            ll[0] = mol;
            s.paa = aminos[i];
            ll[1] = s.pmol;
            ll[2] = nullptr;

            Molecule::conform_molecules(ll, 50);

            aminos[i]->movability = mt;
        }
    }
}

bool AtomGroup::is_bonded_to(Atom* a)
{
    int n = atoms.size();
    if (!n) return false;

    int i;
    for (i=0; i<n; i++)
    {
        if (atoms[i]->is_bonded_to(a)) return true;
        else if (atoms[i]->shares_bonded_with(a)) return true;
    }

    return false;
}

int AtomGroup::heavy_atom_count()
{
    int i, n = atoms.size(), result = 0;

    for (i=0; i<n; i++)
    {
        if (atoms[i]->get_Z() > 2) result++;
    }

    return result;
}

std::vector<std::shared_ptr<AtomGroup>> AtomGroup::get_potential_ligand_groups(Molecule* mol)
{
    std::vector<std::shared_ptr<AtomGroup>> retval;
    if (!mol) return retval;
    int n = mol->get_atom_count();
    if (!n) return retval;

    int i, j, k, l;
    bool dirty[n+4];
    for (i=0; i<n; i++) dirty[i] = false;

    std::vector<Atom*> bd = mol->longest_dimension();
    if (bd.size() < 2) throw 0xbad302;
    float ld = bd[0]->distance_to(bd[1]);

    for (i=0; i<n; i++)
    {
        if (dirty[i]) continue;
        std::shared_ptr<AtomGroup> g(new AtomGroup());
        g->ligand = mol;
        int aliphatic = 0;

        Atom* a = mol->get_atom(i);
        if (!a) continue;
        if (a->get_Z() == 1) continue;
        if (!a->is_polar() && !a->get_charge() && !a->is_pi())
        {
            int bh = a->get_bonded_heavy_atoms_count();
            if (bh > 1) continue;
            aliphatic++;
        }

        Atom* a_ = a;
        g->atoms.push_back(a);
        #if _dbg_groupsel
        cout << "Creating group from " << a->name << "..." << endl;
        #endif

        dirty[i] = true;
        for (j=0; j<n; j++)
        {
            if (j==i) continue;
            if (dirty[j]) continue;

            Atom* b = mol->get_atom(j);
            if (!b) continue;

            // if (g->get_ionic() && !b->get_charge() && !b->is_polar() && !b->is_pi()) continue;

            float r = fmax(0, g->get_center().get_3d_distance(b->get_location()) - 1.5);
            if (r > 2.5 && !a->shares_bonded_with(b))
            {
                continue;
                #if _dbg_groupsel
                //cout << "Rejected " << b->name << " too far away." << endl;
                #endif
            }

            r = a->distance_to(b);
            if (r > ld/1.5) continue;
            if (b->get_Z() == 6 && r > ld/2) continue;

            float simil = a->similarity_to(b);

            // In general, aliphatics and nonpolar aromatics can be regarded as high similarity, since they tend to
            // be mutually highly soluble in one another. But for grouping of ligand atoms, it is important to
            if (simil >= 0.5 || (b->get_Z() == 1 && g->is_bonded_to(b)))
            {
                if (aliphatic < 3 || b->get_Z() == 1)
                {
                    if (g->is_bonded_to(b))
                    {
                        g->atoms.push_back(b);
                        if (b->get_Z() > 1)
                        {
                            if (!b->is_polar() && !b->get_charge() && !b->is_pi()) aliphatic++;
                            else aliphatic--;
                        }

                        #if _dbg_groupsel
                        cout << "Adding " << b->name << " with distance " << r << " and similarity " << simil
                             << " polarity " << b->is_polar() << " aliphatic " << aliphatic << endl;
                        #endif

                        a_ = b;
                        if ((bool)(fabs(a->is_polar()) >= 0.2) == (bool)(fabs(b->is_polar()) >= 0.2)) dirty[j] = true;
                    }
                }
            }
            else
            {
                ;
                #if _dbg_groupsel
                if (b->get_family() == a->get_family()) cout << "Rejected " << b->name << " distance " << r << " similarity " << simil << endl;
                #endif
            }
        }

        bool added = false;
        for (l=0; l < retval.size(); l++)
        {
            if (g->get_sum() > retval[l]->get_sum())
            {
                std::vector<std::shared_ptr<AtomGroup>>::iterator it;
                it = retval.begin();
                retval.insert(it+l, g);
                added = true;
                break;
            }
        }
        if (!added) retval.push_back(g);
        #if _dbg_groupsel
        cout << "Group complete." << endl << endl;
        #endif
    }

    l = retval.size();
    for (i=0; i<l; i++)
    {
        int ni = retval[i]->atoms.size();
        for (j=l-1; j>i; j--)
        {
            if (retval[i]->get_center().get_3d_distance(retval[j]->get_center()) > ld/3) continue;

            int nj = retval[j]->atoms.size();
            int si = retval[i]->intersecting(retval[j].get());

            if (retval[i]->average_similarity(retval[j].get()) >= 0.5 && (si >= nj/2 || si >= ni/2))
            {
                #if _dbg_groupsel
                cout << "Merging groups " << *retval[i] << " and " << *retval[j] << endl << endl;
                #endif
                retval[i]->merge(retval[j].get());
                std::vector<std::shared_ptr<AtomGroup>>::iterator it;
                it = retval.begin();
                retval.erase(it+j);
                l--;
            }
        }
    }

    for (i=0; i<l; i++)
    {
        int ni = retval[i]->atoms.size();
        for (j=l-1; j>i; j--)
        {
            int nj = retval[j]->atoms.size();
            if (retval[i]->intersecting(retval[j].get()))
            {
                for (k=0; k<ni; k++)
                {
                    if (retval[j]->contains_atom(retval[i]->atoms[k]))
                    {
                        if (ni < nj)
                        {
                            #if _dbg_groupsel
                            cout << "Removing atom " << *retval[i]->atoms[k] << " from " << *retval[j] << endl << endl;
                            #endif
                            retval[j]->remove_atom(retval[i]->atoms[k]);
                        }
                        else
                        {
                            #if _dbg_groupsel
                            cout << "Removing atom " << *retval[i]->atoms[k] << " from " << *retval[i] << endl << endl;
                            #endif
                            retval[i]->remove_atom(retval[i]->atoms[k]);
                        }
                    }
                }
            }
        }
    }

    for (i=0; i<n; i++)
    {
        Atom* a = mol->get_atom(i);
        if (a->is_pi())
        {
            int fam = a->get_family();
            if (fam == CHALCOGEN || fam == PNICTOGEN)
            {
                std::shared_ptr<AtomGroup> g(new AtomGroup());
                g->ligand = mol;
                g->atoms.push_back(a);
                
                Bond** b = a->get_bonds();
                for (j=0; b[j] && b[j]->btom; j++)
                {
                    if (b[j]->btom->get_Z() == 1) g->atoms.push_back(b[j]->btom);
                }
                
                retval.push_back(g);

                #if _dbg_groupsel
                cout << "Creating group " << *g << endl << endl;
                #endif
            }
        }
    }

    return retval;
}

bool AtomGroup::contains_atom(Atom* a)
{
    int i, n = atoms.size();
    for (i=0; i<n; i++) if (atoms[i] == a) return true;
    return false;
}

void AtomGroup::remove_atom(Atom* a)
{
    int i, n = atoms.size();
    for (i=0; i<n; i++) if (atoms[i] == a)
    {
        atoms.erase(atoms.begin()+i);
        return;
    }
}

float AtomGroup::average_similarity(AtomGroup* cw)
{
    float totsim = 0;
    int simqty = 0;

    int i, j, m = atoms.size(), n = cw->atoms.size();

    for (i=0; i<n; i++)
    {
        Atom* a = cw->atoms[i];

        for (j=0; j<m; j++)
        {
            if (atoms[j] != a)
            {
                totsim += a->similarity_to(atoms[j]);
                simqty++;
            }
        }
    }

    if (!simqty) return 0;

    return totsim / simqty;
}

int AtomGroup::intersecting(AtomGroup* cw)
{
    int samesies = 0;

    int i, j, m = atoms.size(), n = cw->atoms.size();

    for (i=0; i<n; i++)
    {
        Atom* a = cw->atoms[i];

        for (j=0; j<m; j++)
        {
            if (atoms[j] == a) samesies++;
        }
    }

    return samesies;
}

void AtomGroup::merge(AtomGroup* mw)
{
    if (ligand && mw->ligand && ligand != mw->ligand) throw 0xbad3126;
    if (!ligand) ligand = mw->ligand;

    int i, j, m = atoms.size(), n = mw->atoms.size();

    for (i=0; i<n; i++)
    {
        Atom* a = mw->atoms[i];

        for (j=0; j<m; j++)
        {
            if (atoms[j] == a) goto _already;
        }

        atoms.push_back(a);

        _already:
        ;
    }
}

std::vector<std::shared_ptr<ResidueGroup>> ResidueGroup::get_potential_side_chain_groups(AminoAcid** aalist, Point pcen)
{
    std::vector<std::shared_ptr<ResidueGroup>> retval;
    if (!aalist) return retval;
    int i, j, n;
    for (n=0; aalist[n]; n++);          // Get count.
    bool dirty[n+4];
    for (i=0; i<n; i++) dirty[i] = false;

    for (i=0; i<n; i++)
    {
        if (dirty[i]) continue;

        std::shared_ptr<ResidueGroup> g(new ResidueGroup());
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
        #if _dbg_groupsel
        cout << "Building side chain group from " << aa->get_name() << endl;
        #endif

        for (j=i+1; j<n; j++)
        {
            if (dirty[j]) continue;
            AminoAcid* bb = aalist[j];

            CB = bb->get_atom("CB");
            if (CB)
            {
                float a3d = find_3d_angle(CB->get_location(), pcen, bb->get_CA_location());
                if (a3d > fiftyseventh*120)
                {
                    dirty[j] = true;
                    continue;
                }
                else if (bb->get_num_rings() && bb->ring_is_aromatic(0) && a3d > hexagonal)
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
                #if _dbg_groupsel
                cout << "Adding " << bb->get_name() << " because metal coord." << endl;
                #endif
                continue;
            }
            else if (aa->coordmtl || bb->coordmtl)
            {
                #if _dbg_groupsel
                cout << "Rejected " << bb->get_name() << " metal coord mismatch." << endl;
                #endif
                continue;
            }

            float r = fmax(0, aa->get_CA_location().get_3d_distance(bb->get_CA_location()) - fmax(aa->get_reach(), bb->get_reach()));
            if (r > 2.5)
            {
                #if _dbg_groupsel
                // cout << "Rejected " << bb->get_name() << " distance " << r << endl;
                #endif
                continue;
            }

            float simil = aa->similarity_to(bb);
            if (simil >= hydrophilicity_cutoff)
            {
                g->aminos.push_back(bb);
                dirty[j] = true;
                #if _dbg_groupsel
                cout << "Adding " << bb->get_name() << " distance " << r << " similarity " << simil << endl;
                #endif
            }
            else
            {
                ;
                #if _dbg_groupsel
                cout << "Rejected " << bb->get_name() << " similarity " << simil << endl;
                #endif
            }
        }

        #if _dbg_groupsel
        cout << "Completed group." << endl << endl;
        #endif

        retval.push_back(g);
    }

    return retval;
}

std::ostream& operator<<(std::ostream& os, const AtomGroup& ag)
{
    if (!&ag) return os;
    try
    {
        os << "atom_group[ ";
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

std::ostream& operator<<(std::ostream& os, const ResidueGroup& scg)
{
    if (!&scg) return os;
    try
    {
        os << "residue_group[ ";
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

float GroupPair::get_potential()
{
    if (potential) return potential;
    else
    {
        int m = ag->atoms.size(), n = scg->aminos.size();
        if (!m || !n) return 0;

        bool polar_atoms = (fabs(ag->hydrophilicity()) >= 0.25);
        bool polar_res   = (fabs(scg->hydrophilicity()) >= 0.2);

        if (polar_atoms != polar_res) return 0;

        int i, j, q=0;
        for (i=0; i<m; i++)
        {
            Atom* a = ag->atoms[i];
            if (a->get_Z() == 6 && !a->is_polar()) continue;
            for (j=0; j<n; j++)
            {
                AminoAcid* aa = scg->aminos[j];
                float partial;

                if (aa->coordmtl)
                {
                    partial = InteratomicForce::potential_binding(a, aa->coordmtl);
                }
                else
                {
                    partial = aa->get_atom_mol_bind_potential(a);

                    if ((aa->get_charge() > 1 || aa->conditionally_basic()) && a->is_aldehyde())
                    {
                        partial += protonation(aa->sc_pKa())*60;

                        #if _dbg_groupsel
                        cout << "Aldehyde-base potential for " << *a << "..." << *aa << " = " << partial << endl;
                        #endif
                    }
                    else if (fabs(a->is_polar()) > hydrophilicity_cutoff && aa->is_tyrosine_like())
                    {
                        partial /= 3;
                        #if _dbg_groupsel
                        cout << *a << " is polar and " << *aa << " is tyrosine-like." << endl;
                        #endif
                    }
                    else if (fabs(a->is_polar()) > hydrophilicity_cutoff && fabs(aa->hydrophilicity()) < hydrophilicity_cutoff)
                    {
                        partial /= 3;
                        #if _dbg_groupsel
                        cout << *a << " is polar and " << *aa << " is not." << endl;
                        #endif
                    }

                    #if _dbg_groupsel
                    cout << "Potential for " << *a << "..." << *aa << " = " << partial << endl;
                    #endif
                }

                if (aa->priority)
                {
                    partial *= priority_weight_group;
                    priority = true;

                    #if _dbg_groupsel
                    cout << *aa << " has priority." << endl;
                    #endif
                }

                potential += partial;
                q++;
            }
        }

        if (q) potential /= q;

        float r = pocketcen.get_3d_distance(scg->get_center());
        float r1 = ag->get_center().get_3d_distance(ag->get_ligand()->get_barycenter());
        r = fmax(1.5, r-r1);
        potential /= fmax(1, r-1.5);

        return potential;
    }
}

std::vector<std::shared_ptr<GroupPair>> GroupPair::pair_groups(std::vector<std::shared_ptr<AtomGroup>> ag, std::vector<std::shared_ptr<ResidueGroup>> scg, Point pcen, float rel_stoch)
{
    std::vector<std::shared_ptr<GroupPair>> retval;

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
            GroupPair pair;
            pair.ag = ag[i];
            pair.scg = scg[j];
            pair.pocketcen = pcen;

            // Debug trap.
            if (ag[i]->atoms[0]->get_family() == PNICTOGEN && scg[j]->aminos[0]->get_residue_no() == 6)
            {
                l = n;
            }

            float p1 = pair.get_potential() * frand(1.0-best_binding_stochastic*rel_stoch, 1.0+best_binding_stochastic*rel_stoch);

            int r = retval.size();
            for (l=0; l<r; l++)
            {
                float ra = retval[l]->ag->get_center().get_3d_distance(ag[i]->get_center());
                float rs = retval[l]->scg->get_center().get_3d_distance(scg[j]->get_center());
                float dr = fabs(ra - rs);
                float deff = fmax(1, dr-scg[j]->group_reach());
                #if _dbg_groupsel
                cout << "Group spacing for " << *ag[i] << "-" << *scg[j] << " has delta " << dr << " effectively " << deff << endl;
                #endif
                p1 /= deff;
            }

            #if _dbg_groupsel
            cout << "Potential for " << *ag[i] << "-" << *scg[j] << " = " << p1 << endl << endl;
            #endif

            if (p1 > p)
            {
                j1 = j;
                p = pair.potential;
            }
        }

        if (j1 < 0) continue;

        std::shared_ptr<GroupPair> pair(new GroupPair());
        pair->ag = ag[i];
        pair->scg = scg[j1];

        adirty[i] = true;
        sdirty[j1] = true;

        #if _dbg_groupsel
        cout << "Strongest match for " << *ag[i] << " is " << *scg[j1] << " with potential " << p << endl;
        #endif
        
        bool added = false;
        int r = retval.size();
        if (!r)
        {
            retval.push_back(pair);
            added = true;
            #if _dbg_groupsel
            cout << "Beginning result with " << *pair->ag << "-" << *pair->scg << endl;
            #endif
        }
        else for (l=0; l<r; l++)
        {
            if  (
                    (
                        pair->priority == retval[l]->priority
                        &&
                        pair->get_potential() > retval[l]->get_potential()
                    )
                    ||
                    (pair->priority && !retval[l]->priority)
                )
            {
                std::vector<std::shared_ptr<GroupPair>>::iterator it;
                it = retval.begin();
                retval.insert(it+l, pair);
                added = true;
                
                #if _dbg_groupsel
                cout << "Inserting " << *pair->ag << "-" << *pair->scg << " before " << *retval[l+1]->ag << "-" << *retval[l+1]->scg << endl;
                #endif

                break;
            }
        }
        if (!added)
        {
            retval.push_back(pair);
            added = true;
            #if _dbg_groupsel
            cout << "Appending to result " << *pair->ag << "-" << *pair->scg << endl;
            #endif
        }
    }

    #if _dbg_groupsel
    cout << endl << endl << "Final pair assignments:" << endl;
    for (i=0; i<retval.size(); i++)
    {
        cout << *retval[i]->ag << "-" << *retval[i]->scg;
        if (retval[i]->priority) cout << " *";
        cout << endl;
    }
    #endif

    return retval;
}

void GroupPair::align_groups(Molecule* lig, std::vector<std::shared_ptr<GroupPair>> gp)
{
    GroupPair::align_groups(lig, gp, true);
}

void GroupPair::align_groups_noconform(Molecule* lig, std::vector<std::shared_ptr<GroupPair>> gp)
{
    GroupPair::align_groups(lig, gp, false);
}

void GroupPair::align_groups(Molecule* lig, std::vector<std::shared_ptr<GroupPair>> gp, bool do_conforms)
{
    int n = gp.size();

    if (n < 1) return;
    Rotation rot = align_points_3d(gp[0]->ag->get_center(), gp[0]->scg->get_center(), lig->get_barycenter());
    LocatedVector lv = rot.v;
    lv.origin = lig->get_barycenter();
    #if _dbg_groupsel
    cout << "Rotating " << *gp[0]->ag << " in the direction of " << *gp[0]->scg << endl;
    #endif
    lig->rotate(lv, rot.a);

    // Scooch.
    float r = gp[0]->ag->get_center().get_3d_distance(gp[0]->scg->get_center()) - gp[0]->scg->group_reach();
    if (r > 0)
    {
        Point rel = gp[0]->scg->get_center().subtract(gp[0]->ag->get_center());
        rel.scale(r-1);
        #if _dbg_groupsel
        cout << "Scooching " << *gp[0]->ag << " " << r << "Å into the reach of " << *gp[0]->scg << endl;
        #endif
        lig->move(rel);
    }

    if (do_conforms) gp[0]->scg->conform_to(lig);

    if (n < 2) return;
    rot = align_points_3d(gp[1]->ag->get_center(), gp[1]->scg->get_center(), gp[0]->ag->get_center());
    lv = rot.v;
    lv.origin = gp[0]->ag->get_center();
    #if _dbg_groupsel
    cout << "Rotating " << *gp[1]->ag << " in the direction of " << *gp[1]->scg << endl;
    #endif
    lig->rotate(lv, rot.a);

    if (do_conforms) gp[1]->scg->conform_to(lig);

    if (n < 3) return;
    Point zcen = gp[0]->ag->get_center();
    SCoord axis = gp[1]->ag->get_center().subtract(zcen);
    lv = (SCoord)axis;
    lv.origin = zcen;
    float theta = find_angle_along_vector(gp[2]->ag->get_center(), gp[2]->scg->get_center(), zcen, axis);
    #if _dbg_groupsel
    cout << "\"Rotisserie\" aligning " << *gp[2]->ag << " in the direction of " << *gp[2]->scg << endl;
    #endif
    lig->rotate(lv, theta);

    if (do_conforms) gp[2]->scg->conform_to(lig);
}

