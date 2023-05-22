
#include "group.h"

std::vector<int> extra_wt;
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
        result += fabs(atoms[i]->is_polar());
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

std::vector<std::shared_ptr<AtomGroup>> AtomGroup::get_potential_ligand_groups(Molecule* mol)
{
    std::vector<std::shared_ptr<AtomGroup>> retval;
    if (!mol) return retval;
    int n = mol->get_atom_count();
    if (!n) return retval;

    int i, j, l;
    bool dirty[n+4];
    for (i=0; i<n; i++) dirty[i] = false;

    for (i=0; i<n; i++)
    {
        if (dirty[i]) continue;
        std::shared_ptr<AtomGroup> g(new AtomGroup());
        int aliphatic = 0;

        Atom* a = mol->get_atom(i);
        if (!a) continue;
        if (a->get_Z() == 1) continue;
        if (!a->is_polar() && !a->get_charge() && !a->is_pi())
        {
            if (retval.size() >= 2)
            {
                aliphatic += 10000;
                continue;
            }
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

            float r = fmax(0, g->get_center().get_3d_distance(b->get_location()) - 1.5);
            if (r > 2.5)
            {
                continue;
                #if _dbg_groupsel
                //cout << "Rejected " << b->name << " too far away." << endl;
                #endif
            }
            float simil = 0; // fmax(a->similarity_to(b), a_->similarity_to(b));

            if (a->get_charge() && sgn(a->get_charge()) == sgn(b->get_charge())) simil += 10;
            if ((bool)a->is_polar() == (bool)b->is_polar()) simil += 7;
            if (a->is_pi() == b->is_pi()) simil += 3;
            if (a->is_conjugated_to(b)) simil += 6;

            if (simil >= 8)
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
                        cout << "Adding " << b->name << " with distance " << r << " and similarity " << simil << " " << aliphatic << endl;
                        #endif

                        a_ = b;
                        if ((bool)a->is_polar() == (bool)b->is_polar()) dirty[j] = true;
                    }
                }
            }
            else
            {
                ;
                #if _dbg_groupsel
                //cout << "Rejected " << b->name << " similarity " << simil << endl;
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

    throw 0xbadbeef;
    return retval;
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
            if (simil >= 4)
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
                    if (fabs(a->is_polar()) > 0.333 && aa->is_tyrosine_like()) partial /= 3;

                    if (extra_wt.size()
                            &&
                            std::find(extra_wt.begin(), extra_wt.end(), aa->get_residue_no())!=extra_wt.end()
                    )
                    {
                        partial *= 2.5;		// Extra weight for residues mentioned in a CEN RES or PATH RES parameter.
                    }

                    #if _dbg_groupsel
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

std::vector<std::shared_ptr<GroupPair>> GroupPair::pair_groups(std::vector<std::shared_ptr<AtomGroup>> ag, std::vector<std::shared_ptr<ResidueGroup>> scg, Point pcen)
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
            GroupPair global_pairs;
            global_pairs.ag = ag[i];
            global_pairs.scg = scg[j];
            global_pairs.pocketcen = pcen;

            float p1 = global_pairs.get_potential() * frand(1.0-best_binding_stochastic, 1.0+best_binding_stochastic);

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
                p = global_pairs.potential;
            }
        }

        if (j1 < 0) continue;

        std::shared_ptr<GroupPair> global_pairs(new GroupPair());
        global_pairs->ag = ag[i];
        global_pairs->scg = scg[j1];

        adirty[i] = true;
        sdirty[j1] = true;

        #if _dbg_groupsel
        cout << "Strongest match for " << *ag[i] << " is " << *scg[j1] << endl;
        #endif
        
        bool added = false;
        int r = retval.size();
        if (!r)
        {
            retval.push_back(global_pairs);
            added = true;
            #if _dbg_groupsel
            cout << "Beginning result with " << *global_pairs->ag << "-" << *global_pairs->scg << endl;
            #endif
        }
        else for (l=0; l<r; l++)
        {
            if (global_pairs->get_potential() > retval[l]->get_potential())
            {
                std::vector<std::shared_ptr<GroupPair>>::iterator it;
                it = retval.begin();
                retval.insert(it+l, global_pairs);
                added = true;
                
                #if _dbg_groupsel
                cout << "Inserting " << *global_pairs->ag << "-" << *global_pairs->scg << " before " << *retval[l+1]->ag << "-" << *retval[l+1]->scg << endl;
                #endif

                break;
            }
        }
        if (!added)
        {
            retval.push_back(global_pairs);
            added = true;
            #if _dbg_groupsel
            cout << "Appending to result " << *global_pairs->ag << "-" << *global_pairs->scg << endl;
            #endif
        }
    }

    #if _dbg_groupsel
    cout << endl << endl << "Final pair assignments:" << endl;
    for (i=0; i<retval.size(); i++)
    {
        cout << *retval[i]->ag << "-" << *retval[i]->scg << endl;
    }
    #endif

    return retval;
}

void GroupPair::align_groups(Molecule* lig, std::vector<std::shared_ptr<GroupPair>> global_pairs)
{
    int n = global_pairs.size();

    if (n < 1) return;
    Rotation rot = align_points_3d(global_pairs[0]->ag->get_center(), global_pairs[0]->scg->get_center(), lig->get_barycenter());
    LocatedVector lv = rot.v;
    lv.origin = lig->get_barycenter();
    #if _dbg_groupsel
    cout << "Rotating " << *global_pairs[0]->ag << " in the direction of " << *global_pairs[0]->scg << endl;
    #endif
    lig->rotate(lv, rot.a);

    // Scooch.
    float r = global_pairs[0]->ag->get_center().get_3d_distance(global_pairs[0]->scg->get_center()) - global_pairs[0]->scg->group_reach();
    if (r > 0)
    {
        Point rel = global_pairs[0]->scg->get_center().subtract(global_pairs[0]->ag->get_center());
        rel.scale(r-1);
        #if _dbg_groupsel
        cout << "Scooching " << *global_pairs[0]->ag << " " << r << "Å into the reach of " << *global_pairs[0]->scg << endl;
        #endif
        lig->move(rel);
    }

    global_pairs[0]->scg->conform_to(lig);

    if (n < 2) return;
    rot = align_points_3d(global_pairs[1]->ag->get_center(), global_pairs[1]->scg->get_center(), global_pairs[0]->ag->get_center());
    lv = rot.v;
    lv.origin = global_pairs[0]->ag->get_center();
    #if _dbg_groupsel
    cout << "Rotating " << *global_pairs[1]->ag << " in the direction of " << *global_pairs[1]->scg << endl;
    #endif
    lig->rotate(lv, rot.a);

    global_pairs[1]->scg->conform_to(lig);

    if (n < 3) return;
    Point zcen = global_pairs[0]->ag->get_center();
    SCoord axis = global_pairs[1]->ag->get_center().subtract(zcen);
    lv = (SCoord)axis;
    lv.origin = zcen;
    float theta = find_angle_along_vector(global_pairs[2]->ag->get_center(), global_pairs[2]->scg->get_center(), zcen, axis);
    #if _dbg_groupsel
    cout << "\"Rotisserie\" aligning " << *global_pairs[2]->ag << " in the direction of " << *global_pairs[2]->scg << endl;
    #endif
    lig->rotate(lv, theta);

    global_pairs[2]->scg->conform_to(lig);
}

