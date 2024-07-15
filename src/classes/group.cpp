
#include "group.h"

std::vector<MCoord> mtlcoords;
std::vector<std::shared_ptr<GroupPair>> global_pairs;
std::vector<Moiety> predef_grp;
std::vector<AminoAcid*> ResidueGroup::disqualified_residues;

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
        float m = pt.weight;
        pt.multiply(m);
        result = result.add(pt);
        mass += m;
    }
    if (mass) result.multiply(1.0 / mass);
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

bool AtomGroup::has_hbond_acceptors()
{
    int atct = atoms.size();
    if (!atct) return 0;
    int i;
    bool result = false;

    for (i=0; i<atct; i++) if (atoms[i]->is_polar() <= -hydrophilicity_cutoff && atoms[i]->get_bonded_atoms_count() < 4) result = true;

    return result;
}

bool AtomGroup::has_hbond_donors()
{
    int atct = atoms.size();
    if (!atct) return 0;
    int i;
    bool result = false;

    for (i=0; i<atct; i++) if (atoms[i]->is_polar() >= -hydrophilicity_cutoff) result = true;

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
        if (!aa) continue;
        Atom* a = aa[0];

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

Atom* ResidueGroup::get_nearest_atom(Point pt)
{
    int i, n = aminos.size();

    if (n<1 || n>29) return nullptr;

    Atom* result = nullptr;
    float r = Avogadro;

    for (i=0; i<n; i++)
    {
        Atom* a = aminos[i]->get_nearest_atom(pt);
        float r1 = a->get_location().get_3d_distance(pt);

        // if (r1 > 10000) throw 0xbad;

        if (!i || r1 < r)
        {
            result = a;
            r = r1;
        }
    }

    return result;
}

float ResidueGroup::distance_to(Point pt)
{
    if (!this) return 0;
    Atom* a = get_nearest_atom(pt);
    if (!a) return 0;
    return a->get_location().get_3d_distance(pt);
    // return pt.get_3d_distance(get_center());
}

float ResidueGroup::pi_stackability()
{
    int amsz = aminos.size();
    if (!amsz) return 0;
    float result = 0;
    int i;
    for (i=0; i<amsz; i++)
    {
        result += aminos[i]->pi_stackability(false);
    }

    result /= amsz;
    return result;
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

float ResidueGroup::group_reach(bool ap)
{
    int n = aminos.size();
    if (!n) return 0;
    if (metallic) return 0;
    int i;
    float retval = 0;
    for (i=0; i<n; i++)
    {
        AminoAcid* aa = aminos[i];
        float r = aa->get_reach();
        if (aa->is_tyrosine_like() && !ap) r /= 2;
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

            Molecule::conform_molecules(ll, 20);

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
        else if (atoms[i]->shares_bonded_with(a) && a->get_Z() > 2
            && (atoms[i]->get_family() != TETREL || a->get_family() != TETREL)) return true;
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

std::vector<std::shared_ptr<AtomGroup>> AtomGroup::get_potential_ligand_groups(Molecule* mol, bool sep_mcoord)
{
    std::vector<std::shared_ptr<AtomGroup>> retval;
    if (!mol) return retval;
    int n = mol->get_atom_count();
    if (!n) return retval;

    mol->identify_rings();
    mol->identify_cages();

    if (!predef_grp.size())
    {
        FILE* fp = fopen("data/moieties.dat", "rb");
        if (!fp) throw 0xffff;
        char buffer[1024];
        while (!feof(fp))
        {
            fgets(buffer, 1022, fp);
            if (buffer[0])
            {
                char** words = chop_spaced_words(buffer);

                Moiety m;
                m.pattern = words[0];
                if (words[1]) m.pKa = atof(words[1]);
                predef_grp.push_back(m);
            }
        }
        fclose(fp);
    }

    int i, j, k, l, m;
    bool dirty[n+4];
    for (i=0; i<n; i++) dirty[i] = false;

    int n1 = predef_grp.size();
    for (i=0; i<n1; i++)
    {
        Atom* matches[n*32];
        int times = predef_grp[i].contained_by(mol, matches);
        if (times)
        {
            int n2;
            for (n2=0; matches[n2]; n2++);      // count

            int per_grp = n2 / times;
            for (l=0; l<times; l++)
            {
                bool any_dirty = false;
                for (j=0; j<per_grp; j++)
                {
                    k = mol->atom_idx_from_ptr(matches[l*per_grp+j]);
                    if (dirty[k])
                    {
                        any_dirty = true;
                        break;
                    }

                    Bond* bonds[16];
                    matches[l*per_grp+j]->fetch_bonds(bonds);
                    for (m=0; bonds[m]; m++)
                    {
                        if (!bonds[m]->atom2) continue;
                        if (bonds[m]->atom2->get_Z() > 1) continue;
                        k = mol->atom_idx_from_ptr(bonds[m]->atom2);
                        if (dirty[k])
                        {
                            any_dirty = true;
                            break;
                        }
                    }
                    if (any_dirty) break;
                }
                if (any_dirty) continue;

                #if _dbg_groupsel
                cout << predef_grp[i].pattern << "Matched:";
                for (j=0; matches[j]; j++) cout << " " << *matches[j];
                cout << endl;
                #endif

                std::shared_ptr<AtomGroup> g(new AtomGroup());
                g->ligand = mol;
                for (j=0; j<per_grp; j++)
                {
                    Atom* a = matches[l*per_grp+j];
                    k = mol->atom_idx_from_ptr(a);
                    if (dirty[k]) continue;
                    g->atoms.push_back(a);
                    int fam = a->get_family();

                    Bond* bonds[16];
                    a->fetch_bonds(bonds);
                    for (m=0; bonds[m]; m++)
                    {
                        if (!bonds[m]->atom2) continue;
                        if (bonds[m]->atom2->get_Z() > 1) continue;
                        k = mol->atom_idx_from_ptr(bonds[m]->atom2);
                        if (dirty[k]) continue;
                        g->atoms.push_back(bonds[m]->atom2);
                        dirty[k] = true;
                    }
                }

                g->remove_duplicates();

                j = g->atoms.size();
                for (m=0; m<j; m++)
                {
                    Atom* a = g->atoms[m];
                    if (a->get_Z() == 1) a = a->get_bond_by_idx(0)->get_atom1();
                    if (a->get_family() != TETREL) continue;
                    a->pK = predef_grp[i].pKa;
                    k = mol->atom_idx_from_ptr(a);
                    if (k >= 0) dirty[k] = true;
                }

                // If the group forms a substantial part of the molecule, ignore it.
                if (g->heavy_atom_count() > 0.4 * mol->get_heavy_atom_count()) continue;

                retval.push_back(g);

                if (g->atoms.size() > 2 && g->atoms[0]->num_rings() && g->get_pi() > 0.5*g->atoms.size()
                    && (g->has_hbond_acceptors() || g->has_hbond_donors()))
                {
                    std::vector<std::shared_ptr<AtomGroup>> subg = make_hbond_subgroups(g);
                    retval.insert(std::end(retval), std::begin(subg), std::end(subg));
                }
            }
        }
    }

    int groupcount = retval.size();
    for (j=0; j<groupcount; j++)
    {
        int groupsize = retval[j]->atoms.size();
        for (m=0; m<groupsize; m++)
        {
            int idx = mol->atom_idx_from_ptr(retval[j]->atoms[m]);
            if (idx >= 0) dirty[idx] = true;
        }
    }

    std::vector<Atom*> bd = mol->longest_dimension();
    float ld;
    if (bd.size() < 2) ld = 1;
    else ld = bd[0]->distance_to(bd[1]);

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

        if (a->get_Z() > 1)
        {
            Bond* bb[16];
            a->fetch_bonds(bb);
            if (bb)
            {
                int h, bbhidx;
                for (h=0; bb[h]; h++)
                {
                    if (bb[h]->atom2 && bb[h]->atom2->get_Z() == 1)
                    {
                        bbhidx = mol->atom_idx_from_ptr(bb[h]->atom2);
                        if (bbhidx >= 0 && !dirty[bbhidx])
                        {
                            g->atoms.push_back(bb[h]->atom2);
                            dirty[bbhidx] = true;

                            #if _dbg_groupsel
                            cout << "Adding " << bb[h]->atom2->name << "." << endl;
                            #endif
                        }
                    }
                }
            }
        }

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
                #if _dbg_groupsel
                cout << "Rejected " << b->name << " too far away " << r << " (criterion 516)." << endl;
                #endif
                continue;
            }

            r = a->distance_to(b);
            if (r > ld/1.5)
            {
                #if _dbg_groupsel
                cout << "Rejected " << b->name << " too far away " << r << " (criterion 528)." << endl;
                #endif
                continue;
            }
            if (b->get_Z() == 6 && r > ld/2)
            {
                #if _dbg_groupsel
                cout << "Rejected " << b->name << " too far away " << r << " (criterion 535)." << endl;
                #endif
                continue;
            }

            float simil = a->similarity_to(b);
            if (sep_mcoord)
            {
                int a_fam = a->get_family(), b_fam = b->get_family();
                if (a_fam == TETREL && (b_fam == CHALCOGEN || b_fam == PNICTOGEN)) simil /= 4;
                else if (b_fam == TETREL && (a_fam == CHALCOGEN || a_fam == PNICTOGEN)) simil /= 4;
            }

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

                        if (b->get_Z() > 1)
                        {
                            Bond* bb[16];
                            b->fetch_bonds(bb);
                            if (bb)
                            {
                                int h, bbhidx;
                                for (h=0; bb[h]; h++)
                                {
                                    if (bb[h]->atom2 && bb[h]->atom2->get_Z() == 1)
                                    {
                                        bbhidx = mol->atom_idx_from_ptr(bb[h]->atom2);
                                        if (bbhidx >= 0 && !dirty[bbhidx])
                                        {
                                            g->atoms.push_back(bb[h]->atom2);
                                            dirty[bbhidx] = true;

                                            #if _dbg_groupsel
                                            cout << "Adding " << bb[h]->atom2->name << "." << endl;
                                            #endif
                                        }
                                    }
                                }
                            }
                        }

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

        // If the group forms a substantial part of the molecule, ignore it.
        if (g->heavy_atom_count() > 0.4 * mol->get_heavy_atom_count()) continue;

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
        if (!added)
        {
            g->remove_duplicates();
            retval.push_back(g);
        }
        #if _dbg_groupsel
        cout << "Group complete." << endl << endl;
        #endif

        if (g->atoms.size() > 2 && g->atoms[0]->num_rings() && g->get_pi() > 0.5*g->atoms.size()
            && (g->has_hbond_acceptors() || g->has_hbond_donors()))
        {
            std::vector<std::shared_ptr<AtomGroup>> subg = make_hbond_subgroups(g);
            retval.insert(std::end(retval), std::begin(subg), std::end(subg));
        }
    }

    l = retval.size();
    for (i=0; i<l; i++)
    {
        int ni = retval[i]->atoms.size();
        if (retval[i]->heavy_atom_count() == 1)
        {
            int fam = retval[i]->atoms[0]->get_family();
            if (fam == CHALCOGEN || fam == PNICTOGEN) continue;
        }
        for (j=l-1; j>i; j--)
        {
            if (retval[i]->get_center().get_3d_distance(retval[j]->get_center()) > ld/3) continue;

            int nj = retval[j]->atoms.size();
            if (retval[j]->heavy_atom_count() == 1)
            {
                int fam = retval[j]->atoms[0]->get_family();
                if (fam == CHALCOGEN || fam == PNICTOGEN) continue;
            }
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
        if (retval[i]->heavy_atom_count() == 1)
        {
            int fam = retval[i]->atoms[0]->get_family();
            if (fam == CHALCOGEN || fam == PNICTOGEN) continue;
        }
        for (j=l-1; j>i; j--)
        {
            int nj = retval[j]->atoms.size();
            if (retval[j]->heavy_atom_count() == 1)
            {
                int fam = retval[j]->atoms[0]->get_family();
                if (fam == CHALCOGEN || fam == PNICTOGEN) continue;
            }
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

    return retval;
}

std::vector<std::shared_ptr<AtomGroup>> AtomGroup::make_hbond_subgroups(std::shared_ptr<AtomGroup> g)
{
    std::vector<std::shared_ptr<AtomGroup>> retval;
    int j, l;

    for (l=0; l<g->atoms.size(); l++)
    {
        if (g->atoms[l]->get_family() == CHALCOGEN || g->atoms[l]->get_family() == PNICTOGEN || g->atoms[l]->get_family() == HALOGEN)
        {
            std::shared_ptr<AtomGroup> g1(new AtomGroup());
            g1->atoms.push_back(g->atoms[l]);
            g1->ligand = g->ligand;

            #if _dbg_groupsel
            cout << "Creating subgroup from " << g->atoms[l]->name << "..." << endl;
            #endif

            Bond* lbb[16];
            g->atoms[l]->fetch_bonds(lbb);
            for (j=0; lbb[j]; j++)
            {
                if (lbb[j]->atom2 && lbb[j]->atom2->get_Z() == 1)
                {
                    g1->atoms.push_back(lbb[j]->atom2);

                    #if _dbg_groupsel
                    cout << "Adding " << lbb[j]->atom2->name << "..." << endl;
                    #endif
                }
            }

            retval.push_back(g1);

            #if _dbg_groupsel
            cout << "Group complete." << endl << endl;
            #endif
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
    int i, j, m, n;
    for (n=0; aalist[n]; n++);          // Get count.
    bool dirty[n+4];
    for (i=0; i<n; i++) dirty[i] = false;

    for (i=0; i<n; i++)
    {
        #if _dbg_groupsel
        if (dirty[i]) cout << aalist[i]->get_name() << " already in use." << endl;
        #endif
        if (dirty[i]) continue;

        std::shared_ptr<ResidueGroup> g(new ResidueGroup());
        AminoAcid* aa = aalist[i];

        m = disqualified_residues.size();
        bool dq = false;
        for (j=0; j<m; j++) if (disqualified_residues[j] == aa) dq = true;
        if (dq) continue;

        Atom* CB = aa->get_atom("CB");
        if (CB)
        {
            float a3d = find_3d_angle(CB->get_location(), pcen, aa->get_CA_location());
            if (a3d > fiftyseventh*120)
            {
                #if _dbg_groupsel
                cout << "Rejecting " << aa->get_name() << " for CB angle " << (a3d*fiftyseven) << endl;
                #endif

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
                    // dirty[j] = true;
                    continue;
                }
                else if (bb->get_num_rings() && bb->ring_is_aromatic(0) && a3d > hexagonal)
                {
                    // dirty[j] = true;
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

            Atom *reach1, *reach2;
            reach1 = aa->get_reach_atom();
            if (!reach1) reach1 = aa->get_atom("CA");
            reach2 = bb->get_reach_atom();
            if (!reach2) reach2 = bb->get_atom("CA");
            if ((fabs(reach1->is_polar()) >= hydrophilicity_cutoff) != (fabs(reach2->is_polar()) >= hydrophilicity_cutoff))
            {
                Atom *nearest1, *nearest2;
                nearest1 = aa->get_nearest_atom(reach2->get_location());
                nearest2 = bb->get_nearest_atom(reach1->get_location());

                if ((fabs(reach1->is_polar()) >= hydrophilicity_cutoff) == (fabs(nearest2->is_polar()) >= hydrophilicity_cutoff))
                    reach2 = nearest2;
                else if ((fabs(nearest1->is_polar()) >= hydrophilicity_cutoff) == (fabs(reach2->is_polar()) >= hydrophilicity_cutoff))
                    reach1 = nearest1;
                else
                {
                    #if _dbg_groupsel
                    cout << "Rejected " << bb->get_name() << " incompatible polarities of reach atoms." << endl;
                    #endif
                    continue;
                }
            }

            float r = fmax(0, reach1->distance_to(reach2));
            if (r > bb_group_distance_cutoff)
            {
                #if _dbg_groupsel
                cout << "Rejected " << bb->get_name() << " distance " << r << endl;
                #endif
                continue;
            }

            float theta = 0;
            Atom* CB = bb->get_atom("CB");
            if (CB)
            {
                theta = bb->CB_angle(pcen);
                if (theta > square)
                {
                    #if _dbg_groupsel
                    cout << "Rejected " << bb->get_name() << " pocket center angle " << (theta*fiftyseven) << endl;
                    #endif
                    continue;
                }
            }

            float simil = aa->similarity_to(bb);
            int simil_n = 1, i2;

            for (i2=1; i2 < g->aminos.size(); i2++)
            {
                simil += g->aminos[i2]->similarity_to(bb);
                simil_n++;
            }

            simil /= simil_n;

            if (aa->get_charge() && sgn(aa->get_charge()) == -sgn(bb->get_charge()))
            {
                #if _dbg_groupsel
                cout << "Rejected " << bb->get_name() << " opposite charges." << endl;
                #endif
                continue;
            }

            if (simil >= group_simil_threshold)
            {
                g->aminos.push_back(bb);
                dirty[j] = true;
                #if _dbg_groupsel
                cout << "Adding " << bb->get_name() << " distance " << r << " similarity " << simil << " pocket angle " << (theta*fiftyseven) << endl;
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

float GroupPair::get_weighted_potential()
{
    float f = get_potential();

    if (scg->metallic) f *= 60;
    else if (fabs(ag->hydrophilicity()) >= hydrophilicity_cutoff && fabs(scg->hydrophilicity()) >= hydrophilicity_cutoff) f *= 25;
    else if ((ag->get_pi()/ag->atoms.size()) >= 0.25 && scg->pi_stackability() > 0.25 ) f *= 7;

    // Pair priority never should have overridden binding strength. It was always a means to ensure that specific ligand-protein contacts
    // be made without undue precedence given to other potential contacts *of the same type*. It was never a feature to allow e.g. a van der
    // Waals contact to override e.g. a hydrogen bond.
    if (priority) f *= 2.5;

    return f;
}

float GroupPair::get_potential()
{
    if (potential) return potential;
    else
    {
        int m = ag->atoms.size(), n = scg->aminos.size();
        if (!m || !n) return 0;

        bool polar_atoms = (fabs(ag->hydrophilicity()) >= hydrophilicity_cutoff);
        bool polar_res   = (fabs(scg->hydrophilicity()) >= hydrophilicity_cutoff);

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

                if (polar_atoms && polar_res)
                {
                    if (!aa->has_hbond_acceptors() && a->is_polar() > 0) continue;
                    if (!aa->has_hbond_donors() && a->is_polar() < 0) continue;
                }

                if (aa->coordmtl)
                {
                    partial = InteratomicForce::potential_binding(a, aa->coordmtl);
                }
                else
                {
                    partial = aa->get_atom_mol_bind_potential(a);

                    #if _dbg_groupsel
                    // cout << "Initial partial for " << *a << "..." << *aa << " = " << partial << endl;
                    #endif

                    if (polar_atoms && polar_res && aa->get_charge()) partial *= 1.0 + fabs(aa->get_charge());

                    Moiety amide;
                    amide.pattern = "ocn";
                    Atom* matches[128];

                    if ((aa->get_charge() > 1) && a->is_aldehyde())
                    {
                        partial += protonation(aa->sc_pKa())*60;

                        #if _dbg_groupsel
                        cout << "Aldehyde-base potential for " << *a << "..." << *aa << " = " << partial << endl;
                        #endif
                    }
                    else if (fabs(a->is_polar()) > hydrophilicity_cutoff && a->get_family() != PNICTOGEN && amide.contained_by(aa, matches))
                    {
                        partial *= 2;
                        #if _dbg_groupsel
                        cout << *a << " is polar and " << *aa << " is amide." << endl;
                        #endif
                    }
                    else if (fabs(a->is_polar()) > hydrophilicity_cutoff && aa->is_tyrosine_like())
                    {
                        partial /= 4;
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

                    if (fabs(aa->get_charge()) > hydrophilicity_cutoff && fabs(a->get_charge()) > hydrophilicity_cutoff
                        && sgn(aa->get_charge()) == -sgn(a->get_charge()))
                    {
                        partial += 60.0 * fabs(aa->get_charge()) * fabs(a->get_charge());
                        #if _dbg_groupsel
                        cout << *a << " and " << *aa << " are charged." << endl;
                        #endif
                    }

                    if (!a->get_charge() || !aa->get_charge())
                    {
                        bool apol = fabs(a->is_polar()) > hydrophilicity_cutoff;
                        bool aapol = fabs(aa->hydrophilicity()) > hydrophilicity_cutoff;
                        if (apol != aapol) partial /= 3;
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

                // potential += partial;
                potential = fmax(potential, partial);
                q++;
            }
        }

        // if (q) potential = potential / q + 0.25 * potential;

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

            float p1 = pair.get_potential() * frand(1.0-best_binding_stochastic*rel_stoch, 1.0+best_binding_stochastic*rel_stoch);

            int r = retval.size();
            for (l=0; l<r; l++)
            {
                float ra = retval[l]->ag->get_center().get_3d_distance(ag[i]->get_center());
                float rs = retval[l]->scg->get_center().get_3d_distance(scg[j]->get_center());
                float dr = fabs(ra - rs);
                float deff = fmax(1, dr-scg[j]->group_reach( fabs(ag[i]->get_polarity()) >= hydrophilicity_cutoff ));
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
        pair->get_potential();          // Determines priority.

        adirty[i] = true;
        sdirty[j1] = true;

        #if _dbg_groupsel
        cout << "Strongest match for " << *ag[i] << " is " << *scg[j1] << " with potential " << p;
        if (pair->is_priority()) cout << " and priority";
        cout << "." << endl;
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
            if  (pair->get_weighted_potential() > retval[l]->get_weighted_potential())
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

    #if _dbg_groupsel || _show_final_group_pairs
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

void GroupPair::disqualify()
{
    if (!scg) return;
    if (!scg->aminos.size()) return;
    int i, n = scg->aminos.size();
    for (i=0; i<n; i++)
    {
        if (scg->aminos[i]) ResidueGroup::disqualified_residues.push_back(scg->aminos[i]);
    }
}

void GroupPair::align_groups(Molecule* lig, std::vector<std::shared_ptr<GroupPair>> gp)
{
    GroupPair::align_groups(lig, gp, true);
}

void GroupPair::align_groups_noconform(Molecule* lig, std::vector<std::shared_ptr<GroupPair>> gp)
{
    GroupPair::align_groups(lig, gp, false, bb_realign_amount);
}

void GroupPair::align_groups(Molecule* lig, std::vector<std::shared_ptr<GroupPair>> gp, bool do_conforms, float amount)
{
    int n = min((int)gp.size(), _bb_max_grp);
    if (n < 1) return;
    float r;
    Point origin;
    Rotation rot;
    LocatedVector lv;

    r = gp[0]->scg->distance_to(gp[0]->ag->get_center());
    if (r >= bb_realign_threshold_distance)
    {
        origin = (n > 1) ? gp[1]->ag->get_center() : lig->get_barycenter();
        rot = align_points_3d(gp[0]->ag->get_center(), gp[0]->scg->get_center(), origin);
        lv = rot.v;
        lv.origin = origin;
        #if _dbg_groupsel || _dbg_groupsalign
        cout << "I. Rotating " << *gp[0]->ag << "(" << gp[0]->ag->get_center() << ") "
            << (rot.a * fiftyseven) << "deg about " << rot.v
            << " in the direction of " << *gp[0]->scg << "(" << gp[0]->scg->get_center() << ")." << endl;
        #endif

        #if _dbg_improvements_only_rule
        excuse_deterioration = true;
        #endif

        if (rot.a >= bb_realign_threshold_angle) lig->rotate(lv, rot.a*amount);
    }

    #if enable_bb_scooch
    // Scooch.
    // float r = gp[0]->ag->get_center().get_3d_distance(gp[0]->scg->get_center()) - gp[0]->scg->group_reach();
    float r0 = gp[0]->scg->distance_to(gp[0]->ag->get_center());
    float r1 = 0;
    try
    {
        if (n > 1 && gp.at(1) && gp[1]->scg && gp[1]->ag
            && abs((long)gp[1]->ag.get() - (long)gp[1]->scg.get()) < memsanity)
            r1 = gp[1]->scg->distance_to(gp[1]->ag->get_center());
    }
    catch (...)
    {
        ;
    }

    bool do0 = false, do1 = false;
    Point rel(0,0,0);
    if (r0 > bb_scooch_threshold_distance)
    {
        Point pt = gp[0]->scg->get_center().subtract(gp[0]->ag->get_center());
        pt.scale((r0-bb_scooch_threshold_distance) *amount);
        rel = rel.add(pt);
        do0 = true;
    }
    if (r1 > bb_scooch_threshold_distance)
    {
        Point pt = gp[1]->scg->get_center().subtract(gp[1]->ag->get_center());
        pt.scale((r1-bb_scooch_threshold_distance) *amount);
        rel = rel.add(pt);
        do1 = true;
    }
    if (rel.magnitude())
    {
        if (rel.magnitude() > speed_limit*100) throw 0xbad;
        #if _dbg_groupsel || _dbg_groupsalign
        // cout << "Atom group is " << gp[0]->scg->distance_to(gp[0]->ag->get_center()) << "A from side chain group. ";
        // cout << "Side chain group reach is " << gp[0]->scg->group_reach() << "A." << endl;
        cout << "Ia. Scooching ";
        if (do0) cout << *gp[0]->ag << " ";
        if (do1) cout << *gp[1]->ag << " ";
        cout << rel.magnitude() << "Å into the reach of ";
        if (do0) cout << *gp[0]->scg << " ";
        if (do1) cout << *gp[1]->scg << " ";
        cout << endl;
        #endif

        #if _dbg_improvements_only_rule
        excuse_deterioration = true;
        #endif

        lig->move(rel);
    }
    #endif

    if (do_conforms) gp[0]->scg->conform_to(lig);

    if (n < 2) return;
    r = gp[1]->scg->distance_to(gp[1]->ag->get_center()); // - gp[1]->scg->group_reach();
    if (r >= bb_realign_threshold_distance)
    {
        origin = gp[0]->ag->get_center();
        rot = align_points_3d(gp[1]->ag->get_center(), gp[1]->scg->get_center(), origin);
        lv = rot.v;
        lv.origin = origin;
        #if _dbg_groupsel || _dbg_groupsalign
        cout << "II. Rotating " << *gp[1]->ag << " (" << gp[1]->ag->get_center() << ") "
            << (rot.a * fiftyseven) << "deg about " << rot.v
            << " in the direction of " << *gp[1]->scg << " (" << gp[1]->scg->get_center() << ")." << endl;
        #endif

        #if _dbg_improvements_only_rule
        excuse_deterioration = true;
        #endif

        lig->rotate(lv, rot.a*amount);
        // cout << r << " ~ " << gp[1]->scg->distance_to(gp[1]->ag->get_center()) << endl << endl;

        if (do_conforms) gp[1]->scg->conform_to(lig);
    }

    if (n < 3) return;
    Point zcen = gp[0]->ag->get_center();
    SCoord axis = gp[1]->ag->get_center().subtract(zcen);
    lv = (SCoord)axis;
    lv.origin = zcen;
    float theta = find_angle_along_vector(gp[2]->ag->get_center(), gp[2]->scg->get_center(), zcen, axis);
    #if _dbg_groupsel || _dbg_groupsalign
    cout << "III. \"Rotisserie\" aligning " << *gp[2]->ag << " in the direction of " << *gp[2]->scg << endl;
    #endif

    #if _dbg_improvements_only_rule
    excuse_deterioration = true;
    #endif

    lig->rotate(lv, theta*amount);

    if (do_conforms) gp[2]->scg->conform_to(lig);
}

void AtomGroup::remove_duplicates()
{
    int i, j, n = atoms.size();
    for (i=n-1; i; i--)
    {
        for (j=0; j<i; j++)
        {
            if (atoms[i] == atoms[j])
            {
                atoms.erase(atoms.begin()+i);
                break;
            }
        }
    }
}
