
#include "neighbor.h"
#include "atom.h"
#include "protein.h"

Neighborhood the_neighborhood;

void Neighborhood::add_atom(Atom* a)
{
    Point p = a->get_location();
    Molecule* mol = a->get_molecule();
    Block* pb = get_block_from_location(p);
    if (pb)
    {
        if (pb->atom_count >= block_max_atoms) throw 0x700a7077;
        pb->atoms[pb->atom_count++] = a;

        if (mol)
        {
            int i;
            for (i=0; i<pb->mol_count; i++) if (pb->mols[i] == mol) return;
            pb->mols[pb->mol_count++] = mol;
        }

        /*cout << "Atom located at " << a->get_location() << " belongs to block " << pb->serial
            << " with xyz = " << pb->cx << "," << pb->cy << "," << pb->cz << "." << endl << flush;*/
    }
    else
    {
        int x = round(p.x / block_size), y = round(p.y / block_size), z = round(p.z / block_size);
        Block b;
        b.cx = x;
        b.cy = y;
        b.cz = z;
        b.atoms[0] = a;
        b.atom_count = 1;
        b.mols[0] = mol;
        b.mol_count = 1;
        b.serial = blocks.size();
        blocks.push_back(b);

        /*cout << "Atom located at " << a->get_location()
            << " belongs to new block " << b.serial << " with xyz = " << b.cx << "," << b.cy << "," << b.cz << "." << endl << flush;*/
    }
}

void Neighborhood::update_atom(Atom* a, Point old)
{
    Block* ob = get_block_from_location(old);
    if (ob)
    {
        int i, j, n = ob->atom_count;
        for (i=0; i<n; i++)
        {
            if (ob->atoms[i] == a)
            {
                for (j=i+1; j<n; j++) ob->atoms[j-1] = ob->atoms[j];
                ob->atom_count--;
                break;
            }
        }
    }

    add_atom(a);
}

void Neighborhood::remove_atom(Atom* a)
{
    Block* ob = get_block_from_location(a->get_location());
    if (ob)
    {
        int i, j, n = ob->atom_count;
        for (i=0; i<n; i++)
        {
            if (ob->atoms[i] == a)
            {
                for (j=i+1; j<n; j++) ob->atoms[j-1] = ob->atoms[j];
                ob->atom_count--;
                break;
            }
        }
    }
}

void Neighborhood::set_active_protein(Protein* p)
{
    int i, n = p->get_end_resno();

    for (i=1; i<=n; i++)
    {
        AminoAcid* aa = p->get_residue(i);
        if (!aa) continue;

        set_active_ligand(aa);
    }
}

void Neighborhood::set_active_ligand(Molecule* m)
{
    int i, n = m->get_atom_count();
    for (i=0; i<n; i++)
    {
        Atom* a = m->get_atom(i);
        if (!a) continue;
        a->active_neighbor = true;
        actives.push_back(a);
    }
}

void Neighborhood::clear_active_neighbors()
{
    int i, n = actives.size();
    for (i=0; i<n; i++)
    {
        actives[i]->active_neighbor = false;
    }
    actives.clear();
}


int Neighborhood::fetch_atoms_near(Atom** r, const int m, const Point p, const int d)
{
    int i, j, l, n = blocks.size(), count = 0;
    if (!n) return 0;

    int x = round(p.x / block_size), y = round(p.y / block_size), z = round(p.z / block_size);

    for (i=0; i<n; i++)
    {
        if (abs(blocks[i].cx - x) <= d && abs(blocks[i].cy - y) <= d && abs(blocks[i].cz - z) <= d)
        {
            l = blocks[i].atom_count;
            for (j=0; j<l; j++)
            {
                if (blocks[i].atoms[j]->active_neighbor)
                {
                    r[count++] = blocks[i].atoms[j];
                    if (count >= m) goto _too_many_atoms;
                }
            }
        }
    }

    _too_many_atoms:
    r[count] = nullptr;
    return count;
}

int Neighborhood::fetch_molecules_near(Molecule** r, const int m, const Point p, const int d)
{
    int i, j, n = blocks.size(), count = 0;
    if (!n) return 0;

    int x = round(p.x / block_size), y = round(p.y / block_size), z = round(p.z / block_size);

    for (i=0; i<n; i++)
    {
        if (abs(blocks[i].cx - x) <= d && abs(blocks[i].cy - y) <= d && abs(blocks[i].cz - z) <= d)
        {
            n = blocks[i].mol_count;
            for (j=0; j<n; j++)
            {
                r[count++] = blocks[i].mols[j];
                if (count >= m) goto _too_many_mols;
            }
        }
    }

    _too_many_mols:
    r[count] = nullptr;
    return count;
}

Block* Neighborhood::get_block_from_location(Point p)
{
    int i, n = blocks.size();
    if (!n) return nullptr;

    int x = round(p.x / block_size), y = round(p.y / block_size), z = round(p.z / block_size);

    for (i=0; i<n; i++)
    {
        if (blocks[i].cx == x && blocks[i].cy == y && blocks[i].cz == z) return &blocks[i];
    }

    return nullptr;
}

double Neighborhood::total_system_energy()
{
    int i, j, l, blockct, atomct, interct;
    double result = 0;

    blockct = blocks.size();
    for (i=0; i<blockct; i++)
    {
        atomct = blocks[i].atom_count;
        for (j=0; j<atomct; j++)
        {
            Atom* a = blocks[i].atoms[j];
            if (!a->active_neighbor) continue;
            Atom* nearby[16384];
            interct = fetch_atoms_near(nearby, 16380, a->get_location(), 1);
            for (l=0; l<interct; l++)
            {
                Atom* b = nearby[l];
                if (a->is_bonded_to(b) || a->shares_bonded_with(b)) continue;
                if (!b->active_neighbor) continue;
                if (a->aaletter == 'R' && b->aaletter == 'R' && !strcmp(a->name, "CZ") && !strcmp(b->name, "HE1")) continue;        // KLUDGE!!!
                if (b > a)
                {
                    double d = -InteratomicForce::total_binding(a, b);
                    /*if (d > 1000)
                    {
                        cout << "CLASH " << a->aa3let << a->residue << ":" << a->name
                            << " ! " << b->aa3let << b->residue << ":" << b->name << ": " << d << endl << flush;
                    }*/
                    result += d;
                }
            }
        }
    }

    return result;
}

void Neighborhood::set_initial_energy()
{
    initial_energy = total_system_energy();
}

double Neighborhood::total_energy_delta()
{
    return total_system_energy() - initial_energy;
}
