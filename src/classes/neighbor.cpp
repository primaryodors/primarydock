
#include "neighbor.h"
#include "atom.h"
#include "molecule.h"

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

int Neighborhood::fetch_atoms_near(Atom** r, const int m, const Point p, const int d)
{
    int i, j, n = blocks.size(), count = 0;
    if (!n) return 0;

    int x = round(p.x / block_size), y = round(p.y / block_size), z = round(p.z / block_size);

    for (i=0; i<n; i++)
    {
        if (abs(blocks[i].cx - x) <= d && abs(blocks[i].cy - y) <= d && abs(blocks[i].cz - z) <= d)
        {
            n = blocks[i].atom_count;
            for (j=0; j<n; j++) r[count++] = blocks[i].atoms[j];
        }
    }

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
            for (j=0; j<n; j++) r[count++] = blocks[i].mols[j];
        }
    }

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
        if (blocks[i].cx = x && blocks[i].cy == y && blocks[i].cz == z) return &blocks[i];
    }

    return nullptr;
}
