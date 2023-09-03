
#include <iostream>
#include <cstring>
#include "moiety.h"

int Moiety::contained_by(Molecule* mol, Atom** out_matches)
{
    int i, j, l=0, num_matched = 0;
    molsize = mol->get_atom_count();

    for (i=0; i<molsize; i++)
    {
        Atom* a = mol->get_atom(i);
        if (!a) break;

        if (j = does_atom_match(a, out_matches + l))
        {
            #if _dbg_moieties
            cout << "Matched " << j << " atoms." << endl;
            #endif
            l += j;
            num_matched++;
        }
    }
    out_matches[l] = nullptr;

    return num_matched;
}

bool Moiety::atom_matches_string(Atom* a, char* buffer)
{
    int Z = Atom::Z_from_esym(buffer);
    if (Z)
    {
        return (Z == a->get_Z());
    }

    buffer[0] &= 0xdf;
    Z = Atom::Z_from_esym(buffer);

    if (Z)
    {
        if (Z != a->get_Z()) return false;
        if (!a->is_pi()) return false;
        return true;
    }

    return false;
}

int Moiety::does_atom_match(Atom* a, Atom** out_matches)
{
    if (!molsize) throw 0xffff;

    const char* p = pattern.c_str();
    int i, j, l, n = strlen(p), result = 0, parens = 0, card = 0;
    char c, buffer[10];
    Atom* cursor[20];
    Atom* numbered[10];

    int atoms_used = 0;

    for (i=0; i<10; i++) numbered[i] = nullptr;

    cursor[parens] = a;
    c = p[0];
    buffer[0] = c;
    buffer[1] = 0;
    if (!atom_matches_string(cursor[parens], buffer)) return 0;
    out_matches[atoms_used++] = cursor[parens];

    for (i=1; i<n; i++)
    {
        c = p[i];
        buffer[0] = c;
        buffer[1] = 0;

        if (c == '-')
        {
            card = 1;
            continue;
        }
        else if (c == '=')
        {
            card = 2;
            continue;
        }
        else if (c == '#')
        {
            card = 3;
            continue;
        }

        Bond** b = cursor[parens]->get_bonds();
        int bn = cursor[parens]->get_geometry();

        bool found = 0;
        for (j=0; j<bn; j++)
        {
            if (!b[j]) continue;
            if (!b[j]->btom) continue;
            bool used = 0;
            for (l=0; l<atoms_used; l++) if (out_matches[l] == b[j]->btom) used = true;
            if (used) continue;
            if (atom_matches_string(b[j]->btom, buffer))
            {
                if (card)
                {
                    Bond* between = cursor[parens]->get_bond_between(b[j]->btom);
                    int lcard = between->cardinality;
                    if (between->cardinality > 1 && between->cardinality < 2 && card >= 1 && card <= 2) lcard = card;
                    if (lcard != card) return 0;
                }
                found = true;
                cursor[parens] = b[j]->btom;
                out_matches[atoms_used++] = b[j]->btom;
                break;
            }
        }
        if (!found) return 0;
    }

    out_matches[atoms_used] = nullptr;
    return atoms_used;
}
