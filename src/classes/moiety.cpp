
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
        if (a->is_backbone) continue;

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
    char* charge = strchr(buffer, '+');
    if (!charge) charge = strchr(buffer, '-');
    if (!charge) charge = strchr(buffer, '0');

    if (charge)
    {
        if (*charge == '+' && a->get_charge() <= 0) return false;
        if (*charge == '-' && a->get_charge() >= 0) return false;
        if (*charge == '0' && fabs(a->get_charge()) >= 0.2) return false;
        *charge = 0;
    }

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
    int i, j, l, n = strlen(p), result = 0, parens = 0, card = 0, iter;
    char c, buffer[10];
    Atom* cursor[20];
    Atom* numbered[10];
    Atom* dud_atoms[32];
    int num_duds = 0;
    int atoms_used = 0;

    for (iter=0; iter<10; iter++)
    {
        // cout << endl << iter << endl;
        atoms_used = parens = card = 0;
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

            if (c <= ' ') continue;
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
            else if (c >= '0' && c <= '9')
            {
                if (numbered[c-'0'])
                {
                    Bond* b = cursor[parens]->get_bond_between(numbered[c-'0']);
                    if (!b) return 0;
                    if (card)
                    {
                        int lcard = b->cardinality;
                        if (b->cardinality > 1 && b->cardinality < 2 && card >= 1 && card <= 2) lcard = card;
                        if (lcard != card) return 0;
                        card = 0;
                    }

                    numbered[c-'0'] = nullptr;
                    continue;
                }
                else
                {
                    numbered[c-'0'] = cursor[parens];
                    continue;
                }
            }
            else if (c == '(')
            {
                parens++;
                cursor[parens] = cursor[parens-1];
                continue;
            }
            else if (c == ')')
            {
                parens--;
                continue;
            }
            else if (c == '[')
            {
                const char* bracket = strchr(&p[i], ']');
                if (!bracket) return 0;
                j = bracket - p - i;
                for (l=1; l<j; l++)
                {
                    buffer[l-1] = p[i+l];
                }
                buffer[l-1] = 0;
                i += j;
            }
            else
            {
                buffer[0] = c;
                buffer[1] = 0;
            }

            Bond* b[16];
            cursor[parens]->fetch_bonds(b);
            int bn = cursor[parens]->get_geometry();

            bool found = 0;
            for (j=0; j<bn; j++)
            {
                if (!b[j]) continue;
                if (!b[j]->btom) continue;
                bool used = false;
                for (l=0; l<atoms_used; l++) if (out_matches[l] == b[j]->btom) used = true;
                for (l=0; l<num_duds; l++) if (dud_atoms[l] == b[j]->btom) used = true;
                if (used) continue;
                // cout << "Trying " << *cursor[parens] << "~" << *(b[j]->btom) << " for " << buffer << "..." << endl;
                if (atom_matches_string(b[j]->btom, buffer))
                {
                    if (card)
                    {
                        Bond* between = cursor[parens]->get_bond_between(b[j]->btom);
                        int lcard = between->cardinality;
                        if (between->cardinality > 1 && between->cardinality < 2 && card >= 1 && card <= 2) lcard = card;
                        if (lcard != card) return 0;
                        card = 0;
                    }
                    found = true;
                    // cout << "Good!" << endl;
                    cursor[parens] = b[j]->btom;
                    out_matches[atoms_used++] = b[j]->btom;
                    break;
                }
            }
            if (!found)
            {
                if (cursor[parens] == a) return 0;
                // cout << *cursor[parens] << " is a dud." << endl;
                dud_atoms[num_duds++] = cursor[parens];
                dud_atoms[num_duds] = nullptr;
                atoms_used = 0;
            }
        }
        next_iter:
        if (atoms_used) break;
    }
    out_matches[atoms_used] = nullptr;
    // cout << "Result: " << atoms_used << endl;
    return atoms_used;
}
