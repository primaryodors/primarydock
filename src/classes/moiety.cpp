
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

        for (j=0; j<molsize; j++)
        {
            Atom* b = mol->get_atom(j);
            if (b) b->used = false;
        }

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
    char* charge = nullptr;
    if (buffer[0] == '[')
    {
        charge = strchr(buffer, '+');
        if (!charge) charge = strchr(buffer, '-');
        if (!charge) charge = strchr(buffer, '0');
        buffer++;
        char* c = strchr(buffer, ']');
        if (c) *c = 0;
    }

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

#define max_bracketed_moiety 10
int Moiety::does_atom_match(Atom* a, Atom** out_matches)
{
    if (!molsize) throw 0xffff;

    const char* p = pattern.c_str();
    int i, j, l, n = strlen(p), result = 0, parens = 0, card = 0, iter;
    char c, buffer[max_bracketed_moiety];
    Atom* cursor[20];
    Atom* numbered[10];
    Atom* dud_atoms[32];
    int num_duds = 0;
    int atoms_used = 0;
    int num_numbers = 0;

    for (iter=0; iter<10; iter++)
    {
        // cout << endl << iter << endl;
        atoms_used = parens = card = 0;
        for (i=0; i<10; i++) numbered[i] = nullptr;

        cursor[parens] = a;
        c = p[0];
        if (c == '[')
        {
            const char* c1 = strchr(p, ']');
            if (!c1)
            {
                cout << "Bad moiety string." << endl;
                throw 0xffff;
            }
            i = c1 - p;
            if (i >= max_bracketed_moiety)
            {
                cout << "max_bracketed_moiety too small." << endl;
                throw 0xffff;
            }
            buffer[i+1] = 0;
            for (; i>=0; i--) buffer[i] = p[i];
            #if _dbg_moieties
            cout << "Starting on bracketed expression " << buffer << endl;
            #endif
        }
        else
        {
            buffer[0] = c;
            buffer[1] = 0;
        }
        i = strlen(buffer);
        if (!atom_matches_string(cursor[parens], buffer)) return 0;
        #if _dbg_moieties
        cout << "Matched " << cursor[parens]->name << " to start of expression " << p << endl;
        #endif
        out_matches[atoms_used++] = cursor[parens];
        cursor[parens]->used = true;

        for (; i<n; i++)
        {
            c = p[i];

            if (c <= ' ') continue;
            /*if (c == '-')
            {
                card = 1;
                continue;
            }*/
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
                    num_numbers--;
                    continue;
                }
                else
                {
                    numbered[c-'0'] = cursor[parens];
                    num_numbers++;
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
                for (l=0; l<=j; l++)
                {
                    buffer[l] = p[i+l];
                }
                buffer[l] = 0;
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
                if (b[j]->btom->used)
                {
                    #if _dbg_moieties
                    cout << "Skipping " << b[j]->btom->name << " as 'used'..." << endl;
                    #endif
                    continue;
                }
                bool used = false;
                for (l=0; l<atoms_used; l++) if (out_matches[l] == b[j]->btom)
                {
                    #if _dbg_moieties
                    cout << "Skipping " << b[j]->btom->name << " already matched..." << endl;
                    #endif
                    used = true;
                }
                for (l=0; l<num_duds; l++) if (dud_atoms[l] == b[j]->btom)
                {
                    #if _dbg_moieties
                    cout << "Skipping " << b[j]->btom->name << " as a 'dud'..." << endl;
                    #endif
                    used = true;
                }
                if (used) continue;
                #if _dbg_moieties
                cout << "Trying " << *cursor[parens] << "~" << *(b[j]->btom) << " for " << buffer << "..." << endl;
                #endif
                if (atom_matches_string(b[j]->btom, buffer))
                {
                    if (card)
                    {
                        Bond* between = cursor[parens]->get_bond_between(b[j]->btom);
                        int lcard = between->cardinality;
                        if (between->cardinality > 1 && between->cardinality < 2 && card >= 1 && card <= 2) lcard = card;
                        if (lcard != card)
                        {
                            #if _dbg_moieties
                            cout << "Skipping " << b[j]->btom->name << " for bond cardinality..." << endl;
                            #endif
                            // return 0;
                            continue;
                        }
                        card = 0;
                    }
                    found = true;
                    #if _dbg_moieties
                    cout << "Matched " << b[j]->btom->name << endl;
                    #endif
                    cursor[parens] = b[j]->btom;
                    out_matches[atoms_used++] = b[j]->btom;
                    b[j]->btom->used = true;
                    break;
                }
                #if _dbg_moieties
                cout << b[j]->btom->name << " did not match " << buffer << "." << endl;
                #endif
            }
            if (!found)
            {
                if (cursor[parens] == a) return 0;
                #if _dbg_moieties
                cout << *cursor[parens] << " is a dud." << endl;
                #endif
                dud_atoms[num_duds++] = cursor[parens];
                dud_atoms[num_duds] = nullptr;
                atoms_used = 0;
            }
        }
        next_iter:
        if (atoms_used) break;
    }
    out_matches[atoms_used] = nullptr;
    if (num_numbers) return 0;
    // cout << "Result: " << atoms_used << endl;
    return atoms_used;
}
