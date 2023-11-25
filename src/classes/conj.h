
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "atom.h"

#ifndef _CONJ
#define _CONJ

class Conjugation
{
    public:
    Conjugation();
    Conjugation(Atom* from_atom);
    ~Conjugation();

    void add_atom(Atom* add);
    float get_net_charge() { return net_charge; }
    Atom* get_nearest_atom(Point pt);
    Atom* get_nearest_atom(Conjugation* conj);

    protected:
    Atom** atoms = nullptr;           // Not using a vector for performance reasons.
    float net_charge = 0;

    void add_atom(Atom* a, Atom* prev, Atom* orig);
};


#endif
