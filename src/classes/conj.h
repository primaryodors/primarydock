
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
    float get_net_charge();
    Point get_barycenter();
    Atom* get_nearest_atom(Point pt);
    Atom* get_nearest_atom(Conjugation* conj);

    bool spent = false;

    protected:
    Atom** atoms = nullptr;           // Not using a vector for performance reasons.
    float net_charge = 0;
    bool net_charge_known = false;
    Conjugation* mutual = nullptr;
    Atom* mutual_nearest = nullptr;

    void add_atom(Atom* a, Atom* prev, Atom* orig);
};


#endif
