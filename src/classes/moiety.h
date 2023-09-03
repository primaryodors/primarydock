
#ifndef _MOIETY
#define _MOIETY

#include <string>
#include "molecule.h"

#define _dbg_moieties 1

class Moiety
{
    public:
    int contained_by(Molecule* mol, Atom** out_matches);
    std::string pattern;

    protected:
    int does_atom_match(Atom* a, Atom** out_matches);
    bool atom_matches_string(Atom* a, char* buffer);

    int molsize = 0;
};

#endif


