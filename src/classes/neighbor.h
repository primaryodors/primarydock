
#include "intera.h"
#include <vector>

#ifndef _NEIGHBOR
#define _NEIGHBOR

class Neighborhood
{
    public:
    void add_atom(Atom* neighborly_atom);
    int fetch_atoms_near(Atom* results, const Point location, int depth = 1);

    protected:
    //
};







#endif
