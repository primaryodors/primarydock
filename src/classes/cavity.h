#include "protein.h"

#ifndef _CAVITY
#define _CAVITY

#define min_partial_radius 1.5
#define min_dist_bounding_box 11

struct CPartial
{
    Sphere s;
    bool chargedp = false;
    bool chargedn = false;
    bool metallic = false;
    bool polar = false;
    bool thio = false;
    bool pi = false;
};

class Cavity
{
    public:
    static int scan_in_protein(Protein* p, Cavity* results, int results_max);
    float partial_intersects_cavity(CPartial p);
    void add_partial(CPartial p);
    void output_ngl_js(FILE* fp);
    int count_partials();

    protected:
    CPartial* partials = nullptr;
    int pallocd = 0;
};

#endif
