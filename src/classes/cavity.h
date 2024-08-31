#include "protein.h"

#ifndef _CAVITY
#define _CAVITY

#define min_partial_radius 0.7
#define min_dist_bounding_box 11
#define cav_360_step fiftyseventh*2.5
#define cav_xyz_step 1.6
#define cav_min_partials 4
#define cav_linking_threshold 2.2

struct CPartial
{
    Sphere s;
    bool priority = false;
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
    Point get_center();
    CPartial* point_inside_pocket(Point pt);
    float containment_violations(Molecule* m);
    float find_best_containment(Molecule* m);

    protected:
    CPartial* partials = nullptr;
    int pallocd = 0;
    bool priority = false;
};

#endif
