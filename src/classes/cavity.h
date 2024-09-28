#include "protein.h"

#ifndef _CAVITY
#define _CAVITY

#define min_partial_radius 0.7
#define min_dist_bounding_box 11
#define cav_360_step fiftyseventh*4
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
    // int resno = 0;
    std::string resnos_as_string(Protein* p);
};

class Cavity
{
    public:
    static int scan_in_protein(Protein* p, Cavity* results, int results_max);
    CPartial* get_nearest_partial(Point pt);
    float partial_intersects_cavity(CPartial p);
    void add_partial(CPartial p);
    void output_ngl_js(FILE* fp);
    int count_partials();
    CPartial* get_partial_by_idx(int idx) { return &partials[idx]; }
    Point get_center();
    CPartial* point_inside_pocket(Point pt);
    float containment_violations(Molecule* m, float stop_if_more_than = Avogadro);
    float find_best_containment(Molecule* m);

    protected:
    void compute_vdW_surface(float d);
    Point nearest_surface_vertex(Point pt);
    CPartial* partials = nullptr;
    int pallocd = 0;
    bool priority = false;
    Point* vdw_surface = nullptr;
    CPartial** vdw_vertex_partial = nullptr;
    int vdw_vertex_count = 0;
};

extern float cav_xmax, cav_xmin, cav_ymax, cav_ymin, cav_zmax, cav_zmin, cav_xyrlim, cav_xzrlim, cav_yzrlim;
extern int cav_resmin, cav_resmax;

#endif

