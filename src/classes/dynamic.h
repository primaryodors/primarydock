
#ifndef _DYNCLS
#define _DYNCLS

#include "protein.h"

#define MAX_DYN_CONSTRAINTS 5
#define DYN_RANDOM_RANGE 0.333
#define DYN_BIAS_EFFECT 0.4
#define DYN_MAX_OVERAGE 0.25

#define debug_dyn_motion 1

enum DynamicType
{
    dyn_rock,
    dyn_bend = dyn_rock,
    dyn_move,
    dyn_wind
};

enum DynamicConstraintType
{
    ddep_MIN,
    ddep_MAX
};

class DynamicMotion;

class DynamicConstraint
{
    public:
    DynamicMotion* depends_on = nullptr;
    DynamicConstraintType type;
};

class DynamicMotion
{
    public:
    DynamicType type;
    std::string name;
    BallesterosWeinstein start_resno, end_resno, fulcrum_resno, axis_resno;
    float bias = 0;

    DynamicMotion(Protein* ppro);
    bool add_constraint(DynamicConstraint* new_cons);           // Returns false if unable to add.
    void read_config_line(const char* line, DynamicMotion** all_motions);
    float apply_incremental(float additional_amount);           // Returns total applied.
    float apply_absolute(float target_amount);                  // Returns actual applied, e.g. within constraints etc.
    void make_random_change();
    void undo();

    protected:
    Protein* prot = nullptr;
    float applied = 0, last_change = 0;
    DynamicConstraint* constraints[MAX_DYN_CONSTRAINTS+1];      // Not using a std::vector since those cause double free errors on destruct.
};


#endif
