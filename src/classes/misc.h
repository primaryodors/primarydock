
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>
#include <math.h>

#include "constants.h"

#ifndef _MISC
#define _MISC

struct SCoord;
struct Point;
struct LocatedVector;
struct Rotation;
struct LocRotation;
class Atom;
class Conjugation;
struct Tug;
class Bond;
class InteratomicForce;
class Molecule;
class AminoAcid;
class Protein;
struct Region;
struct Rotation;

union Star
{
    long long n;
    float f;
    char* psz;
    const char* cpsz;
    void* p;
    Point* ppt;
    SCoord* pv;
    Atom* pa;
    Atom** ppa;
    Bond* pb;
    Bond** ppb;
    InteratomicForce* pif;
    Molecule* pmol;
    AminoAcid* paa;
    Protein* pprot;
    Region* preg;
    Rotation* prot;
    Conjugation* pcon;
};

int in_array(void* needle, void** haystack);
int in_array(int needle, int* haystack);
int in_array(Star needle, Star* haystack);
Star* array_unique(Star* input_array);
char** chop_spaced_words(char* line, char separator = ' ');
int greek_from_aname(const char* aname);
int randsgn();
float frand(float min, float max);
float Pearson_correlation(float* xarr, float* yarr, int length);
enum STR_PAD {STR_PAD_RIGHT, STR_PAD_LEFT, STR_PAD_BOTH};
std::string str_pad(const std::string &str, int pad_length, std::string pad_string=" ", STR_PAD pad_type=STR_PAD_RIGHT);
std::string cardinality_printable(float card);
float protonation(float pKa);
float larger(float v1, float v2);
bool file_exists(std::string fname);

#if active_persistence
float residue_binding_multiplier(int resno);
#endif

void colorrgb(int r, int g, int b);
void colorize(float f);
void colorless();

// From here: https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

template <typename T> void init_nulls(T* array, int length)
{
    int i;
    for (i=0; i<length; i++) array[i] = 0;
}

extern float _INTERA_R_CUTOFF;
extern const char* Greek;
extern std::ofstream *debug;
extern bool last_iter;
extern bool differential_dock;
extern float pre_ligand_multimol_radius;
extern float pre_ligand_flex_radius;

extern char asterisk[5];

#if active_persistence_noflex
extern bool allow_ligand_flex;
#endif

#if _dbg_too_slow
extern std::vector<std::string> performance_function_names;
extern std::vector<double> performance_function_times;
extern std::vector<double> performance_function_began;

void performance_begin_clock(std::string fnname);
void performance_end_clock(std::string fnname);
#endif

#endif

