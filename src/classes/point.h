
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>

#include "constants.h"

#ifndef _POINT
#define _POINT

using namespace std;

class SCoord;

class Point
{
	public:
    float x=0, y=0, z=0;
    float weight = 1;

    Point();
    Point(SCoord v);
    Point(SCoord* v);
    Point(float x, float y, float z);
    Point add(Point add_to)    {        return add(&add_to);    }
    Point add(Point* add_to);
    Point add(SCoord* add_to);
    Point subtract(const Point subtracted);
    Point subtract(const Point* subtracted);
    Point subtract(SCoord* subtracted)
    {
        Point pt(subtracted);
        return subtract(pt);
    }
    Point negate();
    float get_3d_distance(const Point reference)    {        return get_3d_distance(&reference);    }
    float get_3d_distance(const Point* reference);
    Point multiply_3d_distance(const Point* reference, float r_mult);
    bool pt_in_bounding_box(const Point* corner1, const Point* corner2);
    float magnitude();
    void scale(float new_magn);

    std::string printable() const;

    Point& operator=(SCoord v);
    friend std::ostream& operator<<(std::ostream& os, const SCoord& v);
};

class SCoord
{
	public:
    float r=0, theta=0, phi=0;

    SCoord()    {        r = theta = phi = 0;    }
    SCoord(const Point from);
    SCoord(const Point* from);
    SCoord(float r, float theta, float phi);
    SCoord negate()
    {
        Point pt(this);
        pt = pt.negate();
        SCoord v(&pt);
        return v;
    }
    SCoord add(SCoord v)    {        return add(&v);    };
    SCoord add(SCoord* v);

    std::string printable() const;

    SCoord& operator=(Point p);
    friend std::ostream& operator<<(std::ostream& os, const Point& p);
};

class LocatedVector : public SCoord
{
	public:
	LocatedVector() { r = theta = phi = origin.x = origin.y = origin.z = 0; }
	LocatedVector(SCoord sc) { copy(sc); }
	void copy(SCoord sc) { r=sc.r; theta=sc.theta; phi=sc.phi; }
    Point origin;
};

class Rotation
{
	public:
    SCoord v;
    float a=0;
    Rotation()	{	a = 0;	};

    Rotation add(Rotation rot)	{	return add(&rot);	};
    Rotation add(Rotation* rot);

    char* printable();
};

class LocRotation : public Rotation
{
	public:
    Point origin;

    LocRotation();
    LocRotation(LocatedVector lv);
    LocatedVector get_lv();
};

class Atom;

class Tug
{
	public:
    SCoord vec;
    Rotation rot;

    Tug() { ; }
    Tug(Atom* atom, Point molcen, SCoord pull);

    Tug add(Tug t);
};

class Bond;
class InteratomicForce;
class Molecule;
class AminoAcid;
class Protein;
class Region;
class Rotation;

union Star
{
    int n;
    char* psz;
    const char* cpsz;
    void* p;
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
};

Point average_of_points(Point* points, int count);

float find_angle(float dx, float dy);
float find_angle_delta(float a1, float a2);
float find_3d_angle(Point* A, Point* B, Point* source);
float find_3d_angle(Point A, Point B, Point source);
float find_angle_along_vector(Point* pt1, Point* pt2, Point* source, SCoord* v);
float find_angle_along_vector(Point pt1, Point pt2, Point source, SCoord v);

Point rotate3D(Point* point, Point* source, SCoord* axis, float theta);
Point rotate3D(Point point, Point source, SCoord axis, float theta);
Point rotate3D(Point* point, Point* source, Rotation* rot);

Rotation align_points_3d(Point* point, Point* align, Point* center);
Rotation* align_2points_3d(Point* point1, Point* align1, Point* point2, Point* align2, Point* center);

SCoord compute_normal(Point* pt1, Point* pt2, Point* pt3);
SCoord compute_normal(Point pt1, Point pt2, Point pt3);

float are_points_planar(Point p1, Point p2, Point p3, Point p4);

float sphere_intersection(float r1, float r2, float d);

SCoord v_from_pt_sub(Point distal, Point reference);

// Misc. functions not Point-related but not enough to warrant their own .h and .cpp files:
int in_array(void* needle, void** haystack);
int in_array(int needle, int* haystack);
int in_array(Star needle, Star* haystack);
Star* array_unique(Star* input_array);
char** chop_spaced_fields(char* line, char separator = ' ');
float polygon_radius(float side_length, int num_sides);
int greek_from_aname(const char* aname);
int randsgn();
float frand(float min, float max);
float Pearson_correlation(float* xarr, float* yarr, int length);
enum STR_PAD {STR_PAD_RIGHT, STR_PAD_LEFT, STR_PAD_BOTH};
std::string str_pad(const std::string &str, int pad_length, std::string pad_string=" ", STR_PAD pad_type=STR_PAD_RIGHT);
std::string cardinality_printable(float card);

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

std::ostream& operator<<(std::ostream& os, const Point& p);
std::ostream& operator<<(std::ostream& os, const SCoord& v);
std::ostream& operator<<(std::ostream& os, const Rotation& r);


const float tetrahedral_angle = acos(-1.0/3);
extern const char* Greek;

extern std::ofstream *debug;

extern bool last_iter;

#endif

