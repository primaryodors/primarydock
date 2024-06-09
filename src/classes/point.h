
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>

#include "misc.h"

#ifndef _POINT
#define _POINT

using namespace std;

struct SCoord;

struct Point
{
    double x=0, y=0, z=0;
    float weight = 1;

    Point();
    Point(SCoord v);
    Point(SCoord* v);
    Point(double x, double y, double z);
    Point add(Point add_to)
    {
        return add(&add_to);
    }
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
    float get_3d_distance(const Point reference)
    {
        return get_3d_distance(&reference);
    }
    float get_3d_distance(const Point* reference);
    float get_distance_to_line(const Point a, const Point b);         // Where a and b are the termini of the line.
    Point multiply_3d_distance(const Point* reference, float r_mult);
    bool pt_in_bounding_box(const Point* corner1, const Point* corner2);
    float magnitude();
    void scale(float new_magn);
    void multiply(float multiplier);
    bool fits_inside(Point container);

    std::string printable() const;

    Point& operator=(SCoord v);
    friend std::ostream& operator<<(std::ostream& os, const SCoord& v);
};

struct SCoord
{
    double r=0, theta=0, phi=0;

    SCoord()
    {
        r = theta = phi = 0;
    }
    SCoord(const Point from);
    SCoord(const Point* from);
    SCoord(double r, double theta, double phi);
    SCoord negate()
    {
        Point pt(this);
        pt = pt.negate();
        SCoord v(&pt);
        return v;
    }
    SCoord add(SCoord v)
    {
        return add(&v);
    };
    SCoord add(SCoord* v);

    std::string printable() const;

    SCoord& operator=(Point p);
    friend std::ostream& operator<<(std::ostream& os, const Point& p);
};

struct LocatedVector : public SCoord
{
    LocatedVector()
    {
        r = theta = phi = origin.x = origin.y = origin.z = 0;
    }
    LocatedVector(SCoord sc)
    {
        copy(sc);
    }
    void copy(SCoord sc)
    {
        r=sc.r;
        theta=sc.theta;
        phi=sc.phi;
    }
    Point to_point();
    Point origin;
};

struct Rotation
{
    SCoord v;
    float a=0;
    Rotation()
    {
        a = 0;
    };

    Rotation add(Rotation rot)
    {
        return add(&rot);
    };
    Rotation add(Rotation* rot);

    char* printable();
};

struct LocRotation : public Rotation
{
    Point origin;

    LocRotation();
    LocRotation(LocatedVector lv);
    LocatedVector get_lv();
};

class Atom;

class LocationProbability
{
    protected:
    int edge_size = 0;
    float cell_size = 1;
    float* probabilities = nullptr;
    float uncertainty = 1e9;
    Point center = Point(0,0,0);

    public:
    LocationProbability();
    LocationProbability(Point pt);              // Creates a probability of 100% at that exact point.
    LocationProbability(float cell_size, float spatial_extent);
    ~LocationProbability();

    float get_cell_size() { return cell_size; }
    float get_extent() { return cell_size * 0.5 * edge_size; }
    float probability_at(Point pt);
    void set_probability(Point pt, float new_prob);
    void resample(float new_cell_size);
    LocationProbability compound(Atom* less_uncertain, Atom* more_uncertain);

    protected:
    int index_from_coord(Point pt);
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
Rotation align_points_3d(Point point, Point align, Point center);
Rotation* align_2points_3d(Point* point1, Point* align1, Point* point2, Point* align2, Point* center);

SCoord compute_normal(Point* pt1, Point* pt2, Point* pt3);
SCoord compute_normal(Point pt1, Point pt2, Point pt3);

float are_points_planar(Point p1, Point p2, Point p3, Point p4);
float polygon_radius(float side_length, int num_sides);

float sphere_intersection(float r1, float r2, float d);         // Volume of the lens composed of the two caps.
float sphere_inter_area(float r1, float r2, float d);           // Area of the circle formed by the spheres' intersection.

SCoord v_from_pt_sub(Point distal, Point reference);

std::ostream& operator<<(std::ostream& os, const Point& p);
std::ostream& operator<<(std::ostream& os, const SCoord& v);
std::ostream& operator<<(std::ostream& os, const Rotation& r);

#endif

