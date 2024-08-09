
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include "point.h"

using namespace std;

Point::Point()
{
    x = y = z = 0;
}

Point::Point(double px, double py, double pz)
{
    x = px;
    y = py;
    z = pz;
}

Point Point::add(Point* add_to)
{
    Point retval;
    retval.x = x + add_to->x;
    retval.y = y + add_to->y;
    retval.z = z + add_to->z;
    if (!retval.x) retval.x = 0;
    if (!retval.y) retval.y = 0;
    if (!retval.z) retval.z = 0;
    return retval;
}

Point Point::add(SCoord* add_to)
{
    Point rel(add_to);
    return add(&rel);
}

Point Point::subtract(const Point subtracted)
{
    Point retval;
    retval.x = x - subtracted.x;
    retval.y = y - subtracted.y;
    retval.z = z - subtracted.z;

    if (isnan(retval.x)) retval.x = 0;
    if (isnan(retval.y)) retval.y = 0;
    if (isnan(retval.z)) retval.z = 0;

    return retval;
}

Point Point::subtract(const Point* subtracted)
{
    Point retval;
    retval.x = x - subtracted->x;
    retval.y = y - subtracted->y;
    retval.z = z - subtracted->z;

    if (isnan(retval.x)) retval.x = 0;
    if (isnan(retval.y)) retval.y = 0;
    if (isnan(retval.z)) retval.z = 0;

    return retval;
}

Point Point::negate()
{
    Point retval;
    retval.x = -x;
    retval.y = -y;
    retval.z = -z;

    if (isnan(retval.x)) retval.x = 0;
    if (isnan(retval.y)) retval.y = 0;
    if (isnan(retval.z)) retval.z = 0;

    return retval;
}

float Point::get_3d_distance(const Point* reference)
{
    float dx = x - reference->x,
          dy = y - reference->y,
          dz = z - reference->z;

    if (isnan(dx)) dx = 0;
    if (isnan(dy)) dy = 0;
    if (isnan(dz)) dz = 0;

    return sqrt(dx*dx + dy*dy + dz*dz);
}

std::string Point::printable() const
{
    std::stringstream buffer;
    buffer << "[" << (0.001 * (int)(x*1000)) << "," << (0.001 * (int)(y*1000)) << "," << (0.001 * (int)(z*1000)) << "]";
    return buffer.str();
}

std::string SCoord::printable() const
{
    std::stringstream buffer;
    buffer << "[φ=" << phi*180/M_PI << ",θ=" << theta*180/M_PI << ",r=" << r << "]";
    return buffer.str();
}

char* Rotation::printable()
{
    char* buffer = new char[256];
    sprintf(buffer, "[φ=%f°, θ=%f°, a=%f°]", v.phi*180/M_PI, v.theta*180/M_PI, a*180/M_PI);
    return buffer;
}

Point Point::multiply_3d_distance(const Point* reference, float r_mult)
{
    Point retval;
    retval.x = (x-reference->x)*r_mult+reference->x;
    retval.y = (y-reference->y)*r_mult+reference->y;
    retval.z = (z-reference->z)*r_mult+reference->z;
    return retval;
}

void Point::scale(float new_magn)
{
    float old_magn = magnitude();
    if (!old_magn) return;
    float multiplier = new_magn / old_magn;
    x *= multiplier;
    y *= multiplier;
    z *= multiplier;
}

void Point::multiply(float m)
{
    x *= m;
    y *= m;
    z *= m;
}

bool Point::pt_in_bounding_box(const Point* corner1, const Point* corner2)
{
    if (x > std::max(corner1->x, corner2->x)) return false;
    if (x < std::min(corner1->x, corner2->x)) return false;
    if (y > std::max(corner1->y, corner2->y)) return false;
    if (y < std::min(corner1->y, corner2->y)) return false;
    if (z > std::max(corner1->z, corner2->z)) return false;
    if (z < std::min(corner1->z, corner2->z)) return false;
    return true;
}

float Point::magnitude()
{
    return sqrt(x*x + y*y + z*z);
}

Point::Point(SCoord* v)
{
    x = v->r * sin(v->phi) *  cos(v->theta);
    z = v->r * cos(v->phi) *  cos(v->theta);
    y = v->r * sin(v->theta);
}

Point::Point(SCoord v)
{
    x = v.r * sin(v.phi) *  cos(v.theta);
    z = v.r * cos(v.phi) *  cos(v.theta);
    y = v.r * sin(v.theta);
}

SCoord::SCoord(const Point* from)
{
    double pxz = from->x*from->x + from->z*from->z;
    r = sqrt(pxz + from->y*from->y);
    pxz = sqrt(pxz);
    phi = find_angle(from->z, from->x);
    theta = find_angle(pxz, from->y);
}

SCoord::SCoord(const Point from)
{
    double pxz = from.x*from.x + from.z*from.z;
    r = sqrt(pxz + from.y*from.y);
    pxz = sqrt(pxz);
    phi = find_angle(from.z, from.x);
    theta = find_angle(pxz, from.y);
}

SCoord::SCoord(double lr, double ltheta, double lphi)
{
    r = lr;
    theta = ltheta;
    phi = lphi;
}

// https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
float Point::get_distance_to_line(Point a, Point b)
{
    float r2 = pow(a.get_3d_distance(b), 2);
    if (!r2) return get_3d_distance(a);
    
    float t = fmax(0, fmin(1,  ((x - a.x) * (b.x - a.x) + (y - a.y) * (b.y - a.y) + (z - a.z) * (b.z - a.z)) / r2));
    Point p(a.x + t * (b.x-a.x), a.y + t * (b.y-a.y), a.z + t * (b.z-a.z));
    
    return get_3d_distance(p);
}


Point average_of_points(Point* points, int count)
{
    int i;
    float x, y, z, sum;

    x=y=z=sum=0;
    for (i=0; i<count; i++)
    {
        x += points[i].weight * points[i].x;
        y += points[i].weight * points[i].y;
        z += points[i].weight * points[i].z;
        sum += points[i].weight;

        #if _dbg_point_avg
        cout << "Adding " << points[i] << " for total [" << x << "," << y << "," << z << "]." << endl;
        #endif
    }

    #if _dbg_point_avg
    cout << "Sum of weights: " << sum << "; return value = [";
    #endif
    x /= sum;
    y /= sum;
    z /= sum;
    #if _dbg_point_avg
    cout << x << "," << y << "," << z << "]." << endl << endl;
    #endif

    Point retval(x, y, z);
    return retval;
}

Point size_of_point_space(Point* points, int count)
{
    int i;
    float x0, y0, z0, x1, y1, z1;

    for (i=0; i<count; i++)
    {
        Point pt = points[i];
        if (!i || pt.x < x0) x0 = pt.x;
        if (!i || pt.y < y0) y0 = pt.y;
        if (!i || pt.z < z0) z0 = pt.z;
        if (!i || pt.x > x1) x1 = pt.x;
        if (!i || pt.y > y1) y1 = pt.y;
        if (!i || pt.z > z1) z1 = pt.z;
    }

    return Point(x1-x0, y1-y0, z1-z0);
}

float find_angle(float dx, float dy)
{
    float angle = atan2(dy,dx);
    if (angle < 0)
    {
        angle += 2 * M_PI;
    }
    return angle;
}

float find_angle_delta(float a1, float a2)
{
    if (a2 > a1)
    {
        float b1 = a1 + M_PI*2;
        if (fabs(b1-a2) < fabs(a2-a1)) return b1-a2;
    }
    else
    {
        float b2 = a2 + M_PI*2;
        if (fabs(b2-a1) < fabs(a1-a2)) return b2-a1;
    }
    return a1 - a2;
}

float find_3d_angle(Point* A, Point* B, Point* source)
{
    if (!source) source = new Point();

    Point lA = A->subtract(source);
    lA.scale(1);
    Point lB = B->subtract(source);
    lB.scale(1);

    // https://stackoverflow.com/questions/1211212/how-to-calculate-an-angle-from-three-points
    float P12 = lA.magnitude();
    float P13 = lB.magnitude();
    float P23 = lA.get_3d_distance(lB);

    float param = (P12*P12 + P13*P13 - P23*P23)/(2 * P12 * P13+.00000000001);
    if (param < -1) param = -1;
    if (param >  1) param =  1;
    float retval = acos(param);
    if (isnan(retval))
    {
        cout << "P12 " << P12 << " P13 " << P13 << " P23 " << P23 << endl;
        throw 0xbad9a9;
    }
    return retval;
}

float find_3d_angle(Point A, Point B, Point source)
{
    Point a=A, b=B, c=source;
    return find_3d_angle(&a, &b, &c);
}

float find_angle_along_vector(Point pt1, Point pt2, Point source, SCoord v)
{
    return find_angle_along_vector(&pt1, &pt2, &source, &v);
}

float find_angle_along_vector(Point* pt1, Point* pt2, Point* source, SCoord* v)
{
    Point vp(v);

    Point lpt1 = pt1->subtract(source);
    Point lpt2 = pt2->subtract(source);

    // Rotate points so v becomes Z axis.
    Point cen;
    Point z(0,0,1);
    Rotation rots = align_points_3d(&vp, &z, &cen);
    Point npt1 = rotate3D(&lpt1, &cen, &rots.v, rots.a);
    Point npt2 = rotate3D(&lpt2, &cen, &rots.v, rots.a);

    // Return the XY angle between the points.
    npt1.z = 0;
    npt2.z = 0;
    // return find_3d_angle(&npt1, &npt2, &z);
    float a1 = find_angle(npt1.x, npt1.y);
    if (a1 > M_PI) a1 -= M_PI*2;
    float a2 = find_angle(npt2.x, npt2.y);
    if (a2 > M_PI) a2 -= M_PI*2;
    return a2 - a1;
}

Point rotate3D(Point* point, Point* source, Rotation* rot)
{
    return rotate3D(point, source, &rot->v, rot->a);
}

Point rotate3D(Point point, Point source, SCoord axis, float theta)
{
    return rotate3D(&point, &source, &axis, theta);
}

Point rotate3D(Point* point, Point* source, SCoord* axis, float theta)
{
    // Originally from http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/

    if (!source)
    {
        Point cen;
        source = &cen;
    }

    if (!axis->r) return *point;
    Point lvec(axis);

    double x = point->x, y = point->y, z = point->z;
    double u = lvec.x / axis->r, v = lvec.y / axis->r, w = lvec.z / axis->r;
    double a, b, c;
    if (source)
    {
        a = source->x;
        b = source->y;
        c = source->z;
    }
    else a = b = c = 0;
    double u2 = u*u;
    double v2 = v*v;
    double w2 = w*w;
    double sint = sin(theta), cost = cos(theta), _1_cost = (1.0 - cost);

    double x1 = (a * (v2+w2) - u * (b*v + c*w - u*x - v*y - w*z)) * _1_cost
               + x * cost
               + (-c*v + b*w - w*y + v*z) * sint;

    double y1 = (b * (u2+w2) - v * (a*u + c*w - u*x - v*y - w*z)) * _1_cost
               + y * cost
               + ( c*u - a*w + w*x - u*z) * sint;

    double z1 = (c * (u2+v2) - w * (a*u + b*v - u*x - v*y - w*z)) * _1_cost
               + z * cost
               + (-b*u + a*v - v*x + u*y) * sint;

    Point pt(x1,y1,z1);
    return pt;
}

Rotation align_points_3d(Point point, Point align, Point center)
{
    return align_points_3d(&point, &align, &center);
}

Rotation align_points_3d(Point* point, Point* align, Point* center)
{
    SCoord n = compute_normal(point, align, center);

    if (n.r < 0.0001)
    {
        Point lpt, lan;

        lpt = *point;
        lan = *align;

        lpt.scale(1);
        lan.scale(1);

        Rotation rot;
        if (lpt.get_3d_distance(&lan) < 0.01)
        {
            rot.v = n;
            rot.a = 0;
            return rot;
        }

        Point pt(0,0,1);
        n = compute_normal(point, align, &pt);
        if (n.r < 0.1)
        {
            pt = Point(0,0,1);
            n = compute_normal(point, align, &pt);
        }

        rot.v = n;
        rot.a = M_PI;
        return rot;
    }

    // Find the 3D angle between pp and pl relative to center.
    float theta = find_3d_angle(point, align, center);
    // cout << " theta = " << theta << " ";

    // Rotate pl positively or negatively along that normal by the found angle, and use the better of the two values.
    Point plus  = rotate3D(point, center, &n,  theta);
    Point minus = rotate3D(point, center, &n, -theta);

    float rplus  = plus.get_3d_distance(align);
    float rminus = minus.get_3d_distance(align);

    float angle;
    if (rplus <= rminus) angle =  theta;
    else                 angle = -theta;

    Rotation rot;
    rot.v = n;
    rot.a = angle;
    return rot;
}

Rotation* align_2points_3d(Point* point1, Point* align1, Point* point2, Point* align2, Point* center)
{
    Rotation* retval = new Rotation[2];
    retval[0] = align_points_3d(point1, align1, center);

    Point point2a = rotate3D(point2, center, &retval[0]);

    SCoord v = v_from_pt_sub(*align1, *center);
    retval[1].v = v;

    // float theta = find_3d_angle(&point2a, align2, center);
    float theta = find_angle_along_vector(&point2a, align2, center, &v);
    if (isnan(theta)) cout << point2a.printable() << ", " << align2->printable() << ", " << center->printable() << endl;

    Point plus  = rotate3D(&point2a, center, &v,  theta);
    Point minus = rotate3D(&point2a, center, &v, -theta);

    float rplus  = align2->get_3d_distance(plus);
    float rminus = align2->get_3d_distance(minus);

    float angle;
    if (rplus <= rminus) angle =  theta;
    else                 angle = -theta;

    retval[1].a = angle;

    return retval;
}

SCoord compute_normal(Point* pt1, Point* pt2, Point* pt3)
{
    Point U = pt2->subtract(pt1);
    Point V = pt3->subtract(pt1);

    Point pt(	U.y * V.z - U.z * V.y,
                U.z * V.x - U.x * V.z,
                U.x * V.y - U.y * V.x
            );
    SCoord v(&pt);
    return v;
}

SCoord compute_normal(Point pt1, Point pt2, Point pt3)
{
    return compute_normal(&pt1, &pt2, &pt3);
}


// https://www.geeksforgeeks.org/program-to-check-whether-4-points-in-a-3-d-plane-are-coplanar/
float are_points_planar(Point p1, Point p2, Point p3, Point p4)
{
    float a1 = p2.x - p1.x ;
    float b1 = p2.y - p1.y ;
    float c1 = p2.z - p1.z ;
    float a2 = p3.x - p1.x ;
    float b2 = p3.y - p1.y ;
    float c2 = p3.z - p1.z ;
    float a  = b1 * c2 - b2 * c1 ;
    float b  = a2 * c1 - a1 * c2 ;
    float c  = a1 * b2 - b1 * a2 ;
    float d  = (- a * p1.x - b * p1.y - c * p1.z) ;

    // Equation of plane is: a*x + b*y + c*z = 0

    // Checking if the 4th point satisfies
    // the above equation.

    // Zero means coplanar; nonzero equals deviation
    return fabs(a * p4.x + b * p4.y + c * p4.z + d);
}

float sphere_intersection(float r1, float r2, float d)
{
    // https://mathworld.wolfram.com/Sphere-SphereIntersection.html
    return M_PI * pow(r1 + r2 - d, 2)
           * (d*d + 2*d*r2 - 3*r2*r2 + 2*d*r1 + 6*r2*r1 - 3*r1*r1)
           / (12*d);
}

float sphere_inter_area(float r1, float r2, float d)
{
    // https://mathworld.wolfram.com/Sphere-SphereIntersection.html
    float a = 1.0 / (d*2) * sqrt(4.0 * d*d * r1*r1 - pow(d*d - r2*r2 + r1*r1, 2));
    return M_PI * a*a;
}

SCoord v_from_pt_sub(Point distal, Point reference)
{
    Point p = distal.subtract(&reference);
    SCoord v(&p);
    return v;
}

float polygon_radius(float side_length, int num_sides)
{
    return side_length/(2.0*sin(M_PI/num_sides));
}

Point& Point::operator=(SCoord v)
{
    Point pt(&v);
    *this = pt;
    return *this;
}

SCoord& SCoord::operator=(Point p)
{
    SCoord v(&p);
    *this = v;
    return *this;
}

SCoord SCoord::add(SCoord* v)
{
    Point pt(v);
    pt.add(this);
    SCoord lv(&pt);
    return lv;
}

Rotation Rotation::add(Rotation* rot)
{
    Point pt(1,0,0), cen;

    float m = 1;
    if ((rot->a + a) > (M_PI/2))
    {
        m = (M_PI/2) / (rot->a + a);
        Rotation r1 = *this;
        Rotation r2 = *rot;

        r1.a *= m;
        r2.a *= m;

        return r1.add(&r2);		// RECURSION!
    }

    Point pt1 = rotate3D(&pt, &cen, this);
    pt1 = rotate3D(&pt1, &cen, rot);

    return align_points_3d(&pt, &pt1, &cen);
}

Tug Tug::add(Tug t)
{
    Tug retval;
    retval.vec = vec.add(t.vec);
    retval.rot = rot.add(t.rot);
    return retval;
}

std::ostream& operator<<(std::ostream& os, const Point& p)
{
    os << p.printable();
    return os;
}

std::ostream& operator<<(std::ostream& os, const SCoord& v)
{
    os << v.printable();
    return os;
}

std::ostream& operator<<(std::ostream& os, const Rotation& r)
{
    os << "[" << (r.v.theta/M_PI*180) << "°, " << (r.v.phi/M_PI*180) << "°, " << (r.a/M_PI*180) << "°]";
    return os;
}

LocRotation::LocRotation()
{
    ;
}

LocRotation::LocRotation(LocatedVector lv)
{
    v.r = lv.r;
    v.theta = lv.theta;
    v.phi = lv.phi;
    a = 0;
    origin = lv.origin;
}

LocatedVector LocRotation::get_lv()
{
    LocatedVector retval;
    retval.r = v.r;
    retval.theta = v.theta;
    retval.phi = v.phi;
    retval.origin = origin;
    return retval;
}

Point LocatedVector::to_point()
{
    return origin.add(this);
}

bool Point::fits_inside(Point c)
{
    if (fabs(x) < fabs(c.x)) return false;
    if (fabs(y) < fabs(c.y)) return false;
    if (fabs(z) < fabs(c.z)) return false;
    return true;
}










