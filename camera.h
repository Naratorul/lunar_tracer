#pragma once
#include "vec3.h"
#include "ray.h"
#include "atmosphere.h"
#include <cmath>

struct Camera {
    pt3  origin;
    vec3 horiz, vert, ll;
    vec3 up;
    double alt;

    Camera(double fov, double aspect, pt3 from, pt3 at, vec3 vup, double moon_alt=40)
        : up(norm(vup)), alt(moon_alt)
    {
        origin = from;
        double h  = std::tan(fov*3.14159265358979323846/180.0*0.5);
        vec3 w = norm(from-at);
        vec3 u = norm(cross(vup,w));
        vec3 v = cross(w,u);
        horiz = 2*h*aspect*u;
        vert  = 2*h*v;
        ll    = origin - horiz*0.5 - vert*0.5 - w;
    }

    ray get(double s, double t) const {
        ray r(origin, norm(ll + s*horiz + t*vert - origin));
        return refract_view(r, up, alt);
    }
};
