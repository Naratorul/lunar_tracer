#pragma once
#include "vec3.h"

struct hit {
    pt3    p;
    vec3   n;
    double t;
    double u, v;
    bool   front;

    void set_normal(const vec3& dir, const vec3& outward) {
        front = dot(dir, outward) < 0;
        n     = front ? outward : -outward;
    }
};

struct ray {
    pt3  o;
    vec3 d;

    ray() {}
    ray(const pt3& origin, const vec3& dir) : o(origin), d(dir) {}

    pt3 at(double t) const { return o + t*d; }
};
