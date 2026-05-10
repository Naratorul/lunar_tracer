#pragma once
#include <cmath>
#include <algorithm>
#include "vec3.h"
#include "ray.h"

static constexpr double D2R_ATM = 3.14159265358979323846/180.0;

ray refract_view(const ray& r, const vec3& zenith, double alt_deg) {
    double a = std::max(alt_deg, 0.0);
    double R = 1.0/std::tan((a + 7.31/(a+4.4))*D2R_ATM) * D2R_ATM/60.0;
    vec3 V  = norm(r.d);
    vec3 Zp = zenith - dot(zenith,V)*V;
    if (Zp.len() < 1e-8) return r;
    Zp = norm(Zp);
    return ray(r.o, norm(V*std::cos(R) + Zp*std::sin(R)));
}

struct Earthshine {
    vec3  dir;
    color light = {0.012, 0.017, 0.030};

    Earthshine(const vec3& toward_earth) : dir(norm(toward_earth)) {}

    color eval(const vec3& N, const vec3& sun, const color& alb) const {
        double dark = std::clamp(-dot(N,sun)/0.087+0.5, 0.0, 1.0);
        if (dark <= 0) return {0,0,0};
        double ec = std::max(0.0, dot(N,dir));
        return light * ec * dark * alb;
    }
};
