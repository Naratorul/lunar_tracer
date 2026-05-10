#pragma once
#include "vec3.h"
#include "ray.h"
#include "brdf.h"
#include "atmosphere.h"
#include "mesh.h"
#include "ephemeris.h"
#include <cmath>
#include <algorithm>

static constexpr double SUN_HALF_ANGLE = 0.00436;

static double rng(unsigned int& s) {
    s ^= s<<13; s ^= s>>17; s ^= s<<5;
    return (s & 0xFFFFFF) / (double)0x1000000;
}

double shadow_vis(const pt3& p, const vec3& N, const vec3& sun_dir,
                  const Mesh& mesh, unsigned int& s)
{
    vec3 up  = (std::fabs(sun_dir.z()) < 0.9) ? vec3(0,0,1) : vec3(1,0,0);
    vec3 t1  = norm(cross(sun_dir, up));
    vec3 t2  = cross(sun_dir, t1);
    pt3 orig = p + N * 0.5;
    int lit = 0;
    for (int i=0;i<4;++i) {
        double u,v,r2;
        do { u=2*rng(s)-1; v=2*rng(s)-1; r2=u*u+v*v; } while(r2>1.0);
        vec3 sd = norm(sun_dir + SUN_HALF_ANGLE*(u*t1 + v*t2));
        hit rec;
        double tmax = 1e9;
        if (!mesh.intersect(ray(orig,sd), 0.01, tmax, rec, pt3(0,0,0), 1.0)) ++lit;
    }
    return lit / 4.0;
}

color trace(const ray& r, const Mesh& mesh, const Ephem& ep,
            const Hapke& mat, unsigned int& s)
{
    hit rec;
    double tmax = 1e9;
    if (!mesh.intersect(r, 0.01, tmax, rec, pt3(0,0,0), 1.0))
        return {0,0,0};

    const vec3 N   = rec.n;
    const vec3 V   = norm(-r.d);
    const vec3 L   = ep.sun_dir;
    color       tex = mesh.tex.sample(rec.u, rec.v);

    color direct{0,0,0};
    double mu0 = dot(N, L);
    if (mu0 > 0) {
        double vis = shadow_vis(rec.p, N, L, mesh, s);
        if (vis > 0)
            direct = mat.eval(N, L, V, tex) * (mu0 * vis);
    }

    Earthshine es(ep.earth_dir);
    color earth = es.eval(N, L, tex);

    return direct + earth;
}
