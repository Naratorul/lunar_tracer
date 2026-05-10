#pragma once
#include <cmath>
#include <algorithm>
#include "vec3.h"

static constexpr double INV4PI = 1.0/(4.0*3.14159265358979323846);

struct Hapke {
    double B0    = 1.1;
    double hs    = 0.06;
    double g     = -0.28;
    double wbase = 1.0;

    double phase(double cos_g) const {
        double g2 = g*g;
        return (1.0-g2) / std::pow(std::max(1e-8, 1.0+g2-2.0*g*cos_g), 1.5);
    }

    double surge(double phase_angle) const {
        return B0 / (1.0 + std::tan(0.5*phase_angle) / hs);
    }

    double H(double x, double w) const {
        if (x <= 0) return 1.0;
        double gamma = std::sqrt(std::max(0.0, 1.0 - w));
        return (1.0 + 2.0*x) / (1.0 + 2.0*x*gamma);
    }

    color eval(const vec3& N, const vec3& L, const vec3& V, const color& tex) const {
        double mu0 = std::max(0.0, dot(N,L));
        double mu  = std::max(0.0, dot(N,V));
        if (mu0 <= 0 || mu <= 0) return {0,0,0};

        double cos_g = std::clamp(dot(-L,-V), -1.0, 1.0);
        double g_ang = std::acos(cos_g);
        double pg    = phase(cos_g);
        double Bg    = surge(g_ang);
        double LS    = mu0 / (mu0 + mu);

        color res;
        for (int ch=0;ch<3;++ch) {
            double w = std::clamp(std::pow(std::max(0.0,tex[ch]), 2.2) * wbase * 2.0, 0.001, 0.990);
            res[ch] = std::max(0.0,
                (w * INV4PI) * LS * ((1.0+Bg)*pg + H(mu0,w)*H(mu,w) - 1.0));
        }
        return res;
    }
};
