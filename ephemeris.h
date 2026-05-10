#pragma once
#include <cmath>
#include <algorithm>
#include <string>
#include <cstdio>
#include "vec3.h"

static constexpr double EP_D2R = 3.14159265358979323846/180.0;
static constexpr double EP_R2D = 180.0/3.14159265358979323846;

static double mod360(double d) { d=std::fmod(d,360); return d<0?d+360:d; }

static double jd(int Y,int M,int D,double h) {
    if (M<=2){Y--;M+=12;}
    int A=Y/100, B=2-A+A/4;
    return std::floor(365.25*(Y+4716))+std::floor(30.6001*(M+1))+D+B-1524.5+h/24.0;
}

struct Ephem {
    vec3   sun_dir, earth_dir;
    double phase_deg, illum;
    double sub_lat, sub_lon;
    int    year, mon, day, hour, min;

    void compute(int Y,int M,int D,double h,double m=0) {
        year=Y; mon=M; day=D; hour=(int)h; min=(int)m;
        double J = jd(Y,M,D,h+m/60.0);
        double T = (J-2451545.0)/36525.0;
        double eps = (23.439291-0.013004*T)*EP_D2R;

        double L0=mod360(280.46646+36000.76983*T);
        double Ms=mod360(357.52911+35999.05029*T);
        double Mr=Ms*EP_D2R;
        double C=(1.914602-0.004817*T)*std::sin(Mr)+(0.019993-0.000101*T)*std::sin(2*Mr)+0.000289*std::sin(3*Mr);
        double slon=(L0+C)*EP_D2R;
        double sdist=1.000001018*(1-0.01671123*0.01671123)/(1+0.01671123*std::cos((Ms+C)*EP_D2R));

        double D_=mod360(297.85036+445267.111480*T), Dr=D_*EP_D2R;
        double Mm=mod360(357.52772+35999.050340*T),  Mmr=Mm*EP_D2R;
        double Mp=mod360(134.96298+477198.867398*T),  Mpr=Mp*EP_D2R;
        double F =mod360(93.27191 +483202.017538*T),  Fr=F*EP_D2R;

        double sl =  6288774*std::sin(Mpr)
                  +  1274027*std::sin(2*Dr-Mpr)
                  +   658314*std::sin(2*Dr)
                  +   213618*std::sin(2*Mpr)
                  -   185116*std::sin(Mmr)
                  -   114332*std::sin(2*Fr)
                  +    58793*std::sin(2*Dr-2*Mpr)
                  +    57066*std::sin(2*Dr-Mmr-Mpr)
                  +    53322*std::sin(2*Dr+Mpr)
                  +    45758*std::sin(2*Dr-Mmr);
        double sr = -20905355*std::cos(Mpr)
                  -  3699111*std::cos(2*Dr-Mpr)
                  -  2955968*std::cos(2*Dr)
                  -   569925*std::cos(2*Mpr)
                  +    48888*std::cos(Mmr)
                  -   246158*std::cos(2*Dr-2*Mpr);
        double sb =  5128122*std::sin(Fr)
                  +   280602*std::sin(Mpr+Fr)
                  +   277693*std::sin(Mpr-Fr)
                  +   173237*std::sin(2*Dr-Fr);

        double mlon = mod360(218.3164477+481267.88123421*T+sl*1e-6);
        double mlat = sb*1e-6;
        double mdist= 385000.56+sr*1e-3;

        auto ecl=[&](double lo,double la,double r)->vec3{
            double l=lo*EP_D2R, b=la*EP_D2R;
            return {r*std::cos(b)*std::cos(l),
                    r*(std::cos(b)*std::sin(l)*std::cos(eps)-std::sin(b)*std::sin(eps)),
                    r*(std::cos(b)*std::sin(l)*std::sin(eps)+std::sin(b)*std::cos(eps))};
        };

        vec3 mp = ecl(mlon,mlat,mdist);
        vec3 sp = ecl(L0+C, 0.0, sdist*1.495978707e8);

        sun_dir   = norm(sp-mp);
        earth_dir = norm(-mp);
        double ca = std::clamp(dot(sun_dir,earth_dir),-1.0,1.0);
        phase_deg = std::acos(ca)*EP_R2D;
        illum     = 0.5*(1+std::cos(phase_deg*EP_D2R));
        sub_lat   = std::asin(std::clamp(sun_dir.z(),-1.0,1.0))*EP_R2D;
        sub_lon   = std::atan2(sun_dir.y(),sun_dir.x())*EP_R2D;
    }

    void print() const {
        printf("Phase: %.1f  Illum: %.1f%%  Sub-solar: (%.1f, %.1f)\n",
               phase_deg, illum*100, sub_lat, sub_lon);
    }
};
