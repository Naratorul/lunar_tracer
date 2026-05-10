#pragma once
#include <cmath>

struct vec3 {
    double e[3];

    vec3() : e{0,0,0} {}
    vec3(double a, double b, double c) : e{a,b,c} {}

    double x() const { return e[0]; }
    double y() const { return e[1]; }
    double z() const { return e[2]; }

    vec3 operator-() const { return {-e[0],-e[1],-e[2]}; }
    double  operator[](int i) const { return e[i]; }
    double& operator[](int i)       { return e[i]; }

    vec3& operator+=(const vec3& v) { e[0]+=v[0]; e[1]+=v[1]; e[2]+=v[2]; return *this; }
    vec3& operator*=(double t)      { e[0]*=t; e[1]*=t; e[2]*=t; return *this; }
    vec3& operator/=(double t)      { return *this *= 1.0/t; }

    double len()  const { return std::sqrt(len2()); }
    double len2() const { return e[0]*e[0]+e[1]*e[1]+e[2]*e[2]; }
};

using pt3   = vec3;
using color = vec3;

inline vec3 operator+(const vec3& a, const vec3& b) { return {a[0]+b[0],a[1]+b[1],a[2]+b[2]}; }
inline vec3 operator-(const vec3& a, const vec3& b) { return {a[0]-b[0],a[1]-b[1],a[2]-b[2]}; }
inline vec3 operator*(const vec3& a, const vec3& b) { return {a[0]*b[0],a[1]*b[1],a[2]*b[2]}; }
inline vec3 operator*(double t, const vec3& v)      { return {t*v[0],t*v[1],t*v[2]}; }
inline vec3 operator*(const vec3& v, double t)      { return t*v; }
inline vec3 operator/(const vec3& v, double t)      { return (1.0/t)*v; }

inline double dot(const vec3& a, const vec3& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
inline vec3 cross(const vec3& a, const vec3& b) {
    return { a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] };
}
inline vec3 norm(const vec3& v) { return v / v.len(); }
