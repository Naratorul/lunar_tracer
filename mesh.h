#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <limits>
#include "vec3.h"
#include "ray.h"

static constexpr double R_MOON = 1737.4;
static constexpr double PI     = 3.14159265358979323846;
static constexpr double D2R    = PI / 180.0;

struct Texture {
    int w=0, h=0;
    std::vector<uint8_t> px;

    void load_raw(const std::string& path, int W, int H) {
        w=W; h=H;
        px.resize((size_t)W*H*3);
        std::ifstream f(path, std::ios::binary);
        if (!f) throw std::runtime_error("Cannot open texture: " + path);
        f.read((char*)px.data(), px.size());
    }

    color sample(double u, double v) const {
        u = u - std::floor(u);
        v = std::clamp(v, 0.0, 1.0);
        double fx = u * (w - 1);
        double fy = (1.0 - v) * (h - 1);
        int x0 = (int)fx,  y0 = (int)fy;
        int x1 = std::min(x0+1, w-1);
        int y1 = std::min(y0+1, h-1);
        double tx = fx - x0, ty = fy - y0;
        auto get = [&](int x, int y) -> color {
            size_t idx = ((size_t)y * w + x) * 3;
            return color(px[idx]/255.0, px[idx+1]/255.0, px[idx+2]/255.0);
        };
        color c00=get(x0,y0), c10=get(x1,y0);
        color c01=get(x0,y1), c11=get(x1,y1);
        return (1-ty)*((1-tx)*c00 + tx*c10) + ty*((1-tx)*c01 + tx*c11);
    }
};

struct Vertex {
    float px,py,pz;
    float nx,ny,nz;
    float u,v;
};

struct AABB {
    vec3 mn{ 1e30, 1e30, 1e30};
    vec3 mx{-1e30,-1e30,-1e30};

    void expand(const vec3& p) {
        for (int i=0;i<3;++i){ mn[i]=std::min(mn[i],p[i]); mx[i]=std::max(mx[i],p[i]); }
    }
    void expand(const AABB& o) { expand(o.mn); expand(o.mx); }

    bool hit(const ray& r, double tmin, double tmax) const {
        for (int i=0;i<3;++i) {
            double inv = 1.0 / r.d[i];
            double t0 = (mn[i] - r.o[i]) * inv;
            double t1 = (mx[i] - r.o[i]) * inv;
            if (inv < 0) std::swap(t0,t1);
            tmin = std::max(tmin, t0);
            tmax = std::min(tmax, t1);
            if (tmax <= tmin) return false;
        }
        return true;
    }

    int longest_axis() const {
        vec3 d = mx - mn;
        if (d[0] > d[1] && d[0] > d[2]) return 0;
        if (d[1] > d[2]) return 1;
        return 2;
    }

    vec3 centroid() const { return (mn + mx) * 0.5; }
};

struct BVHNode {
    AABB box;
    int  left  = -1;
    int  right = -1;
    int  tri_start = 0;
    int  tri_count = 0;
};

struct Mesh {
    std::vector<Vertex>   verts;
    std::vector<uint32_t> idx;
    std::vector<int>      tri_order;
    std::vector<BVHNode>  bvh;
    Texture tex;

    void load_bin(const std::string& path) {
        std::ifstream f(path, std::ios::binary);
        if (!f) throw std::runtime_error("Cannot open mesh: " + path);
        uint32_t n_verts, n_tris;
        f.read((char*)&n_verts, 4);
        f.read((char*)&n_tris,  4);
        verts.resize(n_verts);
        f.read((char*)verts.data(), n_verts * sizeof(Vertex));
        idx.resize((size_t)n_tris * 3);
        f.read((char*)idx.data(), idx.size() * 4);
        printf("Mesh: %u verts, %u tris\n", n_verts, n_tris);
        for (auto& v : verts) {
            v.px *= (float)R_MOON;
            v.py *= (float)R_MOON;
            v.pz *= (float)R_MOON;
        }
        build_bvh();
    }

    void build_bvh() {
        int n_tris = (int)(idx.size() / 3);
        tri_order.resize(n_tris);
        std::iota(tri_order.begin(), tri_order.end(), 0);
        bvh.reserve(n_tris * 2);
        build_node(0, n_tris);
        bvh.shrink_to_fit();
        printf("BVH: %d nodes\n", (int)bvh.size());
    }

    AABB tri_box(int ti) const {
        AABB b;
        for (int k=0;k<3;++k) {
            const Vertex& v = verts[idx[(size_t)ti*3+k]];
            b.expand(vec3(v.px,v.py,v.pz));
        }
        return b;
    }

    vec3 tri_centroid(int ti) const {
        vec3 c(0,0,0);
        for (int k=0;k<3;++k) {
            const Vertex& v = verts[idx[(size_t)ti*3+k]];
            c += vec3(v.px,v.py,v.pz);
        }
        return c / 3.0;
    }

    int build_node(int start, int end) {
        int node_idx = (int)bvh.size();
        bvh.push_back(BVHNode{});
        AABB box;
        for (int i=start;i<end;++i) box.expand(tri_box(tri_order[i]));
        bvh[node_idx].box = box;
        int count = end - start;
        if (count <= 8) {
            bvh[node_idx].tri_start = start;
            bvh[node_idx].tri_count = count;
            return node_idx;
        }
        constexpr int N_BINS = 8;
        int   best_axis = -1;
        int   best_split = 0;
        double best_cost = 1e30;
        for (int axis=0;axis<3;++axis) {
            double lo = box.mn[axis], hi = box.mx[axis];
            if (hi - lo < 1e-6) continue;
            double bin_w = (hi - lo) / N_BINS;
            int   cnt[N_BINS]  = {};
            AABB  bbs[N_BINS];
            for (int i=start;i<end;++i) {
                vec3 c = tri_centroid(tri_order[i]);
                int b = std::min((int)((c[axis]-lo)/bin_w), N_BINS-1);
                cnt[b]++;
                bbs[b].expand(tri_box(tri_order[i]));
            }
            for (int s=1;s<N_BINS;++s) {
                AABB L, R;
                int  cL=0, cR=0;
                for (int i=0;i<s;++i){ if(cnt[i]){ L.expand(bbs[i]); cL+=cnt[i]; }}
                for (int i=s;i<N_BINS;++i){ if(cnt[i]){ R.expand(bbs[i]); cR+=cnt[i]; }}
                if (!cL||!cR) continue;
                auto sa = [](const AABB& b)->double{
                    vec3 d=b.mx-b.mn;
                    return 2.0*(d[0]*d[1]+d[1]*d[2]+d[2]*d[0]);
                };
                double cost = sa(L)*cL + sa(R)*cR;
                if (cost < best_cost) { best_cost=cost; best_axis=axis; best_split=s; }
            }
        }
        if (best_axis < 0) {
            bvh[node_idx].tri_start = start;
            bvh[node_idx].tri_count = count;
            return node_idx;
        }
        double lo = box.mn[best_axis], hi = box.mx[best_axis];
        double bin_w = (hi - lo) / N_BINS;
        auto it = std::partition(tri_order.begin()+start, tri_order.begin()+end,
            [&](int ti) {
                vec3 c = tri_centroid(ti);
                int b = std::min((int)((c[best_axis]-lo)/bin_w), N_BINS-1);
                return b < best_split;
            });
        int mid = (int)(it - tri_order.begin());
        if (mid == start || mid == end) { mid = (start+end)/2; }
        int left  = build_node(start, mid);
        int right = build_node(mid,   end);
        bvh[node_idx].left  = left;
        bvh[node_idx].right = right;
        return node_idx;
    }

    bool intersect_tri(int ti, const ray& r, double tmin, double& tmax, hit& rec) const {
        uint32_t ia = idx[(size_t)ti*3+0];
        uint32_t ib = idx[(size_t)ti*3+1];
        uint32_t ic = idx[(size_t)ti*3+2];
        const Vertex& Va = verts[ia];
        const Vertex& Vb = verts[ib];
        const Vertex& Vc = verts[ic];
        vec3 A(Va.px,Va.py,Va.pz);
        vec3 B(Vb.px,Vb.py,Vb.pz);
        vec3 C(Vc.px,Vc.py,Vc.pz);
        vec3 e1 = B - A;
        vec3 e2 = C - A;
        vec3 h  = cross(r.d, e2);
        double det = dot(e1, h);
        if (std::fabs(det) < 1e-12) return false;
        double inv_det = 1.0 / det;
        vec3 s = r.o - A;
        double bu = dot(s, h) * inv_det;
        if (bu < 0 || bu > 1) return false;
        vec3 q = cross(s, e1);
        double bv = dot(r.d, q) * inv_det;
        if (bv < 0 || bu+bv > 1) return false;
        double t = dot(e2, q) * inv_det;
        if (t < tmin || t >= tmax) return false;
        double bw = 1.0 - bu - bv;
        tmax = t;
        rec.t = t;
        rec.p = r.at(t);
        vec3 N = norm(
            bw * vec3(Va.nx,Va.ny,Va.nz) +
            bu * vec3(Vb.nx,Vb.ny,Vb.nz) +
            bv * vec3(Vc.nx,Vc.ny,Vc.nz)
        );
        rec.set_normal(r.d, N);
        rec.u = bw*Va.u + bu*Vb.u + bv*Vc.u;
        rec.v = bw*Va.v + bu*Vb.v + bv*Vc.v;
        return true;
    }

    bool intersect(const ray& r, double tmin, double tmax, hit& rec,
                   const pt3&, double) const
    {
        bool found = false;
        int stack[64];
        int top = 0;
        stack[top++] = 0;
        while (top > 0) {
            int ni = stack[--top];
            const BVHNode& nd = bvh[ni];
            if (!nd.box.hit(r, tmin, tmax)) continue;
            if (nd.left < 0) {
                for (int i=0;i<nd.tri_count;++i)
                    if (intersect_tri(tri_order[nd.tri_start+i], r, tmin, tmax, rec))
                        found = true;
            } else {
                stack[top++] = nd.left;
                stack[top++] = nd.right;
            }
        }
        return found;
    }
};
