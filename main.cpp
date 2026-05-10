#include <cstring>
#include <cstdio>
#include <cmath>
#include "vec3.h"
#include "ray.h"
#include "image.h"
#include "camera.h"
#include "mesh.h"
#include "brdf.h"
#include "ephemeris.h"
#include "renderer.h"

struct Cfg {
    int    year=2026, mon=5, day=10;
    int    hour=22,   min=0;
    int    W=1920,    H=1080;
    int    spp=8;
    double fov=12.0;
    double alt=40.0;
    std::string mesh_bin = "moon_mesh.bin";
    std::string tex_raw  = "moon_tex.raw";
    std::string out      = "lunar";
};

Cfg parse(int argc, char** argv) {
    Cfg c;
    for (int i=1;i<argc;++i) {
        if      (!strcmp(argv[i],"--date")  &&i+3<argc){c.year=atoi(argv[++i]);c.mon=atoi(argv[++i]);c.day=atoi(argv[++i]);}
        else if (!strcmp(argv[i],"--time")  &&i+2<argc){c.hour=atoi(argv[++i]);c.min=atoi(argv[++i]);}
        else if (!strcmp(argv[i],"--width") &&i+1<argc) c.W=atoi(argv[++i]);
        else if (!strcmp(argv[i],"--height")&&i+1<argc) c.H=atoi(argv[++i]);
        else if (!strcmp(argv[i],"--spp")   &&i+1<argc) c.spp=atoi(argv[++i]);
        else if (!strcmp(argv[i],"--fov")   &&i+1<argc) c.fov=atof(argv[++i]);
        else if (!strcmp(argv[i],"--alt")   &&i+1<argc) c.alt=atof(argv[++i]);
        else if (!strcmp(argv[i],"--mesh")  &&i+1<argc) c.mesh_bin=argv[++i];
        else if (!strcmp(argv[i],"--tex")   &&i+1<argc) c.tex_raw=argv[++i];
        else if (!strcmp(argv[i],"--out")   &&i+1<argc) c.out=argv[++i];
    }
    return c;
}

int main(int argc, char** argv) {
    Cfg cfg = parse(argc, argv);

    Ephem ep;
    ep.compute(cfg.year, cfg.mon, cfg.day, cfg.hour, cfg.min);
    ep.print();

    Mesh mesh;
    mesh.load_bin(cfg.mesh_bin);
    mesh.tex.load_raw(cfg.tex_raw, 4096, 2048);
    printf("Texture: %dx%d\n", mesh.tex.w, mesh.tex.h);

    Hapke mat;
    pt3  center(0,0,0);
    double dist = 10.0 * R_MOON;
    pt3 cam_pos = center + dist * ep.earth_dir;

    vec3 cam_up = {0,0,1};
    if (std::fabs(dot(cam_up, ep.earth_dir)) > 0.99) cam_up = {0,1,0};

    Camera cam(cfg.fov, (double)cfg.W/cfg.H, cam_pos, center, cam_up, cfg.alt);
    Image  img(cfg.W, cfg.H);

    printf("Rendering %dx%d  spp=%d  fov=%.1f  dist=%.0fkm\n",
           cfg.W, cfg.H, cfg.spp, cfg.fov, dist);

    int done = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,4) shared(done)
#endif
    for (int y=0;y<cfg.H;++y) {
        unsigned int s = (unsigned int)(y*1000003u+17u);
        for (int x=0;x<cfg.W;++x) {
            color acc{0,0,0};
            for (int k=0;k<cfg.spp;++k) {
                s^=s<<13; s^=s>>17; s^=s<<5; double jx=(s&0xFFFFFF)/(double)0x1000000;
                s^=s<<13; s^=s>>17; s^=s<<5; double jy=(s&0xFFFFFF)/(double)0x1000000;
                double u = (x+jx)/(cfg.W-1);
                double v = (y+jy)/(cfg.H-1);
                acc += trace(cam.get(u,v), mesh, ep, mat, s);
            }
            img.at(x,y) = acc / (double)cfg.spp;
        }
#ifdef _OPENMP
#pragma omp atomic
#endif
        ++done;
        if (y%40==0){ printf("\r  %.1f%%", 100.0*done/cfg.H); fflush(stdout); }
    }
    printf("\r  Done.         \n");

    std::string bmp = cfg.out + ".bmp";
    printf("Writing %s\n", bmp.c_str());
    img.save_bmp(bmp);

    std::string meta = cfg.out + ".meta";
    FILE* mf = fopen(meta.c_str(), "w");
    if (mf) {
        fprintf(mf, "date %04d-%02d-%02d\n", cfg.year, cfg.mon, cfg.day);
        fprintf(mf, "time %02d:%02d UTC\n", cfg.hour, cfg.min);
        fprintf(mf, "phase_deg %.2f\n", ep.phase_deg);
        fprintf(mf, "illum_pct %.1f\n", ep.illum*100.0);
        fprintf(mf, "sub_solar_lat %.2f\n", ep.sub_lat);
        fprintf(mf, "sub_solar_lon %.2f\n", ep.sub_lon);
        fprintf(mf, "sun_dir %.4f %.4f %.4f\n",
                ep.sun_dir.x(), ep.sun_dir.y(), ep.sun_dir.z());
        fprintf(mf, "earth_dir %.4f %.4f %.4f\n",
                ep.earth_dir.x(), ep.earth_dir.y(), ep.earth_dir.z());
        fprintf(mf, "fov_deg %.1f\n", cfg.fov);
        fprintf(mf, "spp %d\n", cfg.spp);
        fprintf(mf, "resolution %dx%d\n", cfg.W, cfg.H);
        fclose(mf);
    }

    return 0;
}
