#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include "vec3.h"

struct Image {
    int w, h;
    std::vector<color> px;

    Image(int W,int H) : w(W),h(H),px((size_t)W*H,{0,0,0}) {}

    color& at(int x,int y) { return px[(size_t)y*w+x]; }

    void save_bmp(const std::string& path, double gamma=2.2) const {
        double peak=0;
        for (auto& p:px) peak=std::fmax(peak,std::fmax(p[0],std::fmax(p[1],p[2])));
        if (peak<1e-8) peak=1;

        int stride=((w*3+3)/4)*4;
        int pxsz=stride*h;

        std::ofstream f(path,std::ios::binary);
        if (!f) throw std::runtime_error("Can't write: "+path);

        auto u16=[&](uint16_t v){ f.put(v&0xFF); f.put(v>>8); };
        auto u32=[&](uint32_t v){ f.put(v&0xFF);f.put((v>>8)&0xFF);f.put((v>>16)&0xFF);f.put(v>>24); };

        f.put('B'); f.put('M');
        u32(54+pxsz); u16(0); u16(0); u32(54);
        u32(40); u32(w); u32(h); u16(1); u16(24); u32(0); u32(pxsz);
        u32(2835); u32(2835); u32(0); u32(0);

        double ig=1.0/gamma;
        std::vector<uint8_t> pad(stride-w*3,0);
        for (int y=0;y<h;++y) {
            for (int x=0;x<w;++x) {
                const color& c=px[(size_t)y*w+x];
                auto b=[&](double v)->uint8_t{
                    return (uint8_t)(std::pow(std::clamp(v/peak,0.0,1.0),ig)*255+0.5);
                };
                f.put(b(c[2])); f.put(b(c[1])); f.put(b(c[0]));
            }
            f.write((char*)pad.data(),pad.size());
        }
    }
};
