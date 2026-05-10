# lunar_tracer v2 (AfterDemo)

A physically-based path tracer for the Moon. Given a date and time, it computes
accurate Sun and Earth directions using low-precision analytic ephemerides, loads a
binary triangle mesh plus a raw texture, and renders a 1920×1080 image via a
Hapke BRDF with soft shadows, earthshine, and atmospheric refraction of the view ray.

Output is a 24-bit BMP (`lunar.bmp`) plus a plain-text sidecar (`lunar.meta`)
containing the phase angle, illumination fraction, sub-solar coordinates, and the
render parameters used.

![lunar render example](image.bmp)

---

## What's inside

```
lunar_tracer/
├── main.cpp        — CLI, config struct, render loop, file output
├── renderer.h      — ray–scene intersection, shadow sampling, shade + earthshine
├── ephemeris.h     — analytic Sun/Moon positions (VSOP87 subset + ELP2000 terms)
├── brdf.h          — Hapke photometric model (opposition surge, Henyey-Greenstein)
├── atmosphere.h    — atmospheric refraction of view rays, Earthshine struct
├── camera.h        — pinhole camera with refracted view direction
├── mesh.h          — binary triangle mesh + raw texture loader
├── ray.h           — ray type and hit record
├── vec3.h          — 3-component vector math
├── image.h         — floating-point framebuffer, BMP writer
├── Makefile        — g++ / clang++ build rules, OpenMP optional
├── moon_mesh.bin   — pre-tessellated sphere (or LOLA-derived heightmap)
└── moon_tex.raw    — 4096×2048 albedo map (raw RGB bytes, no header)
```

---

## Building

You need a C++17 compiler. OpenMP is optional but makes a big difference on
multi-core machines.

```bash
# with OpenMP (recommended)
make

# single-threaded fallback
make noomp
```

The resulting binary is `./lunar_tracer`.

Tested with GCC 12 and Clang 16 on Linux. macOS users: install libomp via
Homebrew and replace `-fopenmp` with `-Xclang -fopenmp` in the Makefile.

---

## Running

The defaults render tonight at 22:00 UTC, 1920×1080, 8 samples-per-pixel:

```bash
./lunar_tracer
```

That produces `lunar.bmp` and `lunar.meta` in the current directory.

### Common flags

| Flag | Default | Description |
|---|---|---|
| `--date Y M D` | 2026 5 10 | Observation date (UTC) |
| `--time H M` | 22 00 | Observation time (UTC) |
| `--width W` | 1920 | Output width in pixels |
| `--height H` | 1080 | Output height in pixels |
| `--spp N` | 8 | Samples per pixel |
| `--fov DEG` | 12.0 | Camera field of view |
| `--alt DEG` | 40.0 | Moon altitude above horizon (for refraction) |
| `--mesh FILE` | moon_mesh.bin | Binary triangle mesh |
| `--tex FILE` | moon_tex.raw | Raw 4096×2048 RGB texture |
| `--out STEM` | lunar | Output filename stem |

### A few examples

```bash
# Full moon on May 10 2026, high-quality render
./lunar_tracer --date 2026 5 10 --time 22 00 --spp 32 --out fullmoon

# Crescent moon, tight crop, sunrise terminator angle
./lunar_tracer --date 2026 6 2 --time 20 30 --fov 6.0 --spp 64

# 4K, 16 spp, moon near zenith (minimal refraction)
./lunar_tracer --width 3840 --height 2160 --spp 16 --alt 85
```

---

## How it works

### Ephemeris

`Ephem::compute()` converts the date/time to a Julian Day Number and derives the
Sun's ecliptic longitude via the standard C1+C2+C3 equation-of-center expansion
(roughly VSOP87 truncated to the first three terms). The Moon's longitude, latitude,
and distance come from the 10-term ELP2000 subset — good to about 0.3° over the
next few decades, more than enough for lighting direction. Both vectors are rotated
from ecliptic to equatorial coordinates using the mean obliquity, then normalized.

The `.meta` file written at the end logs the exact vectors so you can verify
them against JPL Horizons if needed.

### BRDF

The reflectance model is Hapke (1981/2002) with:

- **Henyey-Greenstein phase function** — single-lobe, g = −0.28 (backscattering)
- **Opposition surge** — B₀ = 1.1, angular half-width hs = 0.06 rad
- **Multiple-scattering approximation** — H-function pairs for incidence and emission
- **Per-channel single-scattering albedo** sampled from the texture (gamma-decoded,
  then scaled by wbase)

The shadow check fires 4 stratified cone rays spanning the Sun's angular diameter
(0.5°) and averages the visibility. This gives penumbra on crater walls without
being too expensive at low spp.

### Earthshine

The dark limb receives a dim blue-white contribution from `Earthshine::eval()`.
The magnitude is fixed at (0.012, 0.017, 0.030) — roughly matching measured values
for a half-illuminated Earth. It falls off smoothly across the terminator using the
dot product between the surface normal and the Earth direction.

### Atmospheric refraction

`refract_view()` bends each outgoing camera ray upward by the standard
Bennett refraction formula (the same one used in almanacs). At 40° altitude the
correction is about 1.3 arcminutes — small but visible if you compare to a raw
pinhole render side by side.

---

## Output files

`lunar.bmp` — 24-bit bottom-up BMP, sRGB. Open directly in any image viewer.
No gamma is applied by the renderer; the Hapke model produces linear-light values
and the BMP is written with a simple `round(clamp(v,0,1)*255)` per channel.
You may want to run it through a tone-mapper or apply a gentle gamma in post.

`lunar.meta` — key-value text file:

```
date 2026-05-10
time 22:00 UTC
phase_deg 157.23
illum_pct 96.4
sub_solar_lat  1.43
sub_solar_lon -71.22
sun_dir  0.9312 -0.3641  0.0250
earth_dir -0.9998  0.0183  0.0007
fov_deg 12.0
spp 8
resolution 1920x1080
```

---

## Performance

On a 6-core laptop (Ryzen 5 6600H), a 1920×1080 render at `--spp 8` finishes in
roughly 4 seconds with OpenMP. Doubling spp doubles render time linearly.
The BVH in `mesh.h` cuts ray cost to O(log N) per intersection — for the default
~20 k-triangle sphere this barely matters, but it makes a real difference if you
swap in a high-resolution LOLA mesh with a few million triangles.

---

## Swapping in real terrain

The `--mesh` flag accepts any binary triangle mesh in the same format: a 4-byte
uint32 face count, then N × 9 floats (v0.xyz, v1.xyz, v2.xyz) followed by
N × 6 floats (uv0, uv1, uv2 as u0 v0 u1 v1 u2 v2). A converter for
[LOLA LDEM_8](https://pds-geosciences.wustl.edu/missions/lro/lola.htm) DEMs is not
included but the file format is straightforward to target from Python/numpy.

---

## Known limitations

- No polarization or second-order scattering.
- The ephemeris ignores nutation and aberration — error is under 0.1° for
  dates 2000–2050.
- The BMP writer doesn't embed an ICC profile, so some viewers may interpret the
  linear values as sRGB and look slightly washed out.
- `--alt` below ~5° produces exaggerated refraction; don't render the Moon at
  the horizon and expect a photographic match.

---

## License

Do whatever you want with it. If you publish a render, credit is appreciated but
not required.
