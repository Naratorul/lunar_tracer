# ============================================================
#  Makefile — Lunar Path-Tracer
# ============================================================

CXX      = g++
CXXFLAGS = -std=c++17 -O3 -march=native -ffast-math -Wall -Wextra
OMP      = -fopenmp
TARGET   = lunar_tracer
SRC      = main.cpp

# Link math library (Linux)
LIBS = -lm

.PHONY: all noomp clean help

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(OMP) -o $@ $^ $(LIBS)
	@echo "Build successful: ./$(TARGET)"

# Build without OpenMP (single-threaded)
noomp: $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $^ $(LIBS)
	@echo "Build successful (no OpenMP): ./$(TARGET)"

clean:
	rm -f $(TARGET) *.ppm

help:
	@echo ""
	@echo "Usage:"
	@echo "  make              Build with OpenMP (parallel)"
	@echo "  make noomp        Build without OpenMP"
	@echo "  make clean        Remove binaries and .ppm files"
	@echo ""
	@echo "Run examples:"
	@echo ""
	@echo "  # Tonight at 23:00 UTC, synthetic craters, default 1920x1080"
	@echo "  ./lunar_tracer --date 2026 5 7 --time 23 00 --preview"
	@echo ""
	@echo "  # With real LOLA heightmap (LDEM_8.IMG, ~2 GB)"
	@echo "  ./lunar_tracer --date 2026 5 7 --time 23 00 \\"
	@echo "                 --lola /data/LDEM_8.IMG \\"
	@echo "                 --width 3840 --height 2160 \\"
	@echo "                 --spp 16 --preview"
	@echo ""
	@echo "  # Tight crop, high-SPP, mare albedo"
	@echo "  ./lunar_tracer --fov 0.6 --spp 64 --alt 55 --preview"
	@echo ""
