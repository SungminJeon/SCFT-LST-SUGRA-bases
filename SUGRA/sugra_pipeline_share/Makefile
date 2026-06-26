CXX      = g++
CXXFLAGS = -std=c++17 -O2 -I.

# OpenMP — try plain -fopenmp first (gcc, Homebrew LLVM); fall back to
# Homebrew libomp on Apple Clang where -fopenmp is not directly supported.
ifeq ($(shell $(CXX) -fopenmp -dM -E -x c++ /dev/null > /dev/null 2>&1 && echo yes),yes)
  OMP_CXXFLAGS = -fopenmp
  OMP_LDFLAGS  = -fopenmp
else ifneq ($(wildcard /opt/homebrew/opt/libomp/lib/libomp.dylib),)
  OMP_CXXFLAGS = -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
  OMP_LDFLAGS  = -L/opt/homebrew/opt/libomp/lib -lomp
else ifneq ($(wildcard /usr/local/opt/libomp/lib/libomp.dylib),)
  OMP_CXXFLAGS = -Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include
  OMP_LDFLAGS  = -L/usr/local/opt/libomp/lib -lomp
else
  $(warning OpenMP not detected — building serial.)
  OMP_CXXFLAGS =
  OMP_LDFLAGS  =
endif
CXXFLAGS += $(OMP_CXXFLAGS)

# Auto-detect Eigen3 path
# macOS (Homebrew Apple Silicon)
ifneq ($(wildcard /opt/homebrew/include/eigen3),)
  CXXFLAGS += -I/opt/homebrew/include/eigen3
# macOS (Homebrew Intel)
else ifneq ($(wildcard /usr/local/include/eigen3),)
  CXXFLAGS += -I/usr/local/include/eigen3
# Linux
else ifneq ($(wildcard /usr/include/eigen3),)
  CXXFLAGS += -I/usr/include/eigen3
endif

# Override: make EIGEN=/your/path
ifdef EIGEN
  CXXFLAGS += -I$(EIGEN)
endif

LDFLAGS  = $(OMP_LDFLAGS)

# Common object files (shared across all phases)
COMMON_SRC = Tensor.C TopoLineCompact_enhanced.cpp Topology_enhanced.cpp
COMMON_OBJ = $(COMMON_SRC:.cpp=.o)

# Fix Tensor.C → Tensor.o
COMMON_OBJ := Tensor.o TopoLineCompact_enhanced.o Topology_enhanced.o

all: gen_sugra_phase1 gen_sugra_phase2 gen_sugra_phase3 gen_sugra_nhc_ext gen_sugra_nhc_ext_phase2 gen_sugra_special_m1 cat2tex

Tensor.o: Tensor.C Tensor.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

TopoLineCompact_enhanced.o: TopoLineCompact_enhanced.cpp TopoLineCompact_enhanced.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

Topology_enhanced.o: Topology_enhanced.cpp Topology_enhanced.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

gen_sugra_phase1: gen_sugra_phase1.cpp $(COMMON_OBJ) sugra_generator.h theory_sugra.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(COMMON_OBJ) $(LDFLAGS)

gen_sugra_phase2: gen_sugra_phase2.cpp $(COMMON_OBJ) sugra_generator.h theory_sugra.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(COMMON_OBJ) $(LDFLAGS)

gen_sugra_phase3: gen_sugra_phase3.cpp $(COMMON_OBJ) sugra_generator.h theory_sugra.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(COMMON_OBJ) $(LDFLAGS)

gen_sugra_nhc_ext: gen_sugra_nhc_ext.cpp $(COMMON_OBJ) sugra_generator.h theory_sugra.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(COMMON_OBJ) $(LDFLAGS)

gen_sugra_nhc_ext_phase2: gen_sugra_nhc_ext_phase2.cpp $(COMMON_OBJ) sugra_generator.h theory_sugra.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(COMMON_OBJ) $(LDFLAGS)

# Special single-external phase1: gauge-less (-1) on kernel-component-1 curves
gen_sugra_special_m1: gen_sugra_special_m1.cpp $(COMMON_OBJ) sugra_generator.h theory_sugra.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(COMMON_OBJ) $(LDFLAGS)

cat2tex: cat2tex.cpp $(COMMON_OBJ) sugra_generator.h theory_sugra.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(COMMON_OBJ) $(LDFLAGS)

clean:
	rm -f *.o gen_sugra_phase1 gen_sugra_phase2 gen_sugra_phase3 gen_sugra_nhc_ext gen_sugra_nhc_ext_phase2 gen_sugra_special_m1 cat2tex

# ── Quick test: T<=5 ──
test: all
	@echo "=== Phase 1: 1-ext catalog (T<=5) ==="
	./gen_sugra_phase1 unified.cat 5 0
	@echo ""
	@echo "=== Phase 2: 2-ext ==="
	./gen_sugra_phase2 cat_1ext
	@echo ""
	@echo "=== Phase 2: 3-ext ==="
	./gen_sugra_phase2 cat_2ext
	@echo ""
	@echo "=== Phase 3: gluing ==="
	./gen_sugra_phase3 cat_1ext 193

# ── Full run: T<=9 ──
run: all
	@echo "=== Phase 1 ==="
	./gen_sugra_phase1 unified.cat 9 0
	@echo ""
	@echo "=== Phase 2: 2-ext ==="
	./gen_sugra_phase2 cat_1ext
	@echo ""
	@echo "=== Phase 3 ==="
	./gen_sugra_phase3 cat_1ext 193

# ── LaTeX compilation ──
pdf: sugra_phase2_T0_9.pdf sugra_phase3.pdf

sugra_phase2_T0_9.pdf: sugra_phase2_T0_9.tex
	pdflatex -interaction=nonstopmode $<

sugra_phase3.pdf: sugra_phase3.tex
	pdflatex -interaction=nonstopmode $<

.PHONY: all clean test run pdf
