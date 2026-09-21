# ---------------------------------------------------------------------------
# Makefile for con2prim_plots
#
# Files:
#   brent_alg.fpp        -> module brent_alg      (needs fypp preprocessing)
#   kastaun.f90      -> module kastaun        (depends on brent_alg)
#   con2prim_plots.f90    -> program con2prim_plots  (depends on kastaun1)
# ---------------------------------------------------------------------------

FC      = gfortran
FYPP    = fypp
FCFLAGS = -O2 -Wall -J$(BUILD_DIR) -I$(BUILD_DIR)
LDFLAGS =

BUILD_DIR = build
TARGET    = con2prim_plots

# ---------------------------------------------------------------------------
# Source -> object mapping
#   brent_alg.fpp is preprocessed into build/brent_alg.f90 first,
#   the other two files are plain Fortran and compiled directly.
# ---------------------------------------------------------------------------

#OBJS = $(BUILD_DIR)/brent_alg.o 
OBJS =       $(BUILD_DIR)/kastaun.o \
       $(BUILD_DIR)/con2prim_plots.o

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $(OBJS) $(LDFLAGS)

# --- fypp preprocessing step -----------------------------------------------
# $(BUILD_DIR)/brent_alg.f90: brent_alg.fpp | $(BUILD_DIR)
# 	$(FYPP) $< $@

$(BUILD_DIR)/kastaun.f90: kastaun.fpp | $(BUILD_DIR)
	$(FYPP) $< $@

# --- compilation rules -------------------------------------------------------
# Order-only dependencies (after the |) enforce compile order for modules,
# since gfortran needs brent_alg.mod before compiling kastaun.f90, etc.

# $(BUILD_DIR)/brent_alg.o: $(BUILD_DIR)/brent_alg.f90 | $(BUILD_DIR)
# 	$(FC) $(FCFLAGS) -c $< -o $@

$(BUILD_DIR)/kastaun.o: $(BUILD_DIR)/kastaun.f90 | $(BUILD_DIR)
	$(FC) $(FCFLAGS) -c $< -o $@

$(BUILD_DIR)/a.o: con2prim_plots.f90 $(BUILD_DIR)/kastaun.o | $(BUILD_DIR)
	$(FC) $(FCFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(TARGET)
