# Builds the PyXSCAT shared library into the Python package.
#
#   make                       gfortran (default)
#   make FC=ifx                Intel oneAPI ifx (also ifort); load the oneAPI environment first,
#                              e.g. `source /opt/intel/oneapi/setvars.sh` or `module load intel`
#   make FC=ifx EXTRA=-xHost   add flags (here: optimise for the build machine's CPU)
#   make FFLAGS="..."          replace the default flags entirely (command line only: FFLAGS and
#                              LDFLAGS set in the environment, e.g. by conda, are ignored)
#
# `pip install .` runs this Makefile with FC and EXTRA from the environment, e.g.
#   FC=ifx pip install .        FC=ifx EXTRA=-xHost pip install .
ifeq ($(origin FC),default)
  FC := gfortran
endif
COMPILER := $(notdir $(firstword $(FC)))

ifneq (,$(filter ifx ifort %-ifx %-ifort,$(COMPILER)))
  # Intel: value-safe floating point (the default fp-model fast=1 reassociates sums and uses
  # approximate divisions, which the 1e-12 tests would notice); Intel runtime linked statically
  # into the library so only the OpenMP runtime (libiomp5) is needed at run time
  DEF_FFLAGS  := -O3 -fPIC -qopenmp -fp-model precise
  DEF_LDFLAGS := -static-intel
  MODDIR       = -module $(BUILD)
else
  DEF_FFLAGS  := -O3 -fPIC -fopenmp
  DEF_LDFLAGS :=
  MODDIR       = -J$(BUILD)
endif
ifneq ($(origin FFLAGS),command line)
  FFLAGS := $(DEF_FFLAGS)
endif
ifneq ($(origin LDFLAGS),command line)
  LDFLAGS := $(DEF_LDFLAGS)
endif
EXTRA ?=

# The bitwise RDM builder counts bits with popcnt; x86-64 compilers emit the hardware
# instruction only when asked (every x86-64 CPU since 2008 has it), otherwise a library call
# that is 7x slower.  Only xs_rdm.f90 gets this flag.
# Decided by the compiler's target, not by `uname -m`: on Apple silicon a shell or Python under
# Rosetta reports x86_64 while Homebrew's gfortran builds for arm64 (which has no -mpopcnt).
FCMACH := $(shell $(FC) -dumpmachine 2>/dev/null)
ifeq ($(strip $(FCMACH)),)
  FCMACH := $(shell uname -m)
endif
ifneq (,$(findstring x86_64,$(FCMACH)))
  ifneq (,$(filter ifx ifort %-ifx %-ifort,$(COMPILER)))
    POPCNT := -march=x86-64-v2
  else
    POPCNT := -mpopcnt
  endif
endif

UNAME   := $(shell uname -s)
ifeq ($(UNAME),Darwin)
  LIBEXT := dylib
else
  LIBEXT := so
endif
SRC   := src/fortran
BUILD := build
OBJS  := $(BUILD)/xs_bessel_j0.o $(BUILD)/xs_p0.o $(BUILD)/xs_hermite.o $(BUILD)/xs_radial.o $(BUILD)/xs_kernel_j0.o \
         $(BUILD)/xs_pairft.o $(BUILD)/xs_nuclear.o $(BUILD)/xs_elastic.o $(BUILD)/xs_rdm.o $(BUILD)/xs_capi.o
LIB   := python/pyxscat/libpyxscat.$(LIBEXT)

all: $(LIB)

# Everything is rebuilt when the compiler, the flags or the content of any source file change.
# (File times alone are not enough: unpacking a newer tarball over an old build gives sources
# that look older than the stale objects, and switching FC would silently reuse the old ones.)
STAMP     := $(BUILD)/.build-stamp
STAMP_NOW := $(shell echo "$(FC) $(FFLAGS) $(EXTRA) $(LDFLAGS) $(POPCNT)"; cat $(SRC)/*.f90 | cksum)
$(STAMP): FORCE
	@mkdir -p $(BUILD)
	@if [ "$$(cat $@ 2>/dev/null)" != "$(STAMP_NOW)" ]; then \
	    rm -f $(BUILD)/*.o $(BUILD)/*.mod $(LIB); echo "$(STAMP_NOW)" > $@; \
	    echo "fresh build with $(COMPILER)"; fi

$(BUILD)/%.o: $(SRC)/%.f90 $(STAMP)
	$(FC) $(FFLAGS) $(EXTRA) $(MODDIR) -c $< -o $@

$(BUILD)/xs_rdm.o: $(SRC)/xs_rdm.f90 $(STAMP)
	$(FC) $(FFLAGS) $(POPCNT) $(EXTRA) $(MODDIR) -c $< -o $@

$(BUILD)/xs_kernel_j0.o: $(BUILD)/xs_bessel_j0.o $(BUILD)/xs_p0.o $(BUILD)/xs_hermite.o $(BUILD)/xs_radial.o
$(BUILD)/xs_nuclear.o: $(BUILD)/xs_pairft.o $(BUILD)/xs_radial.o
$(BUILD)/xs_elastic.o: $(BUILD)/xs_pairft.o $(BUILD)/xs_radial.o
$(BUILD)/xs_capi.o: $(BUILD)/xs_kernel_j0.o $(BUILD)/xs_pairft.o $(BUILD)/xs_nuclear.o $(BUILD)/xs_elastic.o $(BUILD)/xs_rdm.o

$(LIB): $(OBJS) $(STAMP)
	$(FC) $(FFLAGS) $(EXTRA) -shared $(OBJS) $(LDFLAGS) -o $@
	@echo "built $@ with $(COMPILER)"

test: $(LIB)
	python -m pytest -q tests

clean:
	rm -rf $(BUILD) python/pyxscat/libpyxscat.*

.PHONY: all test clean FORCE
FORCE:
