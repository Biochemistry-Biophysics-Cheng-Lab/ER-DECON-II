# erdecon -- ER-Decon I and II command-line programs (core_decon, core2_decon), Priism-free.
#
#   make                      build/bin/core2_decon and build/bin/core_decon
#   make STATIC=1             the same, statically linked (needs libfftw3.a)
#   make test                 unit tests + golden tests (Matlab reference)
#   make check                test + regression suite against the original binaries (~15 s)
#   make install PREFIX=dir   dir/bin/core2_decon, dir/bin/core_decon (e.g. PREFIX=$CONDA_PREFIX)
#   make clean
#
# Needs a C99 compiler with OpenMP and FFTW 3 (double precision, with the
# threads library).  Non-default FFTW location: FFTW_INC=... FFTW_LIB=...
# (FFTW_LIB is also added as the run-time search path, e.g. a conda env's lib/).
# Inside a conda environment with c-compiler and fftw installed:
#   make install PREFIX=$CONDA_PREFIX FFTW_INC=$CONDA_PREFIX/include FFTW_LIB=$CONDA_PREFIX/lib
# macOS/Homebrew (untested): make FFTW_INC=/opt/homebrew/include FFTW_LIB=/opt/homebrew/lib
#
# FFT_THREADED=1 (default) multi-threaded FFTs, as the original Priism build.
# FFTW_DETERMINISTIC=1 (default) plans with FFTW_ESTIMATE so results do not
# depend on machine load; 0 restores the original FFTW_MEASURE planning.

CC       ?= gcc
CSTD     ?= -std=gnu99
OPT      ?= -O2
WARN     ?= -Wall -Wno-unused-variable -Wno-unused-but-set-variable
FFTW_INC ?=
FFTW_LIB ?=
STATIC   ?= 0
FFT_THREADED ?= 1
FFTW_DETERMINISTIC ?= 1
PREFIX   ?= /usr/local

CPPFLAGS = -Isrc -Isrc/recon -Isrc/recon/compat $(if $(FFTW_INC),-I$(FFTW_INC)) -DUSE_OPENMP
CFLAGS   = $(CSTD) $(OPT) $(WARN) -fopenmp
LDFLAGS  = -fopenmp $(if $(FFTW_LIB),-L$(FFTW_LIB) $(if $(filter 1,$(STATIC)),,-Wl,-rpath,$(FFTW_LIB))) \
           $(if $(filter 1,$(STATIC)),-static) $(EXTRA_LDFLAGS)
LDLIBS   = -lfftw3 -lm -lpthread
ifeq ($(FFT_THREADED),1)
  CPPFLAGS += -DFFT_THREADED
  LDLIBS    = -lfftw3_threads -lfftw3 -lm -lpthread
endif
ifeq ($(FFTW_DETERMINISTIC),1)
  CPPFLAGS += -include src/recon/fftw_deterministic.h
endif

B = build

# delivered ER-Decon sources (UCSF), unchanged; see NOTICE.md
SUPPORT   = fft_double.c getmcf.c interp_util.c binary_double_heap.c core_estimate_background.c misc.c misc2d.c
ERDECON2  = core2.c core2_2d.c jacobi2.c jacobi2_2d.c jacobi2_common.c regc.c regc2d.c \
            JacobiNL3Ddu.c JacobiNL2Ddu.c
# ER-Decon I (core, Jacobi, Update, ApplyReg and ApplyComplexTFp* are defined in
# both families, so they are linked into separate programs, as in Priism)
ERDECON1  = core.c jacobi.c applyreg.c
# reconstructions of the files that were never delivered (src/recon)
RECON     = recon/dtests.c recon/volumestruct.c recon/volumestruct2d.c \
            recon/derivatives.c recon/derivatives2d.c recon/solver2014.c \
            recon/Res3D.c recon/Res2D.c recon/LineSearch3D2n.c recon/LineSearch2D2n.c
# new code: Priism-free image I/O (the plugin stub recon/masku_stubs.c is linked into core2_decon only)
IO        = recon/imlite.c
obj = $(patsubst %.c,$(B)/%.o,$(1))

.PHONY: all test unit golden check install clean
all: $(B)/bin/core2_decon $(B)/bin/core_decon

$(B)/%.o: src/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

$(B)/libsupport.a: $(call obj,$(SUPPORT) $(RECON) $(IO))
	ar rcs $@ $^
$(B)/liberdecon2.a: $(call obj,$(ERDECON2))
	ar rcs $@ $^
$(B)/liberdecon1.a: $(call obj,$(ERDECON1))
	ar rcs $@ $^

$(B)/bin/core2_decon: $(B)/core2_decon.o $(B)/recon/masku_stubs.o $(B)/liberdecon2.a $(B)/libsupport.a
	@mkdir -p $(dir $@)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(B)/bin/core_decon: $(B)/core_decon.o $(B)/liberdecon1.a $(B)/libsupport.a
	@mkdir -p $(dir $@)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(B)/run_golden_tests: src/recon/run_golden_tests.c $(B)/JacobiNL3DduTests.o $(B)/JacobiNL2DduTests.o \
                       $(B)/liberdecon2.a $(B)/libsupport.a
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(B)/test_derivatives: src/recon/test_derivatives.c src/recon/derivatives.c src/recon/derivatives2d.c \
                       src/recon/volumestruct.c src/recon/volumestruct2d.c
	@mkdir -p $(B)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

unit: $(B)/test_derivatives
	./$(B)/test_derivatives
golden: $(B)/run_golden_tests
	./$(B)/run_golden_tests
test: unit golden

check: test all
	python3 tests/run_check.py

install: all
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 $(B)/bin/core2_decon $(B)/bin/core_decon $(DESTDIR)$(PREFIX)/bin/

clean:
	rm -rf $(B)
