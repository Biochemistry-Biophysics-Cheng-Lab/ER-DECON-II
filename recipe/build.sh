#!/bin/bash
# conda-build script for erdecon: build against the host environment's FFTW and install into $PREFIX.
set -euxo pipefail
make -j"${CPU_COUNT:-2}" CC="${CC}" FFTW_INC="${PREFIX}/include" FFTW_LIB="${PREFIX}/lib"
make install PREFIX="${PREFIX}"
