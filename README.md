# ER-Decon II

Source code and a Priism-free build of the ER-Decon deconvolution programs (Arigovindan et al., *PNAS* 2013):

- **`core2_decon`**: ER-Decon II (entropy-regularized deconvolution), including the 2014 line search. This is the program [AR-Decon](https://github.com/yifancheng-ucsf/AR-Decon) runs.
- **`core_decon`**: ER-Decon I.

The front ends and solvers are the original ER-Decon source from the Sedat and Agard laboratories (UCSF). The few source files that were missing have been reconstructed, and a small MRC/Priism file reader/writer replaces Priism's image library. The programs take the same command line, read and write the same files and print the same log as the Priism binaries. Their results match the originals to the last bit on most test cases; otherwise they differ by 1–2 float32 units in the last place in a few voxels per million (see *Verification*).

> **Maintenance status.** This repository keeps the original ER-Decon programs available and buildable. It receives bug fixes only; no new features are planned. New development of the method happens in a separate GPU implementation.

## Install

You need a C99 compiler with OpenMP, GNU make, and FFTW 3 (double precision, with `libfftw3_threads`). Python 3 with NumPy is needed only for `make check`.

**Ubuntu / Debian:**

    sudo apt install build-essential libfftw3-dev python3-numpy
    make -j                      # -> build/bin/core2_decon, build/bin/core_decon
    make check                   # a few seconds; must end with "ALL PASSED"
    sudo make install            # -> /usr/local/bin  (or: make install PREFIX=$HOME/.local)

**conda (any Linux):**

    conda create -n erdecon --override-channels -c conda-forge --strict-channel-priority \
        c-compiler make fftw numpy
    conda activate erdecon
    make -j FFTW_INC=$CONDA_PREFIX/include FFTW_LIB=$CONDA_PREFIX/lib
    make check
    make install PREFIX=$CONDA_PREFIX

Use only conda-forge, as above. With the `defaults` channel also configured, the solver can stall or mix incompatible builds. Conda 23.10 and later use the libmamba solver by default (older versions: add `--solver=libmamba`).

**Static build** (a single file that runs without FFTW installed; needs the static `libfftw3.a` and `libfftw3_threads.a`, which Ubuntu's `libfftw3-dev` provides):

    make -j STATIC=1

A statically linked binary contains FFTW, which is licensed under the GPL (version 2 or later); if you redistribute such a binary, the GPL applies to it.

**macOS** (Homebrew, untested): `brew install fftw gcc`, then
`make CC=gcc-14 FFTW_INC=/opt/homebrew/include FFTW_LIB=/opt/homebrew/lib`.

**Build options:**

| Option | Default | Meaning |
|---|---|---|
| `FFTW_DETERMINISTIC` | `1` | Plan FFTs with `FFTW_ESTIMATE`, so repeated runs give identical results. `0` restores the original `FFTW_MEASURE` planning. |
| `FFT_THREADED` | `1` | Multi-threaded FFTs (`-np=N` threads), as in the Priism build. |
| `STATIC` | `0` | Link statically. |
| `PREFIX` | `/usr/local` | Install location. |

A conda-build recipe is in `recipe/` (`conda build recipe/ --override-channels -c conda-forge`).

## Use

The programs are used exactly like the Priism versions:

    core2_decon measurement.mrc result.mrc otf.otf -lamf=0.5 -alpha=10000 -ncycl=50 -np=4
    core_decon  measurement.mrc result.mrc otf.otf -ncycl=10

Run either program without arguments for the list of options. `docs/SPECIFICATION.md` describes what each program computes.

**With AR-Decon.** AR-Decon finds `core2_decon` on `PATH`. Install this `core2_decon` ahead of any other copy, e.g. `make install PREFIX=$CONDA_PREFIX` inside the AR-Decon environment. No AR-Decon file needs to change.

## Verification

| Test | Result |
|---|---|
| `make test`: the golden tests delivered with the source (`JacobiNL{2,3}DduTests.c`, Matlab reference) | 6/6 pass |
| `make check`: AR-Decon's configurations (4 cases on a synthetic 64³ map) against outputs of the original Priism 4.7.2 binaries, plus AR-Decon's probe | bit-identical or within 2 ulp; cost traces identical |
| Full regression suite against the original binaries during development: 29 cases (23 ER-Decon II, 6 ER-Decon I), 87 further synthetic runs, 15 runs on public fluorescence data | bit-identical or within 1–2 ulp; traces identical |
| AR-Decon tutorial (EMD-8731, 256³, 50 iterations) | 14 of 16,777,216 voxels differ, each by 1 ulp (rel RMS 4.8e-11). The original binary differs from itself by a similar amount between two ordinary runs. |
| Image I/O (`imlite.c`) against Priism's image library, 36 invocations | output files, stdout and stderr byte-identical |

The references in `tests/ref` were made by the original binaries with deterministic FFT planning. `docs/SPECIFICATION.md` §10 gives the full verification record.

## Not included

ER-SIM-Decon (`core_si_decon`), the 4D solver, the Motif GUIs, and mask-updater plugins (`-umask` stops with "failed to load mask update plugin").

## Citation

If you use this software, please cite:

M. Arigovindan, J. C. Fung, D. Elnatan, V. Mennella, Y.-H. M. Chan, M. Pollard, E. Branlund, J. W. Sedat and D. A. Agard, "High-resolution restoration of 3D structures from widefield images with extreme low signal-to-noise-ratio", *Proc. Natl. Acad. Sci. USA* **110**, 17344–17349 (2013), doi:10.1073/pnas.1315675110.

## Credits

ER-Decon was developed by Muthuvel Arigovindan and colleagues in the laboratories of John W. Sedat and David A. Agard at UCSF. The reconstruction of the missing source files, the Priism-free build and the test suite were done in the Cheng laboratory, UCSF.

## Licence

BSD 3-Clause; see `LICENSE`. `NOTICE.md` records the origin of every file.

## Files

| Path | Contents |
|---|---|
| `src/*.c`, `src/*.h` | the delivered ER-Decon source, unchanged |
| `src/recon/` | reconstructed files, new image I/O, test drivers |
| `src/recon/compat/` | declarations of the Priism functions the front ends call |
| `tests/` | regression cases, synthetic input data, reference outputs, `run_check.py` |
| `contrib/fftw_estimate_shim.c` | makes the original Priism binaries deterministic, for comparisons |
| `docs/SPECIFICATION.md` | numerical specification of ER-Decon I and II |
| `docs/DELIVERED-MD5.txt` | checksums of the delivered source files |
| `recipe/` | conda-build recipe |
| `NOTICE.md` | origin of every file |
| `CHANGELOG.md` | changes |
