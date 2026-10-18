# NOTICE: provenance of every file

Every file in this package belongs to exactly one of the groups below.

## 1. Delivered ER-Decon source (UCSF), unchanged

**Files:** `src/*.c`, `src/*.h` (46 files).

- **ER-Decon II:** `core2_decon.c`, `core2.c/h`, `core2_2d.c/h`, `core2_params.h`, `jacobi2*.c/h`, `regc*.c/h`, `JacobiNL3Ddu.c/h`, `JacobiNL2Ddu.c/h`.
- **ER-Decon I:** `core_decon.c`, `core.c/h`, `jacobi.c/h`, `applyreg.c/h`.
- **Shared:** `misc*.c/h`, `fft_double.c/h`, `getmcf.c/h`, `interp_util.c/h`, `binary_double_heap.c`, `binary_heap.h`, `core_estimate_background.c/h`, `res_reporter.h`.
- **Golden tests:** `JacobiNL3DduTests.c/h`, `JacobiNL2DduTests.c/h`.

These are the ER-Decon programs and solvers as delivered by the Sedat/Agard laboratory (UCSF) in August–September 2026. Each file is byte-identical to the delivered copy (checked by MD5 against the original deliveries on 2026-09-24; the checksums are in `docs/DELIVERED-MD5.txt`).

The file headers carry the notices they were delivered with: "Copyright 2009 Muthuvel Arigovindan" in 2 files (`regc.c`, `regc2d.c`), the Xcode placeholder "Copyright 2010/2011 \_\_MyCompanyName\_\_" in 20 files, and no notice in the other 24. They were left as delivered, so that the checksums still match; the copyright and licence of the whole repository are stated once, in `LICENSE`. `docs/DELIVERED-MD5.txt` lists the MD5 checksum of every delivered file (check with `cd src && md5sum -c ../docs/DELIVERED-MD5.txt`).

The method is described in M. Arigovindan, J. C. Fung, D. Elnatan, V. Mennella, Y.-H. M. Chan, M. Pollard, E. Branlund, J. W. Sedat and D. A. Agard, *PNAS* 110, 17344 (2013). The 2014 line search is a later refinement.

## 2. Reconstructed files (Cheng laboratory, UCSF, 2026)

**Files:** `src/recon/derivatives.c/h`, `derivatives2d.c/h`, `deriv_kernel.h`, `volumestruct.c/h`, `volumestruct2d.c/h`, `solver2014.c/h`, `Res3D.c/h`, `Res2D.c/h`, `LineSearch3D2n.c/h`, `LineSearch2D2n.c/h`, `dtests.c/h`.

The delivered files call these, but they were not part of the delivery:

- the finite-difference operators (`Derivative`, `Derivative2D`) and the volume geometry (`ConstructVolumeStruct*`);
- the 2014 residual/diagonal step (`Res3D`, `Res2D`) and exact line search (`LineSearch3D2n`, `LineSearch2D2n`);
- the golden-test comparison helper (`dtests`).

They were written from the PNAS paper and its SI, the calling code, the delivered golden tests and black-box comparison with the original binaries. No original code for them was available. Their function signatures follow the calls in the delivered code.

## 3. New code (Cheng laboratory, UCSF, 2026)

- `src/recon/imlite.c` and `src/recon/compat/*.h`: the MRC/Priism image-file functions the two front ends call, replacing Priism's `libimlib`/`libive`.
  - `compat/*.h` declare the subset of Priism's IM/IVE API that the delivered front ends include, with the same names and constant values, so that the front ends compile unchanged.
  - `compat/masku_plugin.h` declares the mask-plugin interface types that `core2_decon.c` refers to.
  - `imlite.c` reproduces libimlib's header-summary text, so the logs are unchanged.
- `src/recon/masku_stubs.c`: mask-updater plugins are not supported.
- `src/recon/fftw_deterministic.h`: forces `FFTW_ESTIMATE` planning (build option `FFTW_DETERMINISTIC`).
- `src/recon/run_golden_tests.c`, `src/recon/test_derivatives.c`: test drivers.
- `contrib/fftw_estimate_shim.c`: an `LD_PRELOAD` shim that forces deterministic FFT planning in the original Priism binaries, for like-for-like comparisons.
- `Makefile`, `tests/run_check.py`, `tests/cases.txt`, `recipe/`, `.github/`, `.gitignore`, and the documentation (`README.md`, `CHANGELOG.md`, `NOTICE.md`, `CITATION.cff`, `docs/SPECIFICATION.md`, `docs/DELIVERED-MD5.txt`).

## 4. Third-party libraries (linked, not included in the source)

- **FFTW 3** (GPL-2.0-or-later). The source in this repository does not include FFTW. A statically linked binary contains it, and the GPL then applies to that binary when it is redistributed.
- A static binary also contains the C runtime and the GNU OpenMP runtime (`libgomp`, GPL-3 with the GCC Runtime Library Exception).

## 5. Test data (`tests/data`, `tests/ref`)

| File | Origin |
|---|---|
| `spa_map64.mrc`, `spa_TF64.otf` | synthetic, Cheng laboratory |
| `tests/ref/*` | outputs of the original Priism 4.7.2 binaries (deterministic FFT planning) on the data above |
