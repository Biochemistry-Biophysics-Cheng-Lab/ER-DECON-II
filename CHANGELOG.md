# Changelog

## 1.0.0

First public release.

- `core2_decon` (ER-Decon II, with the 2014 line search) and `core_decon` (ER-Decon I), built from the original source without Priism.
- Results match the original Priism 4.7.2 binaries to the last bit, or within 1–2 float32 ulp in a few voxels per million.
- FFT planning is deterministic by default (`FFTW_DETERMINISTIC=1`).
