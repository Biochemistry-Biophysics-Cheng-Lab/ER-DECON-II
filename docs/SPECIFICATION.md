# ER-Decon I and II: numerical specification

What the original Priism `core_decon` / `core2_decon` compute, established
either by reading the delivered source or, for the undelivered pieces, by
reconstruction verified against the original binaries. Each statement
carries its evidence tag:

- **[D]** delivered source, file cited; verified end to end against the
  original binaries.
- **[R]** reconstructed; verified by the test cited.
- "Golden tests" = the delivered `JacobiNL{2,3}DduTests.c` with their
  Matlab reference values (`make test`).
- "Original binaries" = comparison of this build with the original Priism
  4.7.2 binaries under deterministic FFT planning (`make check` and the
  wider runs in §10).

Status of the evidence: all six golden tests pass; 87 synthetic runs,
29 regression cases and 15 public-data runs (§10) are bit-identical or
within 1-2 float32 ulp of the original binaries, with identical cost traces.

---

## 0. Conventions

- Volumes are `Nx × Ny × Nz` (x fastest), stored in the FFTW half-complex
  padded layout: rows of `Nx + NxFSkip` doubles, `NxFSkip = 2` if `Nx` even
  else 1; only the first `Nx` entries of a row are data. `N = Nx·Ny·Nz`.
- `F` = unnormalised forward real-to-complex FFT, `F⁻¹` its unnormalised
  inverse; the code folds `1/N` into spectral multiplies. Below, `h * x`
  means the normalised circular convolution `F⁻¹(H · F x)/N`.
- All arithmetic is float64; the output file is float32.
- δ = `delfrac` = dxy/dz (lateral over axial pixel size).
- FFT planning: the delivered code plans with `FFTW_MEASURE` (timing
  dependent, last-bit nondeterministic). A faithful re-implementation should
  use a fixed plan; this build forces `FFTW_ESTIMATE` by default
  (`recon/fftw_deterministic.h`).

## 1. Front end (both tools) [D] `core_decon.c`, `core2_decon.c`

1. **Region** `-x/-y/-z/-w/-t` selected; one (wavelength, time) volume is
   processed at a time.
2. **Background**: `-abs_background -sub=b` subtracts `b`; otherwise the
   background is estimated per wavelength (`core_estimate_background.c`)
   and scaled by the `-sub` factor (default 1). Voxels are `data − bg`.
3. **Z padding** to `Nz_pad` (`-nzpad`, default the smallest power of two
   ≥ 1.3·Nz); the extra sections are filled by `pad_ramp_3d` (a ramp between
   the last and first processed sections, so the volume is smooth under the
   periodic wrap). x and y are not padded.
4. **Rescale** (`rescale()` in `core2_decon.c`): the volume is divided by
   `s` = the intensity above which a fraction `cuth` (default 0.001) of the
   voxels lie. The output is multiplied back by `s`. The per-volume `s` is
   printed in the `scale` column of the log.
5. **OTF**: the on-file cylindrical OTF (`.realctf`: one section per
   wavelength, x = k_z in FFT order with spacing `d[0]`, y = k_r from 0 with
   spacing `d[1]`, in 1/µm; real or complex) is interpolated linearly onto
   the padded 3D grid (`interpolate_cyl_to_rect_rc/_cc`), with frequencies
   scaled by λ_OTF/λ_data, a taper at the last k_r sample and zero beyond.
   Full 3D OTFs are rescaled to H(0) = 1 unless `-nonorm_otf`; radially
   averaged (`.realctf`) ones are never rescaled, so they must be written
   with H(0) = 1.
6. **Output**: the solution × `s`, cropped to the processed z range unless
   `-include_padding`.

## 2. Operators [R] `recon/deriv_kernel.h`; unit test `make unit`; golden tests; original binaries

Periodic (wrap-around) boundaries; with `e_a` the unit step along axis a:

    L_aa X(r) = X(r − e_a) − 2 X(r) + X(r + e_a)                 (a = x, y, z)
    L_ab X(r) = X(r) − X(r − e_a) − X(r − e_b) + X(r − e_a − e_b)   (ab = xy, yz, xz; backward)

`Derivative(X, R, G, dir, mode)`: `'f'` R = L X (overwrite); `'r'`
R += Lᵀ X (**accumulates** — callers zero R first); `'s'` R = (L∘L)ᵀ X, the
squared stencil (1, 4, 1 / 1, 1, 1, 1) applied in the adjoint direction, so
that `(L∘L)ᵀ W = diag(Lᵀ W L)`. Slots `"  "` (identity), `xx yy zz xy yz xz`;
2D: `"  " xx yy xy`.

Evidence: golden costs match to 1e-16 only with periodic + backward
alignment (every other boundary rule or alignment misses by ≥ 1e-5);
adjointness to 1e-13 on nine geometries; ER-Decon I (which uses `'f'` on
full arrays) and ER-Decon II (which uses all three modes) both bit-identical
against the original binaries.

## 3. ER-Decon II setup [D] `core2.c`

    V      = Duchon kernel (AddDuchon3D, misc.c): Σ over the 3-point stencils
             with weights 1, 1, δ⁴ (xx, yy, zz) and 2, 2δ², 2δ² (xy, yz, xz),
             i.e. V̂ = Σ_k w_k |L̂_k|²
    regt   = (1 + laml2) · (lr0 + (1 − lr0) · (lr1 + V̂))
             2014 line search: lr0 = 0, lr1 = lamratio[1] if regtype = log else 0
             2013 line search: lr0 = lamratio[0], lr1 = lamratio[1]
    PCF    = 1 / ( |H|² + lamf · (lampc/lamf · |H_c|² + regt) )    (lampc term only if lampc > 0;
             H_c = missing-cone filter from getmcf.c)
    X0     = h̃ * Imn  with spectrum  conj(H)·Î·PCF   (Tikhonov initial guess; unless -guess)
    P      = sqrt(PCF)                                   (preconditioner)
    B      = hᵀ * Imn  (spectrum conj(H)·Î)
    Diagd  = mean over the full spectrum of |H|²

## 4. ER-Decon II, 2014 line search [R] `recon/solver2014.c`

Per iteration (`JacobiNL3Ddu.c` [D] drives it):

**Cost** (reported and minimised by the line search) (matches the golden-test costs to 1e-16):

    J(X) = Σ_r (h*X − Imn)²
         + (lamf/α) Σ_r [ ρ(α E²(r)) + laml2 · E²(r) ]
         + lamf · lampos  · Σ_r min(X, 0)²
         + lamf · lammask · Σ_r (M · X)²
    ρ(u) = log(1 + u)      regularization_type 1 ("ma", default)
         = sqrt(1 + u)     regularization_type 2 ("tv2"; lamratio is 0 here)
    E²   = lamratio·X² + (L_xx X)² + (L_yy X)² + δ⁴ (L_zz X)² + 4 (L_xy X)² + 4δ² [(L_yz X)² + (L_xz X)²]

**Direction** (`Res3D`):

    E²_g = lamratio·max(X,0)² + (L_xx X)² + (L_yy X)² + δ⁴(L_zz X)² + 2(L_xy X)² + 2δ²[(L_yz X)² + (L_xz X)²]
    W    = ρ'(α E²_g) + laml2            ρ' = 1/(1+u)  or  1/(2 sqrt(1+u))
    R    = hᵀh * X + lamf · [ lamratio·W·X + Σ_k w_k L_kᵀ(W · L_k X) + lampos·min(X,0) + lammask·M·X ]
    Dreg = lamratio·W + Σ_k w_k (L_k∘L_k)ᵀ W                 w = (1, 1, δ⁴, 2, 2δ², 2δ²)
    D    = lamf·Dreg + (max_r Dreg / max_r E²_g) · ( lamf·(lampos·[X<0] + lammask·M) + Diagd )
    U    = P * ( (1/D) · (P * (B − R)) )

**Step** (`LineSearch3D2n`): α_w = argmin_a J(X + aU), by Newton's method on
dJ/da from a = 0 to convergence (J is piecewise smooth in a; every term is
quadratic or ρ of a quadratic per voxel). X ← X + α_w U.

Iterations: fixed `-ncycl` (default 50). The progress line reports J and its
parts (`cost_fid`, `cost_reg` = regulariser incl. laml2, `cost_pos`) after
each step; `rfactor = Σ|h*X − Imn| / Σ|Imn|`.

**The six deliberate inconsistencies** (the update is not the gradient
descent of J; all verified — reproduce them in "reproduce" mode):

| # | cost J | direction R, D | evidence |
|---|---|---|---|
| 1 | mixed weights 4, 4δ² | 2, 2δ² | golden-test costs vs trajectories |
| 2 | laml2 inside (lamf/α)[…] | laml2 added to W without 1/α | golden test 2 |
| 3 | lamratio·X² | lamratio·max(X,0)² in E²_g | golden 3D test 1 (52 % negative voxels) |
| 4 | mask squared, (M X)² | M linear in R and D | fractional mask vs original, bit-identical |
| 5 | positivity from min(X,0)² (Hessian factor 2) | factor 1, and scaled by max Dreg/max E²_g | original binaries: background-subtraction runs, lampos/lamf scans |
| 6 | — | Diagd enters D scaled by max Dreg/max E²_g | plane-wave and two-level initial guesses vs original binaries |

Equivalent reading of the D formula, up to the global scale the exact line
search removes: the regulariser diagonal is rescaled so that its maximum
equals max E²_g, and the positivity, mask and data diagonals are added
unscaled.

2D (`JacobiNL2Ddu.c`, `Res2D`, `LineSearch2D2n`): the same with slots
xx, yy, xy and δ = 1.

## 5. ER-Decon II, 2013 line search [D] `jacobi2.c`, `regc.c`

    Res3D_old:  R = hᵀh * X + lamf · Rt,   D = lamf · Dt + Diagd
    Regsum:     first  Reg, slot "  " only, λ = 1, α₀:
                           W₁ = 1/(1 + α₀ max(X,0)²) + laml2;   R₁ = W₁·X,  D₁ = W₁
                second Reg, slots "  xxyyzzxyyzxz", λ = (lamratio[1], 1, 1, δ⁴, 2, 2δ², 2δ²), α₁:
                           W₂ = 1/(1 + α₁ Σ_k λ_k (L_k X)²) + laml2   (identity slot uses max(X,0))
                           R₂ = Σ_k L_kᵀ(λ_k W₂ L_k X),  D₂ = Σ_k λ_k (L_k∘L_k)ᵀ W₂
                R_t = lr0·R₁ + (1 − lr0)·R₂ + lampos·min(X,0)
                D_t = lr0·D₁ + (1 − lr0)·D₂ + lampos·[X<0]          (lr0 = lamratio[0];
                α₀, α₁ = -alpha values; beta is accepted but unused)
    Direction:  U = P * ((1/D) · (P * (B − R)))
    Step:       JacobiNL3Dstep: ω starts at the previous ω (initially -omega, 0.8);
                X_t = X + ω U; accept if ‖B − R(X_t)‖/‖B‖ decreases, else ω ← 0.7 ω
                (≤ 50 tries). Stops at -tol or -ncycl attempts.

## 6. ER-Decon I [D] `core.c`, `jacobi.c`, `applyreg.c`

Delivered in full; only `Derivative('f')` is ours. Summary (see the source
for exact expressions):

- **Setup** (`core.c`): Poisson-like data weights `Wd = 1 + 10·Ĩ²/max Ĩ²`
  from a smoothed copy Ĩ of the data; Tikhonov guess with
  `PCF = 1/(|H|² + lamf((1−lamratio)V̂ + lamratio) [+ lamf·lampc|H_c|²])`,
  clipped to X ≥ 0; frozen weights `alpham = X0²` and
  `Wm = 1 + 100·tanh(α X0²)/(α X0²)`; per-slot α_k = α / max(1, max_r (L_k X0)²)
  (`RegRMax`, 8 slots: positivity, identity, xx … xz).
- **Iteration** (`jacobi.c`, `ApplyReg`): data term weighted by `Wd`;
  regulariser slots applied as FFT-domain filters built from
  `F(Derivative(PCF_kernel))`: slot 0 positivity (`W = [X<0]`), slot 1
  identity (`W = 1/(1 + α₁ X²)`), slots 2-7 derivatives
  (`W = Wm/(1 + α_k·alpham·(L_k X)²)`); relaxation step ω = 0.8, ω ← 0.7 ω on
  a rejected step. **`-ncycl` counts attempts, including rejected ones.**
- ER-Decon I is the path most sensitive to FFT rounding: its accept/reject
  decisions are thresholds on residuals.

## 7. What is not covered

SIM (`core_si_decon`; `sim_xye` never delivered), the 4D solver (`Res4D`,
`LineSearch4D2n`, `Decon4Ddu.h` not delivered), mask-updater plugins
(`-umask`; stubbed in this build), the Motif GUIs.

## 8. Acceptance criteria for a re-implementation

Against this build (or the original under deterministic FFT planning),
in float64: initial guess, cone filter, J(X0) components ≤ 1e-12 relative;
cost traces ≤ 1e-9 for ≤ 10 iterations, ≤ 1e-7 at 50; images rel RMS
≤ 1e-7 (≤ 1e-6 at 50 iterations), same negative-voxel sets; ER-Decon I runs
with near-tie accept/reject decisions compared up to the branch point.

## 9. Where each piece lives

| Piece | File | Status |
|---|---|---|
| operators | `src/recon/deriv_kernel.h`, `derivatives{,2d}.c` | [R] |
| 2014 R, D, J, line search | `src/recon/solver2014.{c,h}`, `Res{2,3}D.c`, `LineSearch{2,3}D2n.c` | [R] |
| everything else | delivered `src/*.c` | [D] |

The reconstructed solver was also checked against independent NumPy models
(dense, on the golden-test grids: equal to 1e-15; FFT-based, full scale:
equal to 8e-16).

## 10. Verification record

| Suite | Result |
|---|---|
| `make unit` | adjointness / accumulation / `'s'` / in-place, 9 geometries |
| `make test` (the delivered golden tests, Matlab reference) | 6/6; 3D test 1 residual vs Matlab 8.378e-10, equal to the original C's own |
| Regression suite against the original binaries (29 cases, 23 ER-Decon II + 6 ER-Decon I: synthetic crops, public *C. elegans* data with 3 wavelengths, AR-Decon's configuration on a 64³ map; plus AR-Decon's probe) | 27 bit-identical, 2 within 2 ulp in at most 2e-5 of the voxels; traces identical. Same with Ubuntu's FFTW 3.3.10 and conda-forge's FFTW 3.3.11. The public repository ships the 4 AR-Decon cases (synthetic data) and the probe as `make check`. |
| Synthetic runs against the original binaries (87 runs) | 72 bit-identical, the rest within 1 ulp; traces identical |
| EPFL *C. elegans* embryo (3 wavelengths, 256²×64 crop, 15 runs) | 7 bit-identical, 8 within 1 ulp in at most 1e-6 of the voxels; traces identical (one plot-file value differs in the 7th digit) |
| AR-Decon tutorial, EMD-8731 256³ (padded 512), 50 iterations | Against the Priism binary with deterministic planning: 1 ulp in 6.0e-7 of the voxels (rel RMS 4.8e-11), trace identical. The Priism binary run twice as users run it differs from itself by rel RMS 6.4e-12. Time 2:08.5 vs 1:58–2:03, 4 threads |
| AR-Decon pipeline, 0.9 (Priism binary) vs 1.0 (this engine), seeded dFSC | all intermediates identical; final map 9–12 voxels at 1 ulp |
