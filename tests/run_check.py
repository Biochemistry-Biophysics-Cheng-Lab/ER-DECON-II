#!/usr/bin/env python3
"""Regression suite (`make check`).

Runs every case in cases.txt (next to this script) with core_decon /
core2_decon and compares with references produced by the ORIGINAL Priism
4.7.2 binaries (ref/, written with --make-refs under the FFTW_ESTIMATE shim):

  * result.mrc: BIT-identical, or ULP (at most MAX_ULP float32 ulps in at
    most MAX_FRAC of the voxels), else FAIL;
  * mcf.mrc (the cone filter, when written): by absolute error, see
    compare_mcf();
  * the cost / residual trace (numeric lines of stdout and plot.plt): must be
    identical as printed, else FAIL.

  run_check.py                    check build/bin/core*_decon
  run_check.py --make-refs BINDIR write ref/ from BINDIR/core*_decon
  run_check.py --bin-dir DIR      check DIR/core2_decon and DIR/core_decon
  run_check.py -k NAME            only cases whose name contains NAME

Needs Python 3 and NumPy.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
REF = os.path.join(HERE, "ref")
# exactly the options ardecon passes to core2_decon (scripts/ardecon), minus -np/-ncycl
COMMON_SPA = ["-np=1", "-alpha=10000", "-lamratio=0:1", "-lamf=0.5", "-lampc=0", "-lampos=1",
              "-lamsmooth=100", "-laml2=0", "-cuth=0.001", "-na=1.4", "-nimm=1.512", "-omega=0.8",
              "-sub=1:1:1:1:1", "-tol=0.0001", "-guess=", "-linesearch=2014", "-regtype=ma",
              "-oplotfile=plot.plt"]
# AR-Decon's param_search passes only these (plus -lamf -alpha -lampc=0 -ncycl in the case)
COMMON_PSEARCH = ["-np=1", "-oplotfile=plot.plt"]
USAGE = b"usage: core2_decon measurement result CTF"
MAX_ULP = 2
MAX_FRAC = 1e-3
OUTPUTS = ("result.mrc", "mcf.mrc")


def read_cases():
    cases = []
    for line in open(os.path.join(HERE, "cases.txt")):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        name, tool, opts = [f.strip() for f in line.split("|")]
        cases.append((name, tool, opts))
    return cases


def mrc_data(path):
    """Return the voxel data of an MRC file as a raw uint32 array (bit pattern)."""
    raw = open(path, "rb").read()
    for bo in ("<", ">"):
        nx, ny, nz, mode = np.frombuffer(raw[:16], dtype=bo + "i4")
        if 0 < nx < 1 << 16 and 0 < ny < 1 << 16 and 0 < nz < 1 << 16 and mode in (2, 4):
            nxt = int(np.frombuffer(raw[92:96], dtype=bo + "i4")[0])
            n = int(nx) * int(ny) * int(nz) * (2 if mode == 4 else 1)
            a = np.frombuffer(raw[1024 + nxt:1024 + nxt + 4 * n], dtype=bo + "u4")
            return a.astype("<u4")
    raise ValueError("%s: not a float32/complex64 MRC file" % path)


def compare_image(a_path, b_path):
    if os.path.basename(a_path) == "mcf.mrc":
        return compare_mcf(a_path, b_path)
    return compare_result(a_path, b_path)


def compare_mcf(a_path, b_path):
    """The missing-cone filter (-omcf) by absolute error, <= 1e-12 of its
    largest value.  Its imaginary parts are zero in exact arithmetic and hold
    only FFT rounding (~1e-17), which differs between FFTW builds (Ubuntu's
    vs conda-forge's), so a float32-ulp rule does not apply to this file.
    Results (result.mrc) keep the ulp rule."""
    a, b = mrc_data(a_path), mrc_data(b_path)
    if a.shape != b.shape:
        return "FAIL", "shape %s vs %s" % (a.shape, b.shape)
    if np.array_equal(a, b):
        return "BIT", ""
    fa, fb = a.view(np.float32).astype(float), b.view(np.float32).astype(float)
    err = float(np.abs(fa - fb).max() / max(np.abs(fb).max(), 1e-300))
    info = "max |diff| %.1e of max (FFT rounding)" % err
    return ("ULP" if err <= 1e-12 else "FAIL"), info


def compare_result(a_path, b_path):
    a, b = mrc_data(a_path), mrc_data(b_path)
    if a.shape != b.shape:
        return "FAIL", "shape %s vs %s" % (a.shape, b.shape)
    if np.array_equal(a, b):
        return "BIT", ""
    ia = a.view(np.int32).astype(np.int64)
    ib = b.view(np.int32).astype(np.int64)
    # map sign-magnitude float bits onto a monotone integer line
    ia = np.where(ia < 0, -(ia & 0x7FFFFFFF), ia)
    ib = np.where(ib < 0, -(ib & 0x7FFFFFFF), ib)
    d = np.abs(ia - ib)
    frac = float(np.mean(d > 0))
    fa, fb = a.view(np.float32).astype(float), b.view(np.float32).astype(float)
    rel = float(np.linalg.norm(fa - fb) / max(np.linalg.norm(fb), 1e-300))
    info = "max %d ulp in %.2e of voxels, rel RMS %.2e" % (int(d.max()), frac, rel)
    if d.max() <= MAX_ULP and frac <= MAX_FRAC:
        return "ULP", info
    return "FAIL", info


NUMLINE = re.compile(r"^\s*-?\d+(\s+[-+]?\d*\.?\d+(e[-+]?\d+)?)+\s*$", re.I)


def trace(rundir):
    """Numeric lines of stdout.log and plot.plt: the cost / residual trace."""
    out = []
    for f in ("stdout.log", "plot.plt"):
        p = os.path.join(rundir, f)
        if os.path.exists(p):
            out += [l.rstrip() for l in open(p, errors="replace") if NUMLINE.match(l)]
    return out


def run_case(bindir, name, tool, opts, workdir):
    os.makedirs(workdir, exist_ok=True)
    env = dict(os.environ)
    nslots = "1"
    if "@NP4@" in opts:
        opts = opts.replace("@NP4@", "")
        nslots = "4"
    env["NSLOTS"] = nslots
    if "@PSEARCH@" in opts:
        opts = opts.replace("@PSEARCH@", "")
        inp, ctf, common = "spa_map64.mrc", "spa_TF64.otf", COMMON_PSEARCH
    elif "@SPA@" in opts:
        # AR-Decon configuration: synthetic 64^3 map, full 3D complex OTF, wavelength 0
        opts = opts.replace("@SPA@", "")
        inp, ctf, common = "spa_map64.mrc", "spa_TF64.otf", COMMON_SPA
    else:
        raise ValueError("case %s: no input data (use @SPA@ or @PSEARCH@)" % name)
    args = [os.path.join(bindir, tool), os.path.join(DATA, inp), "result.mrc",
            os.path.join(DATA, ctf)] + common + opts.replace("@DATA@", DATA).split()
    with open(os.path.join(workdir, "stdout.log"), "w") as fo:
        rc = subprocess.run(args, cwd=workdir, env=env, stdout=fo, stderr=subprocess.STDOUT).returncode
    return rc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--make-refs", metavar="BINDIR")
    ap.add_argument("--bin-dir", default=None)
    ap.add_argument("-k", default="")
    ap.add_argument("--keep", action="store_true", help="keep the run directories")
    opt = ap.parse_args()

    cases = [c for c in read_cases() if opt.k in c[0]]

    if opt.make_refs:
        for name, tool, opts in cases:
            d = os.path.join(REF, name)
            if os.path.isdir(d):
                shutil.rmtree(d)
            rc = run_case(os.path.abspath(opt.make_refs), name, tool, opts, d)
            if rc != 0:
                print("REF FAILED: %s (rc %d)" % (name, rc))
                return 1
            open(os.path.join(d, "trace.txt"), "w").write("\n".join(trace(d)) + "\n")
            for f in ("stdout.log", "plot.plt"):
                if os.path.exists(os.path.join(d, f)):
                    os.remove(os.path.join(d, f))
            print("ref  %-18s %s" % (name, " ".join(sorted(os.listdir(d)))))
        return 0

    bindir = os.path.abspath(opt.bin_dir or os.path.join(HERE, "..", "build", "bin"))
    tools = {"core_decon": "core_decon", "core2_decon": "core2_decon"}
    tmp = tempfile.mkdtemp(prefix="erdecon_check_")
    t0 = time.time()
    nfail = 0
    counts = {}
    for name, tool, opts in cases:
        exe = tools[tool]
        work = os.path.join(tmp, name)
        rc = run_case(bindir, name, exe, opts, work)
        ref = os.path.join(REF, name)
        verdicts = []
        if rc != 0:
            verdicts.append(("FAIL", "exit status %d" % rc))
        else:
            for f in OUTPUTS:
                if os.path.exists(os.path.join(ref, f)):
                    if not os.path.exists(os.path.join(work, f)):
                        verdicts.append(("FAIL", "%s not written" % f))
                    else:
                        v, info = compare_image(os.path.join(work, f), os.path.join(ref, f))
                        verdicts.append((v, "%s %s" % (f, info)))
            want = open(os.path.join(ref, "trace.txt")).read().split("\n")[:-1]
            got = trace(work)
            if got != want:
                ndiff = sum(1 for a, b in zip(got, want) if a != b) + abs(len(got) - len(want))
                verdicts.append(("FAIL", "trace differs in %d of %d lines" % (ndiff, len(want))))
        worst = "FAIL" if any(v == "FAIL" for v, _ in verdicts) else (
            "ULP" if any(v == "ULP" for v, _ in verdicts) else "BIT")
        counts[worst] = counts.get(worst, 0) + 1
        nfail += worst == "FAIL"
        detail = "; ".join(i for v, i in verdicts if v != "BIT" and i)
        print("%-4s %-18s %s" % ("PASS" if worst != "FAIL" else "FAIL", name,
                                 ("(%s) " % worst.lower()) + detail if worst != "BIT" else "(bit-identical)"))
    # the probe AR-Decon's scripts use to find a working core2_decon: run it
    # without arguments; the usage line must appear (on stderr) and the exit
    # status must be non-zero
    if not opt.k:
        p = subprocess.run([os.path.join(bindir, tools["core2_decon"])], capture_output=True)
        ok = p.returncode != 0 and USAGE in p.stderr
        counts["BIT" if ok else "FAIL"] = counts.get("BIT" if ok else "FAIL", 0) + 1
        nfail += not ok
        cases.append(("probe", "core2_decon", ""))
        print("%-4s %-18s %s" % ("PASS" if ok else "FAIL", "probe",
                                 "(usage on stderr, exit %d)" % p.returncode))
    if not opt.keep:
        shutil.rmtree(tmp)
    else:
        print("run directories kept in", tmp)
    print("check: %d cases, %s, %.1f s -> %s" % (
        len(cases), ", ".join("%d %s" % (n, k.lower()) for k, n in sorted(counts.items())),
        time.time() - t0, "ALL PASSED" if nfail == 0 else "%d FAILED" % nfail))
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
