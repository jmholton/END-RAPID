#!/usr/bin/env phenix.python
"""
Perturb Fobs by rms(sigF):  F_new = F + N(0,1) * sigF
Drop-in Python replacement for kick_data.com

Usage: phenix.python kick_data.py [refme.mtz] [seed=N] [F=label]
Output: kicked.mtz
"""
from __future__ import print_function
import sys
import os
import time
import numpy as np
import iotbx.mtz
from scitbx.array_family import flex


def _looks_like_fobs(label):
    """True if the column label looks like an observed amplitude (not a derived quantity)."""
    s = label.upper().replace('-', '').replace('_', '')
    return any(s == k for k in ('FOBS', 'FP', 'F', 'FMEAS', 'FOBS1')) or \
           s.startswith('FOBS') or s.startswith('FMEAS')


def _pick_f_sigf(mtz_obj, user_F=None):
    """Return (F_label, SIGF_label) for the best F/SIGF pair in the MTZ.

    Preference order:
      1. user-specified label if given
      2. F-type column whose name looks like observed amplitudes (FOBS, FP, …)
         paired with the immediately following Q-type column
      3. Highest-completeness F/Q pair with the highest mean amplitude
         (discriminates FOBS ~60 e- from K_MASK ~0.007 or K_ISO ~1.0)
    """
    candidates = []
    for crystal in mtz_obj.crystals():
        for dataset in crystal.datasets():
            cols   = {c.label(): c for c in dataset.columns()}
            labels = list(cols.keys())
            for idx, fl in enumerate(labels):
                if cols[fl].type() != 'F':
                    continue
                if idx + 1 >= len(labels) or cols[labels[idx+1]].type() != 'Q':
                    continue
                sl     = labels[idx + 1]
                f_vals = np.array(cols[fl].extract_values())
                comp   = float(np.sum(np.isfinite(f_vals))) / max(len(f_vals), 1)
                mean_f = float(np.nanmean(f_vals)) if f_vals.size else 0.0
                named  = _looks_like_fobs(fl)
                candidates.append((named, comp, mean_f, fl, sl))

    if not candidates:
        return None, None
    if user_F:
        for _, _, _, fl, sl in candidates:
            if fl == user_F:
                return fl, sl
    # Sort: name-match first, then completeness, then mean amplitude
    candidates.sort(reverse=True)
    return candidates[0][3], candidates[0][4]


def kick_data(mtzfile="refme.mtz", seed=None, F_label=None, output="kicked.mtz"):
    if seed is None:
        seed = int(time.time() * 1e6) % (2**31)
    np.random.seed(seed)
    print("using seed = %d" % seed)

    mtz_obj = iotbx.mtz.object(file_name=mtzfile)
    F, SIGF = _pick_f_sigf(mtz_obj, user_F=F_label)
    if F is None:
        print("ERROR: no F/SIGF pair found in %s" % mtzfile)
        sys.exit(9)
    print("selected F=%s SIGF=%s" % (F, SIGF))

    for crystal in mtz_obj.crystals():
        for dataset in crystal.datasets():
            cols = {c.label(): c for c in dataset.columns()}
            if F in cols and SIGF in cols:
                f_vals  = np.array(cols[F].extract_values(),    dtype=np.float64)
                sf_vals = np.array(cols[SIGF].extract_values(), dtype=np.float64)

                noise = np.random.normal(0.0, 1.0, len(f_vals)) * sf_vals
                f_new = np.maximum(f_vals + noise, 0.0)

                missing = np.isnan(f_vals) | np.isnan(sf_vals)
                valid = ~missing
                f_new[missing] = 0.0  # placeholder for missing positions

                cols[F].set_values(
                    flex.float(f_new.astype(np.float32)),
                    flex.bool(valid.tolist())
                )

    mtz_obj.write(file_name=output)
    print("kicked.mtz contains %s from %s modified by rms %s" % (F, mtzfile, SIGF))


if __name__ == "__main__":
    mtzfile = "refme.mtz"
    seed = None
    F_label = None
    for arg in sys.argv[1:]:
        if arg.endswith(".mtz"):
            mtzfile = arg
        elif arg.startswith("seed="):
            seed = int(arg.split("=")[1])
        elif arg.startswith("F="):
            F_label = arg.split("=")[1]
    kick_data(mtzfile=mtzfile, seed=seed, F_label=F_label)
