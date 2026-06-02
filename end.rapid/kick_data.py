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


def _pick_f_sigf(mtz_obj, user_F=None):
    """Return (F_label, SIGF_label) for the best F/SIGF pair in the MTZ."""
    candidates = []
    for crystal in mtz_obj.crystals():
        for dataset in crystal.datasets():
            cols = {c.label(): c for c in dataset.columns()}
            f_cols = [l for l, c in cols.items() if c.type() == 'F']
            for fl in f_cols:
                # Look for paired Q (sigma) column immediately after
                labels = list(cols.keys())
                try:
                    idx = labels.index(fl)
                    if idx + 1 < len(labels) and cols[labels[idx+1]].type() == 'Q':
                        sl = labels[idx + 1]
                        f_vals = np.array(cols[fl].extract_values())
                        comp = np.sum(~np.isnan(f_vals)) / max(len(f_vals), 1)
                        reso = dataset.wavelength()  # crude proxy; completeness wins
                        candidates.append((comp, fl, sl))
                except ValueError:
                    pass
    if not candidates:
        return None, None
    if user_F:
        for _, fl, sl in candidates:
            if fl == user_F:
                return fl, sl
    candidates.sort(reverse=True)  # highest completeness first
    return candidates[0][1], candidates[0][2]


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
