#!/usr/bin/env phenix.python
"""
Perturb Fobs by rms(Fobs - Fcalc):  F_new = F + N(0,1) * |F - FC|
Drop-in Python replacement for kick_data_bydiff.com

Usage: phenix.python kick_data_bydiff.py [refined.mtz] [seed=N] [F=label] [FC=label]
Output: kicked.mtz
"""
from __future__ import print_function
import sys
import time
import numpy as np
import iotbx.mtz
from scitbx.array_family import flex


def _find_column(cols, user_label, col_type):
    """Return first column label matching type, or user_label if given."""
    if user_label and user_label in cols:
        return user_label
    for label, col in cols.items():
        if col.type() == col_type:
            return label
    return None


def kick_data_bydiff(mtzfile="refined.mtz", seed=None,
                     F_label=None, FC_label=None, output="kicked.mtz"):
    if seed is None:
        seed = int(time.time() * 1e6) % (2**31)
    np.random.seed(seed)
    print("using seed = %d" % seed)

    mtz_obj = iotbx.mtz.object(file_name=mtzfile)

    # Collect all columns across datasets into a flat dict (label -> column)
    all_cols = {}
    for crystal in mtz_obj.crystals():
        for dataset in crystal.datasets():
            for col in dataset.columns():
                all_cols[col.label()] = col

    F    = _find_column(all_cols, F_label,  'F')
    SIGF = None
    FC   = _find_column(all_cols, FC_label, 'F')  # Fcalc also stored as type F

    # Prefer an F that has a paired Q (sigF)
    labels = list(all_cols.keys())
    for i, lbl in enumerate(labels):
        if all_cols[lbl].type() == 'F' and lbl != FC:
            if i + 1 < len(labels) and all_cols[labels[i+1]].type() == 'Q':
                if F_label is None or lbl == F_label:
                    F    = lbl
                    SIGF = labels[i+1]
                    break

    if F is None:
        print("ERROR: no Fobs in %s" % mtzfile)
        sys.exit(9)
    if FC is None:
        print("ERROR: no Fcalc in %s" % mtzfile)
        sys.exit(9)

    print("selected F=%s SIGF=%s" % (F, SIGF or ""))
    print("using FC=%s" % FC)

    f_vals  = np.array(all_cols[F].extract_values(),  dtype=np.float64)
    fc_vals = np.array(all_cols[FC].extract_values(), dtype=np.float64)
    delta   = np.abs(f_vals - fc_vals)

    noise = np.random.normal(0.0, 1.0, len(f_vals)) * delta
    f_new = np.maximum(f_vals + noise, 0.0)

    missing = np.isnan(f_vals) | np.isnan(fc_vals)
    f_new[missing] = np.nan

    all_cols[F].set_values(flex.double(f_new))

    # Propagate missing: set SIGF absent where F is absent
    if SIGF and SIGF in all_cols:
        sf_vals = np.array(all_cols[SIGF].extract_values(), dtype=np.float64)
        sf_vals[missing] = np.nan
        all_cols[SIGF].set_values(flex.double(sf_vals))

    mtz_obj.write(file_name=output)
    print("kicked.mtz contains %s from %s modified by rms ( %s - %s )" % (F, mtzfile, F, FC))


if __name__ == "__main__":
    mtzfile = "refined.mtz"
    seed = None
    F_label = None
    FC_label = None
    for arg in sys.argv[1:]:
        if arg.endswith(".mtz"):
            mtzfile = arg
        elif arg.startswith("seed="):
            seed = int(arg.split("=")[1])
        elif arg.startswith("F="):
            F_label = arg.split("=")[1]
        elif arg.startswith("FC="):
            FC_label = arg.split("=")[1]
    kick_data_bydiff(mtzfile=mtzfile, seed=seed, F_label=F_label, FC_label=FC_label)
