#!/usr/bin/env phenix.python
"""
Find the "vacuum level" of a CCP4 electron density map using robust statistics.
Drop-in Python replacement for map_vacuum_level.com

Usage: phenix.python map_vacuum_level.py map.ccp4
"""
from __future__ import print_function
import sys
import numpy as np
from iotbx.map_manager import map_manager
from scitbx.array_family import flex


def find_vacuum_level(mapfile, verbose=True):
    """Return (vacuum, mad) for the given CCP4/MRC map."""
    mm = map_manager(file_name=mapfile)
    data = np.array(mm.map_data(), dtype=np.float64)
    n_voxels = data.size
    mean = data.mean()
    cell_volume = mm.crystal_symmetry().unit_cell().volume()

    if verbose:
        print("cell volume is : %.3f    mean electron density: %g electrons/A^3"
              % (cell_volume, mean))

    # Work only on the negative side: start with values below mean
    cut = mean if mean > 0 else 1e-99
    neg = np.sort(data[data < cut])
    median = cut
    mad = 1e-6

    # --- symmetric median rejection ---
    last_median = 1e99
    while median != last_median and len(neg) >= n_voxels * 0.01:
        mn = neg[0]
        med = float(np.median(neg))
        new_cut = mn + 2.0 * (med - mn)
        new_neg = neg[neg < new_cut]
        if len(new_neg) < n_voxels * 0.01:
            break
        frac = int(len(new_neg) / n_voxels * 100)
        if verbose:
            print("%s median= %g   cut= %g   frac= %d%%   min= %g"
                  % (mapfile, mn, new_cut, frac, mn))
        last_median = median
        median = mn
        cut = new_cut
        neg = new_neg

    # --- upper MAD rejection (+sig * MAD cutoff) ---
    for sig in [4, 3]:
        last_median = 1e99
        while median != last_median and len(neg) >= n_voxels * 0.01:
            mad = max(float(np.median(np.abs(neg - median))), 1e-6)
            new_neg = neg[neg <= median + sig * mad]
            if len(new_neg) < n_voxels * 0.01:
                break
            frac = int(len(new_neg) / n_voxels * 100)
            last_median = median
            median = float(np.median(new_neg))
            if verbose:
                print("%s median= %g   mad= %g   frac= %d%%   min= %g"
                      % (mapfile, median, mad, frac, new_neg[0]))
            neg = new_neg

    # --- symmetric MAD rejection (±sig * MAD cutoff) ---
    for sig in [4, 3]:
        last_median = 1e99
        while median != last_median and len(neg) >= n_voxels * 0.01:
            mad = max(float(np.median(np.abs(neg - median))), 1e-6)
            new_neg = neg[np.abs(neg - median) <= sig * mad]
            if len(new_neg) < n_voxels * 0.01:
                break
            frac = int(len(new_neg) / n_voxels * 100)
            last_median = median
            median = float(np.median(new_neg))
            if verbose:
                print("%s median= %g   mad= %g   frac= %d%%   min= %g"
                      % (mapfile, median, mad, frac, new_neg[0]))
            neg = new_neg

    vacuum = float(neg.mean()) if len(neg) > 0 else mean
    frac = int(len(neg) / n_voxels * 100)
    print("%s vacuum level: %g +/- %g  occupies %d%% of map" % (mapfile, vacuum, mad, frac))

    F000 = cell_volume * (mean - vacuum)
    print("estimated F000 = %.1f" % F000)

    # Write offset map (vacuum = 0)
    offset_data = flex.double((data - vacuum).flatten())
    mm_out = mm.customized_copy(map_data=offset_data)
    mm_out.write_map(file_name="vacuum_zero.map")
    print("vacuum_zero.map has a vacuum level of zero")

    return vacuum, mad


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: %s map.ccp4" % sys.argv[0])
        sys.exit(9)
    find_vacuum_level(sys.argv[1])
