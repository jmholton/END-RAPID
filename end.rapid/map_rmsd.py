#!/usr/bin/env phenix.python
"""
Compute per-voxel RMS deviation of a set of CCP4 maps from a reference.
Drop-in Python replacement for map_rmsd.com

Usage: phenix.python map_rmsd.py *.map [ref=avg.map]
Output: sigma.map  (and avg.map if no ref given)
"""
from __future__ import print_function
import sys
import os
import numpy as np
from iotbx.map_manager import map_manager
from scitbx.array_family import flex


def map_rmsd(map_files, ref_file=None, verbose=True):
    if not map_files:
        print("usage: %s *.map [ref=map]" % sys.argv[0])
        sys.exit(9)

    # Load reference template (first map) for metadata
    ref_mm = map_manager(file_name=map_files[0])

    def load_data(path):
        return np.array(map_manager(file_name=path).map_data(), dtype=np.float64)

    def write_map(data, filename, template_mm):
        out = template_mm.customized_copy(map_data=flex.double(data.flatten()))
        out.write_map(file_name=filename)

    # Compute or load reference
    if ref_file and os.path.isfile(ref_file):
        ref = load_data(ref_file)
        print("using reference: %s" % ref_file)
    else:
        print("computing average map from %d inputs" % len(map_files))
        acc = np.zeros_like(load_data(map_files[0]))
        for f in map_files:
            print(f)
            acc += load_data(f)
        ref = acc / len(map_files)
        write_map(ref, "avg.map", ref_mm)
        print("avg.map is the average of all maps")
        ref_file = "avg.map"

    # Accumulate sum of squared deviations
    n = len(map_files)
    sumsq = np.zeros_like(ref)
    for f in map_files:
        diff = load_data(f) - ref
        print("%s - %s" % (f, ref_file))
        sumsq += diff * diff

    scale = 1.0 / (n - 1) if n > 1 else 1.0
    if n == 1:
        print("WARNING: only one map! output will simply be absolute difference.")
    variance = sumsq * scale
    sigma = np.sqrt(np.maximum(variance, 0.0))

    write_map(sigma, "sigma.map", ref_mm)
    print("sigma.map is the rms deviation from %s" % ref_file)

    # Print map statistics
    mm_sigma = map_manager(file_name="sigma.map")
    d = np.array(mm_sigma.map_data())
    print("sigma.map :")
    print("  Minimum density ........ %12g" % d.min())
    print("  Maximum density ........ %12g" % d.max())
    print("  Mean density ........... %12g" % d.mean())
    print("  Rms deviation from mean  %12g" % d.std())


if __name__ == "__main__":
    maps = []
    ref = None
    for arg in sys.argv[1:]:
        if arg.startswith("ref="):
            ref = arg[4:]
        else:
            maps.append(arg)
    map_rmsd(maps, ref_file=ref)
