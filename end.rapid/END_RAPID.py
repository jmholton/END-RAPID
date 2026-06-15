#!/usr/bin/env phenix.python
"""
END/RAPID: Electron Number Density map and Refinement Against Perturbed Input Data.

Python replacement for END_RAPID.com — uses phenix.python/cctbx for all
map and MTZ operations; calls phenix.refine and phenix.fmodel via subprocess.

Usage:
  phenix.python END_RAPID.py phenixrefine.eff [seeds=5] [cycles=5] [cpus=N]
                              [-nofofc] [-nosigf] [-norapid]

                                                           -James Holton 4-14-15
  updated for Phenix 2.x / cctbx Python rewrite            James Holton 6-1-26
"""
from __future__ import print_function
import sys, os, re, math, subprocess, multiprocessing, shutil, glob, time
import numpy as np
import iotbx.mtz, iotbx.pdb
from iotbx.map_manager import map_manager
from scitbx.array_family import flex
from cctbx import crystal, sgtbx, uctbx

# ── helpers ──────────────────────────────────────────────────────────────────

def run(cmd, log=None, check=True):
    """Run a shell command, optionally writing stdout+stderr to a log file."""
    with open(log, 'w') if log else open(os.devnull, 'w') as fh:
        ret = subprocess.call(cmd, shell=True, stdout=fh, stderr=subprocess.STDOUT)
    if check and ret:
        raise RuntimeError("Command failed (exit %d): %s" % (ret, cmd))
    return ret


def mtz_columns(mtzfile):
    """Return {label: type} for every column in an MTZ."""
    obj = iotbx.mtz.object(file_name=mtzfile)
    result = {}
    for c in obj.crystals():
        for d in c.datasets():
            for col in d.columns():
                result[col.label()] = col.type()
    return result


def mtz_arrays(mtzfile):
    return iotbx.mtz.object(file_name=mtzfile).as_miller_arrays()


def find_F_SIGF(mtzfile):
    """Return (F_label, SIGF_label) for the best amplitude pair."""
    obj = iotbx.mtz.object(file_name=mtzfile)
    for arr in obj.as_miller_arrays():
        if arr.is_xray_amplitude_array() and arr.sigmas() is not None:
            labels = arr.info().label_string().split(',')
            return labels[0], labels[1]
    return None, None


def fft_map_to_file(mtzfile, F_label, PHI_label, outfile, grid=None, scale=1.0):
    """FFT a (F, PHI) pair from an MTZ into a CCP4 map (volume-scaled, ASU)."""
    obj = iotbx.mtz.object(file_name=mtzfile)
    arrays = obj.as_miller_arrays()

    # Find the complex (amplitude+phase) pair
    target = None
    for arr in arrays:
        ls = arr.info().label_string()
        if F_label in ls and PHI_label in ls:
            target = arr
            break

    if target is None:
        # amplitude and phase stored in separate arrays — combine
        f_arr = phi_arr = None
        for arr in arrays:
            ls = arr.info().label_string()
            if F_label in ls and arr.is_real_array():
                f_arr = arr
            if PHI_label in ls and arr.is_real_array():
                phi_arr = arr
        if f_arr is None or phi_arr is None:
            raise RuntimeError("Cannot find %s/%s in %s" % (F_label, PHI_label, mtzfile))
        target = f_arr.phase_transfer(phi_arr, deg=True)

    if scale != 1.0:
        target = target.customized_copy(data=target.data() * scale)

    fmap = target.fft_map(resolution_factor=0.25)
    fmap.apply_volume_scaling()

    # real_map_unpadded() gives the full FFT grid; wrap into map_manager
    # so downstream code can use map_manager I/O and customized_copy.
    rmap = fmap.real_map_unpadded()
    mm = map_manager(
        map_data=rmap,
        unit_cell_crystal_symmetry=target.crystal_symmetry(),
        unit_cell_grid=rmap.focus(),
        wrapping=True,
    )
    mm.write_map(file_name=outfile)
    return mm


def _write_map(mm_template, data_np, outfile):
    """Write a numpy array as a CCP4 map using mm_template for metadata."""
    shape = mm_template.map_data().all()
    fd = flex.double(data_np.flatten())
    fd.reshape(flex.grid(shape))
    mm_template.customized_copy(map_data=fd).write_map(file_name=outfile)


def map_add_offset(infile, offset, outfile):
    """Write infile + offset as outfile."""
    mm = map_manager(file_name=infile)
    _write_map(mm, np.array(mm.map_data()) + offset, outfile)


def map_scale(infile, scale, outfile):
    mm = map_manager(file_name=infile)
    _write_map(mm, np.array(mm.map_data()) * scale, outfile)


def map_ratio(num_file, denom_file, outfile):
    """Write per-voxel num/denom, skipping zero-denominator voxels."""
    mm_n = map_manager(file_name=num_file)
    n = np.array(mm_n.map_data())
    d = np.array(map_manager(file_name=denom_file).map_data())
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(d != 0, n / d, 0.0)
    _write_map(mm_n, ratio, outfile)


def map_stats(mapfile):
    d = np.array(map_manager(file_name=mapfile).map_data())
    print("           Minimum density .................................  %12g" % d.min())
    print("           Maximum density .................................  %12g" % d.max())
    print("           Mean density ....................................  %12g" % d.mean())
    print("           Rms deviation from mean density .................  %12g" % d.std())


def scale_mtz(ref_mtz, ref_label, target_mtz, target_label,
              mode='fit_scaleB_use_scale'):
    """Return scale factor(s) so that scale × target ≈ ref.

    Mirrors CCP4 scaleit refine isotropic:
      - Initial k = sqrt(sum(ref²)/sum(tgt²))          [scaleit formula]
      - Residuals: D1=0.5*(ref²-S2²*tgt²), D2=0.5*(tgt²-ref²/S2²)  [NOWT=1]
      - S2 = k * exp(-B * sin²θ/λ²)                    [Debye-Waller]
      - Printed B_scaleit = 4 * B_fit (CCP4 SSQ=sin²θ/λ², so BETA*SSQ*0.25)

    mode:
      'fit_scale_use_scale'   : sqrt(sum(ref²)/sum(tgt²)) only, return scalar k
      'fit_scaleB_use_scale'  : fit k+B (like scaleit refine isotropic),
                                discard B, return scalar k  [DEFAULT]
      'fit_scaleB_use_scaleB' : fit k+B, return per-reflection numpy array
                                k × exp(−B × sin²θ/λ²) for all tgt reflections

    Bugs fixed vs old mean-ratio implementation:
      1. common_sets() instead of separate common_set() calls — the latter
         returns ref and tgt in different index orders (0/16686 matched in test).
      2. sqrt(sum(FP²)/sum(FPH²)) matches scaleit's initial k (not sum-ratio).
      3. Intensity residuals D1/D2 match scaleit NOWT=1 (not amplitude LS).
    """
    def _find(mtzfile, label):
        for arr in mtz_arrays(mtzfile):
            ls = arr.info().label_string()
            if label in ls:
                return arr.amplitudes() if not arr.is_real_array() else arr
        return None

    ref_arr = _find(ref_mtz,    ref_label)
    tgt_arr = _find(target_mtz, target_label)
    if ref_arr is None or tgt_arr is None:
        return 1.0
    # common_sets() returns both arrays with the SAME index order (aligned)
    common_r, common_t = ref_arr.common_sets(tgt_arr)
    if common_r.size() == 0:
        return 1.0
    r  = common_r.data().as_numpy_array()
    t  = common_t.data().as_numpy_array()

    # Initial k: sqrt(sum(ref²)/sum(tgt²)) — scaleit's initial scale formula
    k0 = float(np.sqrt(np.sum(r**2) / max(np.sum(t**2), 1e-12)))

    if mode == 'fit_scale_use_scale':
        return k0

    # ── isotropic k+B fit (scaleit refine isotropic, NOWT=1 unweighted) ──────
    from scipy.optimize import least_squares
    cs     = ref_arr.crystal_symmetry()
    M_frac = np.array(cs.unit_cell().fractionalization_matrix()).reshape(3, 3)
    hkl    = np.array(common_t.indices())
    h_orth = hkl.astype(float) @ M_frac
    s2     = np.einsum('ni,ni->n', h_orth, h_orth) / 4.0   # sin²θ/λ² = 1/(4d²)

    valid = (r > 0.001) & (t > 0.001)
    rv, tv, s2v = r[valid], t[valid], s2[valid]

    p0 = np.array([np.log(max(k0, 1e-10)), 0.0])

    def _residuals(p):
        S2sq = np.exp(2*p[0]) * np.exp(-2*p[1] * s2v)
        D1   = 0.5 * (rv**2 - S2sq * tv**2)
        D2   = 0.5 * (tv**2 - rv**2 / np.maximum(S2sq, 1e-30))
        return np.concatenate([D1, D2])

    try:
        res = least_squares(_residuals, p0, method='lm', max_nfev=400)
        p   = res.x
    except Exception:
        p = p0

    k_fit = float(np.exp(p[0]))
    B_fit = float(p[1])              # effective B; scaleit prints 4×this
    print("scale_mtz: k=%.4f  B_eff=%.2f A^2  B_scaleit=%.2f A^2  (mode=%s)" % (
        k_fit, B_fit, 4*B_fit, mode))

    if mode == 'fit_scaleB_use_scale':
        return k_fit

    # ── fit_scaleB_use_scaleB: per-reflection scale for all target reflections
    hkl_all = np.array(tgt_arr.indices())
    h_all   = hkl_all.astype(float) @ M_frac
    s2_all  = np.einsum('ni,ni->n', h_all, h_all) / 4.0
    return k_fit * np.exp(-B_fit * s2_all)


def _looks_like_fobs(label_string):
    """True if a miller array label looks like observed amplitudes."""
    for part in label_string.split(','):
        s = part.strip().upper().replace('-', '').replace('_', '')
        if s in ('FOBS', 'FP', 'FMEAS') or s.startswith('FOBS') or s.startswith('FMEAS'):
            return True
    return False


def extract_fobs_free(src_mtz, dst_mtz):
    """Write dst_mtz with only FOBS/SIGFOBS and R_FREE_FLAGS from src_mtz.
    Equivalent to: cad labin E1=FOBS E2=SIGFOBS E3=R_FREE_FLAGS

    Selects the amplitude array whose name looks like observed data (FOBS, FP, …)
    rather than relying on column order, which can put K_MASK before FOBS.
    """
    arrays = mtz_arrays(src_mtz)
    # Prefer arrays whose label looks like observed amplitudes
    fobs = next((a for a in arrays
                 if _looks_like_fobs(a.info().label_string())
                 and a.is_xray_amplitude_array() and a.sigmas() is not None), None)
    # Fall back to first amplitude+sigma array
    if fobs is None:
        fobs = next((a for a in arrays
                     if a.is_xray_amplitude_array() and a.sigmas() is not None), None)
    free = next((a for a in arrays
                 if 'FREE' in a.info().label_string().upper()), None)
    if fobs is None:
        raise RuntimeError("No Fobs/SigF found in %s" % src_mtz)
    ds = fobs.as_mtz_dataset(column_root_label='FOBS')
    if free is not None:
        ds.add_miller_array(free, column_root_label='R_FREE_FLAGS')
    ds.mtz_object().write(file_name=dst_mtz)


def apply_mtz_scale(mtzfile, scale, outfile):
    """Multiply all F-type columns in mtzfile by scale and write outfile."""
    obj = iotbx.mtz.object(file_name=mtzfile)
    for c in obj.crystals():
        for d in c.datasets():
            for col in d.columns():
                if col.type() in ('F', 'G', 'K'):
                    v = np.array(col.extract_values(), dtype=np.float64)
                    missing = np.isnan(v)
                    v[missing] = 0.0
                    scaled = flex.float((v * scale).astype(np.float32))
                    col.set_values(scaled, flex.bool((~missing).tolist()))
    obj.write(file_name=outfile)


def cell_volume(cell_params):
    a, b, c, al, be, ga = cell_params
    dtr = math.pi / 180.0
    ca, cb, cg = math.cos(al*dtr), math.cos(be*dtr), math.cos(ga*dtr)
    skew = abs(1 + 2*ca*cb*cg - ca**2 - cb**2 - cg**2)
    return a * b * c * math.sqrt(skew)


# ── vacuum level ─────────────────────────────────────────────────────────────

def find_vacuum_level(mapfile):
    from end.rapid.map_vacuum_level import find_vacuum_level as _fvl
    return _fvl(mapfile)


def _vacuum_level_here(mapfile):
    """Import map_vacuum_level.py from same directory as this script."""
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, here)
    from map_vacuum_level import find_vacuum_level as fvl
    return fvl(mapfile)


# ── map RMSD ─────────────────────────────────────────────────────────────────

def _map_rmsd_here(map_files, ref_file=None):
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, here)
    from map_rmsd import map_rmsd
    map_rmsd(map_files, ref_file=ref_file)


# ── eff file preprocessing ────────────────────────────────────────────────────

def preprocess_eff(eff_file, out_file, strip_data_manager=False):
    """
    Strip unrecognized Phenix 2.x blocks and set volume scaling.

    strip_data_manager=True additionally removes all data_manager {} and
    refinement.input.{pdb,xray_data} sections so the resulting eff can be
    used with different input files on the command line (needed for RAPID seeds).
    """
    with open(eff_file) as fh:
        lines = fh.readlines()

    # Tags to skip entirely (including their brace-delimited block)
    skip_keywords = {'ddr'}
    if strip_data_manager:
        skip_keywords.add('data_manager')

    cleaned = []
    skip_depth = 0
    # Track brace context for refinement.input.pdb / xray_data removal
    context = []    # stack of (keyword, depth_at_open)

    for line in lines:
        opens  = line.count('{')
        closes = line.count('}')
        word   = line.strip().split()[0] if line.strip() else ''

        if skip_depth > 0:
            skip_depth += opens - closes
            if skip_depth <= 0:
                skip_depth = 0
            continue

        if word in skip_keywords and '{' in line:
            skip_depth = opens - closes
            continue

        if strip_data_manager:
            # Also drop refinement { input { pdb { } } } and xray_data sub-blocks
            # by tracking context
            if word == 'pdb' and len(context) >= 2 and context[-1] == 'input' and context[-2] == 'refinement' and '{' in line:
                skip_depth = opens - closes; continue
            if word == 'xray_data' and len(context) >= 2 and context[-1] == 'input' and context[-2] == 'refinement' and '{' in line:
                skip_depth = opens - closes; continue

        # maintain context stack
        if '{' in line:
            context.append(word if word else '?')
        for _ in range(closes):
            if context:
                context.pop()

        cleaned.append(line)

    # Volume scaling fixes
    out = []
    for line in cleaned:
        s = line.strip()
        if re.search(r'\bscale\b.*\bsigma\b', s):
            line = 'scale=volume\n'
        elif re.search(r'\bapply\b.*\bscaling\b', s):
            line = 'apply_volume_scaling=True\napply_sigma_scaling=False\n'
        out.append(line)

    with open(out_file, 'w') as fh:
        fh.writelines(out)


# ── atomic Z-sum ──────────────────────────────────────────────────────────────

# Atomic numbers for elements common in macromolecular models
_ATOMIC_Z = {
    'H':1,'HE':2,'LI':3,'BE':4,'B':5,'C':6,'N':7,'O':8,'F':9,'NE':10,
    'NA':11,'MG':12,'AL':13,'SI':14,'P':15,'S':16,'CL':17,'AR':18,
    'K':19,'CA':20,'SC':21,'TI':22,'V':23,'CR':24,'MN':25,'FE':26,
    'CO':27,'NI':28,'CU':29,'ZN':30,'GA':31,'GE':32,'AS':33,'SE':34,
    'BR':35,'KR':36,'RB':37,'SR':38,'Y':39,'ZR':40,'NB':41,'MO':42,
    'TC':43,'RU':44,'RH':45,'PD':46,'AG':47,'CD':48,'IN':49,'SN':50,
    'SB':51,'TE':52,'I':53,'XE':54,'CS':55,'BA':56,'W':74,'PT':78,
    'AU':79,'HG':80,'PB':82,'U':92,
}


def z_sum_vacuum_level(pdbfile):
    """Compute vacuum level from atomic numbers — equivalent to the shell SFALL Z-sum."""
    pdb_inp = iotbx.pdb.input(file_name=pdbfile)
    xrs     = pdb_inp.xray_structure_simple()

    z_sum = 0.0
    for sc in xrs.scatterers():
        elem = sc.element_symbol().strip().upper()
        z = _ATOMIC_Z.get(elem, 0)
        if z == 0:
            print("WARNING: unknown element %s, Z=0" % elem)
        z_sum += z * sc.occupancy

    symops     = xrs.space_group().order_z()
    vol        = xrs.unit_cell().volume()
    model_vac  = -(z_sum * symops) / vol
    print("Z-sum predicted no-solvent vacuum level: %.6f" % model_vac)
    return model_vac


# ── PDB helpers ───────────────────────────────────────────────────────────────

def pdb_cell_sg(pdbfile):
    """Return (cell_params_tuple, sg_symbol) from a PDB file via cctbx."""
    pdb_inp = iotbx.pdb.input(file_name=pdbfile)
    cs = pdb_inp.crystal_symmetry()
    return tuple(cs.unit_cell().parameters()), str(cs.space_group_info())


def set_B_to_80(pdbfile, outfile):
    """Rewrite a PDB with all B-factors set to 80."""
    with open(pdbfile) as fh, open(outfile, 'w') as fo:
        for line in fh:
            if line.startswith(('ATOM  ', 'HETATM')) and not line.startswith('ANISOU'):
                line = line[:60] + ' 80.00' + line[66:]
            fo.write(line)


def has_hydrogens(pdbfile, threshold=0.5):
    counts = {}
    with open(pdbfile) as fh:
        for line in fh:
            if line.startswith(('ATOM  ', 'HETATM')):
                elem = line[76:78].strip().upper()
                counts[elem] = counts.get(elem, 0) + 1
    C = counts.get('C', 0)
    H = counts.get('H', 0)
    return C > 0 and H / C > threshold


def get_ksol_bsol(logfile):
    """Extract ksol and bsol from a phenix.refine log.

    Phenix 2.x writes a per-resolution-bin table ending with a 'kmask' column;
    the value in the lowest-resolution bin of the last such table is used as ksol.
    Phenix 1.x writes 'k_sol = X  b_sol = Y' lines.
    """
    ksol = bsol = None
    with open(logfile) as fh:
        lines = fh.readlines()

    for i, line in enumerate(lines):
        # Phenix 2.x: header row "Resolution ... kmask" followed by data rows
        if 'Resolution' in line and 'kmask' in line:
            for j in range(i + 1, min(i + 30, len(lines))):
                parts = lines[j].split()
                # data rows have a resolution range "lo-hi" as first field
                if len(parts) >= 6 and '-' in parts[0]:
                    try:
                        ksol = float(parts[-1])   # kmask is last column
                        bsol = 0.0
                    except ValueError:
                        pass
                    break
                elif not parts or parts[0].startswith('='):
                    break  # hit next section header

        # Phenix 1.x: k_sol = X  b_sol = Y
        m = re.search(r'k_sol\s*=\s*([0-9.]+).*b_sol\s*=\s*([0-9.]+)', line, re.I)
        if m:
            ksol, bsol = float(m.group(1)), float(m.group(2))

    return ksol, bsol


def get_rsolv_rshrink(logfile, efffile):
    """Extract solvent mask radii from log or eff file."""
    rsolv = rshrink = None
    for path in [logfile, efffile]:
        if not path or not os.path.isfile(path):
            continue
        with open(path) as fh:
            for line in fh:
                m = re.search(r'solvent_radius\s*=\s*([0-9.]+)', line)
                if m and rsolv is None:
                    rsolv = float(m.group(1))
                m = re.search(r'shrink_truncation_radius\s*=\s*([0-9.]+)', line)
                if m and rshrink is None:
                    rshrink = float(m.group(1))
    return rsolv, rshrink


# ── main ─────────────────────────────────────────────────────────────────────

def main(args):
    # --- parse arguments ---
    eff_file   = None
    cycles     = 5
    rapid_itrs = 5
    cpus       = 'auto'
    no_rapid   = False
    no_fofc    = False
    no_sigf    = False
    no_hadd    = False
    no_scale   = False
    ksol_arg   = None
    bsol_arg   = None

    for arg in args:
        if arg.endswith('.eff') and os.path.isfile(arg):
            eff_file = arg
        elif arg.startswith('cycles='):
            cycles = int(arg.split('=')[1])
        elif arg.startswith('seeds='):
            rapid_itrs = int(arg.split('=')[1])
        elif arg.startswith('cpus='):
            cpus = int(arg.split('=')[1])
        elif arg.startswith('ksol='):
            ksol_arg = float(arg.split('=')[1])
        elif arg.startswith('Bsol='):
            bsol_arg = float(arg.split('=')[1])
        elif arg.startswith('-norapid'):
            no_rapid = True
        elif arg.startswith('-nofofc'):
            no_fofc = True
        elif arg.startswith('-nosig'):
            no_sigf = True
        elif arg.startswith('-nohydro'):
            no_hadd = True
        elif arg.startswith('-noscale'):
            no_scale = True

    if not eff_file:
        print("ERROR: cannot find phenix.refine eff file")
        _usage(); return 9

    # --- CPU count ---
    if cpus == 'auto':
        cpus = multiprocessing.cpu_count()
    needed = min(cpus, rapid_itrs)
    print("%d CPUs available; will run %d seeds of %d cycles each"
          % (cpus, rapid_itrs, cycles))

    prefix    = 'find_F000'
    phenix_log = prefix + '_001.log'

    # --- preprocess eff file ---
    preprocess_eff(eff_file, 'find_F000_input.eff')
    os.environ['PHENIX_OVERWRITE_ALL'] = 'true'

    # --- initial phenix.refine ---
    print("running initial phenix.refine job...")
    run('phenix.refine find_F000_input.eff output.prefix=%s '
        'main.number_of_macro_cycles=%d export_final_f_model=True' % (prefix, cycles),
        log=phenix_log)

    # find output files (handle _001 and _1 naming)
    def find_output(stem, ext):
        for suf in ['_001', '_1']:
            path = stem + suf + ext
            if os.path.isfile(path):
                return path
            path2 = stem + suf + '_map_coeffs' + ext
            if os.path.isfile(path2):
                return path2
        return None

    pdbfile = find_output(prefix, '.pdb')
    mtzfile = find_output(prefix, '.mtz')
    efffile = find_output(prefix, '.eff')
    if not pdbfile or not mtzfile:
        print("ERROR: phenix.refine output not found"); return 9

    # --- traditional 2FoFc and FoFc maps (not absolute scale) ---
    print("calculating traditional 2FoFc.map")
    fft_map_to_file(mtzfile, '2FOFCWT', 'PH2FOFCWT', '2FoFc.map')

    mm_ref = map_manager(file_name='2FoFc.map')
    xyzgrid = mm_ref.unit_cell_grid

    fft_map_to_file(mtzfile, 'FOFCWT', 'PHFOFCWT', 'FoFc.map')

    # --- cell and space group ---
    cell_params, pdbSG = pdb_cell_sg(pdbfile)
    vol = cell_volume(cell_params)
    print("unit cell volume: %.3f A^3" % vol)

    # --- bulk solvent parameters ---
    ksol   = ksol_arg
    bsol   = bsol_arg if bsol_arg is not None else 0.0
    if ksol is None:
        ksol, bsol_log = get_ksol_bsol(phenix_log)
        if bsol_arg is None and bsol_log is not None:
            bsol = bsol_log
        if ksol is None:
            print("ERROR: unable to obtain k_sol"); return 9
        print("got ksol = %g  bsol = %g" % (ksol, bsol))

    rsolv, rshrink = get_rsolv_rshrink(phenix_log, efffile)
    print("bulk solvent params: ksol=%.4f  Bsol=%.1f  Rsolv=%s  Rshrink=%s"
          % (ksol, bsol, rsolv, rshrink))

    phenix_reso = None
    with open(pdbfile) as fh:
        for line in fh:
            if 'RESOLUTION RANGE HIGH' in line:
                phenix_reso = float(line.split()[-1]) * 0.98
    if phenix_reso is None:
        # fall back to cctbx resolution from the MTZ
        phenix_reso = max(a.d_min() for a in mtz_arrays(mtzfile))
        print("WARNING: no RESOLUTION RANGE remark; using d_min=%.4f from MTZ" % phenix_reso)

    # --- set B=80 and add hydrogens ---
    print("setting all B factors to 80")
    set_B_to_80(pdbfile, 'refined.pdb')

    if not no_hadd and not has_hydrogens(pdbfile):
        print("Adding hydrogens...")
        run('phenix.ready_set refined.pdb', log='ready_set.log')
        if os.path.isfile('refined.updated.pdb'):
            shutil.move('refined.updated.pdb', 'refined.pdb')
        else:
            print("ERROR: unable to add hydrogens"); return 9

    # --- phenix.fmodel runs ---
    extra_mask = ''
    if rsolv   is not None: extra_mask += ' mask.solvent_radius=%.3f'          % rsolv
    if rshrink is not None: extra_mask += ' mask.shrink_truncation_radius=%.3f' % rshrink

    print("generating phenix models...")
    print("no solvent")
    run('phenix.fmodel refined.pdb k_sol=0 b_sol=0 high_res=%.4f file_name=nobulk.mtz'
        % phenix_reso, log='phenix_fmodel.log')

    print("sharp solvent")
    run('phenix.fmodel refined.pdb k_sol=%.4f b_sol=0 high_res=%.4f file_name=ss.mtz%s'
        % (ksol, phenix_reso, extra_mask), log='phenix_fmodel.log')

    print("regular")
    run('phenix.fmodel %s k_sol=%.4f b_sol=%.4f high_res=%.4f file_name=fmodel.mtz%s'
        % (pdbfile, ksol, bsol, phenix_reso, extra_mask), log='phenix_fmodel.log')

    # --- sfall-equivalent: SFALL via phenix.fmodel with no B-factor damping ---
    # Use the vacuum-level from the no-bulk model (mean of map = -vacuum)
    print("making sfall map for vacuum level reference")
    fft_map_to_file('nobulk.mtz', 'FMODEL', 'PHIFMODEL', 'nobulk.map')

    # --- compute atomic vacuum level from Z-sum (F000 not in diffraction data) ---
    model_vac = z_sum_vacuum_level('refined.pdb')

    # --- fmodel maps ---
    fft_map_to_file('ss.mtz',     'FMODEL', 'PHIFMODEL', 'fmodel_ss.map')
    fft_map_to_file('fmodel.mtz', 'FMODEL', 'PHIFMODEL', 'fmodel.map')

    # --- vacuum level of sharp-solvent-only region ---
    print("subtracting nobulk map from sharp-solvent map to obtain solvent-only map")
    nb_data = np.array(map_manager(file_name='nobulk.map').map_data())
    ss_data = np.array(map_manager(file_name='fmodel_ss.map').map_data())
    solvent_data = ss_data - nb_data
    _write_map(map_manager(file_name='nobulk.map'), solvent_data, 'sharpsolvent.map')

    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, here)
    from map_vacuum_level import find_vacuum_level as fvl
    print("finding lower median with outlier rejection")
    solvent_vac, solvent_sig = fvl('sharpsolvent.map')
    sharpsolvent_vac_log = open('sharpsolvent_vac.log', 'w')
    sharpsolvent_vac_log.write("vacuum level: %.6f +/- %.6f\n" % (solvent_vac, solvent_sig))
    sharpsolvent_vac_log.close()

    # --- absolute-scale map coefficients ---
    baseline = -model_vac - solvent_vac
    baseline_sig = math.sqrt(0 + solvent_sig**2)
    SIGF000 = vol * baseline_sig

    # Scale map coefficients: scaleit(fmodel vs 2FOFCWT)
    print("putting map coefficients on absolute scale")
    map_scale_factor = scale_mtz('fmodel.mtz', 'FMODEL', mtzfile, '2FOFCWT')
    if no_scale:
        print("but no scaling will be applied.")
        map_scale_factor = 1.0
    print("scale factor: %.4f" % map_scale_factor)

    print("applying %.4f to map coefficients" % map_scale_factor)
    apply_mtz_scale(mtzfile, map_scale_factor, 'scaled_map_coeffs.mtz')

    print("calculating 2FoFc_scaled.map on absolute scale")
    fft_map_to_file('scaled_map_coeffs.mtz', '2FOFCWT', 'PH2FOFCWT', '2FoFc_scaled.map')

    print("calculating FoFc_scaled.map on absolute scale")
    fft_map_to_file('scaled_map_coeffs.mtz', 'FOFCWT', 'PHFOFCWT', 'FoFc_scaled.map')

    print("adding %.5f to 2FoFc_scaled.map to form 2FoFc_END.map" % baseline)
    map_add_offset('2FoFc_scaled.map', baseline, '2FoFc_END.map')

    end_mean = float(np.array(map_manager(file_name='2FoFc_END.map').map_data()).mean())
    F000 = vol * end_mean
    print("overall F000 = %.1f +/- %.2f" % (F000, SIGF000))
    print("2FoFc_END.map is on an absolute electron number-density scale (electrons/A^3)")
    print("FoFc_scaled.map is on the same scale")

    # --- place Fobs on absolute scale (kickme.mtz) ---
    fmodelfile = find_output(prefix, '_f_model.mtz')
    fobs_scale = scale_mtz('fmodel.mtz', 'FMODEL', fmodelfile, 'FOBS')
    apply_mtz_scale(fmodelfile, fobs_scale, 'kickme.mtz')
    print("scale factor (Fobs): %.4f" % fobs_scale)

    # strip eff file for RAPID seeds (remove input file paths)
    if efffile and os.path.isfile(efffile):
        preprocess_eff(efffile, prefix + '.eff', strip_data_manager=True)

    if no_rapid:
        return _finish()

    # ── RAPID: Fo-Fc perturbation ─────────────────────────────────────────────
    if not no_fofc:
        print("\ncomputing RAPID map using Fo-Fc as the error to propagate")
        seeds = list(range(1, rapid_itrs + 1))
        with open('seeds.txt', 'w') as fh:
            fh.write('\n'.join(str(s) for s in seeds) + '\n')

        import shutil as _sh
        _sh.copy(pdbfile.replace('001', '001').replace('_1.pdb', '_001.pdb'), 'refme.pdb')
        # remove hexdigest lines
        with open(pdbfile) as rf, open('refme.pdb', 'w') as wf:
            for line in rf:
                if 'hexdigest' not in line:
                    wf.write(line)

        from kick_data_bydiff import kick_data_bydiff

        procs = []
        for seed in seeds:
            kick_data_bydiff('kickme.mtz', seed=seed, FC_label='FMODEL', output='kicked_%d.mtz' % seed)
            extract_fobs_free('kicked_%d.mtz' % seed, 'refme_FOFC_%d.mtz' % seed)

            cmd = ('phenix.refine %s.eff refme.pdb refme_FOFC_%d.mtz '
                   'export_final_f_model=True write_geo_file=False write_def_file=False '
                   'main.number_of_macro_cycles=%d output.prefix=seed_%d'
                   % (prefix, seed, cycles, seed))
            p = subprocess.Popen(cmd, shell=True,
                                 stdout=open('seed_%d.log' % seed, 'w'),
                                 stderr=subprocess.STDOUT)
            procs.append((seed, p))
            if len(procs) >= needed:
                for s, pr in procs:
                    pr.wait()
                    print("seed %d done" % s)
                procs = []
        for s, pr in procs:
            pr.wait()
            print("seed %d done" % s)

        # scale and FFT each seed map
        for seed in seeds:
            smtz = find_output('seed_%d' % seed, '.mtz')
            if smtz is None:
                raise RuntimeError("seed_%d phenix.refine failed — check seed_%d.log" % (seed, seed))
            sc   = scale_mtz('fmodel.mtz', 'FMODEL', smtz, '2FOFCWT')
            print("seed %d scale factor: %.4f" % (seed, sc))
            fft_map_to_file(smtz, '2FOFCWT', 'PH2FOFCWT',
                            'seed_%d_2FoFc.map' % seed, scale=sc)
            fft_map_to_file(smtz, 'FOFCWT',  'PHFOFCWT',
                            'seed_%d_FoFc.map'  % seed, scale=sc)

        from map_rmsd import map_rmsd
        map_rmsd(sorted(glob.glob('seed_*_2FoFc.map')), ref_file='2FoFc_scaled.map')
        os.rename('sigma.map', '2FoFc_error.map')
        print("2FoFc_error.map is the RAPID map of error bars propagated from Fobs-Fcalc")
        map_stats('2FoFc_error.map')

        map_rmsd(sorted(glob.glob('seed_*_FoFc.map')), ref_file='FoFc_scaled.map')
        os.rename('sigma.map', 'FoFc_error.map')
        print("FoFc_error.map is the RAPID map of error bars propagated from Fobs-Fcalc")
        map_stats('FoFc_error.map')

        map_ratio('2FoFc_END.map', '2FoFc_error.map', '2FoFc_snr.map')
        map_ratio('FoFc_scaled.map', 'FoFc_error.map', 'FoFc_snr.map')

    # ── RAPID: sigF perturbation ──────────────────────────────────────────────
    if not no_sigf:
        print("\ncomputing RAPID map using sigF as the error to propagate")
        from kick_data import kick_data as _kick

        procs = []
        for seed in range(1, rapid_itrs + 1):
            _kick('kickme.mtz', seed=seed, output='kicked_sigF_%d.mtz' % seed)
            extract_fobs_free('kicked_sigF_%d.mtz' % seed, 'refme_sigF_%d.mtz' % seed)

            cmd = ('phenix.refine %s.eff refme.pdb refme_sigF_%d.mtz '
                   'export_final_f_model=True write_geo_file=False write_def_file=False '
                   'main.number_of_macro_cycles=%d output.prefix=sigF_seed_%d'
                   % (prefix, seed, cycles, seed))
            p = subprocess.Popen(cmd, shell=True,
                                 stdout=open('sigF_seed_%d.log' % seed, 'w'),
                                 stderr=subprocess.STDOUT)
            procs.append((seed, p))
            if len(procs) >= needed:
                for s, pr in procs:
                    pr.wait()
                    print("sigF seed %d done" % s)
                procs = []
        for s, pr in procs:
            pr.wait()
            print("sigF seed %d done" % s)

        for seed in range(1, rapid_itrs + 1):
            smtz = find_output('sigF_seed_%d' % seed, '.mtz')
            if smtz is None:
                raise RuntimeError("sigF_seed_%d phenix.refine failed — check sigF_seed_%d.log" % (seed, seed))
            sc   = scale_mtz('fmodel.mtz', 'FMODEL', smtz, '2FOFCWT')
            print("sigF seed %d scale factor: %.4f" % (seed, sc))
            fft_map_to_file(smtz, '2FOFCWT', 'PH2FOFCWT',
                            'sigF_seed_%d_2FoFc.map' % seed, scale=sc)
            fft_map_to_file(smtz, 'FOFCWT',  'PHFOFCWT',
                            'sigF_seed_%d_FoFc.map'  % seed, scale=sc)

        from map_rmsd import map_rmsd
        map_rmsd(sorted(glob.glob('sigF_seed_*_2FoFc.map')), ref_file='2FoFc_scaled.map')
        os.rename('sigma.map', '2FoFc_sigF_error.map')
        print("2FoFc_sigF_error.map is the RAPID map of error bars propagated from sigF")
        map_stats('2FoFc_sigF_error.map')

        map_rmsd(sorted(glob.glob('sigF_seed_*_FoFc.map')), ref_file='FoFc_scaled.map')
        os.rename('sigma.map', 'FoFc_sigF_error.map')
        print("FoFc_sigF_error.map is the RAPID map of error bars propagated from sigF")
        map_stats('FoFc_sigF_error.map')

        map_ratio('2FoFc_END.map', '2FoFc_sigF_error.map', '2FoFc_sigF_snr.map')
        map_ratio('FoFc_scaled.map', 'FoFc_sigF_error.map', 'FoFc_sigF_snr.map')

    return _finish()


def _finish():
    for name, desc in [
        ('2FoFc_END.map',          'absolute electron number-density scale'),
        ('2FoFc_error.map',        'RAPID error bars from Fobs-Fcalc'),
        ('FoFc_error.map',         'RAPID error bars from Fobs-Fcalc'),
        ('2FoFc_sigF_error.map',   'RAPID error bars from sigFobs'),
        ('FoFc_sigF_error.map',    'RAPID error bars from sigFobs'),
        ('2FoFc_snr.map',          'ratio of 2FoFc_END to 2FoFc_error'),
        ('FoFc_snr.map',           'ratio of FoFc_scaled to FoFc_error'),
        ('2FoFc_sigF_snr.map',     'ratio of 2FoFc_END to 2FoFc_sigF_error'),
        ('FoFc_sigF_snr.map',      'ratio of FoFc_scaled to 2FoFc_sigF_error'),
    ]:
        if os.path.isfile(name):
            print("%s is the %s" % (name, desc))
    return 0


def _usage():
    print("""usage:
  phenix.python END_RAPID.py phenixrefine.eff [options]

options:
  seeds=n      RAPID refinements per map (default 5)
  cycles=n     macro cycles per refinement (default 5)
  cpus=n       parallel jobs (default: all available)
  ksol=x       override bulk-solvent k_sol
  Bsol=x       override bulk-solvent B_sol
  -nofofc      skip Fo-Fc RAPID maps
  -nosigf      skip sigF RAPID maps
  -norapid     skip all RAPID maps
""")


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]) or 0)
