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


def fft_map_to_file(mtzfile, F_label, PHI_label, outfile, grid=None):
    """FFT a (F, PHI) pair from an MTZ into a CCP4 map (volume-scaled, ASU)."""
    obj = iotbx.mtz.object(file_name=mtzfile)
    target = None
    for arr in obj.as_miller_arrays():
        ls = arr.info().label_string()
        if F_label in ls and PHI_label in ls:
            target = arr
            break
        # Also accept F,PHI stored as separate columns merged into a complex array
        if F_label in ls:
            # try amplitude+phase pair
            target = arr
    if target is None:
        raise RuntimeError("Cannot find %s/%s in %s" % (F_label, PHI_label, mtzfile))

    if not hasattr(target, 'phases') and target.is_real_array():
        # amplitude only — find companion phase column
        for arr in obj.as_miller_arrays():
            ls = arr.info().label_string()
            if PHI_label in ls:
                target = target.phase_transfer(arr, deg=True)
                break

    fmap = target.fft_map(resolution_factor=0.25)
    fmap.apply_volume_scaling()
    rmap = fmap.real_map_unpadded()

    mm = map_manager(
        map_data=rmap,
        unit_cell_crystal_symmetry=target.crystal_symmetry(),
        unit_cell_grid=rmap.focus(),
        wrapping=True,
    )
    mm.write_map(file_name=outfile)
    return mm


def map_add_offset(infile, offset, outfile):
    """Write infile + offset as outfile."""
    mm = map_manager(file_name=infile)
    d = np.array(mm.map_data()) + offset
    mm.customized_copy(map_data=flex.double(d.flatten())).write_map(file_name=outfile)


def map_scale(infile, scale, outfile):
    mm = map_manager(file_name=infile)
    d = np.array(mm.map_data()) * scale
    mm.customized_copy(map_data=flex.double(d.flatten())).write_map(file_name=outfile)


def map_ratio(num_file, denom_file, outfile):
    """Write per-voxel num/denom, skipping zero-denominator voxels."""
    mm_n = map_manager(file_name=num_file)
    mm_d = map_manager(file_name=denom_file)
    n = np.array(mm_n.map_data())
    d = np.array(mm_d.map_data())
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(d != 0, n / d, 0.0)
    mm_n.customized_copy(map_data=flex.double(ratio.flatten())).write_map(file_name=outfile)


def map_stats(mapfile):
    d = np.array(map_manager(file_name=mapfile).map_data())
    print("           Minimum density .................................  %12g" % d.min())
    print("           Maximum density .................................  %12g" % d.max())
    print("           Mean density ....................................  %12g" % d.mean())
    print("           Rms deviation from mean density .................  %12g" % d.std())


def scale_mtz(ref_mtz, ref_label, target_mtz, target_label):
    """Return least-squares scale s such that s*target ≈ ref (isotropic)."""
    ref_arr = None
    tgt_arr = None
    for arr in mtz_arrays(ref_mtz):
        if ref_label in arr.info().label_string() and arr.is_real_array():
            ref_arr = arr; break
    for arr in mtz_arrays(target_mtz):
        if target_label in arr.info().label_string() and arr.is_real_array():
            tgt_arr = arr; break
    if ref_arr is None or tgt_arr is None:
        return 1.0
    common = ref_arr.common_set(tgt_arr)
    if common.size() == 0:
        return 1.0
    r = common.data().as_numpy_array()
    t = tgt_arr.common_set(ref_arr).data().as_numpy_array()
    scale = np.dot(r, t) / max(np.dot(t, t), 1e-12)
    return float(scale)


def apply_mtz_scale(mtzfile, scale, outfile):
    """Multiply all F-type columns in mtzfile by scale and write outfile."""
    obj = iotbx.mtz.object(file_name=mtzfile)
    for c in obj.crystals():
        for d in c.datasets():
            for col in d.columns():
                if col.type() in ('F', 'G', 'K'):
                    v = np.array(col.extract_values()) * scale
                    col.set_values(flex.double(v))
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

def preprocess_eff(eff_file, out_file):
    """
    Strip unrecognized Phenix 2.x blocks (e.g. ddr{}) and set volume scaling.
    """
    with open(eff_file) as fh:
        lines = fh.readlines()

    cleaned = []
    skip_depth = 0
    brace_depth = 0
    skip_keywords = {'ddr'}

    for line in lines:
        opens  = line.count('{')
        closes = line.count('}')

        if skip_depth == 0:
            word = line.strip().split()[0] if line.strip() else ''
            if word in skip_keywords and '{' in line:
                skip_depth = opens - closes
                continue
            cleaned.append(line)
        else:
            skip_depth += opens - closes
            if skip_depth <= 0:
                skip_depth = 0
            continue

    # Volume scaling fixes
    out = []
    for line in cleaned:
        s = line.strip()
        if re.search(r'\bscale\b.*\bsigma\b', s):
            line = 'scale=volume\n'
        elif re.search(r'\bapply\b.*\bscaling\b', s):
            line = 'apply_volume_scaling=True\napply_sigma_scaling=False\n'
        if 'serial' in s and '=' in s:
            out.append('serial_format = "%03d"\n')
        line = line.replace('%d', '%03d')
        out.append(line)

    with open(out_file, 'w') as fh:
        fh.writelines(out)


# ── PDB helpers ───────────────────────────────────────────────────────────────

def pdb_cell_sg(pdbfile):
    """Return (cell_params_tuple, sg_symbol) from CRYST1 record."""
    with open(pdbfile) as fh:
        for line in fh:
            if line.startswith('CRYST1'):
                cell = tuple(float(line[6+i*9:6+(i+1)*9]) for i in range(6))
                sg   = line[55:66].strip()
                return cell, sg
    raise RuntimeError("No CRYST1 in %s" % pdbfile)


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
    """Extract ksol and bsol from a phenix.refine log."""
    ksol = bsol = None
    with open(logfile) as fh:
        for line in fh:
            if 'kmask' in line.lower():
                # Phenix 1.8+ style: line after "kmask" header has the value
                pass
            m = re.search(r'kmask\s*=?\s*([0-9.]+)', line, re.I)
            if m and ksol is None:
                ksol = float(m.group(1))
                bsol = 0.0
            m2 = re.search(r'k_sol\s*=\s*([0-9.]+).*b_sol\s*=\s*([0-9.]+)', line, re.I)
            if m2:
                ksol, bsol = float(m2.group(1)), float(m2.group(2))
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

    with open(pdbfile) as fh:
        for line in fh:
            if 'RESOLUTION RANGE HIGH' in line:
                phenix_reso = float(line.split()[-1]) * 0.98

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

    mm_nb = map_manager(file_name='nobulk.map')
    model_vac = -float(np.array(mm_nb.map_data()).mean())
    print("no-solvent mean=0 map vacuum level: %.5f" % model_vac)

    # --- fmodel maps ---
    fft_map_to_file('ss.mtz',     'FMODEL', 'PHIFMODEL', 'fmodel_ss.map')
    fft_map_to_file('fmodel.mtz', 'FMODEL', 'PHIFMODEL', 'fmodel.map')

    # --- vacuum level of sharp-solvent-only region ---
    print("subtracting nobulk map from sharp-solvent map to obtain solvent-only map")
    nb_data = np.array(map_manager(file_name='nobulk.map').map_data())
    ss_data = np.array(map_manager(file_name='fmodel_ss.map').map_data())
    solvent_data = ss_data - nb_data
    mm_ref2 = map_manager(file_name='nobulk.map')
    mm_ref2.customized_copy(
        map_data=flex.double(solvent_data.flatten())
    ).write_map(file_name='sharpsolvent.map')

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
        preprocess_eff(efffile, prefix + '.eff')

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
            kick_data_bydiff('kickme.mtz', seed=seed, output='kicked_%d.mtz' % seed)
            import shutil as _sh2
            _sh2.copy('kicked_%d.mtz' % seed, 'refme_FOFC_%d.mtz' % seed)

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
            sc   = scale_mtz('fmodel.mtz', 'FMODEL', smtz, '2FOFCWT')
            print("seed %d scale factor: %.4f" % (seed, sc))
            fft_map_to_file(smtz, '2FOFCWT', 'PH2FOFCWT',
                            'seed_%d_2FoFc.map' % seed)
            fft_map_to_file(smtz, 'FOFCWT',  'PHFOFCWT',
                            'seed_%d_FoFc.map'  % seed)

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
            import shutil as _sh3
            _sh3.copy('kicked_sigF_%d.mtz' % seed, 'refme_sigF_%d.mtz' % seed)

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
            sc   = scale_mtz('fmodel.mtz', 'FMODEL', smtz, '2FOFCWT')
            print("sigF seed %d scale factor: %.4f" % (seed, sc))
            fft_map_to_file(smtz, '2FOFCWT', 'PH2FOFCWT',
                            'sigF_seed_%d_2FoFc.map' % seed)
            fft_map_to_file(smtz, 'FOFCWT',  'PHFOFCWT',
                            'sigF_seed_%d_FoFc.map'  % seed)

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
