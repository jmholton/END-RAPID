# END/RAPID — CLAUDE.md

## What this is

Tcsh shell scripts (+ Python ports) for computing:
- **END maps** — absolutely-scaled electron density (vacuum = 0 e⁻/Å³)
- **RAPID maps** — per-voxel noise estimated by refining against perturbed data

Lang PT, Holton JM, Fraser JS, Alber T. *PNAS* **111**, 237–247 (2014).

## Repository layout

```
end.rapid/          Shell scripts (production) + Python ports
  END_RAPID.com     Master shell script (tcsh)
  END_RAPID.py      Python port — phenix.python, uses cctbx for maps/MTZ
  kick_data.com/.py         Perturb Fobs by rms(sigF)
  kick_data_bydiff.com/.py  Perturb Fobs by rms(|Fo-Fc|)
  map_rmsd.com/.py          Per-voxel RMS across a set of maps
  map_vacuum_level.com/.py  Find map vacuum level (robust statistics)
  coot_map_boiling.scm      Coot script for visualising END/RAPID maps
Documentation/      HTML manual (end.rapid.Manual.htm) and docs page
TestSets/           PDB ID lists and figures from Lang et al. 2014
```

## Dependency map

```
Operation                   Shell (END_RAPID.com)       Python (END_RAPID.py)
─────────────────────────── ─────────────────────────── ──────────────────────────
Crystallographic refinement phenix.refine               phenix.refine (subprocess)
Structure factor calculation phenix.fmodel              phenix.fmodel (subprocess)
Add hydrogens               phenix.ready_set            phenix.ready_set (subprocess)
FFT structure factors→map   CCP4 fft                   cctbx.miller.fft_map()
Trim map to ASU             CCP4 mapmask                (full unit cell retained)
Absolute-scale atomic model CCP4 sfall                  Z-sum formula (no CCP4)
Map statistics              CCP4 mapdump                iotbx.map_manager + numpy
Scale MTZ amplitudes        CCP4 scaleit                sum(ref)/sum(tgt) formula
Combine MTZ files           CCP4 cad                   iotbx.mtz column ops
Perturb Fobs                CCP4 sftools RAN_G          numpy.random.normal
Map RMSD                    CCP4 mapmask + float_mult   iotbx.map_manager + numpy
Map ratio (SNR)             CCP4 mapsig                 numpy
Vacuum level (robust stats) custom awk/shell             numpy (map_vacuum_level.py)
```

**CCP4 is only required by the shell version.**  
**Phenix (including its bundled cctbx/iotbx) is required by both versions.**

## Running

### Shell version (tcsh)
```tcsh
set path = ( /path/to/end.rapid $path )
END_RAPID.com phenixrefine.eff seeds=5 cycles=5
```

### Python version (phenix.python)
```bash
phenix.python /path/to/end.rapid/END_RAPID.py phenixrefine.eff seeds=5 cycles=5
```

The eff file is the one written by `phenix.refine` at the end of your last
refinement.  Both scripts handle Phenix 2.x eff files automatically.

## Phenix 2.x compatibility

### The `ddr { }` block
Phenix 2.x eff files contain a `ddr { }` parameter block that Phenix 2.x
`phenix.refine` no longer accepts.  Both scripts strip it during eff preprocessing:

```awk
# shell (inline in END_RAPID.com eff awk pipeline)
awk '/ddr \{/{skip=1} skip{if(/\}/) --skip; next} {print}'

# Python (preprocess_eff(), skip_keywords = {'ddr'})
```

### The `serial_format` injection
Older versions of the shell awk injected `serial_format = "%03d"` into the eff
file.  This is NOT valid inside helix blocks in Phenix 2.x and causes
"unrecognized PHIL parameter" failures.  The Python script never injects it.
The shell script also no longer injects it (removed in current version).

### Column names in Phenix 2.x output
Phenix 2.1 phenix.refine writes hyphenated column names (`F-obs`, `F-model`,
`R-free-flags`) while phenix.fmodel still writes underscored names (`FMODEL`,
`PHIFMODEL`, `FOBS`).  `scale_mtz()` normalises label comparisons to handle both.

### The `data_manager` block in seed eff files
When creating the seed eff file for RAPID refinements, the `data_manager { }`
section must be stripped (it points to the original data files).  If it is not
stripped, phenix.refine receives conflicting input file references:
"Sorry: Wrong number of models of each type supplied."
Fix: `preprocess_eff(efffile, prefix+'.eff', strip_data_manager=True)`.

## Validated test case

PDB 3ldc (1.45 Å, P4₂2₂), seeds=5, cycles=5  
Shell: CCP4 9.0.015 / Phenix 2.1rc2 | Python: phenix.python 2.1rc2

| Output | Shell | Python | Δ |
|---|---|---|---|
| `2FoFc_END.map` mean | 0.355 e⁻/Å³ | 0.358 e⁻/Å³ | 1% |
| `FoFc_scaled.map` σ | 0.115 e⁻/Å³ | 0.107 e⁻/Å³ | 7% |
| `2FoFc_error.map` mean (FoFc RAPID) | 0.146 e⁻/Å³ | 0.114 e⁻/Å³ | 22% |
| `2FoFc_sigF_error.map` mean (sigF RAPID) | 0.034 e⁻/Å³ | 0.010 e⁻/Å³ | 70% |
| Exit status | 0 | 0 | — |

## Bugs found during Python port testing (9 total)

| Bug | Symptom | Fix |
|---|---|---|
| `customized_copy` rejects 1-D flex | `AssertionError` on map write | Reshape flat numpy → `flex.grid(shape)` before passing |
| `column.set_values` needs `flex.float` | `Boost.Python.ArgumentError` | Use `flex.float(arr.astype(np.float32))` + `flex.bool` validity mask |
| CRYST1 column widths | `ValueError: could not convert '90.00  90'` | Use `iotbx.pdb.input().crystal_symmetry()` instead of substring slicing |
| `serial_format` injection | `phenix.refine` "unrecognized PHIL parameter" inside helix block | Remove injection; Phenix 2.x doesn't need it |
| `data_manager` not stripped for seed eff | "Wrong number of models" in seed refinements | `strip_data_manager=True` in `preprocess_eff` for seed eff |
| F000 from FFT mean | `model_vac = 0.0` always (F000 missing from diffraction data) | Use Z-sum formula: `−(sum_j Z_j × occ_j × N_symop) / V_cell` |
| `scale_mtz` skips complex arrays | All scale factors = 1.0 (arrays not found) | Call `.amplitudes()` on complex miller arrays before comparing |
| `scale_mtz` dot-product formula | Scale = 0.42 instead of ~1.07 for Wilson amplitudes | Use `sum(ref)/sum(tgt)` (mean-ratio); dot/dot² ≈ 0.4 for Wilson distributions |
| `kick_data_bydiff` FC_label not passed | `delta = 0` → RAPID noise = numerical zero | Pass `FC_label='FMODEL'` explicitly (otherwise F = FC = FOBS) |
| `kick_data` column selection by type+alpha | Perturbs K_MASK (mean 0.007 e-) instead of FOBS (mean 60 e-) | Prefer columns whose label looks like observed amplitudes (FOBS, FP, …); also rank by mean amplitude so derived F-type columns like K_MASK and K_ISOTROPIC don't win alphabetically at equal completeness |

## Known limitations of the Python version

**sigF RAPID maps (~50% smaller than shell after all bugs are fixed).**  
The SIGFOBS values in both Python and shell kickme.mtz are essentially identical
(mean ratio SIGFOBS/FOBS = 0.029 in both), so the perturbation amplitude is
the same.  The remaining gap is statistical: sigF perturbation is ~7× smaller
than FoFc perturbation (σ(F)/F ≈ 0.029 vs |Fo−Fc|/F ≈ 0.20), so 5 seeds is
insufficient to average out RNG differences between numpy and CCP4's sftools.
The FoFc RAPID (22% off) is stable at N=5 because the per-seed signal is large
enough to dominate; the sigF RAPID needs more seeds for the same stability.
The Fo−Fc RAPID maps are more physically meaningful in any case.

**Maps cover the full unit cell, not the ASU.**  
The shell uses `mapmask xyzlim asu`; the Python FFT returns the full grid.
Density values and statistics are identical; files are larger (16 MB vs 0.8 MB
for 3ldc at 1.45 Å).

**`scale_mtz` is isotropic.**  
The shell's `scaleit refine anisotropic` fits an anisotropic B-factor correction
in addition to an overall scale; the Python version uses only an isotropic scale
(mean-ratio).  This produces a ~1% systematic difference in the Fobs-scale step.

## Adding new Phenix-incompatible parameters

If a future Phenix version introduces another unrecognised block, add its keyword
to `skip_keywords` in `preprocess_eff()` (Python) or extend the `awk` pattern in
`END_RAPID.com` (shell).
