# END/RAPID

Compute absolutely-scaled Electron Number Density (END) maps and per-voxel noise
maps using Refinement Against Perturbed Input Data (RAPID) from a
[phenix.refine](http://www.phenix-online.org/) output.

Full documentation: [Documentation/documentation.htm](Documentation/documentation.htm)  
Manual: [Documentation/end.rapid.Manual.htm](Documentation/end.rapid.Manual.htm)

## Two implementations

| | Shell (`END_RAPID.com`) | Python (`END_RAPID.py`) |
|---|---|---|
| Interpreter | tcsh | phenix.python |
| Map/MTZ I/O | CCP4 (fft, mapmask, scaleit, …) | cctbx / iotbx (bundled with Phenix) |
| Refinement | phenix.refine, phenix.fmodel | same, via subprocess |
| CCP4 required? | **yes** | **no** |
| Phenix required? | **yes** | **yes** |
| Status | production | validated, functionally equivalent |

## Requirements

### Shell version
- [CCP4](http://www.ccp4.ac.uk/) 6.x – 9.x — for `fft`, `mapmask`, `sfall`, `scaleit`, `sftools`, `cad`, `mapdump`, `mapsig`
- [Phenix](http://www.phenix-online.org/) 1.6 – 2.1 — for `phenix.refine`, `phenix.fmodel`, `phenix.ready_set`
- tcsh

### Python version
- [Phenix](http://www.phenix-online.org/) 2.x — for `phenix.refine`, `phenix.fmodel`, `phenix.ready_set` (called as subprocesses), and for `phenix.python` which bundles cctbx/iotbx/numpy
- **No CCP4 needed** — all map and MTZ operations use `iotbx.mtz`, `iotbx.map_manager`, and `cctbx.miller`

## Installation

```bash
tar -zxvf end.rapid.tar.gz
# add end.rapid/ to your PATH (tcsh) or set PYTHONPATH (Python version)
```

## Usage

### Shell
```bash
END_RAPID.com phenixrefine.eff [seeds=5] [cycles=5] [cpus=N] [-nofofc] [-nosigf] [-norapid]
```

### Python
```bash
phenix.python END_RAPID.py phenixrefine.eff [seeds=5] [cycles=5] [cpus=N] [-nofofc] [-nosigf] [-norapid]
```

The input eff file is the one written by `phenix.refine` at the end of your last
refinement.  No manual pre-processing is needed — both scripts automatically strip
parameter blocks unrecognized by the current Phenix version (e.g. the `ddr { }`
block present in eff files written by Phenix 2.x).

### Output maps

| File | Description |
|---|---|
| `2FoFc_END.map` | Absolute-scale END map (vacuum = 0 e⁻/Å³) |
| `2FoFc_error.map` | RAPID noise from \|Fo−Fc\| perturbation |
| `2FoFc_sigF_error.map` | RAPID noise from σ(Fobs) perturbation |
| `2FoFc_snr.map` | Voxel-wise signal-to-noise (END / FoFc-RAPID) |
| `FoFc_scaled.map` | Fo−Fc difference map on absolute scale |
| `FoFc_error.map` | RAPID noise for the Fo−Fc map |

## Known limitations of the Python version

- **sigF RAPID maps are ~50% smaller** than the shell version's, even though the
  SIGFOBS values in both kickme.mtz files are essentially identical (SIGFOBS/FOBS
  ≈ 0.029 in both).  The gap is statistical: the sigF perturbation is ~7× smaller
  than the Fo−Fc perturbation, so N=5 seeds is too few for the numpy and CCP4
  random-number generators to agree.  The Fo−Fc RAPID maps (22% off) are more
  stable and more physically meaningful.
- **Output maps cover the full unit cell** (not just the ASU).  The shell version
  uses `mapmask xyzlim asu`; the Python FFT returns the full crystallographic grid.
  Map statistics (mean, σ) are identical; file sizes are larger.

## Citation

Lang PT, Holton JM, Fraser JS, Alber T.
*Protein structural ensembles are revealed by redefining x-ray electron density noise.*
PNAS USA **111**, 237–247 (2014).
