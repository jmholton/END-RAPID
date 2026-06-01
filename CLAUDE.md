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

## Dependencies

| Tool | Version tested | Purpose |
|---|---|---|
| CCP4 | 6.x – 9.x | fft, mapmask, sfall, scaleit, sftools (shell scripts) |
| Phenix | 1.6 – 2.1 | phenix.refine, phenix.fmodel, phenix.ready_set |
| phenix.python | 2.1 | Python scripts (numpy, cctbx, iotbx) |

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
refinement.  Both scripts handle Phenix 2.x eff files automatically (the
`ddr { }` block is stripped before calling phenix.refine).

## Phenix 2.x compatibility fix

Phenix 2.x eff files contain a `ddr { }` parameter block that Phenix 2.x's
own `phenix.refine` no longer accepts as input.  Both `END_RAPID.com` and
`END_RAPID.py` strip this block automatically during eff pre-processing:

```awk
# shell version (inline in END_RAPID.com)
awk '/ddr \{/{skip=1} skip{if(/\}/) --skip; next} {print}'

# Python version (inline in END_RAPID.py preprocess_eff())
skip_keywords = {'ddr'}
```

## Validated test case

PDB 3ldc (1.45 Å, P4₂2₂), seeds=5, cycles=5, CCP4 9.0.015 / Phenix 2.1rc2:

| Output | Value |
|---|---|
| F000 | 63,599 ± 1,541 e⁻ |
| RAPID noise (Fo-Fc, mean) | 0.146 e⁻/Å³ |
| RAPID noise (sigF, mean)  | 0.034 e⁻/Å³ |
| 1σ / sigF-RAPID           | ~10× |
| Exit status | 0 |

## Adding new Phenix-incompatible parameters

If a future Phenix version introduces another unrecognised block, add its
keyword to `skip_keywords` in `preprocess_eff()` (Python) or extend the
`awk` pattern in `END_RAPID.com`.

## Python script design notes

- `map_vacuum_level.py`, `map_rmsd.py` — pure numpy on data from
  `iotbx.map_manager`; no CCP4 required
- `kick_data.py`, `kick_data_bydiff.py` — column-level MTZ manipulation
  via `iotbx.mtz`; no CCP4 required
- `END_RAPID.py` — cctbx for FFT, map stats, and scaling; subprocess calls
  for `phenix.refine`, `phenix.fmodel`, `phenix.ready_set`
