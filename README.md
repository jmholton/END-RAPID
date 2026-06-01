# END/RAPID

Compute absolutely-scaled Electron Number Density (END) maps and per-voxel noise
maps using Refinement Against Perturbed Input Data (RAPID) from a
[phenix.refine](http://www.phenix-online.org/) output.

Full documentation: [Documentation/documentation.htm](Documentation/documentation.htm)  
Manual: [Documentation/end.rapid.Manual.htm](Documentation/end.rapid.Manual.htm)

## Requirements

- [CCP4](http://www.ccp4.ac.uk/) (tested with 6.x – 9.x)
- [Phenix](http://www.phenix-online.org/) (tested with 1.6 – 2.1)
- tcsh

## Installation

```bash
tar -zxvf end.rapid.tar.gz
# add end.rapid/ to your PATH
```

## Usage

```bash
END_RAPID.com phenixrefine.eff [seeds=5] [cycles=5] [cpus=N] [-nofofc] [-nosigf] [-norapid]
```

The input eff file is the one written by `phenix.refine` at the end of your
last refinement run.  No pre-processing of the eff file is needed —
`END_RAPID.com` automatically strips any parameter blocks unrecognized by the
current Phenix version (e.g. the `ddr { }` block introduced in Phenix 2.x).

### Output maps

| File | Description |
|---|---|
| `2FoFc_END.map` | Absolute-scale END map (vacuum = 0 e⁻/Å³) |
| `2FoFc_error.map` | RAPID noise from \|Fo−Fc\| perturbation |
| `2FoFc_sigF_error.map` | RAPID noise from σ(Fobs) perturbation |
| `2FoFc_snr.map` | Voxel-wise signal-to-noise ratio |
| `FoFc_scaled.map` | Fo−Fc difference map on absolute scale |
| `FoFc_error.map` | RAPID noise for the Fo−Fc map |

## Citation

Lang PT, Holton JM, Fraser JS, Alber T.
*Protein structural ensembles are revealed by redefining x-ray electron density noise.*
PNAS USA **111**, 237–247 (2014).
