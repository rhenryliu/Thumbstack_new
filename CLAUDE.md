# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

ThumbStack produces stacked CMB maps and radial profiles for thermal and kinematic Sunyaev-Zel'dovich (tSZ/kSZ) measurements. Given a galaxy catalog and a CMB map, it outputs 2D stacked maps and radial profiles (aperture photometry filters at varying radii), along with covariance estimates via bootstrap resampling.

## Running scripts

Scripts live in `scripts/` and add `../src/` to `sys.path` to import library modules. They are run from the `scripts/` directory:

```bash
cd scripts/
python ThumbStack_AllProfiles.py          # direct run
python ThumbStack_NewCatalogues.py '{"catalog_name": "LRG_clustering.dat", "save": true}'  # JSON params
```

For SLURM jobs, use the existing shell scripts in `scripts/`:
- `runCPU.sh` — general CPU job (DESI cosmodesi environment)
- `runCPU_JSON.sh` — CPU job with JSON parameter passing (rhliu_tSZ conda env)
- `runCPU_debug.sh`, `runCPU2.sh`, `runCPU_Martine.sh`, `runCPU_Yifei.sh` — variants

Two environments are in use depending on the script:
- **DESI cosmodesi**: `source /global/common/software/desi/desi_environment.sh 23.1 && source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main`
- **Custom conda**: `conda activate rhliu_tSZ`

SLURM output logs go to `Outputs_Perlmutter/slurm-<jobid>.out`.

## Architecture

### Core library (`src/`)

All modules import shared dependencies via `from headers import *` (`headers.py` centralises numpy, scipy, matplotlib, pixell, healpy, pathos, sharedmem, pyfftw, classylss, etc.).

Key classes and their roles:

| Module | Class/function | Role |
|---|---|---|
| `thumbstack.py` | `ThumbStack` | Main stacking engine. Extracts cutouts, applies aperture filters, computes stacked profiles and covariances. |
| `cmbMap.py` | `cmbMap` | Thin wrapper around a pixell enmap file. Provides `.map()`, `.mask()`, `.hit()` methods; handles tSZ unit conversion from muK to Compton-y. |
| `catalog.py` | `Catalog`, helper functions | Galaxy catalog as a pandas DataFrame. Helpers: `addHaloMass`, `addIntegratedTau`, etc. |
| `universe.py` | `Universe`, `UnivMariana` | Wraps classylss/CLASS for cosmology (distances, growth, power spectra). |
| `mass_conversion.py` | `MassConversionKravtsov14` | Stellar-to-halo mass conversion (Kravtsov 2014 relation). |
| `computeProfiles.py` | `computeProfiles()` | Standalone profile computation extracted from ThumbStack, callable separately. |
| `flat_map.py` | `FlatMap` | Flat-sky 2D map with FFT operations. |
| `cmb.py` | `CMB` | CMB power spectrum templates (reads CAMB output from `input/universe_Planck15/camb/`). |
| `basic_functions.py` | utilities | Misc math/astro helpers. |
| `headers.py` | imports | Shared imports for all `src/` modules. |

### Typical analysis script pattern

```python
sys.path.append('../src/')
from universe import *
from mass_conversion import *
from catalog import *
from thumbstack import *
from cmbMap import *
from computeProfiles import computeProfiles

u = UnivMariana()
massConversion = MassConversionKravtsov14()
cat = Catalog(u, massConversion, name="...", save=False)
cmap = cmbMap(pathMap, pathMask=pathMask, pathHit=None, nu=93.e9, unitLatex=r'y', convert_y=False, name="...")
ts = ThumbStack(u, cat, cmap.map(), cmap.mask(), cmap.hit(),
                name="...", save=True, nProc=64,
                filterTypes='diskring', doBootstrap=True, rApMinArcmin=1.)
# Results: ts.stackedProfile[filterType+"_"+est], ts.sStackedProfile[...], ts.RApArcmin
```

### Filter types (`filterTypes` parameter)

`'diskring'`, `'ringring'`, `'ringring2'`, `'ringring3'`, `'disk'`, `'ring'`, `'all'`

### Estimators

- `'tsz_uniformweight'` — always available
- `'tsz_varweight'` — only when a hit map is provided

### Parallelism

`nProc` controls the number of parallel processes. Uses `pathos.multiprocessing.ProcessingPool` and `sharedmem` (Yu Feng's fork-based shared memory). On a Perlmutter CPU node, `nProc=64` (one full node = 256 hardware threads, but 64 tasks is typical).

### Data paths

Input CMB maps and output profiles are hardcoded to `/pscratch/sd/r/rhliu/projects/ThumbStack/`. This is known technical debt — do not refactor without explicit instruction.

### Catalog preprocessing

`data_cuts/DESI_ThumbStack_cuts.py` handles cuts on DESI LRG catalogs (nobs, EBV, star density, DEC, RA masks) before passing to `Catalog`.

### Older/reference code in `src/`

`driver_*.py`, `generate_*.py`, `coadd_pact_cmb.py`, `reconvolve_*.py` are legacy Cori-era scripts kept for reference. The `*.py.bak` files are old versions of active modules.
