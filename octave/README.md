# Running METIS with GNU Octave

METIS runs headless under [GNU Octave](https://octave.org) (tested with 8.4), so no MATLAB licence is needed to run simulations from scripts. The graphical interface is not supported.

## Quick start

```bash
sudo apt-get install octave octave-dev   # Debian/Ubuntu; octave-dev provides mkoctfile
octave/build_mex.sh                      # optional: compile the C interpolation routines (about 5x faster)
```

```matlab
addpath('octave'); metis_octave_setup;   % sets the path in the same order as MATLAB's metis_path
s = load('certification/metis/ITER_rampup_ECCD.mat');
inp = s.post.z0dinput; inp.option.machine = inp.machine;
[zs, info, profil0d] = zerod(inp.option, inp.cons, inp.geo, inp.exp0d);
```

## What was changed for Octave

Changes to METIS sources, all neutral for MATLAB:

- `zerod/z0j3.m`: a local variable named `do` (a reserved word in Octave) renamed `do_`.
- `import/divers/iscomplex.m`: returns false for non-numeric input. Octave's own library calls `iscomplex` on cell arrays.
- `import/sampling/tsample.c`: adds the missing `#include <float.h>`.

Octave-only files in `octave/` (never on the MATLAB path):

- `metis_octave_setup.m`: calls `zineb_path`, then puts `octave/compat` first and selects the pure-M PDE solver (the compiled one needs MPI/MUMPS).
- `compat/griddedInterpolant.m`, `compat/pwch.m`, `compat/chckxy.m`: minimal versions of MATLAB functions missing in Octave, used by `z0icrh.m` and `import/sampling/pchip.m`.
- `compat/time.m`: restores Octave's built-in `time()`, which `trait/time.m` shadows and which Octave's `clock()` needs.

## Certification against MATLAB

`run_certification_case.m` reruns a case from `certification/metis` the same way `metis_test.m` does and compares with the MATLAB result stored in the file (MATLAB R2019b, pure-M solver). The GitHub Actions workflow `.github/workflows/octave-certification.yml` runs all 28 cases on every push.

```bash
octave --no-gui --quiet octave/run_certification_case.m ITER_rampup_ECCD
```

The metric is the largest relative difference in the flattop means (Ip at least 98% of max) of `vloop`, `pfus`, `li`, `iboot`, `te0` and `betan`.

- **19 cases** reproduce MATLAB to machine precision (below 1e-7; tolerance 1e-6).
- **9 cases** differ by 0.02% to 9%, for a reason not yet identified. They are well conditioned (a 1e-12 input perturbation moves the result by at most 1e-6), so the difference comes from some built-in function behaving differently, not from chaos. They pass against per-case tolerances in `certification_tolerances.csv`, set at about twice the observed difference, so CI catches regressions.

| Case | Octave vs MATLAB |
|---|---|
| iter_2nbi_pfus | 5.3e-4 |
| reference_NTM_ITER | 7.6e-4 |
| JET53521_itb_test_nomorph | 1.7e-4 |
| DIII-D_like_shot_zjrun_test | 1.5e-3 |
| JET53521_itb_bgb | 4.9e-3 |
| JET53521_itb_test | 2.5e-2 |
| JET84746_density_4DoF | 2.7e-2 |
| TS28305_DDS_2 | 4.0e-2 |
| TFTR80539_like_neutron_rate_test | 9.1e-2 |

Ruled out so far: version drift (references regenerated in 2026 with the same solver), path order, `pchip`, Signal Toolbox replacements (`sgolayfilt`, `medfilt1`) and `interp1` default extrapolation. Dumping intermediate outputs from one MATLAB run of a diverging case would locate the remaining difference quickly.
