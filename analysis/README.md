# Analysis Pipeline

Turns raw per-trial MCell output (thousands of `s_<seed>/dat/*.dat` files per
condition) into averaged, per-condition derived quantities: concentration
traces, ion fluxes, paired-pulse/train release probabilities and facilitation.

## Runs on ada.local, not locally

Both scripts here are written to execute **on the PBS/Torque cluster
(ada.local)**. This repo holds a **versioned copy** of the scripts for review and
history, not something you run from this checkout. Concretely:

- `analysisCluster.py` is submitted directly via `qsub` (its shebang and
  `#PBS` directives are read by PBS itself) and hardcodes ada.local paths:
  `/home/surabhi/Output/model_testing/...` for raw input,
  `/home/surabhi/results/...` for output.

## Files

### `modFunc.py`

The reusable `class analysis`, instantiated once per condition/protocol
result directory (`analysis(dataPath, resultPath)`):

| Method | What it does |
|---|---|
| `getDirs` | Lists the `s_<seed>` trial subdirectories |
| `getData` | Manual `.dat` reader, skips `#` header lines |
| `avg_dat` | Averages one `.dat` file across all seeds → `resultPath` |
| `computeCaConc` | Step-windowed finite difference: `(time·Ca)[i+step] − (time·Ca)[i]) / dt` |
| `conc_calc` | Applies `computeCaConc` to `ca.dat` col 1 (`Ca.Conc.az`) → `CaConc` |
| `pre_ca_conc` | Direct count→concentration on col 2 (`Ca.Pre`) → `pre_ca_conc.dat`, **nM** |
| `er_ca_conc` | Direct count→concentration on col 3 (`Ca.ER`) → `er_ca_conc.dat`, **µM** |
| `relppf` | Paired-pulse `Pr1`/`Pr2` via 1000-resample bootstrap; AP times `(i·isi+2.0)/1000` s |
| `relptp` | Same idea for an n-pulse train; AP times `(i·isi)/1000` s (no 2ms offset) |
| `combineReleaseFiles` | Concatenates all seeds' `vdcc.sync_*`/`vdcc.async_*`/`rel.dat` |
| `fluxCurrent` | Windowed finite difference on a cumulative flux file |

### `analysisCluster.py`

A PBS array-job driver (`#PBS -t 1-11`). Each array task picks one function
out of `fl` by index (`fl[i]()`) and runs it against one result directory
(passed in via the `d` environment variable). `dataType` is currently
hardcoded to `'train/control/'` — edit this line before submitting for a
different condition or protocol.

Wrapper functions map each `.dat` file to the right `modFunc.py` call:
`azAvg`, `caAvg` (bundles `conc_calc`+`pre_ca_conc`+`er_ca_conc`),
`rrpAvg`, `vdccFlux`, `pmcaAvg`, `calbAvg`, `sercaFluxAvg`, `sercaMolAvg`,
`ryrFluxAvg`, `ryrMolAvg`, plus `ppf`/`rel`/`ptp` for release-probability
analysis (defined but not all wired into the current `fl`, which is
Train-protocol-oriented).

## `ca.dat` column reference

4 columns, header row `#Seconds Ca.Conc.az Ca.Pre Ca.ER`:

| Col | Name | Meaning | Handled by |
|---|---|---|---|
| 0 | `Seconds` | Time | — |
| 1 | `Ca.Conc.az` | MCell `ESTIMATE_CONCENTRATION` on the `active_zone_plane` region. This is a **cumulative average since t=0** (already µM), not an instantaneous reading — a raw value at time T is the mean over [0,T], not the value at T | `conc_calc` |
| 2 | `Ca.Pre` | Presynaptic/cytosolic Ca²⁺ molecule count | `pre_ca_conc`|
| 3 | `Ca.ER` | ER lumenal Ca²⁺ molecule count | `er_ca_conc` |

If `Ca.Conc.az(T)` is the mean concentration over `[0, T]`, then
`time × Ca.Conc.az ≈ ∫₀ᵀ c(t)dt`, the cumulative integral. A finite
difference of that integral between two nearby output steps —
`(∫₀^(T+Δ) c dt − ∫₀ᵀ c dt) / Δ` — gives the concentration averaged over just
that small window `Δ`, which is what `computeCaConc(time, ca, step)` computes
(`step` sets how many output rows wide that window is). `Ca.Pre`/`Ca.ER`
don't need this because they're plain volume `COUNT`s (instantaneous
molecule counts), not `ESTIMATE_CONCENTRATION` — they only need a
count→concentration unit conversion, not a cumulative→instantaneous one.

## All raw `.dat` files

Written by `core/rxn_outputRS.mdl` into each trial's `dat/` folder (see
[presynaptic_model/README.md](../presynaptic_model/README.md) for where those
trials land). 20 files total;
`vdcc.sync_*`/`vdcc.async_*` are grouped below since they're structurally
identical, one file per release pathway.

| File | Columns | Content | Handled by |
|---|---|---|---|
| `ca.dat` | 4 | `Ca.Conc.az`, `Ca.Pre`, `Ca.ER` — see table above | `conc_calc`, `pre_ca_conc`, `er_ca_conc` |
| `vdcc_pq_ca_flux.dat` | 2 | Cumulative Ca²⁺ flux through P/Q-type VDCCs | `avg_dat` → `fluxCurrent` (rate) |
| `vdcc.sync_{0,1,2}.dat` | event log | Synchronous vesicle-release events, one file per release pathway (0–2) | concatenated into `rel.dat` by `relppf`/`relptp`/`combineReleaseFiles` |
| `vdcc.async_{0..4}.dat` | event log | Asynchronous vesicle-release events, pathways 0–4 | same |
| `vdcc.spont.dat` | event log | Spontaneous vesicle-release events | same |
| `serca_mol.dat` | 7 | `SERCA_X1/X1A/X2/Y2/Y1A/Y1` — 6 pump conformational-state counts (X=cytosol-facing, Y=ER-facing; 1=apo, 1A=1 Ca²⁺, 2=2 Ca²⁺) | `avg_dat` only |
| `serca_ca_flux.dat` | 3 | `SERCA_Ca_out_flux`, `SERCA_Ca_in_flux` | `avg_dat` → `fluxCurrent` (rate) |
| `ryr_mol.dat` | 15 | 14 RyR channel-state counts across the `_L` and `_H1` conformational modes | `avg_dat` only |
| `ryr_ca_flux.dat` | 3 | `out`, `in` — summed RyR Ca²⁺ flux (across `_L`/`_H1` substates) | `avg_dat` → `fluxCurrent` (rate) |
| `az.dat` | 19 | `X0Y0`..`X5Y2` — 18-cell active-zone SNARE-assembly grid occupancy. **Not a calcium quantity** despite living next to `ca.dat` | `avg_dat` only |
| `pmca_mol.dat` | 4 | `PMCA_0/1/2` — pump state counts | `avg_dat` only |
| `pmca&leak_ca_flux.dat` | 4 | `PMCAFlux` (net), `PMCAleakFlux`, `ERleakFlux` | `avg_dat` only|
| `rrp.dat` | 2 | `RRP` — readily releasable pool size | `avg_dat` only |
| `calbindin_mol.dat` | 10 | 9 calbindin Ca²⁺-bound-state counts (`h0m0`..`h2m2`, H/M site occupancy combinations) | `avg_dat` only |

## Requirements

Python 3 with `numpy`/`scipy`; PBS/Torque `qsub` access on ada.local.
