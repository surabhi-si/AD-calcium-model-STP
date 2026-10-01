# Postsynaptic CA1 Model

NEURON-based validation of the presynaptic MCell model's predictions at the
postsynaptic side of the CA3→CA1 Schaffer-collateral synapse. Two stages,
corresponding to Figures 8 and 9 of the paper:

1. **Unitary synapse validation** (Fig. 8) — a single vesicle release calibrated
   against dendritic patch-clamp recordings (Smith, Ellis-Davies & Magee 2003).
2. **Train replay across conditions** (Fig. 9) — MCell-simulated single-bouton
   release trains (20 pulses, 20 Hz) for Control, AD, Control+ER-blocked, and
   AD+ER-blocked, replayed through the same calibrated synapse to produce
   somatic EPSCs.

## Attribution

`mechanisms/` is a **modified, trimmed subset** of ModelDB accession
[20212](https://modeldb.science/20212) — Poirazi P, Brannon T, Mel BW (2003)
*Arithmetic of subthreshold synaptic summation in a model CA1 pyramidal cell*,
Neuron 37:977-987, and the companion paper *Pyramidal Neuron as Two-Layer
Neural Network*, Neuron 37:989-999. Original model, full ModelDB entry, and
online supplement at the link above.


**Placements:** an `Exp2Syn` point process (absolute conductance in pS, not a
density), one `NetCon` event per vesicle, on `apical_dendrite[9]` (stratum
radiatum, ~249 µm from soma), read out by a single somatic `SEClamp`.

## What's in `mechanisms/`, and what's deliberately left out

Verified empirically (`strace`-traced a real run using cached results, then a
full compile + live single-trial test of this exact copy — not just read from
the original `.hoc` files) rather than guessed:

**Included**:
- `mechanism/*.mod` — 17 mechanism definitions (compile with `nrnivmodl`,
  see Setup below)
- `lib/`, `template/` — whole directories; several files here are pulled in
  generically by the model's `ExperimentControl` setup rather than named
  directly in `bpap_1.hoc`
- `morphology/n123/` — whole directory (cell geometry + compartment lists)
- `experiment/cell-setup.hoc` and `experiment/spike-train-attenuation/bpap_1.hoc`
  only — the specific experiment entry point used here


## Setup

```bash
cd postsynaptic_model/mechanisms/mechanism
nrnivmodl .
```

This creates a local `x86_64/` (or platform-equivalent) directory that NEURON
auto-loads. 

Requires the `neuron` Python package (`pip install neuron` or via conda).

## Input and output data (not included in this repo)

Both are large, per-trial, machine-specific, and deliberately **not** vendored
here — `epsc_worker.py`'s `RELEASE_DIRS`/`OUT_DIR` constants point at them and
need to be re-pointed at wherever you keep them.

**Input** — vesicle release-time files, the bridge from the presynaptic MCell
model to this postsynaptic one. Directory layout:
```
<BASE_RELEASE_DIR>/
├── cont_reltimes_single_bouton/V80_reltimes/trial<N>/rel.txt
├── AD_reltimes_single_bouton/V80_reltimes/trial<N>/rel.txt
├── CB_reltimes_single_bouton/V80_reltimes/trial<N>/rel.txt
└── ADB_reltimes_single_bouton/V80_reltimes/trial<N>/rel.txt
```
One `rel.txt` per simulated trial, containing that trial's vesicle release
timestamps, **in seconds, unsorted** (`epsc_worker.load_release_times` sorts
and converts to ms). These are derived from the presynaptic MCell model's
`vdcc.sync_*.dat`/`vdcc.async_*.dat`/`vdcc.spont.dat` trigger outputs (see
`presynaptic_model/README.md` and `analysis/README.md`) for one condition each
— `cont`=Control, `AD`=AD, `CB`=Control+ER-blocked, `ADB`=AD+ER-blocked.

**Output** — one `.npz` checkpoint per condition, `<bouton>_<condition>.npz`
(e.g. `V80_reltimes_AD.npz`), containing:
- `A` — per-trial, per-pulse EPSC amplitude (pA), shape `(n_trials, 20)`
- `V` — per-trial, per-pulse vesicle count, same shape
- `T` — per-trial decimated somatic current trace (pA)
- `cfg` — a fingerprint of the synapse/electrode config used, so a checkpoint
  computed under different parameters is detected and discarded rather than
  silently reused (see `epsc_worker._FINGERPRINT_KEYS`)

Checkpointing means a run can be safely interrupted and resumed — `run_epsc.py`
and `epsc_worker.run_condition` both check existing files in `OUT_DIR` before
launching new simulations.

## How to run

Interactively — see `figures/fig8_poirazi_synapse_validation.ipynb` (run
first, calibrates the synapse) and `figures/fig9_postsynaptic_train_EPSC.ipynb`
(replays all 4 conditions; set `CONDITION` and re-run the cell for each one).

Headless, one condition at a time:
```bash
python run_epsc.py --condition AD --bouton V80_reltimes --trials 1000 --workers 1 --batch 10
```
`--workers` is hard-capped at `epsc_worker.MAX_WORKERS = 30` 

## Key configuration (`epsc_worker.DEFAULT_CFG`)

| Parameter | Value | Note |
|---|---|---|
| Synapse location | `apical_dendrite[9]`, ~249 µm from soma | stratum radiatum |
| Unitary conductance | 383 pS | fit to Smith et al. 2003's 26.8±2.0 pA distal dendritic EPSC |
| τ rise / decay | 0.2 / 4.4 ms | 4.4 ms is Smith et al.'s *distal* (240-280 µm) value |
| Holding potential | −70 mV | isolates AMPA (NMDA largely Mg²⁺-blocked) |
| Series resistance | 10 MΩ | somatic electrode; Smith et al. report 10-30 MΩ dendritic |
| Protocol | 20 pulses, 20 Hz | matches the presynaptic Train protocol |
