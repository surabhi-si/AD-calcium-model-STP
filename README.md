# AD Presynaptic Calcium Model

Code and processed-data archive for **"Reassessing the Calcium Basis of Short-Term Plasticity
Deficit in Alzheimer's Disease."** Models a single hippocampal CA3→CA1 presynaptic bouton with
stochastic MCell simulation (VDCCs, RyR-mediated ER calcium release, calbindin/SERCA/PMCA
buffering, vesicle fusion), validated against a postsynaptic CA1 pyramidal cell model in NEURON.

## Repository structure

```
presynaptic_model/     MCell model definition files (.mdl) -- the presynaptic bouton itself
postsynaptic_model/    NEURON validation of postsynaptic EPSCs (vendored ModelDB #20212)
analysis/               Raw MCell output -> processed .dat conversion pipeline
figures/                One notebook per paper figure; equation reference PDF; small processed
                         data that feeds multiple figures
```

Each of the four has its own README with directory layout, naming conventions, and setup
details. This file is the map between them and the paper.

## Figure → code map

| Figure | Notebook | Contents |
|---|---|---|
| 3 | [`figures/fig3_single_AP.ipynb`](figures/fig3_single_AP.ipynb) | Single-AP presynaptic dynamics: VDCC/RyR flux, cytosolic/AZ/ER calcium, calbindin/SERCA buffering, PMCA extrusion, vesicle release timing |
| 4 | [`figures/fig4_PSblocked_paired_pulse.ipynb`](figures/fig4_PSblocked_paired_pulse.ipynb) | Paired-pulse (AP1/AP2) traces: Control vs presenilin(PS)-blocked |
| 5 | [`figures/fig5_RyR_paired_pulse_facilitation.ipynb`](figures/fig5_RyR_paired_pulse_facilitation.ipynb) | RyR-dependent paired-pulse facilitation across conditions |
| 6 | [`figures/fig6_RyR_train.ipynb`](figures/fig6_RyR_train.ipynb) | RyR flux and ER depletion across a 20-pulse, 20 Hz train |
| 7 | [`figures/fig7_ER_hyperactivity.ipynb`](figures/fig7_ER_hyperactivity.ipynb) | ER hyperactivity, 20-pulse train summary (VDCC = 80) |
| 8 | [`figures/fig8_poirazi_synapse_validation.ipynb`](figures/fig8_poirazi_synapse_validation.ipynb) | Postsynaptic CA1 synapse calibration against Smith, Ellis-Davies & Magee (2003) |
| 9 | [`figures/fig9_postsynaptic_train_EPSC.ipynb`](figures/fig9_postsynaptic_train_EPSC.ipynb) | Postsynaptic somatic EPSC across the train, 4 conditions |
| 10 | [`figures/fig10_sustained_unreliability.ipynb`](figures/fig10_sustained_unreliability.ipynb) | AD synapse reliability breakdown under sustained stimulation |
| 11 | [`figures/fig11_sync_async_reliability.ipynb`](figures/fig11_sync_async_reliability.ipynb) | Synchronous vs. asynchronous release tracking reliability differently |
| S1 | [`figures/figS1_model_fits_timescales.ipynb`](figures/figS1_model_fits_timescales.ipynb) | Kinetics/timescale model fits underlying Figure 3's rise/decay numbers |
| S2 | [`figures/figS2_ER_steadystate.ipynb`](figures/figS2_ER_steadystate.ipynb) | Confirms resting ER calcium is at steady state pre-stimulus |
| S3 | [`figures/figS3_ER_hyperactivity_RyRunder250.ipynb`](figures/figS3_ER_hyperactivity_RyRunder250.ipynb) | RyR under-expression without ER calcium overload |
| S4 | [`figures/figS4_paired_pulse_ISI_sync_coupling.ipynb`](figures/figS4_paired_pulse_ISI_sync_coupling.ipynb) | Paired-pulse ISI dependence and synchronous-release calcium coupling |

Every equation computed inline across these notebooks is documented in
[`figures/equations_reference.pdf`](figures/equations_reference.pdf). Equations in the
raw-to-processed conversion pipeline (`analysis/modFunc.py`) are documented separately in
[`analysis/conversion_equations.pdf`](analysis/conversion_equations.pdf).

## Quickstart

```bash
conda env create -f environment.yml
conda activate ad-presynaptic-calcium-model
```

1. **Run the presynaptic model** (MCell) -- see [`presynaptic_model/README.md`](presynaptic_model/README.md)
   for directory layout, naming conventions, and the exact `mcell` invocation (path resolution
   is working-directory-based, not relative to the included file).
2. **Convert raw output to processed data** -- see [`analysis/README.md`](analysis/README.md)
   for the `.dat` column reference and the PBS array-job pipeline (`analysisCluster.py`).
3. **Validate against the postsynaptic model** (NEURON) -- see
   [`postsynaptic_model/README.md`](postsynaptic_model/README.md); requires compiling the
   vendored mechanisms once (`nrnivmodl`).
4. **Regenerate a figure** -- open the corresponding notebook in `figures/`. Most read directly
   from `results/` directories produced by steps 1-2, which are **not included in this repo**
   (large, per-condition/per-seed simulation output) -- each notebook's own `CONFIG` cell states
   the exact paths it expects. `figures/paired_pulse_summary/` and the two
   `figS1_model_fits_timescales_*` files are small, pre-aggregated exceptions that **are**
   included, since they're already small enough to version directly.

## Data availability

Raw MCell/NEURON simulation output is not included (too large, and per-trial/per-condition).
What's included instead:
- All model definition, analysis, and figure-generation **code**
- The small, pre-aggregated summary files each figure notebook's `CONFIG` cell names explicitly
- Rendered figure outputs (`figures/outputs/`)
- Full equation references for both the analysis pipeline and the figure notebooks

Each of `presynaptic_model/`, `postsynaptic_model/`, and `analysis/` documents exactly what its
un-vendored inputs/outputs are (directory layout, file format) in its own README, so a reader
with access to the underlying simulation output can reproduce any figure end to end.

## Citation

See [`CITATION.cff`](CITATION.cff).
