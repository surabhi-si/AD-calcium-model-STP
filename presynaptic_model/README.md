# Presynaptic MCell Bouton Model

MCell (MDL) model of Ca²⁺ dynamics and neurotransmitter release at a presynaptic CA3
bouton, used to compare control and Alzheimer's-disease-associated RyR/ER
conditions across three stimulus protocols.

## Directory structure

```
presynaptic_model/
├── core/                          # files shared by every condition and protocol
└── per_condition_per_protocol/
    ├── vdcc_dat/                  # VDCC rate-constant .dat files, shared across all conditions
    ├── vdcc_disks/                # VDCC release-site disk geometry, keyed by channel count
    ├── ryr_disks/                 # RyR release-site disk geometry, keyed by channel count
    └── <condition>/
        ├── <condition-specific reaction/parameter files>
        ├── Train/<mainfile>.mdl
        ├── paired-pulse/<mainfile>.mdl
        └── singleAP/<mainfile>.mdl       
```

`core/` holds the geometry, channel/pump/vesicle-release reaction definitions,
and output-reporting logic common to the whole model: `AZ.mdl`, `calbindin.mdl`,
`initialization.mdl`, `instantiate.mdl`, `outputLoc.mdl`, `PMCA.mdl`,
`rxn_outputRS.mdl`, `rxn_rate_constants.mdl`, `surface_classes.mdl`,
`transmitter_trigger_complex_{mol,rxn}_15d7_sudhof_v2.mdl`,
`tripartite_synapse_v7_rrp_15_er_vdcc_junction.geometry.mdl`, `VDCC.mdl`,
`viz_output.mdl`.

Each `<condition>/` folder holds exactly the files that define that disease
state: the RyR reaction-variant file, `SERCA_ss_*.mdl`, `ERleak_*.mdl`, and
`misc_*.mdl`. Everything else needed for a run — the geometry, reaction
network, and stimulus-specific VDCC rate data — is pulled in from `core/` and
the two `per_condition_per_protocol/` root folders (`vdcc_dat/`, `vdcc_disks/`,
`ryr_disks/`).

## Conditions

| Folder | Description | ER [Ca²⁺]|
|---|---|---|
| `Control` | Baseline model | 250 µM |
| `Control-ERblocked` | Control with ER Ca²⁺ stores blocked (SERCA off, RyR off, fixed ER leak) | 250 µM |
| `PSblocked` | Presenilin blocked (ER leak zeroed; distinct mechanism from ER-blocked) | 250 µM |
| `AD-1_2xRyR` | 1/2× RyR overexpression | 750 µM |
| `AD-1_3xRyR` | 1/3× RyR overexpression | 750 µM |
| `AD-3xRyR` | 3× RyR overexpression | 750 µM |
| `AD-5xRyR` | 5× RyR overexpression | 750 µM |
| `AD-ERblocked` | 5× RyR overexpression with ER Ca²⁺ stores blocked | 750 µM |

All conditions currently use 80 P/Q-type VDCC channels (`VDCC_num = 80`).

## Protocols

| Folder | Stimulus | Iterations |
|---|---|---|
| `Train` | 20 pulses at 20 Hz | 1,000,000 |
| `paired-pulse` | 2 pulses, 40 ms inter-stimulus interval | 70,000 |
| `singleAP` | 1 pulse | 50,000 |

## Naming conventions

**Main `.mdl` filenames** (one per condition × protocol):
- Train: `V{vdcc}C{dose}uM_{n}_{freq}hz[_<condition suffix>].mdl`, e.g. `V80C750uM_20_20hz_5xryr.mdl`
- Paired-pulse: `I{isi}V{vdcc}[_<condition suffix>].mdl`, e.g. `I40V80_5xryr_ERblocked.mdl`
- Single AP: `V{vdcc}_single_AP.mdl`

**Output `fname`/`output_folder`** (set inside each main file, drives where results land — see [analysis README] for the resulting directory schema): distinct per condition and protocol, e.g. Train → `train/AD_5xRyR/RSV80C750uM_20_p20hz_5xryr/s_<seed>/`, paired-pulse → `prvsppf_modeltesting/AD_5xRyR/RSI40V80C750uM_5xryr/s_<seed>/`.


## Path resolution

MCell resolves every relative path — `INCLUDE_FILE`, reaction-rate `.dat`
references, disk-geometry includes — relative to the **working directory
`mcell` is invoked from**, not relative to the file containing the reference.
The required convention is to run from inside the protocol folder:

```bash
cd per_condition_per_protocol/<condition>/<protocol>/
mcell <mainfile>.mdl -seed <N>
```

From that working directory, main files reach `core/` via `../../../core/`,
the shared `vdcc_dat/`/`vdcc_disks/`/`ryr_disks/` via `../../<folder>/`, and
their own condition folder via `../`. If this folder is ever moved or
renamed, every relative path across all `core/` and `per_condition_per_protocol/`
files needs re-deriving from the new location and re-verified against this
same CWD convention — don't assume the existing depths still hold.

## Requirements

MCell (version 3.2).
