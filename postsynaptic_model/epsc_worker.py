"""Shared configuration + parallel worker for the CA3->CA1 EPSC pipeline.

Used by both ``poirazi_EPSC_validation.ipynb`` (single unitary release, calibrated
against published measurements) and ``poirazi_EPSC_production.ipynb`` (full release
trains, three conditions) -- in this repo, ``figures/fig8_poirazi_synapse_validation.ipynb``
and ``figures/fig9_postsynaptic_train_EPSC.ipynb``.

``DEFAULT_CFG`` is the single source of truth for the synapse and electrode, so the
validation notebook and the production notebook cannot silently drift apart.

Kept as a module rather than notebook cells because the ``spawn`` start method
cannot pickle functions defined inside a notebook.

IMPORTANT -- never launch more than ``MAX_WORKERS`` processes; see the note on that
constant for why the number is what it is.

Vendored into this repo 2026-09-03 from the original working copy at
``Paper_figure_scripts/epsc_worker.py`` -- only ``CA1_ROOT`` changed (now
repo-relative, pointing at ./mechanisms). ``OUT_DIR``/``RELEASE_DIRS`` are
deliberately left as machine-specific absolute paths -- see the README in this
folder for what they need to contain and why they aren't vendored too.
"""

from __future__ import annotations

import os
from contextlib import contextmanager
from pathlib import Path

import numpy as np

# Cap on concurrent worker processes. Provenance, so it is not over-trusted:
#
#   2026-08-14  user: "use only 90 cores at a time, otherwise the system hangs."
#               Workload, duration and mechanism were not recorded.
#   2026-08-17  a sustained NEURON sweep took the PSU down (the notebook default
#               was 60 workers); three reboots that afternoon.
#
# Whether those are two failure modes or one is NOT established. A power supply
# sagging under sustained load can present as a freeze, so the 90 figure may
# never have been safe -- it may simply be where the symptom was first noticed,
# on a shorter or lighter job. 30 is the only number with a clean observation
# behind it. Do not raise this on the strength of the 90 figure.
MAX_WORKERS = 30

SAFE_WORKERS = MAX_WORKERS  # kept so existing imports keep working

# Vendored ModelDB #20212 (Poirazi, Brannon & Mel 2003) subset, next to this file.
# See postsynaptic_model/README.md for what's included/excluded and why.
CA1_ROOT = Path(__file__).resolve().parent / "mechanisms"
MECH_DIR = CA1_ROOT / "mechanism"
MODEL_FILE = "../experiment/spike-train-attenuation/bpap_1.hoc"

# -----------------------------------------------------------------------------
# NOT vendored into the repo -- machine-specific, large, per-trial data.
# Re-point these at wherever you keep them before running anything in this
# folder. See postsynaptic_model/README.md ("Input and output data") for the
# exact directory layout and file format expected at each of these.
# -----------------------------------------------------------------------------
OUT_DIR = Path("/data/hdrive/fEPSC_results/epsc_corrected")

BASE_RELEASE_DIR = Path("/data/hdrive/fEPSC_results/train_single_bouton")
RELEASE_DIRS = {
    "Control": BASE_RELEASE_DIR / "cont_reltimes_single_bouton",
    "AD":      BASE_RELEASE_DIR / "AD_reltimes_single_bouton",
    "CB":      BASE_RELEASE_DIR / "CB_reltimes_single_bouton",
    "ADB":      BASE_RELEASE_DIR / "ADB_reltimes_single_bouton",
}
# condition key -> plot_config.CONDITIONS key
PLOT_KEY = {"Control": "control", "AD": "ad", "CB": "control_er_blocked", "ADB": "ad_er_blocked"}

# ---------------------------------------------------------------------------
# Single source of truth for the model
# ---------------------------------------------------------------------------
_PRE_PAD, _ISI, _N_PULSES, _POST_PAD = 100.0, 50.0, 20, 250.0

DEFAULT_CFG = dict(
    # numerics
    dt=0.025,                 # ms
    celsius=34.0,             # C, matches the Poirazi cell-setup default

    # synapse -- POINT_PROCESS, absolute conductance, NOT a density
    syn_sec_index=9,          # h.apical_dendrite[9], ~249 um from soma, stratum radiatum
    syn_seg_x=0.5,
    g_unitary_pS=383.0,       # 26.8 pA / 70 mV, Smith, Ellis-Davies & Magee 2003 (240-280 um)
    tau_rise=0.2,             # ms, their rise tau 169-184 us
    tau_decay=4.4,            # ms. Smith et al. report 4.4 ms for their DISTAL group
                              # (240-280 um) and 5.9 ms proximal; apical_dendrite[9] is at
                              # 249 um, so the distal value is the matched one. Was 5.0 ms
                              # (the mid-range value) until 2026-08-18 -- results computed
                              # before that date are NOT comparable to results after it.
    e_ampa=0.0,               # mV

    # recording -- ONE somatic electrode with realistic series resistance
    v_hold=-70.0,             # mV; isolates AMPA (NMDA largely Mg-blocked)
    soma_sec=1,               # h.soma[1], the 562 um^2 main somatic compartment
    series_r=10.0,            # MOhm; None -> ideal h.VClamp (fires a 70 nA transient at t=0)

    # trial timing
    pre_pad=_PRE_PAD,         # quiet baseline before pulse 1
    isi=_ISI,                 # ms (20 Hz)
    n_pulses=_N_PULSES,
    post_pad=_POST_PAD,       # tail so pulse 20 decays fully
    tstop=_PRE_PAD + _N_PULSES * _ISI + _POST_PAD,

    # data
    release_unit="s",         # release-time files are in SECONDS
    trace_decimate=8,         # store traces at dt*8 = 0.2 ms to keep memory sane
)


def load_release_times(path, unit="s") -> np.ndarray:
    """Release times, sorted, in ms. Files on disk are unsorted and in seconds."""
    x = np.atleast_1d(np.loadtxt(path)).ravel()
    if unit == "s":
        x = x * 1000.0
    return np.sort(x)


def list_trial_files(condition, bouton, root=None):
    root = (root or RELEASE_DIRS[condition]) / bouton
    dirs = sorted(d for d in root.iterdir() if d.is_dir() and d.name.startswith("trial"))
    return [str(next(d.glob("*.txt"))) for d in dirs]


# ---------------------------------------------------------------------------
# Worker
# ---------------------------------------------------------------------------
_W: dict = {}


def _mute_stdout() -> None:
    """Permanently point this worker's stdout at /dev/null.

    Loading the CA1 model prints ~45 lines of hoc chatter; with dozens of
    workers that buries the notebook output. Neither ``h.hoc_stdout`` nor a
    scoped fd-1 redirect suppresses it -- NEURON buffers the text in C and
    flushes it after the scope closes -- so the redirect has to be permanent.

    Verified safe: all of that output goes to stdout and none to stderr, and a
    worker never communicates through stdout anyway. Results come back as
    return values and exceptions propagate through the pool, so real failures
    are still visible.
    """
    devnull = os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull, 1)
    os.close(devnull)


def init_worker(cfg: dict) -> None:
    """Load NEURON + the CA1 cell, build one synapse and one somatic electrode."""
    os.chdir(MECH_DIR)
    os.environ["NEURON_MODULE_OPTIONS"] = "-nogui -NSTACK 100000 -NFRAME 20000"
    _mute_stdout()

    from neuron import h

    assert h.load_file(MODEL_FILE), "Failed to load the Poirazi CA1 model."
    h.dt = cfg["dt"]
    h.celsius = cfg["celsius"]

    sec = h.apical_dendrite[cfg["syn_sec_index"]]
    syn = h.Exp2Syn(sec(cfg["syn_seg_x"]))
    syn.tau1, syn.tau2, syn.e = cfg["tau_rise"], cfg["tau_decay"], cfg["e_ampa"]

    # Exp2Syn is normalised so that peak g of a single event == weight (in uS)
    nc = h.NetCon(None, syn)
    nc.weight[0] = cfg["g_unitary_pS"] * 1e-6
    nc.delay = 0.0

    soma = h.soma[cfg["soma_sec"]]
    if cfg["series_r"] is None:
        clamp = h.VClamp(soma(0.5))
        clamp.dur[0], clamp.dur[1], clamp.dur[2] = cfg["tstop"], 0, 0
        for j in range(3):
            clamp.amp[j] = cfg["v_hold"]
    else:
        clamp = h.SEClamp(soma(0.5))
        clamp.dur1, clamp.amp1, clamp.rs = cfg["tstop"], cfg["v_hold"], cfg["series_r"]
        clamp.dur2 = clamp.dur3 = 0

    _W.update(h=h, sec=sec, syn=syn, nc=nc, clamp=clamp, cfg=cfg)


def simulate(times_ms, cfg=None):
    """Run one trial with the given vesicle release times (ms, already sorted).

    Returns the somatic clamp current in pA as a fresh float array.
    """
    cfg = cfg or _W["cfg"]
    h, nc, clamp = _W["h"], _W["nc"], _W["clamp"]

    times = np.asarray(times_ms, dtype=float)

    def _arm():
        for t in times:
            nc.event(float(t))

    _W["fih"] = h.FInitializeHandler(_arm)   # keep the reference alive for the run
    rec = h.Vector().record(clamp._ref_i)

    h.finitialize(cfg["v_hold"])
    h.continuerun(cfg["tstop"])

    # np.asarray() on a NEURON Vector gives a READ-ONLY view into NEURON's memory.
    # np.array() copies, which is what we need before any arithmetic.
    return np.array(rec, dtype=float) * 1e3   # nA -> pA


def per_pulse_amplitudes(i_pA, cfg=None, baseline_ms=5.0):
    """Peak EPSC of each pulse, baselined to the 5 ms immediately before it.

    Returned positive-going (an inward current gives a positive amplitude), which
    is how EPSC amplitudes are conventionally tabulated.
    """
    cfg = cfg or _W["cfg"]
    dt, pre, isi, npulse = cfg["dt"], cfg["pre_pad"], cfg["isi"], cfg["n_pulses"]
    amps = np.empty(npulse)
    for k in range(npulse):
        te = pre + k * isi
        base = i_pA[int((te - baseline_ms) / dt):int(te / dt)].mean()
        seg = i_pA[int(te / dt):int((te + isi) / dt)] - base
        amps[k] = -seg.min() if len(seg) else np.nan
    return amps


def run_trial_file(path):
    """Replay one trial file.

    Returns (per-pulse amplitudes pA, vesicles per pulse, decimated trace pA).
    """
    cfg = _W["cfg"]
    times = load_release_times(path, cfg["release_unit"])

    vesicles = np.zeros(cfg["n_pulses"])
    idx = np.clip((times // cfg["isi"]).astype(int), 0, cfg["n_pulses"] - 1)
    for j in idx:
        vesicles[j] += 1

    i_pA = simulate(times + cfg["pre_pad"], cfg)
    amps = per_pulse_amplitudes(i_pA, cfg)
    trace = i_pA[::cfg["trace_decimate"]].astype(np.float32)
    return amps, vesicles, trace


# ---------------------------------------------------------------------------
# Checkpointed batch runner
#
# Lives here rather than in the notebook so that the notebook and the headless
# CLI (run_epsc.py) execute the same code. A long run on this machine will
# probably be interrupted at some point, so every path through this is
# resumable.
# ---------------------------------------------------------------------------
def cache_path(condition, bouton, out_dir=None):
    """One file per condition per bouton.

    The trial count is deliberately NOT in the filename. It used to be, which
    meant raising n_trials from 100 to 1000 looked for a different file, found
    nothing, and silently re-simulated the 100 trials already on disk. With a
    stable name, raising n_trials simply tops the existing file up.
    """
    return (out_dir or OUT_DIR) / f"{bouton}_{condition}.npz"


# Config keys that change the simulated result. A checkpoint computed under
# different values of these is not comparable and must not be silently reused --
# the resumable cache would otherwise report "already complete" and hand back
# stale data after a parameter change.
_FINGERPRINT_KEYS = ("dt", "celsius", "syn_sec_index", "syn_seg_x", "g_unitary_pS",
                     "tau_rise", "tau_decay", "e_ampa", "v_hold", "soma_sec",
                     "series_r", "pre_pad", "isi", "n_pulses", "tstop",
                     "release_unit", "trace_decimate")


def _fingerprint(cfg):
    cfg = cfg or DEFAULT_CFG
    return np.array([f"{k}={cfg.get(k)}" for k in _FINGERPRINT_KEYS], dtype="U64")


def load_partial(condition, bouton, out_dir=None, log=print, cfg=None):
    """Return (A, V, T) as lists, empty if there is no usable checkpoint.

    A checkpoint whose stored configuration differs from `cfg` is discarded, so
    changing a synapse parameter can never silently resume onto stale trials.
    """
    d = out_dir or OUT_DIR
    p = cache_path(condition, bouton, d)
    if not p.exists():
        # Fall back to the old "<bouton>_<cond>_n<N>.npz" layout so runs made
        # before the rename are not stranded; take the one with most trials.
        legacy = sorted(d.glob(f"{bouton}_{condition}_n*.npz"))
        if not legacy:
            return [], [], []
        p = max(legacy, key=lambda q: q.stat().st_size)
        log(f"  reading legacy checkpoint {p.name}")
    try:
        z = np.load(p)
        if "cfg" in z.files:
            want, got = _fingerprint(cfg), z["cfg"]
            if list(want) != list(got):
                diff = [f"{w.split('=')[0]}: {g.split('=')[1]} -> {w.split('=')[1]}"
                        for w, g in zip(want, got) if w != g]
                log(f"  DISCARDING {p.name}: it was computed with a different "
                    f"configuration ({'; '.join(diff)}). Recomputing from scratch.")
                return [], [], []
        else:
            log(f"  NOTE {p.name} predates configuration stamping -- its parameters "
                f"cannot be verified. Delete it if it may be stale.")
        return list(z["A"]), list(z["V"]), list(z["T"])
    except Exception as e:                       # truncated by a power cut
        log(f"  ignoring unreadable checkpoint {p.name}: {type(e).__name__}")
        return [], [], []


def save_partial(condition, bouton, A, V, T, out_dir=None, cfg=None):
    p = cache_path(condition, bouton, out_dir)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".npz.tmp")
    # Write to a temp file and rename, so a power cut can never leave a
    # half-written .npz behind. The open handle matters: np.savez_compressed
    # appends ".npz" to any *path* lacking it, which would silently produce
    # "....npz.tmp.npz" and break the rename.
    with open(tmp, "wb") as fh:
        np.savez_compressed(fh, A=np.array(A), V=np.array(V), T=np.array(T),
                            cfg=_fingerprint(cfg))
    tmp.replace(p)


def run_condition(condition, bouton="V80", n_trials=1000, n_workers=1,
                  batch=10, out_dir=None, cfg=None, cooldown=0.0, log=print):
    """Simulate one condition, checkpointing every `batch` trials.

    Re-running resumes from the last checkpoint. `n_workers` is capped at
    MAX_WORKERS -- see the note on that constant before raising it.
    """
    import multiprocessing as mp
    import time

    if n_workers > MAX_WORKERS:
        raise ValueError(
            f"{n_workers} workers exceeds MAX_WORKERS={MAX_WORKERS}; "
            "a sustained run above this took the PSU down.")

    cfg = dict(cfg or DEFAULT_CFG)
    files = list_trial_files(condition, bouton)[:n_trials]
    A, V, T = load_partial(condition, bouton, out_dir, log, cfg)
    done = len(A)

    if done >= len(files):
        log(f"{condition}: already complete ({done} trials)")
        return dict(A=np.array(A), V=np.array(V), T=np.array(T))

    todo = files[done:]
    eta = len(todo) * 17.6 / n_workers
    log(f"{condition}/{bouton}: {done} done, {len(todo)} to go, "
        f"{n_workers} worker(s), batch={batch}  (~{eta / 60:.0f} min)")

    t0 = time.time()
    ctx = mp.get_context("spawn")
    with ctx.Pool(n_workers, initializer=init_worker, initargs=(cfg,)) as pool:
        for k, (a, v, t) in enumerate(pool.imap(run_trial_file, todo, chunksize=1), 1):
            A.append(a); V.append(v); T.append(t)
            if k % batch == 0 or k == len(todo):
                save_partial(condition, bouton, A, V, T, out_dir, cfg)
                el = time.time() - t0
                log(f"  {len(A):5d}/{len(files)}  {el / 60:5.1f} min elapsed, "
                    f"~{(el / k) * (len(todo) - k) / 60:5.1f} min left "
                    f"[checkpointed]")
                if cooldown and k < len(todo):
                    time.sleep(cooldown)

    log(f"{condition}: done in {(time.time() - t0) / 60:.1f} min -> "
        f"{cache_path(condition, bouton, out_dir).name}")
    return dict(A=np.array(A), V=np.array(V), T=np.array(T))
