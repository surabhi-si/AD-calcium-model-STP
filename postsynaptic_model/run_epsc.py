#!/usr/bin/env python3
"""Headless runner for one condition of the CA3->CA1 EPSC pipeline.

Runs in batches and checkpoints to disk after each one, so an interrupted run
resumes from the last checkpoint instead of starting over. Safe to re-run.

    python run_epsc.py --condition Control --trials 100 --workers 1 --batch 10

Deliberately a separate file from ``epsc_worker`` rather than a ``__main__``
block inside it: with the ``spawn`` start method, functions referenced from
``__main__`` are pickled by qualified name, and the pool initializer would
resolve differently in the parent and the children. Importing the module by
name avoids that entirely.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import epsc_worker as W  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--condition", required=True, choices=sorted(W.RELEASE_DIRS))
    ap.add_argument("--bouton", default="V92")
    ap.add_argument("--trials", type=int, default=100)
    ap.add_argument("--workers", type=int, default=1,
                    help=f"processes; hard capped at MAX_WORKERS={W.MAX_WORKERS}")
    ap.add_argument("--batch", type=int, default=10,
                    help="checkpoint to disk every N trials")
    ap.add_argument("--cooldown", type=float, default=0.0,
                    help="seconds to idle between batches")
    ap.add_argument("--out-dir", default=str(W.OUT_DIR))
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    def log(msg):
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

    log(f"condition={args.condition} bouton={args.bouton} trials={args.trials} "
        f"workers={args.workers} batch={args.batch}")
    log(f"synapse: apical_dendrite[{W.DEFAULT_CFG['syn_sec_index']}] "
        f"{W.DEFAULT_CFG['g_unitary_pS']:.0f} pS, "
        f"Vhold {W.DEFAULT_CFG['v_hold']} mV, Rs {W.DEFAULT_CFG['series_r']} MOhm")
    log(f"output: {W.cache_path(args.condition, args.bouton, out_dir)}")

    try:
        r = W.run_condition(args.condition, bouton=args.bouton, n_trials=args.trials,
                            n_workers=args.workers, batch=args.batch,
                            out_dir=out_dir, cooldown=args.cooldown, log=log)
    except ValueError as e:            # worker cap exceeded
        log(f"REFUSED: {e}")
        return 2
    except KeyboardInterrupt:
        log("interrupted -- progress up to the last checkpoint is saved; "
            "re-run this command to resume")
        return 130

    A, V = r["A"], r["V"]
    m = A.mean(axis=0)
    sem0 = A[:, 0].std() / max(len(A) ** 0.5, 1.0)

    log(f"n={len(A)}  mean vesicles/trial={V.sum(1).mean():.2f}")
    log(f"pulse 1: EPSC={m[0]:.2f} +/- {sem0:.2f} pA, "
        f"no release in {100 * (V[:, 0] == 0).mean():.0f} % of trials")
    log(f"peak absolute EPSC = {m.max():.2f} pA at pulse {int(m.argmax()) + 1}")

    # Facilitation is normalised to pulse 1, and Control's pulse-1 release
    # probability is low (~0.15), so at small n the denominator is a handful of
    # events and the ratio is meaningless. Refuse to quote it rather than print
    # a number that looks like a result.
    if m[0] < 3 * sem0:
        log(f"peak facilitation vs pulse 1 = {(m / m[0]).max():.1f}  <-- NOT RELIABLE: "
            f"pulse-1 mean ({m[0]:.2f} pA) is within 3 SEM of zero. "
            f"Use the absolute amplitudes, or collect more trials.")
    else:
        log(f"peak facilitation vs pulse 1 = {(m / m[0]).max():.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
