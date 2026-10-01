#!/usr/bin/env python3
"""
Runs ON ada.local. Aggregates per-trial paired-pulse data for Control and
AD+5xRyR into small summary files, following the exact conventions in
~/analysisscripts/modFuncNew.py (relppf) and caStat.py/ca-serca-peak-trials.py:

- AP nominal times: ts[i] = (i*ISI_ms + 2.0)/1000  seconds  (2ms delay, matches relppf)
- "Released" success window: (ts[i], ts[i]+0.02)              (tc=0.02, matches relppf)
- VDCC influx point-sample: flux value at ts[i]+0.001s          (matches caVDCC's 300/2300 index convention)
- AZ calcium: ca.dat col1 (time,col1,col2,col3 = t, Ca.Conc.az, Ca.Pre, Ca.ER),
  transformed via computeCaConc (t*x differenced over step=5), NOT used raw.
- Cytosolic calcium: ca.dat col2 / 602.3 (molecule count -> uM), used directly.
- SERCA occupied sites(t) = X1A(t) + Y1A(t) + 2*(X2(t)+Y2(t)); serca_mol.dat cols
  0..6 = time, X1, X1A, X2, Y2, Y1A, Y1. Total sites = 8678 (established divisor).

Usage:
  python3 remote_pairedpulse_agg.py mode [args...]

Modes:
  sweep_vdcc  <cond> <isi_ms> <vdcc_list_csv> <max_trials> <out_csv>
  sweep_isi   <cond> <vdcc> <isi_list_csv> <max_trials> <out_csv>
  deepdive    <cond> <isi_ms> <vdcc> <max_trials> <out_prefix>
"""
import os
import sys
import numpy as np

BASE = {
    "control": "/data/surabhi/Output/model_testing/prvsppf_modeltesting/nish/control/nish2sur_nishRyR",
    "ad_5x_ryr": "/data/surabhi/Output/model_testing/prvsppf_modeltesting/nish/AD/nish2sur_nishRyR",
}
SUFFIX = {
    "control": "C250uM_nish2sur_nishRyR",
    "ad_5x_ryr": "C750uM_5xryr_nish2sur_nishRyR",
}

AP_DELAY_S = 0.002
SUCCESS_WINDOW_S = 0.02
VDCC_SAMPLE_OFFSET_S = 0.001
DT = 1e-5
SERCA_TOTAL_SITES = 8678.0


def cond_dir(cond, isi_ms, vdcc):
    return os.path.join(BASE[cond], f"RSI{isi_ms}V{vdcc}{SUFFIX[cond]}")


def trial_dirs(cdir, cap):
    dirs = sorted(d for d in os.listdir(cdir) if d.startswith("s_"))
    return dirs[:cap]


def ap_times(isi_ms, n=2):
    return [(i * isi_ms + AP_DELAY_S * 1000) / 1000.0 for i in range(n)]


def p11_and_pr(cdir, isi_ms, cap):
    ts = ap_times(isi_ms)
    dirs = trial_dirs(cdir, cap)
    n0 = n1 = n01 = 0
    ndirs = 0
    for d in dirs:
        relpath = os.path.join(cdir, d, "dat", "rel.dat")
        if not os.path.exists(relpath):
            continue
        ndirs += 1
        p0 = p1 = 0
        if os.path.getsize(relpath) > 0:
            with open(relpath) as f:
                for line in f:
                    parts = line.split(None, 1)
                    if not parts:
                        continue
                    t = float(parts[0])
                    if ts[0] < t < ts[0] + SUCCESS_WINDOW_S:
                        p0 = 1
                    elif ts[1] < t < ts[1] + SUCCESS_WINDOW_S:
                        p1 = 1
        n0 += p0
        n1 += p1
        n01 += 1 if (p0 and p1) else 0
    if ndirs == 0:
        return float("nan"), float("nan"), float("nan"), 0
    pr1 = n0 / ndirs
    pr2 = n1 / ndirs
    p11 = (n01 / ndirs) / pr1 if pr1 > 0 else float("nan")
    return pr1, pr2, p11, ndirs


def read_rows(fpath, row_indices, col):
    wanted = set(row_indices)
    max_idx = max(row_indices)
    out = {}
    with open(fpath) as f:
        for i, line in enumerate(f):
            if i in wanted:
                parts = line.split()
                if len(parts) > col:
                    out[i] = float(parts[col])
            if i >= max_idx:
                break
    return out


def vdcc_influx_cv(cdir, isi_ms, cap):
    ts = ap_times(isi_ms)
    idx = [int(round((t + VDCC_SAMPLE_OFFSET_S) / DT)) + 1 for t in ts]  # +1: header line
    dirs = trial_dirs(cdir, cap)
    vals = {0: [], 1: []}
    for d in dirs:
        fpath = os.path.join(cdir, d, "dat", "vdcc_pq_ca_flux.dat")
        if not os.path.exists(fpath):
            continue
        rows = read_rows(fpath, idx, col=1)
        for k, i in enumerate(idx):
            if i in rows:
                vals[k].append(rows[i])
    out = {}
    for k in (0, 1):
        arr = np.array(vals[k])
        if len(arr) == 0 or arr.mean() == 0:
            out[k] = (float("nan"), float("nan"), len(arr))
        else:
            out[k] = (arr.mean(), arr.std(ddof=1) / arr.mean(), len(arr))
    return out


def compute_ca_conc(time, ca, step=5):
    c_tc = time * ca
    dt = step * (time[1] - time[0])
    n = (len(time) - step - 1) // step
    idx = np.arange(0, n * step, step)
    t_out = time[idx]
    c_out = (c_tc[idx + step] - c_tc[idx]) / dt
    return t_out, c_out


def deepdive(cond, isi_ms, vdcc, cap, out_prefix):
    cdir = cond_dir(cond, isi_ms, vdcc)
    ts = ap_times(isi_ms)
    win = SUCCESS_WINDOW_S
    dirs = trial_dirs(cdir, cap)

    peak_cyto = {0: [], 1: []}
    peak_az = {0: [], 1: []}
    serca_max_pct = {0: [], 1: []}
    serca_delta_pct = {0: [], 1: []}

    az_traces = []
    az_time = None

    for n, d in enumerate(dirs):
        cafile = os.path.join(cdir, d, "dat", "ca.dat")
        sfile = os.path.join(cdir, d, "dat", "serca_mol.dat")
        if not (os.path.exists(cafile) and os.path.exists(sfile)):
            continue

        t, az_raw, cyto_raw, _er = np.genfromtxt(cafile, usecols=(0, 1, 2, 3), unpack=True)
        t_az, az_conc = compute_ca_conc(t, az_raw, step=5)
        cyto_conc = cyto_raw / 602.3

        if az_time is None:
            az_time = t_az
        if len(az_conc) == len(az_time):
            az_traces.append(az_conc)

        s = np.genfromtxt(sfile, usecols=(0, 1, 2, 3, 4, 5, 6), unpack=True)
        st, sx1, sx1a, sx2, sy2, sy1a, sy1 = s
        occ = sx1a + sy1a + 2.0 * (sx2 + sy2)
        occ_pct = occ / SERCA_TOTAL_SITES * 100.0

        for k in (0, 1):
            # Every metric gets exactly one (possibly NaN) value per trial per AP,
            # so the four arrays stay index-aligned to the same trial -- required
            # since panels e/f plot (serca_max_pct, peak_cyto) as paired points.
            m_cyto = (t >= ts[k]) & (t < ts[k] + win)
            peak_cyto[k].append(cyto_conc[m_cyto].max() if m_cyto.any() else np.nan)

            m_az = (t_az >= ts[k]) & (t_az < ts[k] + win)
            peak_az[k].append(az_conc[m_az].max() if m_az.any() else np.nan)

            m_s = (st >= ts[k]) & (st < ts[k] + win)
            if m_s.any():
                init_occ = occ_pct[m_s][0]
                max_occ = occ_pct[m_s].max()
                serca_max_pct[k].append(max_occ)
                serca_delta_pct[k].append(max_occ - init_occ)
            else:
                serca_max_pct[k].append(np.nan)
                serca_delta_pct[k].append(np.nan)

    # Save scalar per-trial arrays (one row per AP, index-aligned across columns)
    with open(out_prefix + "_scalars.csv", "w") as f:
        f.write("ap,peak_cyto,peak_az,serca_max_pct,serca_delta_pct\n")
        for k in (0, 1):
            for i in range(len(peak_cyto[k])):
                f.write(f"{k+1},{peak_cyto[k][i]},{peak_az[k][i]},"
                        f"{serca_max_pct[k][i]},{serca_delta_pct[k][i]}\n")

    # Save AZ CV(t) trace
    if az_traces:
        az_traces = np.array(az_traces)
        mean_t = az_traces.mean(axis=0)
        std_t = az_traces.std(axis=0, ddof=1)
        cv_t = np.divide(std_t, mean_t, out=np.full_like(mean_t, np.nan), where=mean_t != 0)
        np.savetxt(out_prefix + "_az_cv_over_time.csv",
                   np.column_stack([az_time, mean_t, std_t, cv_t]),
                   header="time_s,mean_az_uM,std_az_uM,cv_az", delimiter=",", comments="")

    print(f"deepdive done: {cond} isi={isi_ms} vdcc={vdcc} trials={len(dirs)} "
          f"az_traces={len(az_traces)}")


def main():
    mode = sys.argv[1]
    if mode == "sweep_vdcc":
        cond, isi_ms, vdcc_csv, cap, out_csv = sys.argv[2:7]
        isi_ms = int(isi_ms)
        cap = int(cap)
        vdccs = [int(v) for v in vdcc_csv.split(",")]
        with open(out_csv, "w") as f:
            f.write("vdcc,pr1,pr2,p11,ntrials,vdcc_influx_ap1_mean,vdcc_influx_ap1_cv,vdcc_influx_ap2_mean,vdcc_influx_ap2_cv\n")
            for vdcc in vdccs:
                cdir = cond_dir(cond, isi_ms, vdcc)
                if not os.path.isdir(cdir):
                    print("MISSING", cdir)
                    continue
                pr1, pr2, p11, n = p11_and_pr(cdir, isi_ms, cap)
                vinf = vdcc_influx_cv(cdir, isi_ms, cap)
                f.write(f"{vdcc},{pr1},{pr2},{p11},{n},"
                        f"{vinf[0][0]},{vinf[0][1]},{vinf[1][0]},{vinf[1][1]}\n")
                print(cond, vdcc, pr1, pr2, p11, n)

    elif mode == "sweep_isi":
        cond, vdcc, isi_csv, cap, out_csv = sys.argv[2:7]
        vdcc = int(vdcc)
        cap = int(cap)
        isis = [int(v) for v in isi_csv.split(",")]
        with open(out_csv, "w") as f:
            f.write("isi_ms,pr1,pr2,p11,ntrials\n")
            for isi_ms in isis:
                cdir = cond_dir(cond, isi_ms, vdcc)
                if not os.path.isdir(cdir):
                    print("MISSING", cdir)
                    continue
                pr1, pr2, p11, n = p11_and_pr(cdir, isi_ms, cap)
                f.write(f"{isi_ms},{pr1},{pr2},{p11},{n}\n")
                print(cond, isi_ms, pr1, pr2, p11, n)

    elif mode == "deepdive":
        cond, isi_ms, vdcc, cap, out_prefix = sys.argv[2:7]
        deepdive(cond, int(isi_ms), int(vdcc), int(cap), out_prefix)

    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main()
