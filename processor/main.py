"""
processor-seizure-ei: the Epileptogenicity Index (Bartolomei et al. 2008) per channel and seizure.

Reads the ieeg-prep derivative (or any BIDS-iEEG / EDF input) and, optionally, seizure_onsets.csv
from seizure-dynasd, and writes

  OUTPUT_DIR/ei.csv            recording_id, channel, ei, detection_s, ratio_peak
  OUTPUT_DIR/provenance.json

The analysed window of each recording is [anchor - WINDOW_BEFORE_S, anchor + WINDOW_AFTER_S],
where the anchor is the DynaSD onset (input `onsets`), else recordings.csv `onset_s` (or a
`seizure` event); without either, everything after the baseline.
"""
from __future__ import annotations

import csv
import json
import logging
import os
import platform
import sys
import time
import traceback
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path

import numpy as np
from scipy import signal

from processor import ei as E
from processor import ieeg as I

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("processor-seizure-ei")
for noisy in ("mne", "matplotlib"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "processor-seizure-ei"
REPO_ROOT = Path(__file__).resolve().parent.parent
PARAMETERS = ("LAMBDA", "NU", "H_S", "TAU_S", "WINDOW_S", "STEP_S", "GAMMA_HIGH_HZ",
              "WINDOW_BEFORE_S", "WINDOW_AFTER_S", "BASELINE_S", "USE_ONSETS")


def version_of(root: Path) -> str:
    try:
        return (root / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def tool_versions() -> dict[str, str]:
    out = {"python": platform.python_version()}
    for pkg in ("mne", "numpy", "scipy"):
        try:
            out[pkg] = importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            pass
    return out


def _num(name: str, default: float, lo: float | None = None, hi: float | None = None) -> float:
    raw = (os.environ.get(name) or "").strip()
    if not raw:
        return default
    try:
        v = float(raw)
    except ValueError as exc:
        raise ValueError(f"{name}={raw!r} is not a number") from exc
    if (lo is not None and v < lo) or (hi is not None and v > hi) or not np.isfinite(v):
        raise ValueError(f"{name}={raw!r} out of range [{lo}, {hi}]")
    return v


def env_params() -> tuple[E.Params, dict]:
    p = E.Params(
        lam=_num("LAMBDA", 50.0, lo=0.0), nu=_num("NU", 2.0, lo=0.0), h_s=_num("H_S", 5.0, lo=0.1),
        tau_s=_num("TAU_S", 1.0, lo=0.001), win_s=_num("WINDOW_S", 1.0, lo=0.25, hi=10.0),
        step_s=_num("STEP_S", 0.25, lo=0.01, hi=5.0), gamma_hi=_num("GAMMA_HIGH_HZ", 97.0, lo=30.0),
    )
    use = (os.environ.get("USE_ONSETS") or "true").strip().lower()
    if use not in ("true", "false"):
        raise ValueError(f"USE_ONSETS={use!r}: true or false")
    other = {"window_before_s": _num("WINDOW_BEFORE_S", 20.0, lo=0.0),
             "window_after_s": _num("WINDOW_AFTER_S", 60.0, lo=5.0),
             "baseline_s": _num("BASELINE_S", 60.0, lo=5.0), "use_onsets": use == "true"}
    return p, other


def notch(x: np.ndarray, fs: float, line: float | None) -> np.ndarray:
    lf = line if line in (50.0, 60.0) else 60.0
    for f0 in (lf, 2 * lf):
        if f0 + 2 < fs / 2:
            b, a = signal.butter(4, (f0 - 2, f0 + 2), "bandstop", fs=fs)
            x = signal.filtfilt(b, a, x, axis=-1)
    return x


def process_recording(rec: I.Recording, p: E.Params, o: dict, dynasd_onset: float | None,
                      have_onsets_file: bool) -> dict:
    t_start = time.monotonic()
    raw = I.open_raw(rec)
    fs = float(raw.info["sfreq"])
    if fs < 200:
        raise ValueError(f"sampling rate {fs:g} Hz < 200 Hz")
    duration = raw.n_times / fs
    chans, minfo = I.brain_channels(rec, raw)
    if len(chans) < 2:
        raise ValueError(f"{len(chans)} usable channel(s)")
    if dynasd_onset is not None:
        anchor, src = dynasd_onset, "seizure-dynasd onset"
    elif rec.onset_s is not None:
        anchor = rec.onset_s
        src = "onset_s (no DynaSD onset for this recording)" if have_onsets_file else "onset_s"
    else:
        anchor, src = None, "none"
    if anchor is not None:
        w0, w1 = max(0.0, anchor - o["window_before_s"]), min(duration, anchor + o["window_after_s"])
        window_src = f"{src}: anchor {anchor:.2f} s"
    else:
        base, bsrc = I.baseline_interval(rec, duration, o["baseline_s"])
        w0, w1 = base[1], duration
        window_src = f"after the baseline ({bsrc})"
    if w1 - w0 < p.win_s + p.h_s:
        raise ValueError(f"analysed window {w0:.1f}-{w1:.1f} s is too short")
    start, stop = int(np.floor(w0 * fs)), int(np.ceil(w1 * fs))
    ers, t = [], None
    for _, x in I.iter_blocks(raw, chans, start=start, stop=stop):
        er, t = E.energy_ratio(notch(x, fs, rec.line_freq), fs, p)
        ers.append(er)
    er = np.concatenate(ers)
    res = E.epileptogenicity(er, t + start / fs, p)
    res.update({"names": [c.name for c in chans], "sfreq": fs, "duration_s": round(duration, 3),
                "window": [round(start / fs, 3), round(stop / fs, 3)], "window_source": window_src,
                "montage": minfo["montage"], "channels_dropped": minfo["dropped"],
                "seconds": round(time.monotonic() - t_start, 2)})
    return res


def fmt(v, nd=4) -> str:
    if v is None or not np.isfinite(v):
        return ""
    return f"{v:.{nd}f}"


def main() -> int:
    t_run = time.monotonic()
    input_dir = Path(os.environ.get("INPUT_DIR") or "/data/input")
    output_dir = Path(os.environ.get("OUTPUT_DIR") or "/data/output")
    try:
        p, o = env_params()
    except ValueError as exc:
        log.error("%s", exc)
        return 2
    output_dir.mkdir(parents=True, exist_ok=True)
    recs = I.find_recordings(input_dir) if input_dir.exists() else []
    onsets = I.read_onsets_csv(input_dir) if (o["use_onsets"] and input_dir.exists()) else {}
    log.info("%s %s: %d recording(s) under %s; %d DynaSD onset(s); lambda=%g nu=%g H=%g s", APP_NAME,
             version_of(REPO_ROOT), len(recs), input_dir, sum(v is not None for v in onsets.values()), p.lam, p.nu,
             p.h_s)
    statuses, rows = {}, []
    for rec in recs:
        try:
            res = process_recording(rec, p, o, onsets.get(rec.recording_id), bool(onsets))
        except Exception as exc:  # noqa: BLE001 - a failed recording is logged and skipped
            log.error("%s: failed: %s", rec.recording_id, exc)
            log.debug("%s", traceback.format_exc())
            statuses[rec.recording_id] = {"status": "failed", "error": str(exc), "input": rec.rel}
            continue
        order = np.argsort(-res["ei"])
        for i, name in enumerate(res["names"]):
            rows.append({"recording_id": rec.recording_id, "channel": name, "ei": fmt(res["ei"][i]),
                         "detection_s": fmt(res["detection_s"][i], 3), "ratio_peak": fmt(res["ratio_peak"][i])})
        top = [f"{res['names'][j]} {res['ei'][j]:.2f}" for j in order[:3]]
        log.info("%s: %d channels (%s), window %.1f-%.1f s (%s), %d detected, %d with EI > 0.3, top: %s, %.1f s",
                 rec.recording_id, len(res["names"]), res["montage"], res["window"][0], res["window"][1],
                 res["window_source"], res["n_detected"], int((res["ei"] > 0.3).sum()), ", ".join(top),
                 res["seconds"])
        statuses[rec.recording_id] = {
            "status": "ok", "input": rec.rel, "n_channels": len(res["names"]), "n_detected": res["n_detected"],
            "n_ei_over_0_3": int((res["ei"] > 0.3).sum()), "reference_time_s": res["t0"],
            "window_s": res["window"], "window_source": res["window_source"], "montage": res["montage"],
            "channels_dropped": res["channels_dropped"], "sfreq": res["sfreq"], "duration_s": res["duration_s"],
            "seconds": res["seconds"],
        }
    with open(output_dir / "ei.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["recording_id", "channel", "ei", "detection_s", "ratio_peak"])
        w.writeheader()
        w.writerows(rows)
    ok = sum(1 for s in statuses.values() if s["status"] == "ok")
    params = {"LAMBDA": p.lam, "NU": p.nu, "H_S": p.h_s, "TAU_S": p.tau_s, "WINDOW_S": p.win_s, "STEP_S": p.step_s,
              "GAMMA_HIGH_HZ": p.gamma_hi, "WINDOW_BEFORE_S": o["window_before_s"],
              "WINDOW_AFTER_S": o["window_after_s"], "BASELINE_S": o["baseline_s"], "USE_ONSETS": o["use_onsets"],
              "bands_hz": E.BANDS}
    prov = {
        "app": APP_NAME, "version": version_of(REPO_ROOT), "tools": tool_versions(), "parameters": params,
        "inputs": [r.rel for r in recs] + sorted(
            str(q.relative_to(input_dir)) for q in I._walk(input_dir) if q.name == "seizure_onsets.csv")
        if input_dir.exists() else [],
        "started": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "seconds": round(time.monotonic() - t_run, 2), "recordings": statuses,
    }
    (output_dir / "provenance.json").write_text(json.dumps(prov, indent=2, default=str) + "\n")
    if not recs:
        log.error("no recordings found under %s", input_dir)
        return 1
    if ok == 0:
        log.error("every recording failed")
        return 1
    log.info("done: %d/%d recordings in %.1f s", ok, len(recs), time.monotonic() - t_run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
