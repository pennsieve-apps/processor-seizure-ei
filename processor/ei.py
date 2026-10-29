"""
The Epileptogenicity Index (EI) of Bartolomei, Chauvel & Wendling (2008), Brain 131:1818-1830,
reimplemented in Python from the paper's Methods.

Per channel (bipolar iEEG):

1. Energy ratio. Over a sliding window of D seconds (step S) the periodogram (Hann-windowed
   |FFT|^2) gives the energy in theta 3.5-7.4, alpha 7.4-12.4, beta 12.4-24 and gamma 24-97 Hz
   (the paper's Table 2), and

       ER[n] = (E_beta + E_gamma) / (E_theta + E_alpha).

2. Cumulative sum and Page-Hinkley detection of the rise in ER:

       U_N = sum_{n=1..N} (ER[n] - mean(ER[1..n]) - nu),      u_N = min_{n<=N} U_n,

   an alarm is raised at the first N_a with U_{N_a} - u_{N_a} > lambda, and the detection
   (change-point) time N_d is the last minimum of U before the alarm. nu is the size of a
   change that must not raise an alarm; lambda sets the false-alarm rate.

3. Epileptogenicity Index. With N_0 the earliest detection over all channels, tau = 1 s and
   H = 5 s:

       EI_i = (1 / (N_di - N_0 + tau)) * mean(ER_i[N_di .. N_di + H])

   (the paper sums over the H window; with the normalisation below the two are identical). A
   channel without a detection gets N_di = the end of the analysed window (the paper: "the time
   corresponding to seizure termination"). EI is divided by its maximum over the channels of the
   recording, so it lies in [0, 1]; the paper reads EI > 0.3 as epileptogenic.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

BANDS = {"theta": (3.5, 7.4), "alpha": (7.4, 12.4), "beta": (12.4, 24.0), "gamma": (24.0, 97.0)}


@dataclass
class Params:
    lam: float = 15.0        # lambda, Page-Hinkley threshold (ER units x windows)
    nu: float = 0.5          # bias
    h_s: float = 5.0         # H, seconds of ER averaged after detection
    tau_s: float = 1.0       # tau, seconds added to the delay
    win_s: float = 1.0       # D, periodogram window
    step_s: float = 0.25     # window step
    gamma_hi: float = 97.0   # upper edge of gamma (clipped to 0.45 x sampling rate)


def band_edges(fs: float, gamma_hi: float) -> dict[str, tuple[float, float]]:
    b = dict(BANDS)
    b["gamma"] = (b["gamma"][0], min(gamma_hi, 0.45 * fs))
    return b


def energy_ratio(x: np.ndarray, fs: float, p: Params) -> tuple[np.ndarray, np.ndarray]:
    """
    x: [n_ch, n] (one block of channels). Returns (er [n_ch, n_win], t [n_win]) where t is the
    centre of each window in seconds from the start of x.
    """
    n_win = int(round(p.win_s * fs))
    n_step = max(1, int(round(p.step_s * fs)))
    if x.shape[1] < n_win:
        raise ValueError("analysed span shorter than one ER window")
    starts = np.arange(0, x.shape[1] - n_win + 1, n_step)
    freqs = np.fft.rfftfreq(n_win, 1.0 / fs)
    edges = band_edges(fs, p.gamma_hi)
    masks = {k: (freqs >= lo) & (freqs < hi) for k, (lo, hi) in edges.items()}
    taper = np.hanning(n_win)
    er = np.empty((x.shape[0], len(starts)))
    # Frames in chunks to bound memory: [n_ch, chunk, n_win].
    chunk = max(1, int(2e7 // max(1, x.shape[0] * n_win)))
    for i in range(0, len(starts), chunk):
        idx = starts[i:i + chunk, None] + np.arange(n_win)[None, :]
        frames = x[:, idx]
        frames = frames - frames.mean(axis=-1, keepdims=True)
        spec = np.abs(np.fft.rfft(frames * taper, axis=-1)) ** 2
        e = {k: spec[..., m].sum(axis=-1) for k, m in masks.items()}
        low = e["theta"] + e["alpha"]
        er[:, i:i + chunk] = (e["beta"] + e["gamma"]) / np.maximum(low, np.finfo(float).tiny)
    t = (starts + n_win / 2.0) / fs
    return er, t


def page_hinkley(er: np.ndarray, lam: float, nu: float) -> tuple[int | None, int | None, np.ndarray]:
    """
    One channel's ER [n]. Returns (n_detect, n_alarm, U) with indices into er (None when no
    alarm). U_N = cumsum(ER - running mean - nu); alarm when U_N - min(U_1..U_N) > lambda;
    detection = the last minimum of U before the alarm.
    """
    n = np.arange(1, len(er) + 1)
    running_mean = np.cumsum(er) / n
    u = np.cumsum(er - running_mean - nu)
    umin = np.minimum.accumulate(u)
    over = np.flatnonzero(u - umin > lam)
    if len(over) == 0:
        return None, None, u
    na = int(over[0])
    # last index at which U attained its running minimum up to the alarm
    nd = int(np.flatnonzero(u[: na + 1] == umin[na])[-1])
    return nd, na, u


def epileptogenicity(er: np.ndarray, t: np.ndarray, p: Params) -> dict:
    """
    er [n_ch, n_win] for every channel of one recording, t [n_win] (seconds, any origin).
    Returns detection/alarm times (seconds, NaN when none), EI (normalised) and ER peak.
    """
    n_ch, n_win = er.shape
    step = float(np.median(np.diff(t))) if n_win > 1 else p.step_s
    h = max(1, int(round(p.h_s / step)))
    nd = np.full(n_ch, -1)
    na = np.full(n_ch, -1)
    for i in range(n_ch):
        d, a, _ = page_hinkley(er[i], p.lam, p.nu)
        if d is not None:
            nd[i], na[i] = d, a
    detected = nd >= 0
    end_idx = n_win - 1
    td = np.where(detected, t[np.clip(nd, 0, end_idx)], np.nan)
    ta = np.where(detected, t[np.clip(na, 0, end_idx)], np.nan)
    t_end = t[end_idx]
    t0 = np.nanmin(td) if detected.any() else t_end
    raw = np.zeros(n_ch)
    for i in range(n_ch):
        if detected[i]:
            seg = er[i, nd[i]: nd[i] + h + 1]
            delay = td[i] - t0
        else:
            seg = er[i, max(0, end_idx - h): end_idx + 1]
            delay = t_end - t0
        raw[i] = float(np.mean(seg)) / (delay + p.tau_s)
    top = raw.max() if n_ch and raw.max() > 0 else 1.0
    return {
        "ei": raw / top,
        "ei_raw": raw,
        "detection_s": td,
        "alarm_s": ta,
        "ratio_peak": er.max(axis=1) if n_win else np.full(n_ch, np.nan),
        "t0": t0,
        "n_detected": int(detected.sum()),
    }
