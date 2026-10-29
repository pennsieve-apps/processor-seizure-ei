"""
Scores ei.csv against expert labels (development tool; not in the image).

  python scripts/validate.py <output-dir> <truth> [annotations.tsv]

<truth> is either
  - the Multi-Expert Seizure Annotation dataset folder (Discover 530; has annotations.tsv):
    expert onset = consensus UEO time, onset channels = ueo_consensus contacts (a bipolar channel
    counts when either contact is listed), clip start = approximate_onset - 120 s; or
  - an ieeg-prep derivative (recordings.csv onset_s, _channels.tsv `onset` column).

Per recording: AUC of EI for the expert onset channels, how many of the 5 highest-EI channels
are onset channels, and the low-voltage fast activity flag (annotations.tsv `lvfa`). Prints the
table, median [range], and the split by LVFA.
"""
from __future__ import annotations

import ast
import csv
import re
import sys
from pathlib import Path

import numpy as np


def auc(scores: np.ndarray, labels: np.ndarray) -> float:
    """Mann-Whitney AUC (ties count 1/2); NaN when one class is empty."""
    pos, neg = scores[labels], scores[~labels]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    diff = pos[:, None] - neg[None, :]
    return float(((diff > 0) + 0.5 * (diff == 0)).mean())


def truth_from_annotations(folder: Path) -> dict:
    rows = list(csv.DictReader(open(folder / "annotations.tsv", encoding="utf-8-sig"), delimiter="\t"))
    out = {}
    for edf in folder.rglob("*_ieeg.edf"):
        rid = edf.name[: -len("_ieeg.edf")]
        pat = re.match(r"sub-([A-Za-z0-9]+)_", rid).group(1)
        n = int(re.search(r"task-[a-z]+(\d+)", rid).group(1))
        hit = [r for r in rows if r["patient"] == pat and int(float(r["approximate_onset"])) == n]
        if len(hit) != 1:
            continue
        r = hit[0]
        start = float(r["approximate_onset"]) - 120.0
        out[rid] = {"onset_s": float(r["ueo_time_consensus"]) - start,
                    "contacts": set(ast.literal_eval(r["ueo_consensus"])),
                    "lvfa": r["lvfa"] == "yes", "stim": r["stim"]}
    return out


def truth_from_derivative(folder: Path) -> dict:
    out = {}
    for idx in folder.rglob("recordings.csv"):
        for r in csv.DictReader(open(idx, encoding="utf-8-sig")):
            edf = idx.parent / r["path"]
            ch = edf.parent / edf.name.replace("_ieeg.edf", "_channels.tsv")
            labels = {}
            if ch.exists():
                for c in csv.DictReader(open(ch, encoding="utf-8-sig"), delimiter="\t"):
                    if (c.get("onset") or "").lower() in ("true", "false"):
                        labels[c["name"]] = c["onset"].lower() == "true"
            out[r["recording_id"]] = {"onset_s": float(r["onset_s"]) if r.get("onset_s") else None,
                                      "labels": labels, "lvfa": None, "stim": None}
    return out


def add_flags(truth: dict, annotations: Path) -> None:
    """lvfa / stim flags from the Multi-Expert annotations.tsv for a derivative's recordings."""
    flags = truth_from_annotations_rows(annotations)
    for rid, t in truth.items():
        f = flags.get(rid)
        if f:
            t["lvfa"], t["stim"] = f


def truth_from_annotations_rows(annotations: Path) -> dict:
    rows = list(csv.DictReader(open(annotations, encoding="utf-8-sig"), delimiter="\t"))
    out = {}
    for r in rows:
        out[(r["patient"], int(float(r["approximate_onset"])))] = (r["lvfa"] == "yes", r["stim"])

    class Lookup(dict):
        def get(self, rid, default=None):
            m = re.match(r"sub-([A-Za-z0-9]+)_.*task-[a-z]+(\d+)", rid)
            return out.get((m.group(1), int(m.group(2)))) if m else default
    return Lookup()


def is_onset(t: dict, channel: str) -> bool | None:
    if "labels" in t:
        return t["labels"].get(channel)
    return bool(set(channel.split("-")) & t["contacts"])


def main(out_dir: Path, truth_dir: Path, annotations: Path | None = None) -> None:
    truth = truth_from_annotations(truth_dir) if (truth_dir / "annotations.tsv").exists() else truth_from_derivative(truth_dir)
    if annotations is not None:
        add_flags(truth, annotations)
    rows: dict[str, list[dict]] = {}
    for r in csv.DictReader(open(out_dir / "ei.csv")):
        rows.setdefault(r["recording_id"], []).append(r)
    res = []
    for rid in sorted(rows):
        t = truth.get(rid)
        if t is None:
            continue
        rr = [r for r in rows[rid] if r["ei"] != "" and is_onset(t, r["channel"]) is not None]
        lab = np.array([is_onset(t, r["channel"]) for r in rr], bool)
        ei = np.array([float(r["ei"]) for r in rr])
        top = lab[np.argsort(-ei)[:5]]
        res.append({"recording": rid, "auc": auc(ei, lab), "top5": int(top.sum()), "n_onset": int(lab.sum()),
                    "n_ch": len(rr), "n_over_0_3": int((ei > 0.3).sum()),
                    "onset_over_0_3": int((lab & (ei > 0.3)).sum()), "lvfa": t["lvfa"], "stim": t["stim"]})
    print(f"{'recording':55s} {'auc':>5s} {'top5':>4s} {'n_on':>4s} {'>0.3':>4s} {'on>0.3':>6s} {'n_ch':>4s} lvfa stim")
    for r in res:
        print(f"{r['recording']:55s} {r['auc']:5.2f} {r['top5']:4d} {r['n_onset']:4d} {r['n_over_0_3']:4d} "
              f"{r['onset_over_0_3']:6d} {r['n_ch']:4d} {str(r['lvfa'])[:1]:4s} {str(r['stim'])[:5]}")

    def summ(name, vals):
        v = np.array([x for x in vals if np.isfinite(x)])
        return f"{name}: median {np.median(v):.2f} [{v.min():.2f}, {v.max():.2f}] (n={len(v)})" if len(v) else f"{name}: n=0"

    print(summ("AUC", [r["auc"] for r in res]))
    for flag in (True, False):
        print(summ(f"AUC lvfa={flag}", [r["auc"] for r in res if r["lvfa"] is flag]))
    for kind in sorted({str(r["stim"]) for r in res}):
        print(summ(f"AUC {kind}", [r["auc"] for r in res if str(r["stim"]) == kind]))
    print(summ("onset channels among the top 5", [r["top5"] for r in res]))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]) if len(sys.argv) > 3 else None)
