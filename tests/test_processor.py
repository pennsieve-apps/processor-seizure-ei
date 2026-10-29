"""
Unit tests on a synthetic fixture (tests/make_fixture.py); no Docker, no real data.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from processor import ei as E
from processor import handler
from processor import main as M
from tests import make_fixture


@pytest.fixture(scope="module")
def fixture_dir(tmp_path_factory) -> Path:
    return make_fixture.make(tmp_path_factory.mktemp("fixture"))


def run(monkeypatch, input_dir: Path, output_dir: Path, **env) -> int:
    for k in M.PARAMETERS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    return M.main()


def read_csv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


# --- the method -------------------------------------------------------------------

def test_energy_ratio_bands():
    fs = 512.0
    t = np.arange(int(20 * fs)) / fs
    x = np.vstack([np.sin(2 * np.pi * 30 * t), np.sin(2 * np.pi * 5 * t), np.sin(2 * np.pi * 30 * t) + np.sin(2 * np.pi * 10 * t)])
    er, tt = E.energy_ratio(x, fs, E.Params())
    assert er.shape == (3, len(tt)) and np.allclose(np.diff(tt), 0.25) and tt[0] == 0.5
    assert np.median(er[0]) > 100          # gamma only
    assert np.median(er[1]) < 0.01         # theta only
    assert 0.5 < np.median(er[2]) < 2      # equal gamma and alpha


def test_band_edges_follow_sampling_rate():
    assert E.band_edges(256, 97.0)["gamma"] == (24.0, 97.0)
    assert E.band_edges(200, 97.0)["gamma"] == (24.0, 90.0)


def test_page_hinkley_detects_the_change_point():
    rng = np.random.default_rng(0)
    er = np.concatenate([0.5 + 0.1 * rng.standard_normal(200), 8 + rng.standard_normal(100)])
    nd, na, u = E.page_hinkley(er, lam=50, nu=2)
    assert nd is not None and 198 <= nd <= 201 and nd <= na <= nd + 15
    nd, na, _ = E.page_hinkley(0.5 + 0.1 * rng.standard_normal(300), lam=50, nu=2)
    assert nd is None and na is None


def test_epileptogenicity_weights_early_and_strong():
    t = np.arange(400) * 0.25
    base = np.full(400, 0.5)
    early = base.copy(); early[200:] = 10       # 50 s, strong
    late = base.copy(); late[240:] = 10         # 60 s, as strong but 10 s later
    weak = base.copy(); weak[200:] = 3          # 50 s, weaker
    r = E.epileptogenicity(np.vstack([early, late, weak, base]), t, E.Params())
    assert r["ei"][0] == 1.0 and r["n_detected"] == 3
    assert r["ei"][1] < 0.15 and r["ei"][2] < r["ei"][0] and r["ei"][3] < 0.01
    assert np.isnan(r["detection_s"][3]) and abs(r["detection_s"][0] - 50) <= 0.5 and r["t0"] == r["detection_s"][0]


# --- end to end -----------------------------------------------------------------------

def test_main_on_derivative(fixture_dir, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(monkeypatch, fixture_dir / "linked", out) == 0
    # everything under a folder named after the method (one link of a fan-in port)
    assert [q.name for q in out.iterdir()] == ["ei"]
    rows = read_csv(out / "ei" / "ei.csv")
    assert set(rows[0]) == {"recording_id", "channel", "ei", "detection_s", "ratio_peak"}
    for rid in ("sub-01_task-ictal_run-01", "sub-01_task-ictal_run-02"):
        rr = {r["channel"]: r for r in rows if r["recording_id"] == rid}
        top2 = sorted(rr, key=lambda c: -float(rr[c]["ei"]))[:2]
        assert set(top2) == {"A01-A02", "A02-A03"}
        assert max(float(r["ei"]) for r in rr.values()) == 1.0
        for c in ("B01-B02", "B02-B03", "B03-B04"):
            assert float(rr[c]["ei"]) < 0.1 and rr[c]["detection_s"] == ""
        assert make_fixture.ONSET - 2 <= float(rr["A01-A02"]["detection_s"]) <= make_fixture.ONSET + 5
    prov = json.loads((out / "ei" / "provenance.json").read_text())
    st = prov["recordings"]["sub-01_task-ictal_run-01"]
    assert st["status"] == "ok" and st["window_s"] == [70.0, 150.0] and st["window_source"].startswith("onset_s")


def test_onsets_input_anchors_the_window(fixture_dir, tmp_path, monkeypatch):
    """`ieeg` and `onsets` in subfolders (linked inputs); a DynaSD onset moves the window."""
    inp = tmp_path / "in"
    # as seizure-dynasd writes it: its output folder dynasd/ inside the onsets port
    (inp / "onsets" / "dynasd").mkdir(parents=True)
    (inp / "ieeg").symlink_to(fixture_dir / "linked" / "ieeg", target_is_directory=True)
    (inp / "onsets" / "dynasd" / "seizure_onsets.csv").write_text(
        "recording_id,onset_s,method\nsub-01_task-ictal_run-01,95.0,dynasd-lindda\nsub-01_task-ictal_run-02,,dynasd-lindda\n")
    out = tmp_path / "out"
    assert run(monkeypatch, inp, out) == 0
    prov = json.loads((out / "ei" / "provenance.json").read_text())
    assert prov["recordings"]["sub-01_task-ictal_run-01"]["window_s"] == [75.0, 150.0]
    assert prov["recordings"]["sub-01_task-ictal_run-01"]["window_source"].startswith("seizure-dynasd onset")
    assert prov["recordings"]["sub-01_task-ictal_run-02"]["window_source"].startswith("onset_s (no DynaSD onset")
    assert "onsets/dynasd/seizure_onsets.csv" in prov["inputs"]
    out2 = tmp_path / "out2"
    assert run(monkeypatch, inp, out2, USE_ONSETS="false") == 0
    prov = json.loads((out2 / "ei" / "provenance.json").read_text())
    assert prov["recordings"]["sub-01_task-ictal_run-01"]["window_s"] == [70.0, 150.0]


def test_without_onset_uses_everything_after_the_baseline(tmp_path, monkeypatch):
    fx = make_fixture.make(tmp_path / "fx")
    ev = fx / "plain" / "sub-02" / "ieeg" / "sub-02_task-ictal_events.tsv"
    ev.unlink()
    (fx / "plain" / "sub-02" / "ieeg" / "sub-02_task-broken_ieeg.edf").unlink()
    out = tmp_path / "out"
    assert run(monkeypatch, fx / "plain", out) == 0
    st = json.loads((out / "ei" / "provenance.json").read_text())["recordings"]["sub-02_task-ictal"]
    assert st["window_s"] == [50.0, 150.0] and st["window_source"].startswith("after the baseline")
    rr = {r["channel"]: float(r["ei"]) for r in read_csv(out / "ei" / "ei.csv")}
    assert max(rr, key=rr.get) in ("A01-A02", "A02-A03")


def test_broken_file_is_skipped(fixture_dir, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(monkeypatch, fixture_dir / "plain", out) == 0
    prov = json.loads((out / "ei" / "provenance.json").read_text())
    assert prov["recordings"]["sub-02_task-broken"]["status"] == "failed"
    assert prov["recordings"]["sub-02_task-ictal"]["status"] == "ok"


def test_no_input_fails(tmp_path, monkeypatch):
    (tmp_path / "empty").mkdir()
    assert run(monkeypatch, tmp_path / "empty", tmp_path / "out") == 1


@pytest.mark.parametrize("env", [{"LAMBDA": "-1"}, {"NU": "x"}, {"USE_ONSETS": "maybe"}, {"WINDOW_S": "0"}])
def test_bad_parameters(fixture_dir, tmp_path, monkeypatch, env):
    assert run(monkeypatch, fixture_dir / "linked", tmp_path / "out", **env) == 2


def test_handler(fixture_dir, tmp_path, monkeypatch):
    for k in M.PARAMETERS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("LAMBDA", "999")   # left over from a previous warm invocation: must be cleared
    event = {"inputDir": str(fixture_dir / "linked"), "outputDir": str(tmp_path / "out"), "NU": 2, "USE_ONSETS": True}
    result = handler.handler(event, context=None)
    assert result["status"] == "ok" and (tmp_path / "out" / "ei" / "ei.csv").exists()
    prov = json.loads((tmp_path / "out" / "ei" / "provenance.json").read_text())
    assert prov["parameters"]["LAMBDA"] == 50.0 and prov["parameters"]["USE_ONSETS"] is True
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(RuntimeError):
        handler.handler({"inputDir": str(empty), "outputDir": str(tmp_path / "out2")}, context=None)
