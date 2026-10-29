# processor-seizure-ei

The Epileptogenicity Index (EI) step of the Seizure Onset Localisation workflow
(`ieeg-prep → seizure-dynasd → {seizure-ei, seizure-fragility} → seizure-consensus → ieeg-report`).

**For an epileptologist.** EI is the quantity Bartolomei, Chauvel & Wendling (2008) proposed to
put a number on what a reader of SEEG looks for at seizure onset: which contacts switch to
*fast* activity (beta/gamma, typically low-voltage fast activity), and *how early*. Each
channel gets a value between 0 and 1 for each seizure: 1 for the channel with the strongest and
earliest rise of fast activity, small values for channels recruited late or never. In the
paper, structures with EI > 0.3 were considered epileptogenic. It is a ranking of channels
within one seizure, not a probability, and it says nothing about channels that were not
implanted.

**For a researcher.** A Python reimplementation of the Methods of Bartolomei et al. (2008),
run per bipolar channel and per seizure:

1. **Energy ratio.** A periodogram (Hann-tapered |FFT|²) over a sliding 1 s window (0.25 s step)
   gives the energies in θ 3.5–7.4, α 7.4–12.4, β 12.4–24 and γ 24–97 Hz (the paper's bands),
   and ER[n] = (E_β + E_γ) / (E_θ + E_α).
2. **Detection (Page–Hinkley).** U_N = Σ_{n≤N} (ER[n] − mean(ER[1..n]) − ν), u_N = min_{n≤N} U_n;
   the alarm is the first N with U_N − u_N > λ, and the detection (change-point) time N_d is
   the last minimum of U before the alarm.
3. **Index.** With N_0 the earliest detection of the seizure, τ = 1 s and H = 5 s:
   EI_i = mean(ER_i[N_di … N_di + H]) / (N_di − N_0 + τ), divided by the maximum over the
   channels of the seizure. A channel never detected gets N_di = the end of the analysed
   window (the paper: "seizure termination").

The analysed window is placed around the seizure: 20 s before to 60 s after the onset found by
**seizure-dynasd** (input `onsets`), else the recording's own onset (`recordings.csv onset_s`
or a `seizure` event); with neither, everything after the baseline.

## Ports

| Port | Direction | Placement | What |
|---|---|---|---|
| `ieeg` | input | `INPUT_DIR/ieeg/` | The ieeg-prep derivative (`recordings.csv`, EDF + `_channels.tsv` + `_events.tsv` + `_ieeg.json`), or any BIDS-iEEG dataset or EDF / BrainVision / FIF files, found by name anywhere under `INPUT_DIR` (symlinks followed). |
| `onsets` | input, optional | `INPUT_DIR/onsets/` | `seizure_onsets.csv` from seizure-dynasd (recording_id, onset_s). |
| `ei` | output | `OUTPUT_DIR` | `ei.csv` and `provenance.json`. |

Without links (everything merged at the root of `INPUT_DIR`) the files are found by name too.

## Parameters

| Name | Default | What it does |
|---|---|---|
| `LAMBDA` | `50` | λ, the Page–Hinkley threshold: an alarm when the cumulative excess of ER over its running mean (minus ν), summed over 0.25 s steps, exceeds λ. Higher = fewer, later detections. |
| `NU` | `2` | ν, the bias: a rise in ER of this size per step does not count. |
| `H_S` | `5` | H, seconds of ER averaged after each detection (paper: 5 s). |
| `TAU_S` | `1` | τ, added to each delay (paper: 1). |
| `WINDOW_S` | `1` | Periodogram window length (s). |
| `STEP_S` | `0.25` | Its step (s). λ is in ER × steps, so change both together. |
| `GAMMA_HIGH_HZ` | `97` | Upper edge of γ (paper: 97 Hz, at 256 Hz sampling); clipped to 0.45 × the sampling rate. |
| `WINDOW_BEFORE_S` | `20` | The window starts this long before the anchor onset. |
| `WINDOW_AFTER_S` | `60` | … and ends this long after it (or at the end of the recording). |
| `BASELINE_S` | `60` | Without any onset: the window is everything after the baseline (the input's `baseline`, else the first `BASELINE_S` s). |
| `USE_ONSETS` | `true` | Anchor on the DynaSD onsets when the `onsets` input is given; `false` uses `onset_s`. A recording DynaSD found no onset for falls back to `onset_s`. |

### Parameter values: the paper, EZEI and here

| | Bartolomei et al. 2008 | EZEI (KarasLabRAVE) | this app |
|---|---|---|---|
| Bands (Hz) | θ 3.5–7.4, α 7.4–12.4, β 12.4–24, γ 24–97 | same, γ to 140 | paper's (γ to 97, clipped to 0.45 fs) |
| Spectral estimate | periodogram, sliding window D (not stated) | multitaper, window/step given by the user (1 s / 0.2 s in its example) | Hann periodogram, 1 s / 0.25 s |
| λ, ν | "defined by the epileptologist" per seizure; no values | λ = 15, ν = 0.5 | λ = 50, ν = 2 |
| H, τ | 5 s, 1 | 5 s, 1 | 5 s, 1 s |
| Non-detected channel | N_d = seizure termination | N_d = end − H | N_d = end of the window |
| Normalisation | by the maximum (per patient) | by the maximum | by the maximum over the seizure's channels |

λ and ν were chosen on the demo data from a coarse grid (λ ∈ {15, 50, 150} × ν ∈ {0.5, 1, 2, 5}
× window start 10 or 20 s before the onset): the median AUC is flat at 0.86–0.91 for
λ = 50–150 and ν = 1–5, and drops to 0.81–0.82 with EZEI's λ = 15, ν = 0.5 (which, at a 0.25 s
step, flags half of all channels within the first seconds of the window). This is an in-sample
choice on 24 seizures of 6 patients; change λ/ν with care, and keep `STEP_S` at 0.25 s when
you do.

## Outputs

**`ei.csv`**: one row per channel and recording.

| Column | |
|---|---|
| `recording_id` | the recording (BIDS entity chain, e.g. `sub-HUP257_ses-postimplant_task-ictal387681_run-01`) |
| `channel` | the bipolar channel (`LA01-LA02`) |
| `ei` | Epileptogenicity Index, 0–1 within the recording (1 = most epileptogenic; > 0.3 in the paper's reading) |
| `detection_s` | change-point time N_d (s from the recording start); empty when the channel was never detected |
| `ratio_peak` | the largest energy ratio in the analysed window |

**`provenance.json`**: app and version, tool versions, parameters (and the bands), inputs seen
(recordings and `seizure_onsets.csv`), timings, and per recording: status, analysed window and
where its anchor came from, channels, detections, number of channels with EI > 0.3, the
reference time N_0, montage, dropped channels.

## Validation on the demo data

24 SEEG seizures (6 patients, Multi-Expert Seizure Annotation, doi:10.26275/dely-fzwx), raw
clips re-referenced to neighbouring-contact bipolar pairs (85–208 channels). Truth: the
consensus of 3–5 epileptologists' onset channels (`ueo_consensus`; a bipolar channel counts
when either contact is listed); AUC of EI for those channels, per seizure
(`scripts/validate.py`):

| Window anchored on | AUC, median [range] | LVFA seizures (n = 9) | other (n = 15) | spontaneous (18) | stimulation-induced (6) |
|---|---|---|---|---|---|
| expert consensus onset (as ieeg-prep's `onset_s` will give) | **0.89** [0.26, 0.99] | 0.88 [0.67, 0.97] | 0.89 [0.26, 0.99] | 0.90 [0.34, 0.99] | 0.72 [0.26, 0.93] |
| clinically marked onset (the clip's `ictal` event) | 0.88 [0.27, 0.99] | 0.83 [0.54, 0.97] | 0.90 [0.27, 0.99] | 0.90 [0.39, 0.99] | 0.58 [0.27, 0.96] |
| seizure-dynasd onset (default settings) | 0.87 [0.31, 0.99] | 0.87 [0.56, 0.99] | 0.88 [0.31, 0.98] | 0.88 [0.31, 0.98] | 0.69 [0.44, 0.99] |
| none (everything after a 60 s baseline)¹ | 0.69 [0.24, 0.99] | 0.62 | 0.82 | | |

¹ replayed on cached energy ratios with the same code, not a full app run.

- EI ranks the expert onset channels well (median AUC ≈ 0.88) as long as the window is placed
  around the seizure; searching a whole clip lets pre-ictal bursts of fast activity win
  (0.69). The top of the ranking is less precise than the AUC suggests: of the 5 highest-EI
  channels, a median of 0.5–1 is an expert onset channel, because the experts mark 1–13
  channels among 85–208 and EI often puts a channel recruited a few seconds later, or one
  with a short burst of fast activity, first.
- **LVFA.** EI is known to work best on seizures that start with low-voltage fast activity
  (it measures exactly that). **We do not see that here**: AUC is 0.88 for the 9 LVFA
  seizures and 0.89 for the 15 others (expert anchor; 0.83 vs 0.90 with the clinical mark). With
  bipolar SEEG at 1024 Hz many non-LVFA onsets (e.g. low-amplitude rhythmic or spiky onsets)
  also raise the β/γ share; and the sample is small (6 patients, several of whose LVFA
  seizures are stimulation-induced).
- Stimulation-induced seizures (one per patient) are clearly worse: the stimulation artefacts
  before the onset are broadband and raise ER.

**On the ieeg-prep derivative** (`workflow-epilepsy-ieeg/data/work/prep/ieeg`: same 24 clips,
bipolar, bad channels marked, `onset_s` = the expert consensus onset; truth = its
`_channels.tsv` `onset` column). `onset_s` is the expert onset, so anchoring on it is
optimistic; the onset-blind chain removes onset and baseline from the input and lets DynaSD
find the onset first:

| Window anchored on | AUC, median [range] | LVFA (9) | other (15) | spontaneous (18) | stim-induced (6) |
|---|---|---|---|---|---|
| seizure-dynasd onset, DynaSD itself searching near the expert onset (the workflow default; `data/work/ei`) | **0.90** [0.34, 0.98] | 0.87 | 0.91 | 0.90 | 0.85 |
| the expert onset (`onset_s`, `USE_ONSETS=false`; optimistic) | 0.89 [0.34, 0.99] | 0.89 | 0.93 | 0.92 | 0.75 |
| **onset-blind**: DynaSD onset found without any onset in the input | 0.88 [0.31, 0.98] | 0.87 | 0.89 | 0.90 | 0.74 |
| none (everything after the first 60 s) | 0.68 [0.25, 0.98] | 0.62 | 0.80 | 0.78 | 0.64 |

Same picture as on the raw clips: ≈ 0.88–0.90 whenever the window is placed by an onset (even
DynaSD's own, found blind), 0.68 without one, and no LVFA advantage.

**Cross-check against EZEI.** On one clip (sub-HUP224 ictal519177, 10 channels, 100–180 s), the
R package EZEI (KarasLabRAVE, run unchanged in its own R container, not bundled) and this app
with EZEI's λ = 15, ν = 0.5 and γ to 140 Hz give: detection times equal within 0.25 s (one step)
on 8 of 10 channels (the other two differ by 1.25 and 14 s; EZEI's multitaper ER is smoother);
ER series on the same scale (median ratio 0.93–1.14 per channel, correlation of log ER
0.82–0.93). The EI values differ, because EZEI's EI line indexes the time axis with
`Nd[ie]:Nd[ie]+hspan` (R precedence: averages one sample) and `stimes[Nd[ie]-t0+tau]`
(mixes indices and seconds); this app follows the paper's formula.

## Runtime

`computeTypes: [standard, lambda]`, 1 vCPU / 2 GiB, 15 min timeout. Only the analysed window
(80 s by default) is read, one block of 16 channels at a time, so memory hardly depends on
the recording's length (without any onset the whole recording after the baseline is read, a
block at a time: ~0.5 GB per hour of 1024 Hz data in a block). On the 24 ieeg-prep clips the image takes 18.6 s
for the whole run.

Measured in Docker (linux/amd64 **emulated on an M-series Mac**, `--cpus 1 --memory 3008m`)
on 3 demo clips (208, 116 and 106 channels, 1024 Hz): **8.6, 2.3 and 2.1 s** per clip (22 s for
the run including start-up), cgroup `memory.peak` **0.44 GB**. Natively (arm64, same Mac), the 24 clips take
20–32 s together. The image is 125 MB compressed.

## Running it locally

```sh
python3 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python        # unit tests on a synthetic fixture
make run PYTHON=.venv/bin/python         # Docker over the fixture -> data/output/
make run-demo DEMO=../workflow-epilepsy-ieeg/data/work/prep/ieeg ONSETS=../processor-seizure-dynasd/data/output-demo
make validate TRUTH=../workflow-epilepsy-ieeg/data/seizures/DATA
```

## Limitations

- EI measures one onset pattern: a rise of β/γ over θ/α energy. Onsets without fast activity
  (slow rhythmic, DC shifts, spike-wave) or with very high-frequency activity above 97 Hz
  are seen weakly or not at all; electrodecremental onsets can even lower ER.
- It ranks channels **within one seizure** (normalised to its maximum); values are not
  comparable across seizures or patients, and an EI of 1 exists in every seizure, focal or not.
- The window matters: anchor it on a sensible onset (DynaSD's or the clinical mark). Pre-ictal
  bursts of fast activity, stimulation artefacts and muscle on superficial contacts raise ER.
- λ and ν have no published values; the defaults were picked on the demo (see above).
- A recording that is not bipolar is re-referenced to neighbouring contacts named
  `<lead><number>`; contacts named otherwise are left out.

## Licences

Apache-2.0 for this app's code (`LICENSE`). The image bundles MNE-Python, NumPy, SciPy (BSD)
and their dependencies, and awslambdaric (Apache-2.0); see `THIRD_PARTY.md`. No GPL code: EZEI
(GPL-3.0-or-later) was only used, separately, to cross-check numbers.

## Citations

The method:

> Bartolomei F, Chauvel P, Wendling F. **Epileptogenicity of brain structures in human temporal
> lobe epilepsy: a quantified study from intracerebral EEG.** Brain 131(7):1818–1830 (2008).
> doi:10.1093/brain/awn111, PMID 18556663.

Reading EDF / BrainVision / FIF:

> Gramfort A, Luessi M, Larson E, Engemann DA, Strohmeier D, Brodbeck C, Goj R, Jas M, Brooks T,
> Parkkonen L, Hämäläinen M. **MEG and EEG data analysis with MNE-Python.** Front Neurosci 7:267
> (2013). doi:10.3389/fnins.2013.00267, PMID 24431986.

Demo data:

> Ojemann WK, Zhou D, Armstrong C, Sinha N, Litt B, Conrad E. **Multi-Expert Seizure
> Annotation.** Pennsieve Discover (2026). doi:10.26275/dely-fzwx (CC BY 4.0).

Reference implementation used for the cross-check only: EZEI, https://github.com/KarasLabRAVE/EZEI
(GPL-3.0-or-later).

DOIs were checked against Crossref (articles) and DataCite (dataset): title, authors and year.
Crossref lists only the first author of the MNE-Python paper; the full author list is from
PubMed (PMID 24431986).
