# Third-party software

This app is Apache-2.0 (see `LICENSE`). The Epileptogenicity Index is reimplemented here from
the paper (Bartolomei, Chauvel & Wendling 2008); no code of other EI implementations is
included. KarasLabRAVE/EZEI (R, GPL-3.0-or-later) was only run separately, in its own
container, to cross-check the numbers (see README "Cross-check").

The image contains the following components, installed from PyPI at the versions pinned in
`requirements.txt`. All are permissively licensed; nothing GPL is linked in-process.

| Component | Version | What it does here | Licence | Source |
|---|---|---|---|---|
| MNE-Python | 1.13.2 | Reads EDF / BrainVision / FIF | BSD-3-Clause | https://mne.tools |
| NumPy | 2.4.6 | Periodograms (FFT), the EI arithmetic | BSD-3-Clause (bundled parts 0BSD, MIT, Zlib, CC0-1.0) | https://numpy.org |
| SciPy | 1.18.1 | Notch filter | BSD-3-Clause | https://scipy.org |
| awslambdaric (+ simplejson, snapshot-restore-py) | 3.1.1 (4.2.0, 1.0.0) | Lambda runtime interface client (`computeTypes: lambda`) | Apache-2.0 (MIT / AFL-2.1 for simplejson; Apache-2.0) | https://github.com/aws/aws-lambda-python-runtime-interface-client |
| Matplotlib (+ contourpy, cycler, fonttools, kiwisolver, pillow, pyparsing) | 3.11.2 (1.4.0, 0.12.1, 4.68.0, 1.5.1, 13.0.0, 3.3.3) | MNE dependency (no plots here) | Matplotlib License (PSF-based); BSD-3-Clause; BSD-3-Clause; MIT; BSD-3-Clause; MIT-CMU; MIT | https://matplotlib.org |
| decorator, lazy-loader, Jinja2, MarkupSafe, pooch, platformdirs, packaging, tqdm, requests, urllib3, idna, charset-normalizer, certifi | 5.3.1, 0.6, 3.1.6, 3.0.4, 1.9.0, 4.14.0, 26.3, 4.70.1, 2.34.2, 2.8.0, 3.20, 3.5.2, 2026.10.23 | MNE dependencies (no network is used at run time) | BSD-2-Clause; BSD-3-Clause; BSD-3-Clause; BSD-3-Clause; BSD-3-Clause; MIT; Apache-2.0 OR BSD-2-Clause; MPL-2.0 AND MIT; Apache-2.0; MIT; BSD-3-Clause; MIT; MPL-2.0 | PyPI |
| python-dateutil, six, tzdata | 3.0.0, 1.17.0, 2026.6 | dependencies | Apache-2.0 / BSD dual; MIT; Apache-2.0 | PyPI |
| python:3.12.16-slim-bookworm | — | Base image (Debian 12 + CPython 3.12, pip) | PSF-2.0 (Python); Debian packages under their own free licences; pip MIT | https://hub.docker.com/_/python |

Demo data used for testing (not in the image or the repository): Ojemann WK, Zhou D,
Armstrong C, Sinha N, Litt B, Conrad E. Multi-Expert Seizure Annotation. Pennsieve Discover
2026, doi:10.26275/dely-fzwx, CC BY 4.0.
