# AGENTS.md

Notes for coding agents working in this repository.  `README.md` explains what
is compared, how, and why the tolerances are what they are.  Read it before
you change a check or a tolerance.

## What this repository is

A harness, not a library.  It runs `uwa-channels/python` and
`uwa-channels/matlab` on the same released Zenodo files and reports where they
differ.  A difference is a finding about one of those two repositories.  Fix it
there, not here.

## Running it

```bash
pip install numpy scipy matplotlib h5py
python calib.py                      # fetch the Zenodo files (740 MB)
python make_probe.py                 # write artifacts/probe.mat
PYTHONPATH=../replay_python/src python run_python.py
matlab -batch run_matlab             # uses ../replay_matlab/src, or $UWA_MATLAB_SRC
python compare.py                    # write REPORT.md and figures/*.png
```

The paths assume clones of the Python and MATLAB repositories next to this
one, named `replay_python` and `replay_matlab`.  Set
`UWA_CALIBRATION_CASES=blue_1` to run one case while you iterate.
`compare.py` exits non-zero when a check is outside its tolerance.

## Rules

- `config.json` holds every parameter, and both languages read it.  Do not copy
  a value from it into Python or MATLAB code.
- Do not loosen a tolerance to make a run pass.  A run outside tolerance is a
  finding; report it.
- Do not edit `REPORT.md` or `figures/` by hand.  The workflow
  (`.github/workflows/calibrate.yml`) writes them and commits them to `main`.
  Pull before you push, because that workflow commits too.
- Do not commit `*.mat` files or anything under `artifacts/`.
- To add a channel, add one entry to `config.json`.  No code changes are
  needed.
- Commit messages are one imperative sentence that says what changed and why.

## Citation

`CITATION.cff` cites the arXiv paper (arXiv:2609.03207) and the Zenodo data
(doi:10.5281/zenodo.21287414).
