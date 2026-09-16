# Cross-implementation calibration

Checks that [`uwa-channels/python`](https://github.com/uwa-channels/python) and
[`uwa-channels/matlab`](https://github.com/uwa-channels/matlab) do the same
thing to the same data.  Both replay the same probe through the same channel
file and generate noise from the same mixing coefficients; this repository runs
them side by side, on files fetched from the published Zenodo record, and
reports where they agree and where they do not.

Each suite already tests its own implementation against synthetic fixtures.  A
fixture is written by the same conventions its reader assumes, so a convention
the two packages disagree about cancels inside each of them: both suites pass,
the outputs differ, and neither suite is positioned to notice.  Checking them
against each other on the files a user actually downloads is what closes that
gap.

## What it reports

Both implementations, on the released files, at the numerical floor:

| | `blue_1` | `purple_3` |
|---|---|---|
| resampler alone, no channel | -258 dB | -260 dB |
| replay, tracking removed | -234 dB | -232 dB |
| replay, tracking active | -233 dB | -232 dB |
| unpack, no `f_resamp` | -291 dB | -298 dB |
| unpack, `f_resamp` active | -260 dB | -264 dB |

Noise cannot be compared sample by sample, so it is compared as a
distribution, a spectrum and a spatial covariance:

| | `blue_1` | `purple_3` |
|---|---|---|
| std, max relative difference | 0.75% | 0.67% |
| excess kurtosis vs Gaussian, worst side | 0.030 | 0.020 |
| KS statistic, max | 0.0084 | 0.0068 |
| PSD, rms difference | 0.40 dB | 0.40 dB |
| coherence vs `beta`, rms off-diagonal | 0.0021 | 0.0023 |
| coherence vs *transposed* `beta` | 0.1402 | 0.1075 |

The PSD figure is at the limit of what the estimator can resolve: Welch on
2^19 samples gives two independent estimates 0.38 dB apart on average, so
0.40 dB is agreement, not disagreement.  The last row is the control -- `beta`
is not symmetric in (i, j), so an implementation mixing `sum_j beta_ji z_j`
instead of `sum_j beta_ij z_j` lands on the transposed covariance, two orders
of magnitude away.

Two checks measure the harness rather than the implementations.  One replays
with MATLAB given `start` rather than `start+1`, which is a one-sample
misalignment and reads -58 to -60 dB; the other confirms `f_resamp` changes the
unpacked output at all.  Without them, a comparison too blunt to resolve a
single sample, or one passing because a code path did nothing, would report
success.

One difference is reported rather than flagged: MATLAB's `replay` returns 98
more samples than Python's for the same input, because it allocates
`T + buffer + L` with `buffer = 20` where Python allocates `T + L`.  It
changes no sample they have in common.

## What it covers

Two cases, listed in `config.json`, chosen so that between them they take every
branch:

| case | tracking | array | exercises |
|---|---|---|---|
| `blue_1.mat` + `blue_noise.mat` | `phi_hat` | 12 | replay's delay-tracking branch, all 12 hydrophones of noise |
| `purple_3.mat` + `purple_noise_3.mat` | `theta_hat` | 24 | replay's phase-only branch, and noise on 12 of 24 hydrophones, so the `beta` reference has to be subset the same way |

Each case is run through `replay` (with tracking, with tracking removed, and
with the resampler alone), `unpack` (with and without `f_resamp`), and
`noisegen`.  `f_resamp` is attached by both runners from `config.json`, since the released
files do not carry one and the path would otherwise go untested.

Adding a channel is one entry in `config.json`; nothing in the code changes.

## Design

**The probe is generated once, in Python, and written to `artifacts/probe.mat`
for both runners to read.**  NumPy's Mersenne Twister and MATLAB's produce
different streams from the same seed, so "same seed" is not "same signal".
Nothing in the replay comparison calls a random number generator.

**Noise is compared statistically, because it cannot be compared any other
way.**  The two drivers are different streams by construction.  So the noise
checks are distribution, power spectrum and spatial covariance -- and each is
*also* checked against what `beta` predicts in closed form,

```
C      = sum_k B_k B_k^T
S_i(f) = (2 / Fs) sum_j |sum_k beta(i, j, k) exp(-2 pi i f k / Fs)|^2
```

which is an absolute reference rather than a cross-check.  Two implementations
that made the same mistake would still fail against it.  The transposed
covariance is computed too, as the control: `beta` is not symmetric, so an
implementation mixing `sum_j beta_ji z_j` lands on it.

**Every difference is isolated before it is judged.**  A single end-to-end
number says only that two outputs differ, which is the least useful thing to
know: it cannot separate a resampler from an interpolator from an index.  So
each comparison is run as a ladder -- the resampler alone, then the channel
with tracking removed, then tracking active, then `unpack` with and without
`f_resamp` -- and a difference is attributed to the rung it first appears on.

## Running it

```bash
pip install numpy scipy matplotlib h5py
python calib.py                      # fetch the library files from Zenodo (740 MB)
python make_probe.py                 # writes artifacts/probe.mat

PYTHONPATH=../replay_python/src python run_python.py
matlab -batch run_matlab             # picks up ../replay_matlab/src
python compare.py                    # writes artifacts/report.md, figures/*.png
```

`UWA_CALIBRATION_CASES=blue_1` restricts every step to a subset while
iterating, which both runners and the comparer honour.

`compare.py` exits non-zero if any check exceeds its tolerance.  Override the
locations with `UWA_CHANNELS_CACHE` (data files), `UWA_PYTHON_SRC` via
`PYTHONPATH`, and `UWA_MATLAB_SRC` (MATLAB `src/`).

Everything about the run -- which files, which hydrophones, which start index,
which tolerances -- is in `config.json`, read by both languages.  Neither owns
a second copy.

## Setting the tolerances

The tolerances in `config.json` are not guesses, and they are not fitted to
whatever the last run produced.  The bit-exact ones are set at -200 dB, far
above the measured -234 dB, because they are claims about floating point rather
than about the algorithm.  The statistical ones come from the sampling error of
the estimator:

- **Spectrum.**  Welch on 2^19 samples with `nperseg` 4096 averages 255
  segments, so a single PSD estimate has a standard deviation of
  `10/ln(10) * sqrt(1/255)` = 0.27 dB, and the difference of two independent
  estimates 0.38 dB rms.  Measured: 0.40 dB.  Tolerance 0.8 dB.
- **Kurtosis.**  The mixing is 129 taps long, so the 2^19 samples carry roughly
  `2^19 / 129` independent ones and the sample excess kurtosis has a spread
  near `sqrt(24 * 129 / 2^19)` = 0.08.  Measured: within 0.03.  Tolerance 0.15.
- **Coherence.**  Same effective sample count gives a few per cent on each
  off-diagonal entry.  Measured: 0.003 rms.  Tolerance 0.03.

Replay and unpack have no sampling error to allow for -- same input, same
file, same parameters -- so they are held at -200 dB and measure between -232
and -298 dB.

If a tolerance ever needs loosening to make a run pass, that is a finding, not
maintenance.

## CI

`.github/workflows/calibrate.yml` runs the whole thing weekly and on demand,
pulling both implementations from their default branches so that a change to
either is caught against the other.  It caches the Zenodo download between
runs; the first run fetches 740 MB, later ones nothing.

MATLAB runs via `matlab-actions/setup-matlab`.  On GitHub-hosted runners this
needs no license for a public repository; a private one needs a
`MATLAB_BATCH_LICENSE_TOKEN` secret.  The MATLAB job is allowed to be skipped
rather than fail the workflow when no license is available, so a fork still
gets the Python-side checks against `beta`.

## Files

| file | what it does |
|---|---|
| `config.json` | every parameter, read by both languages |
| `calib.py` | paths, the case list, and the Zenodo fetch with MD5 checks |
| `make_probe.py` | builds the shared probe |
| `run_python.py` | runs the Python implementation |
| `run_matlab.m` | runs the MATLAB implementation |
| `compare.py` | metrics, figures, `artifacts/report.md`, exit status |

Each case in `config.json` names a channel and a noise file and is run through
`replay`, `unpack` (with and without `f_resamp`) and `noisegen`.  Adding a
channel is one entry; nothing in the code changes.
