# Cross-implementation calibration

Checks that [`uwa-channels/python`](https://github.com/uwa-channels/python) and
[`uwa-channels/matlab`](https://github.com/uwa-channels/matlab) do the same
thing to the same data.  Both replay the same probe through the same channel
file and generate noise from the same mixing coefficients; this repository runs
them side by side, on files fetched from the published Zenodo record, and
reports where they agree and where they do not.

Each suite already tests its own implementation against synthetic fixtures.  A
fixture is written by the same conventions its reader assumes, so the two
packages can both pass and still disagree with each other -- and one of them
can be wrong about the released files without any test noticing.  That is not
hypothetical: it is how a transposed `beta` shipped in the Python package.

## What it found

Two off-by-one errors in `uwa-channels/matlab`, both fixed on Sep. 16, 2026,
both invisible to either package's own test suite.

### `replay`: the phase trajectory ran one sample ahead of the impulse response

| what is compared | before | after |
|---|---|---|
| the rational resampler alone, no channel | -258 dB | -258 dB |
| full replay with tracking removed | -234 dB | -234 dB |
| full replay with `phi_hat` active | **-63 dB** | **-233 dB** |

The first two rows said `h_hat`'s axis layout, the spline interpolation, the
time-varying convolution, both resamplings and the up-conversion were identical
to the last bit -- MATLAB's `resample` and SciPy's `resample_poly` turn out to
build the same filter, so there was no resampler floor to hide behind.  The
third put the disagreement in the delay trajectory and nowhere else:

> With `start = s`, both implementations interpolated `h_hat` onto
> `(s + 0:N-1) / fs_delay`.  MATLAB then took the phase from
> `phi_hat(s : s+N-1)` -- 1-based, so elements *s* through *s+N-1*, which sit
> at `(s-1 : s+N-2) / fs_delay`.  Every output sample therefore carried a
> phase one sample of `fs_delay` too early.

### `unpack`: the `f_resamp` ramp started one sample in

| what is compared | before | after |
|---|---|---|
| unpack, no `f_resamp` | -291 dB | -291 dB |
| unpack, `f_resamp` active | **-44 dB** | **-260 dB** |

`unpack` built the ramp on `(1:N_phi)` where `t_orig`, the grid it is then
interpolated from, is `(0:N_phi-1)/fs_delay`.  That put a constant phase
rotation -- measured at +0.00607 rad against a predicted one-sample step of
+0.00636 rad -- and through `phase_drift` a constant delay offset, on every
unpacked tap.

Both bugs were MATLAB's, and in both cases MATLAB's own code contained the
evidence: `unpack.m` already used the right origin for `replay.m`'s bug, and
`t_orig` two lines below the ramp already used it for `unpack.m`'s.

### Why neither suite caught them

`testUnpack.m` has no assertions at all -- it runs `unpack` and draws
pictures -- and its `f_resamp` cases were among them.  `test_unpack.py` is
mostly the same.  The replay suites do assert, but they build their fixtures
by the same convention they then check, so a shared origin error cancels and no
assertion moves.  Both repositories now carry a test that fails on the old
code: `testFResampRampOrigin` and `test_f_resamp_ramp_starts_at_zero`, each
comparing a channel unpacked with and without `f_resamp` at the first output
sample, where the ramp must contribute nothing.

### Everything else agrees

Noise: KS statistic below 0.01, excess kurtosis within 0.03 of Gaussian on both
sides, spectra 0.40 dB rms apart against 0.39 dB of Welch sampling error,
coherence 0.003 rms.  Both land on the covariance `beta` predicts and nowhere
near its transpose.

The harness keeps two probes on itself: a deliberately misaligned replay
(MATLAB given `start` rather than `start+1`, which reads -60 dB) and a check
that `f_resamp` changes the unpacked output at all.  A comparison that could
not see a one-sample offset, or that passed because a code path quietly did
nothing, would pass everything.

Still unfixed, and reported rather than changed: MATLAB's `replay` returns 98
more samples than Python's for the same input -- it allocates `T + buffer + L`
with `buffer = 20` where Python allocates `T + L`.  Harmless, but a caller
porting between the two will see it.

## What it covers

Two cases, listed in `config.json`, chosen so that between them they take every
branch:

| case | tracking | array | exercises |
|---|---|---|---|
| `blue_1.mat` + `blue_noise.mat` | `phi_hat` | 12 | replay's delay-tracking branch, all 12 hydrophones of noise |
| `purple_3.mat` + `purple_noise_3.mat` | `theta_hat` | 24 | replay's phase-only branch, and noise on 12 of 24 hydrophones, so the `beta` reference has to be subset the same way |

Each case is run through `replay` (with tracking, with tracking removed, and
with the resampler alone), `unpack` (with and without `f_resamp`), and
`noisegen`.  `f_resamp` is attached by both runners from `config.json`: the
released files do not carry one, and the path needs exercising -- it is where
the second off-by-one lived.

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

**Every difference is isolated before it is judged.**  A single end-to-end NMSE
would have said "-63 dB, good enough" and hidden a real bug.  Removing
`phi_hat`, then the channel, then emulating the other convention, is what turns
a number into a diagnosis.

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
