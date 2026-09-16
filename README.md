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

## What it finds today

Replay, on `blue_1.mat`, hydrophones 0/3/7, the same probe and start index:

| what is compared | NMSE |
|---|---|
| the rational resampler alone, no channel | -258 dB |
| full replay with `phi_hat` removed | -234 dB |
| full replay with `phi_hat` active | **-63 dB** |
| ditto, with Python delaying `phi_hat` one sample | -233 dB |

The first two rows say that `h_hat`'s axis layout, the spline interpolation,
the time-varying convolution, both resamplings and the up-conversion are
identical to the last bit -- MATLAB's `resample` and SciPy's `resample_poly`
turn out to build the same filter.  The third row says the delay/phase
trajectory is not.  The fourth says exactly why:

> With `start = s`, both implementations interpolate `h_hat` onto
> `(s + 0:N-1) / fs_delay`.  MATLAB then takes the phase from
> `phi_hat(s : s+N-1)` -- 1-based, so elements *s* through *s+N-1* -- while
> Python takes `phi_hat[s : s+N]` -- 0-based, so elements *s+1* onward.
> Element *n* of `phi_hat` belongs at time `(n-1) / fs_delay`, so Python's
> pairing is the self-consistent one and **MATLAB's phase runs one sample of
> `fs_delay` ahead of the impulse response it multiplies**.

Delaying `phi_hat` by one sample on the Python side reproduces the MATLAB
output to the numerical floor, so that offset is the whole of the difference,
not merely most of it.  It is worth about 0.07% rms on this channel, because
`h_hat` changes slowly against `fs_time`.  The fix belongs in `replay.m`.

The MATLAB output is also 98 samples longer for the same input: it allocates
`T + buffer + L` with `buffer = 20` where Python allocates `T + L`.  Harmless,
but a caller porting between the two will see it.

Noise agrees on every measure: distribution (KS statistic < 0.01, excess
kurtosis within 0.03 of Gaussian on both sides), spectrum (0.40 dB rms, against
0.39 dB of Welch sampling error), and spatial coherence (0.003 rms).  Both land
on the covariance `beta` predicts and nowhere near its transpose.

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
python calib.py                      # fetch blue_1.mat and blue_noise.mat from Zenodo
python make_probe.py                 # writes artifacts/probe.mat

PYTHONPATH=../replay_python/src python run_python.py
matlab -batch run_matlab             # picks up ../replay_matlab/src
python compare.py                    # writes artifacts/report.md, figures/*.png
```

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

If a tolerance ever needs loosening to make a run pass, that is a finding, not
maintenance.

## CI

`.github/workflows/calibrate.yml` runs the whole thing weekly and on demand,
pulling both implementations from their default branches so that a change to
either is caught against the other.  It caches the Zenodo download between
runs.

MATLAB runs via `matlab-actions/setup-matlab`.  On GitHub-hosted runners this
needs no license for a public repository; a private one needs a
`MATLAB_BATCH_LICENSE_TOKEN` secret.  The MATLAB job is allowed to be skipped
rather than fail the workflow when no license is available, so a fork still
gets the Python-side checks against `beta`.

## Files

| file | what it does |
|---|---|
| `config.json` | every parameter, read by both languages |
| `calib.py` | paths and the Zenodo fetch, with MD5 checks |
| `make_probe.py` | builds the shared probe |
| `run_python.py` | runs the Python implementation |
| `run_matlab.m` | runs the MATLAB implementation |
| `compare.py` | metrics, figures, `artifacts/report.md`, exit status |
