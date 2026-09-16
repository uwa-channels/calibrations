"""Compare the two implementations and write artifacts/report.md.

Replay is compared sample by sample: the two runners were handed the same
probe, the same channel file, the same hydrophones and the same start index, so
the outputs should agree to within the numerical floor.  Noise cannot be
compared that way -- NumPy's and MATLAB's generators produce different streams
from the same seed -- so it is compared as a distribution, a spectrum and a
spatial covariance, each also checked against what ``beta`` predicts, which is
an absolute reference rather than a cross-check.

Exits non-zero if any comparison exceeds its tolerance in config.json.

Author: Zhengnan Li
Email : uwa-channels@ofdm.link
License: MIT
"""

import sys

import numpy as np
import scipy.signal as sg
from scipy.io import loadmat
from scipy.stats import kurtosis, ks_2samp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from calib import ARTIFACTS, FIGURES, config  # noqa: E402

CHECKS = []      # (name, value, tolerance, ok, units)
LINES = []       # report body


def check(name, value, limit, ok, units=""):
    CHECKS.append((name, value, limit, ok, units))
    return ok


def say(line=""):
    LINES.append(line)
    print(line)


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def nmse_db(a, b):
    """Error of ``b`` relative to ``a``, in dB."""
    return 10 * np.log10(np.sum(np.abs(a - b) ** 2) / np.sum(np.abs(a) ** 2))


def best_lag(a, b, max_lag=64):
    """Integer lag of ``b`` relative to ``a`` that maximizes correlation."""
    a = a - a.mean()
    b = b - b.mean()
    c = sg.correlate(b, a, mode="full")
    lags = sg.correlation_lags(len(b), len(a), mode="full")
    keep = np.abs(lags) <= max_lag
    return int(lags[keep][np.argmax(np.abs(c[keep]))])


def trim(a, b, margin):
    """Common interior of two signals, dropping ``margin`` samples each end.

    The two outputs differ in length (the MATLAB replay carries 20 extra
    samples of extrapolation buffer before the final resampling), and a
    rational resampler's output near either edge depends on that padding, so
    the ends are not a like-for-like comparison.
    """
    n = min(len(a), len(b))
    return a[margin:n - margin], b[margin:n - margin]


# ---------------------------------------------------------------------------
# Replay
# ---------------------------------------------------------------------------

def compare_replay(cfg):
    py = loadmat(ARTIFACTS / "python_replay.mat")
    ml = loadmat(ARTIFACTS / "matlab_replay.mat")
    tol = cfg["tolerances"]
    idx = cfg["replay"]["array_index"]
    fs = cfg["replay"]["fs"]
    margin = 200

    def worst(key_py, key_ml):
        return max(nmse_db(*trim(py[key_py][:, m], ml[key_ml][:, m], margin))
                   for m in range(py["y"].shape[1]))

    say("## Replay")
    say()
    say(f"Probe: {len(py['y_resamp'])} samples at {fs / 1e3:g} kHz through "
        f"`{cfg['data']['channel_file']}`, Python start index "
        f"{int(py['start'].ravel()[0])} and MATLAB start index "
        f"{int(py['start'].ravel()[0]) + 1}, hydrophones {idx} (0-based).  "
        "Both implementations were handed the same probe samples from "
        "`probe.mat`, so nothing here depends on a random number generator.")
    say()
    say(f"Python returns {py['y'].shape[0]} samples, MATLAB "
        f"{ml['y'].shape[0]}: MATLAB allocates `T + buffer + L` with "
        f"`buffer = 20` where Python allocates `T + L`, which after the "
        f"closing resample by q/p is {ml['y'].shape[0] - py['y'].shape[0]} "
        f"extra samples.  The comparison uses the common interior, dropping "
        f"{margin} samples at each end.")
    say()

    say("| what is compared | NMSE | reading |")
    say("|---|---|---|")

    a, b = trim(py["y_resamp"].ravel(), ml["y_resamp"].ravel(), margin)
    e_res = nmse_db(a, b)
    say(f"| resampler alone, no channel | **{e_res:.0f} dB** | "
        "`resample` and `resample_poly` build the same filter |")
    check("resampler control", e_res, tol["replay_resampler_nmse_db"],
          e_res <= tol["replay_resampler_nmse_db"], "dB")

    e_static = worst("y_static", "y_static")
    say(f"| full replay, `phi_hat` removed | **{e_static:.0f} dB** | "
        "`h_hat` layout, spline interpolation, time-varying convolution, both "
        "resamplings and the up-conversion agree to the last bit |")
    check("replay with no tracking", e_static, tol["replay_static_nmse_db"],
          e_static <= tol["replay_static_nmse_db"], "dB")

    e_phi = worst("y", "y")
    say(f"| full replay, `phi_hat` active | **{e_phi:.0f} dB** | "
        "the delay/phase trajectory agrees too |")
    check("replay with delay tracking", e_phi, tol["replay_nmse_db"],
          e_phi <= tol["replay_nmse_db"], "dB")

    e_off = worst("y", "y_same_start")
    say(f"| ditto, MATLAB given `start` instead of `start+1` | "
        f"**{e_off:.0f} dB** | deliberately one sample out, and it shows |")
    check("one-sample offset is detected", -e_off, -tol["replay_offset_probe_db"],
          e_off > tol["replay_offset_probe_db"], "dB")
    say()

    say("The last row is the harness testing itself.  A comparison that cannot "
        "see a one-sample offset in `phi_hat` would pass everything, including "
        "the defect this suite was built to find: before "
        "`replay.m` was corrected on Sep. 16, 2026, the aligned row above read "
        "-63 dB rather than "
        f"{e_phi:.0f} dB, because MATLAB interpolated `h_hat` onto "
        "`(start + 0:N-1)/fs_delay` while taking the phase from "
        "`phi_hat(start : start+N-1)`, whose samples sit one sample of "
        "`fs_delay` earlier.  `unpack.m` had the origin right; `replay.m` did "
        "not.")
    say()

    say("### Per hydrophone")
    say()
    say("| hydrophone | NMSE, aligned | lag | NMSE, one sample out |")
    say("|---|---|---|---|")
    for m in range(py["y"].shape[1]):
        a1, b1 = trim(py["y"][:, m], ml["y"][:, m], margin)
        a2, b2 = trim(py["y"][:, m], ml["y_same_start"][:, m], margin)
        say(f"| {idx[m]} | {nmse_db(a1, b1):.0f} dB | {best_lag(a1, b1):+d} | "
            f"{nmse_db(a2, b2):.1f} dB |")
        lag = abs(best_lag(a1, b1))
        check(f"replay lag, hydrophone {idx[m]}", lag, tol["replay_lag_samples"],
              lag <= tol["replay_lag_samples"], "samples")
    say()

    # Figures
    a, b = trim(py["y"][:, 0], ml["y"][:, 0], margin)
    ao, bo = trim(py["y"][:, 0], ml["y_same_start"][:, 0], margin)
    n0 = len(a) // 2
    seg = slice(n0, n0 + 400)
    fig, ax = plt.subplots(3, 1, figsize=(9, 9), constrained_layout=True)
    ax[0].plot(a[seg], label="Python", lw=1.2)
    ax[0].plot(b[seg], "--", label="MATLAB", lw=1.2)
    ax[0].set_title(f"Replay, hydrophone {idx[0]}, 400 samples mid-record")
    ax[0].set_xlabel("Sample"); ax[0].legend(); ax[0].grid(alpha=0.3)

    ref = 20 * np.log10(np.abs(a).max())
    ax[1].plot(20 * np.log10(np.abs(ao - bo) + 1e-300) - ref, lw=0.6,
               label="MATLAB one sample out")
    ax[1].plot(20 * np.log10(np.abs(a - b) + 1e-300) - ref, lw=0.6,
               label="aligned")
    ax[1].set_ylim(-320, 0)
    ax[1].set_title("Python minus MATLAB, dB relative to peak")
    ax[1].set_xlabel("Sample"); ax[1].legend(); ax[1].grid(alpha=0.3)

    f, Pa = sg.welch(a, fs=fs, nperseg=4096)
    _, Pb = sg.welch(b, fs=fs, nperseg=4096)
    ax[2].semilogy(f / 1e3, Pa, label="Python")
    ax[2].semilogy(f / 1e3, Pb, "--", label="MATLAB")
    ax[2].set_xlim(5, 21); ax[2].set_xlabel("Frequency [kHz]")
    ax[2].set_title("Replay spectrum"); ax[2].legend(); ax[2].grid(alpha=0.3)
    fig.savefig(FIGURES / "replay.png", dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Noise
# ---------------------------------------------------------------------------

def theory(beta, Fs, nfreq):
    """Covariance and one-sided PSD implied by the mixing coefficients.

    ``w(n, i) = sum_j sum_k beta(i, j, k) z(n + k, j)`` with ``z`` iid of unit
    variance gives ``C = sum_k B_k B_k^T`` and, per output channel,
    ``S_i(f) = (2 / Fs) sum_j |sum_k beta(i, j, k) e^(-2 pi i f k / Fs)|^2``.
    """
    M, _, K = beta.shape
    C = np.einsum("ijk,ljk->il", beta, beta)
    f = np.linspace(0, Fs / 2, nfreq)
    phase = np.exp(-2j * np.pi * np.outer(f, np.arange(K)) / Fs)   # (nf, K)
    H = np.einsum("ijk,fk->ijf", beta, phase)                      # (M, M, nf)
    S = (2.0 / Fs) * np.sum(np.abs(H) ** 2, axis=1)                # (M, nf)
    return C, f, S


def coherence(C):
    d = np.sqrt(np.diag(C))
    return np.abs(C / np.outer(d, d))


def compare_noise(cfg):
    py = loadmat(ARTIFACTS / "python_noise.mat")
    ml = loadmat(ARTIFACTS / "matlab_noise.mat")
    tol = cfg["tolerances"]
    nper = cfg["noise"]["welch_nperseg"]

    wp, wm = py["w"], ml["w"]
    beta = np.asarray(py["beta"], dtype=float)
    Fs = float(np.asarray(py["Fs"]).ravel()[0])
    fs = float(cfg["noise"]["fs"])
    assert np.allclose(beta, np.asarray(ml["beta"], dtype=float)), \
        "the two runners disagree about beta itself"

    say("## Noise")
    say()
    say(f"{wp.shape[0]} samples on {wp.shape[1]} hydrophones per implementation "
        f"({wp.shape[0] / fs:.1f} s at {fs / 1e3:g} kHz), alpha = "
        f"{float(np.asarray(py['alpha']).ravel()[0]):g}, beta {beta.shape}.")
    say()
    say("The two draw from different pseudo-random streams, so they are "
        "compared as distributions rather than sample by sample, and each is "
        "also compared against what `beta` predicts.")
    say()

    # --- distribution ---
    sp, sm = wp.std(axis=0), wm.std(axis=0)
    kp, km = kurtosis(wp, axis=0), kurtosis(wm, axis=0)
    rel_std = np.abs(sp - sm) / sp
    say("### Distribution")
    say()
    say("| hydrophone | std (Py) | std (ML) | rel. diff | excess kurtosis (Py) | (ML) | KS stat |")
    say("|---|---|---|---|---|---|---|")
    ks = []
    rng = np.random.default_rng(0)
    for m in range(wp.shape[1]):
        sub = rng.choice(wp.shape[0], size=min(60000, wp.shape[0]), replace=False)
        d = ks_2samp(wp[sub, m], wm[sub, m]).statistic
        ks.append(d)
        say(f"| {m} | {sp[m]:.4g} | {sm[m]:.4g} | {rel_std[m]:.2%} | "
            f"{kp[m]:+.3f} | {km[m]:+.3f} | {d:.4f} |")
    say()
    check("noise std, max relative difference", rel_std.max(),
          tol["noise_std_rel"], rel_std.max() <= tol["noise_std_rel"], "")
    alpha = float(np.asarray(py["alpha"]).ravel()[0])
    if alpha == 2:
        # For alpha = 2 the driver is Gaussian and the mixing is linear, so the
        # excess kurtosis of both outputs must be zero.  That is an absolute
        # reference; the pairwise difference below is only a cross-check.
        check("excess kurtosis vs Gaussian (Python)", np.abs(kp).max(),
              tol["noise_kurtosis_abs"], np.abs(kp).max() <= tol["noise_kurtosis_abs"], "")
        check("excess kurtosis vs Gaussian (MATLAB)", np.abs(km).max(),
              tol["noise_kurtosis_abs"], np.abs(km).max() <= tol["noise_kurtosis_abs"], "")
    check("noise excess kurtosis, max |difference|", np.abs(kp - km).max(),
          tol["noise_kurtosis_abs"], np.abs(kp - km).max() <= tol["noise_kurtosis_abs"], "")
    check("noise KS statistic, max", max(ks), tol["noise_ks_stat"],
          max(ks) <= tol["noise_ks_stat"], "")

    # --- spectrum ---
    f, Pp = sg.welch(wp, fs=fs, nperseg=nper, axis=0)
    _, Pm = sg.welch(wm, fs=fs, nperseg=nper, axis=0)
    Pp, Pm = Pp.T, Pm.T
    band = Pp.mean(axis=0) > Pp.max() * 1e-3          # where there is power
    dpsd = 10 * np.log10(Pp[:, band] / Pm[:, band])
    say("### Spectrum")
    say()
    say(f"Welch, nperseg {nper}, {int(band.sum())} bins inside the 30 dB band.")
    say(f"Python vs MATLAB: rms difference **{np.sqrt((dpsd ** 2).mean()):.3f} dB**, "
        f"max |difference| {np.abs(dpsd).max():.2f} dB.")
    check("noise PSD, rms Python-vs-MATLAB difference", np.sqrt((dpsd ** 2).mean()),
          tol["noise_psd_rms_db"], np.sqrt((dpsd ** 2).mean()) <= tol["noise_psd_rms_db"], "dB")

    C_th, f_th, S_th = theory(beta, Fs, len(f))
    dpy = 10 * np.log10(Pp[:, band] / S_th[:, band])
    dml = 10 * np.log10(Pm[:, band] / S_th[:, band])
    say(f"Against the spectrum `beta` implies: Python {dpy.mean():+.2f} dB mean "
        f"({np.sqrt((dpy ** 2).mean()):.2f} dB rms), MATLAB {dml.mean():+.2f} dB "
        f"({np.sqrt((dml ** 2).mean()):.2f} dB rms).")
    say()

    # --- spatial correlation ---
    Cp = np.cov(wp, rowvar=False)
    Cm = np.cov(wm, rowvar=False)
    gp, gm, gth = coherence(Cp), coherence(Cm), coherence(C_th)
    off = ~np.eye(len(gp), dtype=bool)
    rms_xy = np.sqrt(np.mean((gp[off] - gm[off]) ** 2))
    rms_pth = np.sqrt(np.mean((gp[off] - gth[off]) ** 2))
    rms_mth = np.sqrt(np.mean((gm[off] - gth[off]) ** 2))
    # The transposed convention, the bug this whole exercise exists to catch.
    gT = coherence(np.einsum("jik,ljk->il", beta, beta))
    rms_pT = np.sqrt(np.mean((gp[off] - gT[off]) ** 2))
    rms_mT = np.sqrt(np.mean((gm[off] - gT[off]) ** 2))

    say("### Spatial correlation")
    say()
    say("| comparison | rms difference in coherence |")
    say("|---|---|")
    say(f"| Python vs MATLAB | {rms_xy:.4f} |")
    say(f"| Python vs theory from `beta` | {rms_pth:.4f} |")
    say(f"| MATLAB vs theory from `beta` | {rms_mth:.4f} |")
    say(f"| Python vs *transposed* `beta` | {rms_pT:.4f} |")
    say(f"| MATLAB vs *transposed* `beta` | {rms_mT:.4f} |")
    say()
    say("The last two rows are the control: `beta` is not symmetric in (i, j), "
        "so an implementation that mixes `sum_j beta_ji z_j` instead of "
        "`sum_j beta_ij z_j` lands on the transposed covariance.  Both "
        "implementations must be far closer to theory than to its transpose.")
    say()
    check("noise coherence, Python vs MATLAB", rms_xy, tol["noise_coherence_rms"],
          rms_xy <= tol["noise_coherence_rms"], "")
    check("noise coherence, Python vs theory", rms_pth, tol["noise_theory_coherence_rms"],
          rms_pth <= tol["noise_theory_coherence_rms"], "")
    check("noise coherence, MATLAB vs theory", rms_mth, tol["noise_theory_coherence_rms"],
          rms_mth <= tol["noise_theory_coherence_rms"], "")
    check("noise coherence closer to beta than to beta-transposed (Python)",
          rms_pth, rms_pT, rms_pth < rms_pT, "")
    check("noise coherence closer to beta than to beta-transposed (MATLAB)",
          rms_mth, rms_mT, rms_mth < rms_mT, "")

    # Figures
    fig, ax = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    ax[0, 0].semilogy(f / 1e3, Pp[0], label="Python")
    ax[0, 0].semilogy(f / 1e3, Pm[0], "--", label="MATLAB")
    ax[0, 0].semilogy(f_th / 1e3, S_th[0], ":", color="k", label="from beta")
    ax[0, 0].set_title("Noise PSD, hydrophone 0")
    ax[0, 0].set_xlabel("Frequency [kHz]"); ax[0, 0].legend(); ax[0, 0].grid(alpha=0.3)

    edges = np.linspace(-4 * sp[0], 4 * sp[0], 160)
    ax[0, 1].hist(wp[:, 0], bins=edges, density=True, histtype="step", label="Python")
    ax[0, 1].hist(wm[:, 0], bins=edges, density=True, histtype="step", label="MATLAB")
    ax[0, 1].set_yscale("log"); ax[0, 1].set_title("Amplitude distribution, hydrophone 0")
    ax[0, 1].legend(); ax[0, 1].grid(alpha=0.3)

    for a_, g_, t_ in ((ax[1, 0], gp, "Python coherence"),
                       (ax[1, 1], gp - gth, "Python minus theory")):
        im = a_.imshow(g_, cmap="viridis" if "minus" not in t_ else "coolwarm")
        a_.set_title(t_); fig.colorbar(im, ax=a_)
    fig.savefig(FIGURES / "noise.png", dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------------------

def main():
    cfg = config()
    FIGURES.mkdir(parents=True, exist_ok=True)
    say("# Cross-implementation calibration report")
    say()
    say(f"`{cfg['data']['channel_file']}` and `{cfg['data']['noise_file']}` "
        f"from Zenodo record [{cfg['data']['zenodo_record']}]"
        f"(https://doi.org/{cfg['data']['zenodo_doi']}).")
    say()
    compare_replay(cfg)
    compare_noise(cfg)

    say("## Verdict")
    say()
    say("| check | value | tolerance | |")
    say("|---|---|---|---|")
    failed = 0
    for name, value, limit, ok, units in CHECKS:
        failed += not ok
        say(f"| {name} | {value:.4g} {units} | {limit:.4g} {units} | "
            f"{'pass' if ok else '**FAIL**'} |")
    say()
    say(f"{len(CHECKS) - failed} of {len(CHECKS)} checks pass.")

    (ARTIFACTS / "report.md").write_text("\n".join(LINES) + "\n")
    print(f"\nwrote {ARTIFACTS / 'report.md'} and {FIGURES}/*.png")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

# [EOF]
