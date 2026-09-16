"""Compare the two implementations and write artifacts/report.md.

Replay and unpack are compared sample by sample: both runners were handed the
same probe, the same channel file, the same hydrophones and the same start
index, so the outputs should agree to the numerical floor.  Noise cannot be
compared that way -- NumPy's and MATLAB's generators produce different streams
from the same seed -- so it is compared as a distribution, a spectrum and a
spatial covariance, each also checked against what ``beta`` predicts in closed
form, which is an absolute reference rather than a cross-check.

Exits non-zero if any comparison exceeds its tolerance in config.json.

Author: Zhengnan Li
Email : uwa-channels@ofdm.link
License: MIT
"""

import sys

import h5py
import numpy as np
import scipy.signal as sg
from scipy.io import loadmat
from scipy.stats import kurtosis, ks_2samp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from calib import ARTIFACTS, FIGURES, cases, config  # noqa: E402

CHECKS = []
LINES = []


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
    den = np.sum(np.abs(a) ** 2)
    return 10 * np.log10(np.sum(np.abs(a - b) ** 2) / den) if den else np.nan


def best_lag(a, b, max_lag=64):
    a, b = a - a.mean(), b - b.mean()
    c = sg.correlate(b, a, mode="full")
    lags = sg.correlation_lags(len(b), len(a), mode="full")
    keep = np.abs(lags) <= max_lag
    return int(lags[keep][np.argmax(np.abs(c[keep]))])


def trim(a, b, margin):
    """Common interior of two signals, dropping ``margin`` samples each end."""
    n = min(len(a), len(b))
    return a[margin:n - margin], b[margin:n - margin]


def load_v73_complex(path, keys):
    """Read complex arrays MATLAB wrote with ``-v7.3``, in MATLAB's layout.

    HDF5 reverses every axis and splits complex into a compound dtype, so a
    MATLAB (K, M, T) array arrives as (T, M, K) with 'real'/'imag' fields.
    """
    out = {}
    with h5py.File(path, "r") as f:
        for k in keys:
            d = f[k]
            z = d["real"][...] + 1j * d["imag"][...]
            out[k] = np.transpose(z, tuple(reversed(range(z.ndim))))
    return out


# ---------------------------------------------------------------------------
# Replay
# ---------------------------------------------------------------------------

def compare_replay(case, cfg):
    name = case["name"]
    py = loadmat(ARTIFACTS / f"python_{name}_replay.mat")
    ml = loadmat(ARTIFACTS / f"matlab_{name}_replay.mat")
    tol = cfg["tolerances"]
    idx = case["replay"]["array_index"]
    fs = cfg["probe"]["fs"]
    margin = 200

    def worst(kp, km):
        return max(nmse_db(*trim(py[kp][:, m], ml[km][:, m], margin))
                   for m in range(py["y"].shape[1]))

    say(f"### Replay ({name})")
    say()
    say(f"Python start {int(py['start'].ravel()[0])}, MATLAB start "
        f"{int(py['start'].ravel()[0]) + 1}, hydrophones {idx} (0-based).  "
        f"Python returns {py['y'].shape[0]} samples, MATLAB {ml['y'].shape[0]} "
        f"(MATLAB's `buffer = 20`, resampled by q/p).  The comparison uses the "
        f"common interior, dropping {margin} samples at each end.")
    say()
    say("| what is compared | NMSE | reading |")
    say("|---|---|---|")

    e_res = nmse_db(*trim(py["y_resamp"].ravel(), ml["y_resamp"].ravel(), margin))
    say(f"| resampler alone, no channel | **{e_res:.0f} dB** | "
        "`resample` and `resample_poly` build the same filter |")
    check(f"{name}: resampler control", e_res, tol["replay_resampler_nmse_db"],
          e_res <= tol["replay_resampler_nmse_db"], "dB")

    e_static = worst("y_static", "y_static")
    say(f"| replay, tracking removed | **{e_static:.0f} dB** | "
        "`h_hat` layout, spline interpolation, convolution, both resamplings "
        "and the up-conversion agree to the last bit |")
    check(f"{name}: replay with no tracking", e_static,
          tol["replay_static_nmse_db"], e_static <= tol["replay_static_nmse_db"], "dB")

    e_phi = worst("y", "y")
    say(f"| replay, tracking active | **{e_phi:.0f} dB** | "
        "the delay/phase trajectory agrees too |")
    check(f"{name}: replay with tracking", e_phi, tol["replay_nmse_db"],
          e_phi <= tol["replay_nmse_db"], "dB")

    e_off = worst("y", "y_same_start")
    say(f"| ditto, MATLAB given `start` not `start+1` | **{e_off:.0f} dB** | "
        "a deliberate misalignment, to show the comparison resolves one |")
    check(f"{name}: one-sample offset is detected", -e_off,
          -tol["replay_offset_probe_db"], e_off > tol["replay_offset_probe_db"], "dB")
    say()

    for m in range(py["y"].shape[1]):
        a, b = trim(py["y"][:, m], ml["y"][:, m], margin)
        lag = abs(best_lag(a, b))
        check(f"{name}: replay lag, hydrophone {idx[m]}", lag,
              tol["replay_lag_samples"], lag <= tol["replay_lag_samples"], "samples")

    a, b = trim(py["y"][:, 0], ml["y"][:, 0], margin)
    ao, bo = trim(py["y"][:, 0], ml["y_same_start"][:, 0], margin)
    n0 = len(a) // 2
    fig, ax = plt.subplots(2, 1, figsize=(9, 6), constrained_layout=True)
    ax[0].plot(a[n0:n0 + 400], label="Python", lw=1.2)
    ax[0].plot(b[n0:n0 + 400], "--", label="MATLAB", lw=1.2)
    ax[0].set_title(f"{name}: replay, hydrophone {idx[0]}, 400 samples mid-record")
    ax[0].set_xlabel("Sample"); ax[0].legend(); ax[0].grid(alpha=0.3)
    ref = 20 * np.log10(np.abs(a).max())
    ax[1].plot(20 * np.log10(np.abs(ao - bo) + 1e-300) - ref, lw=0.6,
               label="deliberately misaligned")
    ax[1].plot(20 * np.log10(np.abs(a - b) + 1e-300) - ref, lw=0.6, label="aligned")
    ax[1].set_ylim(-320, 0); ax[1].set_xlabel("Sample")
    ax[1].set_title("Python minus MATLAB, dB relative to peak")
    ax[1].legend(); ax[1].grid(alpha=0.3)
    fig.savefig(FIGURES / f"replay_{name}.png", dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Unpack
# ---------------------------------------------------------------------------

def compare_unpack(case, cfg):
    name = case["name"]
    py = loadmat(ARTIFACTS / f"python_{name}_unpack.mat")
    ml = load_v73_complex(ARTIFACTS / f"matlab_{name}_unpack.mat", ("u", "u_fr"))
    tol = cfg["tolerances"]

    say(f"### Unpack ({name})")
    say()
    say(f"Unpacked to {float(py['fs_out'].ravel()[0]):g} Hz on hydrophones "
        f"{case['unpack']['array_index']} (0-based), shape {py['u'].shape} "
        f"(delay, element, time), with and without "
        f"`f_resamp = {float(py['f_resamp'].ravel()[0]):.10g}`.")
    say()
    say("| what is compared | NMSE |")
    say("|---|---|")

    e_plain = nmse_db(py["u"], ml["u"])
    say(f"| unpack, no `f_resamp` | **{e_plain:.0f} dB** |")
    check(f"{name}: unpack without f_resamp", e_plain, tol["unpack_nmse_db"],
          e_plain <= tol["unpack_nmse_db"], "dB")

    e_fr = nmse_db(py["u_fr"], ml["u_fr"])
    say(f"| unpack, `f_resamp` active | **{e_fr:.0f} dB** |")
    check(f"{name}: unpack with f_resamp", e_fr, tol["unpack_nmse_db"],
          e_fr <= tol["unpack_nmse_db"], "dB")

    # Guard against the f_resamp path quietly doing nothing, which would make
    # the row above pass for the wrong reason.
    e_effect = nmse_db(py["u"], py["u_fr"])
    say(f"| how much `f_resamp` changes the output (Python) | "
        f"**{e_effect:.0f} dB** |")
    check(f"{name}: f_resamp actually does something", -e_effect,
          -tol["unpack_f_resamp_effect_db"],
          e_effect > tol["unpack_f_resamp_effect_db"], "dB")
    say()
    say("The third row is the guard on the second: `f_resamp` has to change the "
        "output substantially, or the row above it would agree for the wrong "
        "reason.")
    say()

    t = py["u"].shape[2] // 2
    peak = np.abs(py["u"]).max()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    ax[0].plot(20 * np.log10(np.abs(py["u"][:, 0, t]) / peak + 1e-12), label="Python")
    ax[0].plot(20 * np.log10(np.abs(ml["u"][:, 0, t]) / peak + 1e-12), "--",
               label="MATLAB")
    ax[0].set_ylim(-80, 5)
    ax[0].set_title(f"{name}: unpacked delay profile, mid-record")
    ax[0].set_xlabel("Delay tap"); ax[0].set_ylabel("dB re peak")
    ax[0].legend(); ax[0].grid(alpha=0.3)

    # Sorted magnitude of the disagreement, relative to the peak tap.  The
    # buffer taps are exactly zero in both, hence the floor.
    for key, lab in (("u", "no f_resamp"), ("u_fr", "f_resamp")):
        d = np.sort(np.abs(py[key] - ml[key]).ravel())
        d = d[:: max(1, len(d) // 2000)]
        ax[1].plot(np.linspace(0, 100, len(d)),
                   20 * np.log10(d / peak + 1e-18), label=lab)
    ax[1].set_ylim(-360, 0)
    ax[1].set_title("Sorted |Python - MATLAB|, dB re peak tap")
    ax[1].set_xlabel("Percentile of samples"); ax[1].legend(); ax[1].grid(alpha=0.3)
    fig.savefig(FIGURES / f"unpack_{name}.png", dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Noise
# ---------------------------------------------------------------------------

def theory(beta, Fs, nfreq):
    """Covariance and one-sided PSD implied by the mixing coefficients.

    ``w(n, i) = sum_j sum_k beta(i, j, k) z(n + k, j)`` with ``z`` iid of unit
    variance gives ``C = sum_k B_k B_k^T`` and
    ``S_i(f) = (2 / Fs) sum_j |sum_k beta(i, j, k) e^(-2 pi i f k / Fs)|^2``.
    """
    M, _, K = beta.shape
    C = np.einsum("ijk,ljk->il", beta, beta)
    f = np.linspace(0, Fs / 2, nfreq)
    phase = np.exp(-2j * np.pi * np.outer(f, np.arange(K)) / Fs)
    H = np.einsum("ijk,fk->ijf", beta, phase)
    return C, f, (2.0 / Fs) * np.sum(np.abs(H) ** 2, axis=1)


def coherence(C):
    d = np.sqrt(np.diag(C))
    return np.abs(C / np.outer(d, d))


def compare_noise(case, cfg):
    name = case["name"]
    py = loadmat(ARTIFACTS / f"python_{name}_noise.mat")
    ml = loadmat(ARTIFACTS / f"matlab_{name}_noise.mat")
    tol = cfg["tolerances"]
    nper = cfg["noise"]["welch_nperseg"]

    wp, wm = py["w"], ml["w"]
    beta = np.asarray(py["beta"], dtype=float)
    Fs = float(np.asarray(py["Fs"]).ravel()[0])
    fs = float(cfg["noise"]["fs"])
    alpha = float(np.asarray(py["alpha"]).ravel()[0])
    assert np.allclose(beta, np.asarray(ml["beta"], dtype=float)), \
        "the two runners disagree about beta itself"

    say(f"### Noise ({name})")
    say()
    say(f"{wp.shape[0]} samples on {wp.shape[1]} hydrophones per "
        f"implementation ({wp.shape[0] / fs:.1f} s at {fs / 1e3:g} kHz), "
        f"alpha = {alpha:g}, beta {beta.shape}.  The two draw from different "
        "pseudo-random streams, so they are compared as distributions, and "
        "each is also compared against what `beta` predicts.")
    say()

    sp, sm = wp.std(axis=0), wm.std(axis=0)
    kp, km = kurtosis(wp, axis=0), kurtosis(wm, axis=0)
    rel_std = np.abs(sp - sm) / sp
    rng = np.random.default_rng(0)
    sub = rng.choice(wp.shape[0], size=min(60000, wp.shape[0]), replace=False)
    ks = [ks_2samp(wp[sub, m], wm[sub, m]).statistic for m in range(wp.shape[1])]

    check(f"{name}: noise std, max relative difference", rel_std.max(),
          tol["noise_std_rel"], rel_std.max() <= tol["noise_std_rel"], "")
    if alpha == 2:
        # A Gaussian driver through a linear mixer: excess kurtosis must be
        # zero on both sides.  That is absolute; the pairwise row is a
        # cross-check only.
        check(f"{name}: excess kurtosis vs Gaussian (Python)", np.abs(kp).max(),
              tol["noise_kurtosis_abs"], np.abs(kp).max() <= tol["noise_kurtosis_abs"], "")
        check(f"{name}: excess kurtosis vs Gaussian (MATLAB)", np.abs(km).max(),
              tol["noise_kurtosis_abs"], np.abs(km).max() <= tol["noise_kurtosis_abs"], "")
    check(f"{name}: noise KS statistic, max", max(ks), tol["noise_ks_stat"],
          max(ks) <= tol["noise_ks_stat"], "")

    f, Pp = sg.welch(wp, fs=fs, nperseg=nper, axis=0)
    _, Pm = sg.welch(wm, fs=fs, nperseg=nper, axis=0)
    Pp, Pm = Pp.T, Pm.T
    band = Pp.mean(axis=0) > Pp.max() * 1e-3
    dpsd = 10 * np.log10(Pp[:, band] / Pm[:, band])
    rms = np.sqrt((dpsd ** 2).mean())
    check(f"{name}: noise PSD, rms Python-vs-MATLAB", rms, tol["noise_psd_rms_db"],
          rms <= tol["noise_psd_rms_db"], "dB")

    # beta describes the whole array; the run may have generated a subset of
    # its hydrophones, so the reference has to be taken on the same ones.
    nidx = np.asarray(py["array_index"]).ravel().astype(int)
    C_full, f_th, S_full = theory(beta, Fs, len(f))
    C_th = C_full[np.ix_(nidx, nidx)]
    S_th = S_full[nidx]
    dpy = 10 * np.log10(Pp[:, band] / S_th[:, band])
    dml = 10 * np.log10(Pm[:, band] / S_th[:, band])

    Cp, Cm = np.cov(wp, rowvar=False), np.cov(wm, rowvar=False)
    gp, gm, gth = coherence(Cp), coherence(Cm), coherence(C_th)
    off = ~np.eye(len(gp), dtype=bool)
    rms_xy = np.sqrt(np.mean((gp[off] - gm[off]) ** 2))
    rms_pth = np.sqrt(np.mean((gp[off] - gth[off]) ** 2))
    rms_mth = np.sqrt(np.mean((gm[off] - gth[off]) ** 2))
    gT = coherence(np.einsum("jik,ljk->il", beta, beta)[np.ix_(nidx, nidx)])
    rms_pT = np.sqrt(np.mean((gp[off] - gT[off]) ** 2))
    rms_mT = np.sqrt(np.mean((gm[off] - gT[off]) ** 2))

    say("| measure | Python vs MATLAB | Python vs `beta` | MATLAB vs `beta` |")
    say("|---|---|---|---|")
    say(f"| std, max rel. difference | {rel_std.max():.2%} | - | - |")
    say(f"| excess kurtosis, max abs | {np.abs(kp - km).max():.3f} | "
        f"{np.abs(kp).max():.3f} | {np.abs(km).max():.3f} |")
    say(f"| KS statistic, max | {max(ks):.4f} | - | - |")
    say(f"| PSD, rms over {int(band.sum())} in-band bins | {rms:.3f} dB | "
        f"{np.sqrt((dpy ** 2).mean()):.3f} dB | {np.sqrt((dml ** 2).mean()):.3f} dB |")
    say(f"| coherence, rms off-diagonal | {rms_xy:.4f} | {rms_pth:.4f} | "
        f"{rms_mth:.4f} |")
    say(f"| coherence vs *transposed* `beta` | - | {rms_pT:.4f} | {rms_mT:.4f} |")
    say()
    say("The last row is the control: `beta` is not symmetric in (i, j), so an "
        "implementation mixing `sum_j beta_ji z_j` instead of "
        "`sum_j beta_ij z_j` lands on the transposed covariance.  Both must be "
        "far closer to `beta` than to its transpose.")
    say()

    check(f"{name}: coherence, Python vs MATLAB", rms_xy, tol["noise_coherence_rms"],
          rms_xy <= tol["noise_coherence_rms"], "")
    check(f"{name}: coherence, Python vs theory", rms_pth,
          tol["noise_theory_coherence_rms"], rms_pth <= tol["noise_theory_coherence_rms"], "")
    check(f"{name}: coherence, MATLAB vs theory", rms_mth,
          tol["noise_theory_coherence_rms"], rms_mth <= tol["noise_theory_coherence_rms"], "")
    check(f"{name}: closer to beta than to its transpose (Python)", rms_pth,
          rms_pT, rms_pth < rms_pT, "")
    check(f"{name}: closer to beta than to its transpose (MATLAB)", rms_mth,
          rms_mT, rms_mth < rms_mT, "")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
    ax[0].semilogy(f / 1e3, Pp[0], label="Python")
    ax[0].semilogy(f / 1e3, Pm[0], "--", label="MATLAB")
    ax[0].semilogy(f_th / 1e3, S_th[0], ":", color="k", label="from beta")
    ax[0].set_title(f"{name}: noise PSD, hydrophone 0")
    ax[0].set_xlabel("Frequency [kHz]"); ax[0].legend(); ax[0].grid(alpha=0.3)
    edges = np.linspace(-4 * sp[0], 4 * sp[0], 160)
    ax[1].hist(wp[:, 0], bins=edges, density=True, histtype="step", label="Python")
    ax[1].hist(wm[:, 0], bins=edges, density=True, histtype="step", label="MATLAB")
    ax[1].set_yscale("log"); ax[1].set_title("Amplitude distribution")
    ax[1].legend(); ax[1].grid(alpha=0.3)
    im = ax[2].imshow(gp - gth, cmap="coolwarm")
    ax[2].set_title("Python coherence minus theory"); fig.colorbar(im, ax=ax[2])
    fig.savefig(FIGURES / f"noise_{name}.png", dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------------------

def main():
    cfg = config()
    FIGURES.mkdir(parents=True, exist_ok=True)
    say("# Cross-implementation calibration report")
    say()
    say(f"Files from Zenodo record [{cfg['data']['zenodo_record']}]"
        f"(https://doi.org/{cfg['data']['zenodo_doi']}).  Both runners were "
        "handed the same probe samples from `probe.mat`, so nothing in the "
        "replay or unpack comparisons depends on a random number generator.")
    say()

    for case in cases():
        say(f"## {case['name']} (`{case['channel_file']}`, "
            f"`{case['noise_file']}`)")
        say()
        compare_replay(case, cfg)
        compare_unpack(case, cfg)
        compare_noise(case, cfg)

    say("## Verdict")
    say()
    say("| check | value | tolerance | |")
    say("|---|---|---|---|")
    failed = 0
    for nm, value, limit, ok, units in CHECKS:
        failed += not ok
        say(f"| {nm} | {value:.4g} {units} | {limit:.4g} {units} | "
            f"{'pass' if ok else '**FAIL**'} |")
    say()
    say(f"{len(CHECKS) - failed} of {len(CHECKS)} checks pass.")

    (ARTIFACTS / "report.md").write_text("\n".join(LINES) + "\n")
    print(f"\nwrote {ARTIFACTS / 'report.md'} and {FIGURES}/*.png")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

# [EOF]
