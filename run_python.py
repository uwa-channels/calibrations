"""Run the Python implementation and write its results to artifacts/.

Produces
    artifacts/python_replay.mat   y            (N, M) replay output
                                  y_resamp     (N2, 1) resampler-only control
    artifacts/python_noise.mat    w            (n, M) noisegen output
                                  beta         (M, M, K) in MATLAB layout

The resampler control matters: ``replay`` puts the probe through a rational
resampler twice, and MATLAB's ``resample`` and SciPy's ``resample_poly`` build
their anti-aliasing filters differently (least-squares versus windowed sinc,
same Kaiser beta and same length).  That difference alone puts a floor under
any cross-language comparison, so measure it separately instead of attributing
it to ``replay``.

Author: Zhengnan Li
Email : uwa-channels@ofdm.link
License: MIT
"""

import subprocess
import sys
from fractions import Fraction

import numpy as np
import scipy.signal as sg
from scipy.io import loadmat, savemat

from calib import ARTIFACTS, config, ensure_data


def git_describe(path):
    try:
        out = subprocess.run(
            ["git", "-C", str(path), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=10,
        )
        return out.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def main():
    import uwa_channels
    from uwa_channels import load_channel, load_noise, noisegen, replay

    cfg = config()
    rcfg, ncfg = cfg["replay"], cfg["noise"]
    channel_path, noise_path = ensure_data()

    probe = loadmat(ARTIFACTS / "probe.mat")
    x = np.asarray(probe["input"]).ravel()

    channel = load_channel(str(channel_path))
    fs_delay = channel["params"]["fs_delay"][0, 0]

    # replay's own resampling ratio, reproduced here for the control.
    frac = Fraction(fs_delay / rcfg["fs"]).limit_denominator()
    down = sg.resample_poly(x, frac.numerator, frac.denominator)
    y_resamp = sg.resample_poly(down, frac.denominator, frac.numerator)

    array_index = np.asarray(rcfg["array_index"], dtype=int)
    y = replay(x, rcfg["fs"], array_index, channel, start=rcfg["start"])

    # The same replay with the phase trajectory removed.  h_hat, the spline
    # interpolation and the two resamplings are then the only things acting, so
    # this separates a disagreement about the signal path from a disagreement
    # about how phi_hat is indexed against it.
    static = {k: channel[k] for k in ("h_hat", "params", "version")}
    y_static = replay(x, rcfg["fs"], array_index, static, start=rcfg["start"])

    savemat(
        ARTIFACTS / "python_replay.mat",
        {
            "y": y,
            "y_static": y_static,
            "y_resamp": y_resamp.reshape(-1, 1),
            "start": float(rcfg["start"]),
            "array_index": array_index.reshape(1, -1).astype(float),
            "fs_delay": float(fs_delay),
            "resamp_p": float(frac.numerator),
            "resamp_q": float(frac.denominator),
            "impl": "python",
            "version": uwa_channels.__version__ if hasattr(uwa_channels, "__version__") else "",
            "commit": git_describe(__import__("pathlib").Path(uwa_channels.__file__).parent),
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
        },
        do_compression=True,
    )
    print(f"python replay: y {y.shape}, control {y_resamp.shape}, p/q = {frac}")

    noise = load_noise(str(noise_path))
    beta = np.asarray(noise["beta"])
    # Store beta in the stored MATLAB layout so the comparer has one reference
    # regardless of which loader read it.  h5py reverses every axis.
    if beta.ndim == 3 and beta.shape[0] != beta.shape[1]:
        beta = np.transpose(beta, (2, 1, 0))

    np.random.seed(ncfg["seed"])
    idx = list(ncfg["array_index"])
    w = noisegen((ncfg["n_samples"], len(idx)), ncfg["fs"], idx, noise)

    savemat(
        ARTIFACTS / "python_noise.mat",
        {
            "w": w,
            "beta": beta,
            "Fs": float(np.asarray(noise["Fs"]).ravel()[0]),
            "fs": float(ncfg["fs"]),
            "alpha": float(np.asarray(noise["alpha"]).ravel()[0]),
            "array_index": np.asarray(idx, dtype=float).reshape(1, -1),
            "impl": "python",
        },
        do_compression=True,
    )
    print(f"python noisegen: w {w.shape}, beta {beta.shape}")


if __name__ == "__main__":
    main()

# [EOF]
