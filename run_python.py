"""Run the Python implementation and write its results to artifacts/.

For each case in config.json, produces
    artifacts/python_<case>_replay.mat   y, y_static, y_resamp
    artifacts/python_<case>_unpack.mat   u, u_fr
    artifacts/python_<case>_noise.mat    w, beta

The resampler control matters: ``replay`` puts the probe through a rational
resampler twice, so if MATLAB's ``resample`` and SciPy's ``resample_poly``
disagreed at all, that would floor every other comparison.  Measuring it
separately is what lets the rest be read as statements about the algorithms.

``u`` and ``u_fr`` are the same channel unpacked without and with ``f_resamp``.
The released files carry no ``f_resamp``, and the path needs exercising: it is
where a second off-by-one lived until Sep. 16, 2026.

Author: Zhengnan Li
Email : uwa-channels@ofdm.link
License: MIT
"""

import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy.signal as sg
from scipy.io import loadmat, savemat

from calib import ARTIFACTS, cases, config, ensure_case


def run_case(case, cfg, x):
    from uwa_channels import load_channel, load_noise, noisegen, replay, unpack

    name = case["name"]
    fs = cfg["probe"]["fs"]
    channel_path, noise_path = ensure_case(case)
    channel = load_channel(str(channel_path))
    fs_delay = channel["params"]["fs_delay"][0, 0]

    # ---- replay -----------------------------------------------------------
    frac = Fraction(fs_delay / fs).limit_denominator()
    down = sg.resample_poly(x, frac.numerator, frac.denominator)
    y_resamp = sg.resample_poly(down, frac.denominator, frac.numerator)

    idx = np.asarray(case["replay"]["array_index"], dtype=int)
    start = case["replay"]["start"]
    y = replay(x, fs, idx, channel, start=start)

    # The same replay with the phase trajectory removed, so that h_hat, the
    # spline interpolation and the two resamplings are the only things acting.
    static = {k: channel[k] for k in ("h_hat", "params", "version")}
    y_static = replay(x, fs, idx, static, start=start)

    savemat(ARTIFACTS / f"python_{name}_replay.mat", {
        "y": y, "y_static": y_static, "y_resamp": y_resamp.reshape(-1, 1),
        "start": float(start), "array_index": idx.reshape(1, -1).astype(float),
        "fs_delay": float(fs_delay),
        "resamp_p": float(frac.numerator), "resamp_q": float(frac.denominator),
    }, do_compression=True)
    print(f"[{name}] replay: y {y.shape}, control {y_resamp.shape}, p/q = {frac}")

    # ---- unpack -----------------------------------------------------------
    ucfg, shared = case["unpack"], cfg["unpack"]
    uidx = list(ucfg["array_index"])
    u = unpack(ucfg["fs_out"], uidx, channel, shared["buffer_left"],
               shared["buffer_right"])
    # Build the f_resamp channel explicitly rather than copying the open
    # file: an h5py.File also carries `meta` and `#refs#`, which unpack has no
    # use for and which would differ from what MATLAB's load() hands over.
    tracked = {k: channel[k] for k in ("h_hat", "params", "version")}
    for k in ("phi_hat", "theta_hat"):
        if k in channel:
            tracked[k] = channel[k]
    tracked["f_resamp"] = np.array([[ucfg["f_resamp"]]])
    u_fr = unpack(ucfg["fs_out"], uidx, tracked,
                  shared["buffer_left"], shared["buffer_right"])
    savemat(ARTIFACTS / f"python_{name}_unpack.mat", {
        "u": u, "u_fr": u_fr, "fs_out": float(ucfg["fs_out"]),
        "f_resamp": float(ucfg["f_resamp"]),
        "array_index": np.asarray(uidx, dtype=float).reshape(1, -1),
    }, do_compression=True)
    print(f"[{name}] unpack: {u.shape} at {ucfg['fs_out']} Hz, "
          f"f_resamp {ucfg['f_resamp']}")

    # ---- noise ------------------------------------------------------------
    ncfg = cfg["noise"]
    noise = load_noise(str(noise_path))
    beta = np.asarray(noise["beta"])
    # Store beta in the stored MATLAB layout so the comparer has one reference
    # whichever loader read it.  h5py reverses every axis.
    if beta.ndim == 3 and beta.shape[0] != beta.shape[1]:
        beta = np.transpose(beta, (2, 1, 0))

    np.random.seed(ncfg["seed"])
    nidx = list(case["noise"]["array_index"])
    w = noisegen((ncfg["n_samples"], len(nidx)), ncfg["fs"], nidx, noise)
    savemat(ARTIFACTS / f"python_{name}_noise.mat", {
        "w": w, "beta": beta,
        "Fs": float(np.asarray(noise["Fs"]).ravel()[0]),
        "alpha": float(np.asarray(noise["alpha"]).ravel()[0]),
        "fs": float(ncfg["fs"]),
        "array_index": np.asarray(nidx, dtype=float).reshape(1, -1),
    }, do_compression=True)
    print(f"[{name}] noisegen: w {w.shape}, beta {beta.shape}")


def main():
    cfg = config()
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    x = np.asarray(loadmat(ARTIFACTS / "probe.mat")["input"]).ravel()

    import uwa_channels
    savemat(ARTIFACTS / "python_env.mat", {
        "impl": "python", "python": sys.version.split()[0],
        "numpy": np.__version__, "scipy": __import__("scipy").__version__,
        "package": str(Path(uwa_channels.__file__).parent),
    })

    for case in cases():
        run_case(case, cfg, x)


if __name__ == "__main__":
    main()

# [EOF]
