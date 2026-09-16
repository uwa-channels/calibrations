"""Build the probe signal both implementations replay.

The one thing a cross-language comparison must not do is draw its own random
numbers on each side: NumPy's Mersenne Twister and MATLAB's are seeded the same
way and still produce different streams, so "same seed" is not "same signal".
The probe is therefore generated once, here, and written to a MAT-file that
both runners read.  Nothing downstream calls a random number generator.

Writes ``artifacts/probe.mat`` with

    input    (N, 1) real passband probe at ``fs``
    symbols  (n_symbols, 1) the +/-1 sequence, for matched filtering
    pulse    (span*sps+1, 1) the RRC pulse, likewise

Author: Zhengnan Li
Email : uwa-channels@ofdm.link
License: MIT
"""

import numpy as np
from scipy.io import savemat

from calib import ARTIFACTS, config


def rrc(rolloff, sps, span):
    """Root-raised-cosine pulse, unit energy, ``span*sps + 1`` taps.

    Written out rather than pulled from a toolbox so that the probe does not
    depend on which signal-processing package happens to be installed.
    """
    n = np.arange(-span * sps / 2, span * sps / 2 + 1, dtype=float)
    t = n / sps
    h = np.zeros_like(t)

    # The three cases are the removable singularities of the closed form.
    at_zero = np.isclose(t, 0.0)
    at_pole = np.isclose(np.abs(4 * rolloff * t), 1.0)
    ok = ~(at_zero | at_pole)

    h[at_zero] = 1.0 - rolloff + 4 * rolloff / np.pi
    if rolloff > 0:
        h[at_pole] = (
            rolloff
            / np.sqrt(2)
            * (
                (1 + 2 / np.pi) * np.sin(np.pi / (4 * rolloff))
                + (1 - 2 / np.pi) * np.cos(np.pi / (4 * rolloff))
            )
        )
    tt = t[ok]
    h[ok] = (
        np.sin(np.pi * tt * (1 - rolloff))
        + 4 * rolloff * tt * np.cos(np.pi * tt * (1 + rolloff))
    ) / (np.pi * tt * (1 - (4 * rolloff * tt) ** 2))

    return h / np.sqrt(np.sum(h**2))


def main():
    cfg = config()["probe"]
    fs, fc = cfg["fs"], cfg["fc"]
    sps = fs / cfg["symbol_rate"]
    if sps != int(sps):
        raise SystemExit(f"fs/symbol_rate = {sps} must be an integer.")
    sps = int(sps)

    rng = np.random.default_rng(cfg["symbol_seed"])
    symbols = rng.choice([-1.0, 1.0], size=cfg["n_symbols"])

    upsampled = np.zeros(len(symbols) * sps)
    upsampled[::sps] = symbols
    pulse = rrc(cfg["rolloff"], sps, cfg["span_symbols"])
    baseband = np.convolve(upsampled, pulse, mode="full")

    carrier = np.exp(2j * np.pi * fc * np.arange(len(baseband)) / fs)
    passband = np.real(baseband * carrier)

    guard = np.zeros(cfg["guard_samples"])
    probe = np.concatenate((guard, passband, guard))
    probe = probe / np.max(np.abs(probe))

    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    savemat(
        ARTIFACTS / "probe.mat",
        {
            "input": probe.reshape(-1, 1),
            "symbols": symbols.reshape(-1, 1),
            "pulse": pulse.reshape(-1, 1),
            "sps": float(sps),
        },
        do_compression=True,
    )
    print(
        f"probe.mat: {len(probe)} samples ({len(probe) / fs:.3f} s at {fs / 1e3:g} kHz), "
        f"{len(symbols)} symbols at {cfg['symbol_rate'] / 1e3:g} kBd, "
        f"rolloff {cfg['rolloff']}, one-sided bandwidth "
        f"{cfg['symbol_rate'] * (1 + cfg['rolloff']) / 2 / 1e3:g} kHz"
    )


if __name__ == "__main__":
    main()

# [EOF]
