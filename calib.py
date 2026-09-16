"""Shared helpers: config, paths, and the Zenodo fetch.

``make_probe.py`` and ``run_python.py`` import this; ``run_matlab.m`` reads
``config.json`` directly with ``jsondecode``, so neither language owns a second
copy of the parameters.

Author: Zhengnan Li
Email : uwa-channels@ofdm.link
License: MIT
"""

import hashlib
import json
import os
import shutil
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ARTIFACTS = ROOT / "artifacts"
FIGURES = ROOT / "figures"

_TIMEOUT = 300.0
_CHUNK = 1 << 20


def config():
    """The parsed ``config.json``."""
    with open(ROOT / "config.json") as f:
        return json.load(f)


def cases():
    """The channel/noise pairs to calibrate on.

    ``UWA_CALIBRATION_CASES``, a comma-separated list of names, restricts the
    run to a subset -- useful when iterating locally, or in CI where the full
    set is three quarters of a gigabyte to fetch the first time.
    """
    all_cases = config()["data"]["cases"]
    want = os.environ.get("UWA_CALIBRATION_CASES", "").strip()
    if not want:
        return all_cases
    names = [n.strip() for n in want.split(",") if n.strip()]
    chosen = [c for c in all_cases if c["name"] in names]
    missing = set(names) - {c["name"] for c in chosen}
    if missing:
        raise SystemExit(f"UWA_CALIBRATION_CASES names no such case: "
                         f"{', '.join(sorted(missing))}")
    return chosen


def data_dir():
    """Where the library MAT-files live: ``UWA_CHANNELS_CACHE``, else here."""
    env = os.environ.get("UWA_CHANNELS_CACHE")
    return Path(env) if env else ROOT


def fetch(name, md5, record):
    """Return a local path to ``name``, downloading it from Zenodo once.

    The digest is Zenodo's own, so this is a transfer check against their
    manifest rather than a security one.
    """
    dest = data_dir() / name
    if dest.is_file() and _md5(dest) == md5:
        return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    url = f"https://zenodo.org/api/records/{record}/files/{name}/content"
    print(f"downloading {name} from {url}")
    try:
        # nosec B310: https literal above with a file name appended.
        with urllib.request.urlopen(url, timeout=_TIMEOUT) as r, open(part, "wb") as h:
            shutil.copyfileobj(r, h, _CHUNK)
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        part.unlink(missing_ok=True)
        raise SystemExit(f"cannot download {name} from Zenodo: {exc}")

    got = _md5(part)
    if got != md5:
        part.unlink(missing_ok=True)
        raise SystemExit(
            f"{name} downloaded with MD5 {got}, expected {md5}; the record may "
            f"have been revised, or the transfer was truncated."
        )
    part.replace(dest)
    return dest


def ensure_case(case):
    """Make sure one case's files are present; return (channel, noise) paths."""
    record = config()["data"]["zenodo_record"]
    return (
        fetch(case["channel_file"], case["md5"][case["channel_file"]], record),
        fetch(case["noise_file"], case["md5"][case["noise_file"]], record),
    )


def _md5(path):
    digest = hashlib.md5(usedforsecurity=False)
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


if __name__ == "__main__":
    for case in cases():
        for p in ensure_case(case):
            print("ok", p)

# [EOF]
