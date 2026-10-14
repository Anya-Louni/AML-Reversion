"""Shared helpers: config, paths, logging, downloads."""
from __future__ import annotations

import logging
import os
import time
from pathlib import Path

import requests
import yaml

ROOT = Path(__file__).resolve().parents[1]

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(name)


def load_config(path: Path | None = None) -> dict:
    with open(path or ROOT / "config.yaml") as f:
        return yaml.safe_load(f)


def path(cfg: dict, key: str, *parts: str) -> Path:
    p = ROOT / cfg["paths"][key]
    p = p.joinpath(*parts)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def download(url: str, dest: Path, retries: int = 4, timeout: int = 120) -> Path:
    """Download url to dest unless dest exists and matches the remote size."""
    log = get_logger("download")
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):
        try:
            with requests.get(url, stream=True, timeout=timeout, allow_redirects=True) as r:
                r.raise_for_status()
                size = int(r.headers.get("Content-Length", -1))
                if dest.exists() and size > 0 and dest.stat().st_size == size:
                    return dest
                tmp = dest.with_suffix(dest.suffix + ".part")
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
                if size > 0 and tmp.stat().st_size != size:
                    raise IOError(f"size mismatch {tmp.stat().st_size} != {size}")
                os.replace(tmp, dest)
                log.info("downloaded %s (%.1f MB)", dest.name, dest.stat().st_size / 1e6)
                return dest
        except Exception as e:  # network errors: retry with backoff
            log.warning("attempt %d failed for %s: %s", attempt + 1, url, e)
            time.sleep(2 ** attempt)
    raise RuntimeError(f"could not download {url}")
