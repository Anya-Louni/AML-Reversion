"""Fetch raw inputs into data/raw/ (skips files that are already complete)."""
from __future__ import annotations

from pipeline.common import download, get_logger, path

log = get_logger("phase0.download")


def fetch_all(cfg: dict) -> None:
    src = cfg["sources"]
    download(src["vangalen"]["geo_tar"], path(cfg, "raw", "vangalen", "GSE116256_RAW.tar"))
    b = src["beataml"]
    for f in b["files"]:
        download(f"{b['base_url']}/{f}", path(cfg, "raw", "beataml", f))
    lc = src["lincs"]
    for f in lc["lincs2020_files"]:
        download(f"{lc['lincs2020_base']}/{f}", path(cfg, "raw", "lincs", "lincs2020", f))
    for url in lc["geo_meta"]:
        download(url, path(cfg, "raw", "lincs", url.rsplit("/", 1)[1]))
    download(src["hgnc"]["url"], path(cfg, "raw", "hgnc", "hgnc_complete_set.txt"))
