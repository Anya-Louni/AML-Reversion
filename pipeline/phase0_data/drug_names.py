"""Match Beat AML inhibitor names to LINCS compound names.

Matching order for each Beat AML name:
  1. a curated synonym table (below), checked by hand against LINCS names
  2. normalized full name, the name without its parenthetical, and each
     parenthetical alias (e.g. 'Axitinib (AG-013736)' -> axitinib, ag013736)
  3. the same keys against LINCS compound aliases
Then LINCS entries listed in LINCS_EXCLUDE are removed (e.g. inactive enantiomers).

Normalization lowercases and strips everything except letters and digits.
"""
from __future__ import annotations

import re

import pandas as pd

# Curated synonyms: Beat AML name -> (LINCS cmap_name, match quality, note).
# "exact" = same active compound under another name or salt form.
# "approximate" = closely related but not identical (prodrug, racemate, ...).
CURATED = {
    "Azacytidine": ("azacitidine", "exact", "spelling variant"),
    "BEZ235": ("NVP-BEZ235", "exact", "same compound (dactolisib)"),
    "Flavopiridol": ("alvocidib", "exact", "INN of flavopiridol"),
    "INK-128": ("MLN-0128", "exact", "INK-128 = MLN0128 = sapanisertib"),
    "Rapamycin": ("sirolimus", "exact", "INN of rapamycin"),
    "Vargetef": ("nintedanib", "exact", "Vargatef is the trade name of nintedanib"),
    "Perhexiline maleate": ("perhexiline", "exact", "salt form"),
    "Venetoclax": ("GDC-0199", "exact", "GDC-0199 = ABT-199 = venetoclax"),
    "AT-101": ("gossypol", "approximate", "AT-101 is R-(-)-gossypol; LINCS entry may be racemic"),
    "R406": ("fostamatinib", "approximate", "fostamatinib is the prodrug of R406"),
    "Nutlin 3a": ("nutlin-3", "approximate", "LINCS nutlin-3 pools two pert_ids (one aliased NUTLIN-3A); enantiomer not confirmed"),
}
# LINCS entries that a name match picks up but that are a different form of the compound.
LINCS_EXCLUDE = {
    "JQ1": {"JQ-1-R": "JQ-1-R is (-)-JQ1, the inactive enantiomer; Beat AML JQ1 is the active (+) form"},
}
# Pairs that look similar but are different compounds; never match them.
EXCLUDED = {
    "Baiclein": "LINCS has baicalin (glycoside), not baicalein",
    "Bay 11-7085": "LINCS has BAY-11-7082, a different compound",
    "AMPK Inhibitor": "likely Compound C (dorsomorphin) but not confirmed in Beat AML metadata",
}


def norm(x) -> str:
    return re.sub(r"[^a-z0-9]", "", str(x).lower())


def name_keys(name: str) -> list[str]:
    parts = [name, re.sub(r"\s*\(.*?\)", "", name)] + re.findall(r"\(([^)]+)\)", name)
    return list(dict.fromkeys(norm(p) for p in parts if norm(p)))


def match_drugs(beat_names, lincs_names, lincs_aliases: dict[str, str] | None = None) -> pd.DataFrame:
    """Return one row per Beat AML name with the matched LINCS cmap_name(s)."""
    lincs_names = set(lincs_names)
    by_key: dict[str, set[str]] = {}
    for n in lincs_names:
        by_key.setdefault(norm(n), set()).add(n)
    alias_key: dict[str, set[str]] = {}
    for alias, n in (lincs_aliases or {}).items():
        if n in lincs_names:
            alias_key.setdefault(norm(alias), set()).add(n)

    rows = []
    for b in sorted(set(beat_names)):
        hits, how, quality, note = set(), None, None, ""
        if b in EXCLUDED:
            how, note = "excluded", EXCLUDED[b]
        else:
            if b in CURATED and CURATED[b][0] in lincs_names:
                hits, how = {CURATED[b][0]}, "curated"
                quality, note = CURATED[b][1], CURATED[b][2]
            if not hits:
                for k in name_keys(b):
                    hits |= by_key.get(k, set())
                how, quality = ("name", "exact") if hits else (None, None)
            if not hits:
                for k in name_keys(b):
                    hits |= alias_key.get(k, set())
                how, quality = ("alias", "exact") if hits else (None, None)
            for bad, why in LINCS_EXCLUDE.get(b, {}).items():
                if bad in hits:
                    hits.discard(bad)
                    note = (note + "; " if note else "") + why
        rows.append({"beat_name": b, "lincs_names": "|".join(sorted(hits)) or None,
                     "match_how": how or "none", "match_quality": quality, "note": note})
    return pd.DataFrame(rows)
