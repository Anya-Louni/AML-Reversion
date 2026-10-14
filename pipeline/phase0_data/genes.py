"""Gene identifier harmonization to current HGNC approved symbols."""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd


@dataclass
class HGNCMapper:
    approved: set[str]
    prev: dict[str, str]          # unambiguous previous symbol -> approved
    alias: dict[str, str]         # unambiguous alias symbol -> approved
    ensembl: dict[str, str]       # ENSG -> approved
    stats: dict = field(default_factory=dict)

    @classmethod
    def from_table(cls, df: pd.DataFrame) -> "HGNCMapper":
        df = df[df["status"] == "Approved"]
        approved = set(df["symbol"])

        def unambiguous(col: str) -> dict[str, str]:
            pairs = (
                df[["symbol", col]].dropna()
                .assign(key=lambda d: d[col].str.split("|"))
                .explode("key")[["key", "symbol"]]
            )
            pairs = pairs[~pairs["key"].isin(approved)]  # approved symbols win
            counts = pairs.groupby("key")["symbol"].nunique()
            keep = counts[counts == 1].index
            return dict(pairs[pairs["key"].isin(keep)].drop_duplicates("key").values)

        ens = df[["ensembl_gene_id", "symbol"]].dropna()
        return cls(approved, unambiguous("prev_symbol"), unambiguous("alias_symbol"),
                   dict(ens.values))

    @classmethod
    def from_file(cls, path) -> "HGNCMapper":
        cols = ["symbol", "status", "prev_symbol", "alias_symbol", "ensembl_gene_id"]
        return cls.from_table(pd.read_csv(path, sep="\t", usecols=cols, dtype=str, low_memory=False))

    def map_symbol(self, s: str) -> tuple[str | None, str]:
        """Return (approved symbol or None, how it was mapped)."""
        if s in self.approved:
            return s, "approved"
        if s in self.prev:
            return self.prev[s], "previous"
        if s in self.alias:
            return self.alias[s], "alias"
        return None, "unmapped"

    def map_symbols(self, symbols) -> pd.DataFrame:
        rows = [(s, *self.map_symbol(s)) for s in symbols]
        out = pd.DataFrame(rows, columns=["original", "hgnc", "how"])
        return out

    def map_ensembl(self, ids) -> pd.Series:
        return pd.Series([self.ensembl.get(str(i).split(".")[0]) for i in ids], index=list(ids))


def summarize_mapping(m: pd.DataFrame) -> dict:
    vc = m["how"].value_counts().to_dict()
    mapped = m.dropna(subset=["hgnc"])
    return {
        "n_input": len(m),
        **{f"n_{k}": int(v) for k, v in vc.items()},
        "n_unique_hgnc": int(mapped["hgnc"].nunique()),
        "n_collapsed": int(len(mapped) - mapped["hgnc"].nunique()),
    }
