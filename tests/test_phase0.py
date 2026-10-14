"""Phase 0 unit tests.

All tables in this file are MOCK DATA built by hand for testing only.
They are never used by the real pipeline.
"""
import gzip

import numpy as np
import pandas as pd
import pytest

from pipeline.phase0_data import drug_names, vangalen
from pipeline.phase0_data.genes import HGNCMapper
from pipeline.phase0_data.lincs import select_signatures


@pytest.fixture
def hgnc():
    # MOCK HGNC table
    df = pd.DataFrame({
        "symbol": ["CD34", "PROM1", "CD14", "LYZ", "GENEX"],
        "status": ["Approved"] * 5,
        "prev_symbol": [None, "OLDPROM", None, None, "SHARED"],
        "alias_symbol": [None, None, "MONO1|SHARED2", "SHARED2", None],
        "ensembl_gene_id": ["ENSG1", "ENSG2", "ENSG3", "ENSG4", None],
    })
    return HGNCMapper.from_table(df)


def test_hgnc_symbol_mapping(hgnc):
    assert hgnc.map_symbol("CD34") == ("CD34", "approved")
    assert hgnc.map_symbol("OLDPROM") == ("PROM1", "previous")
    assert hgnc.map_symbol("MONO1") == ("CD14", "alias")
    assert hgnc.map_symbol("SHARED2") == (None, "unmapped")  # ambiguous alias is not used
    assert hgnc.map_symbol("NOPE") == (None, "unmapped")
    m = hgnc.map_ensembl(["ENSG3.5", "ENSG9"])
    assert m.iloc[0] == "CD14" and pd.isna(m.iloc[1])


def test_drug_name_matching():
    lincs = ["axitinib", "sirolimus", "BRD-K1", "JQ-1-R", "BRD-K2", "baicalin", "venetoclax-x"]
    aliases = {"JQ1-(+)": "BRD-K2", "JQ1-(-)": "JQ-1-R"}
    beat = ["Axitinib (AG-013736)", "Rapamycin", "JQ1", "Baiclein", "Unknown-1"]
    m = drug_names.match_drugs(beat, lincs, aliases).set_index("beat_name")
    assert m.loc["Axitinib (AG-013736)", "lincs_names"] == "axitinib"
    assert m.loc["Rapamycin", "match_how"] == "curated"
    assert m.loc["Baiclein", "match_how"] == "excluded"
    assert m.loc["Unknown-1", "match_how"] == "none"


def test_jq1_inactive_enantiomer_removed():
    # MOCK: name normalization would pick up JQ-1-R ('jq1r' != 'jq1', so via alias only)
    m = drug_names.match_drugs(["JQ1"], ["JQ-1-R", "BRD-K2"], {"JQ1": "JQ-1-R", "jq-1": "BRD-K2"})
    assert m.lincs_names.iat[0] == "BRD-K2"


def test_name_keys():
    assert drug_names.name_keys("Quizartinib (AC220)") == ["quizartinibac220", "quizartinib", "ac220"]


def test_sample_info():
    assert vangalen.sample_info("AML328-D29") == {"sample": "AML328-D29", "patient": "AML328", "day": 29, "sample_type": "AML"}
    assert vangalen.sample_info("BM5-34p")["sample_type"] == "healthy_BM"
    assert vangalen.sample_info("OCI-AML3")["sample_type"] == "cell_line"


def _write_gz(p, df):
    with gzip.open(p, "wt") as f:
        df.to_csv(f, sep="\t")


def test_vangalen_build_mock(tmp_path, hgnc):
    # MOCK count and annotation tables in the GEO layout
    counts = pd.DataFrame({"BM1_A": [1, 0, 3], "BM1_B": [0, 2, 1], "BM1_C": [5, 5, 5]},
                          index=pd.Index(["CD34", "OLDPROM", "PROM1"], name="Gene"))
    anno = pd.DataFrame({"PredictionRefined": ["normal", "normal"], "CellType": ["HSC", "Mono"]},
                        index=pd.Index(["BM1_A", "BM1_B"], name="Cell"))
    _write_gz(tmp_path / "GSM1_BM1.dem.txt.gz", counts)
    _write_gz(tmp_path / "GSM1_BM1.anno.txt.gz", anno)
    a, s = vangalen.build(tmp_path, hgnc)
    assert a.n_obs == 2  # unannotated cell dropped
    assert list(a.var_names) == ["CD34", "PROM1"]  # OLDPROM collapsed into PROM1
    prom1 = a[:, "PROM1"].X.toarray().ravel()
    np.testing.assert_array_equal(prom1, [3, 3])  # 0+3 and 2+1
    assert s["gene_mapping"]["n_collapsed"] == 1


def test_vangalen_missing_cells_raise(tmp_path):
    counts = pd.DataFrame({"X_A": [1]}, index=pd.Index(["CD34"], name="Gene"))
    anno = pd.DataFrame({"CellType": ["HSC", "HSC"]}, index=pd.Index(["X_A", "X_Z"], name="Cell"))
    _write_gz(tmp_path / "a.dem.txt.gz", counts)
    _write_gz(tmp_path / "a.anno.txt.gz", anno)
    with pytest.raises(ValueError):
        vangalen.read_sample(tmp_path / "a.dem.txt.gz", tmp_path / "a.anno.txt.gz")


def test_select_signatures_cap_and_order():
    # MOCK siginfo
    sig = pd.DataFrame({
        "sig_id": [f"s{i}" for i in range(6)],
        "cmap_name": ["d1"] * 4 + ["d2", "d1"],
        "cell_iname": ["THP1"] * 4 + ["THP1", "MCF7"],
        "is_hiq": [0, 1, 1, 0, 1, 1],
        "is_exemplar_sig": [0, 0, 1, 0, 0, 0],
        "tas": [0.9, 0.2, 0.1, 0.5, 0.3, 0.3],
    })
    sel = select_signatures(sig, {"d1": "Drug1"}, ["THP1"], cap=2)
    assert list(sel.sig_id) == ["s2", "s1"]  # hiq first, then exemplar, then TAS
    assert set(sel.drug) == {"Drug1"}
