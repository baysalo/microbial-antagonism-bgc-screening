from pathlib import Path
import pytest

from mine_antagonism import (
    DEFAULT_HYDROPHOBIC_AAS,
    compute_atgc_gc_percent,
    load_gtf_cds_lookup,
    parse_annotation_attributes,
    predict_sec_signal_heuristic,
)


def test_empty_gtf_has_no_strain_column(tmp_path: Path):
    empty_gtf = tmp_path / "empty.gtf"
    empty_gtf.write_text("# header only\n", encoding="utf-8")
    lookup, count = load_gtf_cds_lookup(empty_gtf, "Test_strain")
    assert count == 0
    assert "strain" not in lookup.columns
    assert list(lookup.columns) == ["match_key", "gtf_product"]


def test_proline_excluded_from_hydrophobic_set():
    assert "P" not in DEFAULT_HYDROPHOBIC_AAS
    # Construct sequence with high Proline in h-region (aa 6..25) and low A/V/I/L/M/F/W
    pro_rich = "MKR" + "KK" + ("P" * 20) + "AA" + ("S" * 15)
    res = predict_sec_signal_heuristic(pro_rich)
    assert res["has_sec_sig"] is False


def test_unified_atgc_gc_calculation():
    gc_pct, n_cnt, atgc_cnt = compute_atgc_gc_percent("GGCCNNAATT")
    assert atgc_cnt == 8
    assert n_cnt == 2
    assert gc_pct == pytest.approx(25.0)


def test_gff3_and_repeated_keys():
    attrs, warns = parse_annotation_attributes('ID=cds1;Note=first;Note=second;product="chitinase"')
    assert attrs["ID"] == ["cds1"]
    assert attrs["Note"] == ["first", "second"]
    assert len(warns) == 1