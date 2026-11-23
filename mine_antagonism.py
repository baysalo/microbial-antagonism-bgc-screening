import argparse
import hashlib
import json
import logging
import platform
import re
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import yaml
from Bio import SeqIO
from Bio.Seq import Seq

LOGGER = logging.getLogger("antagonism_miner")

DEFAULT_HYDROPHOBIC_AAS = frozenset({"A", "V", "I", "L", "M", "F", "W"})

DEFAULT_WEIGHTS: Dict[str, float] = {
    "has_sec_sig": 2.5,
    "signal_score_scale": 0.4,
    "signal_score_cap": 5.0,
    "hit_antimicr": 3.5,
    "hit_lytic": 2.0,
    "hit_protease": 1.0,
    "hit_lipase": 1.0,
    "hit_export": 1.0,
}


@dataclass
class StrainConfig:
    strain_name: str
    role: str
    fasta_path: Path
    gtf_path: Path
    cds_path: Path
    assembly_accession: str = "unspecified"
    annotation_pipeline: str = "unspecified"


@dataclass
class PipelineConfig:
    antagonist: StrainConfig
    target: StrainConfig
    output_root: Path = Path("runs")
    run_id: Optional[str] = None
    min_gtf_join_rate: float = 0.25
    hydrophobic_aas: frozenset = DEFAULT_HYDROPHOBIC_AAS
    min_hydrophobic_frac: float = 0.55
    cleavage_regex: str = r"[AGST].[AGST]"
    score_cap: float = 10.0
    high_candidate_threshold: float = 5.0
    weights: Dict[str, float] = field(default_factory=lambda: dict(DEFAULT_WEIGHTS))


def resolve_input_path(raw_path: str, base_dir: Path) -> Path:
    p = Path(raw_path)
    if p.is_absolute():
        return p
    cfg_rel = (base_dir / p).resolve()
    if cfg_rel.exists():
        return cfg_rel
    return (Path.cwd() / p).resolve()


def sha256_file(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def get_pkg_version(dist_name: str) -> str:
    try:
        return version(dist_name)
    except PackageNotFoundError:
        return "not-installed"


def setup_logging(run_dir: Path) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    log_path = run_dir / "pipeline.log"
    LOGGER.setLevel(logging.INFO)
    LOGGER.handlers.clear()

    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setFormatter(fmt)

    LOGGER.addHandler(sh)
    LOGGER.addHandler(fh)


def load_config(config_path: Optional[Path], cli_args: argparse.Namespace) -> PipelineConfig:
    raw: Dict[str, Any] = {}
    resolved_cfg_path: Optional[Path] = None

    if config_path is not None:
        resolved_cfg_path = Path(config_path).resolve()
        if not resolved_cfg_path.exists():
            raise FileNotFoundError(f"Config file not found: {resolved_cfg_path}")
    else:
        default_cfg = Path("config.yaml").resolve()
        if default_cfg.exists():
            resolved_cfg_path = default_cfg

    if resolved_cfg_path is not None:
        with resolved_cfg_path.open("r", encoding="utf-8") as f:
            loaded = yaml.safe_load(f) or {}
            if not isinstance(loaded, dict):
                raise ValueError(f"Config file root must be a mapping: {resolved_cfg_path}")
            raw = loaded

    base_dir = resolved_cfg_path.parent if resolved_cfg_path is not None else Path.cwd()
    strains_cfg = raw.get("strains", {}) if isinstance(raw.get("strains"), dict) else {}
    ant_cfg = strains_cfg.get("antagonist", {}) if isinstance(strains_cfg.get("antagonist"), dict) else {}
    tgt_cfg = strains_cfg.get("target", {}) if isinstance(strains_cfg.get("target"), dict) else {}

    def resolve_strain(
        section: Dict[str, Any], defaults: Dict[str, str], prefix: str
    ) -> StrainConfig:
        fasta_val = (
            getattr(cli_args, f"{prefix}_fasta")
            or section.get("fasta_path")
            or defaults["fasta"]
        )
        gtf_val = (
            getattr(cli_args, f"{prefix}_gtf")
            or section.get("gtf_path")
            or defaults["gtf"]
        )
        cds_val = (
            getattr(cli_args, f"{prefix}_cds")
            or section.get("cds_path")
            or defaults["cds"]
        )
        name_val = (
            getattr(cli_args, f"{prefix}_name")
            or section.get("strain_name")
            or defaults["name"]
        )

        return StrainConfig(
            strain_name=str(name_val),
            role=str(section.get("role", prefix)),
            assembly_accession=str(section.get("assembly_accession", "unspecified")),
            annotation_pipeline=str(section.get("annotation_pipeline", "unspecified")),
            fasta_path=resolve_input_path(str(fasta_val), base_dir),
            gtf_path=resolve_input_path(str(gtf_val), base_dir),
            cds_path=resolve_input_path(str(cds_val), base_dir),
        )

    ant = resolve_strain(
        ant_cfg,
        {
            "name": "Bacillus_pumilus",
            "fasta": "data/BacilluspGCF_000017885.4_ASM1788v4_genomic.fna",
            "gtf": "data/Bacilluspgenomic.gtf",
            "cds": "data/Bacilluspcds_from_genomic.fna",
        },
        "bp",
    )
    tgt = resolve_strain(
        tgt_cfg,
        {
            "name": "Staphylococcus_aureus_NCTC8325",
            "fasta": "data/StaphGCF_000013425.1_ASM1342v1_genomic.fna",
            "gtf": "data/Staphgenomic.gtf",
            "cds": "data/Staphcds_from_genomic.fna",
        },
        "sa",
    )

    sig_cfg = (
        raw.get("heuristic_signal", {})
        if isinstance(raw.get("heuristic_signal"), dict)
        else {}
    )
    hydro_list = sig_cfg.get("hydrophobic_aas", sorted(DEFAULT_HYDROPHOBIC_AAS))
    hydro_set = frozenset(str(aa).strip().upper() for aa in hydro_list if str(aa).strip())
    if not hydro_set:
        raise ValueError("heuristic_signal.hydrophobic_aas cannot be empty.")

    score_cfg = (
        raw.get("scoring", {}) if isinstance(raw.get("scoring"), dict) else {}
    )
    weights = dict(DEFAULT_WEIGHTS)
    if isinstance(score_cfg.get("weights"), dict):
        for k, v in score_cfg["weights"].items():
            if k in weights:
                weights[k] = float(v)

    raw_out = cli_args.output_dir or raw.get("output_dir") or "runs"
    output_root = resolve_input_path(str(raw_out), base_dir)
    run_id = (
        cli_args.run_id
        or raw.get("run_id")
        or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    )
    min_join = float(
        cli_args.min_join_rate
        if cli_args.min_join_rate is not None
        else raw.get("min_gtf_join_rate", 0.25)
    )

    return PipelineConfig(
        antagonist=ant,
        target=tgt,
        output_root=output_root,
        run_id=str(run_id),
        min_gtf_join_rate=min_join,
        hydrophobic_aas=hydro_set,
        min_hydrophobic_frac=float(sig_cfg.get("min_hydrophobic_frac", 0.55)),
        cleavage_regex=str(sig_cfg.get("cleavage_regex", r"[AGST].[AGST]")),
        score_cap=float(score_cfg.get("cap", 10.0)),
        high_candidate_threshold=float(
            score_cfg.get("high_candidate_threshold", 5.0)
        ),
        weights=weights,
    )


def validate_strain_files(strain_cfg: StrainConfig) -> None:
    for label, p in (
        ("FASTA", strain_cfg.fasta_path),
        ("GTF", strain_cfg.gtf_path),
        ("CDS", strain_cfg.cds_path),
    ):
        if not p.exists():
            raise FileNotFoundError(f"[{strain_cfg.strain_name}] Missing {label} file: {p}")
        if p.stat().st_size == 0:
            raise ValueError(f"[{strain_cfg.strain_name}] Empty {label} file: {p}")


def compute_atgc_gc_percent(seq_str: str) -> Tuple[float, int, int]:
    """
    Unified GC metric using valid A/T/G/C bases as denominator so genome and CDS
    GC percentages are directly comparable. Also returns n_count and atgc_count.
    """
    s = seq_str.upper()
    n_count = s.count("N")
    g_count = s.count("G")
    c_count = s.count("C")
    a_count = s.count("A")
    t_count = s.count("T")
    atgc_total = a_count + t_count + g_count + c_count
    if atgc_total == 0:
        return 0.0, n_count, 0
    gc_pct = round(((g_count + c_count) / atgc_total) * 100.0, 2)
    return gc_pct, n_count, atgc_total


def summarize_genome_fasta(fasta_path: Path, strain_label: str) -> Dict[str, Any]:
    LOGGER.info("[%s] Reading genome FASTA: %s", strain_label, fasta_path)
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {fasta_path}")

    lengths = sorted([len(rec.seq) for rec in records], reverse=True)
    total_len = sum(lengths)
    joined_upper = "".join(str(rec.seq).upper() for rec in records)
    gc_pct, n_count, atgc_total = compute_atgc_gc_percent(joined_upper)

    running = 0
    n50 = lengths[0]
    for l in lengths:
        running += l
        if running >= total_len * 0.5:
            n50 = l
            break

    return {
        "strain": strain_label,
        "n_contigs": len(records),
        "total_length_bp": total_len,
        "atgc_bases": atgc_total,
        "n_bases": n_count,
        "n50_bp": n50,
        "gc_percent_atgc": gc_pct,
    }


def parse_annotation_attributes(attr_str: str) -> Tuple[Dict[str, List[str]], List[str]]:
    """
    Parse GTF (`key "val";`) or GFF3 (`key=val;`) attributes while preserving
    repeated keys in lists and logging duplicate/ambiguous values.
    """
    attrs: Dict[str, List[str]] = {}
    warnings: List[str] = []

    token_pattern = re.compile(r'([A-Za-z0-9_.-]+)\s*(?:=\s*([^;]+)|"([^"]*)")')
    for m in token_pattern.finditer(attr_str):
        key = m.group(1)
        val = (m.group(2) if m.group(2) is not None else m.group(3) or "").strip()
        if not key:
            continue
        attrs.setdefault(key, [])
        if val and val not in attrs[key]:
            attrs[key].append(val)

    for k, vals in attrs.items():
        if len(vals) > 1:
            warnings.append(f"Repeated attribute key '{k}' with multiple values: {vals}")

    return attrs, warnings


def normalize_feature_id(raw_id: str) -> str:
    """
    Normalize common CDS/GTF/GFF IDs (e.g., stripping 'lcl|', 'cds_', 'gene-',
    'transcript-', or trailing exon/CDS ordinal suffixes) for cross-file joins.
    """
    s = raw_id.strip()
    if "|" in s:
        s = s.split("|")[-1]
    s = re.sub(r"^(gene-|transcript-|cds-|rna-)", "", s, flags=re.IGNORECASE)
    s = re.sub(r"_cds_[^_]+_\d+$", "", s)
    return s.strip()


def _to_scalar_str(val: Any) -> str:
    if isinstance(val, (list, tuple)):
        return str(val[0]) if len(val) > 0 and val[0] is not None else ""
    if val is None or pd.isna(val):
        return ""
    return str(val)


def load_gtf_cds_lookup(gtf_path: Path, strain_label: str) -> Tuple[pd.DataFrame, int]:
    """
    Return a lookup DataFrame with columns ['match_key', 'gtf_product'] (never
    including 'strain' so it cannot collide during merge) and raw CDS feature count.
    """
    records_map: Dict[str, List[str]] = {}
    cds_row_count = 0
    dup_warnings = 0

    with gtf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2].strip()
            if feature != "CDS":
                continue

            cds_row_count += 1
            attrs, warns = parse_annotation_attributes(parts[8])
            if warns:
                dup_warnings += len(warns)

            prod_candidates: List[str] = []
            for p_key in ("product", "Name", "gene_name", "gene", "ID", "locus_tag"):
                prod_candidates.extend(attrs.get(p_key, []))
            product_str = "; ".join(dict.fromkeys(p for p in prod_candidates if p))

            id_candidates: List[str] = []
            for i_key in ("transcript_id", "protein_id", "ID", "locus_tag", "gene_id", "Parent"):
                id_candidates.extend(attrs.get(i_key, []))

            for raw_id in id_candidates:
                norm_id = normalize_feature_id(raw_id)
                if norm_id:
                    records_map.setdefault(norm_id, [])
                    if product_str and product_str not in records_map[norm_id]:
                        records_map[norm_id].append(product_str)

    if dup_warnings:
        LOGGER.warning(
            "[%s] Encountered %d repeated attribute-key warnings in %s",
            strain_label,
            dup_warnings,
            gtf_path,
        )

    rows = [
        {"match_key": _to_scalar_str(k), "gtf_product": "; ".join(v)}
        for k, v in records_map.items()
        if _to_scalar_str(k) != ""
    ]
    if not rows:
        return pd.DataFrame(columns=["match_key", "gtf_product"]), cds_row_count

    return pd.DataFrame(rows).drop_duplicates(subset=["match_key"], keep="first"), cds_row_count


def predict_sec_signal_heuristic(
    aa_seq: str,
    hydrophobic_aas: frozenset = DEFAULT_HYDROPHOBIC_AAS,
    min_hydrophobic_frac: float = 0.55,
    cleavage_regex: str = r"[AGST].[AGST]",
) -> Dict[str, Any]:
    """
    Heuristic Sec-like signal peptide detector (candidate ranking only; benchmark
    against SignalP 6.0 / DeepTMHMM before drawing cross-strain conclusions).
    """
    if not aa_seq or len(aa_seq) < 25:
        return {"has_sec_sig": False, "signal_score": 0.0, "cleavage_aa": None}

    n_region = aa_seq[0:8]
    h_window = aa_seq[5:25]
    c_window = aa_seq[14:35]

    has_m = aa_seq[0] == "M"
    pos_charge = sum(1 for aa in n_region if aa in {"K", "R", "H"})
    hydrophobic_frac = sum(1 for aa in h_window if aa in hydrophobic_aas) / len(h_window)
    cleavage_match = re.search(cleavage_regex, c_window)
    has_cleavage = cleavage_match is not None

    score = (
        (1.0 if has_m else 0.0)
        + min(pos_charge, 3) * 0.5
        + (hydrophobic_frac * 4.0)
        + (1.5 if has_cleavage else 0.0)
    )
    cleavage_aa = (14 + cleavage_match.end()) if cleavage_match else None
    has_signal = bool(
        has_m and pos_charge >= 1 and hydrophobic_frac >= min_hydrophobic_frac and has_cleavage
    )

    return {
        "has_sec_sig": has_signal,
        "signal_score": round(score, 2),
        "cleavage_aa": cleavage_aa,
    }


ANTAGONISM_PATTERNS = {
    "hit_lytic": re.compile(
        r"\b(?:lysozyme|muramidase|peptidoglycan\s+(?:hydrolase|amidase|endopeptidase)|chitinase|chitosanase|glucanase|lytic\s+(?:transglycosylase|murein))\b",
        re.IGNORECASE,
    ),
    "hit_protease": re.compile(
        r"\b(?:protease|proteinase|peptidase|subtilisin|metalloprotease|serine\s+protease|aminopeptidase|carboxypeptidase)\b",
        re.IGNORECASE,
    ),
    "hit_lipase": re.compile(
        r"\b(?:lipase|phospholipase|esterase|cutinase)\b",
        re.IGNORECASE,
    ),
    "hit_antimicr": re.compile(
        r"\b(?:bacteriocin|lantibiotic|subtilin|bacilysin|surfactin|fengycin|iturin|siderophore|polyketide\s+synthase|non-?ribosomal\s+peptide|antibiotic|antimicrobial|antifungal|antagonis(?:t|tic))\b",
        re.IGNORECASE,
    ),
    "hit_export": re.compile(
        r"\b(?:secA|secY|secE|secG|tatA|tatB|tatC|signal\s+peptidase|abc\s+transporter|efflux\s+pump|exoprotein)\b",
        re.IGNORECASE,
    ),
}


def extract_all_cds_keys(record: Any) -> List[str]:
    desc = record.description or record.id
    keys: List[str] = []
    for tag in ("protein_id", "locus_tag", "transcript_id"):
        m = re.search(rf"\[{tag}=([^\]]+)\]", desc)
        if m:
            norm = normalize_feature_id(m.group(1))
            if norm and norm not in keys:
                keys.append(norm)
    fallback = normalize_feature_id(record.id.split()[0])
    if fallback and fallback not in keys:
        keys.append(fallback)
    return keys


def mine_cds_features(
    strain_cfg: StrainConfig,
    gtf_lookup: pd.DataFrame,
    gtf_cds_count: int,
    cfg: PipelineConfig,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    LOGGER.info("[%s] Reading CDS FASTA: %s", strain_cfg.strain_name, strain_cfg.cds_path)
    records = list(SeqIO.parse(str(strain_cfg.cds_path), "fasta"))
    rows: List[Dict[str, Any]] = []
    trimmed_cds_count = 0
    internal_stop_cds_count = 0

    for rec in records:
        nt_seq = str(rec.seq).upper()
        nt_len = len(nt_seq)
        remainder = nt_len % 3
        was_trimmed = remainder != 0
        if was_trimmed:
            trimmed_cds_count += 1

        trimmed_len = max(nt_len - remainder, 0)
        codon_seq = nt_seq[:trimmed_len]
        aa_raw = str(Seq(codon_seq).translate(table=11, to_stop=False)) if trimmed_len >= 3 else ""
        has_internal_stop = "*" in aa_raw.rstrip("*")
        if has_internal_stop:
            internal_stop_cds_count += 1

        aa_clean = aa_raw.rstrip("*").replace("*", "")
        gc_pct, n_bases, atgc_bases = compute_atgc_gc_percent(codon_seq)
        sig_info = predict_sec_signal_heuristic(
            aa_clean,
            hydrophobic_aas=cfg.hydrophobic_aas,
            min_hydrophobic_frac=cfg.min_hydrophobic_frac,
            cleavage_regex=cfg.cleavage_regex,
        )

        cand_keys = extract_all_cds_keys(rec)
        primary_key = (
            _to_scalar_str(cand_keys[0])
            if cand_keys
            else _to_scalar_str(normalize_feature_id(rec.id.split()[0]))
        )

        rows.append(
            {
                "strain": strain_cfg.strain_name,
                "role": strain_cfg.role,
                "cds_id": _to_scalar_str(rec.id.split()[0]),
                "join_key": primary_key,
                "header": rec.description or rec.id,
                "cds_nt_len": nt_len,
                "trimmed_bp": remainder,
                "was_trimmed": was_trimmed,
                "has_internal_stop": has_internal_stop,
                "prot_aa_len": len(aa_clean),
                "cds_gc_pct_atgc": gc_pct,
                "cds_atgc_bases": atgc_bases,
                "cds_n_bases": n_bases,
                "has_sec_sig": sig_info["has_sec_sig"],
                "signal_score": sig_info["signal_score"],
                "cleavage_aa": sig_info["cleavage_aa"],
            }
        )

    df = pd.DataFrame(rows)
    total_cds = len(df)
    join_rate = 0.0

    if not df.empty:
        df["join_key"] = df["join_key"].map(_to_scalar_str)
        lookup_clean = gtf_lookup.copy()
        lookup_clean["match_key"] = lookup_clean["match_key"].map(_to_scalar_str)

        df = df.merge(lookup_clean, how="left", left_on="join_key", right_on="match_key")
        matched_count = int(df["gtf_product"].notna().sum())
        join_rate = round(matched_count / total_cds, 4) if total_cds > 0 else 0.0

        LOGGER.info(
            "[%s] CDS=%d | GTF_CDS_rows=%d | Matched=%d (%.1f%%) | Trimmed=%d | Internal_stops=%d",
            strain_cfg.strain_name,
            total_cds,
            gtf_cds_count,
            matched_count,
            join_rate * 100.0,
            trimmed_cds_count,
            internal_stop_cds_count,
        )
        if join_rate < cfg.min_gtf_join_rate:
            sample_cds_keys = [k for k in df["join_key"].unique().tolist() if k][:8]
            sample_gtf_keys = [k for k in lookup_clean["match_key"].unique().tolist() if k][:8]
            msg = (
                f"[{strain_cfg.strain_name}] GTF join rate {join_rate:.1%} is below minimum "
                f"threshold {cfg.min_gtf_join_rate:.1%}. "
                f"Sample CDS join_keys={sample_cds_keys}; Sample GTF match_keys={sample_gtf_keys}."
            )
            raise RuntimeError(msg)

        df["gtf_product"] = df["gtf_product"].fillna("")
        df["combined_annot"] = (
            (df["header"] + " " + df["gtf_product"])
            .str.replace(r"\s+", " ", regex=True)
            .str.strip()
        )

        for col_name, pattern in ANTAGONISM_PATTERNS.items():
            df[col_name] = df["combined_annot"].apply(lambda text: bool(pattern.search(text)))

        df["n_antag_hits"] = (
            df["hit_lytic"].astype(int)
            + df["hit_protease"].astype(int)
            + df["hit_lipase"].astype(int)
            + df["hit_antimicr"].astype(int)
        )

        w = cfg.weights
        raw_score = (
            df["has_sec_sig"].map({True: w["has_sec_sig"], False: 0.0})
            + df["signal_score"].clip(upper=w["signal_score_cap"]) * w["signal_score_scale"]
            + df["hit_antimicr"].map({True: w["hit_antimicr"], False: 0.0})
            + df["hit_lytic"].map({True: w["hit_lytic"], False: 0.0})
            + df["hit_protease"].map({True: w["hit_protease"], False: 0.0})
            + df["hit_lipase"].map({True: w["hit_lipase"], False: 0.0})
            + df["hit_export"].map({True: w["hit_export"], False: 0.0})
        )
        df["heuristic_ranking_score"] = raw_score.clip(upper=cfg.score_cap).round(2)

        df = df.sort_values(
            by=["heuristic_ranking_score", "signal_score", "prot_aa_len", "cds_id"],
            ascending=[False, False, False, True],
            kind="mergesort",
        ).reset_index(drop=True)

    stats = {
        "strain": strain_cfg.strain_name,
        "total_cds": total_cds,
        "gtf_cds_rows": gtf_cds_count,
        "gtf_join_rate": join_rate,
        "trimmed_cds_count": trimmed_cds_count,
        "internal_stop_cds_count": internal_stop_cds_count,
    }
    return df, stats


def build_strain_summary(combined_cds: pd.DataFrame, threshold: float) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for strain_label, grp in combined_cds.groupby("strain", sort=False):
        n_cds = len(grp)
        scale = 1000.0 / n_cds if n_cds > 0 else 0.0

        sec_n = int(grp["has_sec_sig"].sum())
        lytic_n = int(grp["hit_lytic"].sum())
        prot_n = int(grp["hit_protease"].sum())
        lip_n = int(grp["hit_lipase"].sum())
        anti_n = int(grp["hit_antimicr"].sum())
        exp_n = int(grp["hit_export"].sum())
        high_n = int((grp["heuristic_ranking_score"] >= threshold).sum())

        rows.append(
            {
                "strain": strain_label,
                "total_cds": n_cds,
                "mean_cds_gc_pct_atgc": round(float(grp["cds_gc_pct_atgc"].mean()), 2) if n_cds else 0.0,
                "internal_stop_cds": int(grp["has_internal_stop"].sum()),
                "trimmed_cds": int(grp["was_trimmed"].sum()),
                "sec_sig_count": sec_n,
                "sec_sig_per_1k_cds": round(sec_n * scale, 2),
                "lytic_count": lytic_n,
                "lytic_per_1k_cds": round(lytic_n * scale, 2),
                "protease_count": prot_n,
                "protease_per_1k_cds": round(prot_n * scale, 2),
                "lipase_count": lip_n,
                "lipase_per_1k_cds": round(lip_n * scale, 2),
                "antimicrobial_count": anti_n,
                "antimicrobial_per_1k_cds": round(anti_n * scale, 2),
                "export_count": exp_n,
                "export_per_1k_cds": round(exp_n * scale, 2),
                "high_rank_candidates": high_n,
                "high_rank_per_1k_cds": round(high_n * scale, 2),
                "max_heuristic_ranking_score": round(float(grp["heuristic_ranking_score"].max()), 2) if n_cds else 0.0,
            }
        )
    return pd.DataFrame(rows)


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Heuristic genomic & secretion candidate-ranking pipeline (Bacillus pumilus vs. target)."
    )
    parser.add_argument("--config", type=Path, default=None, help="Path to YAML config file.")
    parser.add_argument("--output-dir", type=str, default=None, help="Base directory for run outputs.")
    parser.add_argument("--run-id", type=str, default=None, help="Explicit run identifier.")
    parser.add_argument("--min-join-rate", type=float, default=None, help="Minimum required CDS-to-GTF join fraction.")

    parser.add_argument("--bp-name", type=str, default=None)
    parser.add_argument("--bp-fasta", type=str, default=None)
    parser.add_argument("--bp-gtf", type=str, default=None)
    parser.add_argument("--bp-cds", type=str, default=None)

    parser.add_argument("--sa-name", type=str, default=None)
    parser.add_argument("--sa-fasta", type=str, default=None)
    parser.add_argument("--sa-gtf", type=str, default=None)
    parser.add_argument("--sa-cds", type=str, default=None)
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)
    cfg = load_config(args.config, args)
    validate_strain_files(cfg.antagonist)
    validate_strain_files(cfg.target)

    run_dir = cfg.output_root / str(cfg.run_id)
    setup_logging(run_dir)

    LOGGER.info("Starting run_id=%s writing to %s", cfg.run_id, run_dir)

    genome_summary = pd.DataFrame(
        [
            summarize_genome_fasta(cfg.antagonist.fasta_path, cfg.antagonist.strain_name),
            summarize_genome_fasta(cfg.target.fasta_path, cfg.target.strain_name),
        ]
    )

    bp_lookup, bp_gtf_cds = load_gtf_cds_lookup(cfg.antagonist.gtf_path, cfg.antagonist.strain_name)
    sa_lookup, sa_gtf_cds = load_gtf_cds_lookup(cfg.target.gtf_path, cfg.target.strain_name)

    bp_mined, bp_stats = mine_cds_features(cfg.antagonist, bp_lookup, bp_gtf_cds, cfg)
    sa_mined, sa_stats = mine_cds_features(cfg.target, sa_lookup, sa_gtf_cds, cfg)

    combined_cds = pd.concat([bp_mined, sa_mined], ignore_index=True)
    strain_summary = build_strain_summary(combined_cds, cfg.high_candidate_threshold)

    genome_summary.to_csv(run_dir / "genome_summary.csv", index=False)
    strain_summary.to_csv(run_dir / "strain_summary_per_1k_cds.csv", index=False)
    bp_mined.to_csv(run_dir / f"{cfg.antagonist.strain_name}_ranked_candidates.csv", index=False)
    sa_mined.to_csv(run_dir / f"{cfg.target.strain_name}_ranked_candidates.csv", index=False)

    plt.figure(figsize=(8.5, 5))
    for strain_label, grp in combined_cds.groupby("strain", sort=False):
        plt.hist(
            grp["heuristic_ranking_score"],
            bins=20,
            range=(0, cfg.score_cap),
            alpha=0.6,
            label=strain_label,
        )
    plt.title(
        "Heuristic Candidate-Ranking Score (Not a Biological Measurement)\n"
        "Compositional GC differences confound cross-strain secretion counts without SignalP 6.0"
    )
    plt.xlabel(f"Heuristic Candidate-Ranking Score (0-{cfg.score_cap:g}, uncalibrated)")
    plt.ylabel("CDS Count")
    plt.legend()
    plt.tight_layout()
    plt.savefig(run_dir / "heuristic_ranking_distribution.png", dpi=300)
    plt.close()

    manifest = {
        "run_id": cfg.run_id,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "packages": {
            "biopython": get_pkg_version("biopython"),
            "pandas": get_pkg_version("pandas"),
            "matplotlib": get_pkg_version("matplotlib"),
            "PyYAML": get_pkg_version("PyYAML"),
        },
        "parameters": {
            "min_gtf_join_rate": cfg.min_gtf_join_rate,
            "hydrophobic_aas": sorted(cfg.hydrophobic_aas),
            "min_hydrophobic_frac": cfg.min_hydrophobic_frac,
            "cleavage_regex": cfg.cleavage_regex,
            "score_cap": cfg.score_cap,
            "high_candidate_threshold": cfg.high_candidate_threshold,
            "weights": cfg.weights,
        },
        "inputs": {
            "antagonist": {
                "strain_name": cfg.antagonist.strain_name,
                "assembly_accession": cfg.antagonist.assembly_accession,
                "annotation_pipeline": cfg.antagonist.annotation_pipeline,
                "fasta_path": str(cfg.antagonist.fasta_path),
                "fasta_sha256": sha256_file(cfg.antagonist.fasta_path),
                "gtf_path": str(cfg.antagonist.gtf_path),
                "gtf_sha256": sha256_file(cfg.antagonist.gtf_path),
                "cds_path": str(cfg.antagonist.cds_path),
                "cds_sha256": sha256_file(cfg.antagonist.cds_path),
                "mining_stats": bp_stats,
            },
            "target": {
                "strain_name": cfg.target.strain_name,
                "assembly_accession": cfg.target.assembly_accession,
                "annotation_pipeline": cfg.target.annotation_pipeline,
                "fasta_path": str(cfg.target.fasta_path),
                "fasta_sha256": sha256_file(cfg.target.fasta_path),
                "gtf_path": str(cfg.target.gtf_path),
                "gtf_sha256": sha256_file(cfg.target.gtf_path),
                "cds_path": str(cfg.target.cds_path),
                "cds_sha256": sha256_file(cfg.target.cds_path),
                "mining_stats": sa_stats,
            },
        },
    }
    with (run_dir / "run_manifest.json").open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    LOGGER.info("Completed run_id=%s successfully.", cfg.run_id)


if __name__ == "__main__":
    main()