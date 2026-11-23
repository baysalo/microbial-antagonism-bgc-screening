from pathlib import Path
import re
import pandas as pd
import yaml
from Bio import SeqIO
from Bio.Seq import Seq

runs_root = Path("runs")
latest_run = sorted(
    [p for p in runs_root.iterdir() if p.is_dir()],
    key=lambda p: p.name,
    reverse=True,
)[0]

cfg_path = Path("config.yaml")
cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8")) if cfg_path.exists() else {}
bp_cds_raw = (
    cfg.get("strains", {})
    .get("antagonist", {})
    .get("cds_path", "data/Bacilluspcds_from_genomic.fna")
)
bp_cds_path = Path(bp_cds_raw)
if not bp_cds_path.is_absolute() and not bp_cds_path.exists():
  bp_cds_path = (cfg_path.resolve().parent / bp_cds_path).resolve()

ranked_df = pd.read_csv(latest_run / "Bacillus_pumilus_ranked_candidates.csv")

housekeeping_re = re.compile(
    r"\b(?:lepB|sipW|sppA|spoIVB|htpX|ctpB|ddcA|signal peptidase|signal peptide peptidase|DNA damage checkpoint)\b",
    re.IGNORECASE,
)

filtered = ranked_df[
    (~ranked_df["has_internal_stop"].astype(bool))
    & (ranked_df["has_sec_sig"].astype(bool))
    & (
        ranked_df["hit_antimicr"].astype(bool)
        | ranked_df["hit_lytic"].astype(bool)
        | ranked_df["hit_protease"].astype(bool)
        | ranked_df["hit_lipase"].astype(bool)
    )
    & (~ranked_df["combined_annot"].fillna("").str.contains(housekeeping_re))
].copy()

target_ids = set(filtered["cds_id"].astype(str))
cds_records = {}
for rec in SeqIO.parse(str(bp_cds_path), "fasta"):
  cid = rec.id.split()[0]
  if cid in target_ids:
    nt_seq = str(rec.seq).upper()
    trimmed_len = (len(nt_seq) // 3) * 3
    aa_seq = (
        str(Seq(nt_seq[:trimmed_len]).translate(table=11, to_stop=False))
        .rstrip("*")
        .replace("*", "")
    )
    cds_records[cid] = aa_seq

out_faa = latest_run / "Bacillus_pumilus_secreted_antagonism_candidates.faa"
written = 0
with out_faa.open("w", encoding="utf-8") as fh:
  for _, row in filtered.iterrows():
    cid = str(row["cds_id"])
    aa_seq = cds_records.get(cid, "")
    if not aa_seq:
      continue
    desc = (
        f"score={row['heuristic_ranking_score']:.2f} "
        f"sig={row['signal_score']:.2f} "
        f"anti={int(bool(row['hit_antimicr']))} "
        f"lytic={int(bool(row['hit_lytic']))} "
        f"prot={int(bool(row['hit_protease']))} "
        f"lip={int(bool(row['hit_lipase']))}"
    )
    fh.write(f">{cid} {desc}\n")
    for i in range(0, len(aa_seq), 60):
      fh.write(aa_seq[i : i + 60] + "\n")
    written += 1

print(f"Exported {written} protein sequences to: {out_faa.resolve()}")