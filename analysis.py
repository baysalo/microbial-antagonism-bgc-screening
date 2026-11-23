import json
import re
from pathlib import Path
from typing import Dict, Optional
import pandas as pd
import yaml
from Bio import SeqIO
from Bio.Seq import Seq

# 1. Son çalıştırma klasörünü ve SignalP 6.0 output.json dosyasını bul
runs_root = Path("runs")
latest_run = sorted(
    [p for p in runs_root.iterdir() if p.is_dir()],
    key=lambda p: p.name,
    reverse=True,
)[0]

signalp_json_path = Path("output.json")
if not signalp_json_path.exists():
  signalp_json_path = latest_run / "output.json"

sig_data = json.loads(signalp_json_path.read_text(encoding="utf-8"))
seq_entries: Dict[str, dict] = sig_data.get("SEQUENCES", {})

# 2. Config'den Bacillus pumilus CDS FASTA yolunu al
cfg_path = Path("config.yaml")
cfg = (
    yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    if cfg_path.exists()
    else {}
)
bp_cds_raw = (
    cfg.get("strains", {})
    .get("antagonist", {})
    .get("cds_path", "data/Bacilluspcds_from_genomic.fna")
)
bp_cds_path = Path(bp_cds_raw)
if not bp_cds_path.is_absolute() and not bp_cds_path.exists():
  bp_cds_path = (cfg_path.resolve().parent / bp_cds_path).resolve()

# 3. CDS ID -> Çevrilmiş tam AA dizisi sözlüğü oluştur
aa_by_cds: Dict[str, str] = {}
for rec in SeqIO.parse(str(bp_cds_path), "fasta"):
  cid = rec.id.split()[0]
  nt_seq = str(rec.seq).upper()
  trimmed_len = (len(nt_seq) // 3) * 3
  aa_seq = (
      str(Seq(nt_seq[:trimmed_len]).translate(table=11, to_stop=False))
      .rstrip("*")
      .replace("*", "")
  )
  aa_by_cds[cid] = aa_seq


def parse_cleavage_pos(cs_pos_str: str) -> Optional[int]:
  """
    'Cleavage site between pos. 25 and 26. Probability 0.979395' -> 25
    (1-based N-terminal signal peptide length to trim).
    """
  if not cs_pos_str:
    return None
  m = re.search(r"between\s+pos\.\s*(\d+)\s+and\s+(\d+)", cs_pos_str)
  if m:
    return int(m.group(1))
  return None


def extract_cds_id_from_sig_name(sig_name: str) -> Optional[str]:
  """
    SignalP names replace '|' with '_' and append 'score_...'.
    Matches against known CDS IDs cleanly.
    """
  norm_sig = sig_name.replace("|", "_")
  for cid in aa_by_cds:
    if norm_sig.startswith(cid.replace("|", "_")):
      return cid
  m = re.match(r"^(lcl_[^_\s]+_cds_[^_\s]+_\d+)", norm_sig)
  if m:
    candidate = m.group(1)
    for cid in aa_by_cds:
      if cid.replace("|", "_") == candidate:
        return cid
  return None


def write_fasta(records_list, out_path: Path) -> None:
  with out_path.open("w", encoding="utf-8") as fh:
    for header, seq in records_list:
      fh.write(f">{header}\n")
      for i in range(0, len(seq), 60):
        fh.write(seq[i : i + 60] + "\n")


mature_records = []
full_for_dtm_records = []
summary_rows = []

for sig_name, info in seq_entries.items():
  cid = extract_cds_id_from_sig_name(sig_name)
  if not cid or cid not in aa_by_cds:
    continue

  full_aa = aa_by_cds[cid]
  pred_type = info.get("Prediction", "Unknown")
  cs_str = info.get("CS_pos", "")
  cs_len = parse_cleavage_pos(cs_str)

  if cs_len is not None and 0 < cs_len < len(full_aa):
    mature_aa = full_aa[cs_len:]
    cleaved = True
    start_1based = cs_len + 1
  else:
    mature_aa = full_aa
    cleaved = False
    start_1based = 1

  mature_header = (
      f"{cid} pred={pred_type.replace(' ', '_')} "
      f"cleaved={int(cleaved)} mature_aa_len={len(mature_aa)} "
      f"cs_pos=\"{cs_str}\""
  )
  mature_records.append((mature_header, mature_aa))

  dtm_header = (
      f"{cid} signalp={pred_type.replace(' ', '_')} "
      f"mature_start={start_1based} len={len(full_aa)}"
  )
  full_for_dtm_records.append((dtm_header, full_aa))

  summary_rows.append({
      "cds_id": cid,
      "signalp_prediction": pred_type,
      "cleavage_pos": cs_str,
      "cleaved_signal_aa": cs_len if cleaved else 0,
      "full_aa_len": len(full_aa),
      "mature_aa_len": len(mature_aa),
  })

mature_faa = latest_run / "Bacillus_pumilus_signalp6_mature_candidates.faa"
dtm_faa = latest_run / "Bacillus_pumilus_deeptmhmm_input.faa"
summary_csv = latest_run / "signalp6_mature_summary.csv"

write_fasta(mature_records, mature_faa)
write_fasta(full_for_dtm_records, dtm_faa)
pd.DataFrame(summary_rows).to_csv(summary_csv, index=False)

print(f"Saved mature FASTA ({len(mature_records)} seqs): {mature_faa.resolve()}")
print(
    f"Saved DeepTMHMM input FASTA ({len(full_for_dtm_records)} seqs):"
    f" {dtm_faa.resolve()}"
)
print(f"Saved summary CSV: {summary_csv.resolve()}")