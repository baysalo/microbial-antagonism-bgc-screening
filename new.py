from pathlib import Path
import pandas as pd

runs_root = Path("runs")
latest_run = sorted(
    [p for p in runs_root.iterdir() if p.is_dir()],
    key=lambda p: p.name,
    reverse=True,
)[0]

df = pd.read_csv(latest_run / "Bacillus_pumilus_ranked_candidates.csv")

intact_hits = df[
    (~df["has_internal_stop"].astype(bool))
    & (
        df["hit_lytic"].astype(bool)
        | df["hit_protease"].astype(bool)
        | df["hit_lipase"].astype(bool)
        | df["hit_antimicr"].astype(bool)
    )
].sort_values(
    by=["heuristic_ranking_score", "signal_score", "prot_aa_len"],
    ascending=[False, False, False],
)

cols = [
    "cds_id",
    "prot_aa_len",
    "has_sec_sig",
    "signal_score",
    "hit_antimicr",
    "hit_lytic",
    "hit_protease",
    "hit_lipase",
    "heuristic_ranking_score",
    "combined_annot",
]
print(intact_hits[cols].head(25).to_string(index=False))