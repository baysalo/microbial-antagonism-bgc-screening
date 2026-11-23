import re
from pathlib import Path
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

runs_root = Path("runs")
latest_run = sorted(
    [p for p in runs_root.iterdir() if p.is_dir()],
    key=lambda p: p.name,
    reverse=True,
)[0]

ranked_df = pd.read_csv(latest_run / "Bacillus_pumilus_ranked_candidates.csv")
merged_sp_dtm = pd.read_csv(latest_run / "signalp6_deeptmhmm_merged.csv")

antismash_regions = [
    ("Region 1", 322101, 405798, "NRPS"),
    ("Region 2", 650067, 731023, "NRPS, T1PKS"),
    ("Region 3", 917306, 938212, "RRE-containing"),
    ("Region 4", 1086948, 1124600, "terpene, NI-siderophore"),
    ("Region 5", 1815104, 1843517, "betalactone"),
    ("Region 6", 1910959, 1932834, "terpene"),
    ("Region 7", 1972259, 2013360, "T3PKS"),
    ("Region 8", 2182287, 2203178, "terpene-precursor"),
    ("Region 9", 2471059, 2503512, "betalactone"),
    ("Region 10", 2828509, 2849433, "terpene-precursor"),
    ("Region 11", 3336263, 3346554, "RiPP-like"),
    ("Region 12", 3372349, 3413771, "other (bacilysin)"),
]
regions_df = pd.DataFrame(
    antismash_regions, columns=["region", "start", "end", "product"]
)


def extract_wp(val: str) -> str:
  m = re.search(r"(WP_\d+\.\d+)", str(val))
  return m.group(1) if m else str(val)


def parse_loc(annot: str):
  m = re.search(r"\[location=(?:complement\()?(\d+)\.\.(\d+)\)?\]", str(annot))
  if not m:
    return np.nan, np.nan
  a, b = int(m.group(1)), int(m.group(2))
  return min(a, b), max(a, b)


def short_label(row) -> str:
  wp = row["wp_id"]
  annot = str(row["combined_annot"]).lower()
  if "lci" in annot:
    desc = "Antimicrobial peptide LCI"
  elif "beta-glucanase" in annot or "bgls" in annot:
    desc = "Beta-glucanase (bglS)"
  elif "siderophore" in annot:
    desc = "Siderophore ABC-SBP"
  elif "uberolysin" in annot or "circular bacteriocin" in annot:
    desc = "Circular bacteriocin (RiPP)"
  elif "bacilysin" in annot:
    desc = "Bacilysin synthase (BacA)"
  elif "non-ribosomal peptide synthetase" in annot:
    desc = "NRPS mega-synthase"
  elif "type i polyketide synthase" in annot:
    desc = "Type I PKS"
  elif "d-alanyl-d-alanine carboxypeptidase" in annot:
    desc = "D-Ala-D-Ala carboxypeptidase"
  elif "lytic transglycosylase" in annot:
    desc = "Lytic transglycosylase/C40"
  elif "sepm" in annot:
    desc = "SepM serine protease"
  elif "trypsin-like" in annot:
    desc = "Trypsin-like serine peptidase"
  elif "zinc-dependent metalloprotease" in annot:
    desc = "Zn-dependent metalloprotease"
  elif "s8 family" in annot:
    desc = "S8 family peptidase"
  elif "c40 family" in annot:
    desc = "C40 family peptidase"
  else:
    desc = "Candidate effector"
  return f"{desc} ({wp})"


ranked_df["wp_id"] = ranked_df["cds_id"].apply(extract_wp)
merged_sp_dtm["wp_id"] = merged_sp_dtm["seq_key"].apply(extract_wp)

coords = ranked_df["combined_annot"].apply(parse_loc)
ranked_df["nt_start"] = [c[0] for c in coords]
ranked_df["nt_end"] = [c[1] for c in coords]

bgc_core_wps = {
    "WP_012008887.1",
    "WP_012009177.1",
    "WP_012009179.1",
    "WP_012011628.1",
}

sub_df = ranked_df[
    ranked_df["wp_id"].isin(set(merged_sp_dtm["wp_id"]))
    | ranked_df["wp_id"].isin(bgc_core_wps)
].copy()

sub_df = sub_df.drop_duplicates(subset=["wp_id"]).merge(
    merged_sp_dtm[[
        "wp_id",
        "signalp_prediction",
        "cleavage_after_aa",
        "deeptmhmm_type",
        "n_tmr",
        "soluble_or_surface_no_tmr",
    ]].drop_duplicates(subset=["wp_id"]),
    on="wp_id",
    how="left",
)


def map_region(row) -> str:
  s, e = row["nt_start"], row["nt_end"]
  if pd.isna(s) or pd.isna(e):
    return "Non-BGC"
  for _, r in regions_df.iterrows():
    if not (e < r["start"] or s > r["end"]):
      return f"{r['region']} ({r['product']})"
    if abs(s - r["start"]) <= 25000 or abs(e - r["end"]) <= 25000:
      return f"Near {r['region']}"
  return "Non-BGC (Secreted/Soluble)"


sub_df["bgc_context"] = sub_df.apply(map_region, axis=1)
sub_df["label"] = sub_df.apply(short_label, axis=1)
sub_df["signalp_prediction"] = sub_df["signalp_prediction"].fillna("Other (BGC)")
sub_df["deeptmhmm_type"] = sub_df["deeptmhmm_type"].fillna("Intracellular")
sub_df["n_tmr"] = sub_df["n_tmr"].fillna(0).astype(int)
sub_df["soluble_or_surface_no_tmr"] = (
    sub_df["soluble_or_surface_no_tmr"].fillna(False).astype(bool)
)


def assign_class(row) -> str:
  if (
      "Region 1 " in row["bgc_context"]
      or "Region 2" in row["bgc_context"]
      or "Region 12" in row["bgc_context"]
      or row["wp_id"] in bgc_core_wps
  ):
    return "BGC Core (NRPS / T1PKS / Bacilysin)"
  if "Region 11" in row["bgc_context"]:
    return "BGC RiPP-like (Circular Bacteriocin)"
  if row["signalp_prediction"] == "Lipoprotein signal peptide (Sec/SPII)":
    return "Sec/SPII Lipoprotein (0 TMR, Near Reg. 4)"
  if row["soluble_or_surface_no_tmr"]:
    return "Sec/SPI Soluble / Surface (0 TMR)"
  return "Other / Membrane-associated (TM/SP)"


palette = {
    "Sec/SPI Soluble / Surface (0 TMR)": "#1b7837",
    "Sec/SPII Lipoprotein (0 TMR, Near Reg. 4)": "#2166ac",
    "BGC Core (NRPS / T1PKS / Bacilysin)": "#b2182b",
    "BGC RiPP-like (Circular Bacteriocin)": "#762a83",
    "Other / Membrane-associated (TM/SP)": "#878787",
}

sub_df["category"] = sub_df.apply(assign_class, axis=1)
sub_df = sub_df.sort_values(
    by=["heuristic_ranking_score", "nt_start"], ascending=[True, False]
).reset_index(drop=True)

plt.rcParams["font.sans-serif"] = "DejaVu Sans"
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

fig = plt.figure(figsize=(15.5, 8.8), dpi=300)
gs = fig.add_gridspec(1, 2, width_ratios=[1.32, 1.08], wspace=0.34)

# --- PANEL A ---
ax1 = fig.add_subplot(gs[0, 0])
y_pos = np.arange(len(sub_df))
colors = [palette[c] for c in sub_df["category"]]

ax1.barh(
    y_pos,
    sub_df["heuristic_ranking_score"],
    color=colors,
    height=0.66,
    edgecolor="#222222",
    linewidth=0.5,
)

ax1.set_yticks(y_pos)
ax1.set_yticklabels(sub_df["label"], fontsize=8.3)
ax1.set_xlabel(
    "Heuristic Antagonism Ranking Score", fontsize=9.5, fontweight="bold"
)
ax1.set_xlim(0, 12.8)
ax1.grid(axis="x", linestyle="--", alpha=0.35, color="#888888")
ax1.set_axisbelow(True)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

for idx, row in sub_df.iterrows():
  sig_short = (
      "Sec/SPI"
      if "Sec/SPI" in str(row["signalp_prediction"])
      else (
          "Sec/SPII"
          if "Sec/SPII" in str(row["signalp_prediction"])
          else "Other"
      )
  )
  tag = (
      f"{row['heuristic_ranking_score']:.1f} | {sig_short} |"
      f" {row['deeptmhmm_type']} (TMR:{row['n_tmr']})"
  )
  ax1.text(
      row["heuristic_ranking_score"] + 0.15,
      idx,
      tag,
      va="center",
      ha="left",
      fontsize=7.2,
      color="#111111",
  )

ax1.set_title(
    "A. Prioritized Effectors & BGC Candidates (SignalP 6.0 + DeepTMHMM)",
    fontsize=10.5,
    fontweight="bold",
    loc="left",
    pad=10,
)

# --- PANEL B ---
ax2 = fig.add_subplot(gs[0, 1])
genome_len_mb = 3.704465
ax2.plot([0.52, 0.52], [0, genome_len_mb], color="#444444", linewidth=3.0, zorder=1)

# Sol taraftaki bolge etiketleri icin cakismasiz dikey yerlesim
left_items = []
for _, r in regions_df.iterrows():
  s_mb = r["start"] / 1e6
  h_mb = max((r["end"] - r["start"]) / 1e6, 0.022)
  mid_mb = s_mb + h_mb / 2
  is_key = r["region"] in {
      "Region 1",
      "Region 2",
      "Region 4",
      "Region 11",
      "Region 12",
  }
  left_items.append({
      "region": r["region"],
      "product": r["product"],
      "s_mb": s_mb,
      "h_mb": h_mb,
      "mid_mb": mid_mb,
      "is_key": is_key,
  })

left_items.sort(key=lambda x: x["mid_mb"])
min_gap_left = 0.135
last_y_left = -0.2
for item in left_items:
  disp_y = max(item["mid_mb"], last_y_left + min_gap_left)
  last_y_left = disp_y
  item["disp_y"] = disp_y

  box_col = "#d95f02" if item["is_key"] else "#bdbdbd"
  ax2.add_patch(
      mpatches.Rectangle(
          (0.41, item["s_mb"]),
          0.09,
          item["h_mb"],
          facecolor=box_col,
          edgecolor="#222222",
          linewidth=0.5,
          zorder=3,
      )
  )
  ax2.plot(
      [0.39, 0.41],
      [disp_y, item["mid_mb"]],
      color="#666666",
      linewidth=0.6,
      linestyle=":",
  )
  ax2.text(
      0.38,
      disp_y,
      f"{item['region']}: {item['product']}",
      ha="right",
      va="center",
      fontsize=7.2,
      fontweight="bold" if item["is_key"] else "normal",
      color="#111111" if item["is_key"] else "#555555",
  )

# Sag taraftaki anahtar aday isaretcileri
key_markers = sub_df[
    sub_df["wp_id"].isin({
        "WP_012008887.1",
        "WP_012009177.1",
        "WP_012009179.1",
        "WP_003214284.1",
        "WP_012009488.1",
        "WP_041815679.1",
        "WP_012010478.1",
        "WP_025093875.1",
        "WP_012011628.1",
    })
].sort_values("nt_start")

last_y_right = -0.2
min_gap_right = 0.17
for _, row in key_markers.iterrows():
  pos_mb = row["nt_start"] / 1e6
  disp_y = max(pos_mb, last_y_right + min_gap_right)
  last_y_right = disp_y
  c = palette[row["category"]]
  ax2.scatter(
      [0.54],
      [pos_mb],
      color=c,
      s=38,
      zorder=4,
      edgecolor="black",
      linewidth=0.5,
  )
  ax2.plot(
      [0.54, 0.62],
      [pos_mb, disp_y],
      color=c,
      linewidth=0.8,
      linestyle=":",
  )
  ax2.text(
      0.63,
      disp_y,
      f"{row['label']} ({pos_mb:.2f} Mb)",
      ha="left",
      va="center",
      fontsize=7.2,
      color="#111111",
  )

max_y_lim = max(genome_len_mb + 0.15, last_y_left + 0.15, last_y_right + 0.15)
ax2.set_xlim(-0.28, 1.38)
ax2.set_ylim(max_y_lim, -0.15)
ax2.set_xticks([])
ax2.set_ylabel(
    "Chromosome Position (NC_009848.4, Mbp)", fontsize=9.5, fontweight="bold"
)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.spines["bottom"].set_visible(False)
ax2.set_title(
    "B. Chromosomal Map of antiSMASH BGCs & Key Candidates",
    fontsize=10.5,
    fontweight="bold",
    loc="left",
    pad=10,
)

legend_handles = [
    mpatches.Patch(color=col, label=cat, edgecolor="#222222", linewidth=0.5)
    for cat, col in palette.items()
]
fig.legend(
    handles=legend_handles,
    loc="lower center",
    ncol=3,
    fontsize=8.0,
    frameon=True,
    facecolor="white",
    edgecolor="#cccccc",
    bbox_to_anchor=(0.5, -0.02),
)

fig.suptitle(
    "Bacillus pumilus (NC_009848.4) Multi-Layer Antagonistic Potential Against S. aureus NCTC 8325",
    fontsize=12.0,
    fontweight="bold",
    y=0.98,
)

out_png = latest_run / "Bacillus_pumilus_antagonism_summary_300dpi_fixed.png"
out_pdf = latest_run / "Bacillus_pumilus_antagonism_summary_300dpi_fixed.pdf"
fig.savefig(out_png, dpi=300, bbox_inches="tight")
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Fixed 300 DPI PNG saved: {out_png.resolve()}")
print(f"Fixed 300 DPI PDF saved: {out_pdf.resolve()}")