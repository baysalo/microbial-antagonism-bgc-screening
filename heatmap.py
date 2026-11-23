from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

runs_root = Path("runs")
latest_run = sorted(
    [p for p in runs_root.iterdir() if p.is_dir()],
    key=lambda p: p.name,
    reverse=True,
)[0]

strain_df = pd.read_csv(latest_run / "strain_summary_per_1k_cds.csv")

rate_cols = [
    "sec_sig_per_1k_cds",
    "export_per_1k_cds",
    "protease_per_1k_cds",
    "antimicrobial_per_1k_cds",
    "lipase_per_1k_cds",
    "lytic_per_1k_cds",
    "high_rank_per_1k_cds",
]

mat = strain_df.set_index("strain")[rate_cols].T

fig, ax = plt.subplots(figsize=(7.5, 5.2))
im = ax.imshow(mat.values, aspect="auto", cmap="YlGnBu")

ax.set_xticks(np.arange(mat.shape[1]))
ax.set_xticklabels(mat.columns, rotation=15, ha="right")
ax.set_yticks(np.arange(mat.shape[0]))
ax.set_yticklabels(mat.index)

for i in range(mat.shape[0]):
  for j in range(mat.shape[1]):
    val = mat.values[i, j]
    ax.text(
        j,
        i,
        f"{val:.2f}",
        ha="center",
        va="center",
        color="black",
        fontsize=9,
    )

cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Rate per 1,000 CDS (Heuristic)")
ax.set_title("Heuristic Feature Rates per 1,000 CDS (Not a Biological Assay)")
fig.tight_layout()

out_png = latest_run / "strain_feature_heatmap_per_1k_cds.png"
fig.savefig(out_png, dpi=300)
plt.close(fig)
print(f"Saved heatmap to: {out_png.resolve()}")