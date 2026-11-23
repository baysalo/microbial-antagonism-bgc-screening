from pathlib import Path
import pandas as pd

runs_root = Path("runs")
run_dirs = sorted(
    [p for p in runs_root.iterdir() if p.is_dir()],
    key=lambda p: p.name,
    reverse=True,
)
if not run_dirs:
  raise FileNotFoundError(f"No run directories found in {runs_root.resolve()}")

latest_run = run_dirs[0]
print(f"Latest run directory: {latest_run.resolve()}\n")

genome_df = pd.read_csv(latest_run / "genome_summary.csv")
strain_df = pd.read_csv(latest_run / "strain_summary_per_1k_cds.csv")

print("--- genome_summary.csv ---")
print(genome_df.to_string(index=False))

print("\n--- strain_summary_per_1k_cds.csv ---")
print(strain_df.to_string(index=False))