import json
from pathlib import Path

json_path = Path("BacilluspGCF_000017885.4_ASM1788v4_genomic.json")  # kendi JSON dosya adinizi yazin
data = json.loads(json_path.read_text(encoding="utf-8"))

for rec in data.get("records", []):
  rec_id = rec.get("id", "")
  for idx, reg in enumerate(rec.get("areas", rec.get("regions", [])), 1):
    print(
        f"{rec_id} | Region {idx} | {reg.get('start')}..{reg.get('end')} |"
        f" {reg.get('products')}"
    )