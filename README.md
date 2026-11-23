# Microbial Antagonism & BGC Screening Pipeline



Genomic and functional candidate prioritization workflow for *Bacillus pumilus* (`NC_009848.4`) against *Staphylococcus aureus* NCTC 8325, integrating heuristic candidate ranking, **SignalP 6.0** secretion prediction, **DeepTMHMM** topology filtering, **antiSMASH** biosynthetic gene cluster (BGC) mapping, and publication-grade (**300 DPI**) summary visualization.



---



## Key Highlights



* **Multi-Layer Screening:** Combines secreted/surface effector candidates (`SignalP 6.0` + `DeepTMHMM`) with secondary metabolite and bacteriocin BGCs (`antiSMASH` Regions 1–12).

* **Prioritized Candidates:** Includes soluble/secreted effectors such as antimicrobial peptide **LCI** (`WP_012010478.1`), **$\beta$-glucanase (`bglS`, `WP_003214284.1`)**, **D-Ala-D-Ala carboxypeptidase (`WP_041815679.1`)**, and **S8/C40/metallopeptidases**, alongside BGC-associated loci in **Region 1** (`NRPS`), **Region 2** (`NRPS, T1PKS`), **Region 4** (`terpene, NI-siderophore` near `WP_012009488.1`), **Region 11** (`RiPP-like` circular bacteriocin `WP_025093875.1`), and **Region 12** (`bacilysin` / `BacA`, `WP_012011628.1`).

* **300 DPI Reporting:** Generates multi-panel summary figures (`PNG` and `PDF`) and pushes reproducible scripts to `https://github.com/baysalo/microbial-antagonism-bgc-screening`.



---



## Repository & Workflow Structure



```text

microbial-antagonism-bgc-screening/

├── README.md

├── .gitignore

├── push.ps1

├── plot_summary_300dpi_fixed.py

└── runs/

    └── <latest_run>/

        ├── Bacillus_pumilus_ranked_candidates.csv

        ├── signalp6_deeptmhmm_merged.csv

        ├── Bacillus_pumilus_antagonism_summary_300dpi_fixed.png

        └── Bacillus_pumilus_antagonism_summary_300dpi_fixed.pdf

```



---



## Quick Usage



### 1. Generate 300 DPI Summary Figure



```bash

python plot_summary_300dpi_fixed.py

```



### 2. Push to GitHub (`push.ps1`)



```powershell

$env:PATH = "C:\Program Files\Git\cmd;$env:PATH"

git fetch origin main

git merge origin/main --allow-unrelated-histories --no-edit

git push -u origin main

```



---



## Next Validation Priorities



1. **In Silico & Structural Analysis:** Evaluate mature peptide/domain features (e.g., **LCI**, **`bglS`**, **Region 11/12**) and cluster similarity (`KnownClusterBlast` / `MIBiG`).

2. **Experimental Validation:** Test *B. pumilus* cell-free supernatants, fractions, and co-culture expression (`RT-qPCR`) against *S. aureus* NCTC 8325.
