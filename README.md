# 🧬 Protein Disorder-to-Function Linker

![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python&logoColor=white)
![Streamlit](https://img.shields.io/badge/Streamlit-app-FF4B4B?logo=streamlit&logoColor=white)
![License](https://img.shields.io/badge/License-MIT-green)
![Course](https://img.shields.io/badge/CMPS%20297AF-Path%20A-purple)

A modular Python bioinformatics pipeline that integrates **per-residue intrinsic disorder predictions** with **functional annotations from UniProt** — including domains, binding sites, PTM motifs, and disease-associated variants. Results are scored, visualized, and accessible via both a CLI and an interactive **Streamlit** web app.

---

## 🧫 Biological Motivation

Many proteins contain regions that lack a fixed 3D structure — **Intrinsically Disordered Regions (IDRs)**. Far from being non-functional "junk," IDRs are enriched at:

- Protein–protein interaction interfaces
- Post-translational modification (PTM) sites
- Disease-related mutation hotspots

This tool systematically identifies **where disorder and function intersect**, quantifying the relationship through a per-feature **Disorder-Function Linkage Score (DFLS)**.

---

## 📁 Project Structure

```
protein_disorder_linker/
├── app.py            – Streamlit web interface
├── main.py           – CLI entry point; orchestrates the pipeline
├── uniprot_api.py    – Fetches sequence, domains, motifs, binding sites, variants
├── disorder_api.py   – IUPred3 per-residue disorder predictions
├── integrator.py     – Maps IDRs to functional annotations
├── scorer.py         – Computes Disorder-Function Linkage Score (DFLS)
├── visualizer.py     – Disorder profile + annotated linear protein map
├── requirements.txt
├── data/             – (Optional) cached API responses
└── outputs/          – Generated plots and JSON results
```

---

## ⚙️ Installation

```bash
git clone https://github.com/abdulraheem-abbas/Protein-Disorder-to-Function-Linker1.git
cd Protein-Disorder-to-Function-Linker1
pip install -r requirements.txt
```

---

## 🚀 Usage

### Streamlit Web App (recommended)

```bash
streamlit run app.py
```

Enter any UniProt accession ID (e.g., `P04637` for TP53) and click **Analyze** to get an interactive view of disorder profiles, functional maps, and scored annotations.

### CLI

```bash
# Run on TP53 (canonical human tumor suppressor)
python3 main.py --uniprot P04637

# Positional input also works
python3 main.py P04637

# Custom disorder threshold + save JSON output
python3 main.py --uniprot P04637 --disorder-threshold 0.6 --save-json

# Save plots to a specific directory
python3 main.py --uniprot P04637 --output results/TP53/

# Skip visualization (useful for batch runs)
python3 main.py --uniprot P04637 --no-plot

# Skip AlphaFold pLDDT overlay/fetch
python3 main.py --uniprot P04637 --no-plddt
```

---

## 🔬 Pipeline Overview

```
UniProt Accession
      │
      ▼
[1] uniprot_api.py   → sequence, domains, motifs, binding sites, variants
      │
      ▼
[2] disorder_api.py  → per-residue IUPred3 disorder scores
      │
      ▼
[3] integrator.py    → IDR intervals + overlap with functional features
      │
      ▼
[4] scorer.py        → IDR summary, variant-dense windows, enrichment statistics
      │
      ▼
[5] visualizer.py    → disorder profile + annotated linear protein map
```

---

## 📊 Output

| File | Description |
|------|-------------|
| `{ACC}_plot.png` | Per-residue disorder score with IDR shading and site markers |
| `{ACC}_linear_map.png` | 4-track linear map: IDRs, domains/motifs, variants, PTMs |
| `{ACC}_results.json` | Full scored results (with `--save-json`) |
| `{ACC}_residues.csv` | Per-residue integrated table (with `--save-csv`) |

---

## 🏅 Scoring: Disorder-Function Linkage Score (DFLS)

| Component | Weight | Description |
|-----------|--------|-------------|
| Local disorder score | 40% | IUPred3 score at the exact site residue |
| IDR mean score | 25% | Average disorder of the surrounding region |
| Overlap fraction | 20% | Fraction of the feature within the IDR |
| Disease association | 15% | Boost for variants linked to disease |

Scores are classified as **high** (≥ 0.75), **medium** (≥ 0.50), or **low**.

The current pipeline reports residue-level IDR summaries, PTM/variant overlap,
variant-dense UniProt annotation windows, and Fisher's exact enrichment tests.
Legacy DFLS helper methods remain in `scorer.py`, but the active CLI/app output
is intentionally focused on transparent count-based summaries and enrichment.

If IUPred3 is unavailable, accession-based analyses fall back to deterministic
mock disorder scores and the result metadata/UI warning marks that clearly.
Raw sequence mode also uses mock disorder scores because the public IUPred3
endpoint used here requires a UniProt accession.

---

## 🧪 Example Proteins to Try

| UniProt ID | Protein | Why interesting |
|------------|---------|-----------------|
| `P04637` | TP53 (human) | Disordered N/C termini, many disease variants |
| `P42858` | Huntingtin | Polyglutamine disorder, neurodegeneration |
| `P10636` | Tau | Fully disordered IDP, Alzheimer's disease |
| `Q9NR61` | DYRK1A | Kinase with IDR-mediated substrate recognition |

---

## 📚 References

- [IUPred3](https://iupred3.elte.hu/) — Intrinsic disorder prediction
- [UniProt REST API](https://www.uniprot.org/help/api) — Protein functional annotations
- [Streamlit](https://streamlit.io/) — Interactive web app framework
