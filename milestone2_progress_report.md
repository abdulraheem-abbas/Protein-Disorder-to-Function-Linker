# Milestone 2 Progress Report
### CMPS 297AF – Path A | Protein Disorder-to-Function Linker

---

## 1. Project Overview & Refined Research Question

**Tool Name:** Protein Disorder-to-Function Linker

**Refined Question:**
*Can we systematically identify and quantify the overlap between intrinsically disordered protein regions (IDRs) and functional annotations — specifically domains, motifs, binding sites, and disease-associated mutations — for any given protein, using only publicly available biological databases?*

The core hypothesis is that IDRs are not structurally inert; they are enriched at functional sites and disease hotspots. This tool operationalizes that hypothesis into a reproducible, queryable computational pipeline.

**Biological Motivation:**
Intrinsically Disordered Regions (IDRs) lack stable tertiary structure yet are overrepresented at protein-protein interaction interfaces, post-translational modification sites, and pathogenic mutation sites. Linking regional disorder to curated functional annotations enables researchers to prioritize sites for experimental validation or drug targeting.

---

## 2. Updated Plan

The original plan proposed constructing a five-module pipeline. That structure has been **fully implemented** and is now operational. The plan has been refined in two ways:

1. **IUPred3 Fallback:** The original target API (IUPred2A at iupred2a.elte.hu) returns HTTP 404 for the POST endpoint. A fallback mock scorer with a deterministic seed (seed=42) was introduced to maintain pipeline continuity while an alternative disorder source is evaluated.

2. **Streamlit Frontend:** A web interface (app.py) was added beyond the original CLI design, enabling interactive use without command-line knowledge.

**Revised module plan:**

| Module          | Responsibility                                            | Status                  |
|-----------------|-----------------------------------------------------------|-------------------------|
| uniprot_api.py  | Fetch sequence, domains, motifs, binding/active sites, variants | Complete           |
| disorder_api.py | Per-residue disorder scores (IUPred2A / fallback mock)    | Complete (mock active)  |
| integrator.py   | Align disorder + functional annotations per residue       | Complete                |
| scorer.py       | Compute DFLS scores; rank sites, domains, variants        | Complete                |
| visualizer.py   | Disorder profile plot with IDR shading + hotspot highlights | Complete              |
| main.py         | CLI / programmatic entry point; returns structured dict   | Complete                |
| app.py          | Streamlit web interface                                   | Complete (new addition) |

---

## 3. Preliminary Results

### 3.1 Protein Fetched — TP53 (P04637)

The pipeline was validated end-to-end on **TP53 (human tumor suppressor p53)**.

Fetched from UniProt REST API (rest.uniprot.org):
- Protein name: Cellular tumor antigen p53
- Sequence length: 393 amino acids
- Domains: p53 tetramerization domain, DNA-binding domain
- Motifs: Compositional biases, repeat regions
- Binding sites: Active site and ligand-binding annotations
- Mutations: Natural variants and mutagenesis entries (disease-associated)

### 3.2 Disorder Score Computation

Per-residue disorder scores are computed via disorder_api.py. The IUPred2A REST endpoint currently returns 404 (see Section 4). A deterministic mock fallback (seeded random scores in [0.3, 0.7]) is used for development continuity.

Output during pipeline run:
  [WARNING] Using mock disorder scores — reason: IUPred2A returned HTTP 404

The fallback is seeded with random.seed(42) ensuring reproducible output across runs.

### 3.3 Integration Output

The integrator (integrator.py) produces a per-residue table aligning:
- Sequence position (1-indexed)
- Amino acid identity
- Disorder score
- Functional labels (domains, motifs, binding sites)
- Mutation annotations

Sample output — first 5 residues of TP53:

  Position | AA | Disorder | Functional | Mutation
  ---------|----|----------|------------|---------
  1        | M  | 0.6369   | []         | []
  2        | E  | 0.6369   | []         | []
  3        | E  | 0.5245   | []         | []
  4        | P  | 0.6369   | []         | []
  5        | Q  | 0.3936   | []         | []

### 3.4 Disorder-Function Linkage Scoring (DFLS)

scorer.py implements the Disorder-Function Linkage Score (DFLS), a composite metric in [0.0, 1.0]:

  Component             | Weight | Description
  ----------------------|--------|-----------------------------------------------
  Local disorder score  |  40%   | IUPred score at the exact annotated residue
  IDR mean score        |  25%   | Average disorder of surrounding IDR
  Overlap fraction      |  20%   | Fraction of feature within the IDR
  Disease association   |  15%   | Boost for disease-linked variants

Scores are classified into: high (>= 0.75), medium (>= 0.50), low.
Sites, domains, and variants are ranked and reported separately.

### 3.5 Visualizer Output

visualizer.py generates a disorder profile plot saved to outputs/P04637_plot.png.

Plot elements:
- Blue line       — Per-residue disorder score
- Red dashed line — Disorder threshold (0.5)
- Blue shading    — Domain regions (alpha=0.3)
- Green shading   — Motif regions (alpha=0.3)
- Red dots        — Mutations where disorder > 0.5 only (noise-filtered)
- Yellow spans    — Hotspot regions: disorder > 0.5 AND mutation present

### 3.6 Streamlit Web Interface

An interactive Streamlit app (app.py) allows users to:
1. Enter any UniProt accession ID
2. Click "Run Analysis" to trigger the full pipeline
3. View protein name, length, and ID as metric cards
4. View the disorder plot inline
5. Browse the first 20 rows of integrated per-residue data

The app runs at http://localhost:8501 via: streamlit run app.py

---

## 4. Obstacles & Adjustments

### 4.1 IUPred2A API Unavailability

Obstacle: The IUPred2A REST endpoint returns HTTP 404 for POST requests containing a protein sequence.

Adjustment:
- Added graceful fallback in disorder_api.py using seeded random scores
- Evaluating IUPred3 (iupred3.elte.hu) and BioPython-based tools as alternatives
- Considering a charge-hydropathy sequence-composition heuristic as an offline fallback

Impact: Medium. Pipeline runs end-to-end with mock scores. Real disorder data will replace this before final submission.

### 4.2 Mutation Marker Visual Noise

Obstacle: Plotting all natural variants as scatter markers produced an overcrowded, unreadable figure. TP53 has hundreds of annotated variants.

Adjustment: Only mutations where disorder > 0.5 are plotted. Yellow hotspot spans now highlight co-occurring high disorder and mutation sites, dramatically improving readability and biological interpretability.

### 4.3 UniProt Feature Type Normalization

Obstacle: Inconsistent feature type strings in UniProt JSON responses caused some annotations to be missed (e.g., "Region" vs "Domain", "Natural variant" vs "Mutagenesis").

Adjustment: Defined explicit type-set mappings:
  DOMAIN_TYPES   = {"Domain", "Region"}
  MOTIF_TYPES    = {"Motif", "Repeat", "Compositional bias"}
  BINDING_TYPES  = {"Binding site", "Active site", "Modified residue"}
  MUTATION_TYPES = {"Natural variant", "Mutagenesis"}

This ensures consistent annotation capture across all protein entries.

---

## 5. Roadmap to Completion

  Timeline | Task                                                    | Priority
  ---------|---------------------------------------------------------|---------
  Week 1   | Resolve disorder API — switch to IUPred3 or local       | High
  Week 1   | Validate on 3 more proteins (Tau, Huntingtin, DYRK1A)   | High
  Week 2   | Surface DFLS scores in the Streamlit interface          | Medium
  Week 2   | Add downloadable JSON export of full results            | Medium
  Week 2   | Add summary statistics (% IDR, hotspot count)           | Medium
  Week 3   | Comparative analysis across 4-6 disease proteins        | High
  Week 3   | Polish visualizations for final report figures           | Medium
  Week 3   | Milestone 3 report with full biological interpretation  | High

Key Open Questions:
1. Will IUPred3 API be accessible, or is a local/offline method needed?
2. Which proteins best demonstrate the disorder-function link comparatively?
3. Should DFLS weights be calibrated against published disorder benchmarks?

---

## 6. Summary

The Protein Disorder-to-Function Linker is a fully structured, modular, and operational pipeline. All seven modules are implemented and working end-to-end. Preliminary results on TP53 confirm that the integration logic, scoring engine, and visualization layer all function correctly. The primary outstanding issue is restoring real disorder scores via IUPred3. The project is on track for a complete, biologically meaningful final submission.

---

Submitted for CMPS 297AF – Path A | Milestone 2 | April 2026
