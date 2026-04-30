"""Streamlit interface for the Protein Disorder-to-Function Linker."""

import streamlit as st
import pandas as pd
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

from main    import run_pipeline
from compare import compare_proteins, build_comparison_plotly
from visualizer import build_plotly_figure

st.set_page_config(
    page_title="Protein Disorder-to-Function Linker",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');

    html, body, [class*="css"] { font-family: 'Inter', sans-serif; }

    .stApp { background: #0f1117; }

    section[data-testid="stSidebar"] {
        background: linear-gradient(180deg, #1a1d27 0%, #141720 100%);
        border-right: 1px solid #2d3148;
    }

    .metric-card {
        background: linear-gradient(135deg, #1e2235 0%, #252a3d 100%);
        border: 1px solid #2d3148;
        border-radius: 12px;
        padding: 16px 20px;
        text-align: center;
    }
    .metric-card .label {
        color: #8892b0;
        font-size: 11px;
        font-weight: 600;
        letter-spacing: 0.08em;
        text-transform: uppercase;
        margin-bottom: 6px;
    }
    .metric-card .value {
        color: #e2e8f0;
        font-size: 22px;
        font-weight: 700;
    }
    .metric-card .unit {
        color: #64ffda;
        font-size: 12px;
        font-weight: 500;
    }

    .section-header {
        color: #64ffda;
        font-size: 13px;
        font-weight: 600;
        letter-spacing: 0.1em;
        text-transform: uppercase;
        margin-bottom: 12px;
        padding-bottom: 6px;
        border-bottom: 1px solid #2d3148;
    }

    .cluster-badge {
        display: inline-block;
        background: rgba(255,140,0,0.2);
        border: 1px solid #ff8c00;
        color: #ff8c00;
        border-radius: 6px;
        padding: 3px 10px;
        font-size: 12px;
        font-weight: 600;
        margin: 3px;
    }

    .idr-badge {
        display: inline-block;
        background: rgba(100,255,218,0.1);
        border: 1px solid #64ffda;
        color: #64ffda;
        border-radius: 6px;
        padding: 3px 10px;
        font-size: 12px;
        margin: 3px;
    }

    div[data-testid="stTabs"] button {
        font-size: 14px;
        font-weight: 500;
    }
</style>
""", unsafe_allow_html=True)


with st.sidebar:
    st.markdown("## 🧬 Disorder Linker")
    st.markdown("---")

    st.markdown("### ⚙️ Analysis Settings")
    threshold  = st.slider("IDR disorder threshold", 0.0, 1.0, 0.5, 0.05,
                           help="Residues above this score are classified as disordered (IDR).")
    fetch_af   = st.checkbox("Fetch AlphaFold pLDDT", value=True,
                             help="Overlay AlphaFold structural confidence as a second disorder proxy.")
    show_csv   = st.checkbox("Show full residue table", value=False)

    st.markdown("---")
    st.markdown("### 🔬 Example proteins")
    examples = {
        "TP53 (P04637)"       : "P04637",
        "α-Synuclein (P37840)": "P37840",
        "p27 Kip1 (P46527)"   : "P46527",
        "Tau (P10636)"         : "P10636",
        "Huntingtin (P42858)"  : "P42858",
    }
    example_choice = st.selectbox("Quick load", ["— select —"] + list(examples.keys()))

    st.markdown("---")
    st.caption("CMPS 297AF · Path A · Bioinformatics Tool")


tab1, tab2, tab3 = st.tabs(["🔬 Single Protein", "📊 Compare Proteins", "📖 About"])

with tab1:
    st.markdown("### Single Protein Disorder Analysis")
    st.markdown("Enter a UniProt accession ID **or** paste a raw amino-acid sequence.")

    default_input = examples[example_choice] if example_choice != "— select —" else "P04637"
    col_input, col_btn = st.columns([4, 1])
    with col_input:
        accession = st.text_input("UniProt ID or sequence", value=default_input,
                                  placeholder="e.g. P04637", label_visibility="collapsed")
    with col_btn:
        run_btn = st.button("Analyze", type="primary", width="stretch")

    if run_btn:
        if not accession.strip():
            st.error("Please enter a UniProt accession ID or sequence.")
        else:
            with st.spinner(f"Running pipeline for **{accession.strip()[:10]}**…"):
                result = run_pipeline(
                    accession.strip(),
                    verbose=False,
                    fetch_plddt=fetch_af,
                    save_json=False,
                    save_csv=False,
                    disorder_threshold=threshold,
                )

            if result is None:
                st.error("Pipeline failed. Please check the accession ID and try again.")
            else:
                scored   = result["scored"]
                summary  = scored["summary"]
                residues = result["integrated_data"]
                clusters = scored.get("clusters", [])
                disorder_metadata = result.get("disorder_metadata", {})

                st.success(f"✅ Analysis complete for **{result['name']}**")
                if disorder_metadata.get("is_mock"):
                    st.warning(
                        "Disorder scores are mock fallback values, not live IUPred3 output. "
                        f"Reason: {disorder_metadata.get('reason', 'Unavailable')}"
                    )

                if result.get("gene_name") or result.get("organism"):
                    st.caption(
                        f"Gene: **{result.get('gene_name','—')}** · "
                        f"Organism: *{result.get('organism','—')}*"
                    )

                m1, m2, m3, m4, m5 = st.columns(5)
                metrics = [
                    (m1, "Length",                             f"{result['length']:,}", "aa"),
                    (m2, f"% IDR (\u2265{threshold:g})",       f"{summary.get('pct_disordered',0):.1f}", "%"),
                    (m3, "UniProt variant positions in IDRs",  str(summary.get("n_mutations_in_idr", 0)), ""),
                    (m4, "PTM Sites in IDRs",                  str(summary.get("n_ptm_in_idr", 0)), ""),
                    (m5, "Variant-Dense Windows",              str(summary.get("n_variant_dense_windows", 0)), ""),
                ]
                for col, label, value, unit in metrics:
                    col.markdown(
                        f'<div class="metric-card">'
                        f'<div class="label">{label}</div>'
                        f'<div class="value">{value}<span class="unit"> {unit}</span></div>'
                        f'</div>', unsafe_allow_html=True
                    )

                st.markdown("")

                with st.expander("❓ What do these numbers mean?", expanded=False):
                    st.markdown(f"""
| Metric | Meaning |
|--------|---------|
| **% IDR (≥{threshold:g})** | Percentage of residues predicted as intrinsically disordered by IUPred3 at the current threshold. Higher = more disordered protein. |
| **UniProt variant/mutagenesis-annotated IDR residues** | Number of residue positions inside IDRs that carry a UniProt *Natural variant* or *Mutagenesis* annotation. **Source: UniProt only — not COSMIC or gnomAD.** |
| **PTM Sites in IDRs** | Phosphorylation, glycosylation, ubiquitination, and other post-translational modification sites (UniProt) that fall within disordered regions. IDR-enriched PTMs are often regulatory switches. |
| **Variant-Dense UniProt Annotation Windows** | Contiguous local 15-residue windows with ≥3 UniProt variant or mutagenesis-annotated positions. This reflects UniProt annotation density only and does not represent somatic cancer mutation recurrence. Treat as hypothesis-generating regions, not validated mutation hotspots. |

**Disorder score interpretation (IUPred3):**
- Score **≥ {threshold:g}** → predicted Intrinsically Disordered Region (IDR)
- Score **< {threshold:g}** → predicted ordered / structured region

**AlphaFold pLDDT (if shown):**
- Displayed as a normalized disorder proxy: `(100 − pLDDT) / 100`
- pLDDT < 50 → very low structural confidence (likely disordered)
- pLDDT > 90 → high structural confidence (well-folded)
- Concordance between IUPred3 and AlphaFold strengthens IDR calls.
                    """)

                st.markdown('<div class="section-header">Disorder Profile (Interactive)</div>',
                            unsafe_allow_html=True)
                fig = build_plotly_figure(
                    residues,
                    protein_name=result["name"],
                    threshold=threshold,
                    clusters=clusters,
                )
                st.plotly_chart(fig, width="stretch")

                if os.path.exists(result.get("linear_map_file", "")):
                    st.markdown('<div class="section-header">Linear Protein Map</div>',
                                unsafe_allow_html=True)
                    st.image(result["linear_map_file"], width="stretch")

                total_annot    = sum(1 for r in residues if r.get('mutation'))
                idr_annot      = summary.get('n_mutations_in_idr', 0)
                if clusters:
                    st.markdown('<div class="section-header">Variant-Dense UniProt Annotation Windows (sliding window: 15 residues, ≥3 annotated positions, top-5 non-overlapping)</div>',
                                unsafe_allow_html=True)
                    st.markdown(
                        f"**Total UniProt variant-annotated residues:** {total_annot} &nbsp;|&nbsp; "
                        f"**IDR-associated variant-annotated residues:** {idr_annot}",
                        unsafe_allow_html=True
                    )
                    cl_cols = st.columns(min(len(clusters), 4))
                    for i, cl in enumerate(clusters):
                        col = cl_cols[i % len(cl_cols)]
                        idr_tag = "✅ IDR-overlapping" if cl.get("is_idr_window", cl.get("is_idr_cluster")) else "🔵 Structured region"
                        n_annot = cl.get('n_annotations', cl.get('n_mutations', 0))
                        density = cl.get('annotation_density', '')
                        col.markdown(
                            f'<span class="cluster-badge">Residues {cl["start"]}–{cl["end"]}</span>'
                            f'<br><small>{n_annot} UniProt variant/mutagenesis annotations · {idr_tag}'
                            + (f' · density={density}' if density else '') + '</small>',
                            unsafe_allow_html=True
                        )
                    st.caption(
                        "⚠️ These windows reflect **UniProt annotation density only** — not cancer recurrence frequency. "
                        "They should be interpreted as hypothesis-generating regions, not validated mutation hotspots."
                    )
                    st.markdown("")

                ptm_in_idr = scored.get("ptm_in_idr", [])
                if ptm_in_idr:
                    st.markdown('<div class="section-header">PTM Sites in IDRs (UniProt Modified residue, Glycosylation, Lipidation)</div>',
                                unsafe_allow_html=True)
                    ptm_df = pd.DataFrame([
                        {
                            "Position": r["position"],
                            "AA":       r["aa"],
                            "Disorder": f"{r['disorder']:.3f}",
                            "PTM":      "; ".join(r.get("ptm", [])),
                        }
                        for r in ptm_in_idr[:25]
                    ])
                    st.dataframe(ptm_df, width="stretch", hide_index=True)

                n_conc = summary.get("n_concordant_idr", 0)
                n_dis  = summary.get("n_disordered", 0)
                if n_conc > 0:
                    pct = round(100 * n_conc / n_dis, 1) if n_dis else 0
                    st.info(
                        f"🔵 **AlphaFold concordance:** {n_conc} / {n_dis} disordered residues "
                        f"({pct}%) are also classified as low-confidence by AlphaFold pLDDT — "
                        f"high-confidence IDR calls."
                    )

                if summary.get("top_idr_mutation_region"):
                    s, e = summary["top_idr_mutation_region"]
                    st.markdown(
                        f'<span class="idr-badge">📍 Largest IDR-associated UniProt annotation-dense region: residues {s}–{e}</span>',
                        unsafe_allow_html=True
                    )
                    st.markdown("")

                st.caption(
                    "⚠️ **Mutation data source:** UniProt REST API — "
                    "*Natural variant* (germline SNPs, ClinVar-linked) and *Mutagenesis* "
                    "(experimentally validated) entries. "
                    "**Somatic mutations (COSMIC) and population variants (gnomAD) are NOT included.**"
                )

                enrichment = result.get("enrichment", {})
                if enrichment:
                    st.markdown('<div class="section-header">&#x1F4CA; Statistical Validation &#8212; Fisher Exact Test</div>',
                                unsafe_allow_html=True)
                    ve = enrichment.get("variant_enrichment", {})
                    pe = enrichment.get("ptm_enrichment", {})
                    ec_a, ec_b = st.columns(2)
                    with ec_a:
                        st.markdown("**UniProt Variant/Mutagenesis Annotation Enrichment in IDRs**")
                        if ve.get("p_value") is not None:
                            color = "green" if ve["significant"] else "gray"
                            _p_ve = ve['p_value']
                            _p_ve_str = 'p < 0.001' if _p_ve < 0.001 else f'p = {_p_ve:.5f}'
                            st.markdown(
                                f"OR: **{ve['odds_ratio']}** &nbsp;|&nbsp; **{_p_ve_str}** &nbsp;"
                                f"<span style='color:{color}'>{ve['significance']}</span>",
                                unsafe_allow_html=True)
                            st.caption(ve["interpretation"])
                            t = ve["table"]
                            st.dataframe(pd.DataFrame({
                                "": ["Has UniProt variant/mutagenesis annotation", "No annotation"],
                                "IDR residues": [t["IDR_with"], t["IDR_without"]],
                                "Ordered residues": [t["ordered_with"], t["ordered_without"]],
                            }), hide_index=True, use_container_width=True)
                        else:
                            st.info("No variant annotations — test not applicable.")
                    with ec_b:
                        st.markdown("**PTM Site Enrichment in IDRs**")
                        if pe.get("p_value") is not None:
                            color = "green" if pe["significant"] else "gray"
                            _p_pe = pe['p_value']
                            _p_pe_str = 'p < 0.001' if _p_pe < 0.001 else f'p = {_p_pe:.5f}'
                            st.markdown(
                                f"OR: **{pe['odds_ratio']}** &nbsp;|&nbsp; **{_p_pe_str}** &nbsp;"
                                f"<span style='color:{color}'>{pe['significance']}</span>",
                                unsafe_allow_html=True)
                            st.caption(pe["interpretation"])
                            t = pe["table"]
                            st.dataframe(pd.DataFrame({
                                "": ["Has PTM", "No PTM"],
                                "IDR": [t["IDR_with"], t["IDR_without"]],
                                "Ordered": [t["ordered_with"], t["ordered_without"]],
                            }), hide_index=True, use_container_width=True)
                        else:
                            st.info("No PTM annotations — test not applicable.")
                    st.caption(f"⚠️ **DFLS disclaimer:** {enrichment.get('dfls_disclaimer','')}")
                    st.markdown("")

                st.markdown('<div class="section-header">Export Results</div>',
                            unsafe_allow_html=True)
                ec1, ec2, ec3 = st.columns(3)

                import json
                json_data = json.dumps({
                    "accession":   result["accession"],
                    "name":        result["name"],
                    "gene_name":   result.get("gene_name", ""),
                    "organism":    result.get("organism", ""),
                    "length":      result["length"],
                    "disorder_threshold": result.get("disorder_threshold", threshold),
                    "disorder_metadata": result.get("disorder_metadata", {}),
                    "summary":     scored["summary"],
                    "clusters":    scored.get("clusters", []),
                }, indent=2)
                ec1.download_button("⬇ Download JSON",
                                    data=json_data,
                                    file_name=f"{result['accession']}_results.json",
                                    mime="application/json")

                csv_rows = []
                for r in residues:
                    csv_rows.append({
                        "position":   r["position"],
                        "aa":         r["aa"],
                        "disorder":   r["disorder"],
                        "plddt":      r.get("plddt", ""),
                        "functional": "; ".join(r["functional"]),
                        "ptm":        "; ".join(r.get("ptm", [])),
                        "mutation":   "; ".join(r["mutation"]),
                    })
                csv_df  = pd.DataFrame(csv_rows)
                ec2.download_button("⬇ Download CSV",
                                    data=csv_df.to_csv(index=False),
                                    file_name=f"{result['accession']}_residues.csv",
                                    mime="text/csv")

                if os.path.exists(result["plot_file"]):
                    with open(result["plot_file"], "rb") as f:
                        ec3.download_button("⬇ Download Plot",
                                            data=f,
                                            file_name=os.path.basename(result["plot_file"]),
                                            mime="image/png")

                if show_csv:
                    st.markdown('<div class="section-header">Per-Residue Data</div>',
                                unsafe_allow_html=True)
                    st.dataframe(csv_df, width="stretch", hide_index=True)


with tab2:
    st.markdown("### Multi-Protein Comparative Analysis")
    st.markdown(
        "Enter comma-separated UniProt accessions to compare disorder profiles, "
        "compare disorder profiles, PTM enrichment, and variant-dense UniProt annotation windows across proteins."
    )

    default_compare = "P04637, P37840, P46527"
    cmp_input = st.text_input("Accessions (comma-separated)",
                              value=default_compare,
                              placeholder="P04637, P37840, P46527")
    cmp_btn = st.button("Compare", type="primary")

    if cmp_btn:
        raw_acc = [a.strip().upper() for a in cmp_input.split(",") if a.strip()]
        if len(raw_acc) < 2:
            st.error("Please enter at least 2 accession IDs.")
        elif len(raw_acc) > 8:
            st.warning("Comparing more than 8 proteins may be slow. Consider reducing the list.")
        else:
            with st.spinner(f"Analyzing {len(raw_acc)} proteins…"):
                rows = compare_proteins(
                    raw_acc,
                    verbose=False,
                    fetch_plddt=fetch_af,
                    disorder_threshold=threshold,
                )

            if not rows:
                st.error("Comparison failed. Check the accession IDs.")
            else:
                st.success(f"✅ Compared {len(rows)} proteins successfully.")

                st.markdown('<div class="section-header">Comparison Chart</div>',
                            unsafe_allow_html=True)
                cmp_fig = build_comparison_plotly(rows)
                st.plotly_chart(cmp_fig, width="stretch")

                st.markdown('<div class="section-header">Ranked Summary Table</div>',
                            unsafe_allow_html=True)
                table_df = pd.DataFrame([
                    {
                        "Accession":             r["accession"],
                        "Protein":               r["name"][:40],
                        "Gene":                  r.get("gene_name", "—"),
                        "Length (aa)":           r["length"],
                        "% Disordered":          f"{r['pct_disordered']:.1f}%",
                        "UniProt variants/IDR":  r["n_mutations_in_idr"],
                        "PTM/IDR":               r["n_ptm_in_idr"],
                        "Annotation Windows":    r.get("n_variant_dense_windows", 0),
                    }
                    for r in rows
                ])
                st.dataframe(table_df, width="stretch", hide_index=True)

                st.markdown('<div class="section-header">Individual Disorder Profiles</div>',
                            unsafe_allow_html=True)
                for r in rows:
                    with st.expander(f"{r.get('gene_name') or r['accession']} — {r['name'][:50]}"):
                        if os.path.exists(r.get("plot_file", "")):
                            st.image(r["plot_file"], width="stretch")
                        if os.path.exists(r.get("linear_map_file", "")):
                            st.image(r["linear_map_file"], width="stretch")

                cmp_df = pd.DataFrame(rows)
                st.download_button(
                    "⬇ Download Comparison CSV",
                    data=cmp_df.to_csv(index=False),
                    file_name="comparison_table.csv",
                    mime="text/csv",
                )


with tab3:
    st.markdown("## About This Tool")
    st.markdown("""
    **Protein Disorder-to-Function Linker** is a bioinformatics pipeline that integrates
    multiple layers of protein information to identify regions where **structural disorder,
    structural disorder, biological annotations, PTM sites, and UniProt variant or mutagenesis annotations co-occur**.

    ---

    ### 🔬 Pipeline Steps

    | Step | Module | Description |
    |------|--------|-------------|
    | 1 | `uniprot_api.py` | Fetch sequence, domains, motifs, binding/PTM sites, variants |
    | 2 | `disorder_api.py` | Per-residue disorder via **IUPred3** REST API |
    | 2b | `disorder_api.py` | Per-residue structural confidence via **AlphaFold** pLDDT |
    | 3 | `integrator.py` | Align all layers into one per-residue table |
    | 4 | `scorer.py` | DFLS scoring + variant-dense window detection + Fisher's exact enrichment test |
    | 5 | `visualizer.py` | Interactive Plotly chart + linear protein map |

    ---

    ### 🏅 Disorder-Function Linkage Score (DFLS)

    | Component | Weight | Description |
    |-----------|--------|-------------|
    | Local disorder score | 40% | IUPred3 score at the exact annotated residue |
    | IDR mean score | 25% | Average disorder of the surrounding region |
    | Overlap fraction | 20% | Fraction of the feature within the IDR |
    | Disease association | 15% | Boost for disease-linked variants |

    Tiers: **High** ≥ 0.75 · **Medium** ≥ 0.50 · **Low** < 0.50

    > ⚠️ **Disclaimer:** DFLS weights are **heuristic and literature-derived**, not optimised
    > against a labelled benchmark. They are transparent and interpretable but not statistically
    > calibrated. Use the Fisher's exact test results (shown after each analysis) as independent
    > statistical evidence for IDR enrichment.

    ---

    ### 📊 Statistical Validation

    For every protein, a **two-sided Fisher's exact test** (scipy) is computed:

    | Test | Null hypothesis | Significant result |
    |------|-----------------|--------------------|
    | Variant enrichment | UniProt variants are uniformly distributed | OR > 1, p < 0.05 → enriched in IDRs |
    | PTM enrichment | PTM sites are uniformly distributed | OR > 1, p < 0.05 → enriched in IDRs |

    **Contingency table:**
    ```
                  IDR (≥0.5)    Ordered (<0.5)
    Has feature       a               b
    No  feature       c               d
    ```
    Odds ratio = (a·d) / (b·c). p-value from exact combinatorial calculation.

    ---

    ### ⚠️ "Variant-dense IDR window" — NOT a mutation cluster

    The sliding-window analysis (window=15, min=3 variants) detects windows
    of **high UniProt annotation density** — it is **not** equivalent to somatic
    mutation hotspot detection (e.g. COSMIC recurrence). It reflects how many
    neighbouring residues have UniProt variant or mutagenesis entries.

    ---

    ### 🔗 Data Sources
    - [UniProt REST API](https://www.uniprot.org/help/api) — Protein annotations
    - [IUPred3](https://iupred3.elte.hu/) — Intrinsic disorder prediction
    - [AlphaFold EBI](https://alphafold.ebi.ac.uk/) — Structure confidence (pLDDT)

    ---
    *CMPS 297AF · Path A · Protein Disorder-to-Function Linker*
    """)
