"""Multi-protein comparison helpers."""

import os
import sys
import csv

sys.path.insert(0, os.path.dirname(__file__))
from main import run_pipeline


def compare_proteins(
    accessions: list[str],
    output_dir: str = "outputs",
    verbose: bool = True,
    fetch_plddt: bool = True,
    disorder_threshold: float = 0.5,
) -> list[dict]:
    """Run the pipeline for several proteins and return summary rows."""
    os.makedirs(output_dir, exist_ok=True)
    rows = []

    for acc in accessions:
        if verbose:
            print(f"\n{'='*55}")
            print(f"  Comparing: {acc}")
            print(f"{'='*55}")

        result = run_pipeline(
            acc,
            verbose=verbose,
            save_json=False,
            save_csv=False,
            fetch_plddt=fetch_plddt,
            disorder_threshold=disorder_threshold,
            output_dir=output_dir,
        )
        if result is None:
            print(f"[WARNING] Pipeline failed for {acc} — skipping.")
            continue

        s = result["scored"]["summary"]
        rows.append({
            "accession":          acc,
            "name":               result["name"],
            "gene_name":          result.get("gene_name", ""),
            "organism":           result.get("organism", ""),
            "length":             result["length"],
            "pct_disordered":     s.get("pct_disordered", 0),
            "n_disordered":       s.get("n_disordered", 0),
            "n_functional_in_idr":s.get("n_functional_in_idr", 0),
            "n_ptm_in_idr":       s.get("n_ptm_in_idr", 0),
            "n_mutations_in_idr": s.get("n_mutations_in_idr", 0),
            "n_idr_mutation_residues": s.get("n_idr_mutation_residues", 0),
            "n_variant_dense_windows": s.get("n_variant_dense_windows", 0),
            "n_concordant_idr":   s.get("n_concordant_idr", 0),
            "plot_file":          result["plot_file"],
            "linear_map_file":    result.get("linear_map_file", ""),
            "disorder_threshold": result.get("disorder_threshold", disorder_threshold),
            "disorder_metadata":  result.get("disorder_metadata", {}),
        })

    if not rows:
        print("[ERROR] No proteins were successfully analyzed.")
        return []

    rows.sort(key=lambda x: x["pct_disordered"], reverse=True)

    _save_comparison_csv(rows, output_dir, verbose)

    _plot_comparison(rows, output_dir, verbose)

    _plot_idr_vs_mutations(rows, output_dir, verbose)

    if verbose:
        _print_comparison_table(rows)

    return rows


def build_comparison_plotly(rows: list[dict]):
    """Build the comparison chart for Streamlit."""
    import plotly.graph_objects as go

    labels = [f"{r['gene_name'] or r['accession']}" for r in rows]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        name="% Disordered",
        x=labels,
        y=[r["pct_disordered"] for r in rows],
        marker_color="#2c72b5",
        text=[f"{r['pct_disordered']:.1f}%" for r in rows],
        textposition="outside",
    ))
    fig.add_trace(go.Bar(
        name="UniProt variants in IDRs",
        x=labels,
        y=[r["n_mutations_in_idr"] for r in rows],
        marker_color="#e05c5c",
        text=[str(r["n_mutations_in_idr"]) for r in rows],
        textposition="outside",
    ))
    fig.add_trace(go.Bar(
        name="PTM sites in IDRs",
        x=labels,
        y=[r["n_ptm_in_idr"] for r in rows],
        marker_color="#9b59b6",
        text=[str(r["n_ptm_in_idr"]) for r in rows],
        textposition="outside",
    ))
    fig.add_trace(go.Bar(
        name="Variant-dense annotation windows",
        x=labels,
        y=[r["n_variant_dense_windows"] for r in rows],
        marker_color="#ff8c00",
        text=[str(r["n_variant_dense_windows"]) for r in rows],
        textposition="outside",
    ))

    fig.update_layout(
        title="Multi-Protein Disorder & Function Comparison",
        barmode="group",
        xaxis_title="Protein",
        yaxis_title="Count / %",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        plot_bgcolor="white",
        height=450,
        margin=dict(l=50, r=20, t=80, b=50),
    )
    fig.update_yaxes(showgrid=True, gridcolor="rgba(0,0,0,0.06)")
    return fig


def _save_comparison_csv(rows: list[dict], output_dir: str, verbose: bool) -> None:
    path = os.path.join(output_dir, "comparison_table.csv")
    fields = [
        "accession", "name", "gene_name", "organism", "length",
        "pct_disordered", "n_disordered", "n_functional_in_idr",
        "n_ptm_in_idr", "n_mutations_in_idr", "n_variant_dense_windows", "n_concordant_idr",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    if verbose:
        print(f"\n[INFO] Comparison table saved → {path}")


def _plot_comparison(rows: list[dict], output_dir: str, verbose: bool) -> None:
    """Grouped bar chart comparing key metrics across proteins."""
    import matplotlib.pyplot as plt
    import numpy as np

    labels = [r.get("gene_name") or r["accession"] for r in rows]
    x      = np.arange(len(labels))
    w      = 0.20

    fig, ax = plt.subplots(figsize=(max(8, len(rows) * 2), 5))

    ax.bar(x - 1.5*w, [r["pct_disordered"]           for r in rows], w,
           label="% Disordered",                      color="#2c72b5", alpha=0.85)
    ax.bar(x - 0.5*w, [r["n_mutations_in_idr"]        for r in rows], w,
           label="UniProt variants in IDRs",          color="#e05c5c", alpha=0.85)
    ax.bar(x + 0.5*w, [r["n_ptm_in_idr"]              for r in rows], w,
           label="PTM sites in IDRs",                 color="#9b59b6", alpha=0.85)
    ax.bar(x + 1.5*w, [r["n_variant_dense_windows"]   for r in rows], w,
           label="Variant-dense annotation windows",  color="#ff8c00", alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("Count / %", fontsize=11)
    ax.set_title("Multi-Protein Disorder & Function Comparison", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    path = os.path.join(output_dir, "comparison_bar.png")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    if verbose:
        print(f"[INFO] Comparison bar chart saved → {path}")


def _plot_idr_vs_mutations(rows: list[dict], output_dir: str, verbose: bool) -> None:
    """Plot % IDR against UniProt variant annotations in IDRs."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 5))

    for r in rows:
        label = r.get("gene_name") or r["accession"]
        ax.scatter(r["pct_disordered"], r["n_mutations_in_idr"],
                   s=80, zorder=3, color="#2c72b5")
        ax.annotate(label,
                    (r["pct_disordered"], r["n_mutations_in_idr"]),
                    textcoords="offset points", xytext=(5, 4), fontsize=8.5)

    ax.set_xlabel("% Disordered residues",                         fontsize=11)
    ax.set_ylabel("UniProt variant/mutagenesis annotations in IDRs", fontsize=11)
    ax.set_title("Disorder vs. UniProt Variant Annotations in IDRs",  fontsize=12, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    path = os.path.join(output_dir, "comparison_scatter.png")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    if verbose:
        print(f"[INFO] Scatter plot saved → {path}")


def _print_comparison_table(rows: list[dict]) -> None:
    print()
    print("━" * 80)
    print(f"  {'Gene':<10} {'Accession':<10} {'Length':>6}  "
          f"{'% IDR':>7}  {'Var/IDR':>7}  {'PTM/IDR':>7}  {'Ann.Win':>8}")
    print("━" * 80)
    for r in rows:
        g = r.get("gene_name") or "—"
        print(f"  {g:<10} {r['accession']:<10} {r['length']:>6}  "
              f"{r['pct_disordered']:>6.1f}%  "
              f"{r['n_mutations_in_idr']:>7}  "
              f"{r['n_ptm_in_idr']:>7}  "
              f"{r.get('n_variant_dense_windows', 0):>8}")
    print("━" * 80)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args:
        print("Usage: python compare.py P04637 P37840 P46527 [--no-plddt]")
        sys.exit(1)

    no_plddt = "--no-plddt" in sys.argv
    compare_proteins(args, fetch_plddt=not no_plddt)
