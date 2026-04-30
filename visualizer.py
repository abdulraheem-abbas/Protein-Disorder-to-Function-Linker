"""Static and interactive plots for residue-level results."""

import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "protein_linker_mpl"))

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec


def plot_disorder(
    residues: list[dict],
    output_path: str = "disorder_plot.png",
    protein_name: str = "Protein",
    threshold: float = 0.5,
    clusters: list[dict] | None = None,
) -> None:
    """Save the disorder profile plot."""
    positions = [r["position"] for r in residues]
    scores    = [r["disorder"]  for r in residues]
    has_plddt = any(r.get("plddt") is not None for r in residues)

    fig, ax = plt.subplots(figsize=(14, 4.5))

    for cl in (clusters or []):
        ax.axvspan(cl["start"], cl["end"], color="#ff8c00", alpha=0.15, zorder=1,
                   label="_cluster")

    hotspot_ranges = _find_hotspot_ranges(residues, threshold)
    for start, end in hotspot_ranges:
        ax.axvspan(start, end, color="yellow", alpha=0.20, zorder=1)

    _shade_functional_regions(ax, residues)

    if has_plddt:
        plddt_vals = [r["plddt"] if r["plddt"] is not None else 0.0 for r in residues]
        ax.plot(positions, plddt_vals, color="#e8822a", linewidth=0.9,
                linestyle="--", alpha=0.75, label="AlphaFold pLDDT proxy (1−pLDDT/100; high = low confidence)", zorder=3)

    ax.plot(positions, scores, color="#2c72b5", linewidth=1.3,
            label="IUPred3 disorder score (0=ordered, 1=disordered)", zorder=4)

    ax.axhline(threshold, color="#e05c5c", linewidth=1.0,
               linestyle=":", label=f"IDR threshold ({threshold}) — above = disordered", zorder=2)

    mut_positions = [r["position"] for r in residues if r["mutation"] and r["disorder"] > threshold]
    mut_scores    = [r["disorder"]  for r in residues if r["mutation"] and r["disorder"] > threshold]
    if mut_positions:
        ax.scatter(mut_positions, mut_scores, color="red", s=18, zorder=5,
                   label="UniProt variant in IDR (germline/experimental)")

    ptm_positions = [r["position"] for r in residues if r.get("ptm") and r["disorder"] > threshold]
    ptm_scores    = [r["disorder"]  for r in residues if r.get("ptm") and r["disorder"] > threshold]
    if ptm_positions:
        ax.scatter(ptm_positions, ptm_scores, color="#9b59b6", s=22,
                   marker="^", zorder=5, label="PTM site in IDR (UniProt: phospho/glyco/ubiq)")

    ax.set_xlabel("Residue Position", fontsize=11)
    ax.set_ylabel("Disorder Score",   fontsize=11)
    ax.set_title(f"{protein_name} — Disorder & Functional Map", fontsize=13, fontweight="bold")
    ax.set_xlim(min(positions), max(positions))
    ax.set_ylim(0, 1.05)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    _add_legend(ax, threshold,
                has_hotspots=bool(hotspot_ranges),
                has_mutations=bool(mut_positions),
                has_plddt=has_plddt,
                has_ptm=bool(ptm_positions),
                has_clusters=bool(clusters))

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Disorder plot saved → {output_path}")


def plot_linear_map(
    residues: list[dict],
    output_path: str = "linear_map.png",
    protein_name: str = "Protein",
    threshold: float = 0.5,
) -> None:
    """Save the linear protein map."""
    length = residues[-1]["position"] if residues else 1

    fig = plt.figure(figsize=(14, 3.5))
    fig.suptitle(f"{protein_name} — Linear Protein Map", fontsize=12, fontweight="bold")

    gs  = gridspec.GridSpec(4, 1, hspace=0.08, figure=fig)
    axes = [fig.add_subplot(gs[i]) for i in range(4)]
    track_labels = ["IDRs", "Domains & Motifs", "Mutations (IDR)", "PTM Sites (IDR)"]

    for ax, label in zip(axes, track_labels):
        ax.set_xlim(0, length)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_ylabel(label, fontsize=7.5, rotation=0, ha="right",
                      va="center", labelpad=60)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        if ax is not axes[-1]:
            ax.set_xticks([])

    axes[-1].set_xlabel("Residue Position", fontsize=9)

    ax_idr = axes[0]
    in_idr, idr_start = False, 0
    for r in residues:
        pos = r["position"]
        if r["disorder"] >= threshold and not in_idr:
            in_idr, idr_start = True, pos
        elif r["disorder"] < threshold and in_idr:
            ax_idr.broken_barh([(idr_start, pos - idr_start)], (0.1, 0.8),
                                facecolors="#4a90d9", alpha=0.7, edgecolors="none")
            in_idr = False
    if in_idr:
        ax_idr.broken_barh([(idr_start, length - idr_start + 1)], (0.1, 0.8),
                            facecolors="#4a90d9", alpha=0.7, edgecolors="none")

    ax_dom = axes[1]
    domain_colors = ["#e74c3c", "#2ecc71", "#f39c12", "#9b59b6", "#1abc9c", "#e67e22"]
    seen_names: dict[str, str] = {}
    color_idx = 0

    in_reg, reg_start, reg_name = False, 0, ""
    for r in residues:
        pos    = r["position"]
        labels = r["functional"]
        name   = labels[0] if labels else None

        if name and not in_reg:
            in_reg, reg_start, reg_name = True, pos, name
        elif (not name or name != reg_name) and in_reg:
            if reg_name not in seen_names:
                seen_names[reg_name] = domain_colors[color_idx % len(domain_colors)]
                color_idx += 1
            col = seen_names[reg_name]
            ax_dom.broken_barh([(reg_start, pos - reg_start)], (0.1, 0.8),
                                facecolors=col, alpha=0.75, edgecolors="white", linewidth=0.5)
            mid = reg_start + (pos - reg_start) / 2
            ax_dom.text(mid, 0.5, reg_name[:18], ha="center", va="center",
                        fontsize=5.5, color="white", fontweight="bold", clip_on=True)
            in_reg = bool(name)
            reg_start, reg_name = pos, name or ""

    if in_reg:
        if reg_name not in seen_names:
            seen_names[reg_name] = domain_colors[color_idx % len(domain_colors)]
        col = seen_names[reg_name]
        width = length - reg_start + 1
        ax_dom.broken_barh([(reg_start, width)], (0.1, 0.8),
                            facecolors=col, alpha=0.75, edgecolors="white", linewidth=0.5)
        mid = reg_start + width / 2
        ax_dom.text(mid, 0.5, reg_name[:18], ha="center", va="center",
                    fontsize=5.5, color="white", fontweight="bold", clip_on=True)

    ax_mut = axes[2]
    ax_mut.set_facecolor("#fff5f5")
    for r in residues:
        if r["mutation"] and r["disorder"] >= threshold:
            ax_mut.axvline(r["position"], color="#e74c3c", linewidth=1.0, alpha=0.8)

    ax_ptm = axes[3]
    ax_ptm.set_facecolor("#f9f5ff")
    for r in residues:
        if r.get("ptm") and r["disorder"] >= threshold:
            ax_ptm.axvline(r["position"], color="#9b59b6", linewidth=1.0, alpha=0.8)

    for ax in axes:
        ax.axhline(0.05, color="#cccccc", linewidth=0.5, zorder=0)

    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Linear map saved → {output_path}")


def build_plotly_figure(
    residues: list[dict],
    protein_name: str = "Protein",
    threshold: float = 0.5,
    clusters: list[dict] | None = None,
) -> "go.Figure":
    """Build the interactive Plotly figure used by Streamlit."""
    import plotly.graph_objects as go

    positions = [r["position"] for r in residues]
    scores    = [r["disorder"]  for r in residues]
    has_plddt = any(r.get("plddt") is not None for r in residues)

    def _disorder_label(score):
        if score >= 0.75: return "High disorder"
        if score >= 0.5:  return "Disordered (IDR)"
        if score >= 0.25: return "Low disorder"
        return "Ordered"

    hover_text = [
        f"<b>Position {r['position']}</b> ({r['aa']})<br>"
        f"<b>IUPred3:</b> {r['disorder']:.3f} → {_disorder_label(r['disorder'])}"
        + (f"<br><b>AlphaFold pLDDT proxy:</b> {r['plddt']:.3f} "
           f"({'likely disordered' if r['plddt']>=0.5 else 'likely ordered'})"
           if r.get("plddt") is not None else "")
        + (f"<br><b>Structural annotations (UniProt):</b> {', '.join(r['functional'])}"
           if r["functional"] else "")
        + (f"<br><b>PTM sites (UniProt):</b> {', '.join(r['ptm'])}"
           if r.get("ptm") else "")
        + (f"<br><b>Variants (UniProt):</b> {', '.join(r['mutation'])}"
           if r["mutation"] else "")
        for r in residues
    ]

    fig = go.Figure()

    for cl in (clusters or []):
        fig.add_vrect(
            x0=cl["start"], x1=cl["end"],
            fillcolor="rgba(255,140,0,0.12)", line_width=0,
            annotation_text="Variant-dense window (UniProt)",
            annotation_position="top left",
            annotation_font_size=9,
        )

    for start, end in _find_hotspot_ranges(residues, threshold):
        fig.add_vrect(
            x0=start, x1=end,
            fillcolor="rgba(255,235,0,0.18)", line_width=0,
        )

    _plotly_shade_regions(fig, residues, keyword="domain",
                          color="rgba(74,144,217,0.15)", label="Domain")

    _plotly_shade_regions(fig, residues, keyword="motif",
                          color="rgba(76,175,80,0.15)", label="Motif")

    fig.add_hline(
        y=threshold, line_dash="dot", line_color="#e05c5c",
        line_width=1.5,
        annotation_text=f"IDR threshold ({threshold})",
        annotation_position="bottom right",
        annotation_font_color="#e05c5c",
    )

    if has_plddt:
        plddt_vals = [r["plddt"] if r.get("plddt") is not None else None for r in residues]
        fig.add_trace(go.Scatter(
            x=positions, y=plddt_vals,
            mode="lines",
            name="AlphaFold pLDDT proxy (1−pLDDT/100)",
            line=dict(color="#e8822a", width=1.2, dash="dash"),
            opacity=0.8,
            hovertemplate="Pos %{x}<br>pLDDT proxy: %{y:.3f}<extra></extra>",
        ))

    fig.add_trace(go.Scatter(
        x=positions, y=scores,
        mode="lines",
        name=f"IUPred3 disorder (≥{threshold:g} = IDR)",
        line=dict(color="#2c72b5", width=2),
        hovertext=hover_text,
        hovertemplate="%{hovertext}<extra></extra>",
    ))

    mut_res = [r for r in residues if r["mutation"] and r["disorder"] > threshold]
    if mut_res:
        fig.add_trace(go.Scatter(
            x=[r["position"] for r in mut_res],
            y=[r["disorder"]  for r in mut_res],
            mode="markers",
            name="UniProt variant in IDR (germline/experimental)",
            marker=dict(color="red", size=7, symbol="circle"),
            hovertemplate="Mutation @ pos %{x}<extra></extra>",
        ))

    ptm_res = [r for r in residues if r.get("ptm") and r["disorder"] > threshold]
    if ptm_res:
        fig.add_trace(go.Scatter(
            x=[r["position"] for r in ptm_res],
            y=[r["disorder"]  for r in ptm_res],
            mode="markers",
            name="PTM site in IDR (phospho/glyco/ubiq)",
            marker=dict(color="#9b59b6", size=8, symbol="triangle-up"),
            hovertemplate="PTM @ pos %{x}<extra></extra>",
        ))

    total_variants = sum(1 for r in residues if r["mutation"])
    idr_variants   = sum(1 for r in residues if r["mutation"] and r["disorder"] >= threshold)
    variant_subtitle = (
        f"UniProt variant-annotated residues, total = {total_variants}; "
        f"IDR-associated = {idr_variants}"
    )

    fig.update_layout(
        title=dict(text=(
            f"<b>{protein_name}</b> — Disorder & Functional Map<br>"
            f"<span style='font-size:12px;color:#8892b0'>{variant_subtitle}</span>"
        ),
                   font=dict(size=15)),
        xaxis=dict(title="Residue Position", showgrid=False),
        yaxis=dict(title="Disorder Score", range=[0, 1.05], showgrid=True,
                   gridcolor="rgba(0,0,0,0.06)"),
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="right", x=1),
        plot_bgcolor="white",
        hovermode="x unified",
        height=420,
        margin=dict(l=50, r=20, t=70, b=50),
    )

    return fig


def _find_hotspot_ranges(
    residues: list[dict],
    threshold: float,
    gap: int = 5,
) -> list[tuple[int, int]]:
    hotspot_positions = [
        r["position"] for r in residues
        if r["disorder"] > threshold and r["mutation"]
    ]
    if not hotspot_positions:
        return []

    ranges = []
    start = hotspot_positions[0]
    prev  = hotspot_positions[0]

    for pos in hotspot_positions[1:]:
        if pos - prev <= gap:
            prev = pos
        else:
            ranges.append((start, prev))
            start = pos
            prev  = pos
    ranges.append((start, prev))
    return ranges


def _shade_functional_regions(ax: plt.Axes, residues: list[dict]) -> None:
    _shade_by_keyword(ax, residues, keyword="domain", color="#4a90d9", alpha=0.25)
    _shade_by_keyword(ax, residues, keyword="motif",  color="#4caf50", alpha=0.25)


def _shade_by_keyword(ax, residues, keyword, color, alpha):
    in_region, start_pos = False, None
    for r in residues:
        pos     = r["position"]
        matches = any(keyword in label.lower() for label in r["functional"])
        if matches and not in_region:
            in_region, start_pos = True, pos
        elif not matches and in_region:
            ax.axvspan(start_pos, pos - 1, color=color, alpha=alpha, zorder=1)
            in_region = False
    if in_region:
        ax.axvspan(start_pos, residues[-1]["position"], color=color, alpha=alpha, zorder=1)


def _plotly_shade_regions(fig, residues, keyword, color, label):
    in_region, start_pos = False, None
    for r in residues:
        pos     = r["position"]
        matches = any(keyword in lbl.lower() for lbl in r["functional"])
        if matches and not in_region:
            in_region, start_pos = True, pos
        elif not matches and in_region:
            fig.add_vrect(x0=start_pos, x1=pos - 1, fillcolor=color,
                          line_width=0, layer="below")
            in_region = False
    if in_region:
        end = residues[-1]["position"]
        fig.add_vrect(x0=start_pos, x1=end, fillcolor=color,
                      line_width=0, layer="below")


def _add_legend(ax, threshold, has_hotspots, has_mutations,
                has_plddt=False, has_ptm=False, has_clusters=False):
    handles = [
        ax.lines[0] if not has_plddt else ax.lines[1],
        mpatches.Patch(color="#4a90d9", alpha=0.6, label="UniProt domain/region"),
        mpatches.Patch(color="#4caf50", alpha=0.6, label="UniProt motif/repeat"),
    ]
    labels = ["IUPred3 disorder score", "UniProt domain/region", "UniProt motif/repeat"]

    if has_plddt:
        handles.insert(1, ax.lines[0])
        labels.insert(1, "AlphaFold proxy")

    if has_mutations:
        handles.append(plt.Line2D([0], [0], marker="o", color="w",
                                  markerfacecolor="red", markersize=7))
        labels.append("Mutation (in IDR)")

    if has_ptm:
        handles.append(plt.Line2D([0], [0], marker="^", color="w",
                                  markerfacecolor="#9b59b6", markersize=7))
        labels.append("PTM site (in IDR)")

    if has_hotspots:
        handles.append(mpatches.Patch(color="yellow", alpha=0.5, label="IDR mutation cluster"))
        labels.append("IDR mutation cluster")

    if has_clusters:
        handles.append(mpatches.Patch(color="#ff8c00", alpha=0.4))
        labels.append("Mutation cluster")

    ax.legend(handles, labels, loc="upper right", fontsize=8.5, framealpha=0.75)
