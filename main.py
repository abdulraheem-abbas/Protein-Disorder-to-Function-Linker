"""Command-line entry point for the disorder/function pipeline."""

import sys
import os
import re
import json
import argparse

from uniprot_api  import fetch_protein, build_protein_from_sequence
from disorder_api import get_disorder_scores, fetch_alphafold_plddt
from integrator   import integrate
from visualizer   import plot_disorder, plot_linear_map
from scorer       import Scorer, compute_enrichment_statistics

_ACCESSION_RE = re.compile(
    r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$',
    re.IGNORECASE
)
_AA_RE = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY]+$', re.IGNORECASE)


def _detect_input(user_input: str) -> str:
    """Return 'accession' or 'sequence' based on what the user typed."""
    cleaned = user_input.strip()
    if cleaned.startswith('>'):
        return 'sequence'
    cleaned = cleaned.replace(' ', '').replace('\n', '')
    if _ACCESSION_RE.match(cleaned) and len(cleaned) <= 10:
        return 'accession'
    if _AA_RE.match(cleaned) and len(cleaned) > 10:
        return 'sequence'
    return 'accession'


def run_pipeline(
    user_input: str,
    verbose: bool = True,
    save_json: bool = False,
    save_csv: bool = False,
    fetch_plddt: bool = True,
    disorder_threshold: float = 0.5,
    output_dir: str = "outputs",
    make_plots: bool = True,
) -> dict | None:
    """Run one protein analysis and return the collected results."""
    try:
        if not 0 <= disorder_threshold <= 1:
            raise ValueError("disorder_threshold must be between 0.0 and 1.0")

        input_type = _detect_input(user_input)
        label      = user_input[:10] + '...' if len(user_input) > 10 else user_input

        if verbose:
            print(f"[INFO] Input detected as: {input_type}")
            print(f"[INFO] Running pipeline for: {label}")

        if verbose:
            print("[INFO] Step 1: Fetching protein data...")
        if input_type == 'accession':
            protein    = fetch_protein(user_input.strip(), verbose=verbose)
            file_label = user_input.strip().upper()
            accession  = file_label
        else:
            sequence_clean = user_input.strip().replace(' ', '').replace('\n', '').upper()
            protein    = build_protein_from_sequence(sequence_clean)
            file_label = f"seq_{len(sequence_clean)}aa"
            accession  = None

        name          = protein["name"]
        gene_name     = protein.get("gene_name", "")
        organism      = protein.get("organism", "")
        sequence      = protein["sequence"]
        length        = protein["length"]
        domains       = protein["domains"]
        motifs        = protein["motifs"]
        binding_sites = protein["binding_sites"]
        ptm_sites     = protein.get("ptm_sites", [])
        mutations     = protein["mutations"]

        if verbose:
            print(f"[INFO] Protein : {name}")
            print(f"[INFO] Gene    : {gene_name}  |  Organism: {organism}")
            print(f"[INFO] Length  : {length} aa  |  "
                  f"PTM sites: {len(ptm_sites)}  |  Mutations: {len(mutations)}")

        if verbose:
            print("[INFO] Step 2: Computing disorder scores (IUPred3)...")
        disorder_scores, disorder_metadata = get_disorder_scores(
            sequence, verbose=verbose,
            accession=accession if input_type == "accession" else None,
            return_metadata=True,
        )
        if verbose and disorder_metadata.get("is_mock"):
            print("[WARNING] Disorder scores are mock fallback values; interpret results cautiously.")

        plddt_scores = None
        if fetch_plddt and accession:
            if verbose:
                print("[INFO] Step 2b: Fetching AlphaFold pLDDT scores...")
            plddt_scores = fetch_alphafold_plddt(accession, verbose=verbose)

        if verbose:
            print("[INFO] Step 3: Integrating annotations...")
        residues = integrate(
            sequence, disorder_scores,
            domains, motifs, binding_sites, mutations,
            ptm_sites=ptm_sites,
            plddt_scores=plddt_scores,
            verbose=verbose,
        )

        if verbose:
            print("[INFO] Step 4: Scoring disorder-function linkages...")
        scorer = Scorer()
        scored = scorer.score_residues(residues, threshold=disorder_threshold)

        if verbose:
            print("[INFO] Step 4b: Running Fisher's exact test for IDR enrichment...")
        enrichment = compute_enrichment_statistics(residues, threshold=disorder_threshold)
        if verbose:
            ve = enrichment["variant_enrichment"]
            pe = enrichment["ptm_enrichment"]
            if ve.get("p_value") is not None:
                pv_s = "p < 0.001" if ve["p_value"] < 0.001 else f"p = {ve['p_value']:.5f}"
                print(f"[INFO] Variant enrichment in IDR: OR={ve['odds_ratio']}, "
                      f"{pv_s} {ve['significance']}")
            if pe.get("p_value") is not None:
                pp_s = "p < 0.001" if pe["p_value"] < 0.001 else f"p = {pe['p_value']:.5f}"
                print(f"[INFO] PTM enrichment in IDR:     OR={pe['odds_ratio']}, "
                      f"{pp_s} {pe['significance']}")

        _print_summary(
            scored["summary"],
            name,
            gene_name,
            organism,
            sequence,
            verbose,
            enrichment=enrichment,
            threshold=disorder_threshold,
            disorder_metadata=disorder_metadata,
        )

        os.makedirs(output_dir, exist_ok=True)
        plot_file   = os.path.join(output_dir, f"{file_label}_plot.png")
        linear_file = os.path.join(output_dir, f"{file_label}_linear_map.png")

        if make_plots:
            if verbose:
                print("[INFO] Step 5: Generating disorder profile plot...")
            plot_disorder(
                residues,
                output_path=plot_file,
                protein_name=name,
                threshold=disorder_threshold,
                clusters=scored.get("clusters"),
            )

            if verbose:
                print("[INFO] Step 6: Generating linear protein map...")
            plot_linear_map(
                residues,
                output_path=linear_file,
                protein_name=name,
                threshold=disorder_threshold,
            )
        else:
            plot_file = ""
            linear_file = ""

        result = {
            "name":            name,
            "gene_name":       gene_name,
            "organism":        organism,
            "accession":       accession or file_label,
            "length":          length,
            "plot_file":       plot_file,
            "linear_map_file": linear_file,
            "disorder_threshold": disorder_threshold,
            "disorder_metadata": disorder_metadata,
            "integrated_data": residues,
            "scored":          scored,
            "enrichment":      enrichment,
        }

        if save_json:
            _export_json(result, file_label, verbose, output_dir=output_dir)
        if save_csv:
            _export_csv(residues, file_label, verbose, output_dir=output_dir)

        return result

    except Exception as e:
        print(f"[ERROR] {e}")
        import traceback; traceback.print_exc()
        return None


def _print_summary(
    summary: dict, protein_name: str, gene_name: str,
    organism: str, sequence: str, verbose: bool,
    enrichment: dict | None = None,
    threshold: float = 0.5,
    disorder_metadata: dict | None = None,
) -> None:
    if not verbose:
        return

    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        analysis = ProteinAnalysis(sequence.upper().replace('*', ''))
        mw   = analysis.molecular_weight()
        pi   = analysis.isoelectric_point()
        props = (f"  Molecular weight           : {mw:>10,.1f} Da\n"
                 f"  Isoelectric point (pI)     : {pi:>10.2f}")
    except Exception:
        props = "  Protein properties         : unavailable"

    clusters = summary.get("n_variant_dense_windows", summary.get("n_clusters", 0))

    print()
    print("━" * 60)
    print(f"  DISORDER-FUNCTION SUMMARY  │  {protein_name}")
    if gene_name:
        print(f"  Gene: {gene_name}  │  Organism: {organism}")
    print("━" * 60)
    print(props)
    print(f"  Disordered residues (>={threshold:g})  : {summary.get('n_disordered', 0)}"
          f"  ({summary.get('pct_disordered', 0):.1f}%)")
    print(f"  UniProt variant-annotated IDR residues : {summary.get('n_idr_mutation_residues', 0)}")
    print(f"  Functional sites in IDRs               : {summary.get('n_functional_in_idr', 0)}")
    print(f"  PTM sites in IDRs                      : {summary.get('n_ptm_in_idr', 0)}")
    print(f"  Variant-dense annotation windows       : {clusters}")
    if summary.get("n_concordant_idr", 0) > 0:
        print(f"  IUPred3+AlphaFold agree      : {summary['n_concordant_idr']} residues")
    if summary.get('top_idr_mutation_region'):
        s, e = summary['top_idr_mutation_region']
        print(f"  Largest IDR annotation-dense region    : residues {s}–{e}")
    print("  Mutation source              : UniProt (Natural variant + Mutagenesis)")
    if disorder_metadata:
        source = disorder_metadata.get("source", "unknown")
        print(f"  Disorder score source        : {source}")
        if disorder_metadata.get("is_mock"):
            print("  WARNING: Disorder scores are mock fallback values, not IUPred3 output.")
            reason = disorder_metadata.get("reason")
            if reason:
                print(f"           Reason: {reason}")
    print("  NOTE: Somatic mutations (COSMIC) and population variants (gnomAD)")
    print("        are NOT included in this analysis.")
    if enrichment and verbose:
        ve = enrichment.get("variant_enrichment", {})
        pe = enrichment.get("ptm_enrichment", {})
        print("━" * 60)
        print("  STATISTICAL ENRICHMENT  (Fisher's exact test, two-sided)")
        print("━" * 60)
        if ve.get("p_value") is not None:
            _pv = ve['p_value']
            _pv_str = 'p < 0.001' if _pv < 0.001 else f'p = {_pv:.5f}'
            print(f"  UniProt variants in IDRs : OR={ve['odds_ratio']:<8} {_pv_str:<14} {ve['significance']}")
            print(f"                      → {ve['interpretation']}")
        if pe.get("p_value") is not None:
            _pp = pe['p_value']
            _pp_str = 'p < 0.001' if _pp < 0.001 else f'p = {_pp:.5f}'
            print(f"  PTM sites in IDRs : OR={pe['odds_ratio']:<8} {_pp_str:<14} {pe['significance']}")
            print(f"                      → {pe['interpretation']}")
        print("  ⚠ DFLS weights are heuristic (not benchmark-validated).")
    print("━" * 60)
    print()


def _export_json(result: dict, file_label: str, verbose: bool, output_dir: str = "outputs") -> None:
    """Write full pipeline results to JSON (serializes residue list)."""
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"{file_label}_results.json")
    exportable = {k: v for k, v in result.items() if k not in ("integrated_data",)}
    exportable["residues"] = result["integrated_data"]
    with open(out_path, "w") as f:
        json.dump(exportable, f, indent=2, default=str)
    if verbose:
        print(f"[INFO] JSON results saved → {out_path}")


def _export_csv(residues: list[dict], file_label: str, verbose: bool, output_dir: str = "outputs") -> None:
    """Write per-residue table to CSV."""
    import csv
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"{file_label}_residues.csv")
    fieldnames = ["position", "aa", "disorder", "plddt", "functional", "ptm", "mutation"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in residues:
            writer.writerow({
                "position":   r["position"],
                "aa":         r["aa"],
                "disorder":   r["disorder"],
                "plddt":      r.get("plddt", ""),
                "functional": "; ".join(r["functional"]),
                "ptm":        "; ".join(r.get("ptm", [])),
                "mutation":   "; ".join(r["mutation"]),
            })
    if verbose:
        print(f"[INFO] CSV residues saved  → {out_path}")


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Link protein intrinsic disorder predictions to UniProt functional annotations."
    )
    parser.add_argument(
        "input",
        nargs="?",
        help="UniProt accession ID or raw amino-acid sequence.",
    )
    parser.add_argument(
        "--uniprot",
        dest="uniprot",
        help="UniProt accession ID. Equivalent to the positional input.",
    )
    parser.add_argument(
        "--disorder-threshold",
        type=float,
        default=0.5,
        help="IDR disorder cutoff used for scoring, enrichment, and plots. Default: 0.5.",
    )
    parser.add_argument(
        "--output",
        default="outputs",
        help="Directory for plots and exported files. Default: outputs.",
    )
    parser.add_argument("--quiet", action="store_true", help="Suppress INFO logs.")
    parser.add_argument("--save-json", action="store_true", help="Export full results as JSON.")
    parser.add_argument("--save-csv", action="store_true", help="Export per-residue data as CSV.")
    parser.add_argument("--no-plddt", action="store_true", help="Skip AlphaFold pLDDT fetch.")
    parser.add_argument("--no-plot", action="store_true", help="Skip PNG plot generation.")
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = _parse_args(sys.argv[1:])

    user_input = args.uniprot or args.input
    if not user_input:
        print("Enter a UniProt accession ID  (e.g. P04637)")
        print("  OR a raw amino-acid sequence (e.g. MAST...)")
        user_input = input("→ ").strip()

    if not user_input:
        print("[ERROR] No input provided.")
        sys.exit(1)

    run_pipeline(
        user_input,
        verbose=not args.quiet,
        save_json=args.save_json,
        save_csv=args.save_csv,
        fetch_plddt=not args.no_plddt,
        disorder_threshold=args.disorder_threshold,
        output_dir=args.output,
        make_plots=not args.no_plot,
    )
