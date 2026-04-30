"""Generate example analyses for the demo proteins."""

import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.dirname(__file__))
from main import run_pipeline

EXAMPLES = [
    {
        "accession":   "P04637",
        "short_name":  "TP53",
        "description": "Tumour suppressor p53; hallmark of cancer biology with "
                       "well-characterised disordered N- and C-terminal regions.",
    },
    {
        "accession":   "P37840",
        "short_name":  "α-Synuclein (SNCA)",
        "description": "Intrinsically disordered neuronal protein; aggregation "
                       "drives Parkinson's and Lewy body diseases.",
    },
    {
        "accession":   "P46527",
        "short_name":  "p27 Kip1 (CDKN1B)",
        "description": "Cell-cycle inhibitor (CDK2/Cyclin A); extensively disordered "
                       "yet tightly regulated through binding-induced folding.",
    },
]


def main() -> None:
    os.makedirs("outputs", exist_ok=True)
    summary_path = os.path.join("outputs", "example_analyses.txt")
    lines = []

    header = (
        "=" * 60 + "\n"
        "  PROTEIN DISORDER-TO-FUNCTION LINKER\n"
        "  Example Analyses\n"
        f"  Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n"
        + "=" * 60
    )
    print(header)
    lines.append(header)

    for i, ex in enumerate(EXAMPLES, 1):
        acc  = ex["accession"]
        name = ex["short_name"]
        desc = ex["description"]

        section = f"\n[{i}] {name}  ({acc})\n{desc}\n{'-' * 60}"
        print(section)
        lines.append(section)

        result = run_pipeline(acc, verbose=True)

        if result is None:
            msg = f"  ✗ Pipeline failed for {acc}.\n"
            print(msg)
            lines.append(msg)
            continue

        scored   = result["scored"]["summary"]
        residues = result["integrated_data"]

        try:
            from Bio.SeqUtils.ProtParam import ProteinAnalysis
            analysis = ProteinAnalysis(result["integrated_data"][0]["aa"]
                                        and "".join(r["aa"] for r in residues))
            mw = f"{analysis.molecular_weight():,.1f} Da"
            pi = f"{analysis.isoelectric_point():.2f}"
        except Exception:
            mw = pi = "N/A"

        block = (
            f"  Protein name        : {result['name']}\n"
            f"  Sequence length     : {result['length']} aa\n"
            f"  Molecular weight    : {mw}\n"
            f"  Isoelectric point   : {pi}\n"
            f"  Disordered (≥0.5)  : {scored['n_disordered']} residues\n"
            f"  UniProt variant-annotated IDR residues: {scored.get('n_idr_mutation_residues', scored.get('n_hotspots',0))}\n"
            f"  Func. sites in IDRs : {scored['n_functional_in_idr']}\n"
        )
        if scored.get("top_idr_mutation_region"):
            s, e = scored["top_idr_mutation_region"]
            block += f"  Largest IDR annotation-dense region: residues {s}–{e}\n"

        block += f"  Visualization       : {result['plot_file']}\n"

        print(block)
        lines.append(block)

    with open(summary_path, "w") as f:
        f.write("\n".join(lines))

    print(f"\n✓ All example analyses complete.")
    print(f"✓ Summary saved → {summary_path}")
    print(f"✓ Plots saved   → outputs/")


if __name__ == "__main__":
    main()
