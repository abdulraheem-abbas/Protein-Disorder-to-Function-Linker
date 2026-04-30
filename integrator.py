"""Per-residue integration of disorder scores and UniProt annotations."""

from collections import defaultdict


def integrate(
    sequence: str,
    disorder: list[float],
    domains: list[tuple],
    motifs: list[tuple],
    binding_sites: list[tuple],
    mutations: list[tuple],
    ptm_sites: list[tuple] | None = None,
    plddt_scores: list[float] | None = None,
    verbose: bool = True,
) -> list[dict]:
    """Align disorder scores and annotations by residue position."""
    _validate(sequence, disorder)

    if plddt_scores is not None and len(plddt_scores) != len(sequence):
        if verbose:
            print(f"[WARNING] pLDDT length ({len(plddt_scores)}) != sequence length "
                  f"({len(sequence)}) — pLDDT will be None for all residues")
        plddt_scores = None

    maps = build_feature_maps(domains, motifs, binding_sites, mutations, ptm_sites or [])

    return integrate_residues(sequence, disorder, maps, plddt_scores, verbose=verbose)


def build_feature_maps(
    domains: list[tuple],
    motifs: list[tuple],
    binding_sites: list[tuple],
    mutations: list[tuple],
    ptm_sites: list[tuple] | None = None,
) -> dict:
    """Build position-keyed annotation maps."""
    return {
        "domain":       _build_range_map(domains),
        "motif":        _build_range_map(motifs),
        "binding_site": _build_point_map(binding_sites),
        "ptm":          _build_point_map(ptm_sites or []),
        "mutation":     _build_point_map(mutations),
    }


def _build_range_map(features: list[tuple]) -> dict[int, list[str]]:
    pos_map: dict[int, list[str]] = defaultdict(list)
    for start, end, name in features:
        if start is None or end is None:
            continue
        for pos in range(start, end + 1):
            pos_map[pos].append(name)
    return dict(pos_map)


def _build_point_map(features: list[tuple]) -> dict[int, list[str]]:
    pos_map: dict[int, list[str]] = defaultdict(list)
    for pos, description in features:
        if pos is not None:
            pos_map[pos].append(description)
    return dict(pos_map)


def integrate_residues(
    sequence: str,
    disorder: list[float],
    maps: dict,
    plddt_scores: list[float] | None = None,
    verbose: bool = True,
) -> list[dict]:
    """Build the per-residue output list using precomputed lookup maps."""
    domain_map   = maps["domain"]
    motif_map    = maps["motif"]
    binding_map  = maps["binding_site"]
    ptm_map      = maps["ptm"]
    mutation_map = maps["mutation"]

    result = []
    for i, (aa, score) in enumerate(zip(sequence, disorder)):
        pos = i + 1

        functional = sorted(
            domain_map.get(pos, []) +
            motif_map.get(pos, []) +
            binding_map.get(pos, [])
        )

        result.append({
            "position":   pos,
            "aa":         aa,
            "disorder":   round(score, 4),
            "plddt":      round(plddt_scores[i], 4) if plddt_scores else None,
            "functional": functional,
            "ptm":        ptm_map.get(pos, []),
            "mutation":   mutation_map.get(pos, []),
        })

    if verbose:
        n_ptm = sum(1 for r in result if r["ptm"])
        print(f"[INFO] Integrated {len(sequence)} residues "
              f"({n_ptm} with PTM annotations)")

    return result


def get_disordered(residues: list[dict], threshold: float = 0.5) -> list[dict]:
    """Return only residues with disorder score >= threshold."""
    return [r for r in residues if r["disorder"] >= threshold]


def get_functional_in_disorder(residues: list[dict], threshold: float = 0.5) -> list[dict]:
    """Return disordered residues that also carry at least one functional annotation."""
    return [r for r in residues if r["disorder"] >= threshold and r["functional"]]


def get_mutations_in_disorder(residues: list[dict], threshold: float = 0.5) -> list[dict]:
    """Return disordered residues that have at least one mutation annotation."""
    return [r for r in residues if r["disorder"] >= threshold and r["mutation"]]


def get_ptm_in_disorder(residues: list[dict], threshold: float = 0.5) -> list[dict]:
    """Return disordered residues that carry at least one PTM annotation."""
    return [r for r in residues if r["disorder"] >= threshold and r["ptm"]]


def _validate(sequence: str, disorder: list[float]) -> None:
    """Raise ValueError if inputs are empty or lengths don't match."""
    if not sequence:
        raise ValueError("Sequence is empty.")
    if not disorder:
        raise ValueError("Disorder scores list is empty.")
    if len(sequence) != len(disorder):
        raise ValueError(
            f"Length mismatch: sequence has {len(sequence)} residues "
            f"but disorder has {len(disorder)} scores."
        )
