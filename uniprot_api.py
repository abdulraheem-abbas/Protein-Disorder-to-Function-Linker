"""UniProt and raw-sequence input helpers."""

import json
import os
import requests

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
CACHE_DIR   = "data"

FIELDS = ",".join([
    "accession",
    "protein_name",
    "gene_names",
    "organism_name",
    "sequence",
    "ft_domain",
    "ft_region",
    "ft_motif",
    "ft_repeat",
    "ft_compbias",
    "ft_binding",
    "ft_act_site",
    "ft_mod_res",
    "ft_carbohyd",
    "ft_lipid",
    "ft_crosslnk",
    "ft_variant",
    "ft_mutagen",
])

DOMAIN_TYPES   = {"Domain", "Region"}
MOTIF_TYPES    = {"Motif", "Repeat", "Compositional bias"}
BINDING_TYPES  = {"Binding site", "Active site"}
PTM_TYPES      = {"Modified residue", "Glycosylation", "Lipidation", "Cross-link",
                   "Carbohyd", "Lipid"}
MUTATION_TYPES = {"Natural variant", "Mutagenesis"}


def fetch_protein(accession: str, verbose: bool = True) -> dict:
    """Fetch and parse a UniProt entry."""
    accession = accession.strip().upper()
    raw = _fetch_raw(accession, verbose=verbose)

    sequence = raw.get("sequence", {}).get("value", "")
    if not sequence:
        raise ValueError(f"No sequence found for {accession}")

    features = raw.get("features", [])

    if verbose:
        print(f"[INFO] Loaded protein {accession} (length={len(sequence)})")

    return {
        "name":          _extract_name(raw),
        "gene_name":     _extract_gene(raw),
        "organism":      _extract_organism(raw),
        "sequence":      sequence,
        "length":        len(sequence),
        "domains":       _extract_domains(features),
        "motifs":        _extract_motifs(features),
        "binding_sites": _extract_binding_sites(features),
        "ptm_sites":     _extract_ptm_sites(features),
        "mutations":     _extract_mutations(features),
    }


def build_protein_from_sequence(raw_input: str, name: str = "User-provided sequence") -> dict:
    """Build a protein dict from raw sequence or FASTA input."""
    from io import StringIO
    from Bio import SeqIO
    from Bio.Seq import Seq

    raw_input = raw_input.strip()
    if not raw_input:
        raise ValueError("Input is empty.")

    if raw_input.startswith(">"):
        record  = SeqIO.read(StringIO(raw_input), "fasta")
        seq_str = str(record.seq).upper().replace(" ", "").replace("\n", "")
        name    = record.description if record.description else name
    else:
        seq_str = raw_input.upper().replace(" ", "").replace("\n", "")

    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    invalid  = set(seq_str) - valid_aa
    if invalid:
        raise ValueError(
            f"Sequence contains invalid characters: {', '.join(sorted(invalid))}. "
            "Only standard single-letter amino-acid codes are accepted."
        )
    if not seq_str:
        raise ValueError("Sequence is empty after parsing.")

    return {
        "name":          name,
        "gene_name":     "",
        "organism":      "",
        "sequence":      str(Seq(seq_str)),
        "length":        len(seq_str),
        "domains":       [],
        "motifs":        [],
        "binding_sites": [],
        "ptm_sites":     [],
        "mutations":     [],
    }


def _fetch_raw(accession: str, verbose: bool = True) -> dict:
    """Fetch the UniProt JSON entry, using disk cache if available."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    cache_path = os.path.join(CACHE_DIR, f"{accession}.json")

    if os.path.exists(cache_path):
        if verbose:
            print(f"[INFO] Loading {accession} from cache ({cache_path})")
        with open(cache_path, "r") as f:
            return json.load(f)

    url = f"{UNIPROT_API}/{accession}"
    try:
        resp = requests.get(
            url, params={"fields": FIELDS},
            headers={"Accept": "application/json"}, timeout=30
        )
    except requests.exceptions.RequestException as exc:
        raise ValueError(f"Failed to fetch data for {accession}: {exc}") from exc

    if resp.status_code == 404:
        raise ValueError(f"UniProt accession not found: '{accession}'")
    if not resp.ok:
        raise ValueError(f"Failed to fetch data for {accession}: HTTP {resp.status_code}")

    raw = resp.json()

    with open(cache_path, "w") as f:
        json.dump(raw, f)
    if verbose:
        print(f"[INFO] Cached {accession} → {cache_path}")

    return raw


def _extract_name(raw: dict) -> str:
    """Extract the protein name."""
    try:
        return raw["proteinDescription"]["recommendedName"]["fullName"]["value"]
    except KeyError:
        try:
            return raw["proteinDescription"]["submissionNames"][0]["fullName"]["value"]
        except (KeyError, IndexError):
            return "Unknown protein"


def _extract_gene(raw: dict) -> str:
    """Extract the primary gene name."""
    try:
        genes = raw.get("genes", [])
        if genes:
            return genes[0].get("geneName", {}).get("value", "")
    except Exception:
        pass
    return ""


def _extract_organism(raw: dict) -> str:
    """Extract the organism scientific name."""
    try:
        return raw["organism"]["scientificName"]
    except (KeyError, TypeError):
        return ""


def _location(feat: dict) -> tuple[int | None, int | None]:
    """Return (start, end) 1-indexed positions from a feature dict."""
    loc   = feat.get("location", {})
    start = loc.get("start", {}).get("value")
    end   = loc.get("end",   {}).get("value")
    return start, end


def _extract_domains(features: list[dict]) -> list[tuple]:
    """Return [(start, end, name)] for Domain and Region annotations."""
    result = []
    for feat in features:
        if feat.get("type") not in DOMAIN_TYPES:
            continue
        start, end = _location(feat)
        name = feat.get("description", feat.get("type", ""))
        if start and end:
            result.append((start, end, name))
    return result


def _extract_motifs(features: list[dict]) -> list[tuple]:
    """Return [(start, end, name)] for Motif, Repeat, and Compositional bias."""
    result = []
    for feat in features:
        if feat.get("type") not in MOTIF_TYPES:
            continue
        start, end = _location(feat)
        name = feat.get("description", feat.get("type", ""))
        if start and end:
            result.append((start, end, name))
    return result


def _extract_binding_sites(features: list[dict]) -> list[tuple]:
    """Return [(position, description)] for binding sites and active sites."""
    result = []
    for feat in features:
        if feat.get("type") not in BINDING_TYPES:
            continue
        start, _ = _location(feat)
        description = feat.get("description", feat.get("type", ""))
        ligand = feat.get("ligand", {})
        if ligand and ligand.get("name"):
            description = f"{description} [{ligand['name']}]".strip(" []")
        if start:
            result.append((start, description))
    return result


def _extract_ptm_sites(features: list[dict]) -> list[tuple]:
    """Return [(position, description)] for PTM annotations."""
    result = []
    for feat in features:
        if feat.get("type") not in PTM_TYPES:
            continue
        start, _ = _location(feat)
        description = feat.get("description", feat.get("type", "PTM"))
        if start:
            result.append((start, description))
    return result


def _extract_mutations(features: list[dict]) -> list[tuple]:
    """Return one UniProt variant or mutagenesis annotation per residue."""
    seen = set()
    result = []
    for feat in features:
        if feat.get("type") not in MUTATION_TYPES:
            continue
        start, _ = _location(feat)
        if not start or start in seen:
            continue
        seen.add(start)
        description = feat.get("description", "") or "Variant"
        result.append((start, description))
    return sorted(result, key=lambda x: x[0])
