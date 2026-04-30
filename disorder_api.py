"""Disorder and AlphaFold confidence fetchers."""

import random
import requests

IUPRED3_BASE  = "https://iupred3.elte.hu/iupred3"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"


def get_disorder_scores(
    sequence: str,
    pred_type: str = "long",
    verbose: bool = True,
    accession: str | None = None,
    return_metadata: bool = False,
) -> list[float] | tuple[list[float], dict]:
    """Return IUPred3 scores when possible, otherwise deterministic mock scores."""
    sequence = sequence.strip().upper()
    if not sequence:
        raise ValueError("Sequence is empty.")

    fallback_reason = None
    if accession:
        try:
            scores = _call_iupred3_by_accession(accession, pred_type)
            _validate_length(scores, sequence)
            if verbose:
                print(f"[INFO] IUPred3 scores fetched for {accession} "
                      f"(length={len(scores)})")
            metadata = {
                "source": "iupred3",
                "is_mock": False,
                "accession": accession,
                "pred_type": pred_type,
                "reason": None,
            }
            return (scores, metadata) if return_metadata else scores
        except Exception as exc:
            fallback_reason = str(exc)
            if verbose:
                print(f"[WARNING] IUPred3 unavailable — reason: {exc}")

    if verbose:
        if not accession:
            print("[INFO] No accession — using mock disorder scores (raw sequence mode)")
        else:
            print("[WARNING] Using mock disorder scores (seeded fallback)")
    scores = mock_disorder(sequence)
    metadata = {
        "source": "mock",
        "is_mock": True,
        "accession": accession,
        "pred_type": pred_type,
        "reason": fallback_reason or "No UniProt accession provided; raw sequences use mock disorder scores.",
    }
    return (scores, metadata) if return_metadata else scores


def fetch_alphafold_plddt(accession: str, verbose: bool = True) -> list[float] | None:
    """Return AlphaFold pLDDT as a disorder-style proxy, or None if unavailable."""
    accession = accession.strip().upper()
    try:
        meta_url = f"{ALPHAFOLD_API}/{accession}"
        resp = requests.get(meta_url, headers={"Accept": "application/json"}, timeout=30)
        if not resp.ok:
            if verbose:
                print(f"[INFO] AlphaFold: no entry for {accession} (HTTP {resp.status_code})")
            return None
        data = resp.json()
        if not data:
            return None

        pdb_url = data[0].get("pdbUrl")
        if not pdb_url:
            return None

        if verbose:
            print(f"[INFO] Fetching AlphaFold structure for {accession}...")
        pdb_resp = requests.get(pdb_url, timeout=60)
        if not pdb_resp.ok:
            return None

        plddt_raw = _parse_plddt_from_pdb(pdb_resp.text)
        if not plddt_raw:
            return None

        normalized = [round((100.0 - v) / 100.0, 4) for v in plddt_raw]
        if verbose:
            print(f"[INFO] AlphaFold pLDDT loaded ({len(normalized)} residues)")
        return normalized

    except Exception as exc:
        if verbose:
            print(f"[WARNING] AlphaFold pLDDT unavailable: {exc}")
        return None


def _call_iupred3_by_accession(accession: str, pred_type: str) -> list[float]:
    """Fetch IUPred3 disorder scores for a UniProt accession."""
    url  = f"{IUPRED3_BASE}/{pred_type}/{accession}"
    resp = requests.get(url, timeout=30)
    if not resp.ok:
        raise RuntimeError(
            f"IUPred3 returned HTTP {resp.status_code} for {accession}"
        )
    if "not valid IUPred type" in resp.text or "not found" in resp.text.lower():
        raise RuntimeError(f"IUPred3 rejected request: {resp.text[:120].strip()}")
    scores = _parse_text_response(resp.text)
    if not scores:
        raise RuntimeError("IUPred3 response parsed to empty list.")
    return scores


def _parse_plddt_from_pdb(pdb_text: str) -> list[float]:
    """Extract CA-atom pLDDT values from an AlphaFold PDB file."""
    try:
        import io
        from Bio.PDB import PDBParser

        parser   = PDBParser(QUIET=True)
        structure = parser.get_structure("AF", io.StringIO(pdb_text))

        plddt = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] != " ":
                        continue
                    if "CA" in residue:
                        plddt.append(residue["CA"].get_bfactor())
            break

        return plddt

    except Exception:
        plddt = {}
        for line in pdb_text.splitlines():
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            try:
                res_seq = int(line[22:26].strip())
                bfactor = float(line[60:66].strip())
                plddt[res_seq] = bfactor
            except (ValueError, IndexError):
                continue
        return [plddt[k] for k in sorted(plddt)] if plddt else []


def _parse_text_response(text: str) -> list[float]:
    """Parse plain-text IUPred output into disorder scores."""
    import re
    text = re.sub(r'<[^>]+>', '', text)

    scores = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            int(parts[0])
            scores.append(float(parts[2]))
        except ValueError:
            continue
    return scores


def _validate_length(scores: list[float], sequence: str) -> None:
    """Raise ValueError if score count doesn't match sequence length."""
    if len(scores) != len(sequence):
        raise ValueError(
            f"Disorder length mismatch: expected {len(sequence)}, got {len(scores)}"
        )


def mock_disorder(sequence: str) -> list[float]:
    """Return deterministic fallback disorder scores."""
    random.seed(42)
    return [round(random.uniform(0.3, 0.7), 4) for _ in sequence]
