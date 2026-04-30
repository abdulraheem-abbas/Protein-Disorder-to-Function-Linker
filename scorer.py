"""Scoring and enrichment summaries for integrated residue annotations."""

from __future__ import annotations
from typing import Any


_W_LOCAL_SCORE = 0.40
_W_IDR_MEAN    = 0.25
_W_OVERLAP     = 0.20
_W_DISEASE     = 0.15


class Scorer:
    """Summarize disorder/function overlap."""

    def score_residues(
        self,
        residues: list[dict],
        threshold: float = 0.5,
    ) -> dict[str, Any]:
        """Compute IDR overlap counts, PTM/variant lists, and dense windows."""
        disordered           = [r for r in residues if r["disorder"] >= threshold]
        idr_mut_residues     = [r for r in disordered if r["mutation"]]
        functional_idr       = [r for r in disordered if r["functional"]]
        ptm_in_idr           = [r for r in disordered if r.get("ptm")]

        concordant_idr = [
            r for r in disordered
            if r.get("plddt") is not None and r["plddt"] >= threshold
        ]

        top_cluster_region = self._largest_idr_mutation_region(idr_mut_residues)
        clusters           = self.detect_variant_dense_windows(residues, threshold=threshold)

        summary = {
            "n_disordered":              len(disordered),
            "pct_disordered":            round(100 * len(disordered) / len(residues), 1) if residues else 0,
            "n_hotspots":                len(idr_mut_residues),
            "n_idr_mutation_residues":   len(idr_mut_residues),
            "n_mutations_in_idr":        len(idr_mut_residues),
            "n_functional_in_idr":       len(functional_idr),
            "n_ptm_in_idr":              len(ptm_in_idr),
            "n_concordant_idr":          len(concordant_idr),
            "n_variant_dense_windows":                len(clusters),
            "top_hotspot_region":        top_cluster_region,
            "top_idr_mutation_region":   top_cluster_region,
            "mutation_source":           "UniProt REST API (Natural variant + Mutagenesis)",
        }

        return {
            "summary":            summary,
            "idr_mutation_residues": idr_mut_residues,
            "hotspots":           idr_mut_residues,
            "disordered":         disordered,
            "clusters":           clusters,
            "ptm_in_idr":         ptm_in_idr,
            "concordant_idr":     concordant_idr,
        }

    @staticmethod
    def detect_variant_dense_windows(
        residues: list[dict],
        window: int = 15,
        min_mutations: int = 3,
        threshold: float = 0.5,
        top_n: int = 5,
    ) -> list[dict]:
        """Find top non-overlapping UniProt variant-density windows."""
        if not residues or len(residues) < window:
            return []

        n           = len(residues)

        scored_windows = []
        for i in range(n - window + 1):
            win_res  = residues[i : i + window]
            n_annot  = sum(1 for r in win_res if r["mutation"])
            if n_annot < min_mutations:
                continue
            n_dis    = sum(1 for r in win_res if r["disorder"] >= threshold)
            is_idr   = (n_dis / len(win_res)) >= 0.5
            scored_windows.append((n_annot, i, win_res))

        if not scored_windows:
            return []

        scored_windows.sort(key=lambda x: x[0], reverse=True)
        selected: list[tuple[int, int]] = []
        result   = []

        for n_annot, i, win_res in scored_windows:
            if len(result) >= top_n:
                break
            end_i = i + window - 1
            if any(not (end_i < s or i > e) for s, e in selected):
                continue

            selected.append((i, end_i))
            start_pos = win_res[0]["position"]
            end_pos   = win_res[-1]["position"]
            n_dis     = sum(1 for r in win_res if r["disorder"] >= threshold)
            is_idr    = (n_dis / len(win_res)) >= 0.5
            annot_pos = [r["position"] for r in win_res if r["mutation"]]

            result.append({
                "start":               start_pos,
                "end":                 end_pos,
                "n_annotations":       n_annot,
                "n_mutations":         n_annot,
                "n_disordered":        n_dis,
                "is_idr_window":       is_idr,
                "is_idr_cluster":      is_idr,
                "annotation_density":  round(n_annot / window, 3),
                "positions":           annot_pos,
                "window_size":         window,
                "note": (
                    "UniProt annotation density window — NOT a somatic cancer hotspot. "
                    "Treat as hypothesis-generating only."
                ),
            })

        result.sort(key=lambda x: x["start"])
        return result

    def score(self, integrated: dict[str, Any]) -> dict[str, Any]:
        scored_sites    = self._score_sites(integrated.get("linked_sites", []))
        scored_domains  = self._score_domains(integrated.get("linked_domains", []))
        scored_variants = self._score_variants(integrated.get("linked_variants", []))

        return {
            "accession":       integrated.get("accession", ""),
            "name":            integrated.get("name", ""),
            "scored_sites":    scored_sites,
            "scored_domains":  scored_domains,
            "scored_variants": scored_variants,
            "summary":         self._summary(scored_sites, scored_domains, scored_variants),
        }

    def _score_sites(self, sites: list[dict]) -> list[dict]:
        scored = []
        for site in sites:
            local    = site.get("local_score") or site.get("idr_mean_score") or 0.0
            idr_mean = site.get("idr_mean_score") or 0.0
            overlap  = site.get("overlap_fraction") or 1.0
            dfls = (
                _W_LOCAL_SCORE * local + _W_IDR_MEAN * idr_mean + _W_OVERLAP * overlap
            ) / (_W_LOCAL_SCORE + _W_IDR_MEAN + _W_OVERLAP)
            scored.append({**site, "dfls": round(dfls, 4), "score_tier": self._tier(dfls)})
        return sorted(scored, key=lambda x: x["dfls"], reverse=True)

    def _score_domains(self, domains: list[dict]) -> list[dict]:
        scored = []
        for dom in domains:
            df   = dom.get("disorder_fraction", 0.0)
            scored.append({**dom, "dfls": round(df, 4), "score_tier": self._tier(df)})
        return sorted(scored, key=lambda x: x["dfls"], reverse=True)

    def _score_variants(self, variants: list[dict]) -> list[dict]:
        scored = []
        for var in variants:
            local         = var.get("local_score") or var.get("idr_mean_score") or 0.0
            idr_mean      = var.get("idr_mean_score") or 0.0
            disease_boost = 1.0 if var.get("disease_associated") else 0.0
            dfls = (
                _W_LOCAL_SCORE * local + _W_IDR_MEAN * idr_mean + _W_DISEASE * disease_boost
            ) / (_W_LOCAL_SCORE + _W_IDR_MEAN + _W_DISEASE)
            scored.append({**var, "dfls": round(dfls, 4), "score_tier": self._tier(dfls)})
        return sorted(scored, key=lambda x: x["dfls"], reverse=True)

    @staticmethod
    def _summary(sites, domains, variants) -> dict:
        return {
            "n_scored_sites":               len(sites),
            "n_high_confidence_sites":      len([s for s in sites if s["score_tier"] == "high"]),
            "n_scored_domains":             len(domains),
            "n_high_confidence_domains":    len([d for d in domains if d["score_tier"] == "high"]),
            "n_scored_variants":            len(variants),
            "n_disease_variants_in_idr":    len([v for v in variants if v.get("disease_associated")]),
        }

    @staticmethod
    def _largest_idr_mutation_region(
        idr_mut_residues: list[dict],
        gap: int = 5,
    ) -> tuple[int, int] | None:
        """Find the largest nearby span of IDR residues carrying variants."""
        if not idr_mut_residues:
            return None

        positions = sorted(r["position"] for r in idr_mut_residues)
        best  = (positions[0], positions[0])
        start = positions[0]
        prev  = positions[0]

        for pos in positions[1:]:
            if pos - prev <= gap:
                prev = pos
            else:
                if prev - start > best[1] - best[0]:
                    best = (start, prev)
                start = pos
                prev  = pos

        if prev - start > best[1] - best[0]:
            best = (start, prev)

        return best

    @staticmethod
    def _tier(score: float) -> str:
        if score >= 0.75:
            return "high"
        if score >= 0.50:
            return "medium"
        return "low"

def compute_enrichment_statistics(
    residues: list[dict],
    threshold: float = 0.5,
) -> dict:
    """Run Fisher's exact tests for variant and PTM enrichment in IDRs."""
    from scipy.stats import fisher_exact

    idr     = [r for r in residues if r["disorder"] >= threshold]
    ordered = [r for r in residues if r["disorder"] <  threshold]

    def _run_test(feature_key: str, label: str) -> dict:
        a = sum(1 for r in idr     if r.get(feature_key))
        c = len(idr)     - a
        b = sum(1 for r in ordered if r.get(feature_key))
        d = len(ordered) - b

        if a + b == 0:
            return {
                "odds_ratio":    None,
                "p_value":       None,
                "significant":   False,
                "table":         [[a, b], [c, d]],
                "interpretation": f"No {label} annotations found — test not applicable.",
            }

        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="two-sided")

        if p_value < 0.001:
            sig_label = "*** (p < 0.001)"
        elif p_value < 0.01:
            sig_label = "** (p < 0.01)"
        elif p_value < 0.05:
            sig_label = "* (p < 0.05)"
        else:
            sig_label = "ns (not significant)"

        if odds_ratio > 1 and p_value < 0.05:
            direction = f"UniProt {label} are significantly enriched in IDRs (OR > 1, p < 0.05)"
        elif odds_ratio < 1 and p_value < 0.05:
            direction = f"UniProt {label} are significantly depleted from IDRs (OR < 1, p < 0.05)"
        else:
            direction = "No statistically significant enrichment or depletion detected (p ≥ 0.05)"

        return {
            "odds_ratio":    round(float(odds_ratio), 3),
            "p_value":       round(float(p_value), 5),
            "significant":   p_value < 0.05,
            "significance":  sig_label,
            "table":         {
                "IDR_with":     a,
                "ordered_with": b,
                "IDR_without":  c,
                "ordered_without": d,
            },
            "interpretation": direction,
        }

    variant_result = _run_test("mutation", "variant/mutagenesis annotations")
    ptm_result     = _run_test("ptm",      "PTM sites")

    return {
        "variant_enrichment": variant_result,
        "ptm_enrichment":     ptm_result,
        "test_method":        "Fisher's exact test (two-sided), scipy.stats.fisher_exact",
        "idr_threshold":      threshold,
        "n_idr":              len(idr),
        "n_ordered":          len(ordered),
        "dfls_disclaimer": (
            "DFLS weights (local=0.40, IDR_mean=0.25, overlap=0.20, disease=0.15) "
            "are heuristic and literature-derived, NOT optimised against a benchmark. "
            "The Fisher's exact test above provides statistically grounded enrichment evidence "
            "independent of the DFLS score."
        ),
    }
