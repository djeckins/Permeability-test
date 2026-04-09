"""Screen molecule records against follicle-oriented topical-delivery criteria.

All descriptor and pH-dependent value calculations are unchanged from the
original pipeline.  Only the threshold boundaries, class assignment, weighted
scoring and the PASS / BORDERLINE / FAIL decision have been updated.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import pandas as pd

from epidermal_barrier_screen.descriptors import calculate
from epidermal_barrier_screen.ionization import (
    PH_SC,
    _hhb_acid,
    _hhb_base,
    predict_pka,
)

# ---------------------------------------------------------------------------
# PAINS and Brenk structural alert filters  (RDKit FilterCatalog)
# ---------------------------------------------------------------------------
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

_pains_catalog: FilterCatalog | None = None
_brenk_catalog: FilterCatalog | None = None


def _get_pains_catalog() -> FilterCatalog:
    global _pains_catalog
    if _pains_catalog is None:
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
        _pains_catalog = FilterCatalog(params)
    return _pains_catalog


def _get_brenk_catalog() -> FilterCatalog:
    global _brenk_catalog
    if _brenk_catalog is None:
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        _brenk_catalog = FilterCatalog(params)
    return _brenk_catalog


def _check_pains(mol) -> str:
    """Return the PAINS alert name if matched, or empty string."""
    entry = _get_pains_catalog().GetFirstMatch(mol)
    return entry.GetDescription() if entry else ""


def _check_brenk(mol) -> str:
    """Return the Brenk alert name if matched, or empty string."""
    entry = _get_brenk_catalog().GetFirstMatch(mol)
    return entry.GetDescription() if entry else ""

# ---------------------------------------------------------------------------
# Scoring configuration
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class _CriterionCfg:
    weight: float   # contribution to the 100-point maximum raw score
    is_core: bool   # True → counted in CorePoorCount


# Edit weights / thresholds here; do NOT touch screen_records().
_CRITERIA: dict[str, _CriterionCfg] = {
    "MW":                _CriterionCfg(weight=17, is_core=True),
    "LogD":              _CriterionCfg(weight=17, is_core=True),
    "TPSA":              _CriterionCfg(weight=14, is_core=True),
    "FormalCharge":      _CriterionCfg(weight=14, is_core=True),
    "UnionizedFraction": _CriterionCfg(weight=14, is_core=True),
    "HBD":               _CriterionCfg(weight=7,  is_core=False),
    "HBA":               _CriterionCfg(weight=7,  is_core=False),
    "RotB":              _CriterionCfg(weight=10, is_core=False),
}

_MAX_RAW_SCORE: float = sum(c.weight for c in _CRITERIA.values())  # 100


# ---------------------------------------------------------------------------
# Per-criterion classification  →  "optimal" / "acceptable" / "poor"
# ---------------------------------------------------------------------------

def _classify_mw(v: float) -> str:
    if v <= 400:
        return "optimal"
    if v <= 500:
        return "acceptable"
    return "poor"


def _classify_logd(v: float | None) -> str:
    if v is None:
        return "poor"
    if 1.0 <= v <= 3.5:
        return "optimal"
    if (0.5 <= v < 1.0) or (3.5 < v <= 4.5):
        return "acceptable"
    return "poor"


def _classify_tpsa(v: float) -> str:
    if v <= 90:
        return "optimal"
    if v <= 120:
        return "acceptable"
    return "poor"


def _classify_formal_charge(v: int) -> str:
    if v == 0:
        return "optimal"
    if abs(v) == 1:
        return "acceptable"
    return "poor"


def _classify_unionized(v: float | None) -> str:
    """Classify the unionized percentage (0–100+%).

    Thresholds are in percent:
    - optimal:    >= 40 %
    - acceptable: >= 10 %
    - poor:       <  10 %
    """
    if v is None:
        return "poor"
    if v >= 40.0:
        return "optimal"
    if v >= 10.0:
        return "acceptable"
    return "poor"


def _classify_hbd(v: int) -> str:
    if v <= 2:
        return "optimal"
    if v == 3:
        return "acceptable"
    return "poor"


def _classify_hba(v: int) -> str:
    if 2 <= v <= 8:
        return "optimal"
    if v <= 1 or (8 < v <= 10):
        return "acceptable"
    return "poor"


def _compute_unionized_pct(logd: float | None, logp: float | None) -> float | None:
    """Compute the unionized percentage from logD and logP.

    Derivation:
        logD = logP + log10(f_neutral)   [standard pH-partition relationship]
        → f_neutral = 10 ** (logD - logP)
        → unionized_pct = f_neutral * 100

    - logD (row["logd"]): pH-adjusted, may come from pKa-corrected calculation,
      neutral default (logD = logP), or experimental SDF input.
    - logP (desc["clogp"]): RDKit cLogP, pH-independent.

    Result is in PERCENT (e.g. 27.54, 100.0).  NOT capped — values slightly
    above 100 % can arise from experimental logD/logP rounding and are kept as-is
    to preserve the numeric information.  Returns None when either input is
    missing, NaN, or non-numeric.
    """
    if logd is None or logp is None:
        return None
    try:
        ld = float(logd)
        lp = float(logp)
    except (TypeError, ValueError):
        return None
    import math as _math
    if _math.isnan(ld) or _math.isnan(lp):
        return None
    # fraction neutral = 10^(logD - logP); convert to percent, no cap
    return round(10.0 ** (ld - lp) * 100.0, 2)


def _classify_rotb(v: int) -> str:
    if v <= 7:
        return "optimal"
    if v <= 10:
        return "acceptable"
    return "poor"


# ---------------------------------------------------------------------------
# Legacy per-criterion status functions
# Kept for backward-compatible colour-coding in app.py.
# Thresholds mirror the new classification above; terminology kept as
# "suboptimal" so existing CSS keys in app.py continue to work.
# ---------------------------------------------------------------------------

def _mw_status(v: float) -> str:
    if v <= 400:   return "optimal"
    if v <= 500:   return "suboptimal"
    return "poor"


def _logd_status(v: float | None) -> str:
    if v is None:  return "poor"
    if 1.0 <= v <= 3.5:  return "optimal"
    if (0.5 <= v < 1.0) or (3.5 < v <= 4.5):  return "suboptimal"
    return "poor"


def _tpsa_status(v: float) -> str:
    if v <= 90:    return "optimal"
    if v <= 120:   return "suboptimal"
    return "poor"


def _hbd_status(v: int) -> str:
    if v <= 2:     return "optimal"
    if v == 3:     return "suboptimal"
    return "poor"


def _hba_status(v: int) -> str:
    if 2 <= v <= 8:        return "optimal"
    if v <= 1 or v <= 10:  return "suboptimal"
    return "poor"


def _rotb_status(v: int) -> str:
    if v <= 7:     return "optimal"
    if v <= 10:    return "suboptimal"
    return "poor"


def _hac_status(v: int) -> str:
    """Informational only — not used in weighted scoring."""
    if v < 30:     return "optimal"
    if v <= 50:    return "suboptimal"
    return "poor"


def _charge_status(v: int) -> str:
    if v == 0:         return "optimal"
    if abs(v) == 1:    return "suboptimal"
    return "poor"


def _ionization_status(pct_unionized: float | None) -> str:
    """Legacy status from unionized percentage (0–100+%)."""
    if pct_unionized is None:  return "poor"
    if pct_unionized >= 40.0:  return "optimal"
    if pct_unionized >= 10.0:  return "suboptimal"
    return "poor"


# ---------------------------------------------------------------------------
# Weighted scoring helpers
# ---------------------------------------------------------------------------

def _criterion_score(cls: str, weight: float) -> float:
    """Convert a classification to its numeric score contribution."""
    if cls == "optimal":    return weight
    if cls == "acceptable": return weight * 0.5
    return 0.0


def _compute_weighted_score(classes: dict[str, str]) -> float:
    """Return the weighted score on a 0–100 scale (weights already sum to 100).

    WeightedScore = raw_score  (no normalisation needed).
    """
    raw = sum(_criterion_score(classes[k], _CRITERIA[k].weight) for k in _CRITERIA)
    return round(raw, 1)


def _count_core_poor(classes: dict[str, str]) -> int:
    return sum(
        1 for k, cfg in _CRITERIA.items()
        if cfg.is_core and classes[k] == "poor"
    )


def _final_decision(weighted_score: float, core_poor: int) -> str:
    """
    PASS        WeightedScore >= 75  AND  CorePoorCount == 0
    BORDERLINE  55 <= score < 75  OR  (score >= 75 and CorePoorCount == 1)
    FAIL        score < 55  OR  CorePoorCount >= 2
    """
    if core_poor >= 2 or weighted_score < 55:
        return "FAIL"
    if weighted_score >= 75 and core_poor == 0:
        return "PASS"
    return "BORDERLINE"


# ---------------------------------------------------------------------------
# Main screening function
# ---------------------------------------------------------------------------

def screen_records(records: list[dict[str, Any]], ph: float = PH_SC) -> pd.DataFrame:
    """Apply follicle-oriented topical-delivery screening criteria to *records*.

    Parameters
    ----------
    records:
        List of molecule records from :func:`epidermal_barrier_screen.io.parse_input`.
    ph:
        Target pH for ionization and logD calculations.  Defaults to
        :data:`~epidermal_barrier_screen.ionization.PH_SC` (5.5).

    Returns
    -------
    pandas.DataFrame containing raw descriptors, legacy status columns,
    new classification / score columns, WeightedScore, CorePoorCount,
    FinalDecision, and a backward-compat ``final_result`` alias.
    """
    rows = []
    for rec in records:
        row: dict[str, Any] = {
            "name":             rec.get("name"),
            "input_smiles":     rec.get("input_smiles"),
            "canonical_smiles": rec.get("canonical_smiles"),
            "parse_status":     rec.get("parse_status", "invalid"),
        }

        if rec.get("parse_status") != "ok" or rec.get("mol") is None:
            row["FinalDecision"] = "invalid_input"
            row["final_result"]  = "invalid_input"
            rows.append(row)
            continue

        mol  = rec["mol"]
        desc = calculate(mol)
        row.update(desc)

        # ── PAINS and Brenk structural alerts ────────────────────────────────
        pains_hit = _check_pains(mol)
        brenk_hit = _check_brenk(mol)
        # Show "alert_name 🚩" when flagged, empty string when clean
        row["PAINS"]            = f"{pains_hit} \U0001f6a9" if pains_hit else ""
        row["Toxicity (BRENK)"] = f"{brenk_hit} \U0001f6a9" if brenk_hit else ""

        # ── pKa / ionization (unchanged calculation engine) ───────────────────
        smiles    = rec.get("canonical_smiles") or ""
        input_pka = rec.get("input_pka")

        if input_pka is not None:
            _, ion_type = predict_pka(smiles, name=rec.get("name", ""), ph=ph)
            if ion_type == "non_ionizable":
                ion_type = "acid"
            row["ionization_class"] = ion_type
            pka_val: float = input_pka
            if ion_type == "base":
                f_neutral, charge = _hhb_base(pka_val, ph)
            else:
                f_neutral, charge = _hhb_acid(pka_val, ph)
            row["predicted_pka"]       = round(pka_val, 2)
            # Store as percent (0–100); raw f_neutral kept for logD calc below
            row["fraction_unionized"]  = round(f_neutral * 100.0, 2)
            row["fraction_ionized"]    = round((1.0 - f_neutral) * 100.0, 2)
            row["expected_net_charge"] = round(charge, 4)
            row["logd"]        = round(desc["clogp"] + math.log10(max(f_neutral, 1e-10)), 4)
            row["logd_method"] = "pKa-corrected (input)"

        else:
            pka_val, ion_type = predict_pka(smiles, name=rec.get("name", ""), ph=ph)
            row["ionization_class"] = ion_type

            if ion_type == "non_ionizable" or pka_val is None:
                row["predicted_pka"]       = None
                row["fraction_unionized"]  = 100.0   # 100 % neutral
                row["fraction_ionized"]    = 0.0
                row["expected_net_charge"] = 0.0
                row["logd"]        = desc["clogp"]
                row["logd_method"] = "neutral (= cLogP)"
            else:
                if ion_type == "base":
                    f_neutral, charge = _hhb_base(pka_val, ph)
                else:
                    f_neutral, charge = _hhb_acid(pka_val, ph)
                row["predicted_pka"]       = round(pka_val, 2)
                # Store as percent; raw f_neutral used for logD calc
                row["fraction_unionized"]  = round(f_neutral * 100.0, 2)
                row["fraction_ionized"]    = round((1.0 - f_neutral) * 100.0, 2)
                row["expected_net_charge"] = round(charge, 4)
                row["logd"]        = round(desc["clogp"] + math.log10(max(f_neutral, 1e-10)), 4)
                row["logd_method"] = "pKa-corrected (pKaLearn)"

        # ── LogD override from SDF input ─────────────────────────────────────
        input_logd = rec.get("input_logd_7_4")
        if input_logd is not None:
            row["logd"]        = input_logd
            row["logd_method"] = "input (experimental)"

        # ── unionized %: percent of neutral species from logD and logP ──────────
        # logD is row["logd"] (pH-adjusted, fully resolved above)
        # logP is desc["clogp"] (RDKit cLogP, pH-independent)
        # Formula: unionized_pct = 10 ** (logD - logP) * 100   — percent, NOT capped
        row["unionized"] = _compute_unionized_pct(row.get("logd"), desc.get("clogp"))

        # ── Legacy status columns (backward-compat / table colour-coding) ─────
        # ionization_status uses the primary ionization value: unionized (%)
        # which is computed from the final logD (after any SDF override).
        row["mw_status"]            = _mw_status(desc["mw"])
        row["logd_status"]          = _logd_status(row["logd"])
        row["tpsa_status"]          = _tpsa_status(desc["tpsa"])
        row["hbd_status"]           = _hbd_status(desc["hbd"])
        row["hba_status"]           = _hba_status(desc["hba"])
        row["rotb_status"]          = _rotb_status(desc["rotb"])
        row["hac_status"]           = _hac_status(desc["hac"])        # informational
        row["formal_charge_status"] = _charge_status(desc["formal_charge"])
        row["ionization_status"]    = _ionization_status(row["unionized"])

        # ── New classification columns  (optimal / acceptable / poor) ─────────
        # UnionizedFraction is scored from row["unionized"] (logD-logP based, %)
        # because that column already reflects any experimental logD override.
        # When unionized is None (pathological: logD or logP absent), fall back
        # to fraction_unionized (HHB-based %) so scoring is never silently empty.
        _u_pct = row["unionized"] if row["unionized"] is not None else row.get("fraction_unionized")
        classes: dict[str, str] = {
            "MW":                _classify_mw(desc["mw"]),
            "LogD":              _classify_logd(row["logd"]),
            "TPSA":              _classify_tpsa(desc["tpsa"]),
            "FormalCharge":      _classify_formal_charge(desc["formal_charge"]),
            "UnionizedFraction": _classify_unionized(_u_pct),
            "HBD":               _classify_hbd(desc["hbd"]),
            "HBA":               _classify_hba(desc["hba"]),
            "RotB":              _classify_rotb(desc["rotb"]),
        }

        row["MW_class"]                = classes["MW"]
        row["LogD_class"]              = classes["LogD"]
        row["TPSA_class"]              = classes["TPSA"]
        row["FormalCharge_class"]      = classes["FormalCharge"]
        row["UnionizedFraction_class"] = classes["UnionizedFraction"]
        row["HBD_class"]               = classes["HBD"]
        row["HBA_class"]               = classes["HBA"]
        row["RotB_class"]              = classes["RotB"]

        # ── Per-criterion score contributions ─────────────────────────────────
        row["MW_score"]                = _criterion_score(classes["MW"],               _CRITERIA["MW"].weight)
        row["LogD_score"]              = _criterion_score(classes["LogD"],             _CRITERIA["LogD"].weight)
        row["TPSA_score"]              = _criterion_score(classes["TPSA"],             _CRITERIA["TPSA"].weight)
        row["FormalCharge_score"]      = _criterion_score(classes["FormalCharge"],     _CRITERIA["FormalCharge"].weight)
        row["UnionizedFraction_score"] = _criterion_score(classes["UnionizedFraction"],_CRITERIA["UnionizedFraction"].weight)
        row["HBD_score"]               = _criterion_score(classes["HBD"],              _CRITERIA["HBD"].weight)
        row["HBA_score"]               = _criterion_score(classes["HBA"],              _CRITERIA["HBA"].weight)
        row["RotB_score"]              = _criterion_score(classes["RotB"],              _CRITERIA["RotB"].weight)

        # ── Final weighted score and decision ─────────────────────────────────
        weighted_score = _compute_weighted_score(classes)
        core_poor      = _count_core_poor(classes)

        row["WeightedScore"]  = weighted_score
        row["CorePoorCount"]  = core_poor
        row["FinalDecision"]  = _final_decision(weighted_score, core_poor)
        row["final_result"]   = row["FinalDecision"]   # backward-compat alias

        rows.append(row)

    col_order = [
        "name",
        "parse_status",
        # ── Raw descriptors ──────────────────────────────────────────────────
        "mw",
        "tpsa",
        "hbd",
        "hba",
        "rotb",
        "hac",
        "predicted_pka",
        "ionization_class",
        "clogp",
        "logd",
        # unionized: continuous float [0,1] derived from 10**(logD-logP)
        # NOT a binary flag — represents the neutral species fraction at pH
        "unionized",
        "logd_method",
        "fraction_unionized",
        "fraction_ionized",
        "expected_net_charge",
        "formal_charge",
        # ── Legacy status columns ─────────────────────────────────────────────
        "mw_status",
        "logd_status",
        "tpsa_status",
        "hbd_status",
        "hba_status",
        "rotb_status",
        "hac_status",
        "formal_charge_status",
        "ionization_status",
        # ── New classification columns ────────────────────────────────────────
        "MW_class",
        "LogD_class",
        "TPSA_class",
        "FormalCharge_class",
        "UnionizedFraction_class",
        "HBD_class",
        "HBA_class",
        "RotB_class",
        # ── Per-criterion score contributions ─────────────────────────────────
        "MW_score",
        "LogD_score",
        "TPSA_score",
        "FormalCharge_score",
        "UnionizedFraction_score",
        "HBD_score",
        "HBA_score",
        "RotB_score",
        # ── Summary ───────────────────────────────────────────────────────────
        "WeightedScore",
        "CorePoorCount",
        "FinalDecision",
        "final_result",
        # ── Structural alerts ─────────────────────────────────────────────────
        "PAINS",
        "Toxicity (BRENK)",
        # ── SMILES ────────────────────────────────────────────────────────────
        "input_smiles",
        "canonical_smiles",
    ]

    df = pd.DataFrame(rows)
    for col in col_order:
        if col not in df.columns:
            df[col] = None
    # Ensure unionized column is numeric float dtype
    if "unionized" in df.columns:
        df["unionized"] = pd.to_numeric(df["unionized"], errors="coerce")
    return df[col_order]


# ---------------------------------------------------------------------------
# Self-test / validation for _compute_unionized
# Run:  python -m epidermal_barrier_screen.screen
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Self-test for _compute_unionized_pct
    # Values are NOW in PERCENT (multiply old expected by 100)
    _CASES = [
        # (logD,  logP,  expected_%,  description)
        (0.94,   1.50,  27.54,  "partial → 27.54 %"),
        (1.00,   1.20,  63.10,  "partial → 63.10 %"),
        (-1.00,  1.00,  1.00,   "mostly ionized → 1.00 %"),
        (1.60,   1.50,  125.89, "logD > logP → > 100 %, not capped"),
        (1.50,   1.50,  100.0,  "logD = logP → 100 %"),
        (None,   1.00,  None,   "missing logD → None"),
        (1.00,   None,  None,   "missing logP → None"),
        (float("nan"), 1.00, None, "NaN logD → None"),
    ]

    print("=" * 60)
    print("_compute_unionized_pct validation  (values in PERCENT)")
    print("=" * 60)
    all_ok = True
    for logd, logp, expected, desc in _CASES:
        result = _compute_unionized_pct(logd, logp)
        if expected is None:
            ok = result is None
        else:
            ok = result is not None and abs(result - expected) < 0.02
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_ok = False
        print(f"  [{status}]  logD={logd!r:>6}, logP={logp!r:>6}  →  {str(result):>10} %  (expected {expected} %)  | {desc}")

    print("=" * 60)
    print("All tests PASSED" if all_ok else "SOME TESTS FAILED")
    print("=" * 60)
