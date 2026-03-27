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
# Scoring configuration
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class _CriterionCfg:
    weight: float   # contribution to the 95-point maximum raw score
    is_core: bool   # True → counted in CorePoorCount


# Edit weights / thresholds here; do NOT touch screen_records().
_CRITERIA: dict[str, _CriterionCfg] = {
    "MW":                _CriterionCfg(weight=20, is_core=True),
    "LogD":              _CriterionCfg(weight=20, is_core=True),
    "TPSA":              _CriterionCfg(weight=15, is_core=True),
    "FormalCharge":      _CriterionCfg(weight=15, is_core=True),
    "UnionizedFraction": _CriterionCfg(weight=15, is_core=True),
    "HBD":               _CriterionCfg(weight=5,  is_core=False),
    "RotB":              _CriterionCfg(weight=5,  is_core=False),
}

_MAX_RAW_SCORE: float = sum(c.weight for c in _CRITERIA.values())  # 95


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
    if v is None:
        return "poor"
    if v >= 0.40:
        return "optimal"
    if v >= 0.10:
        return "acceptable"
    return "poor"


def _classify_hbd(v: int) -> str:
    if v <= 2:
        return "optimal"
    if v == 3:
        return "acceptable"
    return "poor"


def _classify_rotb(v: int) -> str:
    if v <= 8:
        return "optimal"
    if v <= 12:
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
    """Informational only — not used in weighted scoring."""
    if 2 <= v <= 8:    return "optimal"
    if 8 < v <= 10:    return "suboptimal"
    return "poor"


def _rotb_status(v: int) -> str:
    if v <= 8:     return "optimal"
    if v <= 12:    return "suboptimal"
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


def _ionization_status(fraction_unionized: float | None) -> str:
    if fraction_unionized is None:   return "poor"
    if fraction_unionized >= 0.40:   return "optimal"
    if fraction_unionized >= 0.10:   return "suboptimal"
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
    """Return the normalised weighted score on a 0–100 scale (1 d.p.)."""
    raw = sum(_criterion_score(classes[k], _CRITERIA[k].weight) for k in _CRITERIA)
    return round((raw / _MAX_RAW_SCORE) * 100, 1)


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
            row["fraction_unionized"]  = round(f_neutral, 4)
            row["fraction_ionized"]    = round(1.0 - f_neutral, 4)
            row["expected_net_charge"] = round(charge, 4)
            row["logd"]        = round(desc["clogp"] + math.log10(max(f_neutral, 1e-10)), 4)
            row["logd_method"] = "pKa-corrected (input)"

        else:
            pka_val, ion_type = predict_pka(smiles, name=rec.get("name", ""), ph=ph)
            row["ionization_class"] = ion_type

            if ion_type == "non_ionizable" or pka_val is None:
                row["predicted_pka"]       = None
                row["fraction_unionized"]  = 1.0
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
                row["fraction_unionized"]  = round(f_neutral, 4)
                row["fraction_ionized"]    = round(1.0 - f_neutral, 4)
                row["expected_net_charge"] = round(charge, 4)
                row["logd"]        = round(desc["clogp"] + math.log10(max(f_neutral, 1e-10)), 4)
                row["logd_method"] = "pKa-corrected (pKaLearn)"

        # ── LogD override from SDF input ─────────────────────────────────────
        input_logd = rec.get("input_logd_7_4")
        if input_logd is not None:
            row["logd"]        = input_logd
            row["logd_method"] = "input (experimental)"

        # ── Legacy status columns (backward-compat / table colour-coding) ─────
        row["mw_status"]            = _mw_status(desc["mw"])
        row["logd_status"]          = _logd_status(row["logd"])
        row["tpsa_status"]          = _tpsa_status(desc["tpsa"])
        row["hbd_status"]           = _hbd_status(desc["hbd"])
        row["hba_status"]           = _hba_status(desc["hba"])        # informational
        row["rotb_status"]          = _rotb_status(desc["rotb"])
        row["hac_status"]           = _hac_status(desc["hac"])        # informational
        row["formal_charge_status"] = _charge_status(desc["formal_charge"])
        row["ionization_status"]    = _ionization_status(row["fraction_unionized"])

        # ── New classification columns  (optimal / acceptable / poor) ─────────
        classes: dict[str, str] = {
            "MW":                _classify_mw(desc["mw"]),
            "LogD":              _classify_logd(row["logd"]),
            "TPSA":              _classify_tpsa(desc["tpsa"]),
            "FormalCharge":      _classify_formal_charge(desc["formal_charge"]),
            "UnionizedFraction": _classify_unionized(row["fraction_unionized"]),
            "HBD":               _classify_hbd(desc["hbd"]),
            "RotB":              _classify_rotb(desc["rotb"]),
        }

        row["MW_class"]                = classes["MW"]
        row["LogD_class"]              = classes["LogD"]
        row["TPSA_class"]              = classes["TPSA"]
        row["FormalCharge_class"]      = classes["FormalCharge"]
        row["UnionizedFraction_class"] = classes["UnionizedFraction"]
        row["HBD_class"]               = classes["HBD"]
        row["RotB_class"]              = classes["RotB"]

        # ── Per-criterion score contributions ─────────────────────────────────
        row["MW_score"]                = _criterion_score(classes["MW"],               _CRITERIA["MW"].weight)
        row["LogD_score"]              = _criterion_score(classes["LogD"],             _CRITERIA["LogD"].weight)
        row["TPSA_score"]              = _criterion_score(classes["TPSA"],             _CRITERIA["TPSA"].weight)
        row["FormalCharge_score"]      = _criterion_score(classes["FormalCharge"],     _CRITERIA["FormalCharge"].weight)
        row["UnionizedFraction_score"] = _criterion_score(classes["UnionizedFraction"],_CRITERIA["UnionizedFraction"].weight)
        row["HBD_score"]               = _criterion_score(classes["HBD"],              _CRITERIA["HBD"].weight)
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
        "RotB_class",
        # ── Per-criterion score contributions ─────────────────────────────────
        "MW_score",
        "LogD_score",
        "TPSA_score",
        "FormalCharge_score",
        "UnionizedFraction_score",
        "HBD_score",
        "RotB_score",
        # ── Summary ───────────────────────────────────────────────────────────
        "WeightedScore",
        "CorePoorCount",
        "FinalDecision",
        "final_result",
        # ── SMILES ────────────────────────────────────────────────────────────
        "input_smiles",
        "canonical_smiles",
    ]

    df = pd.DataFrame(rows)
    for col in col_order:
        if col not in df.columns:
            df[col] = None
    return df[col_order]
