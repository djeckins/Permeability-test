"""Screen molecule records against epidermal barrier passage criteria."""
from __future__ import annotations

import math
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
# Per-criterion status functions
# Each returns one of: "optimal", "suboptimal", "poor"
# ---------------------------------------------------------------------------


def _mw_status(v: float) -> str:
    if v < 300:
        return "optimal"
    if v <= 500:
        return "suboptimal"
    return "poor"


def _logd_status(v: float | None) -> str:
    if v is None:
        return "poor"
    if 1.0 <= v <= 3.0:
        return "optimal"
    if 0.5 <= v < 1.0 or 3.0 < v <= 5.0:
        return "suboptimal"
    return "poor"


def _tpsa_status(v: float) -> str:
    if v < 60:
        return "optimal"
    if v <= 130:
        return "suboptimal"
    return "poor"


def _hbd_status(v: int) -> str:
    if 0 <= v <= 3:
        return "optimal"
    if 4 <= v <= 5:
        return "suboptimal"
    return "poor"


def _hba_status(v: int) -> str:
    if 2 <= v <= 8:
        return "optimal"
    if 8 < v <= 10:
        return "suboptimal"
    return "poor"


def _rotb_status(v: int) -> str:
    if v < 10:
        return "optimal"
    if v <= 15:
        return "suboptimal"
    return "poor"


def _hac_status(v: int) -> str:
    if v < 30:
        return "optimal"
    if v <= 50:
        return "suboptimal"
    return "poor"


def _charge_status(v: int) -> str:
    if v == 0:
        return "optimal"
    if v in (-1, 1):
        return "suboptimal"
    return "poor"


def _ionization_status(fraction_unionized: float | None) -> str:
    if fraction_unionized is None:
        return "poor"
    if fraction_unionized >= 0.8:
        return "optimal"
    if fraction_unionized >= 0.5:
        return "suboptimal"
    return "poor"


# ---------------------------------------------------------------------------
# Overall result
# ---------------------------------------------------------------------------


def _final_result(statuses: list[str]) -> str:
    """Derive overall PASS / BORDERLINE / FAIL from per-criterion statuses.

    PASS       – at most 1 non-optimal criterion (suboptimal or poor)
    BORDERLINE – multiple suboptimal with no poor, or up to 2 poor
    FAIL       – 3 or more poor criteria
    """
    n_poor = statuses.count("poor")
    n_suboptimal = statuses.count("suboptimal")
    n_non_optimal = n_poor + n_suboptimal

    if n_non_optimal <= 1:
        return "PASS"
    if n_poor <= 2:
        return "BORDERLINE"
    return "FAIL"


def screen_records(records: list[dict[str, Any]], ph: float = PH_SC) -> pd.DataFrame:
    """Apply epidermal barrier criteria to *records* and return a results DataFrame.

    Parameters
    ----------
    records:
        List of molecule records as produced by :func:`epidermal_barrier_screen.io.parse_input`.
    ph:
        Target pH for ionization and logD calculations.  Defaults to
        :data:`~epidermal_barrier_screen.ionization.PH_SC` (5.5, stratum corneum).

    Returns
    -------
    pandas.DataFrame with descriptor columns, per-criterion status columns, and
    a ``final_result`` column (PASS / BORDERLINE / FAIL).
    """
    rows = []
    for rec in records:
        row: dict[str, Any] = {
            "name": rec.get("name"),
            "input_smiles": rec.get("input_smiles"),
            "canonical_smiles": rec.get("canonical_smiles"),
            "parse_status": rec.get("parse_status", "invalid"),
        }

        if rec.get("parse_status") != "ok" or rec.get("mol") is None:
            row["final_result"] = "invalid_input"
            rows.append(row)
            continue

        mol = rec["mol"]
        desc = calculate(mol)
        row.update(desc)

        # ── pKa prediction via pKaLearn GNN ──────────────────────────────────
        smiles = rec.get("canonical_smiles") or ""
        input_pka = rec.get("input_pka")

        if input_pka is not None:
            # User-provided experimental pKa — use pKaLearn only for ion-type detection
            _, ion_type = predict_pka(smiles, name=rec.get("name", ""), ph=ph)
            if ion_type == "non_ionizable":
                ion_type = "acid"  # conservative fallback when type is ambiguous
            row["ionization_class"] = ion_type
            pka_val: float = input_pka
            if ion_type == "base":
                f_neutral, charge = _hhb_base(pka_val, ph)
            else:
                f_neutral, charge = _hhb_acid(pka_val, ph)

            row["predicted_pka"] = round(pka_val, 2)
            row["fraction_unionized"] = round(f_neutral, 4)
            row["fraction_ionized"] = round(1.0 - f_neutral, 4)
            row["expected_net_charge"] = round(charge, 4)
            row["logd"] = round(desc["clogp"] + math.log10(max(f_neutral, 1e-10)), 4)
            row["logd_method"] = "pKa-corrected (input)"

        else:
            # Predict pKa and ion type with pKaLearn
            pka_val, ion_type = predict_pka(smiles, name=rec.get("name", ""), ph=ph)
            row["ionization_class"] = ion_type

            if ion_type == "non_ionizable" or pka_val is None:
                # Non-ionizable — neutral form at all relevant pH values
                row["predicted_pka"] = None
                row["fraction_unionized"] = 1.0
                row["fraction_ionized"] = 0.0
                row["expected_net_charge"] = 0.0
                row["logd"] = desc["clogp"]
                row["logd_method"] = "neutral (= cLogP)"
            else:
                if ion_type == "base":
                    f_neutral, charge = _hhb_base(pka_val, ph)
                else:
                    f_neutral, charge = _hhb_acid(pka_val, ph)

                row["predicted_pka"] = round(pka_val, 2)
                row["fraction_unionized"] = round(f_neutral, 4)
                row["fraction_ionized"] = round(1.0 - f_neutral, 4)
                row["expected_net_charge"] = round(charge, 4)
                row["logd"] = round(desc["clogp"] + math.log10(max(f_neutral, 1e-10)), 4)
                row["logd_method"] = "pKa-corrected (pKaLearn)"

        # ── logD override from input ─────────────────────────────────────────
        input_logd = rec.get("input_logd_7_4")
        if input_logd is not None:
            row["logd"] = input_logd
            row["logd_method"] = "input (experimental)"

        # ── Per-criterion statuses ───────────────────────────────────────────
        row["mw_status"] = _mw_status(desc["mw"])
        row["logd_status"] = _logd_status(row["logd"])
        row["tpsa_status"] = _tpsa_status(desc["tpsa"])
        row["hbd_status"] = _hbd_status(desc["hbd"])
        row["hba_status"] = _hba_status(desc["hba"])
        row["rotb_status"] = _rotb_status(desc["rotb"])
        row["hac_status"] = _hac_status(desc["hac"])
        row["formal_charge_status"] = _charge_status(desc["formal_charge"])
        row["ionization_status"] = _ionization_status(row["fraction_unionized"])

        statuses = [
            row["mw_status"],
            row["logd_status"],
            row["tpsa_status"],
            row["hbd_status"],
            row["hba_status"],
            row["rotb_status"],
            row["hac_status"],
            row["formal_charge_status"],
            row["ionization_status"],
        ]
        row["final_result"] = _final_result(statuses)
        rows.append(row)

    col_order = [
        "name",
        "parse_status",
        "mw",
        "hba",
        "hbd",
        "rotb",
        "hac",
        "predicted_pka",
        "ionization_class",
        "clogp",
        "logd",
        "logd_method",
        "tpsa",
        "fraction_unionized",
        "fraction_ionized",
        "expected_net_charge",
        "formal_charge",
        # ── Per-criterion statuses ───────────────────────────────────────────
        "mw_status",
        "logd_status",
        "tpsa_status",
        "hbd_status",
        "hba_status",
        "rotb_status",
        "hac_status",
        "formal_charge_status",
        "ionization_status",
        "final_result",
        # ── SMILES (moved to end) ────────────────────────────────────────────
        "input_smiles",
        "canonical_smiles",
    ]

    df = pd.DataFrame(rows)
    for col in col_order:
        if col not in df.columns:
            df[col] = None
    return df[col_order]
