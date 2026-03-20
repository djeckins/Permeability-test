"""Screen molecule records against epidermal barrier passage criteria."""
from __future__ import annotations

import math
from typing import Any

import pandas as pd

from epidermal_barrier_screen.descriptors import calculate

# ---------------------------------------------------------------------------
# Criterion definitions
# Each criterion is (column, (opt_min, opt_max), (sub_min, sub_max))
# None means unbounded in that direction.
# ---------------------------------------------------------------------------
_CRITERIA: list[tuple[str, tuple, tuple]] = [
    # (descriptor_col, optimal_range, suboptimal_range)
    # ranges: (low_inclusive, high_exclusive) or None for no bound
    ("mw",           (None, 300),    (300,  500)),
    ("logd",         (1,    3),      (None, None)),   # handled specially below
    ("tpsa",         (None, 60),     (60,   130)),
    ("hbd",          (0,    3),      (3,    5)),       # 0-3 optimal, 4-5 suboptimal
    ("hba",          (2,    8),      (8,    10)),
    ("rotb",         (None, 10),     (10,   15)),
    ("hac",          (None, 30),     (30,   50)),
    ("formal_charge",(0,    0),      (-1,   1)),       # handled specially below
]

_FAIL_THRESHOLDS = {
    "mw":           500,
    "logd_lo":      0.5,  # below 0.5 is fail
    "logd_hi":      5,    # above 5 is fail
    "tpsa":         130,
    "hbd":          5,
    "hba":          10,
    "rotb":         15,
    "hac":          50,
}

_PH_SC = 5.5  # stratum corneum pH used for fraction-unionized estimate


def _fraction_unionized(pka: float, ph: float = _PH_SC) -> float:
    """Fraction in neutral (non-ionized) form assuming the pKa is for an acid.

    Uses Henderson-Hasselbalch: neutral_fraction = 1 / (1 + 10^(pH - pKa))
    """
    try:
        return 1.0 / (1.0 + 10 ** (ph - pka))
    except (OverflowError, ZeroDivisionError):
        return 0.0


def _mw_status(v: float) -> str:
    if v < 300:
        return "optimal"
    if v <= 500:
        return "suboptimal"
    return "fail"


def _logd_status(v: float) -> str:
    if 1.0 <= v <= 3.0:
        return "optimal"
    if 0.5 <= v < 1.0 or 3.0 < v <= 5.0:
        return "suboptimal"
    return "fail"


def _tpsa_status(v: float) -> str:
    if v < 60:
        return "optimal"
    if v <= 130:
        return "suboptimal"
    return "fail"


def _hbd_status(v: int) -> str:
    if 0 <= v <= 3:
        return "optimal"
    if 4 <= v <= 5:
        return "suboptimal"
    return "fail"


def _hba_status(v: int) -> str:
    if 2 <= v <= 8:
        return "optimal"
    if 8 < v <= 10:
        return "suboptimal"
    return "fail"


def _rotb_status(v: int) -> str:
    if v < 10:
        return "optimal"
    if v <= 15:
        return "suboptimal"
    return "fail"


def _hac_status(v: int) -> str:
    if v < 30:
        return "optimal"
    if v <= 50:
        return "suboptimal"
    return "fail"


def _charge_status(v: int) -> str:
    if v == 0:
        return "optimal"
    if v in (-1, 1):
        return "suboptimal"
    return "fail"


def _final_result(statuses: list[str]) -> str:
    fail_count = statuses.count("fail")
    sub_count = statuses.count("suboptimal")
    if fail_count >= 2:
        return "enhancers_likely_needed"
    if fail_count == 1:
        return "borderline"
    if sub_count > 0:
        return "pass_with_manual_review"
    return "pass"


def screen_records(records: list[dict[str, Any]]) -> pd.DataFrame:
    """Apply epidermal barrier criteria to *records* and return a results DataFrame.

    Parameters
    ----------
    records:
        List of molecule records as produced by :func:`epidermal_barrier_screen.io.parse_input`.

    Returns
    -------
    pandas.DataFrame with descriptor columns, per-criterion status columns, and
    a ``final_result`` column.
    """
    rows = []
    for rec in records:
        row: dict[str, Any] = {
            "name": rec.get("name"),
            "input_smiles": rec.get("input_smiles"),
            "canonical_smiles": rec.get("canonical_smiles"),
            "parse_status": rec.get("parse_status", "invalid"),
            "input_pka": rec.get("input_pka"),
            "input_logd_7_4": rec.get("input_logd_7_4"),
        }

        if rec.get("parse_status") != "ok" or rec.get("mol") is None:
            row["final_result"] = "invalid_input"
            rows.append(row)
            continue

        mol = rec["mol"]
        desc = calculate(mol)
        row.update(desc)

        # Determine the logD value to use
        input_logd = rec.get("input_logd_7_4")
        if input_logd is not None:
            row["logd"] = input_logd
            row["logd_source"] = "input_logd_7_4"
        else:
            row["logd"] = desc["clogp"]
            row["logd_source"] = "clogp_proxy"

        # Fraction unionized at pH 5.5
        pka = rec.get("input_pka")
        row["fraction_unionized_pH5_5"] = (
            round(_fraction_unionized(pka), 4) if pka is not None else None
        )

        # Per-criterion statuses
        row["mw_status"] = _mw_status(desc["mw"])
        row["logd_status"] = _logd_status(row["logd"])
        row["tpsa_status"] = _tpsa_status(desc["tpsa"])
        row["hbd_status"] = _hbd_status(desc["hbd"])
        row["hba_status"] = _hba_status(desc["hba"])
        row["rotb_status"] = _rotb_status(desc["rotb"])
        row["hac_status"] = _hac_status(desc["hac"])
        row["formal_charge_status"] = _charge_status(desc["formal_charge"])

        statuses = [
            row["mw_status"],
            row["logd_status"],
            row["tpsa_status"],
            row["hbd_status"],
            row["hba_status"],
            row["rotb_status"],
            row["hac_status"],
            row["formal_charge_status"],
        ]
        row["final_result"] = _final_result(statuses)
        rows.append(row)

    col_order = [
        "name",
        "input_smiles",
        "canonical_smiles",
        "parse_status",
        "mw",
        "clogp",
        "logd",
        "logd_source",
        "tpsa",
        "hbd",
        "hba",
        "rotb",
        "hac",
        "formal_charge",
        "input_pka",
        "input_logd_7_4",
        "fraction_unionized_pH5_5",
        "mw_status",
        "logd_status",
        "tpsa_status",
        "hbd_status",
        "hba_status",
        "rotb_status",
        "hac_status",
        "formal_charge_status",
        "final_result",
    ]

    df = pd.DataFrame(rows)
    # Reorder columns; add missing ones as NaN
    for col in col_order:
        if col not in df.columns:
            df[col] = None
    return df[col_order]
