"""Screen molecule records against epidermal barrier passage criteria."""
from __future__ import annotations

import math
from typing import Any

import pandas as pd

from epidermal_barrier_screen.descriptors import calculate
from epidermal_barrier_screen.ionization import PH_SC, ionize

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


def _logd_status(v: float) -> str:
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


def _ionization_status(fraction_unionized: float) -> str:
    """Classify ionization based on fraction unionized at pH 5.5.

    A high fraction unionized means the molecule is mostly neutral at the
    stratum corneum surface pH, which is favourable for passive permeation.

    - optimal:   >= 0.8  (mostly neutral)
    - suboptimal: 0.5 – 0.8  (partially ionized)
    - poor:       < 0.5  (majority ionized)
    """
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

    PASS       – all criteria optimal (no flags)
    BORDERLINE – one or more suboptimal, but no poor
    FAIL       – one or more poor (strong failure in a core parameter)
    """
    if "poor" in statuses:
        return "FAIL"
    if "suboptimal" in statuses:
        return "BORDERLINE"
    return "PASS"


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

        # ── Ionization at user-specified pH ─────────────────────────────────
        ion = ionize(mol, ph=ph)

        input_pka = rec.get("input_pka")
        if input_pka is not None:
            # Use the experimentally supplied pKa; acid vs base from predicted ionization
            ion_type = ion.dominant_type if ion.dominant_group is not None else "acid"
            if ion_type == "base":
                from epidermal_barrier_screen.ionization import _hhb_base
                f_neutral, _ = _hhb_base(input_pka, ph)
            else:
                from epidermal_barrier_screen.ionization import _hhb_acid
                f_neutral, _ = _hhb_acid(input_pka, ph)
            row["pka_source"] = "input"
            row["predicted_pka"] = None
            row["predicted_pka_type"] = ion_type
            row["fraction_unionized_pH5_5"] = round(f_neutral, 4)
        elif ion.dominant_group is not None:
            row["pka_source"] = "predicted"
            row["predicted_pka"] = ion.dominant_pka
            row["predicted_pka_type"] = ion.dominant_type
            row["fraction_unionized_pH5_5"] = round(ion.fraction_neutral_total, 4)
        else:
            row["pka_source"] = "none"
            row["predicted_pka"] = None
            row["predicted_pka_type"] = None
            row["fraction_unionized_pH5_5"] = 1.0  # non-ionizable → always neutral

        # ── logD at user-specified pH ────────────────────────────────────────
        # Priority: experimentally supplied logD > pH-corrected clogP
        input_logd = rec.get("input_logd_7_4")
        if input_logd is not None:
            row["logd"] = input_logd
            row["logd_source"] = "input_logd_7_4"
        else:
            f_total = row["fraction_unionized_pH5_5"]
            # logD(pH) = clogP + log10(fraction_neutral)
            logd_correction = math.log10(max(f_total, 1e-10))
            row["logd"] = round(desc["clogp"] + logd_correction, 4)
            row["logd_source"] = f"clogp_pH{ph:.1f}"

        row["mean_charge_pH5_5"] = round(ion.mean_charge, 4)
        row["ionization_class"] = ion.ionization_class

        # ── Per-criterion statuses ───────────────────────────────────────────
        row["mw_status"] = _mw_status(desc["mw"])
        row["logd_status"] = _logd_status(row["logd"])
        row["tpsa_status"] = _tpsa_status(desc["tpsa"])
        row["hbd_status"] = _hbd_status(desc["hbd"])
        row["hba_status"] = _hba_status(desc["hba"])
        row["rotb_status"] = _rotb_status(desc["rotb"])
        row["hac_status"] = _hac_status(desc["hac"])
        row["formal_charge_status"] = _charge_status(desc["formal_charge"])
        row["ionization_status"] = _ionization_status(row["fraction_unionized_pH5_5"])

        statuses = [
            row["mw_status"],
            row["logd_status"],
            row["tpsa_status"],
            row["hbd_status"],
            row["hba_status"],
            row["rotb_status"],
            row["hac_status"],
            row["formal_charge_status"],
            row["ionization_status"],  # mandatory criterion
        ]
        row["final_result"] = _final_result(statuses)
        rows.append(row)

    col_order = [
        "name",
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
        # ── pKa / ionization ────────────────────────────────────────────────
        "input_pka",
        "input_logd_7_4",
        "pka_source",
        "predicted_pka",
        "predicted_pka_type",
        "fraction_unionized_pH5_5",
        "mean_charge_pH5_5",
        "ionization_class",
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
        # ── SMILES (wide columns, moved to end) ─────────────────────────────
        "input_smiles",
        "canonical_smiles",
    ]

    df = pd.DataFrame(rows)
    # Reorder columns; add missing ones as NaN
    for col in col_order:
        if col not in df.columns:
            df[col] = None
    return df[col_order]
