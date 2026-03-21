"""Ionization state prediction.

Detects ionizable functional groups via SMARTS patterns, predicts pKa
using the pKaPredict LGBMRegressor model, and computes ionization state
via Henderson-Hasselbalch.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from rdkit import Chem

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PH_SC: float = 5.5  # stratum corneum surface pH

IonType = Literal["acid", "base"]


# ---------------------------------------------------------------------------
# SMARTS-based ionizable group detection
# ---------------------------------------------------------------------------

_ACID_SMARTS: list[tuple[str, str]] = [
    ("[CX3](=O)[OX2H1]", "carboxylic_acid"),
    ("[c][OX2H1]",        "phenol"),
    ("[SX2H1]",           "thiol"),
    ("[NH1]S(=O)(=O)",    "sulfonamide"),
    ("[NH2]S(=O)(=O)",    "sulfonamide"),
    ("[PX4](=O)[OX2H1]", "phosphate"),
]

_BASE_SMARTS: list[tuple[str, str]] = [
    ("[NX3H2;!$(NC=O);!$(NS(=O)(=O))]",          "primary_amine"),
    ("[NX3H1;!$(NC=O);!$(NS(=O)(=O));!$([nH])]",  "secondary_amine"),
    ("[NX3H0;!$(NC=O);!$(NS(=O)(=O))]",           "tertiary_amine"),
    ("[nX2H0]",                                     "pyridine_like"),
    ("[$([CX3](=[NX2])([NX3])[NX3])]",            "guanidine"),
    ("[$([CX3](=[NX2])[NX3])]",                   "amidine"),
]


@dataclass
class DetectedGroup:
    """An ionizable functional group detected via SMARTS."""
    name: str
    ion_type: IonType


def detect_ionizable_groups(mol: Chem.Mol) -> list[DetectedGroup]:
    """Detect ionizable functional groups in *mol* using SMARTS patterns.

    Returns a deduplicated list of :class:`DetectedGroup` instances.
    """
    groups: list[DetectedGroup] = []
    claimed: set[int] = set()

    for patterns, ion_type in [(_ACID_SMARTS, "acid"), (_BASE_SMARTS, "base")]:
        for smarts_str, name in patterns:
            pattern = Chem.MolFromSmarts(smarts_str)
            if pattern is None:
                continue
            for match in mol.GetSubstructMatches(pattern):
                key = match[0]
                if key in claimed:
                    continue
                claimed.add(key)
                groups.append(DetectedGroup(name=name, ion_type=ion_type))

    return groups


def classify_ionization(groups: list[DetectedGroup]) -> str:
    """Classify a molecule as non_ionizable / acid / base / ampholyte."""
    if not groups:
        return "non_ionizable"
    has_acid = any(g.ion_type == "acid" for g in groups)
    has_base = any(g.ion_type == "base" for g in groups)
    if has_acid and has_base:
        return "ampholyte"
    if has_acid:
        return "acid"
    return "base"


# ---------------------------------------------------------------------------
# Henderson-Hasselbalch helpers
# ---------------------------------------------------------------------------


def _hhb_acid(pka: float, ph: float) -> tuple[float, float]:
    """Henderson-Hasselbalch for an *acid* group (loses proton).

    Returns (fraction_neutral, mean_charge).
    """
    try:
        ratio = 10 ** (ph - pka)
        f_neutral = 1.0 / (1.0 + ratio)
        mean_charge = -ratio / (1.0 + ratio)
    except OverflowError:
        f_neutral, mean_charge = 0.0, -1.0
    return round(f_neutral, 6), round(mean_charge, 6)


def _hhb_base(pka: float, ph: float) -> tuple[float, float]:
    """Henderson-Hasselbalch for a *base* group (gains proton).

    Returns (fraction_neutral, mean_charge).
    """
    try:
        ratio = 10 ** (pka - ph)
        f_neutral = 1.0 / (1.0 + ratio)
        mean_charge = ratio / (1.0 + ratio)
    except OverflowError:
        f_neutral, mean_charge = 0.0, 1.0
    return round(f_neutral, 6), round(mean_charge, 6)


# ---------------------------------------------------------------------------
# pKaPredict integration  (molecule-specific ML pKa)
# ---------------------------------------------------------------------------


def _get_pka_model():
    """Lazily load (and cache) the pKaPredict LGBMRegressor model."""
    if _get_pka_model._cache is None:
        from pkapredict import load_model  # type: ignore[import]

        model = load_model()
        _get_pka_model._cache = (model, model.feature_name_)
    return _get_pka_model._cache


_get_pka_model._cache = None  # type: ignore[attr-defined]


def _ml_pka(smiles: str) -> float | None:
    """Return a molecule-specific pKa prediction from pKaPredict, or *None* on failure."""
    import sys

    try:
        from pkapredict import predict_pKa  # type: ignore[import]

        model, feat = _get_pka_model()
        val = float(predict_pKa(smiles, model, feat))
        return val
    except Exception as exc:
        print(f"[pKaPredict] WARNING: prediction failed for {smiles!r}: {exc}", file=sys.stderr)
        return None
