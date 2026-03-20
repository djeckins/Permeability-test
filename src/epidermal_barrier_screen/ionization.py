"""Rule-based pKa estimation and ionization state prediction.

Predicts the protonation state of a molecule at a target pH (default 5.5,
stratum corneum surface pH) using SMARTS-based functional group matching and
the Henderson-Hasselbalch equation.

Limitations
-----------
This is a *rule-based* estimator, not a quantum-mechanical or ML model.
pKa values are typical literature means for isolated functional groups;
electronic effects, ring strain, intramolecular hydrogen bonds, and
conjugation can shift real pKas by ±2 units.  For high-accuracy work,
validated computational tools (e.g. Epik, Jaguar pKa, pkasolver) should
be used.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

from rdkit import Chem

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PH_SC: float = 5.5  # stratum corneum surface pH

IonType = Literal["acid", "base"]

# ---------------------------------------------------------------------------
# SMARTS rules
# ---------------------------------------------------------------------------
# Each entry: (SMARTS_pattern, typical_pKa, ion_type, short_description)
#
# 'acid' groups lose a proton at pH > pKa (become negatively charged).
# 'base' groups gain a proton at pH < pKa (become positively charged).
#
# Rules are ordered from most specific to least specific.  During matching
# the first rule that claims a given atom wins, so more specific patterns
# must come first.
#
# SMARTS atom notation used:
#   [OH]  – oxygen with exactly one H (hydroxyl)
#   [SH]  – sulfur with one H (thiol)
#   [NH1] – nitrogen with exactly one H
#   [NH2] – nitrogen with two H's
#   X3    – 3 total bonds (valence)
#   X4    – 4 total bonds
#   $()   – recursive SMARTS (substructure requirement)
#   c     – aromatic carbon  /  n – aromatic nitrogen
# ---------------------------------------------------------------------------
_RAW_RULES: list[tuple[str, float, IonType, str]] = [
    # ── Strong acids ────────────────────────────────────────────────────────
    ("[OH][SX4](=O)(=O)",          -1.0, "acid", "sulfonic_acid"),
    ("[OH][PX4](=O)",               2.1, "acid", "phosphoric_acid"),

    # ── Carboxylic acids ─────────────────────────────────────────────────────
    ("[OH][CX3](=O)",               4.5, "acid", "carboxylic_acid"),

    # ── Aryl thiol ───────────────────────────────────────────────────────────
    ("[SH][cX3]",                   6.5, "acid", "aryl_thiol"),

    # ── Hydroxamic acid  (-C(=O)-N-OH) ────────────────────────────────────────
    ("[OH][NX3][CX3](=O)",          8.5, "acid", "hydroxamic_acid"),

    # ── Imide N-H  (-C(=O)-NH-C(=O)-) ────────────────────────────────────────
    ("[NH1]([CX3]=O)[CX3]=O",       9.6, "acid", "imide"),

    # ── Sulfonamide N-H ────────────────────────────────────────────────────────
    ("[NH1][SX4](=O)(=O)",         10.0, "acid", "sulfonamide"),

    # ── Phenol / enol OH on aromatic ring ─────────────────────────────────────
    ("[OH][cX3]",                  10.0, "acid", "phenol"),

    # ── Aliphatic thiol ────────────────────────────────────────────────────────
    ("[SH][CX4]",                  10.5, "acid", "aliphatic_thiol"),

    # ── Guanidine  (most basic; must precede amidine and amine rules) ───────────
    (
        "[NX3;H0,H1,H2][CX3](=[NX2;H0,H1])[NX3;H0,H1,H2]",
        12.5, "base", "guanidine",
    ),

    # ── Amidine  (C(=N)N, not guanidine) ──────────────────────────────────────
    (
        "[NX3;H1,H2;!$(NC(=N)N)][CX3;!$(CC=N)](=[NX2])",
        11.0, "base", "amidine",
    ),

    # ── Primary aliphatic amine ─────────────────────────────────────────────────
    (
        "[NX3;H2;"
        "!$(NC=O);!$(NS(=O)=O);!$(Nc);"
        "!$(N=*);!$(N#*)]",
        10.5, "base", "primary_aliphatic_amine",
    ),

    # ── Secondary aliphatic amine ───────────────────────────────────────────────
    (
        "[NX3;H1;"
        "!$(NC=O);!$(NS(=O)=O);!$(Nc);"
        "!$(N=*);!$(N#*)]",
        10.0, "base", "secondary_aliphatic_amine",
    ),

    # ── Tertiary aliphatic amine ────────────────────────────────────────────────
    (
        "[NX3;H0;"
        "!$(NC=O);!$(NS(=O)=O);!$(Nc);"
        "!$(N=*);!$(N#*)]",
        9.0, "base", "tertiary_aliphatic_amine",
    ),

    # ── Aromatic primary amine ──────────────────────────────────────────────────
    ("[NH2][cX3]",                  4.5, "base", "aromatic_primary_amine"),

    # ── Aromatic secondary amine (not sulfonamide, not amide) ───────────────────
    (
        "[NH1;!$(NS(=O)=O);!$(NC=O)][cX3]",
        3.5, "base", "aromatic_secondary_amine",
    ),

    # ── Pyridine-type N in 6-membered aromatic ring ──────────────────────────────
    ("[nH0;r6]",                    5.2, "base", "pyridine_N"),

    # ── Imidazole / benzimidazole basic N in 5-membered ring adjacent to [nH] ──
    (
        "[nH0;r5;$([nH0]1cc[nH]c1),$([nH0]1cnc2ccccc12),$([nH0]1cnc2cccc12)]",
        6.8, "base", "imidazole_N",
    ),

    # ── Other 5-membered aromatic N (oxazole, thiazole, pyrazole …) ────────────
    ("[nH0;r5]",                    2.5, "base", "heteroaromatic_5ring_N"),
]


# ---------------------------------------------------------------------------
# Pre-compile SMARTS patterns once at import time
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class _Rule:
    pattern: Chem.Mol
    pka: float
    ion_type: IonType
    description: str


def _compile_rules() -> list[_Rule]:
    compiled: list[_Rule] = []
    for smarts, pka, ion_type, desc in _RAW_RULES:
        pat = Chem.MolFromSmarts(smarts)
        if pat is None:
            raise ValueError(f"Invalid SMARTS in ionization rules: {smarts!r}")
        compiled.append(_Rule(pat, pka, ion_type, desc))
    return compiled


_RULES: list[_Rule] = _compile_rules()


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
@dataclass
class IonizableGroup:
    """One ionizable functional group found in a molecule."""

    description: str
    pka: float
    ion_type: IonType
    #: fraction of this group in its neutral form at the target pH
    fraction_neutral: float
    #: mean charge contribution of this group at the target pH
    mean_charge: float


@dataclass
class IonizationResult:
    """Full ionization analysis for one molecule at a given pH."""

    ph: float
    #: all ionizable groups detected
    groups: list[IonizableGroup] = field(default_factory=list)

    # ── Summary properties ─────────────────────────────────────────────────

    @property
    def non_ionizable(self) -> bool:
        return not self.groups

    @property
    def mean_charge(self) -> float:
        """Net mean charge at *ph* summed across all ionizable groups."""
        return sum(g.mean_charge for g in self.groups)

    @property
    def dominant_group(self) -> IonizableGroup | None:
        """Group whose pKa is closest to the target pH (most influential)."""
        if not self.groups:
            return None
        return min(self.groups, key=lambda g: abs(g.pka - self.ph))

    @property
    def dominant_pka(self) -> float | None:
        g = self.dominant_group
        return g.pka if g else None

    @property
    def dominant_type(self) -> IonType | None:
        g = self.dominant_group
        return g.ion_type if g else None

    @property
    def fraction_neutral_dominant(self) -> float | None:
        """Fraction neutral of the dominant (closest-to-pH) group."""
        g = self.dominant_group
        return g.fraction_neutral if g else None

    @property
    def ionization_class(self) -> str:
        """Broad ionization class at the target pH.

        Returns one of:
        ``'non_ionizable'``, ``'neutral'``, ``'acid'``, ``'base'``,
        ``'zwitterion'``.
        """
        if self.non_ionizable:
            return "non_ionizable"

        acid_charge = sum(g.mean_charge for g in self.groups if g.ion_type == "acid")
        base_charge = sum(g.mean_charge for g in self.groups if g.ion_type == "base")
        net = self.mean_charge

        # Zwitterion: both acid and base groups carry non-trivial charge
        if acid_charge < -0.1 and base_charge > 0.1:
            return "zwitterion"
        if net < -0.1:
            return "acid"
        if net > 0.1:
            return "base"
        return "neutral"


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def _hhb_acid(pka: float, ph: float) -> tuple[float, float]:
    """Henderson-Hasselbalch for an *acid* group (loses proton).

    Returns (fraction_neutral, mean_charge).
    """
    try:
        ratio = 10 ** (ph - pka)          # [A-] / [HA]
        f_neutral = 1.0 / (1.0 + ratio)   # fraction as HA
        mean_charge = -ratio / (1.0 + ratio)  # mean charge (0 to -1)
    except OverflowError:
        f_neutral, mean_charge = 0.0, -1.0
    return round(f_neutral, 6), round(mean_charge, 6)


def _hhb_base(pka: float, ph: float) -> tuple[float, float]:
    """Henderson-Hasselbalch for a *base* group (gains proton).

    Returns (fraction_neutral, mean_charge).
    """
    try:
        ratio = 10 ** (pka - ph)           # [BH+] / [B]
        f_neutral = 1.0 / (1.0 + ratio)   # fraction as B
        mean_charge = ratio / (1.0 + ratio)   # mean charge (0 to +1)
    except OverflowError:
        f_neutral, mean_charge = 0.0, 1.0
    return round(f_neutral, 6), round(mean_charge, 6)


def ionize(mol: Chem.Mol, ph: float = PH_SC) -> IonizationResult:
    """Predict the ionization state of *mol* at *ph*.

    Uses functional-group SMARTS rules and Henderson-Hasselbalch.  Each
    heavy atom is assigned to at most one rule (first match wins, ordered
    by rule specificity) so the same atom is never double-counted.

    Parameters
    ----------
    mol:
        A valid RDKit molecule (must have implicit Hs computed).
    ph:
        Target pH.  Defaults to 5.5 (stratum corneum surface).

    Returns
    -------
    :class:`IonizationResult` with all matched groups and summary properties.
    """
    result = IonizationResult(ph=ph)
    claimed: set[int] = set()  # atom indices already assigned to a rule

    for rule in _RULES:
        matches = mol.GetSubstructMatches(rule.pattern)
        for match in matches:
            # Skip if any atom in this match is already claimed
            if any(idx in claimed for idx in match):
                continue

            # Use the first (anchor) atom of the match as the "key" atom
            anchor = match[0]
            claimed.update(match)

            if rule.ion_type == "acid":
                f_neutral, m_charge = _hhb_acid(rule.pka, ph)
            else:
                f_neutral, m_charge = _hhb_base(rule.pka, ph)

            result.groups.append(
                IonizableGroup(
                    description=rule.description,
                    pka=rule.pka,
                    ion_type=rule.ion_type,
                    fraction_neutral=f_neutral,
                    mean_charge=m_charge,
                )
            )

    return result
