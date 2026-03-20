"""Ionization state prediction using Dimorphite-DL pKa database.

Predicts the protonation state of a molecule at a target pH (default 5.5,
stratum corneum surface pH) using Dimorphite-DL's experimentally-derived pKa
database and the Henderson-Hasselbalch equation.

Dimorphite-DL (Ropp et al., 2019, J. Cheminformatics 11:14) contains 40+
SMARTS-based functional-group patterns with mean pKa values and standard
deviations derived from curated experimental literature data.  This is
substantially more accurate than a hand-crafted rule table.

Reference
---------
Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019). "Dimorphite-DL: an
open-source program for enumerating the ionization states of drug-like small
molecules." J Cheminform 11, 14. https://doi.org/10.1186/s13321-019-0336-9
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
# Dimorphite-DL group → acid / base classification
#
# Dimorphite-DL's 40+ functional-group patterns are classified into two
# categories:
#   'base'  – neutral nitrogen that gains H⁺ at pH < pKa  (+1 charge)
#   'acid'  – neutral X-H group that loses H⁺ at pH > pKa (-1 charge)
#
# Groups are taken directly from PKaData.get_substructures() names.
# ---------------------------------------------------------------------------
_BASE_GROUPS: frozenset[str] = frozenset(
    {
        "AmidineGuanidine1",       # R-C(=NH)-NR2,  pKa ~12 (strong base)
        "AmidineGuanidine2",       # R2N-C=NR,       pKa ~10
        "Anilines_primary",        # Ar-NH2,          pKa ~3.9 (weak base)
        "Anilines_secondary",      # Ar-NHR,          pKa ~4.3
        "Anilines_tertiary",       # Ar-NR2,          pKa ~4.2
        "Aromatic_nitrogen_unprotonated",   # pyridine-type,  pKa ~4.4
        "*Aromatic_nitrogen_protonated",    # protonated aromatic N, pKa ~7.2 (dimorphite_dl v2 name has * prefix)
        "Amines_primary_secondary_tertiary",  # aliphatic N, pKa ~8.2
        "Primary_hydroxyl_amine",  # NH2OH,           pKa ~4.0
    }
)

# Groups to skip outright (sentinel / non-ionizable at any physiological pH)
_SKIP_GROUPS: frozenset[str] = frozenset(
    {
        "Nitro",  # pKa = -1000, never ionizes at physiological pH
    }
)

# ---------------------------------------------------------------------------
# Public data-classes  (unchanged interface from previous implementation)
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
    #: True when the pKaPredict ML model successfully refined at least one group pKa
    pka_ml_used: bool = False

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
    def fraction_neutral_total(self) -> float:
        """Fraction of molecules where ALL ionizable groups are simultaneously neutral.

        Computed as the product of per-group fraction_neutral values, which
        is the correct approximation for independent equilibria.  For a single
        ionizable group this equals fraction_neutral_dominant.
        """
        result = 1.0
        for g in self.groups:
            result *= g.fraction_neutral
        return result

    @property
    def ionization_class(self) -> str:
        """Broad ionization class at the target pH.

        Returns one of:
        ``'non_ionizable'``, ``'neutral'``, ``'acid'``, ``'base'``,
        ``'zwitterion'``.
        """
        if self.non_ionizable:
            return "non_ionizable"

        acid_charge = sum(
            g.mean_charge for g in self.groups if g.ion_type == "acid"
        )
        base_charge = sum(
            g.mean_charge for g in self.groups if g.ion_type == "base"
        )
        net = self.mean_charge

        if acid_charge < -0.1 and base_charge > 0.1:
            return "zwitterion"
        if net < -0.1:
            return "acid"
        if net > 0.1:
            return "base"
        return "neutral"


# ---------------------------------------------------------------------------
# Henderson-Hasselbalch helpers  (unchanged)
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
# Dimorphite-DL integration
# ---------------------------------------------------------------------------


def _get_detector():
    """Lazily instantiate (and cache) the Dimorphite-DL detector."""
    if _get_detector._instance is None:
        from dimorphite_dl.protonate.detect import ProtonationSiteDetector

        _get_detector._instance = ProtonationSiteDetector()
    return _get_detector._instance


_get_detector._instance = None  # type: ignore[attr-defined]


# ---------------------------------------------------------------------------
# pKaPredict integration  (molecule-specific ML pKa refinement)
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
    """Return a molecule-specific pKa prediction from pKaPredict, or *None* on failure.

    Failures are logged to stderr so they are visible in the Streamlit terminal
    instead of being silently swallowed.
    """
    import sys

    try:
        from pkapredict import predict_pKa  # type: ignore[import]

        model, feat = _get_pka_model()
        val = float(predict_pKa(smiles, model, feat))
        return val
    except Exception as exc:
        print(f"[pKaPredict] WARNING: prediction failed for {smiles!r}: {exc}", file=sys.stderr)
        return None


def ionize(mol: Chem.Mol, ph: float = PH_SC) -> IonizationResult:
    """Predict the ionization state of *mol* at *ph*.

    Uses Dimorphite-DL's experimentally-derived pKa database to detect
    ionizable functional groups, then applies Henderson-Hasselbalch to
    compute the fraction-neutral and mean charge of each group.

    Parameters
    ----------
    mol:
        A valid RDKit molecule.
    ph:
        Target pH.  Defaults to 5.5 (stratum corneum surface).

    Returns
    -------
    :class:`IonizationResult` with all matched groups and summary properties.
    """
    result = IonizationResult(ph=ph)

    try:
        from dimorphite_dl.mol import MoleculeRecord

        smiles = Chem.MolToSmiles(mol)
        if not smiles:
            return result

        detector = _get_detector()
        rec = MoleculeRecord(smiles)
        _, sites = detector.find_sites(rec)
    except Exception:
        # If Dimorphite-DL fails for any reason, return empty (non-ionizable)
        return result

    # Deduplicate: the same heavy atom can appear in multiple SMARTS matches
    # (e.g. diethylamine N matches twice from each C-N bond).
    # Use the resolved atom index from PKaDatum.idx_site as the dedup key.
    claimed_atoms: set[int] = set()

    for site in sites:
        # Skip sentinel / non-physiological groups
        if site.name in _SKIP_GROUPS:
            continue

        ion_type: IonType = (
            "base" if site.name in _BASE_GROUPS else "acid"
        )

        for pka_datum in site.pkas:
            pka = pka_datum.mean

            # Skip groups that are effectively always neutral (contribute ≈1.0 to
            # the fraction_neutral product and carry no ionization information):
            #   • ACID  with pKa > 12 – never deprotonated at physiological pH
            #     (e.g. Alcohol pKa≈14.8, *Amide pKa≈12.0)
            #   • BASE  with pKa < 1  – never protonated at physiological pH
            # Note: BASE groups with pKa > 12 (e.g. guanidinium pKa≈12.0) are
            # kept because they are always charged, which DOES affect ionization.
            if ion_type == "acid" and pka > 12.0:
                continue
            if ion_type == "base" and pka < 1.0:
                continue

            # Resolve the actual heavy-atom index this pKa datum modifies
            atom_idx = site.idxs_match[pka_datum.idx_site]
            if atom_idx in claimed_atoms:
                continue
            claimed_atoms.add(atom_idx)

            if ion_type == "acid":
                f_neutral, m_charge = _hhb_acid(pka, ph)
            else:
                f_neutral, m_charge = _hhb_base(pka, ph)

            result.groups.append(
                IonizableGroup(
                    description=site.name,
                    pka=round(pka, 3),
                    ion_type=ion_type,
                    fraction_neutral=f_neutral,
                    mean_charge=m_charge,
                )
            )

    # ── pKaPredict refinement ────────────────────────────────────────────────
    # Dimorphite-DL assigns category-level mean pKa values (all carboxylic
    # acids → 3.46, all aliphatic amines → 8.16, etc.).  Where possible, we
    # replace the pKa of the *most influential* ionizable group with a
    # molecule-specific prediction from the pKaPredict LGBMRegressor model,
    # which is trained on RDKit descriptors and gives distinct values per
    # molecule.
    #
    # Strategy: pKaPredict outputs a single pKa per molecule.  We assign that
    # value to the group whose Dimorphite-DL pKa is numerically closest to the
    # ML prediction (most likely to refer to the same ionisation event), then
    # recompute its fraction_neutral and mean_charge via Henderson-Hasselbalch.
    if result.groups:
        ml = _ml_pka(smiles)
        if ml is not None:
            target = min(result.groups, key=lambda g: abs(g.pka - ml))
            target.pka = round(ml, 3)
            if target.ion_type == "acid":
                target.fraction_neutral, target.mean_charge = _hhb_acid(ml, ph)
            else:
                target.fraction_neutral, target.mean_charge = _hhb_base(ml, ph)
            result.pka_ml_used = True

    return result
