"""Ionization state prediction using the pKaLearn GNN model.

Reference: https://github.com/MoitessierLab/pKaLearn
"""
from __future__ import annotations

import os
import sys
import tempfile

import pandas as pd
from rdkit import Chem

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PH_SC: float = 5.5  # stratum corneum surface pH

# ---------------------------------------------------------------------------
# Paths to the vendored pKaLearn submodule
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
# src/epidermal_barrier_screen/ → ../../  = repo root
_REPO_ROOT = os.path.normpath(os.path.join(_HERE, "..", ".."))
_PKALEARN_GNN = os.path.join(_REPO_ROOT, "third_party", "pKaLearn", "GNN")
_PKALEARN_MODEL = os.path.join(_REPO_ROOT, "third_party", "pKaLearn", "Model")


# ---------------------------------------------------------------------------
# Henderson-Hasselbalch helpers
# ---------------------------------------------------------------------------


def _hhb_acid(pka: float, ph: float) -> tuple[float, float]:
    """Henderson-Hasselbalch for an acid group (loses proton).

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
    """Henderson-Hasselbalch for a base group (gains proton).

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
# Ion type detection from pKaLearn output
# ---------------------------------------------------------------------------


def _detect_ion_type(original_smiles: str, protonated_smiles: str) -> str:
    """Infer acid / base / non_ionizable by comparing implicit-H counts.

    pKaLearn returns the dominant protonation state at the target pH.
    If the dominant state has more H atoms than the input SMILES the
    ionizable centre is a *base* (it gains a proton below its pKa).
    If it has fewer H atoms the centre is an *acid* (it loses a proton
    above its pKa).
    """
    mol_orig = Chem.MolFromSmiles(original_smiles)
    mol_prot = Chem.MolFromSmiles(protonated_smiles)
    if mol_orig is None or mol_prot is None:
        return "non_ionizable"
    h_orig = sum(a.GetTotalNumHs() for a in mol_orig.GetAtoms())
    h_prot = sum(a.GetTotalNumHs() for a in mol_prot.GetAtoms())
    if h_prot > h_orig:
        return "base"
    if h_prot < h_orig:
        return "acid"
    return "non_ionizable"


# ---------------------------------------------------------------------------
# pKaLearn integration
# ---------------------------------------------------------------------------


def _ensure_pkalearn_importable() -> None:
    """Add the pKaLearn GNN directory to *sys.path* if needed."""
    gnn_dir = os.path.normpath(_PKALEARN_GNN)
    if not os.path.isdir(gnn_dir):
        raise RuntimeError(
            f"pKaLearn GNN directory not found at {gnn_dir!r}.\n"
            "Initialise the submodule with:\n"
            "  git submodule update --init --recursive"
        )
    if gnn_dir not in sys.path:
        sys.path.insert(0, gnn_dir)


def predict_pka(smiles: str, name: str = "", ph: float = PH_SC) -> tuple[float | None, str]:
    """Predict pKa for *smiles* using the pKaLearn GNN model.

    Parameters
    ----------
    smiles:
        Canonical SMILES string for the molecule.
    name:
        Optional molecule name (used only as a label inside pKaLearn).
    ph:
        Target pH for protonation-state determination. Defaults to
        :data:`PH_SC` (5.5, stratum corneum surface pH).

    Returns
    -------
    (pka_value, ion_type)
    - *pka_value*: predicted pKa as a float, or ``None`` for
      non-ionizable molecules or when the prediction fails.
    - *ion_type*: one of ``'acid'``, ``'base'``, or ``'non_ionizable'``.
    """
    _ensure_pkalearn_importable()

    try:
        from predict import predict as _pk_predict  # pKaLearn's predict module
    except ImportError as exc:
        raise RuntimeError(
            f"Cannot import pKaLearn predict module: {exc}\n"
            "Ensure pKaLearn dependencies are installed:\n"
            "  pip install torch torch_geometric"
        ) from exc

    df = pd.DataFrame({"Name": [name or smiles], "Smiles": [smiles]})

    # pKaLearn calls argsParser() which reads sys.argv; temporarily neutralise it
    # so Streamlit / CLI arguments do not confuse argparse.
    old_argv = sys.argv[:]
    sys.argv = ["predict"]

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            # dump_datasets writes the pickled dataset here; using an absolute
            # path makes os.path.join(dirname(args.csv_path), ..., infer_pickled)
            # resolve to just pkl_path on POSIX systems.
            pkl_path = os.path.join(tmpdir, "infer.pkl")
            pka_list, smiles_list = _pk_predict(
                csv_file=df,
                pH=ph,
                model_dir=os.path.normpath(_PKALEARN_MODEL),
                data_path=tmpdir + os.sep,
                infer_pickled=pkl_path,
            )
    except Exception as exc:
        print(
            f"[pKaLearn] WARNING: prediction failed for {smiles!r}: {exc}",
            file=sys.stderr,
        )
        return None, "non_ionizable"
    finally:
        sys.argv = old_argv

    if not pka_list:
        return None, "non_ionizable"

    pka_raw = pka_list[0]
    prot_smiles = smiles_list[0] if smiles_list else smiles

    # pKaLearn returns the string 'NaN' for non-ionizable molecules
    if pka_raw == "NaN" or (isinstance(pka_raw, float) and pka_raw != pka_raw):
        return None, "non_ionizable"

    try:
        pka_val = float(pka_raw)
    except (ValueError, TypeError):
        return None, "non_ionizable"

    ion_type = _detect_ion_type(smiles, prot_smiles)
    return pka_val, ion_type
