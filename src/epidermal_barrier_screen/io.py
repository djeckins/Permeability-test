"""Parse various input formats into a list of molecule records."""
from __future__ import annotations

import io
import zipfile
from typing import Any

from rdkit import Chem

# SDF property name aliases
_PKA_KEYS = ("pKa", "PKA", "predicted_pKa", "input_pka")
_LOGD_KEYS = ("LogD", "LOGD", "logD_7_4", "input_logd_7_4")
_NAME_KEYS = ("_Name", "Name", "TITLE", "ID", "Compound")


def _try_float(value: str | None) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def _sdf_prop(mol: Chem.Mol, keys: tuple[str, ...]) -> str | None:
    """Return the first matching property value from *mol*, or None."""
    for key in keys:
        if mol.HasProp(key):
            return mol.GetProp(key)
    return None


def _record_from_mol(mol: Chem.Mol | None, input_smiles: str | None = None) -> dict[str, Any]:
    """Build a molecule record from an RDKit Mol (may be None for invalid input)."""
    if mol is None:
        return {
            "name": None,
            "input_smiles": input_smiles,
            "canonical_smiles": None,
            "mol": None,
            "parse_status": "invalid",
            "input_pka": None,
            "input_logd_7_4": None,
        }

    name = _sdf_prop(mol, _NAME_KEYS) or input_smiles
    canonical_smiles = Chem.MolToSmiles(mol)
    return {
        "name": name,
        "input_smiles": input_smiles or canonical_smiles,
        "canonical_smiles": canonical_smiles,
        "mol": mol,
        "parse_status": "ok",
        "input_pka": _try_float(_sdf_prop(mol, _PKA_KEYS)),
        "input_logd_7_4": _try_float(_sdf_prop(mol, _LOGD_KEYS)),
    }


def _parse_smiles(smiles: str) -> dict[str, Any]:
    smiles = smiles.strip()
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    return _record_from_mol(mol, input_smiles=smiles or None)


def _records_from_sdf_supplier(supplier: Chem.ForwardSDMolSupplier) -> list[dict[str, Any]]:
    records = []
    for mol in supplier:
        records.append(_record_from_mol(mol))
    return records


def parse_input(
    mode: str,
    payload: str | bytes,
    *,
    filename: str | None = None,
) -> list[dict[str, Any]]:
    """Parse *payload* according to *mode* and return a list of molecule records.

    Parameters
    ----------
    mode:
        One of ``'smiles'``, ``'smiles_list'``, ``'sdf'``, ``'sdf_zip'``.
    payload:
        The raw input (string for SMILES modes, bytes for file modes).
    filename:
        Original file name (used for SDF/ZIP modes; ignored for SMILES modes).

    Returns
    -------
    List of record dicts.  Each record has keys:
        name, input_smiles, canonical_smiles, mol, parse_status,
        input_pka, input_logd_7_4.
    """
    if mode == "smiles":
        return [_parse_smiles(str(payload))]

    if mode == "smiles_list":
        records = []
        for raw_line in str(payload).splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            smiles = parts[0]
            name_hint = parts[1] if len(parts) > 1 else None
            record = _parse_smiles(smiles)
            if name_hint and record["name"] == smiles:
                record["name"] = name_hint
            records.append(record)
        return records

    if mode == "sdf":
        if not isinstance(payload, (bytes, bytearray)):
            raise TypeError("SDF mode expects bytes payload")
        supplier = Chem.ForwardSDMolSupplier(io.BytesIO(payload), removeHs=False)
        return _records_from_sdf_supplier(supplier)

    if mode == "sdf_zip":
        if not isinstance(payload, (bytes, bytearray)):
            raise TypeError("sdf_zip mode expects bytes payload")
        records = []
        with zipfile.ZipFile(io.BytesIO(payload)) as zf:
            sdf_names = [n for n in zf.namelist() if n.lower().endswith(".sdf")]
            if not sdf_names:
                raise ValueError("ZIP archive contains no .sdf files")
            for sdf_name in sorted(sdf_names):
                sdf_bytes = zf.read(sdf_name)
                supplier = Chem.ForwardSDMolSupplier(io.BytesIO(sdf_bytes), removeHs=False)
                records.extend(_records_from_sdf_supplier(supplier))
        return records

    raise ValueError(f"Unknown input mode: {mode!r}")
