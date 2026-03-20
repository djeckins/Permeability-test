"""Command-line interface for the epidermal barrier screener."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from epidermal_barrier_screen.io import parse_input
from epidermal_barrier_screen.screen import screen_records


def _read_payload(mode: str, input_path: str) -> tuple[str | bytes, str]:
    """Read the input file/string and return (payload, filename)."""
    if mode in ("smiles", "smiles_list"):
        p = Path(input_path)
        if p.exists():
            return p.read_text(encoding="utf-8"), p.name
        # treat raw string as payload
        return input_path, ""
    # sdf / sdf_zip — always read as bytes
    p = Path(input_path)
    return p.read_bytes(), p.name


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Screen molecules against epidermal barrier passage criteria."
    )
    parser.add_argument(
        "--mode",
        required=True,
        choices=["smiles", "smiles_list", "sdf", "sdf_zip"],
        help="Input format.",
    )
    parser.add_argument(
        "--input",
        required=True,
        metavar="INPUT",
        help="SMILES string, path to a SMILES file, SDF file, or ZIP archive.",
    )
    parser.add_argument(
        "--output-prefix",
        default="results",
        metavar="PREFIX",
        help="Output path prefix (without extension). "
             "Both .csv and .xlsx will be written.",
    )
    args = parser.parse_args(argv)

    payload, filename = _read_payload(args.mode, args.input)
    records = parse_input(args.mode, payload, filename=filename or None)
    df = screen_records(records)

    out_prefix = Path(args.output_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    csv_path = out_prefix.with_suffix(".csv")
    xlsx_path = out_prefix.with_suffix(".xlsx")

    df.to_csv(csv_path, index=False)
    df.to_excel(xlsx_path, index=False)

    print(f"Screened {len(df)} record(s).")
    print(f"  CSV  -> {csv_path}")
    print(f"  XLSX -> {xlsx_path}")

    pass_count = (df["final_result"] == "PASS").sum()
    borderline_count = (df["final_result"] == "BORDERLINE").sum()
    fail_count = (df["final_result"] == "FAIL").sum()
    print(f"  PASS: {pass_count}  BORDERLINE: {borderline_count}  FAIL: {fail_count}")
    invalid_count = (df["parse_status"] != "ok").sum()
    if invalid_count:
        print(f"  Invalid (could not parse): {invalid_count}", file=sys.stderr)


if __name__ == "__main__":
    main()
