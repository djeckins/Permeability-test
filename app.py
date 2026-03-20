from __future__ import annotations

import io
from pathlib import Path

import pandas as pd
import streamlit as st
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from epidermal_barrier_screen.io import parse_input
from epidermal_barrier_screen.screen import screen_records

# ── Page config ──────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Epidermal Barrier Screen",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ── Custom CSS ────────────────────────────────────────────────────────────────
st.markdown(
    """
    <style>
    [data-testid="stAppViewContainer"] { background: #f0f4f9; }
    [data-testid="stHeader"] { background: transparent; }

    /* Hide GitHub / Edit-source / Deploy toolbar buttons */
    [data-testid="stToolbar"]          { display: none !important; }
    [data-testid="stDecoration"]       { display: none !important; }
    #MainMenu                          { visibility: hidden !important; }
    footer                             { visibility: hidden !important; }

    .eb-hero {
        background: linear-gradient(135deg, #1a3a5c 0%, #2563a8 100%);
        border-radius: 14px;
        padding: 2rem 2.5rem 1.8rem;
        margin-bottom: 1.6rem;
        color: white;
    }
    .eb-hero h1 { font-size: 2rem; font-weight: 700; margin: 0 0 0.4rem; color: white; }
    .eb-hero p  { font-size: 0.97rem; color: #b8d4f0; margin: 0; }

    .eb-card {
        background: white;
        border-radius: 12px;
        padding: 1.4rem 1.6rem;
        box-shadow: 0 1px 6px rgba(0,0,0,.08);
        margin-bottom: 1.2rem;
    }
    .eb-card h3 { margin: 0 0 1rem; font-size: 1rem; font-weight: 600; color: #1a3a5c; }

    .eb-metric {
        border-radius: 10px;
        padding: 1rem 1.2rem;
        text-align: center;
    }
    .eb-metric .val { font-size: 2rem; font-weight: 700; }
    .eb-metric .lbl { font-size: 0.8rem; text-transform: uppercase; letter-spacing: .05em; }
    .m-total      { background: #e3f2fd; color: #1565c0; }
    .m-pass       { background: #e8f5e9; color: #2e7d32; }
    .m-borderline { background: #fff8e1; color: #e65100; }
    .m-fail       { background: #ffebee; color: #c62828; }

    div[data-testid="stDownloadButton"] button {
        background: #2563a8 !important;
        color: white !important;
        border-radius: 8px !important;
        font-weight: 600 !important;
        font-size: 1rem !important;
        padding: .6rem 1.2rem !important;
    }
    div[data-testid="stDownloadButton"] button:hover {
        background: #1a3a5c !important;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# ── Constants ─────────────────────────────────────────────────────────────────
_STATUS_COLS = [
    "mw_status", "logd_status", "tpsa_status", "hbd_status",
    "hba_status", "rotb_status", "hac_status", "formal_charge_status",
    "ionization_status",
]

# Numeric column → its status column (for direct cell colouring)
_NUMERIC_STATUS = {
    "mw":                    "mw_status",
    "logd":                  "logd_status",
    "tpsa":                  "tpsa_status",
    "hbd":                   "hbd_status",
    "hba":                   "hba_status",
    "rotb":                  "rotb_status",
    "hac":                   "hac_status",
    "formal_charge":         "formal_charge_status",
    "fraction_unionized_pH5_5": "ionization_status",
}

# Columns shown in the Streamlit table (status cols hidden; SMILES at end and narrow)
_DISPLAY_COLS = [
    "name", "parse_status",
    "mw", "clogp", "logd", "logd_source",
    "tpsa", "hbd", "hba", "rotb", "hac", "formal_charge",
    "pka_source", "predicted_pka", "predicted_pka_type",
    "fraction_unionized_pH5_5", "mean_charge_pH5_5", "ionization_class",
    "final_result",
    "input_smiles", "canonical_smiles",
]

_CELL_CSS = {
    "optimal":    "background-color:#c8e6c9;color:#1b5e20",
    "suboptimal": "background-color:#fff9c4;color:#e65100",
    "poor":       "background-color:#ffcdd2;color:#b71c1c",
    "PASS":       "background-color:#a5d6a7;color:#1b5e20;font-weight:bold",
    "BORDERLINE": "background-color:#ffe082;color:#e65100;font-weight:bold",
    "FAIL":       "background-color:#ef9a9a;color:#b71c1c;font-weight:bold",
    "invalid_input": "background-color:#e0e0e0;color:#616161;font-weight:bold",
}

_XLSX_FILLS = {
    "optimal":    PatternFill("solid", fgColor="C8E6C9"),
    "suboptimal": PatternFill("solid", fgColor="FFF9C4"),
    "poor":       PatternFill("solid", fgColor="FFCDD2"),
    "PASS":       PatternFill("solid", fgColor="A5D6A7"),
    "BORDERLINE": PatternFill("solid", fgColor="FFE082"),
    "FAIL":       PatternFill("solid", fgColor="EF9A9A"),
    "invalid_input": PatternFill("solid", fgColor="E0E0E0"),
}


# ── Helpers ───────────────────────────────────────────────────────────────────

def _style_df(df: pd.DataFrame):
    def _color(val):
        return _CELL_CSS.get(str(val), "")

    styler = df.style

    # Colour final_result cell
    if "final_result" in df.columns:
        styler = styler.map(_color, subset=["final_result"])

    # Colour numeric cells directly based on their status column
    for num_col, status_col in _NUMERIC_STATUS.items():
        if num_col in df.columns and status_col in df.columns:
            styler = styler.apply(
                lambda col, sc=status_col: [
                    _CELL_CSS.get(str(df.at[i, sc]), "") for i in col.index
                ],
                subset=[num_col],
            )

    return styler


def _build_xlsx(df: pd.DataFrame) -> bytes:
    buf = io.BytesIO()
    df.to_excel(buf, index=False, engine="openpyxl")
    buf.seek(0)
    wb = load_workbook(buf)
    ws = wb.active

    # Header styling
    hdr_fill = PatternFill("solid", fgColor="1A3A5C")
    hdr_font = Font(color="FFFFFF", bold=True, size=10)
    for cell in ws[1]:
        cell.fill = hdr_fill
        cell.font = hdr_font
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    ws.row_dimensions[1].height = 36

    # Resolve which columns need colour
    headers = [c.value for c in ws[1]]
    colour_cols = {
        col: headers.index(col) + 1
        for col in _STATUS_COLS + ["final_result"]
        if col in headers
    }

    # Apply row colours
    for row in ws.iter_rows(min_row=2):
        for col_name, col_idx in colour_cols.items():
            cell = row[col_idx - 1]
            fill = _XLSX_FILLS.get(str(cell.value))
            if fill:
                cell.fill = fill
        # Zebra stripe on non-coloured cells
        for cell in row:
            if cell.column not in colour_cols.values():
                if cell.row % 2 == 0:
                    cell.fill = PatternFill("solid", fgColor="F5F8FF")
            cell.alignment = Alignment(horizontal="center", vertical="center")

    # Column widths
    for col_idx, col_cells in enumerate(ws.columns, 1):
        max_len = max((len(str(c.value or "")) for c in col_cells), default=8)
        ws.column_dimensions[get_column_letter(col_idx)].width = min(max_len + 3, 32)

    # Freeze header
    ws.freeze_panes = "A2"

    out = io.BytesIO()
    wb.save(out)
    return out.getvalue()


def _detect_mode(filename: str) -> str:
    suffix = Path(filename).suffix.lower()
    if suffix == ".zip":
        return "sdf_zip"
    if suffix in (".sdf", ".mol"):
        return "sdf"
    return "smiles_list"   # .txt or anything else


# ── UI ────────────────────────────────────────────────────────────────────────

st.markdown(
    """
    <div class="eb-hero">
      <h1>🧪 Epidermal Barrier Screen</h1>
      <p>Assess whether a molecule fits the physicochemical window for passive epidermal barrier permeation.<br>
         Supports single SMILES, SMILES lists, SDF files, and ZIP archives of SDFs.</p>
    </div>
    """,
    unsafe_allow_html=True,
)

# ── Input section ─────────────────────────────────────────────────────────────
tab_smiles, tab_file = st.tabs(["✏️  Enter SMILES", "📂  Upload File"])

smiles_payload: str | None = None
file_payload: bytes | None = None
file_name: str | None = None

with tab_smiles:
    st.markdown(
        "Enter a **single SMILES** or a **list of SMILES** (one per line, optional name after a space):"
    )
    smiles_text = st.text_area(
        label="SMILES input",
        placeholder="CC(=O)Oc1ccccc1C(=O)O  Aspirin\nCCO  Ethanol\nCN1CCC[C@H]1c2cccnc2  Nicotine",
        height=160,
        label_visibility="collapsed",
    )
    if smiles_text and smiles_text.strip():
        smiles_payload = smiles_text.strip()

with tab_file:
    st.markdown("Upload an **SDF file**, a **ZIP archive** containing SDFs, or a **TXT file** with SMILES:")
    uploaded = st.file_uploader(
        label="File upload",
        type=["sdf", "zip", "txt"],
        label_visibility="collapsed",
    )
    if uploaded is not None:
        file_payload = uploaded.getvalue()
        file_name = uploaded.name

# ── pH input (required) ───────────────────────────────────────────────────────
st.markdown("<div style='height:0.6rem'/>", unsafe_allow_html=True)
_, ph_col, _ = st.columns([2, 3, 2])
with ph_col:
    ph_input = st.number_input(
        "🧪  pH для расчёта (pKa и logD)",
        min_value=0.0,
        max_value=14.0,
        value=None,
        step=0.1,
        format="%.1f",
        placeholder="Введите pH (0–14), например 5.5",
    )

# ── Run button ────────────────────────────────────────────────────────────────
st.markdown("<div style='height:0.4rem'/>", unsafe_allow_html=True)
_, btn_col, _ = st.columns([2, 3, 2])
with btn_col:
    run = st.button("🔬  Run Screening", type="primary", use_container_width=True)

# ── Processing ────────────────────────────────────────────────────────────────
if run:
    if ph_input is None:
        st.error("Введите pH перед запуском расчёта.")
        st.stop()

    # File takes priority over SMILES text
    if file_payload is not None:
        mode = _detect_mode(file_name or "upload.sdf")
        payload = file_payload
        fname = file_name
    elif smiles_payload:
        mode = "smiles_list" if "\n" in smiles_payload else "smiles"
        payload = smiles_payload
        fname = None
    else:
        st.error("Please enter a SMILES string or upload a file before running.")
        st.stop()

    with st.spinner("Calculating descriptors and screening…"):
        records = parse_input(mode, payload, filename=fname)
        df = screen_records(records, ph=float(ph_input))

    # ── Metrics ───────────────────────────────────────────────────────────────
    total      = len(df)
    n_pass      = int((df["final_result"] == "PASS").sum())
    n_border    = int((df["final_result"] == "BORDERLINE").sum())
    n_fail      = int((df["final_result"] == "FAIL").sum())
    n_invalid   = int((df["parse_status"] != "ok").sum())

    st.markdown("<div style='height:0.4rem'/>", unsafe_allow_html=True)
    m1, m2, m3, m4 = st.columns(4)
    for col, cls, label, value in [
        (m1, "m-total",      "Total",      total),
        (m2, "m-pass",       "PASS",        n_pass),
        (m3, "m-borderline", "BORDERLINE",  n_border),
        (m4, "m-fail",       "FAIL",        n_fail),
    ]:
        col.markdown(
            f'<div class="eb-metric {cls}">'
            f'<div class="val">{value}</div>'
            f'<div class="lbl">{label}</div>'
            f'</div>',
            unsafe_allow_html=True,
        )
    if n_invalid:
        st.warning(f"{n_invalid} record(s) could not be parsed and are shown as `invalid_input`.")

    # ── Results table ─────────────────────────────────────────────────────────
    st.markdown("<div style='height:0.8rem'/>", unsafe_allow_html=True)
    st.markdown("### Results")
    display_order = [c for c in _DISPLAY_COLS if c in df.columns]
    st.dataframe(
        _style_df(df),
        use_container_width=True,
        hide_index=True,
        height=min(40 + len(df) * 36, 600),
        column_order=display_order,
        column_config={
            "input_smiles": st.column_config.TextColumn("SMILES", width="small"),
            "canonical_smiles": st.column_config.TextColumn("Canonical SMILES", width="small"),
            "fraction_unionized_pH5_5": st.column_config.NumberColumn(
                f"f_unionized pH {ph_input:.1f}"
            ),
            "logd": st.column_config.NumberColumn(f"logD pH {ph_input:.1f}"),
        },
    )

    # ── Download ──────────────────────────────────────────────────────────────
    st.markdown("<div style='height:0.6rem'/>", unsafe_allow_html=True)
    xlsx_bytes = _build_xlsx(df)
    _, dl_col, _ = st.columns([2, 3, 2])
    with dl_col:
        st.download_button(
            label="⬇️  Download Results (Excel)",
            data=xlsx_bytes,
            file_name="epidermal_barrier_results.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            use_container_width=True,
        )

# ── Criteria reference ────────────────────────────────────────────────────────
st.markdown("<div style='height:1rem'/>", unsafe_allow_html=True)
with st.expander("📋  Screening criteria reference"):
    st.markdown(
        """
        Per-criterion classification: 🟢 **optimal** · 🟡 **suboptimal** · 🔴 **poor**

        Overall result: **PASS** (all optimal) · **BORDERLINE** (≥ 1 suboptimal, no poor) · **FAIL** (≥ 1 poor)

        | Criterion | Optimal | Suboptimal | Poor |
        |:---|:---:|:---:|:---:|
        | **MW** | < 300 Da | 300–500 Da | > 500 Da |
        | **LogD / cLogP** | 1–3 | 0.5–1 or 3–5 | < 0.5 or > 5 |
        | **TPSA** | < 60 Å² | 60–130 Å² | > 130 Å² |
        | **HBD** | 0–3 | 4–5 | > 5 |
        | **HBA** | 2–8 | 8–10 | > 10 |
        | **RotB** | < 10 | 10–15 | > 15 |
        | **HAC** | < 30 | 30–50 | > 50 |
        | **Formal charge** | 0 | ±1 | ≥ ±2 |
        | **Ionization (pH 5.5)** | f_unionized ≥ 0.8 | 0.5–0.8 | < 0.5 |

        > **Ionization is mandatory** — a molecule that is mostly ionised at the stratum corneum
        > surface pH (5.5) automatically results in FAIL regardless of other properties.

        pKa predicted using **Dimorphite-DL** (Ropp et al., 2019, *J. Cheminformatics* 11:14).
        Overridden by `pKa` / `input_pka` SDF property when present.
        """
    )
