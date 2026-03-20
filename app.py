from __future__ import annotations

import io
import os
import sys

# ── Fix import path for Streamlit Cloud (package lives in src/) ──────────────
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

import pandas as pd
import streamlit as st
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from epidermal_barrier_screen.io import parse_input
from epidermal_barrier_screen.screen import screen_records

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Epidermal Barrier Screening Tool",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ── CSS ───────────────────────────────────────────────────────────────────────
st.markdown(
    """
<style>
/* ---------- page background ---------- */
[data-testid="stAppViewContainer"] > .main {
    background: linear-gradient(155deg, #cfe4f8 0%, #e4f0fb 45%, #d0e6f7 100%);
    min-height: 100vh;
}
[data-testid="stHeader"],
[data-testid="stToolbar"] { display: none !important; }
.block-container { padding: 1.2rem 1.8rem !important; max-width: 100% !important; }

/* ---------- hero banner ---------- */
.eb-hero {
    background: linear-gradient(120deg, #b8d8f5 0%, #d4ebff 60%, #bcd9f0 100%);
    border-radius: 18px;
    padding: 1.8rem 2.2rem;
    margin-bottom: 1.4rem;
    position: relative;
    overflow: hidden;
    box-shadow: 0 4px 20px rgba(37,99,168,.12);
}
.eb-hero h1 {
    font-size: 1.85rem;
    font-weight: 800;
    color: #0d2d52;
    margin: 0 0 .3rem;
    letter-spacing: -.02em;
}
.eb-hero .sub { color: #2e6ca0; font-size: .95rem; margin-bottom: 1rem; }
.eb-pills { display: flex; flex-wrap: wrap; gap: .4rem; }
.eb-pill {
    background: rgba(255,255,255,.82);
    border: 1px solid rgba(37,99,168,.18);
    border-radius: 20px;
    padding: .22rem .75rem;
    font-size: .78rem;
    color: #1a5fa8;
    font-weight: 500;
    backdrop-filter: blur(4px);
}

/* ---------- white cards ---------- */
.eb-card {
    background: white;
    border-radius: 16px;
    padding: 1.4rem 1.5rem;
    box-shadow: 0 2px 16px rgba(37,99,168,.1);
    margin-bottom: .8rem;
}
.eb-card-title {
    font-size: 1rem;
    font-weight: 700;
    color: #0d2d52;
    margin: 0 0 1rem;
}

/* ---------- tabs ---------- */
[data-testid="stTabs"] [role="tablist"] {
    background: #edf4fb;
    border-radius: 10px;
    padding: 3px;
    gap: 2px;
    border-bottom: none !important;
}
[data-testid="stTabs"] [role="tab"] {
    border-radius: 8px !important;
    font-weight: 600 !important;
    font-size: .83rem !important;
    color: #5a8ab8 !important;
    border: none !important;
    padding: .35rem .9rem !important;
}
[data-testid="stTabs"] [role="tab"][aria-selected="true"] {
    background: #2563a8 !important;
    color: white !important;
}
[data-testid="stTabs"] [role="tabpanel"] { padding-top: .8rem !important; }

/* ---------- file uploader ---------- */
[data-testid="stFileUploader"] {
    background: #f0f7ff;
    border: 1.5px dashed #90b8e0;
    border-radius: 12px;
    padding: .5rem !important;
}
[data-testid="stFileUploader"] label { color: #2563a8 !important; font-weight: 500 !important; }

/* ---------- primary button (Run Screening) ---------- */
div[data-testid="stButton"] > button[kind="primary"] {
    background: linear-gradient(135deg, #1a5fa8 0%, #2196d4 100%) !important;
    color: white !important;
    border: none !important;
    border-radius: 12px !important;
    padding: .75rem 1rem !important;
    font-size: 1rem !important;
    font-weight: 700 !important;
    letter-spacing: .03em !important;
    box-shadow: 0 4px 18px rgba(37,99,168,.35) !important;
    transition: transform .15s, box-shadow .15s !important;
    width: 100% !important;
}
div[data-testid="stButton"] > button[kind="primary"]:hover {
    transform: translateY(-2px) !important;
    box-shadow: 0 7px 24px rgba(37,99,168,.45) !important;
}

/* ---------- download button ---------- */
div[data-testid="stDownloadButton"] button {
    background: linear-gradient(135deg, #1a5fa8 0%, #2196d4 100%) !important;
    color: white !important;
    border: none !important;
    border-radius: 10px !important;
    font-weight: 600 !important;
    font-size: .9rem !important;
    box-shadow: 0 2px 10px rgba(37,99,168,.3) !important;
    width: 100% !important;
}

/* ---------- results table ---------- */
.rt-wrap { overflow-x: auto; }
table.rt {
    width: 100%;
    border-collapse: collapse;
    font-size: .86rem;
}
table.rt thead tr {
    border-bottom: 2px solid #e0ecfa;
}
table.rt th {
    padding: .55rem .9rem;
    text-align: left;
    color: #4a7aaa;
    font-weight: 600;
    font-size: .78rem;
    white-space: nowrap;
}
table.rt td {
    padding: .5rem .9rem;
    border-bottom: 1px solid #eef3fb;
    color: #0d2d52;
    white-space: nowrap;
}
table.rt tr:last-child td { border-bottom: none; }
table.rt tr:hover td { background: #f4f8ff; }

/* badge */
.bdg {
    display: inline-block;
    padding: .18rem .7rem;
    border-radius: 20px;
    font-size: .75rem;
    font-weight: 700;
    letter-spacing: .02em;
}
.bdg-opt  { background: #00BCD4; color: white; }
.bdg-sub  { background: #FFB300; color: white; }
.bdg-poor { background: #EF5350; color: white; }
.bdg-pass { background: #43A047; color: white; }
.bdg-bord { background: #FB8C00; color: white; }
.bdg-fail { background: #E53935; color: white; }
.bdg-inv  { background: #90A4AE; color: white; }

/* placeholder */
.rt-placeholder {
    display: flex; flex-direction: column;
    align-items: center; justify-content: center;
    min-height: 280px; color: #90b4d8;
    gap: .6rem; text-align: center;
}
.rt-placeholder .ico { font-size: 3rem; }
.rt-placeholder p { font-size: .9rem; margin: 0; }

/* metric mini-cards */
.mini-metrics { display: flex; gap: .7rem; margin-bottom: 1rem; }
.mini-m {
    flex: 1; border-radius: 10px; padding: .6rem .8rem; text-align: center;
}
.mini-m .mv { font-size: 1.5rem; font-weight: 800; }
.mini-m .ml { font-size: .7rem; text-transform: uppercase; letter-spacing: .06em; font-weight: 600; }
.mm-total { background: #e3f2fd; color: #1565c0; }
.mm-pass  { background: #e8f5e9; color: #2e7d32; }
.mm-bord  { background: #fff8e1; color: #e65100; }
.mm-fail  { background: #ffebee; color: #c62828; }
</style>
""",
    unsafe_allow_html=True,
)

# ── Session state ─────────────────────────────────────────────────────────────
if "df" not in st.session_state:
    st.session_state.df = None

# ── Helpers ───────────────────────────────────────────────────────────────────
_BADGE_MAP = {
    "optimal":     ("Optimal",    "bdg-opt"),
    "suboptimal":  ("Suboptimal", "bdg-sub"),
    "poor":        ("Poor",       "bdg-poor"),
    "PASS":        ("PASS",       "bdg-pass"),
    "BORDERLINE":  ("BORDERLINE", "bdg-bord"),
    "FAIL":        ("FAIL",       "bdg-fail"),
    "invalid_input": ("Invalid",  "bdg-inv"),
}

def _badge(val: str) -> str:
    if val in _BADGE_MAP:
        label, cls = _BADGE_MAP[val]
        return f'<span class="bdg {cls}">{label}</span>'
    return val if val else "—"

def _fmt(v, decimals: int = 2) -> str:
    if v is None or (isinstance(v, float) and str(v) == "nan"):
        return "—"
    if isinstance(v, float):
        return f"{v:.{decimals}f}"
    return str(v)

def _build_results_html(df: pd.DataFrame) -> str:
    rows_html = ""
    for _, r in df.iterrows():
        name = r.get("name") or r.get("input_smiles", "—")
        if name and len(str(name)) > 22:
            name = str(name)[:20] + "…"
        rows_html += f"""<tr>
          <td title="{r.get('input_smiles','')}">{name}</td>
          <td>{_fmt(r.get('mw'), 2)}</td>
          <td>{_fmt(r.get('clogp'), 2)}</td>
          <td>{_badge(r.get('logd_status',''))}</td>
          <td>{_fmt(r.get('tpsa'), 1)}</td>
          <td>{_badge(r.get('tpsa_status',''))}</td>
          <td>{_fmt(r.get('hbd'), 0)}</td>
          <td>{_fmt(r.get('hba'), 0)}</td>
          <td>{_fmt(r.get('fraction_unionized_pH5_5'), 3)}</td>
          <td>{_badge(r.get('ionization_status',''))}</td>
          <td>{_badge(r.get('final_result',''))}</td>
        </tr>"""
    return f"""<div class="rt-wrap"><table class="rt">
      <thead><tr>
        <th>Compound</th><th>MW</th><th>LogP</th>
        <th>Lipophilicity</th><th>TPSA</th><th>TPSA status</th>
        <th>HBD</th><th>HBA</th><th>f_unionized</th>
        <th>Ionization</th><th>Result</th>
      </tr></thead>
      <tbody>{rows_html}</tbody>
    </table></div>"""

def _detect_mode(filename: str) -> str:
    from pathlib import Path
    s = Path(filename).suffix.lower()
    if s == ".zip":
        return "sdf_zip"
    if s in (".sdf", ".mol"):
        return "sdf"
    return "smiles_list"

_XLSX_FILLS = {
    "optimal":    PatternFill("solid", fgColor="B2EBF2"),
    "suboptimal": PatternFill("solid", fgColor="FFE082"),
    "poor":       PatternFill("solid", fgColor="FFCDD2"),
    "PASS":       PatternFill("solid", fgColor="C8E6C9"),
    "BORDERLINE": PatternFill("solid", fgColor="FFE0B2"),
    "FAIL":       PatternFill("solid", fgColor="FFCDD2"),
    "invalid_input": PatternFill("solid", fgColor="ECEFF1"),
}
_STATUS_COLS = [
    "mw_status","logd_status","tpsa_status","hbd_status","hba_status",
    "rotb_status","hac_status","formal_charge_status","ionization_status","final_result",
]

def _build_xlsx(df: pd.DataFrame) -> bytes:
    buf = io.BytesIO()
    df.to_excel(buf, index=False, engine="openpyxl")
    buf.seek(0)
    wb = load_workbook(buf)
    ws = wb.active
    hdr_fill = PatternFill("solid", fgColor="0D2D52")
    hdr_font = Font(color="FFFFFF", bold=True, size=9)
    for cell in ws[1]:
        cell.fill = hdr_fill
        cell.font = hdr_font
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    ws.row_dimensions[1].height = 34
    headers = [c.value for c in ws[1]]
    colour_idx = {col: headers.index(col)+1 for col in _STATUS_COLS if col in headers}
    for row in ws.iter_rows(min_row=2):
        for col_name, col_i in colour_idx.items():
            cell = row[col_i-1]
            fill = _XLSX_FILLS.get(str(cell.value))
            if fill:
                cell.fill = fill
        for cell in row:
            if cell.column not in colour_idx.values() and row[0].row % 2 == 0:
                cell.fill = PatternFill("solid", fgColor="F5F8FF")
            cell.alignment = Alignment(horizontal="center", vertical="center")
    for col_i, col_cells in enumerate(ws.columns, 1):
        w = max((len(str(c.value or "")) for c in col_cells), default=8)
        ws.column_dimensions[get_column_letter(col_i)].width = min(w+3, 30)
    ws.freeze_panes = "A2"
    out = io.BytesIO()
    wb.save(out)
    return out.getvalue()

# ─────────────────────────────────────────────────────────────────────────────
# HERO BANNER
# ─────────────────────────────────────────────────────────────────────────────
_MOL_SVG = """
<svg width="260" height="130" viewBox="0 0 260 130"
     style="position:absolute;right:16px;top:10px;opacity:.55"
     xmlns="http://www.w3.org/2000/svg">
  <polygon points="210,8 228,19 228,41 210,52 192,41 192,19"
    fill="none" stroke="#1a5fa8" stroke-width="1.8"/>
  <polygon points="244,28 262,39 262,61 244,72 226,61 226,39"
    fill="none" stroke="#1a5fa8" stroke-width="1.2" opacity=".5"/>
  <polygon points="218,56 236,67 236,89 218,100 200,89 200,67"
    fill="none" stroke="#1a5fa8" stroke-width="1.2" opacity=".4"/>
  <circle cx="160" cy="65" r="20" fill="rgba(37,99,168,.18)"
    stroke="#2563a8" stroke-width="1.8"/>
  <circle cx="210" cy="35" r="13" fill="rgba(26,159,212,.25)"
    stroke="#1a9fd4" stroke-width="1.5"/>
  <circle cx="248" cy="95" r="10" fill="rgba(37,99,168,.18)"
    stroke="#2563a8" stroke-width="1.2"/>
  <circle cx="175" cy="108" r="8" fill="rgba(26,159,212,.2)"
    stroke="#1a9fd4" stroke-width="1.2"/>
  <circle cx="130" cy="38" r="7" fill="rgba(37,99,168,.15)"
    stroke="#2563a8" stroke-width="1"/>
  <line x1="180" y1="57" x2="197" y2="43" stroke="#2563a8" stroke-width="1.6" opacity=".55"/>
  <line x1="223" y1="42" x2="238" y2="85" stroke="#2563a8" stroke-width="1.4" opacity=".4"/>
  <line x1="168" y1="78" x2="172" y2="100" stroke="#2563a8" stroke-width="1.4" opacity=".4"/>
  <line x1="145" y1="60" x2="137" y2="45" stroke="#2563a8" stroke-width="1.2" opacity=".35"/>
</svg>"""

st.markdown(
    f"""<div class="eb-hero">
  {_MOL_SVG}
  <h1>🧪 Epidermal Barrier Screening Tool</h1>
  <p class="sub">Evaluate compounds for epidermal barrier permeation potential</p>
  <div class="eb-pills">
    <span class="eb-pill">⚖️ MW – Dalton</span>
    <span class="eb-pill">💧 LogP</span>
    <span class="eb-pill">⬤ TPSA</span>
    <span class="eb-pill">🔗 HBD &amp; HBA</span>
    <span class="eb-pill">🫧 Lipophilicity</span>
    <span class="eb-pill">🔄 RotB &amp; HAC</span>
    <span class="eb-pill">⚡ Ionization at pH 5.5</span>
    <span class="eb-pill">🔋 Formal Charge</span>
  </div>
</div>""",
    unsafe_allow_html=True,
)

# ─────────────────────────────────────────────────────────────────────────────
# TWO COLUMNS
# ─────────────────────────────────────────────────────────────────────────────
col_in, col_out = st.columns([1, 1.7], gap="large")

# ── LEFT: Input ───────────────────────────────────────────────────────────────
with col_in:
    st.markdown('<div class="eb-card">', unsafe_allow_html=True)
    st.markdown('<p class="eb-card-title">📋 Upload Compound Data</p>', unsafe_allow_html=True)

    tab_s, tab_l, tab_f = st.tabs(["SMILES", "SMILES List", "SDF Archive"])

    smiles_payload = None
    file_payload   = None
    file_name      = None

    with tab_s:
        smiles_single = st.text_input(
            "Enter SMILES",
            placeholder="CC(=O)Oc1ccccc1C(=O)O",
            label_visibility="collapsed",
        )
        if smiles_single.strip():
            smiles_payload = smiles_single.strip()

    with tab_l:
        smiles_list = st.text_area(
            "SMILES list (one per line, optional name after space)",
            placeholder="CCO  Ethanol\nCC(=O)Oc1ccccc1C(=O)O  Aspirin\nCN1CCC[C@H]1c2cccnc2  Nicotine",
            height=130,
            label_visibility="collapsed",
        )
        if smiles_list.strip():
            smiles_payload = smiles_list.strip()

    with tab_f:
        uploaded = st.file_uploader(
            "Upload SDF, ZIP archive with SDFs, or TXT with SMILES",
            type=["sdf", "zip", "txt"],
            label_visibility="collapsed",
        )
        if uploaded is not None:
            file_payload = uploaded.getvalue()
            file_name    = uploaded.name
            smiles_payload = None   # file takes priority

    st.markdown("<div style='height:.6rem'/>", unsafe_allow_html=True)
    run = st.button("🔬  Run Screening", type="primary", use_container_width=True)
    st.markdown("</div>", unsafe_allow_html=True)

# ── RIGHT: Results ────────────────────────────────────────────────────────────
with col_out:
    st.markdown('<div class="eb-card">', unsafe_allow_html=True)

    # Header row with title + download
    h1, h2 = st.columns([2, 1])
    h1.markdown('<p class="eb-card-title">📊 Screening Results</p>', unsafe_allow_html=True)

    df = st.session_state.df

    if df is not None and len(df):
        xlsx_bytes = _build_xlsx(df)
        h2.download_button(
            "⬇️ Download Excel",
            data=xlsx_bytes,
            file_name="epidermal_barrier_results.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            use_container_width=True,
        )

        # Mini metrics
        total   = len(df)
        n_pass  = int((df["final_result"] == "PASS").sum())
        n_bord  = int((df["final_result"] == "BORDERLINE").sum())
        n_fail  = int((df["final_result"] == "FAIL").sum())
        st.markdown(
            f"""<div class="mini-metrics">
              <div class="mini-m mm-total"><div class="mv">{total}</div><div class="ml">Total</div></div>
              <div class="mini-m mm-pass"><div class="mv">{n_pass}</div><div class="ml">Pass</div></div>
              <div class="mini-m mm-bord"><div class="mv">{n_bord}</div><div class="ml">Borderline</div></div>
              <div class="mini-m mm-fail"><div class="mv">{n_fail}</div><div class="ml">Fail</div></div>
            </div>""",
            unsafe_allow_html=True,
        )

        st.markdown(_build_results_html(df), unsafe_allow_html=True)
    else:
        st.markdown(
            """<div class="rt-placeholder">
              <div class="ico">🔬</div>
              <p><strong>No results yet</strong></p>
              <p>Enter a SMILES or upload a file,<br>then click <strong>Run Screening</strong></p>
            </div>""",
            unsafe_allow_html=True,
        )

    st.markdown("</div>", unsafe_allow_html=True)

# ── Process run ───────────────────────────────────────────────────────────────
if run:
    if file_payload is not None:
        mode    = _detect_mode(file_name or "upload.sdf")
        payload = file_payload
        fname   = file_name
    elif smiles_payload:
        mode    = "smiles_list" if "\n" in smiles_payload else "smiles"
        payload = smiles_payload
        fname   = None
    else:
        st.error("Please enter a SMILES or upload a file before running.")
        st.stop()

    with st.spinner("Calculating descriptors and screening molecules…"):
        records = parse_input(mode, payload, filename=fname)
        st.session_state.df = screen_records(records)
    st.rerun()

# ── Criteria reference ────────────────────────────────────────────────────────
st.markdown("<div style='height:.4rem'/>", unsafe_allow_html=True)
with st.expander("📋 Screening criteria reference"):
    st.markdown(
        """
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
| **Ionization pH 5.5** | f_unionized ≥ 0.8 | 0.5–0.8 | < 0.5 |

**Overall**: PASS = all optimal · BORDERLINE = ≥1 suboptimal, no poor · FAIL = ≥1 poor

> Ionization is **mandatory** — majority-ionised molecules fail automatically.
> pKa predicted via **Dimorphite-DL** (Ropp et al., 2019); overridden by SDF property `pKa`.
        """
    )
