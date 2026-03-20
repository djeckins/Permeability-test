from __future__ import annotations

import io
from pathlib import Path

import pandas as pd
import streamlit as st

from epidermal_barrier_screen.io import parse_input
from epidermal_barrier_screen.screen import screen_records

st.set_page_config(page_title='Epidermal Barrier Screen', page_icon='🧪', layout='wide')

st.title('🧪 Epidermal Barrier Screen')
st.caption('GitHub-ready Streamlit app for screening molecules against epidermal barrier passage criteria.')

with st.sidebar:
    st.header('Input mode')
    mode = st.radio('Choose one', ['smiles', 'smiles_list', 'sdf', 'sdf_zip'])
    st.markdown(
        """
        **Supported inputs**
        - `smiles`: one SMILES
        - `smiles_list`: one SMILES per line, optional name after the SMILES
        - `sdf`: single SDF file
        - `sdf_zip`: ZIP archive containing one or more SDF files
        """
    )

payload = None
filename = None

if mode == 'smiles':
    payload = st.text_input('SMILES', placeholder='CC(=O)Oc1ccccc1C(=O)O')
elif mode == 'smiles_list':
    payload = st.text_area(
        'SMILES list',
        height=240,
        placeholder='CCO Ethanol\nCC(=O)Oc1ccccc1C(=O)O Aspirin',
    )
else:
    uploaded = st.file_uploader('Upload file', type=['sdf', 'zip'])
    if uploaded is not None:
        payload = uploaded.getvalue()
        filename = uploaded.name

run = st.button('Run screening', type='primary', use_container_width=True)

if run:
    if not payload:
        st.error('Please provide input first.')
        st.stop()

    records = parse_input(mode, payload, filename=filename)
    df = screen_records(records)

    st.subheader('Results')
    st.dataframe(df, use_container_width=True, hide_index=True)

    c1, c2, c3 = st.columns(3)
    c1.metric('Total records', len(df))
    c2.metric('Pass / pass with review', int(df['final_result'].isin(['pass', 'pass_with_manual_review']).sum()))
    c3.metric('Invalid', int((df['parse_status'] != 'ok').sum()))

    csv_bytes = df.to_csv(index=False).encode('utf-8')
    xlsx_buffer = io.BytesIO()
    df.to_excel(xlsx_buffer, index=False)
    xlsx_buffer.seek(0)

    dl1, dl2 = st.columns(2)
    dl1.download_button('Download CSV', data=csv_bytes, file_name='epidermal_barrier_results.csv', mime='text/csv', use_container_width=True)
    dl2.download_button('Download XLSX', data=xlsx_buffer.getvalue(), file_name='epidermal_barrier_results.xlsx', mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet', use_container_width=True)

with st.expander('Criteria used by the app'):
    st.markdown(
        """
        - **MW**: optimal `< 300`, suboptimal `300–500`
        - **logP / logD 7.4**: optimal `1–3`, suboptimal `0.5–1` or `3–5`
        - **TPSA**: optimal `< 60`, suboptimal `60–130`
        - **HBD**: optimal `0–3`, suboptimal `4–5`
        - **HBA**: optimal `2–8`, suboptimal `8–10`
        - **RotB**: optimal `< 10`, suboptimal `10–15`
        - **HAC**: optimal `< 30`, suboptimal `30–50`
        - **Formal charge**: optimal `0`, suboptimal `±1`
        - **pKa / ionization**: optional if supplied via SDF properties
        """
    )
