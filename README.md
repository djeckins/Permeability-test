# Epidermal Barrier Screen

A GitHub-ready molecular screening tool for estimating whether compounds fit a physicochemical window associated with epidermal barrier passage.

It supports four input modes:

- **single SMILES**
- **SMILES list**
- **single SDF**
- **ZIP archive with SDF files**

The app always **canonicalizes SMILES** and exports a result table with the original input, canonical SMILES, calculated descriptors, per-criterion statuses, and a final result.

## What it calculates

For each molecule the tool reports:

- molecular weight (MW)
- cLogP
- LogD 7.4 if provided in SDF properties, otherwise cLogP is used as a proxy
- TPSA
- HBD
- HBA
- Rotatable bonds (RotB)
- Heavy atom count (HAC)
- formal charge
- pKa if provided in SDF properties
- estimated fraction unionized at pH 5.5 when pKa is available
- HBD:HBA ratio
- criterion-by-criterion pass state
- final screening result

## Repository layout

```text
.
├── app.py
├── requirements.txt
├── runtime.txt
├── pyproject.toml
├── .gitignore
├── examples/
│   ├── example_smiles.txt
│   └── example_smiles_list.txt
└── src/
    └── epidermal_barrier_screen/
        ├── __init__.py
        ├── descriptors.py
        ├── io.py
        ├── screen.py
        └── cli.py
```

## Local run

Create a virtual environment and install the dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Run the Streamlit interface:

```bash
streamlit run app.py
```

Run the command-line interface:

```bash
python -m epidermal_barrier_screen.cli --mode smiles --input "CC(=O)Oc1ccccc1C(=O)O" --output-prefix results/aspirin
python -m epidermal_barrier_screen.cli --mode smiles_list --input examples/example_smiles_list.txt --output-prefix results/list_run
python -m epidermal_barrier_screen.cli --mode sdf --input path/to/file.sdf --output-prefix results/sdf_run
python -m epidermal_barrier_screen.cli --mode sdf_zip --input path/to/archive.zip --output-prefix results/zip_run
```

The CLI writes **both CSV and XLSX**.

## Streamlit Community Cloud deployment

1. Push this repository to GitHub.
2. In Streamlit Community Cloud, create a new app from the repo.
3. Use **`app.py`** as the main file.
4. Keep **Python 3.12** from `runtime.txt` for broad package compatibility.

## Input notes

### 1) single SMILES
Paste one SMILES string.

### 2) SMILES list
One molecule per line. You can optionally place a name after the SMILES.

Example:

```text
CCO Ethanol
CC(=O)Oc1ccccc1C(=O)O Aspirin
```

### 3) SDF
Upload one SDF file.

### 4) ZIP with SDF files
Upload a ZIP archive containing one or more `.sdf` files.

## Optional SDF properties

If your SDF contains these properties, the app will use them automatically:

- `pKa`, `PKA`, `predicted_pKa`, `input_pka`
- `LogD`, `LOGD`, `logD_7_4`, `input_logd_7_4`
- `_Name`, `Name`, `TITLE`, `ID`, `Compound`

## Example outputs

The exported result table includes:

- `input_smiles`
- `canonical_smiles`
- `mw`, `clogp`, `tpsa`, `hbd`, `hba`, `rotb`, `hac`, `formal_charge`
- `input_pka`, `input_logd_7_4`
- `fraction_unionized_pH5_5`
- one status column per criterion
- `final_result`

## Final result labels

- `pass`
- `pass_with_manual_review`
- `borderline`
- `enhancers_likely_needed`
- `invalid_input`

## Notes

- This is a **screening utility**, not a biological guarantee.
- `pKa` and experimental `LogD` are optional because they are not reliably derivable from plain SMILES using base RDKit alone.
- When `LogD 7.4` is missing, the app uses **cLogP as a proxy** and marks that choice in the exported table.
