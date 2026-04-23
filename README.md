<p align="center">
  <h1 align="center">BIOPREDICT</h1>
  <p align="center">
    <b>A QSAR Machine Learning Framework for Bioactivity Prediction of Query Compounds</b><br>
    <i>Cheminformatics · QSAR modelling · Random Forest · PubChem fingerprints · ChEMBL</i>
  </p>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Python-3.x-3776AB?style=flat-square&logo=python&logoColor=white" alt="Python">
  <img src="https://img.shields.io/badge/Jupyter-Notebook-F37626?style=flat-square&logo=jupyter&logoColor=white" alt="Jupyter">
  <img src="https://img.shields.io/badge/scikit--learn-ML-F7931E?style=flat-square&logo=scikit-learn&logoColor=white" alt="scikit-learn">
  <img src="https://img.shields.io/badge/RDKit-Cheminformatics-blue?style=flat-square" alt="RDKit">
  <img src="https://img.shields.io/badge/ChEMBL-Data%20Source-red?style=flat-square" alt="ChEMBL">
  <img src="https://img.shields.io/badge/PaDEL-Descriptor-green?style=flat-square" alt="PaDEL">
  <img src="https://img.shields.io/badge/License-MIT-lightgrey?style=flat-square" alt="License">
</p>

---

## Abstract

BIOPREDICT is an end-to-end quantitative structure–activity relationship (QSAR) modelling framework for predicting the bioactivity (pIC50) of small molecules against a target of interest. The pipeline queries bioactivity data directly from the ChEMBL database, performs exploratory data analysis with rigorous statistical validation, computes PubChem fingerprint descriptors using PaDEL-Descriptor, and trains a Random Forest regression model to predict pIC50 values for query compounds. The framework includes a comprehensive model comparison module using LazyPredict to benchmark the Random Forest against the full suite of available scikit-learn regressors. The trained model is exported as a serialised pickle object for downstream deployment and prediction on novel compounds.

The current implementation demonstrates the pipeline against a **Receptor Tyrosine Kinase (RTK)** target retrieved from ChEMBL, with IC50-based bioactivity classification and pIC50 prediction.

---

## Table of Contents

- [Key Features](#key-features)
- [Pipeline Overview](#pipeline-overview)
- [Repository Structure](#repository-structure)
- [Data Description](#data-description)
- [Bioactivity Classification Thresholds](#bioactivity-classification-thresholds)
- [Descriptors and Feature Engineering](#descriptors-and-feature-engineering)
- [Machine Learning Model](#machine-learning-model)
- [Statistical Analysis](#statistical-analysis)
- [Output Files](#output-files)
- [Installation and Dependencies](#installation-and-dependencies)
- [Usage](#usage)
- [Results](#results)
- [Author](#author)

---

## Key Features

| Feature | Description |
|---|---|
| **Automated data retrieval** | Queries IC50 bioactivity data from ChEMBL via REST API |
| **Activity classification** | Labels compounds as active / inactive / intermediate based on IC50 thresholds |
| **Lipinski descriptor calculation** | MW, LogP, NumHDonors, NumHAcceptors computed via RDKit |
| **PubChem fingerprints** | 881-bit PubChem fingerprint descriptors via PaDEL-Descriptor |
| **Statistical validation** | Mann-Whitney U test comparing active vs inactive compound properties |
| **Low-variance feature removal** | VarianceThreshold-based feature selection before model training |
| **Random Forest QSAR model** | Trained on pIC50 with 80/20 train/test split |
| **Comprehensive model comparison** | LazyPredict benchmarks 40+ regression algorithms side-by-side |
| **Model serialisation** | Trained model exported as `.pkl` for reuse and deployment |
| **Diagnostic visualisations** | Scatter plots (experimental vs predicted pIC50), box plots, R² and RMSE bar charts |

---

## Pipeline Overview

The workflow is implemented as six sequential Jupyter notebooks:

```
ChEMBL Database
      │
      ▼
┌─────────────────────────────────────────────────────────────────────┐
│  Notebook 1 — Data Collection                                       │
│  • Query ChEMBL by target name                                      │
│  • Retrieve IC50 bioactivity records                                │
│  • Remove duplicates and missing values                             │
│  • Classify compounds: active / inactive / intermediate             │
│  Output: 1_bioactivity_data_raw.csv → 3_bioactivity_data_curated.csv│
└──────────────────────────┬──────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│  Notebook 2 — Exploratory Data Analysis                             │
│  • Lipinski descriptors (RDKit): MW, LogP, HBD, HBA                │
│  • IC50 → pIC50 conversion (cap at 100 μM)                         │
│  • Distribution plots, scatter plots, box plots                     │
│  • Mann-Whitney U test: active vs inactive (5 descriptors)          │
│  Output: 4_bioactivity_data_3class_pIC50.csv                        │
│          5_bioactivity_data_2class_pIC50.csv                        │
│          mannwhitneyu_*.csv                                         │
└──────────────────────────┬──────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│  Notebook 3 — Descriptor Calculation                                │
│  • Export SMILES to molecule.smi                                    │
│  • PaDEL-Descriptor: PubChem fingerprints (881 bits, via padel.sh)  │
│  • Combine fingerprints with pIC50 labels                           │
│  Output: descriptors_output.csv                                     │
│          6_bioactivity_data_3class_pIC50_pubchem_fp.csv             │
└──────────────────────────┬──────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│  Notebook 4 — Model Building                                        │
│  • VarianceThreshold feature selection (threshold = 0.16)           │
│  • 80/20 train/test split                                           │
│  • Random Forest Regressor (n_estimators=100)                       │
│  • R² evaluation on test set                                        │
│  • Experimental vs predicted pIC50 scatter plot                     │
└──────────────────────────┬──────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│  Notebook 5 — Model Deployment                                      │
│  • Final model: Random Forest (n_estimators=500, random_state=42)   │
│  • VarianceThreshold feature selection (threshold = 0.1)            │
│  • Trained on full dataset                                          │
│  • MSE and R² evaluation                                            │
│  • Serialised model export: RTK_model.pkl                           │
│  • descriptor_list.csv: selected feature set for prediction         │
└──────────────────────────┬──────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│  Notebook 6 — Model Comparison                                      │
│  • LazyPredict: benchmarks 40+ sklearn regressors automatically     │
│  • Reports R², RMSE, and runtime for each algorithm                 │
│  • Bar chart visualisations (R², RMSE, computation time)            │
│  Output: model comparison PDFs                                      │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Repository Structure

```
BIOPREDICT/
│
├── 1-Data-Collection.ipynb                      # ChEMBL query, IC50 retrieval, activity labelling
├── 2-Exploratory-Data-Analysis.ipynb            # Lipinski descriptors, pIC50, EDA, Mann-Whitney U
├── 3-Descriptor-Calculation.ipynb               # PaDEL PubChem fingerprints
├── 4-Model_Building.ipynb                       # Random Forest model training and evaluation
├── 5-Model-deployment.ipynb                     # Final model, pickle export, deployment
├── 6-Model-Comparison.ipynb                     # LazyPredict multi-algorithm benchmark
│
├── padel.sh                                     # Shell script to run PaDEL-Descriptor via Java
│
├── 1_bioactivity_data_raw.csv                   # Raw IC50 records from ChEMBL
├── 2_bioactivity_data_preprocessed.csv          # Deduplicated: ChEMBL ID, SMILES, IC50
├── 3_bioactivity_data_curated.csv               # + bioactivity class labels
├── 4_bioactivity_data_3class_pIC50.csv          # + Lipinski descriptors + pIC50 (3 classes)
├── 5_bioactivity_data_2class_pIC50.csv          # Active / inactive only (intermediates removed)
├── 6_bioactivity_data_3class_pIC50_pubchem_fp.csv # PubChem fingerprints + pIC50
│
├── molecule.smi                                 # SMILES input file for PaDEL-Descriptor
├── descriptors_output.csv                       # Full 881-bit PubChem fingerprint matrix
├── descriptor_list.csv                          # Variance-filtered feature subset used for prediction
│
├── RTK_model.pkl                                # Serialised trained Random Forest model
│
├── mannwhitneyu_pIC50.csv                       # Mann-Whitney U: pIC50 (active vs inactive)
├── mannwhitneyu_MW.csv                          # Mann-Whitney U: Molecular Weight
├── mannwhitneyu_LogP.csv                        # Mann-Whitney U: LogP
├── mannwhitneyu_NumHDonors.csv                  # Mann-Whitney U: H-bond donors
├── mannwhitneyu_NumHAcceptors.csv               # Mann-Whitney U: H-bond acceptors
│
└── User_Manual.pdf                              # User manual
```

---

## Data Description

Bioactivity data is retrieved from the **ChEMBL database** using the `chembl_webresource_client` Python library.

| Parameter | Value |
|---|---|
| **Target** | Receptor Tyrosine Kinase (RTK) |
| **Activity type** | IC50 (half-maximal inhibitory concentration) |
| **Data source** | ChEMBL database |
| **Query field** | `standard_type = "IC50"` |
| **Columns retained** | `molecule_chembl_id`, `canonical_smiles`, `standard_value` |
| **Preprocessing** | Rows with missing `standard_value` or `canonical_smiles` removed; duplicate SMILES dropped |

---

## Bioactivity Classification Thresholds

Compounds are classified into three bioactivity categories based on IC50 values (in nM):

| Class | IC50 threshold | Interpretation |
|---|---|---|
| **Active** | ≤ 1,000 nM | Potent inhibitors |
| **Intermediate** | 1,000 – 10,000 nM | Moderate activity; excluded from binary model training |
| **Inactive** | ≥ 10,000 nM | Non-inhibitors |

IC50 values are capped at **100,000,000 nM** (100 μM) before conversion to pIC50 to prevent infinite or extreme values.

**IC50 → pIC50 conversion:**

```
pIC50 = −log₁₀(IC50 in M)
```

---

## Descriptors and Feature Engineering

### Lipinski Descriptors (RDKit)

Computed in Notebook 2 for exploratory analysis:

| Descriptor | Definition |
|---|---|
| `MW` | Molecular weight (Da) |
| `LogP` | Octanol–water partition coefficient |
| `NumHDonors` | Number of hydrogen bond donors |
| `NumHAcceptors` | Number of hydrogen bond acceptors |

### PubChem Fingerprints (PaDEL-Descriptor)

Computed in Notebook 3 via `padel.sh`:

- **Fingerprint type:** PubChem fingerprints (`PubchemFingerprinter.xml`)
- **Dimensionality:** 881 binary bits per compound
- **Preprocessing:** Salt removal (`-removesalt`), nitro group standardisation (`-standardizenitro`)
- **Tool:** [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) (Java, ≥ JRE 1.8)

```bash
# padel.sh — executed internally by Notebook 3
java -Xms1G -Xmx1G -Djava.awt.headless=true \
  -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar \
  -removesalt -standardizenitro -fingerprints \
  -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml \
  -dir ./ -file descriptors_output.csv
```

### Feature Selection

Low-variance features are removed using `sklearn.feature_selection.VarianceThreshold`:

| Notebook | Threshold | Use |
|---|---|---|
| Notebook 4 (training) | 0.16 (≈ 80/20 variance) | Model building and test evaluation |
| Notebook 5 (deployment) | 0.10 | Final production model and `descriptor_list.csv` |

The selected feature names are saved to `descriptor_list.csv` and must be used consistently when predicting on new compounds.

---

## Machine Learning Model

### Algorithm

**Random Forest Regressor** (`sklearn.ensemble.RandomForestRegressor`)

| Parameter | Notebook 4 (exploration) | Notebook 5 (deployment) |
|---|---|---|
| `n_estimators` | 100 | 500 |
| `random_state` | — | 42 |
| Training set | 80% | Full dataset |
| Test set | 20% | Full dataset (training evaluation) |

### Prediction target

`pIC50` — continuous regression output. Higher pIC50 indicates greater predicted potency.

### Model file

The final trained model is serialised using Python's `pickle` module:

```python
import pickle
model = pickle.load(open('RTK_model.pkl', 'rb'))
predictions = model.predict(X_new)
```

> **Note:** Input features for prediction must match the descriptor columns in `descriptor_list.csv` and be computed using the same PaDEL PubChem fingerprint settings as the training data.

---

## Statistical Analysis

### Mann-Whitney U Test

Non-parametric statistical comparison between **active** and **inactive** compound classes is performed in Notebook 2 for five descriptors. Results are saved individually as CSV files:

| File | Descriptor tested | Null hypothesis |
|---|---|---|
| `mannwhitneyu_pIC50.csv` | pIC50 | Active and inactive have the same pIC50 distribution |
| `mannwhitneyu_MW.csv` | Molecular weight | Same MW distribution |
| `mannwhitneyu_LogP.csv` | LogP | Same LogP distribution |
| `mannwhitneyu_NumHDonors.csv` | H-bond donors | Same HBD distribution |
| `mannwhitneyu_NumHAcceptors.csv` | H-bond acceptors | Same HBA distribution |

Significance threshold: **α = 0.05**. Results report the test statistic, p-value, and interpretation (reject / fail to reject H₀).

---

## Output Files

| File | Description |
|---|---|
| `1_bioactivity_data_raw.csv` | Raw ChEMBL IC50 records |
| `2_bioactivity_data_preprocessed.csv` | Cleaned: ChEMBL ID, SMILES, IC50 |
| `3_bioactivity_data_curated.csv` | + activity class labels (active/inactive/intermediate) |
| `4_bioactivity_data_3class_pIC50.csv` | + Lipinski descriptors + pIC50 (3 classes) |
| `5_bioactivity_data_2class_pIC50.csv` | Active/inactive only; intermediates removed |
| `6_bioactivity_data_3class_pIC50_pubchem_fp.csv` | PubChem fingerprints (881-bit) + pIC50 |
| `descriptors_output.csv` | Full PaDEL output before feature selection |
| `descriptor_list.csv` | Variance-filtered feature names for deployment |
| `molecule.smi` | SMILES file used as PaDEL input |
| `RTK_model.pkl` | Serialised Random Forest model (pickle) |
| `mannwhitneyu_*.csv` | Mann-Whitney U test results (5 files) |

---

## Installation and Dependencies

### Python packages

Install all required packages via pip:

```bash
pip install chembl_webresource_client rdkit pandas numpy scikit-learn seaborn matplotlib lazypredict
```

Or using conda:

```bash
conda install -c conda-forge rdkit
pip install chembl_webresource_client lazypredict scikit-learn seaborn matplotlib pandas numpy
```

| Package | Purpose |
|---|---|
| `chembl_webresource_client` | ChEMBL database query |
| `rdkit` | Lipinski descriptor calculation, SMILES processing |
| `pandas` / `numpy` | Data manipulation |
| `scikit-learn` | Random Forest, VarianceThreshold, train/test split, metrics |
| `seaborn` / `matplotlib` | Visualisation |
| `lazypredict` | Multi-algorithm model comparison (Notebook 6) |
| `pickle` | Model serialisation |

### PaDEL-Descriptor

PaDEL-Descriptor requires **Java Runtime Environment (JRE) ≥ 1.8**.

1. Download [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) and place the `PaDEL-Descriptor/` directory in the repository root.
2. Install Java:

   ```bash
   # Ubuntu / Debian
   sudo apt install -y default-jre

   # macOS (Homebrew)
   brew install --cask temurin
   ```

3. Verify Java installation:

   ```bash
   java -version
   ```

4. Update the Java path in `padel.sh` if your JRE is not at the default location:

   ```bash
   # Example path — adjust to your system
   /home/rashid/Downloads/jre-8u333-linux-x64/jre1.8.0_333/bin/java \
     -Xms1G -Xmx1G ...
   ```

---

## Usage

Run the notebooks in sequential order:

```bash
jupyter notebook
```

| Step | Notebook | Action |
|---|---|---|
| 1 | `1-Data-Collection.ipynb` | Query ChEMBL, retrieve IC50, classify compounds |
| 2 | `2-Exploratory-Data-Analysis.ipynb` | Lipinski EDA, pIC50 conversion, Mann-Whitney U test |
| 3 | `3-Descriptor-Calculation.ipynb` | Run PaDEL fingerprints via `padel.sh` |
| 4 | `4-Model_Building.ipynb` | Train and evaluate Random Forest on 80/20 split |
| 5 | `5-Model-deployment.ipynb` | Train final model, export `RTK_model.pkl` |
| 6 | `6-Model-Comparison.ipynb` | Benchmark all sklearn regressors with LazyPredict |

### Adapting to a different target

To run BIOPREDICT against a different ChEMBL target:

1. In `1-Data-Collection.ipynb`, modify the target search query:

   ```python
   target_query = target.search('YOUR_TARGET_NAME')
   ```

2. Select the appropriate target from the results dataframe and update the index accordingly.
3. Re-run all notebooks in order. All intermediate CSV files will be regenerated automatically.

### Predicting on new compounds

To predict pIC50 for new query compounds using the trained model:

1. Prepare a SMILES file of your query compounds.
2. Run PaDEL-Descriptor using `padel.sh` to generate PubChem fingerprints.
3. Load `descriptor_list.csv` to align the feature columns:

```python
import pickle, pandas as pd

# Load model and descriptor list
model = pickle.load(open('RTK_model.pkl', 'rb'))
desc_list = pd.read_csv('descriptor_list.csv')

# Load your query fingerprints and align columns
query_fp = pd.read_csv('your_query_descriptors.csv')
query_fp = query_fp[desc_list.columns]

# Predict pIC50
predictions = model.predict(query_fp)
print(predictions)
```

---

## Results

### Model performance summary

| Metric | Value |
|---|---|
| Algorithm | Random Forest Regressor |
| Fingerprint type | PubChem (881 bits) |
| Feature selection | VarianceThreshold (threshold = 0.1) |
| Training set | Full dataset (Notebook 5) |
| Evaluation metric | R², MSE |

Detailed performance metrics and comparison plots are generated within Notebooks 4, 5, and 6.

---

## Author

**Rashid Hussain**, PhD, RSci, MRSC  
Postdoctoral Researcher in Computational Pathology  
Humanitas Research Hospital (IRCCS), Milan, Italy  
[rashid.bioinfo@gmail.com](mailto:rashid.bioinfo@gmail.com) · [https://rashid-bioinfo.github.io](https://rashid-bioinfo.github.io)

---

<p align="center">
  <i>QSAR · Cheminformatics · Machine Learning for Drug Discovery</i>
</p>
