# QSAR_Tools

**Universal QSAR Pipeline for Virtual Screening of Molecular Libraries**

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Sbolivar16/QSAR_Tools/blob/main/QSAR_Universal_Interactive_EN.ipynb)
[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)
[![RDKit](https://img.shields.io/badge/RDKit-2023+-green.svg)](https://www.rdkit.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Overview

**QSAR_Tools** is an interactive Google Colab notebook that implements a complete **Quantitative Structure-Activity Relationship (QSAR)** pipeline ‚Äî from raw ChEMBL bioactivity data to a ranked list of molecular candidates ready for molecular docking.

The notebook is fully generalizable: by changing just **5 variables**, any researcher can adapt it to any therapeutic target available in ChEMBL, screen any molecular library in CSV format, and obtain publication-ready results.

> **Developed as part of an in silico drug discovery study targeting the NS2B-NS3 protease of Dengue virus (PDB: 2FOM) using Colombian natural products (PNDBCOL database).**

---

## What does this notebook do?

```
ChEMBL bioactivity data (IC50 / Ki / EC50)
        |
    Data curation & standardization
        |
    Molecular descriptor calculation (RDKit 2D + ECFP4)
        |
    QSAR model training (Random Forest / XGBoost / SVM)
        |
    5-fold cross-validation + external test set evaluation
        |
    Virtual screening of user-provided molecular library
        |
    Applicability Domain (Williams Plot)
        |
    SAR analysis with SHAP values
        |
    ADMET profiling (Lipinski, Veber, PAINS, QED)
        |
    Analog generation by scaffold modification
        |
    TOP N candidates + SMILES ready for AutoDock Vina
```

---

## Repository Structure

```
QSAR_Tools/
|
|-- QSAR_Universal_Interactive_EN.ipynb   # Main notebook (English, full pipeline)
|-- PASO0_Convertir_PNDBCOL_a_CSV.ipynb  # Utility: convert raw text DB to CSV
|-- README.md                             # This file
|-- LICENSE
```

---

## Quick Start

### Option 1 ‚Äî Google Colab (recommended, no installation needed)

1. Click the **Open in Colab** badge at the top
2. Run the installation cell (Module 0) ‚Äî takes ~4 min on first run
3. Edit the 5 configuration variables in the config cell
4. Run all cells sequentially
5. Download the results ZIP at the end

### Option 2 ‚Äî Local installation

```bash
git clone https://github.com/Sbolivar16/QSAR_Tools.git
cd QSAR_Tools

pip install rdkit chembl-webresource-client xgboost shap mols2grid \
            scikit-learn pandas numpy matplotlib seaborn

jupyter notebook QSAR_Universal_Interactive_EN.ipynb
```

---

## Configuration ‚Äî 5 Variables

The only section you need to edit before running:

```python
TARGET_ID     = 'CHEMBL613966'      # ChEMBL target ID
TARGET_NAME   = 'NS2B-NS3 Dengue'  # Name for figures
ACTIVITY_TYPE = 'IC50'              # IC50, Ki, EC50, or Kd
THRESHOLD_nM  = 10000               # Activity cutoff (nM)
TOP_N         = 10                  # Number of final candidates
```

### Ready-to-use targets

| Target | ChEMBL ID | Activity | Cutoff (nM) |
|--------|-----------|----------|-------------|
| Dengue NS2B-NS3 protease | CHEMBL613966 | IC50 | 10000 |
| HIV-1 Protease | CHEMBL612545 | IC50 | 1000 |
| SARS-CoV-2 Mpro | CHEMBL4523582 | IC50 | 10000 |
| P. falciparum DHFR | CHEMBL3718 | IC50 | 1000 |
| EGFR Kinase | CHEMBL203 | IC50 | 100 |
| COX-2 | CHEMBL230 | IC50 | 1000 |
| Leishmania TDR1 | CHEMBL1741203 | IC50 | 10000 |
| Trypanosoma TbCK1 | CHEMBL3885796 | IC50 | 10000 |

---

## Outputs

After running the full pipeline, a ZIP file is generated containing:

| File | Description |
|------|-------------|
| `01_dataset.png` | IC50 distribution and class balance |
| `02_models_comparison.png` | Cross-validation model comparison |
| `03_test_evaluation.png` | ROC curve and confusion matrix |
| `04_williams_plot.png` | Applicability Domain (Williams Plot) |
| `05_top_structures.png` | 2D structures of top candidates |
| `06_shap_SAR.png` | SAR analysis ‚Äî SHAP summary plot |
| `07_admet_qed.png` | Drug-likeness profile (QED) |
| `08_analogs.png` | Generated analogs from top hits |
| `TOP_candidates_docking.csv` | Ranked candidates with all properties |
| `full_screening_results.csv` | Complete screening results |
| `SMILES_for_Vina.smi` | SMILES ready for OpenBabel -> AutoDock Vina |
| `ADMET_profile.csv` | Full ADMET table |
| `analogs_generated.csv` | Analogs with predicted activities |
| `STUDY_SUMMARY.txt` | Parameters, metrics, and references |

---

## Methods

### Data Source
Bioactivity data retrieved from [ChEMBL](https://www.ebi.ac.uk/chembl/) using the official Python API (`chembl-webresource-client`).

### Curation Pipeline
1. Filter exact measurements (relation `=` only)
2. Remove NaN values
3. Validate SMILES with RDKit
4. Standardize structures (remove salts, largest fragment)
5. Deduplicate by molecule (median activity)
6. Assign binary labels: ACTIVE (<= threshold nM) / INACTIVE (> threshold nM)

### Descriptors
- **RDKit 2D descriptors** (~200 physicochemical properties)
- **Feature selection**: NaN filter (>20%), zero variance, Pearson correlation (r > 0.95)
- **Scaling**: StandardScaler (mean=0, std=1)

### Machine Learning Models

| Algorithm | Hyperparameters |
|-----------|----------------|
| Random Forest | n_estimators=300, class_weight='balanced' |
| XGBoost | n_estimators=300, max_depth=6, lr=0.05, scale_pos_weight=auto |
| SVM (RBF) | kernel='rbf', class_weight='balanced', probability=True |

### Validation
- Stratified 5-fold cross-validation
- Independent external test set (20% holdout)
- Metrics: ROC-AUC, MCC, Balanced Accuracy, Sensitivity, Specificity

### Applicability Domain
Leverage method (Williams Plot):
- Hat value: h = x(X'X)^-1 x'
- Threshold: h* = 3(k+1)/n

### SAR Analysis
SHAP (SHapley Additive exPlanations) using `TreeExplainer` (Random Forest / XGBoost) or `KernelExplainer` (SVM).

### ADMET Profile
- Lipinski Rule-of-Five (MW <= 500, logP <= 5, HBD <= 5, HBA <= 10)
- Veber rules (TPSA <= 140 A^2, RotBonds <= 10)
- PAINS structural alerts (RDKit FilterCatalog)
- QED ‚Äî Quantitative Estimate of Drug-likeness

### Analog Generation
Simplified Matched Molecular Pairs (MMP) strategy using RDKit reaction SMARTS:
- OH -> F (bioisostere)
- OH -> OMe (increased lipophilicity)
- Ph -> Pyridine (reduced logP)
- NH2 -> NHMe (improved cell permeability)
- Me -> CF3 (metabolic stability)

---

## Next Steps After This Notebook

This notebook is **Step 1** of a complete multi-scale in silico pipeline:

```
QSAR (this notebook)
    |
    v  Top N candidates selected
    |
Molecular Docking (AutoDock Vina / Glide)
    |  obabel -ismi SMILES_for_Vina.smi -O candidates.sdf --gen3d -h
    |  vina --receptor receptor.pdbqt --ligand ligand.pdbqt --out result.pdbqt
    |
    v  Top 3-5 docking poses
    |
Molecular Dynamics (AMBER / GROMACS)
    |  100-200 ns production run
    |  RMSD, RMSF, H-bond occupancy analysis
    |
    v
MM-PBSA / MM-GBSA
    |  Binding free energy calculation
    |  Per-residue energy decomposition
    |
    v  Final 1-2 candidates for publication
```

---

## Validation Metrics ‚Äî Interpretation Guide

| Metric | Good value | Interpretation |
|--------|-----------|----------------|
| ROC-AUC | >= 0.75 | Discriminative power (0.5 = random) |
| MCC | >= 0.40 | Best metric for imbalanced datasets |
| Balanced Accuracy | >= 0.70 | Unbiased accuracy |
| Sensitivity | >= 0.65 | Fraction of actives correctly identified |
| Specificity | >= 0.65 | Fraction of inactives correctly rejected |

---

## Dependencies

| Library | Version | Purpose |
|---------|---------|---------|
| `rdkit` | >= 2023.03 | Molecular descriptors, SMILES handling |
| `chembl-webresource-client` | >= 0.10 | ChEMBL API access |
| `scikit-learn` | >= 1.3 | ML models, validation |
| `xgboost` | >= 2.0 | Gradient boosting classifier |
| `shap` | >= 0.44 | SHAP explainability |
| `pandas` | >= 2.0 | Data manipulation |
| `numpy` | >= 1.24 | Numerical computing |
| `matplotlib` | >= 3.7 | Visualization |
| `seaborn` | >= 0.12 | Statistical plots |
| `mols2grid` | >= 2.0 | Interactive molecular grids |

---

## Citation

If you use this notebook in your research, please cite:

```bibtex
@software{QSAR_Tools_2025,
  author  = {Bolivar, S.},
  title   = {QSAR_Tools: Universal QSAR Pipeline for Virtual Screening},
  year    = {2025},
  url     = {https://github.com/Sbolivar16/QSAR_Tools},
  note    = {Google Colab interactive notebook}
}
```

**Key references:**

- Mendez D, et al. (2019) ChEMBL: towards direct deposition of bioassay data. *Nucleic Acids Research* 47:D930-D940
- Pedregosa F, et al. (2011) Scikit-learn: Machine Learning in Python. *JMLR* 12:2825-2830
- Chen T & Guestrin C (2016) XGBoost: A Scalable Tree Boosting System. *KDD*
- Lundberg SM & Lee SI (2017) A unified approach to interpreting model predictions. *NeurIPS*
- Lipinski CA, et al. (1997) Experimental and computational approaches to estimate solubility and permeability. *Adv Drug Deliv Rev* 23:3-25
- Veber DF, et al. (2002) Molecular properties that influence oral bioavailability. *J Med Chem* 45:2615-2623
- Tropsha A, et al. (2003) The importance of being earnest: validation is the absolute essential for successful QSAR. *QSAR Comb Sci* 22:69-77
- Bickerton GR, et al. (2012) Quantifying the chemical beauty of drugs. *Nature Chemistry* 4:90-98

---

## License

This project is licensed under the MIT License ‚Äî see the [LICENSE](LICENSE) file for details.

---

## Contributing

Contributions are welcome. Please open an issue first to discuss what you would like to change.

---

## Author

**Sbolivar16**
*In silico drug discovery | Molecular docking | Molecular dynamics | QSAR*

> *This tool was developed to support research on anti-dengue compounds from Colombian natural products, but is designed to be universally applicable to any ChEMBL target.*

