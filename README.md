# HIV Infant Reservoir Prediction — Immune Cell Correlates of Viral Persistence

> Longitudinal flow cytometry analysis of innate and adaptive immune cell populations in perinatally HIV-infected (HEI) and HIV-exposed uninfected (HEU) infants from Mozambique, using limma/voom and linear mixed-effects models to identify immune signatures of HIV reservoir size.

---

## Overview

This repository contains R-based pipelines for analyzing multi-parameter flow cytometry data from the **Towards AIDS Remission Approaches (TARA) cohort** — a longitudinal study of HIV-exposed infected (HEI) and HIV-exposed uninfected (HEU) infants from Maputo, Mozambique. Infants were enrolled at the first postnatal visit (age 1–2 months), prior to antiretroviral therapy (ART) initiation, and followed longitudinally through 18 months of age.

The central goal is to identify **immune cell phenotypic correlates of HIV reservoir size** (quantified by viral titre) and to characterize how specific innate and adaptive immune compartments differ between HEI and HEU infants over time. Each immune cell compartment is analyzed in a parallel, modular pipeline.

---

## Scientific Background

In perinatal HIV infection, early ART is recommended to limit reservoir seeding, but questions remain about which immune cell populations shape reservoir dynamics. This repository investigates:

1. **HEI vs. HEU differences** — Which cell subsets are differentially abundant or activated between HIV-infected and HIV-exposed-uninfected infants?
2. **Viral titre as a predictor** — How do individual immune cell subset frequencies associate with HIV reservoir size (viral titre)?
3. **Age effects** — How do immune populations evolve over the first 18 months of life, and does infection status modify this trajectory?
4. **Sex effects** — Are there sex-based differences in immune subset frequencies?
5. **Group effects (HEI vs. HEU)** — Which cell subsets remain persistently altered in HEI after ART?
6. **Reservoir prediction** — Can change in immune phenotype between time points predict change in reservoir size?

---

## Cohort

| Feature | Details |
|---|---|
| **Cohort** | TARA (Towards AIDS Remission Approaches) |
| **Location** | Maputo, Mozambique |
| **Groups** | HEI (HIV exposed, infected) and HEU (HIV exposed, uninfected) |
| **Enrollment** | First postnatal visit, age 1–2 months (pre-ART) |
| **Follow-up** | Longitudinal, up to 18 months of age |
| **Reservoir measure** | HIV viral titre (copies/mL), normalized and scaled |
| **Immunophenotyping** | Multi-parameter flow cytometry (FlowJo gating, normalized to CD45+ cells) |

---

## Repository Structure

```
HIV_Infants_Reservoir_Prediction/
│
├── Cleaned_Datasets/                          # Normalized, CD45-scaled flow cytometry data + viral titres
├── Clinical_Dataset/                          # Demographic and clinical metadata
│
├── NK_Mixed_Model_Output/                     # NK cell model outputs and figures
├── Monocyte_Mixed_Model_Output/               # Monocyte model outputs and figures
├── DC_Mixed_Model_Output/                     # Dendritic cell model outputs and figures
├── T_B_Mixed_Model_Output/                    # T and B cell model outputs and figures
├── Total_Mixed_Model_Output/                  # Combined/total immune panel outputs
│
├── NK_Data_Normalization.R                    # Normalize NK flow data to CD45+ counts
├── NK_Limma.R                                 # Limma/voom + mixed-effects modeling for NK cells
│
├── Monocyte_Data_Normalization.R              # Normalize monocyte data
├── Monocytes_Limma.R                          # Limma/voom + mixed-effects modeling for monocytes
│
├── DC_Data_Normalization.R                    # Normalize dendritic cell data
├── DC_Limma.R                                 # Limma/voom + mixed-effects modeling for DCs
│
├── T_B_Cell_Data_Normalization.R              # Normalize T and B cell data
├── T_B_Cell_Limma.R                           # Limma/voom + mixed-effects modeling for T/B cells
│
├── Total_Limma.R                              # Combined analysis across all cell types
├── Demographics_Database_Cleaning.R           # Clean and merge demographic/clinical data
│
├── HIV_Infants_Mixed_Model_Analysis_Report.Rmd  # R Markdown analysis report
├── HIV_Infants_Mixed_Model_Analysis_Report.html # Rendered HTML report
│
├── HIV_Infants_Reservoir_Prediction.Rproj    # RStudio project file
└── README.md
```

---

## Analysis Pipeline

Each immune compartment follows the same modular two-step pipeline:

### Step 1 — Data Normalization (`*_Data_Normalization.R`)
- Reads raw flow cytometry Excel sheets
- Normalizes cell subset frequencies to total CD45⁺ live lymphocytes
- Exports cleaned, normalized CSVs to `Cleaned_Datasets/`

### Step 2 — Statistical Modeling (`*_Limma.R`)
Each Limma script performs two complementary analyses:

**A. Differential Abundance (limma/voom)**
- Transposes normalized flow data into a features × samples matrix
- Applies `voom` transformation to stabilize variance
- Accounts for repeated measures using `duplicateCorrelation` (blocking on `PID`)
- Design matrix: `~ Group + Viral Titre + Age + Sex`
- Compares HEI vs. HEU (`GroupHEI` coefficient)
- Outputs: volcano plots of top differentially abundant subsets (FDR < 0.05 highlighted)

**B. Linear Mixed-Effects Models (lmerTest)**
- Fits per-subset models: `subset ~ Viral Titre + Age + Sex + Group + (1 | PID)`
- Extracts fixed-effect coefficients for Viral Titre, Age, Sex, and Group
- Generates coefficient bar plots with error bars and significance stars for each effect
- Outputs: four coefficient plots per cell type (Viral Titre, Age, Sex, HEI group effect)

---

## Cell Compartments Analyzed

| Compartment | Script Pair | Key Subsets |
|---|---|---|
| **NK cells** | `NK_Data_Normalization.R` / `NK_Limma.R` | CD56bright, CD56dimCD16−, CD56dimCD16+, CD56−CD16+ (dysfunctional) |
| **Monocytes** | `Monocyte_Data_Normalization.R` / `Monocytes_Limma.R` | Classical (CM), Intermediate (IM), Non-classical (NCM) |
| **Dendritic Cells** | `DC_Data_Normalization.R` / `DC_Limma.R` | mDC, pDC, cDC subsets |
| **T & B Cells** | `T_B_Cell_Data_Normalization.R` / `T_B_Cell_Limma.R` | CD4, CD8 T cell subsets, B cell populations |
| **Total Panel** | `Total_Limma.R` | All subsets combined |

---

## Outputs Per Cell Type

Each compartment produces the following output files in its respective `*_Mixed_Model_Output/` folder:

| File | Description |
|---|---|
| `HEIvsHEU_<type>_Volcano.png` | Volcano plot — HEI vs HEU differential abundance |
| `Viral_Effect_<type>_Coefficient_Plot.png` | Effect of viral titre on each subset |
| `Age_Effect_<type>_Coefficient_Plot.png` | Effect of age on each subset |
| `Sex_Effect_<type>_Coefficient_Plot.png` | Effect of sex (Male vs. Female) on each subset |
| `HIV_Effect_<type>_Coefficient_Plot.png` | HEI vs. HEU group effect on each subset |

---

## Dependencies

All scripts are written in **R**. Required packages:

| Package | Purpose |
|---|---|
| `readxl` / `readr` | Data ingestion |
| `limma` | Differential analysis with voom transformation |
| `lmerTest` | Linear mixed-effects models |
| `broom.mixed` | Tidy extraction of mixed model results |
| `dplyr` | Data manipulation |
| `ggplot2` | Visualization |
| `ggrepel` | Non-overlapping labels on volcano plots |
| `pheatmap` | Heatmap visualizations |

Install via:

```r
install.packages(c("readxl", "readr", "dplyr", "ggplot2", "ggrepel", "pheatmap",
                   "lmerTest", "broom.mixed"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("limma")
```

---

## Usage

1. Clone the repository and open `HIV_Infants_Reservoir_Prediction.Rproj` in RStudio.
2. Update the input/output paths in each script:

```r
in.path  <- 'path/to/Cleaned_Datasets/'
out.path <- 'path/to/NK_Mixed_Model_Output/'
```

3. Run normalization scripts first to generate cleaned CSVs, then run the corresponding Limma scripts.
4. Render `HIV_Infants_Mixed_Model_Analysis_Report.Rmd` for a full interactive report.

**Recommended run order:**
```
Demographics_Database_Cleaning.R
→ *_Data_Normalization.R (one per cell type)
→ *_Limma.R (one per cell type)
→ Total_Limma.R
→ HIV_Infants_Mixed_Model_Analysis_Report.Rmd
```

> ⚠️ **Note:** Raw flow cytometry data and identifiable clinical data are not included in this repository due to participant privacy. Contact the study team for data access.

---

## Citation

If you use this code, please cite this repository and acknowledge the TARA cohort and the University of Miami Pahwa Laboratory.

---

## Contact

For questions or collaboration, please open an [issue](https://github.com/codeneeded/HIV_Infants_Reservoir_Prediction/issues) in this repository.
