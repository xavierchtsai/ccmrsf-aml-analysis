# CCM-RSF: analysis code

Code for the development and external validation of the **Clinical-Cytogenetic-Molecular Random Survival Forest (CCM-RSF)**, a survival model for adults with de novo, non-APL acute myeloid leukemia treated with intensive induction.

Web calculator: https://xavierchtsai.github.io/ccmrsf-calculator/

## Contents

| Path | Purpose |
|---|---|
| `python/eln2022_derive.py` | Derives ELN-2022 for AMLSG and checks the derivation against the NCRI ELN-2022 labels |
| `python/beataml_harmonize.py` | Harmonizes the Beat AML cohort |
| `R/00_setup.R` | Shared settings: paths, the 53 predictors, cohort loading and exclusions, model loading |
| `R/01_build_model.R` | Fits the random survival forest on the development cohort (needs NTUH data) |
| `R/02_table_s1_baseline.R` | Supplementary Table 1: baseline characteristics |
| `R/03_table1_ntuh.R` | Table 1: development cohort by ELN-2022 and CCM-RSF group |
| `R/04_table2a_discrimination.R` | Table 2A: c-index with 2,000 paired bootstrap resamples, hazard ratios, comparator models |
| `R/05_table2b_decomposition.R` | Table 2B: source of the discrimination gain; multivariable Cox model |
| `R/06_figures_table3.R` | Figures 1–2, calibration, permutation importance, Table 3 |
| `R/07_decision_curves.R` | Decision curve analysis |
| `R/08_sankey.R` | ELN-2022 to CCM-RSF reclassification |
| `R/09_transplant_table4.R` | Allogeneic HSCT landmark analysis, Figure 5, Supplementary Figure 3, Table 4 |
| `R/10_partial_dependence.R` | Partial dependence plots |
| `R/11_export_web_model.R` | Exports the forest for the web calculator |
| `R/run_all.R`, `R/run_all_windows.R` | Run scripts 01–10 in order |
| `model/rsf_model_NTUH_53feat.rds` | The frozen model (randomForestSRC object) |

## Model

- 400 trees, minimum terminal node size 10, `mtry = floor(sqrt(53))`, log-rank splitting with 10 random split points, fitted with randomForestSRC 3.6.2.
- 53 predictors at diagnosis: age, sex, WBC count, hemoglobin, platelet count, bone marrow blast percentage, 13 ELN-defined cytogenetic abnormalities, complex karyotype, 31 gene mutations (including *KMT2A* partial tandem duplication), *FLT3*-ITD and *CEBPA* bZIP in-frame mutation.
- Risk score: ensemble mortality (`predict()$predicted`). Risk groups use cutoffs locked at the tertiles of the development cohort: 52.86 and 97.52.

## Running the analysis

1. Place the input files in `data/` (see `data/README.md`).
2. From the repository root:
   ```
   python3 python/eln2022_derive.py
   python3 python/beataml_harmonize.py
   ```
3. From `R/`: `Rscript run_all.R` (on Windows, source `run_all_windows.R`).

Tables are written to `outputs/` and figures to `figures/`. Without the NTUH data, skip `01_build_model.R`; the validation scripts use the frozen model in `model/`, but `00_setup.R` reads the NTUH cohort for the development-cohort summaries, so those scripts need the NTUH files.

Software: R 4.5 with randomForestSRC, survival, survminer, ggplot2, ggtext, ggalluvial, glmnet, gbm, gridExtra, cowplot, dplyr, readxl and jsonlite; Python 3 with pandas, numpy and openpyxl.

## Data availability

- **NTUH development cohort**: not publicly available; available from the corresponding author on reasonable request with institutional review board approval.
- **UK-NCRI and AMLSG**: Tazi et al., Nat Commun 2022 (https://doi.org/10.5281/zenodo.6878209).
- **Beat AML**: Bottomly et al., Cancer Cell 2022 (Supplementary Table S1).

## Contact

Xavier Cheng-Hong Tsai, Division of Hematology, Department of Internal Medicine, National Taiwan University Hospital (xavierchtsai@ntu.edu.tw).

## License

Code: MIT License (see `LICENSE`).
