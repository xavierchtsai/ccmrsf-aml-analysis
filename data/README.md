# Input data

Data files are not included in this repository. Place them in this folder with the names below.

| File | Cohort | Source | Used by |
|---|---|---|---|
| `NTUH_preprocessing_df_36.csv` | NTUH (n=801) | Institutional; not publicly available | 00_setup.R |
| `Cohort_20250527_with_RSF.xlsx` | NTUH | Institutional; FLT3-ITD, CEBPA bZIP, bone marrow blasts, complex karyotype, CR/relapse dates, HSCT in CR1 | 00_setup.R, 03, 09 |
| `AML raw data (NTUH+AMLSG+UKNCRI).xlsx` (sheet `total`) | NCRI, AMLSG | Compiled from Tazi et al. 2022; FLT3, ICC 2022, intensive treatment, HSCT in CR1, CEBPA bZIP, bone marrow blasts, complex karyotype | 00_setup.R |
| `NCRI_preprocessing_df.csv` | UK-NCRI (n=2,113 before exclusions) | Harmonized from Tazi et al. 2022 | 00_setup.R, eln2022_derive.py |
| `AMLSG_preprocessing_df.csv` | AMLSG (n=1,540 before exclusions) | Harmonized from Tazi et al. 2022 | eln2022_derive.py |
| `NCRI_AMLSG/paper_aml_prognosis_updated.tsv` | UK-NCRI | Tazi et al. 2022, Zenodo 10.5281/zenodo.6878209 | eln2022_derive.py |
| `NCRI_AMLSG/paper_full_data_validation.tsv` | AMLSG | Tazi et al. 2022, Zenodo 10.5281/zenodo.6878209 | eln2022_derive.py |
| `beat_aml_lab.xlsx` (sheet `summary`) | Beat AML | Bottomly et al., Cancer Cell 2022, Supplementary Table S1 | beataml_harmonize.py, 00_setup.R |

Notes:

- In the NCRI and AMLSG files, rows are matched to the `total` sheet by order within each source; `00_setup.R` checks the alignment on age.
- NCRI and AMLSG record bone marrow blasts but not peripheral blood blasts; the model therefore uses bone marrow blast percentage in all cohorts.
- Overall survival is measured from diagnosis and is taken from `NTUH_preprocessing_df_36.csv`; only HSCT status is taken from the NTUH master file.
