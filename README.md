# CCM-RSF analysis code

Code for the development and external validation of the CCM-RSF survival model in acute myeloid leukemia.

Web calculator: https://xavierchtsai.github.io/ccmrsf-calculator/

## Structure

- `python/`: data preparation for the AMLSG and Beat AML cohorts
- `R/`: model development, validation, tables and figures (numbered in run order)
- `model/`: the frozen model (`randomForestSRC` object)
- `data/`: input files (not included; see `data/README.md`)

## Usage

```
python3 python/eln2022_derive.py
python3 python/beataml_harmonize.py
cd R && Rscript run_all.R
```

## Data availability

The NTUH cohort is available from the corresponding author on reasonable request. The UK-NCRI and AMLSG data are from Tazi et al. (Nat Commun 2022), and the Beat AML data from Bottomly et al. (Cancer Cell 2022).

## License

The code in this repository is released under the [MIT License](LICENSE). The license covers the code and the frozen model file; it does not cover the patient data, which are not distributed here. If you use this code or the model, please cite this repository until the accompanying article is published.
