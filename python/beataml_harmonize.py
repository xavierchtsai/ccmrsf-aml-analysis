"""
Harmonize the Beat AML cohort (Bottomly et al., Cancer Cell 2022) for external validation.

Input:  data/beat_aml_lab.xlsx (sheet "summary", Supplementary Table S1 of Bottomly et al.)
Output: outputs/beataml_harmonized_51feat.csv (one row per patient)

Cohort: acute myeloid leukemia at inclusion, de novo, standard intensive induction
chemotherapy, valid overall survival. One specimen per patient: the first row passing the
filter in sheet order (proxy for the diagnostic specimen); mutations are not pooled across
specimens because later specimens can carry mutations acquired after diagnosis.

Features: 31 gene mutations and FLT3-TKD from variantSummary, FLT3-ITD from the FLT3-ITD
column, 13 cytogenetic abnormalities and complex karyotype from the karyotype string, and
clinical variables. Missing hematologic values are imputed with the cohort median.
ELN-2022 is derived with the rules in eln2022_row().

R/00_setup.R then (i) excludes APL, (ii) restricts KMT2A to partial tandem duplications,
(iii) calls CEBPA bZIP in-frame mutations from variantSummary and re-derives ELN-2022 with
them, and (iv) adds bone marrow blasts and complex karyotype from the same specimen.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import re
import numpy as np
import pandas as pd


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.environ.get("AML_DATA") or os.path.join(ROOT, "data")
OUT = os.environ.get("AML_OUT") or os.path.join(ROOT, "outputs")
os.makedirs(OUT, exist_ok=True)

# ---------- 51-feature constants (model feature set) ----------
KARYOTYPE_FEATURES = ['inv(3)', '-7', '-5', 't(6;9)', 't(3;3)', '-17', 't(16;16)', 'inv(16)',
                       't(9;22)', 't(8;21)', 'del(5q)', 't(9;11)', 't(8;16)']
MUTATION_FEATURES = ['PTPN11', 'NRAS', 'KRAS', 'NPM1', 'FLT_TKD', 'MLL', 'KIT', 'RUNX1', 'WT1',
                      'ASXL1', 'IDH1', 'IDH2', 'TET2', 'DNMT3A', 'TP53', 'SF3B1', 'U2AF1', 'SRSF2',
                      'ZRSR2', 'GATA2', 'STAG1', 'STAG2', 'Rad21', 'SMC1A', 'SMC3', 'PHF6', 'CBL-b',
                      'c-CBL', 'ETV6', 'EZH2', 'BCOR']
BASE_FEATURES = ['Age', 'Gender']
HEME_FEATURES = ['WBC', 'Hb', 'Plt', 'PB_Blast']
FEATURES_51 = BASE_FEATURES + HEME_FEATURES + KARYOTYPE_FEATURES + MUTATION_FEATURES + ['FLT3_ITD']

# genes present (as plain gene symbols) directly in variantSummary tokens, no remap needed
SAME_GENE = ['PTPN11', 'NRAS', 'KRAS', 'NPM1', 'KIT', 'RUNX1', 'WT1', 'ASXL1', 'IDH1', 'IDH2',
             'TET2', 'DNMT3A', 'TP53', 'SF3B1', 'U2AF1', 'SRSF2', 'ZRSR2', 'GATA2', 'STAG1',
             'STAG2', 'SMC1A', 'SMC3', 'PHF6', 'ETV6', 'EZH2', 'BCOR']
# variantSummary token -> our feature name, for symbols that differ
GENE_MAP = {'MLL': 'MLL', 'KMT2A': 'MLL', 'RAD21': 'Rad21', 'CBL': 'c-CBL', 'CBLB': 'CBL-b'}
FLT3_TKD_TOKENS = {'FLT3-TKD', 'FLT3-D835', 'FLT3-N676K', 'FLT3-V491L'}


def variant_tokens(variant_summary):
    """Gene-level tokens from a pipe-separated variantSummary cell, e.g.
    'ASXL1 (p.G646fs*12; MAF 30%)|NPM1|FLT3-ITD' -> {'ASXL1','NPM1','FLT3-ITD'}."""
    if pd.isna(variant_summary):
        return set()
    return {re.split(r'[\s(]', tok.strip())[0] for tok in str(variant_summary).split('|') if tok.strip()}


def mutation_flags(variant_summary):
    toks = variant_tokens(variant_summary)
    rec = {g: (1 if g in toks else 0) for g in SAME_GENE}
    for g in ['MLL', 'Rad21', 'c-CBL', 'CBL-b']:
        rec[g] = 0
    for tok in toks:
        if tok in GENE_MAP:
            rec[GENE_MAP[tok]] = 1
    rec['FLT_TKD'] = 1 if (toks & FLT3_TKD_TOKENS) else 0
    return rec


def kary_flags(karyotype):
    """13 ELN-2022 cytogenetic lesion flags + complex, from free-text karyotype.
    Unchanged from the earlier version of this script."""
    k = str(karyotype) if pd.notna(karyotype) else ''

    def has(pat):
        return 1 if re.search(pat, k) else 0

    f = {
        't(8;21)': has(r't\(8;21'),
        'inv(16)': has(r'inv\(16'),
        't(16;16)': has(r't\(16;16'),
        't(9;11)': has(r't\(9;11'),
        't(6;9)': has(r't\(6;9'),
        't(9;22)': has(r't\(9;22'),
        't(8;16)': has(r't\(8;16'),
        'inv(3)': has(r'inv\(3\b|inv\(3\)'),
        't(3;3)': has(r't\(3;3'),
        '-5': has(r'(^|[,\s\]])-5(?![0-9])'),
        'del(5q)': has(r'del\(5'),
        '-7': has(r'(^|[,\s\]])-7(?![0-9])'),
        '-17': has(r'(^|[,\s\]])-17|del\(17|add\(17|abn\(17|i\(17'),
    }
    body = re.sub(r'\[\d+\]', '', k)
    abn = len(re.findall(r't\(|inv\(|del\(|add\(|der\(|dup\(|[+\-]\d', body))
    f['complex'] = 1 if abn >= 3 else 0
    return f


def eln2022_row(rec):
    """ELN-2022 risk (0=adverse,1=intermediate,2=favorable), ported from
    eln2022_derive.py's eln2022_row(), adapted to this script's
    rec dict (values are already 0/1 ints, no row.get() needed)."""
    def g(col):
        return bool(rec.get(col, 0))

    MDS = ['ASXL1', 'BCOR', 'EZH2', 'RUNX1', 'SF3B1', 'SRSF2', 'STAG2', 'U2AF1', 'ZRSR2']
    cbf = g('t(8;21)') or g('inv(16)') or g('t(16;16)')
    if cbf:
        return 2
    adverse_cyto = (g('t(6;9)') or g('t(9;22)') or g('t(8;16)') or g('inv(3)') or g('t(3;3)') or
                    g('t(v;11)') or g('-5') or g('del(5q)') or g('-7') or g('-17') or
                    g('complex') or g('TP53'))
    if adverse_cyto:
        return 0
    if g('t(9;11)'):
        return 1
    if g('NPM1'):
        return 1 if g('ITD') else 2
    if g('CEBPA_bi'):
        return 2
    if g('ITD'):
        return 1
    if any(g(m) for m in MDS):
        return 0
    return 1


def main():
    s1 = pd.read_excel(f"{DATA}/beat_aml_lab.xlsx", sheet_name='summary', header=0)
    if 'dbgap_subject_id' not in s1.columns:          # some copies have a title row above the header
        s1 = pd.read_excel(f"{DATA}/beat_aml_lab.xlsx", sheet_name='summary', header=1)
    s1 = s1.dropna(subset=['dbgap_dnaseq_sample']).drop_duplicates('dbgap_dnaseq_sample')

    rows = []
    for _, r in s1.iterrows():
        rec = mutation_flags(r['variantSummary'])
        rec.update(kary_flags(r['karyotype']))
        rec['Age'] = r['ageAtDiagnosis']
        rec['Gender'] = 1 if str(r['consensus_sex']).strip().lower().startswith('m') else 0
        rec['WBC'] = r['wbcCount'] * 1000 if pd.notna(r['wbcCount']) else np.nan  # 10^9/L -> /uL
        rec['Hb'] = r['hemoglobin']
        rec['Plt'] = r['plateletCount']
        rec['PB_Blast'] = r['%.Blasts.in.PB']
        rec['ITD'] = 1 if str(r['FLT3-ITD']).strip().lower() in ('positive', 'yes', 'true', '1') else 0
        rec['CEBPA_bi'] = 1 if str(r['CEBPA_Biallelic']).strip().lower() in ('positive', 'yes', 'true', '1') else 0
        rec['OS_days'] = r['overallSurvival']
        vs = str(r['vitalStatus']).strip().lower()
        rec['event'] = 1 if vs == 'dead' else (0 if vs == 'alive' else np.nan)
        rec['isDenovo'] = str(r['isDenovo']).strip().upper()
        rec['dx'] = str(r['dxAtInclusion'])
        rec['induction'] = str(r['typeInductionTx'])
        rec['ELN2017'] = str(r['ELN2017'])
        rec['subject'] = r['dbgap_subject_id']
        rows.append(rec)
    bf = pd.DataFrame(rows)

    # ---------- cohort filter (confirmed 2026-07-12): initial-Dx AML, de novo, standard
    # intensive induction chemo, valid OS -- consistent with NTUH/NCRI/AMLSG eligibility ----------
    keep = (bf['dx'].str.contains('ACUTE MYELOID', case=False, na=False) &
            (bf['induction'] == 'Standard Chemotherapy') &
            (bf['isDenovo'] == 'TRUE') &
            (~bf['ELN2017'].isin(['NonInitial', 'NonAML'])) &
            (bf['OS_days'] > 0) & bf['event'].notna())
    bf = bf[keep].copy()
    # one specimen per patient: first row surviving the filter, in original sheet order
    # (proxy for the initial-diagnosis specimen) -- NOT a union across specimens
    bf = bf.drop_duplicates('subject').reset_index(drop=True)

    # impute remaining missing labs with this (post-filter, post-dedup) cohort's own median
    for c in ['WBC', 'Hb', 'Plt', 'PB_Blast', 'Age']:
        bf[c] = pd.to_numeric(bf[c], errors='coerce')
        bf[c] = bf[c].fillna(bf[c].median())

    bf['FLT3_ITD'] = bf['ITD'].astype(float)
    bf['ELN'] = bf.apply(eln2022_row, axis=1).astype(int)
    bf['T'] = bf['OS_days'] / 365.25
    bf['E'] = bf['event'].astype(int)

    print(f"Beat AML usable n = {len(bf)}")
    print("ELN-2022 derived dist (0=adv,1=int,2=fav):", bf['ELN'].value_counts().sort_index().to_dict())
    print(f"deaths: {int(bf['E'].sum())} ({100 * bf['E'].mean():.0f}%) | median OS y {bf['T'].median():.2f}")
    missing = [c for c in FEATURES_51 if c not in bf.columns]
    print("missing 51-feature columns:", missing)
    assert not missing, "harmonized frame is missing model features"

    out_cols = FEATURES_51 + ['ELN', 'T', 'E', 'subject']
    bf[out_cols].to_csv(f"{OUT}/beataml_harmonized_51feat.csv", index=False)
    print(f"Wrote {OUT}/beataml_harmonized.csv")


if __name__ == '__main__':
    main()
