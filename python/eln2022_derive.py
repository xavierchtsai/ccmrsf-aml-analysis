"""
Phase B helper: derive ELN-2022 risk (0=adverse,1=intermediate,2=favorable) from
genetic/cytogenetic columns, VALIDATE the derivation against the author-provided
NCRI ELN-2022, then apply it to the SG cohort (which lacks an ELN-2022 column).

Hierarchical ELN-2022 logic (Döhner et al., Blood 2022), implemented with
precedence:
 1. Core-binding-factor [t(8;21) or inv(16)/t(16;16)]            -> Favorable
 2. Adverse cytogenetics or TP53 or complex karyotype           -> Adverse
    [t(6;9), t(9;22), t(8;16), inv(3)/t(3;3), KMT2A-r t(v;11),
     -5/del(5q), -7, -17, complex, TP53mut]
 3. t(9;11)/MLLT3::KMT2A  (precedence over adverse gene muts)   -> Intermediate
 4. NPM1 mutated:  FLT3-ITD -> Intermediate ; else              -> Favorable
 5. CEBPA bZIP (CEBPA_bi)                                       -> Favorable
 6. FLT3-ITD (NPM1 wild-type)                                   -> Intermediate
 7. Myelodysplasia-related mutation                             -> Adverse
    [ASXL1,BCOR,EZH2,RUNX1,SF3B1,SRSF2,STAG2,U2AF1,ZRSR2]
 8. otherwise                                                   -> Intermediate
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "outputs")
os.makedirs(OUT, exist_ok=True)
MDS = ['ASXL1', 'BCOR', 'EZH2', 'RUNX1', 'SF3B1', 'SRSF2', 'STAG2', 'U2AF1', 'ZRSR2']


def g(row, col):
    """Safe 0/1 getter."""
    v = row.get(col, 0)
    try:
        return int(v) == 1
    except Exception:
        return False


def eln2022_row(r):
    cbf = g(r, 't(8;21)') or g(r, 'inv(16)') or g(r, 't(16;16)')
    if cbf:
        return 2
    adverse_cyto = (g(r, 't(6;9)') or g(r, 't(9;22)') or g(r, 't(8;16)') or
                    g(r, 'inv(3)') or g(r, 't(3;3)') or g(r, 't(v;11)') or
                    g(r, '-5') or g(r, 'del(5q)') or g(r, '-7') or g(r, '-17') or
                    g(r, 'complex') or g(r, 'TP53'))
    if adverse_cyto:
        return 0
    if g(r, 't(9;11)'):
        return 1
    if g(r, 'NPM1'):
        return 1 if g(r, 'ITD') else 2
    if g(r, 'CEBPA_bi'):
        return 2
    if g(r, 'ITD'):
        return 1
    if any(g(r, m) for m in MDS):
        return 0
    return 1


def derive(df):
    return df.apply(eln2022_row, axis=1).astype(int).values


def main():
    log = []

    # --- NCRI validation: derive from raw, compare to author-provided ---
    ncri_raw = pd.read_csv(f"{ROOT}/data/NCRI_AMLSG/paper_aml_prognosis_updated.tsv",
                           sep=' ', quotechar='"', engine='python')
    ncri_pre = pd.read_csv(f"{ROOT}/data/NCRI_preprocessing_df.csv")
    # NCRI raw lacks single U2AF1; reconstruct from the two hotspot columns
    if 'U2AF1' not in ncri_raw.columns:
        u = np.zeros(len(ncri_raw), dtype=int)
        for c in ['U2AF1_p.S34', 'U2AF1_p.Q157']:
            if c in ncri_raw.columns:
                u = u | (ncri_raw[c].fillna(0).astype(int).values == 1)
        ncri_raw['U2AF1'] = u.astype(int)
    derived = derive(ncri_raw)
    provided = ncri_pre['ELN_2022'].astype(int).values
    agree = float(np.mean(derived == provided))
    log.append(f"NCRI ELN-2022 derivation vs author-provided: agreement = {agree:.3f}")
    # confusion
    cm = pd.crosstab(pd.Series(provided, name='provided'),
                     pd.Series(derived, name='derived'))
    log.append("NCRI confusion (rows=provided, cols=derived):\n" + cm.to_string())
    log.append(f"provided dist {pd.Series(provided).value_counts().sort_index().to_dict()}")
    log.append(f"derived  dist {pd.Series(derived).value_counts().sort_index().to_dict()}")

    # --- SG: derive ELN-2022 from raw, attach to preprocessed SG ---
    sg_raw = pd.read_csv(f"{ROOT}/data/NCRI_AMLSG/paper_full_data_validation.tsv",
                         sep=' ', quotechar='"', engine='python')
    sg_pre = pd.read_csv(f"{ROOT}/data/AMLSG_preprocessing_df.csv")
    assert len(sg_raw) == len(sg_pre)
    sg_eln = derive(sg_raw)
    sg_pre['ELN_2022_derived'] = sg_eln
    outp = f"{OUT}/SG_with_ELN2022.csv"
    sg_pre.to_csv(outp, index=False)
    log.append(f"SG derived ELN-2022 dist "
               f"{pd.Series(sg_eln).value_counts().sort_index().to_dict()} "
               f"(0=adv,1=int,2=fav)")
    log.append(f"SG ELN-2017 (provided) dist "
               f"{sg_pre['ELN_2017'].value_counts().sort_index().to_dict()}")
    # sanity: SG ELN-2022 vs ELN-2017 cross-tab
    cm2 = pd.crosstab(sg_pre['ELN_2017'].astype(int),
                      pd.Series(sg_eln, name='ELN2022_derived'))
    log.append("SG ELN-2017 vs derived ELN-2022:\n" + cm2.to_string())
    log.append(f"saved SG+ELN2022 -> {outp}")

    text = "\n".join(log)
    print(text)
    with open(f"{OUT}/eln2022_derive.log", 'w') as f:
        f.write(text + "\n")


if __name__ == '__main__':
    main()
