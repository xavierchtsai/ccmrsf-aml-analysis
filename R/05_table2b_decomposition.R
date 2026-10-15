## Table 2B: where the discrimination gain comes from (ELN, age, ELN + age, RSF 3 groups, RSF),
## plus the multivariable Cox model (RSF risk per SD + ELN-2022 + age) described in Methods 2.3.
## External cohorts restricted to non-APL, intensively treated patients (see setup).
## Output: outputs/R_table2b_53feat.csv, outputs/R_multivariable_53feat.csv

source("00_setup.R")

ext <- load_external()
COH <- analysis_frames(ext)

## ---- reference Cox models fit on NTUH ----
nt <- load_ntuh()
NT <- ntuh_frame(nt)
NT$eln <- nt$ELN_2022 - 1            # NTUH raw 1=fav,2=int,3=adv -> 0..2, higher = worse
cox_age  <- coxph(Surv(time, status) ~ v1, data = NT)          # v1 = Age
cox_eln  <- coxph(Surv(time, status) ~ eln, data = NT)
cox_both <- coxph(Surv(time, status) ~ v1 + eln, data = NT)

rows <- list(); mv_rows <- list()
for (nm in names(COH)) {
  D <- COH[[nm]]
  r_rsf <- risk_of(D)
  Dx <- data.frame(v1 = D$v1, eln = eln_risk(D))
  c_eln3 <- cidx(D, eln_risk(D))
  stopifnot(abs(c_eln3 - cidx(D, predict(cox_eln, newdata = Dx, type = "lp"))) < 1e-9)
  rows[[nm]] <- data.frame(
    cohort = nm, n = nrow(D),
    ELN = c_eln3,
    Age = cidx(D, predict(cox_age, newdata = Dx, type = "lp")),
    ELN_plus_age = cidx(D, predict(cox_both, newdata = Dx, type = "lp")),
    RSF_3grp = cidx(D, 2 - rsf_grp(r_rsf)),
    RSF = cidx(D, r_rsf))

  ## multivariable Cox in the validation cohort: RSF per SD + ELN-2022 (adverse / intermediate
  ## vs favorable) + age per 10 years
  mv <- data.frame(time = D$time, status = D$status,
                   rsf_sd = (r_rsf - mean(r_rsf)) / sd(r_rsf),
                   eln = factor(GLAB[as.character(D$eln)], levels = c("favorable", "intermediate", "adverse")),
                   age10 = D$v1 / 10)
  s <- summary(coxph(Surv(time, status) ~ rsf_sd + eln + age10, data = mv))
  mv_rows[[nm]] <- data.frame(cohort = nm, term = rownames(s$coefficients),
                              HR = s$conf.int[, "exp(coef)"], lo = s$conf.int[, "lower .95"],
                              hi = s$conf.int[, "upper .95"], p = s$coefficients[, "Pr(>|z|)"],
                              row.names = NULL)
}
res <- bind_rows(rows); mv_res <- bind_rows(mv_rows)
write.csv(res, out_path("table2b"), row.names = FALSE)
write.csv(mv_res, out_path("multivariable"), row.names = FALSE)

cat("\n=== Table 2B (c-index) ===\n")
print(res, digits = 3, row.names = FALSE)
cat("\n=== Multivariable Cox: RSF per SD + ELN-2022 + age (per 10 y) ===\n")
for (nm in names(COH)) {
  r <- mv_res[mv_res$cohort == nm, ]
  cat(sprintf("\n[%s]\n", nm))
  for (j in seq_len(nrow(r)))
    cat(sprintf("  %-16s HR %.2f (%.2f-%.2f) p=%.2e\n", r$term[j], r$HR[j], r$lo[j], r$hi[j], r$p[j]))
}
cat("\nWrote", out_path("table2b"), "and", out_path("multivariable"), "\n")
