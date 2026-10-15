## Table 2A: c-index (RSF vs ELN-2022) with 2,000-resample paired bootstrap, pairwise HRs,
## categorical NRI at 2 years, and the CoxNet / GBSurv comparators.
## External cohorts restricted to non-APL, intensively treated patients (see setup).
## Output: outputs/R_table2a_53feat.csv, outputs/R_cohort_flow_53feat.csv

suppressPackageStartupMessages({ library(glmnet); library(gbm) })
source("00_setup.R")

ext <- load_external()
write.csv(ext$flow, out_path("cohort_flow"), row.names = FALSE)
COH <- analysis_frames(ext)

## ---- comparators fit on NTUH ----
NT <- ntuh_frame(load_ntuh())
## CoxNet: elastic-net Cox (alpha 0.5), predicted at the smallest lambda (sksurv default)
coxnet_fit <- glmnet(as.matrix(NT[, VMAP[FEATS]]), Surv(NT$time, NT$status),
                     family = "cox", alpha = 0.5, cox.ties = "breslow")
coxnet_risk <- function(D) as.numeric(predict(coxnet_fit, newx = as.matrix(D[, VMAP[FEATS]]),
                                              s = min(coxnet_fit$lambda), type = "link"))
## GBSurv: gradient-boosted Cox (300 trees, depth 3, shrinkage 0.05, subsample 0.8)
set.seed(0)
gb_fml <- as.formula(paste("Surv(time,status) ~", paste(VMAP[FEATS], collapse = " + ")))
gbsurv_fit <- gbm(gb_fml, data = NT, distribution = "coxph", n.trees = 300,
                  interaction.depth = 3, shrinkage = 0.05, bag.fraction = 0.8)
gbsurv_risk <- function(D) as.numeric(predict(gbsurv_fit, newdata = D, n.trees = 300, type = "link"))

## ---- paired bootstrap: both c-indices and their difference from the same resamples ----
## p = two-sided bootstrap p for delta (0 means < 1/N_BOOT)
boot_paired <- function(D, risk_a, risk_b, n = N_BOOT) {
  v <- t(replicate(n, { i <- sample(nrow(D), replace = TRUE); Di <- D[i, ]
                        c(cidx(Di, risk_a[i]), cidx(Di, risk_b[i])) }))
  d <- v[, 1] - v[, 2]
  list(a = quantile(v[, 1], c(.025, .975)), b = quantile(v[, 2], c(.025, .975)),
       d = quantile(d, c(.025, .975)), p = min(1, 2 * min(mean(d <= 0), mean(d >= 0))))
}

## ---- pairwise HR (Cox) with log-rank p ----
pairwise_hr <- function(T, E, grp, hi, lo) {
  m <- grp %in% c(hi, lo)
  if (length(unique(grp[m])) < 2) return(c(HR = NA, p = NA))
  x <- as.integer(grp[m] == hi)
  cox_fit <- tryCatch(coxph(Surv(T[m], E[m]) ~ x), error = function(e) NULL)
  hr <- if (is.null(cox_fit)) NA else exp(coef(cox_fit)[[1]])
  sd <- tryCatch(survdiff(Surv(T[m], E[m]) ~ x), error = function(e) NULL)
  p <- if (is.null(sd)) NA else 1 - pchisq(sd$chisq, length(sd$n) - 1)
  c(HR = hr, p = p)
}
hr_triplet <- function(T, E, grp) {
  list(adv_int = pairwise_hr(T, E, grp, 0, 1),
       int_fav = pairwise_hr(T, E, grp, 1, 2),
       adv_fav = pairwise_hr(T, E, grp, 0, 2))
}

## ---- categorical NRI at 2 years (drops patients censored before 2 years) ----
nri_categorical <- function(T, E, old, new, h = 2) {
  event <- (T <= h) & (E == 1)
  usable <- ((E == 1) & (T <= h)) | (T >= h)
  ev <- usable & event; nev <- usable & !event
  up <- new < old; dn <- new > old             # lower group number = higher risk
  nri_ev  <- if (sum(ev))  mean(up[ev])  - mean(dn[ev])  else 0
  nri_nev <- if (sum(nev)) mean(dn[nev]) - mean(up[nev]) else 0
  c(events = nri_ev, nonevents = nri_nev, total = nri_ev + nri_nev)
}

## ---- assemble ----
rows <- list()
for (k in seq_along(COH)) {
  nm <- names(COH)[k]; D <- COH[[nm]]
  T <- D$time; E <- D$status
  r_rsf <- risk_of(D); r_eln <- eln_risk(D)
  g_rsf <- rsf_grp(r_rsf); g_eln <- eln_grp(D)

  set.seed(7 + k)
  b <- boot_paired(D, r_rsf, r_eln)
  ci_rsf <- cidx(D, r_rsf); ci_eln <- cidx(D, r_eln)
  hr_eln <- hr_triplet(T, E, g_eln); hr_rsf <- hr_triplet(T, E, g_rsf)
  nri <- nri_categorical(T, E, g_eln, g_rsf)

  rows[[nm]] <- data.frame(
    cohort = nm, n = nrow(D), deaths = sum(E),
    c_RSF = ci_rsf, RSF_lo = b$a[1], RSF_hi = b$a[2],
    c_ELN = ci_eln, ELN_lo = b$b[1], ELN_hi = b$b[2],
    delta = ci_rsf - ci_eln, delta_lo = b$d[1], delta_hi = b$d[2], delta_p = b$p,
    ELN_HR_adv_int = hr_eln$adv_int["HR"], ELN_p_adv_int = hr_eln$adv_int["p"],
    ELN_HR_int_fav = hr_eln$int_fav["HR"], ELN_p_int_fav = hr_eln$int_fav["p"],
    ELN_HR_adv_fav = hr_eln$adv_fav["HR"], ELN_p_adv_fav = hr_eln$adv_fav["p"],
    RSF_HR_adv_int = hr_rsf$adv_int["HR"], RSF_p_adv_int = hr_rsf$adv_int["p"],
    RSF_HR_int_fav = hr_rsf$int_fav["HR"], RSF_p_int_fav = hr_rsf$int_fav["p"],
    RSF_HR_adv_fav = hr_rsf$adv_fav["HR"], RSF_p_adv_fav = hr_rsf$adv_fav["p"],
    NRI_events = nri["events"], NRI_nonevents = nri["nonevents"], NRI_total = nri["total"],
    c_CoxNet = cidx(D, coxnet_risk(D)), c_GBSurv = cidx(D, gbsurv_risk(D)),
    row.names = NULL)
}
res <- bind_rows(rows)
write.csv(res, out_path("table2a"), row.names = FALSE)

cat(sprintf("\n=== Table 2A (53-feature, %d bootstrap resamples) ===\n", N_BOOT))
for (nm in names(COH)) {
  r <- res[res$cohort == nm, ]
  cat(sprintf("\n[%s] n=%d, deaths=%d\n", nm, r$n, r$deaths))
  cat(sprintf("  ELN-2022  c=%.3f (%.3f-%.3f)  HR adv/int=%.2f int/fav=%.2f adv/fav=%.2f\n",
              r$c_ELN, r$ELN_lo, r$ELN_hi, r$ELN_HR_adv_int, r$ELN_HR_int_fav, r$ELN_HR_adv_fav))
  cat(sprintf("  RSF       c=%.3f (%.3f-%.3f)  HR adv/int=%.2f int/fav=%.2f adv/fav=%.2f\n",
              r$c_RSF, r$RSF_lo, r$RSF_hi, r$RSF_HR_adv_int, r$RSF_HR_int_fav, r$RSF_HR_adv_fav))
  cat(sprintf("  delta     %+.3f (%+.3f to %+.3f), bootstrap p=%s\n", r$delta, r$delta_lo, r$delta_hi,
              if (r$delta_p == 0) sprintf("<%.4f", 1 / N_BOOT) else sprintf("%.4f", r$delta_p)))
  cat(sprintf("  NRI(ELN->RSF, 2y) total=%+.3f (events %+.3f, nonevents %+.3f)\n",
              r$NRI_total, r$NRI_events, r$NRI_nonevents))
  cat(sprintf("  CoxNet c=%.3f | GBSurv c=%.3f\n", r$c_CoxNet, r$c_GBSurv))
}
cat("\nWrote", out_path("table2a"), "\n")
