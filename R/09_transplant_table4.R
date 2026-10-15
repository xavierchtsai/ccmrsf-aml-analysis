## Transplant analyses and Table 4.
##
## Primary (manuscript Results 3.5, Figure 5, Supplementary Figure 3):
##   patients reclassified from ELN-2022 favorable/intermediate to CCM-RSF adverse,
##   UK-NCRI + AMLSG (the cohorts with HSCT-in-CR1 data), 3-month landmark from diagnosis,
##   Cox model for HSCT in CR1 adjusted for age and stratified by cohort; 2-yr OS after the landmark.
## Sensitivity (text output only): same analysis with NTUH added, 4-month landmark,
##   all CCM-RSF-adverse patients, and a pooled multivariable landmark model.
## Table 4: Cox model for OS with age, CCM-RSF group and HSCT in CR1, separately in NTUH,
##   UK-NCRI and AMLSG (no landmark).
## AMLSG / UK-NCRI restricted to non-APL, intensively treated patients (see setup);
## UK-NCRI then keeps only patients with HSCT_CR1 recorded. Beat AML has no transplant timing.
## Output: outputs/R_transplant_53feat.csv / .txt, outputs/R_table4_53feat.csv
##         figures/Fig5_HSCT_reclassified_*.png, figures/FigS3_HSCT_reclassified_by_cohort_*.png

suppressPackageStartupMessages({ library(survminer); library(cowplot) })
source("00_setup.R")

## ELN grouping by median OS, so each cohort's own ELN coding is handled the same way
eln_group_by_median_os <- function(eln_code, os_time) {
  med <- tapply(os_time, eln_code, median, na.rm = TRUE)
  med <- med[!is.na(names(med))]
  out <- rep("intermediate", length(eln_code))
  out[as.character(eln_code) == names(med)[which.min(med)]] <- "adverse"
  out[as.character(eln_code) == names(med)[which.max(med)]] <- "favorable"
  out[is.na(eln_code)] <- NA
  out
}

## ---- NTUH (months) ----
## OS, age and ELN come from the training data (OS measured from diagnosis). Only HSCT_CR1 is
## taken from Cohort_20250527_with_RSF.xlsx: in that file OSD_no_censor_2021 of the 158 patients
## transplanted in CR1 is measured from the transplant date, so it must not be used here.
nt <- load_ntuh()
src <- read_excel(NTUH_SRC, guess_max = 5000)
src[[ID_COL]] <- as.character(src[[ID_COL]])
src <- src[match(nt[[ID_COL]], src[[ID_COL]]), ]
stopifnot(identical(src[[ID_COL]], nt[[ID_COL]]))
ntuh_tx <- data.frame(cohort = "NTUH", grp = rsf_grp(rsf$predicted), risk = rsf$predicted,
                      eln3 = GLAB[as.character(3 - nt$ELN_2022)],      # raw 1=fav,2=int,3=adv
                      cr1 = as.numeric(as.numeric(src$HSCT_CR1) == 1), T = nt$OSD_no_censor_2021,
                      event = nt$OS_coding_no_censor_2021, age = nt$Age)

## ---- AMLSG / UK-NCRI (preprocessing files store OS in years) ----
ext <- load_external()
ext_frame <- function(df, source_name, eln_column) {
  X <- df[, FEATS]; colnames(X) <- VMAP[FEATS]
  X <- data.frame(lapply(X, as.numeric), check.names = FALSE)
  r <- predict(rsf, newdata = X)$predicted
  os_m <- df$OSD_no_censor_2021 * 12
  data.frame(cohort = source_name, grp = rsf_grp(r), risk = r,
             eln3 = eln_group_by_median_os(df[[eln_column]], os_m),
             cr1 = df$HSCT_CR1, T = os_m, event = df$OS_coding_no_censor_2021, age = df$Age)
}
pooled <- bind_rows(ntuh_tx,
                    ext_frame(ext$AMLSG, "AMLSG", "ELN_2022_derived"),
                    ext_frame(ext$NCRI, "NCRI", "ELN_2022"))
pooled <- pooled[!is.na(pooled$T) & pooled$T > 0 & !is.na(pooled$event) & !is.na(pooled$cr1), ]
cat("\nPatients with HSCT_CR1 recorded and OS > 0:\n"); print(table(pooled$cohort))
pooled$rsf <- GLAB[as.character(pooled$grp)]
pooled$eln_adv <- as.integer(pooled$eln3 == "adverse")
pooled$reclass <- as.integer(pooled$rsf == "adverse" & pooled$eln3 %in% c("favorable", "intermediate"))

landmark <- function(df, lm) { ld <- df[df$T >= lm, ]; ld$Tt <- ld$T - lm; ld }
landmark_hr <- function(df, lm) {
  ld <- landmark(df, lm)
  if (sum(ld$cr1 == 1) < 5 || sum(ld$cr1 == 0) < 5) return(NULL)
  rhs <- c("cr1", "age")
  if (length(unique(ld$eln_adv)) > 1) rhs <- c(rhs, "eln_adv")
  fml <- paste("Surv(Tt, event) ~", paste(rhs, collapse = " + "),
               if (length(unique(ld$cohort)) > 1) "+ strata(cohort)" else "")
  s <- summary(coxph(as.formula(fml), data = ld))
  os <- sapply(c(1, 0), function(v)
    summary(survfit(Surv(Tt, event) ~ 1, data = ld[ld$cr1 == v, ]), times = 24, extend = TRUE)$surv)
  list(hr = s$coefficients["cr1", "exp(coef)"], lo = s$conf.int["cr1", "lower .95"],
       hi = s$conf.int["cr1", "upper .95"], p = s$coefficients["cr1", "Pr(>|z|)"],
       n_exp = sum(ld$cr1 == 1), n_noexp = sum(ld$cr1 == 0), os_exp = os[1], os_noexp = os[2])
}

lines <- c("=== Transplant (HSCT in CR1), CCM-RSF 53-feature ===")
csv_rows <- list()
report <- function(label, df, lm) {
  r <- landmark_hr(df, lm)
  if (is.null(r)) return(invisible())
  lines <<- c(lines, sprintf("%-62s HR %.2f (%.2f-%.2f) p=%.2g | HSCT %d vs no HSCT %d | 2-yr OS %.1f%% vs %.1f%%",
                             label, r$hr, r$lo, r$hi, r$p, r$n_exp, r$n_noexp, 100 * r$os_exp, 100 * r$os_noexp))
  csv_rows[[length(csv_rows) + 1]] <<- data.frame(label = label, hr = r$hr, lo = r$lo, hi = r$hi, p = r$p,
                                                  n_exp = r$n_exp, n_noexp = r$n_noexp,
                                                  os2_exp = r$os_exp, os2_noexp = r$os_noexp)
}

ext2 <- pooled[pooled$cohort %in% c("NCRI", "AMLSG"), ]
q2 <- ext2[ext2$reclass == 1, ]
lines <- c(lines, "\n[PRIMARY] ELN fav/int reclassified to CCM-RSF adverse, NCRI + AMLSG")
report("Reclassified, NCRI+AMLSG, 3-month landmark", q2, 3)
report("Reclassified, NCRI only, 3-month landmark", q2[q2$cohort == "NCRI", ], 3)
report("Reclassified, AMLSG only, 3-month landmark", q2[q2$cohort == "AMLSG", ], 3)

lines <- c(lines, "\n[SENSITIVITY]")
report("Reclassified, NCRI+AMLSG, 4-month landmark", q2, 4)
report("Reclassified, NTUH+NCRI+AMLSG, 3-month landmark", pooled[pooled$reclass == 1, ], 3)
report("All CCM-RSF adverse, NCRI+AMLSG, 3-month landmark", ext2[ext2$rsf == "adverse", ], 3)
report("All CCM-RSF adverse, NTUH+NCRI+AMLSG, 3-month landmark", pooled[pooled$rsf == "adverse", ], 3)

## pooled multivariable landmark model (NTUH + NCRI + AMLSG)
mv <- landmark(pooled, 3)
mv$rsf_z <- (mv$risk - mean(mv$risk)) / sd(mv$risk)
s <- summary(coxph(Surv(Tt, event) ~ cr1 + rsf_z + age + eln_adv + strata(cohort), data = mv))
lines <- c(lines, sprintf("\n[Multivariable, 3-month landmark, NTUH+NCRI+AMLSG, n=%d]", nrow(mv)))
for (v in c("cr1", "rsf_z", "age", "eln_adv"))
  lines <- c(lines, sprintf("  %-8s HR %.3f (%.3f-%.3f) p=%.2e", v, s$conf.int[v, 1],
                            s$conf.int[v, 3], s$conf.int[v, 4], s$coefficients[v, 5]))

txt <- paste(lines, collapse = "\n"); cat(txt, "\n")
writeLines(txt, sub("\\.csv$", ".txt", out_path("transplant")))
write.csv(bind_rows(csv_rows), out_path("transplant"), row.names = FALSE)

## ============================ Figure 5 / Supplementary Figure 3 ============================
km_tx <- function(ld, strata_var, labels, palette, linetypes, ttl, p_text) {
  ld$s <- factor(strata_var, levels = names(labels))
  fit <- survfit(Surv(Tt, event) ~ s, data = ld)
  gg <- ggsurvplot(fit, data = ld, palette = palette, linetype = linetypes, censor = FALSE,
                   legend.title = "", legend.labs = unname(labels), legend = c(0.75, 0.88),
                   xlab = "Time since 3-month landmark (months)", ylab = "Overall survival",
                   break.time.by = 40, title = ttl, risk.table = TRUE, risk.table.y.text = FALSE,
                   fontsize = 3.5, size = 0.9)
  gg$plot <- gg$plot + annotate("text", x = 0, y = 0.04 + 0.07 * (seq_along(p_text) - 1),
                                label = rev(p_text), hjust = 0, size = 4) +
    theme(legend.background = element_blank(), plot.title = element_text(size = 11))
  gg$table <- gg$table + labs(title = NULL, x = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  cowplot::plot_grid(gg$plot, gg$table, ncol = 1, align = "v", axis = "lr", rel_heights = c(3, 1))
}
logrank_p <- function(ld) {
  sd <- survdiff(Surv(Tt, event) ~ cr1, data = ld)
  p <- 1 - pchisq(sd$chisq, 1)
  if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
}

ld <- landmark(q2, 3)
f5 <- km_tx(ld, ifelse(ld$cr1 == 1, "tx", "no"), c(tx = "HSCT in CR1", no = "No HSCT"),
            c("#1f77b4", "#e31a1c"), c("solid", "solid"),
            "ELN-2022 favorable/intermediate reclassified to CCM-RSF adverse\n(3-month landmark survivors, NCRI + AMLSG)",
            paste("Log-rank", logrank_p(ld)))
save_fig(f5, fig_path("Fig5_HSCT_reclassified"), width = 6.5, height = 6)

key <- paste(ld$cohort, ifelse(ld$cr1 == 1, "tx", "no"))
s3 <- km_tx(ld, key,
            c("NCRI tx" = "NCRI - HSCT in CR1", "NCRI no" = "NCRI - No HSCT",
              "AMLSG tx" = "AMLSG - HSCT in CR1", "AMLSG no" = "AMLSG - No HSCT"),
            c("#1f77b4", "#a6cee3", "#e31a1c", "#fb9a99"), c("solid", "dashed", "solid", "dashed"),
            "ELN-2022 favorable/intermediate reclassified to CCM-RSF adverse\n(3-month landmark survivors)",
            c(paste("NCRI log-rank", logrank_p(ld[ld$cohort == "NCRI", ])),
              paste("AMLSG log-rank", logrank_p(ld[ld$cohort == "AMLSG", ]))))
save_fig(s3, fig_path("FigS3_HSCT_reclassified_by_cohort"), width = 6.5, height = 7)

## ============================ Table 4 ============================
## OS ~ age + CCM-RSF group (favorable reference) + HSCT in CR1, per cohort, no landmark
t4 <- bind_rows(lapply(c("NTUH", "NCRI", "AMLSG"), function(cc) {
  d <- pooled[pooled$cohort == cc, ]
  d$rsf_f <- factor(d$rsf, levels = c("favorable", "intermediate", "adverse"))
  s <- summary(coxph(Surv(T, event) ~ age + rsf_f + cr1, data = d))
  data.frame(cohort = cc, n = nrow(d),
             term = c("Age (year)", "CCM-RSF intermediate vs favorable", "CCM-RSF adverse vs favorable", "HSCT at first CR"),
             HR = s$conf.int[, 1], lo = s$conf.int[, 3], hi = s$conf.int[, 4], p = s$coefficients[, 5],
             row.names = NULL)
}))
cat("\n=== Table 4: Cox OS ~ age + CCM-RSF group + HSCT in CR1 ===\n")
for (j in seq_len(nrow(t4)))
  cat(sprintf("%-6s n=%-5d %-34s HR %.2f (%.2f-%.2f) p=%s\n", t4$cohort[j], t4$n[j], t4$term[j],
              t4$HR[j], t4$lo[j], t4$hi[j], if (t4$p[j] < 0.001) "<0.001" else sprintf("%.3f", t4$p[j])))
write.csv(t4, out_path("table4"), row.names = FALSE)
cat("\nWrote", out_path("transplant"), ",", out_path("table4"), "and Figure 5 / S3\n")
