## Table 1: baseline characteristics, NTUH (training) and the three external cohorts.
## Variables and tests:
##   continuous -> Kruskal-Wallis; binary / ELN distribution -> chi-square; OS -> log-rank.
## External cohorts restricted to non-APL, intensively treated patients (see setup).
## Change from the Python version: "Median OS" is the Kaplan-Meier median (the Python
## script took the median of raw follow-up times, which counts censored patients as deaths);
## median follow-up is added (reverse Kaplan-Meier).
## Output: outputs/R_table1_53feat.csv

source("00_setup.R")

ext <- load_external()
nt  <- load_ntuh()

## one frame per cohort: Age, Gender, WBC, Hb, Plt, features, eln (0 adv,1 int,2 fav), time (y), status
cohort_df <- function(df, eln, time, status) {
  out <- df[, FEATS]
  out$eln <- as.integer(eln); out$time <- as.numeric(time); out$status <- as.integer(status)
  out
}
CO <- list(
  NTUH    = cohort_df(nt, 3 - nt$ELN_2022, nt$OSD_no_censor_2021 / 12, nt$OS_coding_no_censor_2021),
  NCRI    = cohort_df(ext$NCRI, ext$NCRI$ELN_2022, ext$NCRI$OSD_no_censor_2021, ext$NCRI$OS_coding_no_censor_2021),
  AMLSG   = cohort_df(ext$AMLSG, ext$AMLSG$ELN_2022_derived, ext$AMLSG$OSD_no_censor_2021, ext$AMLSG$OS_coding_no_censor_2021),
  BeatAML = cohort_df(ext$BeatAML, ext$BeatAML$ELN, ext$BeatAML$T, ext$BeatAML$E))
grp <- factor(rep(names(CO), sapply(CO, nrow)), levels = names(CO))
pool <- function(v) unlist(lapply(CO, function(d) as.numeric(d[[v]])), use.names = FALSE)

fmt_p   <- function(p) if (p < 0.001) "<0.001" else sprintf("%.3f", p)
med_iqr <- function(x) { q <- quantile(x, c(.5, .25, .75), na.rm = TRUE); sprintf("%.1f (%.1f-%.1f)", q[1], q[2], q[3]) }
n_pct   <- function(x) sprintf("%d (%.1f%%)", sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))

rows <- list()
add_row <- function(label, vals, p) rows[[length(rows) + 1]] <<- c(Characteristic = label, vals, p_value = p)

for (v in list(c("Age", "Age, median (IQR), y"), c("WBC", "WBC, median (IQR), /uL"),
               c("Hb", "Hb, median (IQR)"), c("Plt", "Platelets, median (IQR)"), c("BM_Blast", "BM blast %, median (IQR)"))) {
  x <- pool(v[1])
  add_row(v[2], sapply(CO, function(d) med_iqr(as.numeric(d[[v[1]]]))), fmt_p(kruskal.test(x ~ grp)$p.value))
}
bin_row <- function(label, var, level = 1) {
  x <- pool(var) == level
  add_row(label, sapply(CO, function(d) n_pct(as.numeric(d[[var]]) == level)),
          fmt_p(suppressWarnings(chisq.test(table(grp, x)))$p.value))
}
bin_row("Male sex, n (%)", "Gender")
for (g in c("NPM1", "FLT3_ITD", "FLT_TKD", "TP53", "DNMT3A", "RUNX1", "IDH2", "SRSF2", "ASXL1", "KIT"))
  bin_row(sprintf("%s mutated, n (%%)", g), g)
for (a in c("-5", "-7", "-17", "inv(16)", "t(8;21)")) bin_row(sprintf("%s, n (%%)", a), a)
bin_row("Complex karyotype, n (%)", "Cplx")
bin_row("CEBPA bZIP in-frame, n (%)", "CEBPA_bZIP")

## ELN-2022: one cohort x 3-level chi-square, same p on all three rows
p_eln <- fmt_p(chisq.test(table(grp, pool("eln")))$p.value)
for (e in 0:2)
  add_row(sprintf("ELN-2022 %s, n (%%)", GLAB[as.character(e)]),
          sapply(CO, function(d) n_pct(d$eln == e)), p_eln)

## outcome
st <- pool("status"); tm <- pool("time")
add_row("Deaths, n (%)", sapply(CO, function(d) n_pct(d$status == 1)),
        fmt_p(chisq.test(table(grp, st))$p.value))
km_med <- function(time, status) summary(survfit(Surv(time, status) ~ 1))$table[["median"]]
lr <- survdiff(Surv(tm, st) ~ grp)
add_row("Median OS (Kaplan-Meier), y", sapply(CO, function(d) sprintf("%.2f", km_med(d$time, d$status))),
        fmt_p(1 - pchisq(lr$chisq, length(lr$n) - 1)))
add_row("Median follow-up (reverse Kaplan-Meier), y",
        sapply(CO, function(d) sprintf("%.2f", km_med(d$time, 1 - d$status))), "")

tab <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(tab) <- c("Characteristic", sprintf("%s (n=%d)", names(CO), sapply(CO, nrow)), "p value (4-cohort)")
write.csv(tab, out_path("table1"), row.names = FALSE)
print(tab, row.names = FALSE, right = FALSE)
cat("\nWrote", out_path("table1"), "\n")
