## Manuscript Table 1 (NTUH, n=801): features and outcomes by ELN-2022 group, split within each
## ELN group by CCM-RSF group. using the 53-feature model groups
## (the Python version used the scikit-survival model, whose groups differ).
## P value compares the three CCM-RSF subgroups within each ELN group:
##   continuous -> Kruskal-Wallis; categorical -> chi-square; 2-yr OS / PFS -> log-rank.
## CR = CR Date recorded; Relapse = Relapse Date recorded, among CR achievers;
## 2-yr OS = Kaplan-Meier, all patients; 2-yr PFS = Kaplan-Meier of RFS, CR achievers only.
## OS / RFS come from the training data (OS from diagnosis); only CR / Relapse dates are taken
## from Cohort_20250527_with_RSF.xlsx.
## Rows follow the order of the Table 1 in "{Table} AML_ML" docx.
## Output: outputs/R_table1_ntuh_53feat.csv

source("00_setup.R")

nt <- load_ntuh()
src <- read_excel(NTUH_SRC, guess_max = 5000)
src[[ID_COL]] <- as.character(src[[ID_COL]])
src <- src[match(nt[[ID_COL]], src[[ID_COL]]), ]
stopifnot(identical(src[[ID_COL]], nt[[ID_COL]]))

nt$CR      <- !is.na(src[["CR Date"]])
nt$Relapse <- ifelse(nt$CR, !is.na(src[["Relapse Date"]]), NA)
nt$OS_y    <- nt$OSD_no_censor_2021 / 12
nt$RFS_y   <- nt$RFS_no_censor_2021 / 12
nt$ELN_g   <- factor(c("Fav", "Int", "Adv")[nt$ELN_2022], levels = c("Fav", "Int", "Adv"))
nt$RSF_g   <- factor(c("Adv", "Int", "Fav")[rsf_grp(rsf$predicted) + 1], levels = c("Fav", "Int", "Adv"))
G <- c("Fav", "Int", "Adv")
cat("ELN x CCM-RSF (NTUH):\n"); print(table(ELN = nt$ELN_g, RSF = nt$RSF_g))

fmt_p   <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "<0.001" else sprintf("%.3f", p)
med_iqr <- function(x, d = 1) { x <- x[!is.na(x)]; if (!length(x)) return("NA")
  q <- quantile(x, c(.5, .25, .75)); sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"), q[1], q[2], q[3]) }
n_pct   <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return("NA")
  sprintf("%d (%.1f%%)", sum(x), 100 * mean(x)) }
p_kw    <- function(v, g) { ok <- !is.na(v); if (length(unique(g[ok])) < 2) NA else kruskal.test(v[ok], g[ok])$p.value }
p_chi   <- function(v, g) { ok <- !is.na(v); tb <- table(g[ok], factor(v[ok], levels = c(FALSE, TRUE)))
  tb <- tb[rowSums(tb) > 0, , drop = FALSE]
  if (nrow(tb) < 2 || any(colSums(tb) == 0)) NA else suppressWarnings(chisq.test(tb)$p.value) }
km2     <- function(t, e) { if (length(t) < 2) return("NA")
  sprintf("%.1f%%", 100 * summary(survfit(Surv(t, e) ~ 1), times = 2, extend = TRUE)$surv) }
p_lr    <- function(t, e, g) { g <- droplevels(g); if (nlevels(g) < 2) return(NA)
  sd <- survdiff(Surv(t, e) ~ g); 1 - pchisq(sd$chisq, length(sd$n) - 1) }

rows <- list()
add <- function(label, cells) rows[[length(rows) + 1]] <<- c(Characteristic = label, cells)
by_eln <- function(fun_cell, fun_p, data = nt) {
  out <- c()
  for (e in G) {
    s <- data[data$ELN_g == e, ]
    out <- c(out, setNames(sapply(G, function(r) fun_cell(s[s$RSF_g == r, ])), paste(e, G, sep = "_")),
             setNames(fmt_p(fun_p(s)), paste0(e, "_P")))
  }
  out
}
cont <- function(label, col, d = 1) add(label, by_eln(function(s) med_iqr(s[[col]], d),
                                                    function(s) p_kw(s[[col]], s$RSF_g)))
cat_ <- function(label, col, data = nt) add(label, by_eln(function(s) n_pct(s[[col]] == 1),
                                                          function(s) p_chi(s[[col]] == 1, s$RSF_g), data))

## clinical
cont("Age, median (IQR), y", "Age"); cont("WBC, median (IQR), /uL", "WBC", 0)
cat_("Male, n (%)", "Gender"); cont("Hb, median (IQR)", "Hb"); cont("PLT, median (IQR)", "Plt")
cont("PB Blast %, median (IQR)", "PB_Blast")
cont("BM Blast %, median (IQR)", "BM_Blast")
## cytogenetics
for (k in c("inv(3)", "-7", "-5", "t(6;9)", "t(3;3)", "-17", "t(16;16)", "inv(16)", "t(9;22)",
            "t(8;21)", "del(5q)", "t(9;11)", "t(8;16)")) cat_(paste0(k, ", n (%)"), k)
cat_("Complex karyotype, n (%)", "Cplx")
## mutations, grouped as in the docx
mut_rows <- c(ASXL1 = "ASXL1", RUNX1 = "RUNX1", BCOR = "BCOR", SRSF2 = "SRSF2", EZH2 = "EZH2",
              SF3B1 = "SF3B1", U2AF1 = "U2AF1", ZRSR2 = "ZRSR2", STAG2 = "STAG2",
              "FLT3-ITD" = "FLT3_ITD", "FLT3-TKD" = "FLT_TKD", KIT = "KIT", PTPN11 = "PTPN11",
              NRAS = "NRAS", KRAS = "KRAS", TP53 = "TP53", WT1 = "WT1", PHF6 = "PHF6", NPM1 = "NPM1",
              IDH1 = "IDH1", IDH2 = "IDH2", TET2 = "TET2", DNMT3A = "DNMT3A", ETV6 = "ETV6",
              "CEBPA bZIP in-frame" = "CEBPA_bZIP", GATA2 = "GATA2", RAD21 = "Rad21", SMC1A = "SMC1A", SMC3 = "SMC3",
              STAG1 = "STAG1", "KMT2A-PTD" = "MLL", CBL = "c-CBL", CBLB = "CBL-b")
for (k in names(mut_rows)) cat_(paste0(k, ", n (%)"), mut_rows[[k]])
## outcomes
add("CR, n (%)", by_eln(function(s) n_pct(s$CR), function(s) p_chi(s$CR, s$RSF_g)))
cr <- nt[nt$CR, ]
add("Relapse, n (%)", by_eln(function(s) n_pct(s$Relapse), function(s) p_chi(s$Relapse, s$RSF_g), cr))
add("2-yr OS, %", by_eln(function(s) km2(s$OS_y, s$OS_coding_no_censor_2021),
                         function(s) p_lr(s$OS_y, s$OS_coding_no_censor_2021, s$RSF_g)))
add("2-yr PFS, %", by_eln(function(s) km2(s$RFS_y, s$RFS_coding_no_censor_2021),
                          function(s) p_lr(s$RFS_y, s$RFS_coding_no_censor_2021, s$RSF_g), cr))

tab <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
hdr_n <- table(nt$ELN_g)
write.csv(tab, out_path("table1_ntuh"), row.names = FALSE)
cat(sprintf("\nELN groups: Fav n=%d, Int n=%d, Adv n=%d\n", hdr_n["Fav"], hdr_n["Int"], hdr_n["Adv"]))
print(tab, row.names = FALSE, right = FALSE)
cat("\nWrote", out_path("table1_ntuh"), "\n")
