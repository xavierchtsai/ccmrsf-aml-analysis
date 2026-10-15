## Shared setup for the R scripts.
## Change log: in Beat AML the "MLL" feature counts KMT2A-PTD only
## (14/305; previously 20/305 including 2 KMT2A fusions, 2 point mutations and 2 untyped).
##   - CEBPA bZIP in-frame added as the 52nd feature, and Beat AML ELN-2022
## re-derived with it (beataml_harmonize.py reads CEBPA_Biallelic as "positive"/"yes" while the
## sheet codes "bi"/"mono", so no Beat AML patient was CEBPA-positive there).
##   - peripheral-blood blast % replaced by bone-marrow blast % (NCRI/AMLSG have
## no PB blasts; their "PB_Blast" column was BM blasts), and complex karyotype (Cplx) added.
##   - 53-feature RSF, trained on NTUH n=801: age, sex, WBC, Hb, platelets, BM blast %,
##     13 cytogenetic abnormalities, complex karyotype, 31 gene mutations, FLT3-ITD, CEBPA bZIP
##   - 2,000 patient-bootstrap resamples
##   - external cohorts restricted to the same population as NTUH:
##     non-APL AND intensively treated
##       UK-NCRI : drop APL (ICC_2022 contains "APL") and Intensive == 0
##       AMLSG   : drop APL (all AMLSG patients are Intensive == 1)
##       Beat AML: drop APL (specificDx contains "promyelocytic"); the harmonized
##                 file is already restricted to de novo + standard chemotherapy
##
## Every script does source("00_setup.R") and must be run with the working
## directory set to R/.

suppressPackageStartupMessages({
  library(randomForestSRC); library(survival); library(dplyr); library(readxl); library(ggplot2)
})

env_or <- function(var, default) { v <- Sys.getenv(var, unset = ""); if (nzchar(v)) v else default }

ROOT <- env_or("AML_ROOT", normalizePath(".."))                     # repository root
DATA <- env_or("AML_DATA", file.path(ROOT, "data"))                 # input files, see data/README.md
IN   <- file.path(ROOT, "outputs")                                   # derived inputs written by python/
OUT  <- env_or("AML_OUT", file.path(ROOT, "outputs"))
FIG  <- env_or("AML_FIG", file.path(ROOT, "figures"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

TAG        <- "53feat"
N_BOOT     <- 2000
MODEL_RDS  <- file.path(ROOT, "model", "rsf_model_NTUH_53feat.rds")
NTUH_SRC   <- file.path(DATA, "Cohort_20250527_with_RSF.xlsx")       # FLT3-ITD, ELN, HSCT_CR1 for NTUH
BEATAML_CSV <- file.path(IN, "beataml_harmonized_51feat.csv")        # python/beataml_harmonize.py output, n=319

out_path <- function(name) file.path(OUT, sprintf("R_%s_%s.csv", name, TAG))
fig_path <- function(name) file.path(FIG, sprintf("%s_%s.png", name, TAG))
save_fig <- function(obj, path, width, height, dpi = 600) {
  ggsave(path, obj, width = width, height = height, dpi = dpi)
  cat("Saved ->", path, "\n")
}

## ---- feature set (53) ----
KARY <- c("inv(3)","-7","-5","t(6;9)","t(3;3)","-17","t(16;16)","inv(16)",
          "t(9;22)","t(8;21)","del(5q)","t(9;11)","t(8;16)")
MUTS <- c("PTPN11","NRAS","KRAS","NPM1","FLT_TKD","MLL","KIT","RUNX1","WT1","ASXL1",
          "IDH1","IDH2","TET2","DNMT3A","TP53","SF3B1","U2AF1","SRSF2","ZRSR2","GATA2",
          "STAG1","STAG2","Rad21","SMC1A","SMC3","PHF6","CBL-b","c-CBL","ETV6","EZH2","BCOR")
FEATS <- c("Age","Gender","WBC","Hb","Plt","BM_Blast",KARY,"Cplx",MUTS,"FLT3_ITD","CEBPA_bZIP")
VMAP  <- setNames(paste0("v", seq_along(FEATS)), FEATS)
GLAB  <- c("0" = "adverse", "1" = "intermediate", "2" = "favorable")

## Display labels for figures (ggtext markdown): gene names in italics, no underscores,
## monosomies with a true minus sign
feat_label <- function(f) {
  fixed <- c(Age = "Age", Gender = "Sex", WBC = "WBC", Hb = "Hemoglobin", Plt = "Platelets",
             BM_Blast = "BM blast %", Cplx = "Complex karyotype",
             MLL = "*KMT2A*-PTD", FLT_TKD = "*FLT3*-TKD", FLT3_ITD = "*FLT3*-ITD",
             CEBPA_bZIP = "*CEBPA* bZIP", Rad21 = "*RAD21*", "CBL-b" = "*CBLB*", "c-CBL" = "*CBL*")
  vapply(f, function(x) {
    if (x %in% names(fixed)) fixed[[x]]
    else if (x %in% KARY) sub("^-", "\u2212", x)
    else paste0("*", x, "*")
  }, character(1), USE.NAMES = FALSE)
}
LV3   <- c("adverse", "intermediate", "favorable")

## Chinese chart-number column, written as a Unicode escape so the scripts parse the same
## way on Windows (non-UTF-8 locales) and macOS
ID_COL <- "\u75c5\u6b77\u865f"

## ---- NTUH training cohort (n=801) ----
## merge() sorts by chart number; the model is trained on exactly this row order, so
## rsf$predicted lines up with load_ntuh() rows in every script.
load_ntuh <- function() {
  nt <- read.csv(file.path(DATA, "NTUH_preprocessing_df_36.csv"), check.names = FALSE,
                 fileEncoding = "UTF-8")
  nt[[ID_COL]] <- as.character(nt[[ID_COL]])
  src <- read_excel(NTUH_SRC, guess_max = 5000)
  src[[ID_COL]] <- as.character(src[[ID_COL]])
  nt <- merge(nt, src[, c(ID_COL, "ITD NGS+PCR", "CEBPA bZIP", "BM_blasts", "Cplx")], by = ID_COL, all.x = TRUE)
  nt$FLT3_ITD <- as.numeric(nt[["ITD NGS+PCR"]]); nt$FLT3_ITD[is.na(nt$FLT3_ITD)] <- 0
  nt$CEBPA_bZIP <- suppressWarnings(as.numeric(nt[["CEBPA bZIP"]]))   # one "?" -> 0
  nt$CEBPA_bZIP[is.na(nt$CEBPA_bZIP)] <- 0
  nt$BM_Blast <- suppressWarnings(as.numeric(nt$BM_blasts))           # 1 missing -> cohort median
  nt$BM_Blast[is.na(nt$BM_Blast)] <- median(nt$BM_Blast, na.rm = TRUE)
  nt$Cplx <- suppressWarnings(as.numeric(nt$Cplx)); nt$Cplx[is.na(nt$Cplx)] <- 0   # 10 missing -> 0
  stopifnot(nrow(nt) == 801)
  nt
}
ntuh_frame <- function(nt) {
  X <- nt[, FEATS]; colnames(X) <- VMAP[FEATS]
  data.frame(time = nt$OSD_no_censor_2021 / 12, status = nt$OS_coding_no_censor_2021,
             X, check.names = FALSE)
}

## ---- external cohorts, with exclusion ----
## NCRI / AMLSG preprocessing files are row-aligned with the "total" sheet within each
## Source; FLT3-ITD, ICC_2022, Intensive and HSCT_CR1 are taken by that row order.
## Alignment is checked on Age.
raw_total <- read_excel(file.path(DATA, "AML raw data (NTUH+AMLSG+UKNCRI).xlsx"),
                        sheet = "total", guess_max = 10000)
ext_meta <- function(source_name, df) {
  m <- raw_total[raw_total$Source == source_name, ]
  stopifnot(nrow(m) == nrow(df), isTRUE(all.equal(as.numeric(m$Age), as.numeric(df$Age))))
  data.frame(FLT3_ITD = as.numeric(grepl("ITD", m$FLT3)),
             APL = grepl("APL", m$ICC_2022),
             non_intensive = !is.na(m$Intensive) & m$Intensive == 0,
             HSCT_CR1 = as.numeric(m$HSCT_CR1),
             CEBPA_bZIP = as.numeric(m$CEBPA_bZIP),
             BM_Blast = as.numeric(m$BM_blasts),
             Cplx = as.numeric(m$Cplx))
}
beataml_apl_subjects <- function() {
  f <- file.path(DATA, "beat_aml_lab.xlsx")
  ## some copies of the sheet have a title row above the header
  hdr <- names(suppressMessages(read_excel(f, sheet = "summary", n_max = 0)))
  s <- read_excel(f, sheet = "summary", skip = if ("dbgap_subject_id" %in% hdr) 0 else 1, guess_max = 2000)
  dx <- paste(s$specificDxAtInclusion, s$specificDxAtAcquisition)
  unique(s$dbgap_subject_id[grepl("promyelocytic", dx, ignore.case = TRUE)])
}

## The model's "MLL" feature is KMT2A-PTD (NTUH: identical to the KMT2A-PTD column of the master
## file; NCRI/AMLSG: KMT2A-PTD in the Tazi resource). beataml_harmonize.py sets MLL = 1 for any
## "MLL" or "KMT2A" token in variantSummary, which also catches KMT2A fusions and point mutations.
## Keep MLL = 1 in Beat AML only when the KMT2A/MLL token is labelled as a PTD.
beataml_kmt2a_ptd <- function(subject) {
  f <- file.path(DATA, "beat_aml_lab.xlsx")
  hdr <- names(suppressMessages(read_excel(f, sheet = "summary", n_max = 0)))
  s <- read_excel(f, sheet = "summary", skip = if ("dbgap_subject_id" %in% hdr) 0 else 1, guess_max = 2000)
  tok <- strsplit(as.character(s$variantSummary), "\\|")
  ptd <- vapply(tok, function(x) {
    x <- trimws(x); x <- x[grepl("^(MLL|KMT2A)\\b", x) & !grepl("^MLL[0-9]", x)]
    any(grepl("PTD|partial tandem", x, ignore.case = TRUE))
  }, logical(1))
  as.numeric(subject %in% unique(s$dbgap_subject_id[ptd]))
}

## CEBPA bZIP in-frame call from a Beat AML variantSummary token (ELN-2022: in-frame mutation in the
## bZIP region, amino acids 278-345: in-frame insertion/deletion/duplication/delins or missense).
cebpa_token_bzip <- function(tok) {
  x <- sub("^CEBPA\\s*", "", tok)
  parts <- unlist(strsplit(x, "\\band\\b"))
  any(vapply(parts, function(p) {
    if (grepl("fs|\\*", p, ignore.case = TRUE)) return(FALSE)               # frameshift / nonsense
    cpos <- regmatches(p, regexpr("c\\.(\\d+)", p))
    if (length(cpos)) {                                                       # cDNA notation
      ins <- regmatches(p, regexpr("(ins|del|dup)[ACGT]+", sub(";.*", "", p)))
      seq <- sub("^(ins|del|dup)", "", ins)
      if (!length(seq) || nchar(seq) %% 3 != 0) return(FALSE)                # not demonstrably in-frame
      aa <- ceiling(as.numeric(sub("c\\.", "", cpos)) / 3)
    } else {
      m <- regmatches(p, gregexpr("[A-Z][a-z]{0,2}(\\d{2,3})", p))[[1]]
      if (!length(m)) return(NA)                                              # no position given
      aa <- as.numeric(gsub("\\D", "", m))
      if (!grepl("ins|del|dup|[A-Z][a-z]{0,2}\\d+[A-Z]", p)) return(FALSE)
    }
    any(aa >= 278 & aa <= 345)
  }, logical(1)), na.rm = TRUE)
}

## Beat AML: re-select the specimen used by beataml_harmonize.py (same filter, first row per patient
## in sheet order), call CEBPA bZIP from that specimen and re-derive ELN-2022 with it.
## The re-derivation uses the same rules as beataml_harmonize.py; with CEBPA set to negative it
## must reproduce the original ELN exactly (checked below).
beataml_cebpa_eln <- function(ba) {
  f <- file.path(DATA, "beat_aml_lab.xlsx")
  hdr <- names(suppressMessages(read_excel(f, sheet = "summary", n_max = 0)))
  s <- read_excel(f, sheet = "summary", skip = if ("dbgap_subject_id" %in% hdr) 0 else 1, guess_max = 2000)
  vs <- tolower(trimws(as.character(s$vitalStatus)))
  keep <- grepl("ACUTE MYELOID", s$dxAtInclusion, ignore.case = TRUE) &
    s$typeInductionTx %in% "Standard Chemotherapy" & toupper(as.character(s$isDenovo)) %in% "TRUE" &
    !(s$ELN2017 %in% c("NonInitial", "NonAML")) & !is.na(s$overallSurvival) & s$overallSurvival > 0 &
    vs %in% c("dead", "alive")
  s <- s[keep, ]; s <- s[!duplicated(s$dbgap_subject_id), ]
  s <- s[match(ba$subject, s$dbgap_subject_id), ]
  stopifnot(identical(as.numeric(s$dbgap_subject_id), as.numeric(ba$subject)))
  tok <- strsplit(as.character(s$variantSummary), "\\|")
  bi  <- tolower(trimws(as.character(s$CEBPA_Biallelic))) %in% "bi"
  ## bZIP in-frame called from the variants listed in variantSummary; one patient (2036) is
  ## recorded only as biallelic without variant detail and is counted as bZIP-positive
  detail <- vapply(tok, function(x) any(grepl("^CEBPA\\s*\\(", trimws(x))), logical(1))
  bzip <- vapply(tok, function(x) { x <- trimws(x); x <- x[startsWith(x, "CEBPA")]
    length(x) > 0 && any(vapply(x, cebpa_token_bzip, logical(1))) }, logical(1)) | (bi & !detail)
  k <- gsub("\\[\\d+\\]", "", ifelse(is.na(s$karyotype), "", s$karyotype))
  complex <- lengths(regmatches(k, gregexpr("t\\(|inv\\(|del\\(|add\\(|der\\(|dup\\(|[+-][0-9]", k))) >= 3
  g <- function(v) ba[[v]] == 1
  mds <- rowSums(ba[, c("ASXL1","BCOR","EZH2","RUNX1","SF3B1","SRSF2","STAG2","U2AF1","ZRSR2")]) > 0
  eln_rule <- function(cebpa) {
    adv_c <- g("t(6;9)") | g("t(9;22)") | g("t(8;16)") | g("inv(3)") | g("t(3;3)") | g("-5") |
             g("del(5q)") | g("-7") | g("-17") | complex | g("TP53")
    ifelse(g("t(8;21)") | g("inv(16)") | g("t(16;16)"), 2L,
    ifelse(adv_c, 0L, ifelse(g("t(9;11)"), 1L,
    ifelse(g("NPM1"), ifelse(g("FLT3_ITD"), 1L, 2L),
    ifelse(cebpa, 2L, ifelse(g("FLT3_ITD"), 1L, ifelse(mds, 0L, 1L)))))))
  }
  stopifnot(all(eln_rule(rep(FALSE, nrow(ba))) == ba$ELN))   # reproduces beataml_harmonize.py exactly
  bm <- suppressWarnings(as.numeric(s[["%.Blasts.in.BM"]]))             # missing -> cohort median
  bm[is.na(bm)] <- median(bm, na.rm = TRUE)
  list(CEBPA_bZIP = as.numeric(bzip), ELN = eln_rule(bzip), BM_Blast = bm, Cplx = as.numeric(complex))
}

## Returns the filtered raw data frames (NCRI / AMLSG / BeatAML) and prints the
## exclusion flow. HSCT_CR1 is attached for NCRI/AMLSG (transplant analysis).
load_external <- function(exclude = TRUE, verbose = TRUE) {
  rd <- function(p) read.csv(p, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  nc <- rd(file.path(DATA, "NCRI_preprocessing_df.csv"))
  sg <- rd(file.path(IN, "SG_with_ELN2022.csv"))
  ba <- rd(BEATAML_CSV)
  mn <- ext_meta("UKNCRI", nc); ms <- ext_meta("AMLSG", sg)
  nc$FLT3_ITD <- mn$FLT3_ITD; nc$HSCT_CR1 <- mn$HSCT_CR1; nc$CEBPA_bZIP <- mn$CEBPA_bZIP
  sg$FLT3_ITD <- ms$FLT3_ITD; sg$HSCT_CR1 <- ms$HSCT_CR1; sg$CEBPA_bZIP <- ms$CEBPA_bZIP
  nc$BM_Blast <- mn$BM_Blast; nc$Cplx <- mn$Cplx; sg$BM_Blast <- ms$BM_Blast; sg$Cplx <- ms$Cplx
  ce <- beataml_cebpa_eln(ba); ba$CEBPA_bZIP <- ce$CEBPA_bZIP; ba$ELN <- ce$ELN
  ba$BM_Blast <- ce$BM_Blast; ba$Cplx <- ce$Cplx
  ba_apl <- ba$subject %in% beataml_apl_subjects()
  ba$MLL <- ba$MLL * beataml_kmt2a_ptd(ba$subject)       # KMT2A-PTD only (see above)

  flow <- data.frame(
    cohort        = c("NCRI", "AMLSG", "BeatAML"),
    n_total       = c(nrow(nc), nrow(sg), nrow(ba)),
    APL           = c(sum(mn$APL), sum(ms$APL), sum(ba_apl)),
    non_intensive = c(sum(mn$non_intensive & !mn$APL), sum(ms$non_intensive & !ms$APL), 0),
    n_analysed    = NA)
  if (exclude) {
    nc <- nc[!(mn$APL | mn$non_intensive), ]
    sg <- sg[!(ms$APL | ms$non_intensive), ]
    ba <- ba[!ba_apl, ]
  }
  flow$n_analysed <- c(nrow(nc), nrow(sg), nrow(ba))
  if (verbose) { cat("\n=== External cohort flow (exclude =", exclude, ") ===\n"); print(flow, row.names = FALSE) }
  list(NCRI = nc, AMLSG = sg, BeatAML = ba, flow = flow)
}

## time (years), status, eln (0=adv,1=int,2=fav), v1..v53
prep <- function(df, time, status, eln) {
  X <- df[, FEATS]; colnames(X) <- VMAP[FEATS]
  X <- data.frame(lapply(X, function(z) as.numeric(z)), check.names = FALSE)
  data.frame(time = as.numeric(time), status = as.integer(status), eln = as.integer(eln),
             X, check.names = FALSE)
}
analysis_frames <- function(ext) {
  list(NCRI    = prep(ext$NCRI,  ext$NCRI$OSD_no_censor_2021,  ext$NCRI$OS_coding_no_censor_2021,  ext$NCRI$ELN_2022),
       AMLSG   = prep(ext$AMLSG, ext$AMLSG$OSD_no_censor_2021, ext$AMLSG$OS_coding_no_censor_2021, ext$AMLSG$ELN_2022_derived),
       BeatAML = prep(ext$BeatAML, ext$BeatAML$T, ext$BeatAML$E, ext$BeatAML$ELN))
}
eln_risk <- function(D) 2 - D$eln   # higher = worse
eln_grp  <- function(D) D$eln       # 0 adv, 1 int, 2 fav

## ---- frozen 53-feature RSF and locked NTUH tertile cutpoints ----
load_model <- function() {
  if (!file.exists(MODEL_RDS)) stop("Model not found: run '01_build_model.R' first")
  rsf <- readRDS(MODEL_RDS)
  stopifnot(length(rsf$xvar.names) == length(FEATS))   # guards against older 50/51-feature models
  rsf
}
if (file.exists(MODEL_RDS)) {
  rsf <- load_model()
  cut <- quantile(rsf$predicted, c(1/3, 2/3))
  risk_of <- function(D) predict(rsf, newdata = D)$predicted          # mortality, higher = worse
  rsf_grp <- function(r) { g <- rep(1L, length(r)); g[r >= cut[2]] <- 0L; g[r < cut[1]] <- 2L; g }
}

## ---- c-index ----
cidx <- function(D, risk) 1 - survival::concordance(Surv(time, status) ~ risk, data = D)$concordance
