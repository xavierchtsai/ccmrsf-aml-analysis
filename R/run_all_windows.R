## Windows launcher for the R pipeline.
##
## How to run (either way):
##   RStudio : open this file and click "Source" (Ctrl+Shift+S)
##   R GUI   : File > Source R code... and pick this file
##   Command : "C:\Program Files\R\R-4.x.x\bin\Rscript.exe" "run_all_windows.R"
##             (run from inside the analysis\R folder)
##
## If the working directory cannot be detected automatically, fill in R_DIR below,
## using forward slashes, e.g. "C:/projects/ccmrsf-aml-analysis/R".
R_DIR <- ""

## ---- R version: 4.2+ runs in UTF-8 on Windows, needed for the Chinese folder name ----
if (getRversion() < "4.2.0")
  stop("R ", getRversion(), " found. Please install R 4.2 or later (UTF-8 support on Windows).")

## ---- locate the R/ folder ----
find_script_dir <- function() {
  if (nzchar(R_DIR)) return(R_DIR)
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))      # Rscript
  if (length(f)) return(dirname(normalizePath(f)))
  for (i in rev(seq_len(sys.nframe()))) {                                          # source()
    o <- sys.frame(i)$ofile
    if (!is.null(o)) return(dirname(normalizePath(o)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) # RStudio
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  stop("Cannot find the R/ folder: set R_DIR at the top of this file.")
}
setwd(find_script_dir())
if (!file.exists("00_setup.R")) stop("00_setup.R not found in ", getwd())
cat("Working directory:", getwd(), "\n")

## ---- packages ----
pkgs <- c("randomForestSRC", "survival", "dplyr", "readxl", "ggplot2", "glmnet", "gbm",
          "survminer", "gridExtra", "cowplot", "ggalluvial", "scales")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  cat("Installing:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, repos = "https://cloud.r-project.org")
}

## ---- run ----
source("run_all.R", encoding = "UTF-8")
cat("\nDone. Tables in outputs/, figures in figures/.\n")
