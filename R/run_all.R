## Run the whole R pipeline in order. Working directory must be R/.
##   macOS:   Rscript run_all.R
##   Windows: use "run_all_windows.R" (sets the working directory and checks packages)
for (f in c("01_build_model", "02_table_s1_baseline", "03_table1_ntuh", "04_table2a_discrimination", "05_table2b_decomposition", "06_figures_table3",
            "07_decision_curves", "08_sankey", "09_transplant_table4", "10_partial_dependence")) {
  cat("\n\n########################", f, "########################\n")
  source(sprintf("%s.R", f), local = new.env(), echo = FALSE, encoding = "UTF-8")
}
