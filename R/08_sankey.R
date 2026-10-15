## ELN-2022 -> RSF reclassification alluvial plot, with cross-tabs and % reclassified.
## External cohorts restricted to non-APL, intensively treated patients (see setup).
## Output: figures/Fig_Sankey_ELN_to_RSF_53feat.png

suppressPackageStartupMessages({ library(ggalluvial); library(gridExtra) })
source("00_setup.R")

ext <- load_external()
COH <- analysis_frames(ext)

make_alluvial_plot <- function(nm, D) {
  plot_data <- data.frame(ELN = factor(GLAB[as.character(eln_grp(D))], levels = LV3),
                          RSF = factor(GLAB[as.character(rsf_grp(risk_of(D)))], levels = LV3)) |>
    count(ELN, RSF, name = "freq")
  ggplot(plot_data, aes(axis1 = ELN, axis2 = RSF, y = freq)) +
    geom_alluvium(aes(fill = ELN), alpha = .7, width = 1/6) +
    geom_stratum(width = 1/6, fill = "grey90", color = "grey40") +
    geom_text(stat = "stratum", size = 2.6,
              aes(label = paste0(after_stat(stratum), "\n(n=", after_stat(count), ")"))) +
    scale_x_discrete(limits = c("ELN-2022", "CCM-RSF"), expand = c(.05, .05)) +
    scale_fill_manual(values = c(adverse = "#d62728", intermediate = "#ff7f0e", favorable = "#2ca02c"), guide = "none") +
    labs(title = sprintf("%s (n=%d)", nm, nrow(D)), y = "Patients", x = NULL) +
    theme_minimal() + theme(panel.grid = element_blank())
}
plots <- lapply(names(COH), function(nm) make_alluvial_plot(nm, COH[[nm]]))
save_fig(arrangeGrob(grobs = plots, ncol = length(plots),
                     top = "Reclassification of ELN-2022 categories by CCM-RSF"),
         fig_path("Fig_Sankey_ELN_to_RSF"), width = 4.5 * length(plots), height = 5)

for (nm in names(COH)) {
  D <- COH[[nm]]
  e <- factor(GLAB[as.character(eln_grp(D))], levels = LV3)
  r <- factor(GLAB[as.character(rsf_grp(risk_of(D)))], levels = LV3)
  cat("\n===", nm, "=== (rows ELN-2022, cols RSF)\n"); print(table(ELN = e, RSF = r))
  k <- sum(as.character(e) != as.character(r))
  cat(sprintf("Reclassified (RSF != ELN-2022): %d / %d (%.1f%%)\n", k, length(e), 100 * k / length(e)))
}
