## Decision curve analysis (death by 2 years; "positive" = classified adverse by RSF or ELN-2022).
## Net benefit = TP/n - FP/n * pt/(1-pt); patients censored before 2 years without an event are dropped.
## External cohorts restricted to non-APL, intensively treated patients (see setup).
## Output: figures/Fig_DCA_<cohort>_53feat.png, outputs/R_dca_53feat.csv

source("00_setup.R")

ext <- load_external()
COH <- analysis_frames(ext)

decision_curve <- function(time, status, eln_adv, rsf_adv, h = 2) {
  event  <- as.integer(time <= h & status == 1)
  usable <- (status == 1 & time <= h) | (time >= h)
  n <- sum(usable); y <- event[usable]; prev <- mean(y)
  eln_pred <- eln_adv[usable]; rsf_pred <- rsf_adv[usable]
  nb <- function(pred, w) sum(pred & y == 1) / n - sum(pred & y == 0) / n * w
  do.call(rbind, lapply(seq(0.05, 0.75, by = 0.05), function(pt) {
    w <- pt / (1 - pt)
    data.frame(threshold = pt, RSF_adverse = nb(rsf_pred, w), ELN_adverse = nb(eln_pred, w),
               all = prev - (1 - prev) * w, none = 0)
  }))
}
dc <- lapply(COH, function(D) decision_curve(D$time, D$status, D$eln == 0, rsf_grp(risk_of(D)) == 0))

plot_dca <- function(d, ttl) {
  long <- bind_rows(
    data.frame(threshold = d$threshold, series = "CCM-RSF adverse", value = d$RSF_adverse),
    data.frame(threshold = d$threshold, series = "ELN adverse",     value = d$ELN_adverse),
    data.frame(threshold = d$threshold, series = "Treat all",       value = d$all),
    data.frame(threshold = d$threshold, series = "Treat none",      value = d$none))
  lv <- c("CCM-RSF adverse", "ELN adverse", "Treat all", "Treat none")
  long$series <- factor(long$series, levels = lv)
  ggplot(long, aes(threshold, value, color = series, linetype = series, shape = series)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
    scale_color_manual(values = setNames(c("#1f77b4", "#ff7f0e", "grey50", "black"), lv)) +
    scale_linetype_manual(values = setNames(c("solid", "solid", "dotted", "solid"), lv)) +
    scale_shape_manual(values = setNames(c(16, 15, NA, NA), lv)) +
    scale_x_continuous(breaks = seq(0, 0.8, 0.1)) +
    labs(x = "Threshold probability", y = "Net benefit",
         title = paste0(ttl, " decision curve (death by 2 yr)"), color = NULL, linetype = NULL, shape = NULL) +
    theme_minimal() + theme(legend.position = "right")
}
for (nm in names(dc))
  save_fig(plot_dca(dc[[nm]], sprintf("%s (n=%d)", nm, nrow(COH[[nm]]))),
           fig_path(paste0("Fig_DCA_", nm)), width = 6.6, height = 4.4)

dc_out <- bind_rows(lapply(names(dc), function(nm) cbind(cohort = nm, dc[[nm]])))
cat("\nNet benefit at the 50% threshold:\n")
print(dc_out[dc_out$threshold == 0.5, c("cohort", "RSF_adverse", "ELN_adverse")], digits = 3, row.names = FALSE)
write.csv(dc_out, out_path("dca"), row.names = FALSE)
cat("\nWrote", out_path("dca"), "\n")
