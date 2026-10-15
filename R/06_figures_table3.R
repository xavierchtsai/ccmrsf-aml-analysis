## Figures 1-4, ELN x RSF cross-tables, and Table 4 (within-ELN restratification).
## External cohorts restricted to non-APL, intensively treated patients (see setup).
## Figure 2 reads the c-index CIs from 04_table2a_discrimination.R so the figure and Table 2A show the same
## 2,000-resample bootstrap; run "04_table2a_discrimination.R" first.
## Output: figures/Fig1_KM_*, Fig2_cindex_forest_*, Fig3_calibration_*, Fig4_importance_*
##         outputs/R_within_eln_53feat.csv, outputs/R_crosstab_53feat.csv

suppressPackageStartupMessages({ library(survminer); library(gridExtra); library(cowplot) })
source("00_setup.R")
set.seed(0)

ext <- load_external()
COH <- analysis_frames(ext)
nt  <- load_ntuh()

## ---- ELN x RSF cross-tables (NTUH training + external) ----
ntuh_eln <- factor(GLAB[as.character(3 - nt$ELN_2022)], levels = LV3)   # NTUH raw 1=fav,2=int,3=adv
ntuh_rsf <- factor(GLAB[as.character(rsf_grp(rsf$predicted))], levels = LV3)
xt_list <- list(NTUH = table(ELN = ntuh_eln, RSF = ntuh_rsf))
for (nm in names(COH)) {
  D <- COH[[nm]]
  xt_list[[nm]] <- table(ELN = factor(GLAB[as.character(eln_grp(D))], levels = LV3),
                         RSF = factor(GLAB[as.character(rsf_grp(risk_of(D)))], levels = LV3))
}
for (nm in names(xt_list)) {
  cat("\n===", nm, "=== (rows ELN-2022, cols RSF)\n"); print(addmargins(xt_list[[nm]]))
  cat(sprintf("Reclassified: %.1f%%\n", 100 * (1 - sum(diag(xt_list[[nm]])) / sum(xt_list[[nm]]))))
}
write.csv(bind_rows(lapply(names(xt_list), function(nm)
  data.frame(cohort = nm, as.data.frame(xt_list[[nm]])))), out_path("crosstab"), row.names = FALSE)

## ============================ FIGURE 1 : KM by group =========================
km_panel <- function(D, grp, ttl) {
  d <- data.frame(time = D$time, status = D$status,
                  grp = factor(GLAB[as.character(grp)], levels = LV3))
  fit <- survfit(Surv(time, status) ~ grp, data = d)
  strata_labs <- sub("^grp=", "", names(fit$strata))
  strata_labs <- paste0(toupper(substr(strata_labs, 1, 1)), substr(strata_labs, 2, nchar(strata_labs)))
  gg <- ggsurvplot(fit, data = d, palette = c("#d62728", "#ff7f0e", "#2ca02c"),
                   conf.int = TRUE, censor = FALSE, legend.title = "", legend.labs = strata_labs,
                   xlab = "Time (years)", ylab = "Overall survival", title = ttl,
                   legend = c(0.82, 0.93), risk.table = TRUE, risk.table.y.text = FALSE,
                   fontsize = 3.2)
  gg$plot <- gg$plot +
    theme(legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
          legend.key.size = unit(0.8, "lines"), legend.text = element_text(size = 9),
          plot.margin = margin(5.5, 5.5, 0, 5.5))
  gg$table <- gg$table + labs(title = NULL, x = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(0, 5.5, 5.5, 5.5))
  cowplot::plot_grid(gg$plot, gg$table, ncol = 1, align = "v", axis = "lr", rel_heights = c(3, 0.9))
}
for (nm in names(COH)) {
  D <- COH[[nm]]
  save_fig(km_panel(D, eln_grp(D), sprintf("%s (n=%d) -- ELN-2022", nm, nrow(D))),
           fig_path(paste0("Fig1_KM_", nm, "_ELN")), width = 5.5, height = 5)
  save_fig(km_panel(D, rsf_grp(risk_of(D)), sprintf("%s (n=%d) -- CCM-RSF", nm, nrow(D))),
           fig_path(paste0("Fig1_KM_", nm, "_RSF")), width = 5.5, height = 5)
}

## ============================ FIGURE 2 : c-index forest ======================
if (!file.exists(out_path("table2a"))) stop("Run '04_table2a_discrimination.R' first (Figure 2 uses its CIs)")
t2a <- read.csv(out_path("table2a"))
fdat <- bind_rows(lapply(names(COH), function(nm) {
  r <- t2a[t2a$cohort == nm, ]
  data.frame(cohort = sprintf("%s (n=%d)", nm, r$n), model = c("ELN-2022", "CCM-RSF"),
             c = c(r$c_ELN, r$c_RSF), lo = c(r$ELN_lo, r$RSF_lo), hi = c(r$ELN_hi, r$RSF_hi))
})) %>% mutate(row = paste(cohort, model, sep = ": "), row = factor(row, levels = rev(unique(row))))
pf <- ggplot(fdat, aes(c, row, color = model)) +
  geom_point(size = 2.5) + geom_errorbar(aes(xmin = lo, xmax = hi), width = .2, orientation = "y") +
  geom_vline(xintercept = .5, linetype = 3, color = "grey") +
  scale_color_manual(values = c("ELN-2022" = "#ff7f0e", "CCM-RSF" = "#1f77b4")) +
  labs(x = sprintf("Harrell c-index (95%% CI, %d bootstrap resamples)", N_BOOT), y = NULL, color = NULL) +
  theme_minimal()
save_fig(pf, fig_path("Fig2_cindex_forest"), width = 6.4, height = 1.1 * length(COH))

## ============================ FIGURE 3 : calibration (deciles) ===============
calib_panel <- function(D, ttl, h = 2, q = 10) {
  sf <- predict(rsf, newdata = D)$survival
  j  <- which.min(abs(rsf$time.interest - h))
  pred <- sf[, j]
  bin <- cut(pred, breaks = quantile(pred, seq(0, 1, 1/q), na.rm = TRUE), include.lowest = TRUE, labels = FALSE)
  rows <- bind_rows(lapply(sort(unique(bin)), function(b) {
    m <- bin == b
    obs <- summary(survfit(Surv(time, status) ~ 1, data = D[m, ]), times = h, extend = TRUE)$surv
    data.frame(pred = mean(pred[m]), obs = obs)
  }))
  ggplot(rows, aes(pred, obs)) + geom_abline(linetype = 2) +
    geom_point(color = "#1f77b4", size = 2) + xlim(0, 1) + ylim(0, 1) +
    labs(x = "Predicted 2-yr survival", y = "Observed (KM) 2-yr survival", title = ttl) + theme_minimal()
}
pc <- lapply(names(COH), function(nm) calib_panel(COH[[nm]], sprintf("%s (n=%d)", nm, nrow(COH[[nm]]))))
save_fig(arrangeGrob(grobs = pc, ncol = length(pc)), fig_path("Fig3_calibration"), width = 4 * length(pc), height = 4)

## ============================ FIGURE 4 : importance (NTUH, unchanged) ========
imp <- sort(rsf$importance, decreasing = TRUE)[1:15]
names(imp) <- names(VMAP)[match(names(imp), VMAP)]
idat <- data.frame(feature = factor(names(imp), levels = rev(names(imp))), vimp = as.numeric(imp))
pi <- ggplot(idat, aes(vimp, feature)) + geom_col(fill = "#1f77b4") +
  scale_y_discrete(labels = function(x) feat_label(x)) +
  labs(x = "Permutation importance", y = NULL, title = "CCM-RSF feature importance (NTUH)") +
  theme_minimal() + theme(plot.title = element_text(size = 10), axis.text.y = ggtext::element_markdown(), axis.text.y.left = ggtext::element_markdown())
save_fig(pi, fig_path("Fig4_importance"), width = 5.4, height = 5)

## ============================ TABLE 4 : within-ELN ===========================
## Per ELN-2022 group: n, 2-yr OS and median OS (years) for each RSF group, and the
## RSF-adverse vs RSF-favorable Cox HR (Wald p) and log-rank p, plus the HR per SD of the continuous
## CCM-RSF score within that ELN group (SD of the whole cohort; Wald p). "Overall" row = all ELN groups.
mk_row <- function(D, g, m, eln_label, nm, z) {
  hi <- m & g == 0; mid <- m & g == 1; lo <- m & g == 2
  ## adverse vs favorable is only estimated when both groups have >= 3 patients;
  ## the HR per SD of the continuous score is estimated regardless
  sub <- D[m & g != 1, ]; x <- as.integer(g[m & g != 1] == 0)
  cox_fit <- if (sum(hi) < 3 || sum(lo) < 3) NULL else
    tryCatch(coxph(Surv(time, status) ~ x, data = sub), error = function(e) NULL)
  hr   <- if (is.null(cox_fit)) NA else exp(coef(cox_fit)[[1]])
  hr_p <- if (is.null(cox_fit)) NA else summary(cox_fit)$coefficients[1, "Pr(>|z|)"]
  lr   <- if (is.null(cox_fit)) NULL else tryCatch(survdiff(Surv(time, status) ~ x, data = sub), error = function(e) NULL)
  lr_p <- if (is.null(lr)) NA else 1 - pchisq(lr$chisq, 1)
  sz   <- summary(coxph(Surv(time, status) ~ z, data = data.frame(D[m, c("time", "status")], z = z[m])))
  os2 <- function(idx) if (sum(idx) == 0) NA else
    summary(survfit(Surv(time, status) ~ 1, data = D[idx, ]), times = 2, extend = TRUE)$surv
  med_os <- function(idx) if (sum(idx) == 0) NA else
    summary(survfit(Surv(time, status) ~ 1, data = D[idx, ]))$table[["median"]]
  data.frame(cohort = nm, ELN = eln_label,
             RSFadv_n = sum(hi), RSFadv_OS2 = round(os2(hi), 2), RSFadv_medOS = round(med_os(hi), 2),
             RSFint_n = sum(mid), RSFint_OS2 = round(os2(mid), 2), RSFint_medOS = round(med_os(mid), 2),
             RSFfav_n = sum(lo), RSFfav_OS2 = round(os2(lo), 2), RSFfav_medOS = round(med_os(lo), 2),
             adv_to_fav_HR = round(hr, 2), adv_to_fav_p = signif(hr_p, 3), logrank_p = signif(lr_p, 3),
             HR_perSD = round(sz$conf.int[1, 1], 2), HR_perSD_lo = round(sz$conf.int[1, 3], 2),
             HR_perSD_hi = round(sz$conf.int[1, 4], 2), HR_perSD_p = signif(sz$coefficients[1, 5], 3))
}
build_wln <- function(D, g, eg, nm, r) {
  z <- (r - mean(r)) / sd(r)
  rbind(do.call(rbind, lapply(0:2, function(e) mk_row(D, g, eg == e, GLAB[as.character(e)], nm, z))),
        mk_row(D, g, rep(TRUE, nrow(D)), "Overall (all ELN)", nm, z))
}
D_ntuh <- data.frame(time = nt$OSD_no_censor_2021 / 12, status = nt$OS_coding_no_censor_2021)
wln <- bind_rows(
  build_wln(D_ntuh, rsf_grp(rsf$predicted), 3 - nt$ELN_2022, "NTUH", rsf$predicted),
  bind_rows(lapply(names(COH), function(nm) {
    D <- COH[[nm]]; r <- risk_of(D); build_wln(D, rsf_grp(r), eln_grp(D), nm, r)
  })))
print(wln, row.names = FALSE)
write.csv(wln, out_path("within_eln"), row.names = FALSE)
cat("\nWrote", out_path("within_eln"), "and", out_path("crosstab"), "\n")
