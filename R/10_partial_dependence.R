## Partial dependence of predicted 2-year OS on the leading predictors: age, WBC, BM blast %,
## TP53 mutation and CEBPA bZIP mutation
## (NTUH training cohort only, so the external-cohort exclusion does not affect it).
## Output: figures/Fig_PDP_53feat.png

suppressPackageStartupMessages({ library(gridExtra) })
source("00_setup.R")
set.seed(0)

nt <- load_ntuh()
NT <- ntuh_frame(nt)
jt <- which.min(abs(rsf$time.interest - 2))
surv2 <- function(D) predict(rsf, D)$survival[, jt]

pdp_cont <- function(varname, grid) {
  v <- VMAP[varname]
  data.frame(x = grid, surv = sapply(grid, function(val) { D <- NT; D[[v]] <- val; mean(surv2(D)) }))
}
g_age <- round(quantile(nt$Age, seq(.05, .95, length = 25)))
g_wbc <- quantile(nt$WBC, seq(.05, .95, length = 25))
g_bm  <- quantile(nt$BM_Blast, seq(.05, .95, length = 25))
d_age <- pdp_cont("Age", g_age); d_wbc <- pdp_cont("WBC", g_wbc); d_bm <- pdp_cont("BM_Blast", g_bm)
mkline <- function(d, xlab) ggplot(d, aes(x, surv)) + geom_line(color = "#1f77b4", linewidth = 1) +
  ylim(0, 1) + labs(x = xlab, y = "Predicted 2-yr survival") + theme_minimal()

binary_pdp <- function(varname, lo_lab, hi_lab) {
  v <- VMAP[varname]; D0 <- NT; D0[[v]] <- 0; D1 <- NT; D1[[v]] <- 1
  data.frame(grp = factor(c(lo_lab, hi_lab), levels = c(lo_lab, hi_lab)), surv = c(mean(surv2(D0)), mean(surv2(D1))))
}
mk_bar <- function(d, varname, ttl) ggplot(d, aes(grp, surv, fill = grp)) + geom_col(width = .6) +
  scale_fill_manual(values = setNames(if (d$surv[1] >= d$surv[2]) c("#2ca02c", "#d62728") else c("#d62728", "#2ca02c"),
                                      levels(d$grp)), guide = "none") +   # green = better predicted survival
  ylim(0, 1) + labs(x = varname, y = "Predicted 2-yr survival", title = ttl) +
  geom_text(aes(label = sprintf("%.2f", surv)), vjust = -0.4) + theme_minimal()
## titles use plotmath: ggtext drops the space after an italic word ending in a capital letter
d_tp    <- binary_pdp("TP53", "wild-type", "mutated")
d_cebpa <- binary_pdp("CEBPA_bZIP", "wild-type", "bZIP in-frame")

save_fig(arrangeGrob(mkline(d_age, "Age (years)"),
                     mkline(d_wbc, "WBC (/uL)") + scale_x_log10(labels = scales::label_comma()),
                     mkline(d_bm, "Bone marrow blasts (%)"),
                     mk_bar(d_tp, quote(italic("TP53")), quote(italic("TP53") ~ "mutation")),
                     mk_bar(d_cebpa, quote(italic("CEBPA")), quote(italic("CEBPA") ~ "bZIP in-frame mutation")), ncol = 3,
                     top = "Partial dependence: predicted 2-year overall survival (CCM-RSF)"),
         fig_path("Fig_PDP"), width = 14, height = 8)
cat(sprintf("Age %d->%d: 2-yr surv %.2f -> %.2f\n", min(g_age), max(g_age), d_age$surv[1], tail(d_age$surv, 1)))
cat(sprintf("TP53 wt vs mut: %.2f vs %.2f\n", d_tp$surv[1], d_tp$surv[2]))
cat(sprintf("CEBPA bZIP wt vs mut: %.2f vs %.2f\n", d_cebpa$surv[1], d_cebpa$surv[2]))
cat(sprintf("BM blast %.0f%% -> %.0f%%: 2-yr surv %.2f -> %.2f\n", min(g_bm), max(g_bm), d_bm$surv[1], tail(d_bm$surv, 1)))
cat("WBC grid (/uL) and predicted 2-yr survival:\n"); print(round(data.frame(WBC = d_wbc$x, surv = d_wbc$surv), 3), row.names = FALSE)
