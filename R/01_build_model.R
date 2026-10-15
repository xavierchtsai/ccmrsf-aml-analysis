## Fit and save the frozen 53-feature RSF on NTUH (n=801).
## Specification: 400 trees, nodesize 10, mtry = floor(sqrt(53)), log-rank splitting with
## 10 random split points, set.seed(0) + seed = -1 (reproducible).
##
## Output: model/rsf_model_NTUH_53feat.rds

source("00_setup.R")
set.seed(0)

nt <- load_ntuh()
NT <- ntuh_frame(nt)
fml <- as.formula(paste("Surv(time,status) ~", paste(VMAP[FEATS], collapse = " + ")))
rsf <- rfsrc(fml, NT, ntree = 400, nodesize = 10, mtry = floor(sqrt(length(FEATS))),
             splitrule = "logrank", importance = "permute", seed = -1)
print(rsf)

cut <- quantile(rsf$predicted, c(1/3, 2/3))
cat(sprintf("\nNTUH tertile cutpoints (locked): %.2f / %.2f\n", cut[1], cut[2]))
cat(sprintf("NTUH OOB c-index: %.3f\n", 1 - rsf$err.rate[rsf$ntree]))

imp <- sort(rsf$importance, decreasing = TRUE)[1:10]
names(imp) <- names(VMAP)[match(names(imp), VMAP)]
cat("\nTop-10 permutation importance:\n"); print(round(imp, 5))

saveRDS(rsf, MODEL_RDS)
cat("\nSaved model ->", MODEL_RDS, "\n")
