## Export the frozen 53-feature CCM-RSF for the static web calculator
## (https://github.com/xavierchtsai/ccmrsf-calculator).
## The forest is written as flat typed arrays so the browser reproduces predict() exactly:
##   - tree topology in preorder (split variable, split value, index of the right child);
##     a case goes left when x <= split value
##   - terminal-node mortality (tnMORT) -> ensemble risk score = mean over the 400 trees
##   - terminal-node survival (tnSURV) on a quarterly grid 0-10 years, stored as uint16
## Also writes a test file (external patients with R predictions) used to check the JavaScript.
## Output: web/model.js, web/tests/test_cases.json

source("00_setup.R")
library(jsonlite)

WEB <- file.path(ROOT, "web")
dir.create(file.path(WEB, "tests"), showWarnings = FALSE, recursive = TRUE)

f  <- rsf$forest; na <- f$nativeArray; tn <- f$nativeArrayTNDS
ti <- rsf$time.interest; nT <- length(ti); leaf <- f$leafCount
grid <- seq(0, 10, by = 0.25)
gidx <- vapply(grid, function(t) sum(ti <= t), integer(1))          # step function: last time <= t

## ---- tree topology ----
n <- nrow(na)
var <- as.integer(na$parmID); split <- ifelse(is.na(na$contPT), 0, na$contPT)
right <- integer(n); leafid <- integer(n) - 1L; start <- integer(rsf$ntree)
leaf_off <- c(0L, cumsum(leaf))[seq_len(rsf$ntree)]
pos <- 0L
for (b in seq_len(rsf$ntree)) {
  rows <- which(na$treeID == b); p <- var[rows]; m <- length(rows); s <- integer(m)
  size <- function(i) { if (p[i] == 0) { s[i] <<- 1L; return(1L) }
    l <- size(i + 1L); r <- size(i + 1L + l); s[i] <<- 1L + l + r; s[i] }
  size(1L)
  for (i in seq_len(m)) {
    g <- rows[i] - 1L                                                   # 0-based global index
    if (p[i] > 0) right[rows[i]] <- g + 1L + s[i + 1L]
    else leafid[rows[i]] <- leaf_off[b] + na$nodeID[rows[i]] - 1L        # 0-based global leaf
  }
  start[b] <- rows[1] - 1L
}

## ---- terminal-node estimates ----
L <- length(tn$tnMORT)
surv_all <- matrix(tn$tnSURV, nrow = nT, ncol = L)                      # column = leaf
surv_grid <- rbind(1, surv_all[gidx[gidx > 0], , drop = FALSE])         # grid 0 -> survival 1
stopifnot(nrow(surv_grid) == length(grid))
q16 <- as.integer(round(pmin(pmax(surv_grid, 0), 1) * 65535))

b64 <- function(x, type, size) base64_enc(writeBin(x, raw(), size = size, endian = "little"))
b64i <- function(x, size) base64_enc(writeBin(as.integer(x), raw(), size = size, endian = "little"))

## ---- reference curves: NTUH Kaplan-Meier by CCM-RSF group ----
nt <- load_ntuh(); g <- rsf_grp(rsf$predicted)
ref <- lapply(c(adverse = 0, intermediate = 1, favorable = 2), function(k) {
  fit <- survfit(Surv(nt$OSD_no_censor_2021 / 12, nt$OS_coding_no_censor_2021) ~ 1, subset = g == k)
  round(summary(fit, times = grid, extend = TRUE)$surv, 4) })

model <- list(
  name = "CCM-RSF", version = "53feat", ntree = rsf$ntree, nleaf = L,
  features = FEATS, cutoffs = unname(round(cut, 6)), grid = grid, ref = ref,
  var = b64i(var, 1), split = b64(split, "double", 8), right = b64i(right, 4),
  leaf = b64i(leafid, 4), start = b64i(start, 4),
  mort = b64(tn$tnMORT, "double", 8), surv = b64i(q16, 2))
writeLines(paste0("window.CCMRSF_MODEL = ", toJSON(model, auto_unbox = TRUE, digits = NA), ";"),
           file.path(WEB, "model.js"))

## ---- test cases: 100 patients per external cohort, with R predictions ----
COH <- analysis_frames(load_external(verbose = FALSE))
set.seed(2)
D <- do.call(rbind, lapply(COH, function(d) d[sample(nrow(d), 100), ]))
pr <- predict(rsf, newdata = D)
at <- function(t) pr$survival[, max(which(ti <= t))]
X <- D[, VMAP]; colnames(X) <- FEATS
tests <- list(x = unname(as.matrix(X)), mort = pr$predicted, group = rsf_grp(pr$predicted),
              s1 = at(1), s2 = at(2), s5 = at(5))
writeLines(toJSON(tests, digits = NA), file.path(WEB, "tests", "test_cases.json"))
cat("Wrote", file.path(WEB, "model.js"), sprintf("(%.1f MB)", file.size(file.path(WEB, "model.js")) / 1e6), "\n")
cat("cutoffs:", cut, "\n")
