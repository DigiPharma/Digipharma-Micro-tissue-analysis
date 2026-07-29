######################################################################
## Figure S3 - Permutation test of the worked example 1 PLS-R model
##
## Tests whether the cross-validated Q2 of worked example 1 (PLS-R of
## CD8A on immune-rich GBM ROIs, R/Figure 1.R) is better than chance.
## The CD8A labels are shuffled 1,000 times, the model is refit on each
## shuffle, and the observed Q2 is compared with the resulting null.
##
## The model is reused from R/Figure 1.R without modification: same
## data, same >=25% CD45+ filter, same Y = CD8A, same X (all measured
## proteins except CD8A and CD3), same ncomp_use = 8, same 10-fold CV.
##
## (A) Null distribution of cross-validated Q2 under permuted labels,
##     with the observed Q2 marked.
##
## Data: Raw data.xlsx, sheet "GBM" (Lu et al., 2021), filtered to
##   Ratio of CD45+ cells >= 0.25.
##
## `Y ~ .` requires the response column in the data frame to be named Y.
## A mismatched name leaves the permuted response among the predictors
## and returns Q2 near 1 on every permutation.
##
## set.seed(42) fixes the CV fold draw and the permutation draws.
######################################################################

suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(dplyr)
  library(pls)
  library(ggplot2)
})

## ---- Paths (relative to repo root, resolved via here::here) --------
DATA    <- here::here("Nat commun data", "Raw data.xlsx")
OUT_PDF <- here::here("figures", "Figure S3.pdf")
OUT_JPG <- here::here("figures", "Figure S3.jpg")
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

## ---- 1. Load and filter to immune-rich (>=25% CD45+) ----------------
## Identical to R/Figure 1.R.
gbm <- read_excel(DATA, sheet = "GBM")
gbm_ir <- gbm %>% filter(`Ratio of CD45+ cells` >= 0.25)

X <- gbm_ir %>% select(-c(1:5)) %>% select(-CD8A, -CD3)
Y <- gbm_ir$CD8A

ncomp_max <- 20
ncomp_use <- 8

## ---- 2. Observed model ---------------------------------------------
set.seed(42)
plsr.mod <- plsr(Y ~ ., data = data.frame(Y = Y, X), ncomp = ncomp_max,
                 validation = "CV", scale = TRUE)

preds    <- as.numeric(predict(plsr.mod, ncomp = ncomp_use))
cv_preds <- as.numeric(plsr.mod$validation$pred[, 1, ncomp_use])

R2X <- sum(explvar(plsr.mod)[1:ncomp_use]) / 100
R2Y <- 1 - sum((Y - preds)^2)    / sum((Y - mean(Y))^2)
Q2  <- 1 - sum((Y - cv_preds)^2) / sum((Y - mean(Y))^2)

cat(sprintf("Observed model:  R2X=%.3f  R2Y=%.3f  Q2=%.3f  (ncomp=%d)\n",
            R2X, R2Y, Q2, ncomp_use))
cat(sprintf("  nROI = %d,  nPred = %d\n", nrow(X), ncol(X)))

## ---- 3. Permutation null -------------------------------------------
set.seed(42)
n_perm  <- 1000
obs_Q2  <- Q2
perm_Q2 <- replicate(n_perm, {
  Y_perm <- sample(Y)
  m   <- plsr(Y ~ ., data = data.frame(Y = Y_perm, X),
              ncomp = ncomp_use, validation = "CV", scale = TRUE)
  cvp <- as.numeric(m$validation$pred[, 1, ncomp_use])
  1 - sum((Y_perm - cvp)^2) / sum((Y_perm - mean(Y_perm))^2)
})
p_perm <- (1 + sum(perm_Q2 >= obs_Q2)) / (1 + n_perm)

n_ge <- sum(perm_Q2 >= obs_Q2)
cat("\n--- Permutation test ---\n")
cat(sprintf("observed Q2      = %.4f\n", obs_Q2))
cat(sprintf("mean(perm_Q2)    = %.4f\n", mean(perm_Q2)))
cat(sprintf("median(perm_Q2)  = %.4f\n", median(perm_Q2)))
cat(sprintf("max(perm_Q2)     = %.4f\n", max(perm_Q2)))
cat(sprintf("# perms >= obs   = %d of %d\n", n_ge, n_perm))
cat(sprintf("p_perm           = %.4f\n", p_perm))

## ---- 4. Figure: null distribution with the observed Q2 --------------
perm_df <- data.frame(Q2 = perm_Q2)

## Plain ">=" rather than the glyph: the Helvetica PDF device cannot
## encode U+2265 and silently substitutes.
metrics_label <- sprintf(
  "observed Q² = %.3f\nnull mean = %.3f\nnull max = %.3f\np = %.3f  (%d/%d >= observed)",
  obs_Q2, mean(perm_Q2), max(perm_Q2), p_perm, n_ge, n_perm)

fig <- ggplot(perm_df, aes(Q2)) +
  geom_histogram(bins = 40, fill = "grey80", colour = "white",
                 linewidth = 0.2) +
  geom_vline(xintercept = obs_Q2, colour = "red", linewidth = 0.9) +
  annotate("text", x = -Inf, y = Inf, label = metrics_label,
           hjust = -0.06, vjust = 1.2, size = 3.2) +
  ## coord_cartesian, not scale_x_continuous(limits): the latter filters
  ## the data before binning.
  coord_cartesian(xlim = range(c(perm_Q2, obs_Q2)) + c(-0.15, 0.15)) +
  labs(x = "Q²", y = "Count") +
  theme_classic(base_size = 11)

## ---- 5. Write vector PDF and JPG -----------------------------------
pdf(OUT_PDF, width = 6.7, height = 4.5, family = "Helvetica")
print(fig)
invisible(dev.off())
cat("\nWrote ", OUT_PDF, "\n", sep = "")
ggsave(OUT_JPG, fig, width = 6.7, height = 4.5, dpi = 300)
cat("Wrote ", OUT_JPG, "\n", sep = "")
