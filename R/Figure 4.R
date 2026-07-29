######################################################################
## Figure 4 - Organ-specific immune-cell analysis (lung vs. gut)
##   (A) Stacked bar of immune-cell fractions per ROI, faceted by tissue
##       (immune-high regions only).
##   (B) Box plot of Macrophage M1 fraction by tissue.
##   (C) PLS-R observed vs. predicted M1 macrophage fraction across all
##       53 ROIs, using the other 21 immune cell types as predictors.
##   (D) VIP score vs. PLS coefficient for the M1 PLS-R model.
##
## Data:
##   absolute cell fraction.csv   - immune-high ROIs only (panels A, B)
##   whole cell fraction.xlsx     - all 53 ROIs           (panels C, D)
######################################################################

suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(readr)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(pls)
  library(patchwork)
})

## CV folds are random; fix for a reproducible figure.
set.seed(42)

## ---- Standard PLS VIP helper ---------------------------------------
## VIP_j = sqrt(p * sum_a (SSY_a * w_aj^2 / sum_k w_ak^2) / sum_a SSY_a)
## (Wold et al., textbook PLS VIP). mean(VIP^2) = 1 so VIP > 1 is the
## conventional "above-average importance" cutoff.
##
## The same definition appears in Figures 1, 2, 3, 4 and S2 - one VIP
## computation across the whole protocol.
pls_vip <- function(model, ncomp) {
  W <- model$loading.weights[, 1:ncomp, drop = FALSE]
  T <- model$scores[,         1:ncomp, drop = FALSE]
  Q <- model$Yloadings[,      1:ncomp, drop = FALSE]
  p <- nrow(W)
  SSY <- as.numeric(Q)^2 * colSums(T * T)
  Wn  <- sweep(W^2, 2, colSums(W^2), "/")
  as.numeric(sqrt(p * (Wn %*% SSY) / sum(SSY)))
}

## ---- Paths (relative to repo root, resolved via here::here) --------
ABS_FILE <- here::here("Cell Rep data", "absolute cell fraction.csv")
ALL_FILE <- here::here("Cell Rep data", "whole cell fraction.xlsx")
OUT_PDF  <- here::here("figures", "Figure 4.pdf")
OUT_JPG  <- here::here("figures", "Figure 4.jpg")

abs_d <- read_csv(ABS_FILE, show_col_types = FALSE)
all_d <- read_excel(ALL_FILE)

## ---- "GI" to "gut" -------------------------------------------------
## "GI" is a value inside the Tissue column of both source files, not a
## display string, and the manuscript says "gut" in every one of its
## mentions. Recode both data frames here, before any select, pivot,
## facet or factor call.
##
## Renaming only the hard-coded colour, level and shape vectors further
## down would leave the data still saying "GI", and then ggplot finds no
## match: the panel A facet strip would still read "GI", panel B's boxes
## would fall back to grey, and scale_shape_manual would hand panel C NA
## shapes, so every point would disappear.
##
## The source CSV and XLSX are not modified.
abs_d$Tissue <- ifelse(abs_d$Tissue == "GI", "gut", abs_d$Tissue)
all_d$Tissue <- ifelse(all_d$Tissue == "GI", "gut", all_d$Tissue)

cell_cols <- setdiff(
  colnames(abs_d),
  c("Tissue", "ROI", "Mixture", "P-value", "Correlation",
    "RMSE", "Absolute score (sig.score)")
)

## ---- Panels A and B only: row-normalize to true fractions ----------
## CIBERSORTx absolute mode returns scores whose row sums vary, 1.259 to
## 2.058 across these ROIs, and equal the "Absolute score (sig.score)"
## column exactly. The raw values are therefore not fractions. Dividing
## each ROI by its own row total makes panel A's stacked bars sum to 1
## and makes panel B's M1 a true fraction of that ROI's immune cells,
## which is what both axis labels already say.
##
## This applies to abs_d ONLY. all_d, which feeds the PLS-R in panels C
## and D, is left on the absolute scale. Once the 22 cell
## types sum to 1, M1 equals 1 minus the other 21 exactly, so the
## regression would be fitting an identity rather than a relationship.
abs_d[cell_cols] <- sweep(as.matrix(abs_d[cell_cols]), 1,
                          rowSums(abs_d[cell_cols]), "/")

## ---- Panel A: stacked bar of immune-cell fractions -----------------
long_d <- abs_d %>%
  select(Tissue, ROI, all_of(cell_cols)) %>%
  pivot_longer(-c(Tissue, ROI),
               names_to = "Cell_Type", values_to = "Fraction")

palette22 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
               "#1a55ff", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
               "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
               "#6b6ecf", "#b5cf6b")

p_A <- ggplot(long_d,
              aes(factor(ROI), Fraction, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette22) +
  facet_wrap(~ Tissue, scales = "free_x") +
  labs(x = "ROI #", y = "Fraction of Immune Cells",
       title = "Cell distribution per ROI (immune-high)") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_blank(),
        legend.text  = element_text(size = 5.5),
        legend.key.size = unit(6, "pt"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 10)) +
  guides(fill = guide_legend(nrow = 5))

## ---- Panel B: M1 macrophage box plot -------------------------------
p_B <- ggplot(abs_d, aes(Tissue, `Macrophages M1`, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85) +
  geom_jitter(width = 0.15, color = "black", size = 1.6, alpha = 0.7) +
  scale_fill_manual(values = c("lung" = "#FFC0CB", "gut" = "#ADD8E6")) +
  labs(title = "M1 fraction by tissue",
       x = "Tissue", y = "M1 Fraction") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = "none")

## ---- Fit PLS-R for panels C, D -------------------------------------
X_all <- all_d %>%
  select(-c("ROI", "Tissue", "Mixture",
            "P-value", "Correlation", "RMSE",
            "Absolute score (sig.score)",
            "Macrophages M1"))
Y_all <- all_d[["Macrophages M1"]]

ncomp_max <- 10
pls_m1 <- plsr(Y_all ~ ., data = data.frame(Y_all = Y_all, X_all),
               ncomp = ncomp_max, validation = "CV", scale = TRUE)

ncomp_use <- 8
preds    <- as.numeric(predict(pls_m1, ncomp = ncomp_use))
cv_preds <- as.numeric(pls_m1$validation$pred[, 1, ncomp_use])

## PLS model-quality metrics (Wold et al., 2001)
R2X_m1 <- sum(explvar(pls_m1)[1:ncomp_use]) / 100
R2Y_m1 <- 1 - sum((Y_all - preds)^2)    / sum((Y_all - mean(Y_all))^2)
Q2_m1  <- 1 - sum((Y_all - cv_preds)^2) / sum((Y_all - mean(Y_all))^2)

cat(sprintf("M1 PLS-R: R2X=%.3f  R2Y=%.3f  Q2=%.3f  (ncomp=%d)\n",
            R2X_m1, R2Y_m1, Q2_m1, ncomp_use))

## ---- Panel C: observed vs. predicted M1 ----------------------------
metrics_label_m1 <- sprintf("R²X=%.3f  R²Y=%.3f\nQ²=%.3f",
                            R2X_m1, R2Y_m1, Q2_m1)

df_C <- data.frame(Observed  = Y_all,
                   Predicted = preds,
                   Tissue    = factor(all_d$Tissue,
                                      levels = c("lung", "gut")))

## Panels C and D run on the absolute scale while panels A and B are
## normalized, so the axes say so. Without this, the same ROI appears as
## two different numbers in panel B and panel C with nothing telling the
## reader they are different quantities.
p_C <- ggplot(df_C, aes(Observed, Predicted, shape = Tissue)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.6) +
  geom_point(size = 1.5, stroke = 0.5) +
  scale_shape_manual(values = c("lung" = 3, "gut" = 1)) +
  annotate("text", x = Inf, y = -Inf, label = metrics_label_m1,
           hjust = 1.1, vjust = -0.5, size = 3) +
  labs(title = "M1: actual vs. predicted",
       x = "Observed M1 (a.u.)",
       y = "Predicted M1 (a.u.)") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = c(0.12, 0.88),
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA, color = NA))

## ---- Panel D: VIP-coefficient plot ---------------------------------
## VIP via the shared `pls_vip` helper at the top of this file.
vip_vals <- pls_vip(pls_m1, ncomp_use)
coef_v   <- coef(pls_m1, ncomp = ncomp_use)

dat_D <- data.frame(
  Cell       = colnames(X_all),
  Importance = vip_vals,
  Coef       = as.numeric(coef_v)
)
dat_D$above <- dat_D$Importance > 1.0    # standard VIP cutoff (1.0 in PLS literature)

p_D <- ggplot(dat_D, aes(Coef, Importance)) +
  geom_hline(yintercept = 1.0,  linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0,    linetype = "dashed", color = "black") +
  geom_point(data = subset(dat_D, !above), color = "grey60", size = 1.6) +
  geom_point(data = subset(dat_D,  above), color = "red",    size = 2.0) +
  geom_text_repel(data = subset(dat_D, above),
                  aes(label = Cell),
                  size = 2.5, box.padding = 0.25,
                  max.overlaps = Inf,
                  segment.color = "grey50", segment.size = 0.3) +
  labs(title = "VIP vs. Coefficient (M1 PLS-R)",
       x = "Coefficient", y = "VIP Score") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

## ---- Compose layout and write outputs ------------------------------
mid4 <- p_B | p_C
fig4 <- (p_A / mid4 / p_D) +
  plot_layout(heights = c(1.3, 1, 1)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) &
  theme(plot.tag = element_text(face = "bold", size = 13))
pdf(OUT_PDF, width = 6.7, height = 8, family = "Helvetica")
print(fig4)
dev.off()

ggsave(OUT_JPG, fig4, width = 6.7, height = 8, dpi = 300)

cat("Wrote ", OUT_PDF, "\n", sep = "")
cat("Wrote ", OUT_JPG, "\n", sep = "")
