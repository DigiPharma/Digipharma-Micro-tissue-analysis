######################################################################
## Figure 2 - PLS-DA prediction model for treatment response
##
## Derives the 5-protein prediction formula directly from the
## PLS-DA model. All panels are generated from the fitted models.
##
## (A) VIP score vs. PLS-DA coefficient (11-protein model).
##     The 5 proteins above VIP > 1 are identified for the
##     final prediction model.
## (B) ROI-level box plot of 5-protein model predictions, with
##     the model-derived prediction formula displayed.
## (C) ROI-level ROC curve and AUC.
## (D) Tissue-level box plot: per-ROI predictions averaged within
##     each Sample ID (n = 8 mNR and 4 mR).
## (E) Tissue-level ROC curve and AUC.
##
## Data:
##   Raw data_plsda round 1.xlsx - 11 proteins (input to Panel A)
##   Raw data_plsda round 2.xlsx - 5 proteins + metadata (input to
##     the 5-protein PLS-DA; the pre-computed Predictions column
##     is NOT used)
######################################################################

suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(dplyr)
  library(pls)
  library(ggplot2)
  library(ggrepel)
  library(pROC)
  library(patchwork)
})

## ---- Standard PLS VIP helper ---------------------------------------
pls_vip <- function(model, ncomp) {
  W <- model$loading.weights[, 1:ncomp, drop = FALSE]
  T <- model$scores[,         1:ncomp, drop = FALSE]
  Q <- model$Yloadings[,      1:ncomp, drop = FALSE]
  p <- nrow(W)
  SSY <- as.numeric(Q)^2 * colSums(T * T)
  Wn  <- sweep(W^2, 2, colSums(W^2), "/")
  as.numeric(sqrt(p * (Wn %*% SSY) / sum(SSY)))
}

## ---- Helper: ROC ggplot panel --------------------------------------
roc_panel <- function(roc_obj, title, auc_x = 0.40, auc_y = 0.08) {
  roc_df <- data.frame(spec = roc_obj$specificities,
                       sens = roc_obj$sensitivities) %>%
    arrange(desc(spec), sens)
  auc_v <- as.numeric(auc(roc_obj))
  ggplot(roc_df, aes(spec, sens)) +
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1),
                 inherit.aes = FALSE,
                 color = "grey60", linewidth = 0.5) +
    geom_step(color = "red", linewidth = 1.4, direction = "vh") +
    scale_x_reverse(limits = c(1, 0), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    annotate("text", x = auc_x, y = auc_y,
             label = sprintf("AUC = %.3f", auc_v)) +
    labs(title = title, x = "Specificity", y = "Sensitivity") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
}

## ---- Paths ---------------------------------------------------------
ROUND1_FILE <- here::here("Nat commun data", "Raw data_plsda round 1.xlsx")
ROUND2_FILE <- here::here("Nat commun data", "Raw data_plsda round 2.xlsx")
OUT_PDF     <- here::here("figures", "Figure 2.pdf")
OUT_JPG     <- here::here("figures", "Figure 2.jpg")

## ====================================================================
## Panel A: PLS-DA on 11 proteins → identifies 5 with VIP > 1
## ====================================================================
r1 <- read_excel(ROUND1_FILE, sheet = "Sheet1")
prots_r1 <- c("AKT", "CD11c", "CD163", "CD44", "CD66B", "FoxP3",
              "HLA-DR", "pAKT", "PTEN", "STAT3", "pSTAT3")
X_r1 <- as.data.frame(r1[, prots_r1])
Y_r1 <- r1$pls_code

plsda_r1 <- plsr(Y_r1 ~ ., data = data.frame(Y_r1 = Y_r1, X_r1),
                 ncomp = 3, validation = "CV", scale = TRUE)

vip_vals <- pls_vip(plsda_r1, 3)
coefs_r1 <- coef(plsda_r1, ncomp = 3)

vip_dat <- data.frame(
  Gene       = colnames(X_r1),
  Importance = vip_vals,
  Coef       = as.numeric(coefs_r1)
)
selected <- c("CD11c", "CD163", "CD44", "CD66B", "PTEN")
vip_dat$above <- vip_dat$Gene %in% selected

p_A <- ggplot(vip_dat, aes(Coef, Importance)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0,   linetype = "dashed", color = "black") +
  geom_point(data = subset(vip_dat, !above), color = "grey60", size = 1.8) +
  geom_point(data = subset(vip_dat,  above), color = "red",    size = 2.4) +
  geom_text_repel(data = subset(vip_dat, above),
                  aes(label = Gene),
                  size = 3, box.padding = 0.35,
                  segment.color = "grey40", segment.size = 0.3,
                  max.overlaps = Inf) +
  labs(title = "VIP Scores vs. Coefficients (11-protein model)",
       x = "Coefficient", y = "VIP Score") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0, size = 10))

## ====================================================================
## 5-protein PLS-DA: fit model and derive prediction formula
## ====================================================================
r2 <- read_excel(ROUND2_FILE, sheet = "Sheet1")
prots_5 <- c("CD11c", "CD163", "CD44", "CD66B", "PTEN")
X_5 <- as.data.frame(r2[, prots_5])
Y_5 <- r2$pls_code

## Fit PLS-DA on the 5 proteins with scaling
plsda_5 <- plsr(Y_5 ~ ., data = data.frame(Y_5 = Y_5, X_5),
                ncomp = 3, validation = "CV", scale = TRUE)

## Extract scaled coefficients (these operate on the original
## unscaled X, so predictions = X %*% coefs + intercept)
coefs_5_raw <- coef(plsda_5, ncomp = 3, intercept = TRUE)
intercept   <- coefs_5_raw[1, , ]
slopes      <- coefs_5_raw[-1, , ]

## Print the derived prediction formula
cat("\n5-protein prediction formula (model-derived):\n")
cat(sprintf("Prediction = %.4f*CD11c %+.4f*CD163 %+.4f*CD44 %+.4f*CD66B %+.4f*PTEN %+.4f\n",
            slopes["CD11c"], slopes["CD163"], slopes["CD44"],
            slopes["CD66B"], slopes["PTEN"], intercept))

## Compute predictions from the model
r2$Predictions_model <- as.numeric(predict(plsda_5, newdata = X_5, ncomp = 3))

## Build the formula label from the model coefficients
formula_lbl <- sprintf(
  "Prediction = %.4f*CD11c %+.4f*CD163\n%+.4f*CD44 %+.4f*CD66B %+.4f*PTEN %+.4f",
  slopes["CD11c"], slopes["CD163"], slopes["CD44"],
  slopes["CD66B"], slopes["PTEN"], intercept)

## ====================================================================
## Panels B-E: use MODEL-derived predictions
## ====================================================================

## Use model predictions for all downstream panels
roi_data <- r2
roi_data$Predictions <- roi_data$Predictions_model

## ---- Panel B: ROI-level box plot -----------------------------------
p_B <- ggplot(roi_data,
              aes(Molecular_Response, Predictions,
                  fill = Molecular_Response)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85) +
  geom_jitter(width = 0.18, color = "black", size = 1.4, alpha = 0.6) +
  scale_fill_manual(values = c("Nonresponder" = "lightblue",
                               "Responder"    = "lightpink")) +
  labs(title = ">25% CD45+ ROIs",
       subtitle = formula_lbl,
       x = "Molecular Response", y = "Predictions") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(size = 9, hjust = 0),
        plot.subtitle = element_text(size = 7, hjust = 0),
        legend.position = "none")

## ---- Panel C: ROI-level ROC ----------------------------------------
roc_roi <- roc(roi_data$pls_code, roi_data$Predictions, quiet = TRUE)
cat(sprintf("\nROI-level AUC (model-derived): %.3f\n", as.numeric(auc(roc_roi))))
p_C <- roc_panel(roc_roi, "ROC Curve (ROI level)")

## ---- Panel D: tissue-level box plot --------------------------------
tissue_data <- roi_data %>%
  group_by(`Sample ID`, Molecular_Response, pls_code) %>%
  summarise(mean_pred = mean(Predictions, na.rm = TRUE), .groups = "drop")

mw <- wilcox.test(mean_pred ~ Molecular_Response, data = tissue_data)
roc_tis <- roc(tissue_data$pls_code, tissue_data$mean_pred, quiet = TRUE)
auc_tis <- as.numeric(auc(roc_tis))

cat(sprintf("Tissue-level: U = %.1f, p = %.3f, AUC = %.3f\n",
            mw$statistic, mw$p.value, auc_tis))

n_nr_tis <- sum(tissue_data$Molecular_Response == "Nonresponder")
n_r_tis  <- sum(tissue_data$Molecular_Response == "Responder")
ymax_tis <- max(tissue_data$mean_pred) * 1.07

p_D <- ggplot(tissue_data,
              aes(Molecular_Response, mean_pred,
                  fill = Molecular_Response)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85) +
  geom_jitter(width = 0.12, color = "black", size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("Nonresponder" = "lightblue",
                               "Responder"    = "lightpink")) +
  annotate("segment", x = 1, xend = 2, y = ymax_tis, yend = ymax_tis) +
  annotate("text", x = 1.5, y = ymax_tis * 1.06,
           label = sprintf("p = %.3f", mw$p.value), size = 3.5) +
  scale_x_discrete(labels = c(sprintf("Nonresponder\n(n = %d)", n_nr_tis),
                               sprintf("Responder\n(n = %d)", n_r_tis))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(title = "Tissue-level prediction",
       x = NULL, y = "Mean Prediction") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = "none")

## ---- Panel E: tissue-level ROC -------------------------------------
p_E <- roc_panel(roc_tis, "ROC Curve (tissue level)")

## ---- Compose and write ---------------------------------------------
top    <- p_A
middle <- p_B | p_C
bottom <- p_D | p_E
fig2 <- (top / middle / bottom) +
  plot_layout(heights = c(1.1, 1, 1)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E"))) &
  theme(plot.tag = element_text(face = "bold", size = 13))
pdf(OUT_PDF, width = 6.7, height = 8, family = "Helvetica")
print(fig2)
dev.off()
cat("\nWrote ", OUT_PDF, "\n", sep = "")
ggsave(OUT_JPG, fig2, width = 6.7, height = 8, dpi = 300)
cat("Wrote ", OUT_JPG, "\n", sep = "")
