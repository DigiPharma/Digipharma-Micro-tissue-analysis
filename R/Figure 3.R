######################################################################
## Figure 3 - Micro-tissue analysis surfaces hidden druggable targets
## Panels:
##   (A) Tissue-level mean Ki67 expression in immune-poor ROIs,
##       mNR (n = 8) vs mR (n = 6), two-sided Mann-Whitney U test.
##       Format follows Lu et al., Nat Commun 2021, Fig. 4a.
##   (B) VIP score vs PLS coefficient from PLS-R with Ki67 as output
##       on the immune-poor ROIs. ncomp = 3 (matches Nat Commun
##       Fig. 4c). B7-H3 and IDO-1 are highlighted in red.
##   (C) Scatter of B7-H3 vs IDO-1 across all immune-poor ROIs,
##       colored by Ki67 (viridis), shape by molecular response.
##
## Data: Raw data.xlsx, sheet "GBM" filtered to Ratio of CD45+ < 0.25.
##
######################################################################

suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(dplyr)
  library(pls)
  library(ggplot2)
  library(ggrepel)
  library(viridis)
  library(patchwork)
})

## ---- Standard PLS VIP helper ---------------------------------------
## VIP_j = sqrt(p * sum_a (SSY_a * w_aj^2 / sum_k w_ak^2) / sum_a SSY_a)
## (Wold et al., textbook PLS VIP). mean(VIP^2) = 1 so VIP > 1 is the
## conventional "above-average importance" cutoff.
##
## The SAME function definition appears in Figure 1.R, Figure 2.R,
## Figure 3.R, and Figure 4.R - one VIP computation across the whole
## protocol. If you ever change it, change it in all four files.
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
DATA     <- here::here("Nat commun data", "Raw data.xlsx")
OUT_PDF  <- here::here("figures", "Figure 3.pdf")
OUT_JPG  <- here::here("figures", "Figure 3.jpg")

## ---- Load and restrict to immune-poor (<25% CD45+) -----------------
gbm <- read_excel(DATA, sheet = "GBM")
ip  <- gbm %>% filter(`Ratio of CD45+ cells` < 0.25)
ip$Molecular_Response <- factor(ip$Molecular_Response,
                                levels = c("Nonresponder", "Responder"))

cat(sprintf("Immune-poor ROIs: n = %d across %d tissues\n",
            nrow(ip), length(unique(ip$`Sample ID`))))

## ---- Panel A: tissue-level Ki67 box plot ---------------------------
## IMPORTANT: the Ki67 values stored in Raw data.xlsx are
## ln(count + 1)-transformed (per Lu et al., Nat Commun 2021 methods).
## The paper's Fig. 4a Mann-Whitney test is computed on the LINEAR
## count density, not on the log values, so we back-transform each
## ROI before averaging to tissue level. This reproduces the
## published statistic (U = 7, two-sided p = 0.029, n = 8 mNR + 6 mR).
tissue <- ip %>%
  mutate(Ki67_linear = exp(Ki67) - 1) %>%
  group_by(`Sample ID`, `Patient ID`, Molecular_Response) %>%
  summarise(mean_Ki67 = mean(Ki67_linear, na.rm = TRUE),
            n_ROIs    = n(),
            .groups   = "drop")

mw <- wilcox.test(mean_Ki67 ~ Molecular_Response,
                  data = tissue)
n_nr <- sum(tissue$Molecular_Response == "Nonresponder")
n_r  <- sum(tissue$Molecular_Response == "Responder")
cat(sprintf("Tissue-level: n_NR=%d, n_R=%d, U=%.1f, p=%.4f\n",
            n_nr, n_r, mw$statistic, mw$p.value))

ymax <- max(tissue$mean_Ki67) * 1.07

p_A <- ggplot(tissue,
              aes(Molecular_Response, mean_Ki67,
                  fill = Molecular_Response)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85) +
  geom_jitter(width = 0.12, color = "black", size = 2.2, alpha = 0.85) +
  scale_fill_manual(values = c("Nonresponder" = "lightblue",
                               "Responder"    = "lightpink")) +
  annotate("segment", x = 1, xend = 2, y = ymax,        yend = ymax) +
  annotate("text", x = 1.5, y = ymax * 1.03,
           label = sprintf("p = %.3f", mw$p.value), size = 3.6) +
  scale_x_discrete(labels = c(sprintf("Nonresponder\n(n = %d)", n_nr),
                              sprintf("Responder\n(n = %d)",   n_r))) +
  labs(title = "Ki67 (immune-poor, tissue-level)",
       x = NULL, y = "Mean Ki67 count density") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = "none")

## ---- Panel B: PLS-R Y=Ki67 on immune-poor proteins -----------------
exclude_cols <- c("Patient ID", "Sample ID", "Molecular_Response",
                  "Average No. of Neighbors", "Ratio of CD45+ cells",
                  "Ki67")
X <- as.data.frame(ip[, !(names(ip) %in% exclude_cols)])
Y <- ip$Ki67

ncomp_use <- 3
pls_ki67  <- plsr(Y ~ ., data = data.frame(Y = Y, X),
                  ncomp = ncomp_use, validation = "CV", scale = TRUE)

## VIP via the shared `pls_vip` helper at the top of this file
## (identical helper in Figures 1-4) at ncomp_use = 3.
vip_vals <- pls_vip(pls_ki67, ncomp_use)
coefs    <- coef(pls_ki67, ncomp = ncomp_use)

dat_B <- data.frame(
  Protein    = colnames(X),
  Importance = vip_vals,
  Coef       = as.numeric(coefs)
)

## Highlight B7-H3 and IDO-1 (the two checkpoints called out in the
## Nat Commun paper as positive correlates of Ki67).
spotlight   <- c("B7-H3", "IDO-1")
dat_B$class <- ifelse(dat_B$Protein %in% spotlight, "spotlight",
                ifelse(dat_B$Importance > 1.0,      "important",
                                                    "other"))

p_B <- ggplot(dat_B, aes(Coef, Importance)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_point(data = subset(dat_B, class == "other"),
             color = "grey70", size = 1.6) +
  geom_point(data = subset(dat_B, class == "important"),
             color = "grey25", size = 1.8) +
  geom_point(data = subset(dat_B, class == "spotlight"),
             color = "red",   size = 2.6) +
  geom_text_repel(data = subset(dat_B, class != "other"),
                  aes(label = Protein,
                      color = class),
                  size = 2.5, box.padding = 0.25,
                  segment.color = "grey40", segment.size = 0.3,
                  max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = c("spotlight" = "red",
                                "important" = "grey25")) +
  coord_cartesian(ylim = c(0.8, NA)) +
  labs(title = "VIP vs. Coefficient (Ki67 PLS-R)",
       x = "Coefficient", y = "VIP Score") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

## ---- Panel C: B7-H3 vs IDO-1 scatter, viridis(Ki67) ----------------
p_C <- ggplot(ip, aes(`IDO-1`, `B7-H3`,
                      color = Ki67,
                      shape = Molecular_Response)) +
  geom_point(size = 2.4, alpha = 0.9) +
  scale_color_viridis(option = "viridis", name = "Ki67") +
  scale_shape_manual(values = c("Nonresponder" = 17,
                                "Responder"    = 16),
                     name = "Molecular Response") +
  labs(title = "B7-H3 and IDO-1 co-vary with Ki67 in nonresponders",
       x = "IDO-1 expression",
       y = "B7-H3 expression") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = c(0.88, 0.32),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.key.size = unit(8, "pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8))

## ---- Compose and write outputs -------------------------------------
top3    <- (p_A | p_B) + plot_layout(widths = c(1, 1.5))
bottom3 <- p_C
fig3 <- (top3 / bottom3) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = list(c("A", "B", "C"))) &
  theme(plot.tag = element_text(face = "bold", size = 13))
pdf(OUT_PDF, width = 6.7, height = 8, family = "Helvetica")
print(fig3)
dev.off()

ggsave(OUT_JPG, fig3, width = 6.7, height = 8, dpi = 300)

cat("Wrote ", OUT_PDF, "\n", sep = "")
cat("Wrote ", OUT_JPG, "\n", sep = "")
