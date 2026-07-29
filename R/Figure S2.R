######################################################################
## Figure S2 - PLS-R of CD8A on 10x Visium breast-cancer pseudo-ROIs
##
## (A) Pseudo-ROI grid on the tissue section, colored by mean CD8A
## (B) One pseudo-ROI magnified, with its Visium spots marked
## (C) MSEP curve vs. number of components (cross-validated)
## (D) Observed vs. predicted CD8A
## (E) VIP score vs. PLS coefficient, cytotoxicity module highlighted
##
## Same PLS-R workflow as R/Figure 1.R, applied to Visium pseudo-ROIs.
##
## Data: visium_demo/breast_roi_immune_matrix.csv
##       (205 pseudo-ROIs x 1,009 genes, written by
##        visium_demo/4_format_breast_regression.py)
## Y: CD8A. X: all genes except CD8A, CD8B, CD3D and CD3E.
######################################################################

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(pls)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(png)        # reads tissue_hires_image.png for panels A and B
  library(grid)       # rasterGrob
  library(jsonlite)   # scalefactors_json.json
})

set.seed(42)   # CV folds are random; fix for a reproducible figure

## ---- Standard PLS VIP helper (verbatim from R/Figure 1.R) ----------
pls_vip <- function(model, ncomp) {
  W <- model$loading.weights[, 1:ncomp, drop = FALSE]
  T <- model$scores[,         1:ncomp, drop = FALSE]
  Q <- model$Yloadings[,      1:ncomp, drop = FALSE]
  p <- nrow(W)
  SSY <- as.numeric(Q)^2 * colSums(T * T)
  Wn  <- sweep(W^2, 2, colSums(W^2), "/")
  as.numeric(sqrt(p * (Wn %*% SSY) / sum(SSY)))
}

## ---- Paths ---------------------------------------------------------
IN_CSV  <- here::here("visium_demo", "breast_roi_immune_matrix.csv")
OUT_PDF <- here::here("figures", "Figure S2.pdf")
OUT_JPG <- here::here("figures", "Figure S2.jpg")
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

META_COLS <- c("roi_id", "label", "n_spots", "x_um", "y_um", "epi_score")
RESPONSE  <- "CD8A"
REMOVE_FROM_X <- c("CD8A", "CD8B", "CD3D", "CD3E")  # response + co-markers
## Cytotoxicity module highlighted in panel E.
CYTO_MODULE <- c("GZMA", "GZMK", "PRF1", "NKG7", "CCL5", "PTPRC")

## The luminal epithelial differentiation programme, highlighted in panel E.
EPI_MODULE <- c("KRT8", "KRT18", "TACSTD2", "ESR1")

## Genes labelled in panel E without being highlighted. TRBC2 is named in
## the figure legend but falls outside the top-10-per-side rule.
FORCE_LABEL <- c("TRBC2")

## Deep purple separates from red by 26 greyscale levels and by a large
## simulated deuteranope distance. The epithelial module also uses a
## distinct shape, which survives greyscale and colour-blind reproduction.
EPI_COLOR <- "#4A148C"   # deep purple
## A filled triangle reads lighter than a filled circle at the same
## nominal size: at 2.4 the triangle is 0.82 the circle's area, so 2.70
## matches it.
EPI_SHAPE <- 17          # solid triangle
EPI_SIZE  <- 2.70        # vs 2.4 for the circular points; 0.99 area match

## ---- Load and assemble X, Y ----------------------------------------
dat <- read_csv(IN_CSV, show_col_types = FALSE)
Y <- dat[[RESPONSE]]
gene_cols <- setdiff(names(dat), META_COLS)
X <- as.data.frame(dat[, setdiff(gene_cols, REMOVE_FROM_X)])

## Guard: drop genes constant across ROIs (scale = TRUE divides by their
## zero SD -> NaN). This does not remove any cytotoxicity-module gene.
sds <- vapply(X, sd, numeric(1))
if (any(sds == 0)) {
  cat(sprintf("Dropping %d zero-variance genes.\n", sum(sds == 0)))
  X <- X[, sds > 0, drop = FALSE]
}
cat(sprintf("PLS-R: Y = %s, %d ROIs x %d predictor genes "
            , RESPONSE, nrow(X), ncol(X)))
cat(sprintf("(removed from X: %s)\n", paste(REMOVE_FROM_X, collapse = ", ")))
present_mod <- intersect(CYTO_MODULE, colnames(X))
cat(sprintf("Cytotoxicity-module genes present in X: %s\n",
            paste(present_mod, collapse = ", ")))

## ---- Fit PLS-R and choose ncomp from the MSEP curve ----------------
NCOMP_MAX <- 15
mod <- plsr(Y ~ ., data = data.frame(Y = Y, X), ncomp = NCOMP_MAX,
            validation = "CV", scale = TRUE)

## selectNcomp with the one-sigma rule picks the most parsimonious
## component count whose CV error is within one SE of the minimum.
## Fall back to the CV-MSEP minimum if it returns 0.
ncomp_use <- selectNcomp(mod, method = "onesigma", validation = "CV")
if (ncomp_use < 1) {
  cv_msep <- as.numeric(MSEP(mod, estimate = "CV")$val[1, 1, -1])
  ncomp_use <- which.min(cv_msep)
}
cat(sprintf("Chosen ncomp (one-sigma rule): %d\n", ncomp_use))

preds    <- as.numeric(predict(mod, ncomp = ncomp_use))
cv_preds <- as.numeric(mod$validation$pred[, 1, ncomp_use])
ok <- !is.na(cv_preds)
if (sum(!ok) > 0) cat(sprintf("WARNING: %d/%d CV predictions NaN, excluded.\n",
                              sum(!ok), length(cv_preds)))

R2X <- sum(explvar(mod)[1:ncomp_use]) / 100
R2Y <- 1 - sum((Y - preds)^2)            / sum((Y - mean(Y))^2)
Q2  <- 1 - sum((Y[ok] - cv_preds[ok])^2) / sum((Y[ok] - mean(Y[ok]))^2)
cat(sprintf("R2X = %.3f, R2Y = %.3f, Q2 = %.3f (ncomp = %d)\n",
            R2X, R2Y, Q2, ncomp_use))

## ---- VIP and reporting ---------------------------------------------
vip_dat <- data.frame(
  Gene       = colnames(X),
  Importance = pls_vip(mod, ncomp_use),
  Coef       = as.numeric(coef(mod, ncomp = ncomp_use))
)
vip_dat$is_module <- vip_dat$Gene %in% CYTO_MODULE
vip_dat$is_epi    <- vip_dat$Gene %in% EPI_MODULE
## Three drawing tiers for panel C. A gene cannot be in both modules.
vip_dat$tier <- ifelse(vip_dat$is_module, "cyto",
                ifelse(vip_dat$is_epi,    "epi", "other"))

cat("\nCytotoxicity-module genes, by VIP (Coef > 0 = tracks CD8A):\n")
print(vip_dat %>% filter(is_module) %>% arrange(desc(Importance)) %>%
        mutate(tracks_CD8A = Coef > 0), row.names = FALSE)
cat("\nEpithelial-module genes, by VIP (Coef < 0 = runs against CD8A):\n")
print(vip_dat %>% filter(is_epi) %>% arrange(desc(Importance)) %>%
        mutate(against_CD8A = Coef < 0, above_VIP1 = Importance > 1),
      row.names = FALSE)
missing_epi <- setdiff(EPI_MODULE, vip_dat$Gene)
if (length(missing_epi))
  cat(sprintf("NOTE: epithelial markers absent from the predictors: %s\n",
              paste(missing_epi, collapse = ", ")))
cat("\nTop 15 VIP genes overall:\n")
print(vip_dat %>% arrange(desc(Importance)) %>% head(15) %>%
        mutate(side = ifelse(Coef > 0, "with CD8A", "against")),
      row.names = FALSE)

## ---- Panel C: MSEP curve -------------------------------------------
msep_obj <- MSEP(mod, estimate = c("train", "CV"))
msep_df  <- data.frame(
  ncomp = rep(0:NCOMP_MAX, 2),
  MSEP  = c(as.numeric(msep_obj$val["train", , ]),
            as.numeric(msep_obj$val["CV",    , ])),
  type  = rep(c("train", "CV"), each = NCOMP_MAX + 1))

p_A <- ggplot(msep_df, aes(ncomp, MSEP, color = type, linetype = type)) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = ncomp_use, linetype = "dotted",
             color = "blue", linewidth = 0.5) +
  annotate("text", x = ncomp_use + 0.3, y = max(msep_df$MSEP) * 0.95,
           label = paste0("ncomp = ", ncomp_use), hjust = 0, size = 3,
           color = "blue") +
  scale_color_manual(breaks = c("train", "CV"),
                     values = c("train" = "black", "CV" = "red"),
                     labels = c("Training", "Cross-validation")) +
  scale_linetype_manual(breaks = c("train", "CV"),
                        values = c("train" = "solid", "CV" = "dashed"),
                        labels = c("Training", "Cross-validation")) +
  labs(title = "MSEP vs. components", x = "Number of components", y = "MSEP") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0, size = 10),
        legend.position = c(0.7, 0.85), legend.title = element_blank(),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.text = element_text(size = 8))

## ---- Panel D: observed vs predicted --------------------------------
df_pred <- data.frame(Observed = Y, Predicted = preds)
metrics_label <- sprintf("R²X=%.3f  R²Y=%.3f\nQ²=%.3f",
                         R2X, R2Y, Q2)
p_B <- ggplot(df_pred, aes(Observed, Predicted)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.6) +
  geom_point(size = 1.8, alpha = 0.6, shape = 1, stroke = 0.7) +
  annotate("text", x = Inf, y = -Inf, label = metrics_label,
           hjust = 1.1, vjust = -0.5, size = 3.2) +
  labs(title = "Observed vs. predicted CD8A",
       x = "Observed CD8A (log-norm)", y = "Predicted") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

## ---- Panel E: VIP vs coefficient, module highlighted ---------------
## Every gene in both highlighted modules is force-labelled. The
## top-N-per-side rule supplies the surrounding context genes.
n_label <- 10
top <- vip_dat %>%
  filter(Importance > 1) %>%
  group_by(side = Coef > 0) %>%
  slice_max(Importance, n = n_label) %>%
  ungroup() %>%
  bind_rows(subset(vip_dat, is_module | is_epi | Gene %in% FORCE_LABEL)) %>%
  distinct(Gene, .keep_all = TRUE)

cat(sprintf("\nPanel C: %d points (cyto %d, epi %d, other %d); %d labelled:\n",
            nrow(vip_dat), sum(vip_dat$tier == "cyto"),
            sum(vip_dat$tier == "epi"), sum(vip_dat$tier == "other"),
            nrow(top)))
cat("  ", paste(sort(top$Gene), collapse = ", "), "\n", sep = "")

## Draw order: grey background, then epithelial, then cytotoxicity on top.
p_C <- ggplot(vip_dat, aes(Coef, Importance)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0,   linetype = "dashed", color = "black") +
  geom_point(data = subset(vip_dat, tier == "other"),
             color = "grey70", size = 1.0, alpha = 0.6) +
  geom_point(data = subset(vip_dat, tier == "epi"),
             color = EPI_COLOR, shape = EPI_SHAPE, size = EPI_SIZE) +
  geom_point(data = subset(vip_dat, tier == "cyto"),
             color = "red", size = 2.4) +
  geom_text_repel(data = top,
                  aes(label = Gene, color = tier),
                  size = 2.5, box.padding = 0.3,
                  segment.color = "grey40", segment.size = 0.3,
                  max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = c(cyto = "red", epi = EPI_COLOR,
                                other = "grey20")) +
  labs(title = "VIP vs. coefficient",
       x = "Coefficient  (positive = co-varies with CD8A)",
       y = "VIP Score") +
  theme_classic(base_size = 10) +
  ## legend.position = "none": a legend box would change the panel's
  ## proportions.
  theme(plot.title = element_text(hjust = 0, size = 9),
        legend.position = "none")

## ====================================================================
## Panels A and B: the pseudo-ROI grid drawn on the tissue section
## ====================================================================
## Placement and membership are read from spot_to_roi.csv; this script
## does no coordinate or bin arithmetic.
##
## These three inputs live directly in visium_demo/ so the figure can be
## reproduced from a fresh clone. The scale factors and the tissue image
## are small derived files from the 10x Genomics public slide (CC BY
## 4.0), credited in the Figure S2 legend. The slide itself is downloaded
## on first run and stays out of the repository.
SPOT_CSV <- here::here("visium_demo", "spot_to_roi.csv")
SCALE_JSON <- here::here("visium_demo", "scalefactors_json.json")
HIRES_PNG  <- here::here("visium_demo", "tissue_hires_image.png")
PIPE_PY  <- here::here("visium_demo", "4_format_breast_regression.py")

## The pixels-per-micrometre scale is read from the pipeline rather than
## restated here, so the two cannot drift apart.
pipe_src <- readLines(PIPE_PY, warn = FALSE)
pipe_num <- function(nm) {
  ln <- grep(paste0("^", nm, "[ ]*=[ ]*[0-9.]+"), pipe_src, value = TRUE)[1]
  if (is.na(ln)) stop("could not read ", nm, " from the pipeline")
  as.numeric(sub(paste0("^", nm, "[ ]*=[ ]*([0-9.]+).*$"), "\\1", ln))
}
BIN_UM_PIPE <- pipe_num("BIN_UM")
BIN_PX      <- pipe_num("BIN_PX_COMMITTED")
PX_PER_UM   <- BIN_PX / BIN_UM_PIPE

sf    <- jsonlite::fromJSON(SCALE_JSON)
HIRES <- sf$tissue_hires_scalef
img   <- png::readPNG(HIRES_PNG)
img_h <- dim(img)[1]; img_w <- dim(img)[2]
img_grob  <- grid::rasterGrob(img, interpolate = TRUE)
BIN_HIRES <- BIN_PX * HIRES

sp <- read_csv(SPOT_CSV, show_col_types = FALSE)
sp$keep_lgl <- as.character(sp$kept) %in% c("True", "TRUE")
## column is x, row is y; the image origin is top-left, so y is flipped
## once here and every later coordinate is in plot space.
sp$x  <- sp$pxl_col_in_fullres * HIRES
sp$yp <- img_h - sp$pxl_row_in_fullres * HIRES
kept_sp <- sp[sp$keep_lgl, ]
drop_sp <- sp[!sp$keep_lgl, ]

roi_ids <- sort(unique(kept_sp$roi_id), method = "radix")

## Panel A fills each kept bin with its mean CD8A, read from the
## committed matrix.
cells <- kept_sp %>%
  distinct(roi_id, bin_i, bin_j) %>%
  left_join(dat %>% select(roi_id, CD8A), by = "roi_id") %>%
  mutate(xmin = bin_i * BIN_HIRES, xmax = (bin_i + 1) * BIN_HIRES,
         ymin = img_h - (bin_j + 1) * BIN_HIRES,
         ymax = img_h - bin_j * BIN_HIRES)

## CD8A across pseudo-ROIs is strongly right-skewed, so a linear colour
## scale leaves the mid range almost uniformly dark.
##
## A square-root transform on the fill. It is monotone, so ranking is
## preserved; no value is clipped. Panel D keeps linear axes.
CD8A_MAX <- max(cells$CD8A)
stopifnot(min(cells$CD8A) >= 0)   # the transform requires non-negative values

## Fill alpha: keeps the H and E architecture visible while the fill
## stays judgeable.
FILL_ALPHA <- 0.62
CRIMSON    <- "#DC143C"   # locator box, RGB 220 20 60
ZOOM_FILL  <- "#aec7e8"   # neutral light fill for the panel B spots

## Crop to the occupied bins plus half a bin of margin.
i_rng <- range(sp$bin_i); j_rng <- range(sp$bin_j)
xlim_A <- c(i_rng[1] - 0.5, i_rng[2] + 1.5) * BIN_HIRES
yimg_A <- c(j_rng[1] - 0.5, j_rng[2] + 1.5) * BIN_HIRES
ylim_A <- c(img_h - yimg_A[2], img_h - yimg_A[1])

grid_x <- seq(i_rng[1], i_rng[2] + 1) * BIN_HIRES
grid_y <- img_h - seq(j_rng[1], j_rng[2] + 1) * BIN_HIRES

## Zoom target, chosen deterministically: a median-sized pseudo-ROI
## nearest the centroid of all pooled spots.
cnt  <- kept_sp %>% count(roi_id, name = "n")
medn <- as.integer(median(cnt$n))
cen  <- kept_sp %>% group_by(roi_id) %>%
  summarise(n = n(), cx = mean(x), cy = mean(yp), .groups = "drop") %>%
  filter(n == medn) %>%
  mutate(d = sqrt((cx - mean(kept_sp$x))^2 + (cy - mean(kept_sp$yp))^2)) %>%
  arrange(d, n)
zoom_id <- cen$roi_id[1]
zi <- kept_sp$bin_i[kept_sp$roi_id == zoom_id][1]
zj <- kept_sp$bin_j[kept_sp$roi_id == zoom_id][1]
zx <- c(zi, zi + 1) * BIN_HIRES
zy <- img_h - c(zj + 1, zj) * BIN_HIRES
zpad   <- 0.55 * BIN_HIRES
xlim_B <- c(zx[1] - zpad, zx[2] + zpad)
ylim_B <- c(zy[1] - zpad, zy[2] + zpad)

## Scale bars, both derived from the measured scale. Sizes are stated in
## the figure legend, so the bars carry no text.
bar_A_px <- 1000 * PX_PER_UM * HIRES    # 1 mm
bar_B_px <-  100 * PX_PER_UM * HIRES    # 100 um, exactly one spot pitch
bar_rect <- function(xlim, ylim, len, frac_x = 0.06, frac_y = 0.07, thick = 0.016) {
  w <- diff(xlim); h <- diff(ylim)
  x2 <- xlim[2] - frac_x * w; x1 <- x2 - len
  y1 <- ylim[1] + frac_y * h
  data.frame(xmin = x1, xmax = x2, ymin = y1, ymax = y1 + thick * h)
}
barA <- bar_rect(xlim_A, ylim_A, bar_A_px)
barB <- bar_rect(xlim_B, ylim_B, bar_B_px)

micrograph_theme <- function() {
  list(theme_void(base_size = 10),
       theme(legend.position = "none", plot.margin = margin(1, 1, 1, 1)))
}

## Colour bar breaks, spaced for the square root scale so the low end
## does not bunch up. The maximum is included as the last break.
cb_breaks <- c(0, 0.05, 0.10, 0.20, 0.35, CD8A_MAX)
cb_labels <- format(round(cb_breaks, 2), nsmall = 2)

p_grid <- ggplot() +
  annotation_custom(img_grob, 0, img_w, 0, img_h) +
  geom_rect(data = cells,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = CD8A), alpha = FILL_ALPHA) +
  geom_segment(data = data.frame(x = grid_x),
               aes(x = x, xend = x, y = ylim_A[1], yend = ylim_A[2]),
               colour = "black", linewidth = 0.18, alpha = 0.55) +
  geom_segment(data = data.frame(y = grid_y),
               aes(x = xlim_A[1], xend = xlim_A[2], y = y, yend = y),
               colour = "black", linewidth = 0.18, alpha = 0.55) +
  ## Spots are drawn faintly; a heavier weight would compete with the
  ## H and E architecture. Panel B shows individual spots.
  geom_point(data = drop_sp, aes(x, yp), shape = 21, fill = NA,
             colour = "grey35", stroke = 0.10, size = 0.35) +
  geom_point(data = kept_sp, aes(x, yp), colour = "grey20", size = 0.12,
             alpha = 0.25) +
  ## Locator box for panel B, drawn from the same bin bounds panel B uses.
  annotate("rect", xmin = zx[1], xmax = zx[2], ymin = zy[1], ymax = zy[2],
           fill = NA, colour = CRIMSON, linewidth = 0.75) +
  geom_rect(data = barA, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "white", colour = "black", linewidth = 0.25) +
  scale_fill_viridis_c(name = "mean CD8A", transform = "sqrt",
                       breaks = cb_breaks, labels = cb_labels,
                       limits = c(0, CD8A_MAX)) +
  coord_fixed(xlim = xlim_A, ylim = ylim_A, expand = FALSE) +
  theme_void(base_size = 10) +
  ## The colour bar sits below the panel, not inside it.
  theme(plot.margin = margin(1, 1, 1, 1),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width = unit(34, "pt"),
        legend.key.height = unit(6, "pt"),
        legend.title = element_text(size = 7, vjust = 0.9),
        legend.text = element_text(size = 6),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"))

p_zoom <- ggplot() +
  annotation_custom(img_grob, 0, img_w, 0, img_h) +
  annotate("rect", xmin = zx[1], xmax = zx[2], ymin = zy[1], ymax = zy[2],
           fill = NA, colour = "black", linewidth = 0.7) +
  geom_point(data = kept_sp[kept_sp$roi_id != zoom_id, ], aes(x, yp),
             shape = 21, fill = NA, colour = "grey40", stroke = 0.3, size = 2.6) +
  geom_point(data = kept_sp[kept_sp$roi_id == zoom_id, ], aes(x, yp),
             shape = 21, fill = ZOOM_FILL,
             colour = "black", stroke = 0.3, size = 3.4) +
  geom_rect(data = barB, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "white", colour = "black", linewidth = 0.25) +
  coord_fixed(xlim = xlim_B, ylim = ylim_B, expand = FALSE) +
  micrograph_theme()

cat(sprintf("\nPanel A: %d spots drawn (%d pooled, %d below the %d-spot rule), %d pseudo-ROIs\n",
            nrow(sp), nrow(kept_sp), nrow(drop_sp),
            pipe_num("MIN_SPOTS_PER_ROI"), length(roi_ids)))
cat(sprintf("  scale read from pipeline: %.4f px/um; bin %.2f um = %.2f full-res px = %.2f hires px\n",
            PX_PER_UM, BIN_UM_PIPE, BIN_PX, BIN_HIRES))
cat(sprintf("  scale bars: A 1 mm = %.2f hires px; B 100 um = %.2f hires px\n",
            bar_A_px, bar_B_px))
n_drop_bins <- length(unique(paste(drop_sp$bin_i, drop_sp$bin_j)))
cat(sprintf("  bins filled by mean CD8A: %d; bins left unfilled by the %d-spot rule: %d\n",
            nrow(cells), pipe_num("MIN_SPOTS_PER_ROI"), n_drop_bins))
cat(sprintf("  CD8A fill: min %.4f, median %.4f, max %.4f; square-root scale, nothing clipped; alpha %.2f\n",
            min(cells$CD8A), median(cells$CD8A), max(cells$CD8A), FILL_ALPHA))
top_roi <- cells$roi_id[which.max(cells$CD8A)]
cat(sprintf("  highest pseudo-ROI: %s at %.4f (unique maximum: %s; %.2fx the runner-up)\n",
            top_roi, CD8A_MAX, sum(cells$CD8A == CD8A_MAX) == 1,
            CD8A_MAX / sort(cells$CD8A, decreasing = TRUE)[2]))
cat(sprintf("  range use under sqrt: median %.1f%%, p95 %.1f%%, p98 %.1f%% of the colour range\n",
            sqrt(median(cells$CD8A)) / sqrt(CD8A_MAX) * 100,
            sqrt(quantile(cells$CD8A, .95)) / sqrt(CD8A_MAX) * 100,
            sqrt(quantile(cells$CD8A, .98)) / sqrt(CD8A_MAX) * 100))
## Panel A and panel D plot the same quantity, so the brightest cell in A
## must be the extreme observed point in D.
cat(sprintf("  same ROI is the extreme observed point in panel D: %s\n",
            identical(top_roi, dat$roi_id[which.max(Y)])))
## The locator box is exactly one bin, so its corners must fall on grid
## intersections.
on_x <- all(sapply(zx, function(v) any(abs(grid_x - v) < 1e-6)))
on_y <- all(sapply(zy, function(v) any(abs(grid_y - v) < 1e-6)))
cat(sprintf("  locator box corners on grid intersections: x %s, y %s\n",
            ifelse(on_x, "yes", "NO"), ifelse(on_y, "yes", "NO")))
if (!(on_x && on_y)) stop("locator box is not aligned to the grid; stopping.")
cat(sprintf("Panel B: %s, %d spots\n", zoom_id, sum(kept_sp$roi_id == zoom_id)))

## ---- Compose -------------------------------------------------------
## Five panels: the two micrographs first, then model, then result.
fig <- (p_grid | p_zoom) / (p_A | p_B | p_C) +
  plot_layout(heights = c(1.55, 1))
FIG_W <- 10; FIG_H <- 9.2
fig <- fig +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E"))) &
  theme(plot.tag = element_text(face = "bold", size = 13))

pdf(OUT_PDF, width = FIG_W, height = FIG_H, family = "Helvetica")
print(fig)
invisible(dev.off())
cat("\nWrote ", OUT_PDF, "\n", sep = "")
ggsave(OUT_JPG, fig, width = FIG_W, height = FIG_H, dpi = 300)
cat(sprintf("Wrote %s  (%.1f x %.1f in at 300 dpi = %.0f x %.0f px)\n",
            OUT_JPG, FIG_W, FIG_H, FIG_W * 300, FIG_H * 300))
