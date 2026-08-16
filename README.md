# Protocol for identifying cellular and molecular signatures from spatial profiling data using micro-tissue analysis

Reproducible R code, raw demonstration data, and pre-rendered figures for the
four worked examples in the STAR Protocols article *Protocol for identifying
cellular and molecular signatures from spatial profiling data using micro-tissue
analysis*.

## Authors

Huiqian Hu,<sup>1,3,\*</sup> Alphonsus H. C. Ng,<sup>1,2,3,\*</sup> and
Yue Lu<sup>1,4,\*\*</sup>

<sup>1</sup>Department of Molecular Pharmaceutics, University of Utah, Salt Lake City, UT 84112, USA
<sup>2</sup>Department of Biomedical Engineering, University of Utah, Salt Lake City, UT 84112, USA
<sup>3</sup>Technical contact
<sup>4</sup>Lead contact

\* Correspondence: alphonsus.ng@pharm.utah.edu; huiqian.hu@utah.edu
\*\* Correspondence: yue.lu@pharm.utah.edu

- **Lead contact** — Yue Lu (`yue.lu@pharm.utah.edu`).
  Further information and requests for resources should be directed to and
  will be fulfilled by the lead contact.
- **Technical contacts** — Huiqian Hu (`huiqian.hu@utah.edu`) and
  Alphonsus H. C. Ng (`alphonsus.ng@pharm.utah.edu`). Technical questions
  on executing this protocol should be directed to either technical contact.

## Summary

Spatial tissue profiling measures hundreds to thousands of genes or proteins
across regions of a tissue section, but most analyses examine each feature one
at a time. This protocol provides a step-by-step workflow for *micro-tissue
analysis*, which treats each profiled region as its own window into the local
tissue biology. By grouping regions by microenvironment and analysing all
features together, the workflow identifies which genes, proteins, or cell
types collectively drive a biological outcome. All code and worked examples
are provided using two published GeoMx Digital Spatial Profiler datasets from
tumour and mucosal tissues.

For complete details on the use and execution of this protocol, please refer
to Lu et al.<sup>1</sup> and Ng et al.<sup>2</sup>

## What this repository contains

```
.
├── LICENSE                 # MIT for the code; demonstration data carved out
├── CITATION.cff            # how to cite this software
├── R/
│   ├── Figure 1.R   — PLS-R of CD8A in immune-rich GBM ROIs
│   ├── Figure 2.R   — PLS-DA prediction model for treatment response (5 panels A–E)
│   ├── Figure 3.R   — PLS-R of Ki67 in immune-poor ROIs (highlights B7-H3 / IDO-1)
│   ├── Figure 4.R   — Lung-vs-gut PLS-R of Macrophage M1 fraction
│   ├── Figure S2.R  — platform generalization: PLS-R of CD8A on 10x Visium pseudo-ROIs
│   └── Figure S3.R  — permutation null for the worked example 1 cross-validated Q²
├── Nat commun data/        # demonstration dataset 1 (Lu et al., 2021)
│   ├── Raw data.xlsx
│   ├── Raw data_plsda round 1.xlsx
│   └── Raw data_plsda round 2.xlsx
├── Cell Rep data/          # demonstration dataset 2 (Ng et al., 2023)
│   ├── absolute cell fraction.csv
│   └── whole cell fraction.xlsx
├── visium_demo/            # Python preprocessing for the Visium demonstration
│   ├── 4_format_breast_regression.py
│   ├── breast_roi_immune_matrix.csv    # derived input to R/Figure S2.R
│   ├── spot_to_roi.csv                 # per-spot pseudo-ROI assignment
│   ├── scalefactors_json.json          # 10x slide scale factors
│   └── tissue_hires_image.png          # 10x tissue image, panels A and B
└── figures/                # pre-rendered, editable PDFs plus 300 dpi JPGs
    ├── Figure 1.pdf   … Figure 4.pdf
    ├── Figure S2.pdf
    └── Figure S3.pdf
```

## How to reproduce the figures

1. Install **R ≥ 4.0** from CRAN (https://cran.r-project.org) and, optionally,
   **RStudio Desktop** from Posit (https://posit.co/download/rstudio-desktop).
2. Clone this repository:

   ```bash
   git clone https://github.com/DigiPharma/Digipharma-Micro-tissue-analysis.git
   ```

   then open `Digipharma-Micro-tissue-analysis.Rproj` in RStudio. (If you do
   not use RStudio, `setwd()` to the cloned repository's root before sourcing
   any script.)
3. Install the required R packages once:

   ```r
   install.packages(c(
     "here", "readxl", "readr", "tidyr", "dplyr",
     "pls", "ggplot2", "ggrepel", "patchwork",
     "pROC", "viridis",
     "jsonlite", "png"          # required by R/Figure S2.R only
   ))
   ```

4. Run any of the scripts:

   ```r
   source("R/Figure 1.R")
   source("R/Figure 2.R")
   source("R/Figure 3.R")
   source("R/Figure 4.R")
   source("R/Figure S2.R")
   source("R/Figure S3.R")
   ```

   Each script reads its inputs from `Nat commun data/`, `Cell Rep data/` or
   `visium_demo/` and writes its panels into `figures/` as a vector PDF and a
   300 dpi JPG.

   `R/Figure S2.R` reads only the derived files committed in `visium_demo/`, so
   it runs from a fresh clone without downloading the 10x slide. Re-running the
   Python preprocessing is optional and is described in `visium_demo/README.md`.

## Notes on the analysis

- **Multivariate model.** All four worked examples use partial least squares
  (PLS).<sup>5</sup> Examples 1, 3, and 4 use PLS regression (PLS-R, continuous
  outcome); example 2 uses PLS discriminant analysis (PLS-DA, categorical
  outcome).
- **VIP scores** are computed via the inline helper `pls_vip()`, which
  implements the standard formula
  `VIP_j = sqrt(p · Σ_a (SSY_a · w_aj² / Σ_k w_ak²) / Σ_a SSY_a)`.
  The same helper definition appears in all five scripts so the protocol
  applies one consistent VIP method throughout. With this convention
  `mean(VIP²) = 1`, so the conventional cutoff is `VIP > 1`.
- **Mann–Whitney p-values** (Figure 2D, Figure 3A) use R's exact small-sample
  test (default behaviour of `wilcox.test()`).
- **ROC curves** (Figure 2C, 2E) are drawn by the `roc_panel()` helper at the
  top of `R/Figure 2.R`.
- **Cross-validation folds are fixed.** Every modelling script calls
  `set.seed(42)` before `plsr()`, so R²X, R²Y, and Q² reproduce exactly from a
  fresh clone.

## Demonstration datasets

| Folder | Source | Tissue | Profiling |
|---|---|---|---|
| `Nat commun data/` | Lu et al., 2021<sup>1</sup> | GBM tumours treated with neoadjuvant anti-PD-1 (168 ROIs across 14 tumours) | GeoMx DSP, 35-protein immuno-oncology panel |
| `Cell Rep data/`   | Ng et al., 2023<sup>2</sup> | Lung and gut autopsy tissue from a COVID-19 case (53 ROIs) | GeoMx DSP, 1,833-transcript Cancer Transcriptome Atlas |

Both demonstration datasets were generated on the NanoString (Bruker) GeoMx
Digital Spatial Profiler.<sup>3</sup> The cell-fraction inputs in
`Cell Rep data/` were derived from the spatial transcript data using
CIBERSORTx.<sup>4</sup>

Please cite the original publications when reusing these data.

## Citation

An archived version of record is deposited at Zenodo:

> Hu, H., Ng, A. H. C., and Lu, Y. (2026). Micro-tissue analysis: R code for
> identifying cellular and molecular signatures from spatial profiling data.
> Zenodo. https://doi.org/10.5281/zenodo.19967582

That DOI is the **concept DOI**, which always resolves to the most recent
version of the archive. Each individual release also has its own version DOI,
listed in the Versions panel on the Zenodo record. `CITATION.cff` carries the
same metadata in machine-readable form.

## License

The code in `R/` and `visium_demo/` is released under the **MIT license**; see
`LICENSE`.

The demonstration data are not covered by it. The files in `Nat commun data/`
and `Cell Rep data/` are taken from Lu et al.<sup>1</sup> and Ng et al.<sup>2</sup>
and remain subject to the terms of those publications. The files in
`visium_demo/` derived from the 10x Genomics
*V1_Breast_Cancer_Block_A_Section_1* dataset are subject to the 10x Genomics
terms of use for that dataset. Cite the original source when reusing either.

## References

1. Lu, Y., Ng, A. H. C., Chow, F. E., Everson, R. G., Helmink, B. A.,
   Tetzlaff, M. T., Thakur, R., Wargo, J. A., Cloughesy, T. F., Prins, R. M.,
   and Heath, J. R. (2021). Resolution of tissue signatures of therapy
   response in patients with recurrent GBM treated with neoadjuvant anti-PD1.
   *Nature Communications* **12**, 4031.
   https://doi.org/10.1038/s41467-021-24293-4
2. Ng, A. H. C., Hu, H., Wang, K., Scherler, K., Warren, S. E.,
   Zollinger, D. R., McKay-Fleisch, J., Sorg, K., Beechem, J. M.,
   Ragaglia, E., et al. (2023). Organ-specific immunity: A tissue analysis
   framework for investigating local immune responses to SARS-CoV-2.
   *Cell Reports* **42**, 113212.
   https://doi.org/10.1016/j.celrep.2023.113212
3. Merritt, C. R., Ong, G. T., Church, S. E., Barker, K., Danaher, P.,
   Geiss, G., Hoang, M., Jung, J., Liang, Y., McKay-Fleisch, J., et al.
   (2020). Multiplex digital spatial profiling of proteins and RNA in fixed
   tissue. *Nature Biotechnology* **38**, 586–599.
   https://doi.org/10.1038/s41587-020-0472-9
4. Newman, A. M., Steen, C. B., Liu, C. L., Gentles, A. J.,
   Chaudhuri, A. A., Scherer, F., Khodadoust, M. S., Esfahani, M. S.,
   Luca, B. A., Steiner, D., et al. (2019). Determining cell type abundance
   and expression from bulk tissues with digital cytometry.
   *Nature Biotechnology* **37**, 773–782.
   https://doi.org/10.1038/s41587-019-0114-2
5. Wold, S., Sjöström, M., and Eriksson, L. (2001). PLS-regression: a basic
   tool of chemometrics. *Chemometrics and Intelligent Laboratory Systems*
   **58**, 109–130.
   https://doi.org/10.1016/S0169-7439(01)00155-1
