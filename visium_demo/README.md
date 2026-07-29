# visium_demo: platform generalization to 10x Visium

Supporting files for **Figure S2**, which shows that the micro-tissue analysis
workflow is not specific to GeoMx DSP. The same PLS-R and VIP code used in
worked example 1 is applied to a public 10x Genomics Visium slide.

This folder holds the Python preprocessing and its derived inputs. The
modelling script lives in `R/Figure S2.R`, alongside the other analysis
scripts.

| Step | File | What it does |
|---|---|---|
| Format | `4_format_breast_regression.py` | Bins Visium spots into pseudo-ROIs and writes the ROI-by-gene matrix |
| Model | `../R/Figure S2.R` | PLS-R with Y = CD8A, mirroring worked example 1; writes `figures/Figure S2.pdf` and `.jpg` |

## What is committed here

| File | Purpose |
|---|---|
| `4_format_breast_regression.py` | the preprocessing step, the only new code |
| `breast_roi_immune_matrix.csv` | the 205 pseudo-ROI by 1,009 gene matrix that `R/Figure S2.R` fits |
| `spot_to_roi.csv` | per-spot pseudo-ROI assignment, used to draw panels A and B |
| `scalefactors_json.json` | 10x scale factors, gives `tissue_hires_scalef` |
| `tissue_hires_image.png` | the tissue image behind panels A and B |

These four derived files are committed so that **`R/Figure S2.R` runs from a
fresh clone with no download**. The raw slide itself is not committed; the
script downloads it on first run and `data/` is gitignored.

Bins containing fewer than eight spots are dropped, which is why some grid
squares in panel A are unfilled.

The tissue image and scale factors are small derived files from the 10x
Genomics *V1_Breast_Cancer_Block_A_Section_1* dataset, released under CC BY
4.0 and credited in the Figure S2 legend.

## Re-running the preprocessing (optional)

Only needed if you want to regenerate the matrix rather than use the committed
one. The script downloads the slide into `data/` on first run.

```bash
python 4_format_breast_regression.py   # -> breast_roi_immune_matrix.csv
```

Python needs `scanpy`, `squidpy`, `anndata`, `pandas`, `numpy` and
`scikit-misc` on Python 3.10 or newer.

## Then run the model step from the repository root

```bash
Rscript "R/Figure S2.R"
```

It writes `figures/Figure S2.pdf` and `figures/Figure S2.jpg`.
