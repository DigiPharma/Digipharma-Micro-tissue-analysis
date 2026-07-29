######################################################################
## Data formatting for the breast-cancer regression example (Y = CD8A).
##
## Bins the Visium spots of one slide into square pseudo-ROIs and writes
## an ROI-by-gene matrix. The gene set is the seurat_v3 highly variable
## panel plus a forced immune / cytotoxicity panel, so the model can see
## the T-cell module.
##
## Reports the per-ROI distribution of CD8A and MKI67; both columns are
## written to the output matrix.
##
## Usage:  python 4_format_breast_regression.py
######################################################################

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")

# ---- Parameters -----------------------------------------------------
SAMPLE_ID = "V1_Breast_Cancer_Block_A_Section_1"
## Grid square in micrometres. The pixels-per-micrometre scale is derived
## from the 100 um centre-to-centre spot pitch, not from
## spot_diameter_fullres, which 10x documents as a visualisation value
## that varies by slide design.
BIN_UM = 413.71
MIN_SPOTS_PER_ROI = 8
## Centre-to-centre spot pitch, a fixed Visium v1 specification.
SPOT_PITCH_UM = 100.0
## Bin width in full-resolution pixels that the committed
## breast_roi_immune_matrix.csv was built with. The scale computed below
## must reproduce this, otherwise spots would move between pseudo-ROIs and the
## analysis would silently change. Asserted below.
BIN_PX_COMMITTED = 1129.4290
N_HVG = 1000
HVG_FLAVOR = "seurat_v3"

# Immune / cytotoxicity panel, force-included so the model can see the
# T-cell module even if these genes are not top HVGs. CD8A is the
# response; CD8B/CD3D/CD3E are its direct co-markers (removed from X in
# the R step); the rest are the cytotoxicity module that co-varies with
# CD8A.
IMMUNE_PANEL = ["CD8A", "CD8B", "CD3D", "CD3E", "GZMB", "GZMK",
                "PRF1", "NKG7", "CCL5", "PTPRC"]
# Secondary response reported alongside CD8A, plus epithelial markers
# kept so this matrix can reproduce the same reference label as Figure S2.
SECOND_Y = ["MKI67"]
EPI_MARKERS = ["EPCAM", "KRT8", "KRT18", "ERBB2"]
FORCE = IMMUNE_PANEL + SECOND_Y + EPI_MARKERS

HERE = Path(__file__).parent
DATA_DIR = HERE / "data"
OUT_CSV = HERE / "breast_roi_immune_matrix.csv"


def load_slide():
    local = DATA_DIR / SAMPLE_ID
    if (local / "filtered_feature_bc_matrix.h5").exists():
        print(f"Loading local slide from {local}")
        return sc.read_visium(local)
    print(f"Downloading {SAMPLE_ID} ...")
    return sc.datasets.visium_sge(sample_id=SAMPLE_ID)


def dense(m):
    return np.asarray(m.todense()) if hasattr(m, "todense") else np.asarray(m)


def describe(series, name):
    q = series.quantile([0, .1, .25, .5, .75, .9, 1.0]).round(3)
    print(f"  {name:6s} mean={series.mean():.3f} sd={series.std():.3f} "
          f"CV={series.std()/max(series.mean(),1e-9):.2f} "
          f"frac_detected(>0)={ (series>0).mean():.2f}")
    print(f"         quantiles 0/10/25/50/75/90/100: {list(q.values)}")


def main():
    adata = load_slide()
    adata.var_names_make_unique()
    adata = adata[adata.obs["in_tissue"] == 1].copy()
    sc.pp.filter_genes(adata, min_cells=10)

    missing = [g for g in FORCE if g not in adata.var_names]
    if missing:
        print(f"NOTE: forced genes absent from slide (dropped): {missing}")

    # seurat_v3 HVGs on raw counts, before normalizing.
    counts = adata.copy()
    sc.pp.highly_variable_genes(counts, n_top_genes=N_HVG, flavor=HVG_FLAVOR)
    hvg = counts.var_names[counts.var["highly_variable"]].tolist()
    del counts

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    force_present = [g for g in FORCE if g in adata.var_names]
    forced_in = [g for g in IMMUNE_PANEL + SECOND_Y if g not in hvg and g in adata.var_names]
    print(f"Forced-in genes not among the top HVGs: {forced_in}")
    genes = sorted(set(hvg) | set(force_present))
    print(f"Genes kept: {len(genes)}")

    # ---- Grid aggregation ------------------------------------------
    ## Pixels per micrometre, derived from the 100 um centre-to-centre spot
    ## pitch. The pitch is measured from the data as the median
    ## nearest-neighbour distance between in-tissue spot centres, so this
    ## stays correct on any slide.
    xy = adata.obsm["spatial"].astype(float)   # column 0 is x, column 1 is y
    nn = cKDTree(xy).query(xy, k=2)[0][:, 1]
    px_per_um = float(np.median(nn)) / SPOT_PITCH_UM
    bin_px = BIN_UM * px_per_um
    print(f"Scale: median spot pitch {np.median(nn):.1f} px = {SPOT_PITCH_UM:.0f} um "
          f"-> {px_per_um:.4f} px/um; bin {BIN_UM:.2f} um = {bin_px:.4f} px")

    ## The committed matrix is only valid if this bin width groups spots
    ## into exactly the same pseudo-ROIs. A loose tolerance is not enough: at BIN_UM = 413.6 the bin
    ## width is only 2.7e-04 off, which would pass a 1e-3 check, yet it
    ## reassigns 34 spots and yields 206 pseudo-ROIs instead of 205. The
    ## tolerance below is tight enough that only a bin width preserving
    ## every spot's assignment can pass.
    rel = abs(bin_px - BIN_PX_COMMITTED) / BIN_PX_COMMITTED
    if rel > 1e-5:
        raise SystemExit(
            f"ABORT: bin width moved. Got {bin_px:.4f} px, expected "
            f"{BIN_PX_COMMITTED:.4f} px (relative difference {rel:.2e} > 1e-5). "
            "Spots would move between pseudo-ROIs and the committed matrix "
            "would no longer describe this grid. Nothing was written.")
    print(f"  bin width check passed: {bin_px:.4f} px vs committed "
          f"{BIN_PX_COMMITTED:.4f} px (relative difference {rel:.1e})")
    grid = np.floor(xy / bin_px).astype(int)
    roi_key = np.array([f"{i}_{j}" for i, j in grid])

    expr = pd.DataFrame(dense(adata[:, genes].X), columns=genes)
    expr["roi_key"] = roi_key
    expr["x_um"] = xy[:, 0] / px_per_um
    expr["y_um"] = xy[:, 1] / px_per_um
    n_spots = expr.groupby("roi_key").size().rename("n_spots")
    roi = expr.groupby("roi_key")[genes + ["x_um", "y_um"]].mean().join(n_spots)
    roi = roi[roi["n_spots"] >= MIN_SPOTS_PER_ROI].copy()

    # Reference epithelial label, same rule as Figure S2 (for context).
    roi["epi_score"] = roi[EPI_MARKERS].mean(axis=1)
    roi["label"] = np.where(roi["epi_score"] >= roi["epi_score"].median(),
                            "epi_high", "epi_low")
    roi = roi.reset_index().rename(columns={"roi_key": "roi_id"})
    roi["roi_id"] = "ROI_" + roi["roi_id"]
    print(f"\nPseudo-ROIs: {len(roi)} (mean {roi['n_spots'].mean():.1f} spots/ROI)")

    # ---- Response distribution: CD8A and MKI67 ----------------------
    print("\n=== Response distribution (per-ROI mean, log-normalized) ===")
    describe(roi["CD8A"], "CD8A")
    describe(roi["MKI67"], "MKI67")
    # A response is usable if it is detected in a reasonable fraction of
    # ROIs and has real spread.
    cd8a_ok = (roi["CD8A"] > 0).mean() >= 0.5 and roi["CD8A"].std() > 0.05
    verdict = "USABLE" if cd8a_ok else "LOW VARIANCE"
    print(f"\nCD8A verdict: {verdict}")

    meta_cols = ["roi_id", "label", "n_spots", "x_um", "y_um",
                 "epi_score"]
    out = roi[meta_cols + genes]
    out.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV}  ({out.shape[0]} ROIs x {len(genes)} genes)")


if __name__ == "__main__":
    main()
