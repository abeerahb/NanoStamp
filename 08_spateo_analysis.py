"""
08_spateo_analysis.py

Purpose
-------
Global Spateo morpho alignment of longitudinal spatial lipidomics data.

Strategy:
    - Align Day0 and Day5 whole tissue morphology
    - Use lipid PCA representation
    - Retain anatomical region labels for interpretation

Outputs
-------
alignment_overlay_global.png
transport_lines.png
morphofield_quiver.png
velocity_by_region.png
morphofield_streamlines.png
Day0_morpho_aligned.h5ad
Day5_morpho_aligned.h5ad
aligned_combined.h5ad
"""

# Load modules

import sys
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import spateo as st

from config import OUTPUTS

# ========================
# ARGUMENT
# ========================

if len(sys.argv) != 2:
    raise ValueError("Usage: python 08_spateo_analysis.py B1_T0")

dataset_id = sys.argv[1]

# ========================
# LOAD DATA
# ========================

adata_d0 = sc.read_h5ad(
    OUTPUTS /
    "anndata_raw" /
    dataset_id /
    f"{dataset_id}_Day0.h5ad"
)

adata_d5 = sc.read_h5ad(
    OUTPUTS /
    "anndata_raw" /
    dataset_id /
    f"{dataset_id}_Day5.h5ad"
)

adata_d0.obs_names_make_unique()
adata_d5.obs_names_make_unique()

adata_d0.var_names_make_unique()
adata_d5.var_names_make_unique()

print("\nDay0 regions:")
print(adata_d0.obs["region"].value_counts())

print("\nDay5 regions:")
print(adata_d5.obs["region"].value_counts())

# ========================
# CHECK REGION LABELS
# ========================

if "region" not in adata_d0.obs.columns:
    raise ValueError("adata_d0.obs['region] missing")

if "region" not in adata_d5.obs.columns:
    raise ValueError("adata_d5.obs['region] missing")

regions = [1, 2, 3]

# ========================
# MATCH COMMON FEATURES
# ========================

common = adata_d0.var_names.intersection(adata_d5.var_names)

adata_d0 = adata_d0[:, common].copy()
adata_d5 = adata_d5[:, common].copy()

# Remove zero-variance lipids

var_mask = (
    np.var(adata_d0.X, axis=0) > 1e-8
) & (
    np.var(adata_d5.X, axis=0) > 1e-8
)

adata_d0 = adata_d0[:, var_mask].copy()
adata_d5 = adata_d5[:, var_mask].copy()

print(f"Common lipid features: {adata_d0.n_vars}")

# ========================
# PCA
# ========================

st.align.group_pca(
    adatas=[adata_d0, adata_d5],
    pca_key="X_pca",
    use_hvg=False
)

# ========================
# OUTPUT DIRECTORY
# ========================

output_dir = OUTPUTS / "spateo_viewer" / dataset_id
output_dir.mkdir(parents=True, exist_ok=True)

# ========================
# GLOBAL MORPHO ALIGNMENT
# FIXED COORDINATE MODE
# ========================

print("\nRunning Spateo alignment with fixed spatial coordinates")

# Preserve original coordinates
adata_d0.obsm["spatial_original"] = (
    adata_d0.obsm["spatial"].copy()
)

adata_d5.obsm["spatial_original"] = (
    adata_d5.obsm["spatial"].copy()
)


aligned_models, pis = st.align.morpho_align(
    models=[
        adata_d0,
        adata_d5
    ],
    rep_layer="X_pca",
    rep_field="obsm",
    spatial_key="spatial",
    key_added="align_spatial",
    mode="SN-S",
    dissimilarity="euclidean",
    max_iter=100,
    device="cpu",
    verbose=True
)


adata_d0_aligned = aligned_models[0]
adata_d5_aligned = aligned_models[1]


# ==================================================
# REMOVE SPATEO ROTATION / FLIPPING EFFECTS
# Keep original tissue coordinates
# ==================================================

adata_d0_aligned.obsm["spatial"] = (
    adata_d0.obsm["spatial"].copy()
)

adata_d5_aligned.obsm["spatial"] = (
    adata_d5.obsm["spatial"].copy()
)

print("Using original spatial coordinates")

# ==================================================
# FORCE ORIGINAL SPATIAL POSITIONS
# Ignore Spateo coordinate transformations
# ==================================================

adata_d0_aligned.obsm["spatial"] = (
    adata_d0.obsm["spatial_original"].copy()
)

adata_d5_aligned.obsm["spatial"] = (
    adata_d5.obsm["spatial_original"].copy()
)

print("Spatial coordinates restored - no Spateo rotation/flipping")

# ========================
# OPTIMAL TRANSPORT
# ========================

coords0 = adata_d0_aligned.obsm["spatial"]
coords5 = adata_d5_aligned.obsm["spatial"]

pi = np.asarray(pis[0])

# normalise rows
row_sum = pi.sum(axis=1, keepdims=True)
row_sum[row_sum == 0] = 1
pi = pi / row_sum

# dimensions can differ slightly because of Spateo sampling
# Ensure OT matrix and coordinates match

n0 = min(coords0.shape[0], pi.shape[0])
n5 = min(coords5.shape[0], pi.shape[1])

coords0_field = coords0[:n0]
coords5_field = coords5[:n5]

pi = pi[:n0, :n5]

future = pi @ coords5_field

vectors = future - coords0_field

velocity = np.linalg.norm(
    vectors,
    axis=1
)

region_vector_labels = (
    adata_d0_aligned.obs["region"]
    .to_numpy()[:n0]
)

velocity = np.linalg.norm(vectors, axis=1)

# ========================
# GLOBAL ALIGNMENT PLOT
# ========================

region_names = {
    1: "Tumour",
    2: "Grey Matter",
    3: "White Matter"
}

region_colors = {
    1: "red",
    2: "royalblue",
    3: "forestgreen"
}

plt.figure(figsize=(8, 8))

for region in [1, 2, 3]:

    d0 = adata_d0_aligned.obs["region"] == region
    d5 = adata_d5_aligned.obs["region"] == region

    plt.scatter(
        adata_d0_aligned.obsm["spatial"][d0, 0],
        adata_d0_aligned.obsm["spatial"][d0, 1],
        s=5,
        color=region_colors[region],
        alpha=0.35,
        label=f"{region_names[region]} Day0"
    )

    plt.scatter(
        adata_d5_aligned.obsm["spatial"][d5,0],
        adata_d5_aligned.obsm["spatial"][d5,1],
        s=5,
        marker="x",
        color=region_colors[region],
        alpha=0.9,
        label=f"{region_names[region]} Day5"
    )

plt.gca().invert_yaxis()
plt.axis("equal")
plt.title(f"{dataset_id}: Global Morpho Alignment")
plt.xlabel("X coordinate")
plt.ylabel("Y coordinate")
plt.legend(markerscale=3, fontsize=8, ncol=2)
plt.tight_layout()

plt.savefig(
    output_dir / "alignment_overlay_global.png",
    dpi=300
)

plt.close()

plt.figure(figsize=(8,8))

best_match = np.argmax(pi, axis=1)

for region in regions:

    idx = np.where(
        region_vector_labels == region
    )[0]

    step = max(1, len(idx)//200)

    for i in idx[::step]:

        j = best_match[i]

        plt.plot(
            [
                coords0_field[i,0],
                coords5_field[j,0]
            ],
            [
                coords0_field[i,1],
                coords5_field[j,1]
            ],
            color=region_colors[region],
            alpha=0.2,
            linewidth=0.5
        )


plt.gca().invert_yaxis()
plt.axis("equal")

plt.title(
    f"{dataset_id}: Optimal Transport Correspondence"
)

plt.xlabel("X coordinate")
plt.ylabel("Y coordinate")

plt.legend(
    handles=[
        plt.Line2D(
            [0],[0],
            color=region_colors[r],
            label=region_names[r]
        )
        for r in regions
    ]
)

plt.tight_layout()

plt.savefig(
    output_dir/"transport_lines.png",
    dpi=300
)

plt.close()

plt.figure(figsize=(8,8))

for region in regions:

    idx = np.where(
        region_vector_labels == region
    )[0]

    step = max(1, len(idx)//100)

    plt.quiver(
        coords0_field[idx[::step],0],
        coords0_field[idx[::step],1],

        vectors[idx[::step],0],
        vectors[idx[::step],1],

        color=region_colors[region],

        angles="xy",
        scale_units="xy",
        scale=1,
        width=0.002,
        alpha=0.8
    )


plt.gca().invert_yaxis()
plt.axis("equal")

plt.title(
    f"{dataset_id}: Morphometric Displacement Field"
)

plt.xlabel("X coordinate")
plt.ylabel("Y coordinate")

plt.tight_layout()

plt.savefig(
    output_dir/"morphofield_quiver.png",
    dpi=300
)

plt.close()

velocity_df = adata_d0_aligned.obs.iloc[:len(velocity)].copy()

velocity_df["velocity"] = velocity

summary = (
    velocity_df
    .groupby("region")["velocity"]
    .mean()
)

labels = [
    region_names[r]
    for r in summary.index
]

colors = [
    region_colors[r]
    for r in summary.index
]

plt.figure(figsize=(5,4))

plt.bar(
    labels,
    summary.values,
    color=colors
)

plt.ylabel("Mean displacement magnitude")
plt.xlabel("Region")
plt.title(f"{dataset_id}: Morphometric Velocity by Region")

plt.tight_layout()

plt.savefig(
    output_dir/"velocity_by_region.png",
    dpi=300
)

plt.close()

# ========================
# STREAMLINES
# ========================

adata_d0_aligned.obsm["spatial_aligned"] = (
    adata_d0_aligned.obsm["spatial"]
)

adata_d5_aligned.obsm["spatial_aligned"] = (
    adata_d5_aligned.obsm["spatial"]
)

try:

    field = st.tdr.morphofield(
        source=adata_d0_aligned,
        target=adata_d5_aligned,
        spatial_key="spatial_aligned",
        vecfld_key_added="VecFld_morpho"
    )

    st.tdr.construct_field(
        field,
        vf_key="VecFld_morpho",
        n_sampling=200,
        grid_num=[60, 60]
    )

    plt.savefig(
        output_dir / "morphofield_streamlines.png",
        dpi=300
    )

    plt.close()

except Exception as e:
    print(f"Streamline construction skipped: {e}")


# ========================
# CONCATENATE & SAVE
# ========================

print("\nSaving AnnData files...")
print(output_dir)

# Remove Spateo alignment history (cannot be written to h5ad)

for adata in [adata_d0_aligned, adata_d5_aligned]:

    if "iter_spatial" in adata.uns:
        del adata.uns["iter_spatial"]

adata_d0_aligned.write_h5ad(
    output_dir / "Day0_morpho_aligned.h5ad"
)

if "align_spatial" in adata_d5_aligned.obsm:
    del adata_d5_aligned.obsm["align_spatial"]

adata_d5_aligned.write_h5ad(
    output_dir / "Day5_morpho_aligned.h5ad"
)


adata_all = ad.concat(
    [
        adata_d0_aligned,
        adata_d5_aligned
    ],
    label="timepoint",
    keys=["Day0", "Day5"],
    index_unique="-"
)


adata_all.write_h5ad(
    output_dir / "aligned_combined.h5ad"
)

print("\nFinished Spateo analysis")