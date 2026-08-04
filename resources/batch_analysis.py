# Load packages
import sys
import numpy as np

# HPC adaptation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pyimzml.ImzMLParser import ImzMLParser
from scipy.signal import find_peaks

from config import (
    PROCESSED_DATA,
    OUTPUTS
)

# =================
# ARGUMENT
# =================

if len(sys.argv) != 2:
    raise ValueError(
        "Usage: python <script_name>.py B1_T0"
    )

dataset_id = sys.argv[1]

# ====================
# PROCESSING FUNCTIONS
# ====================

def process_segmentation(folder, sample_output_dir, sample_name):
    """Handles parsing and drawing segmentation file."""
    segmentation_path = folder / "segmentation.npy"
    if segmentation_path.exists():
        print("--> Processing Segmentation Map...")
        try:
            mask = np.load(segmentation_path)

            plt.figure(figsize=(6, 6))
            plt.imshow(mask, cmap="tab20", origin="lower")
            plt.colorbar(label="Cluster ID")
            plt.title(f"Segmentation Inspection - {sample_name}")

            out_file = sample_output_dir / f"segmentation_{sample_name.split('_')[-1]}.png"
            plt.savefig(out_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"   [Success] Saved: {out_file.name}")
        except Exception as e:
            print(f"  [ERROR] Failed during segmentation plotting: {e}")
            plt.close()
    else:
        print(f"[Warning] segmentation.npy missing in {folder}")

def process_mass_spec(folder, sample_output_dir, sample_name):
    """Loops over representation formats to build spectra and images."""
    representations = ["common_representation", "meaningful_signal"]

    for rep in representations:

        imzml_path_upper = folder / f"{rep}.imzML"
        imzml_path_lower = folder / f"{rep}.imzml"

        if imzml_path_upper.exists():
            imzml_path = imzml_path_upper
        elif imzml_path_lower.exists():
            imzml_path = imzml_path_lower
        else:
            print(f"[Warning] No imzML file found for {rep} in {sample_name}")
            continue

        print(f"--> Processing Representation: {rep} ({imzml_path.name})")

        try:
            p = ImzMLParser(str(imzml_path))
            total_pixels = len(p.coordinates)

            # Get first spectrum to calibrate m/z axis lengths
            mz_axis, _ = p.getspectrum(0)

            # Memory safe mean spectrum calculation
            sum_intensities = np.zeros_like(mz_axis, dtype=np.float64)

            for i in range(total_pixels):
                _, intensity = p.getspectrum(i)
                sum_intensities += intensity

            mean_spectrum = sum_intensities / total_pixels

            # Plot mean spectrum
            plt.figure(figsize=(12, 4))
            plt.plot(mz_axis, mean_spectrum, color="black", lw=0.5)
            plt.xlabel("m/z")
            plt.ylabel("Mean Intensity")
            plt.title(f"Mean Spectrum ({rep}) - {sample_name}")
            plt.grid(True, alpha=0.3)

            spec_out = sample_output_dir / f"mean_spectrum_{rep}.png"
            plt.savefig(spec_out, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  [Success] Saved: {spec_out.name}")

            # Automated target peak diagnostics
            peaks, props = find_peaks(mean_spectrum, height=1000)
            if len(peaks) > 0:
                top_peaks = mz_axis[peaks]
                top_heights = props["peak_heights"]
                idx = np.argsort(top_heights)[::-1]
                print(f" Top 5 detected peaks for {rep}:")
                for i in idx[:5]:
                    print(f"   m/z: {top_peaks[i]:.4f} | Intensity: {top_heights[i]:.2f}")
            else:
                print("  No peaks detected above the height threshold of 1000.")
            
            # Target ion image generation (m/z 885.5)
            target_mz = 885.5
            tol = 0.3

            coords = p.coordinates
            max_x = max(c[0] for c in coords)
            max_y = max(c[1] for c in coords)

            img = np.zeros((max_y, max_x))

            for i, (x, y, *z) in enumerate(coords):
                mz, intensity = p.getspectrum(i)
                mask_mz = (mz > target_mz - tol) & (mz < target_mz + tol)
                img[y-1, x-1] = intensity[mask_mz].sum()

            # Plot and save individual representation ion image
            plt.figure(figsize=(6, 6))
            plt.imshow(img, cmap="inferno", origin="lower")
            plt.colorbar(label="Intensity")
            plt.title(f"Ion image m/z {target_mz} ({rep})\n{sample_name}")
            plt.axis('off')

            ion_out = sample_output_dir / f"ion_885_{rep}.png"
            plt.savefig(ion_out, dpi=300, bbox_inches='tight')
            plt.close()
            print(f" [Success] Saved: {ion_out.name}")

        except Exception as ms_error:
            print(f"  [ERROR] Failed inspection for {rep}: {ms_error}")
            plt.close()

# =================
# MAIN EXECUTION
# =================

def main():
    dataset_folder = PROCESSED_DATA / dataset_id

    print(f"PROCESSING DATASET: {dataset_id}")

    # Loop through each timepoint folder
    for day in ["Day0", "Day5"]:

        folder = dataset_folder / day

        if not folder.exists() or not folder.is_dir():
            print(f"[Skipping] Folder missing or not a directory: {folder}")
            continue

        imzml_path = folder / "common_representation.imzML"
        mask_path = folder / "segmentation.npy"

        if not imzml_path.exists() or not mask_path.exists():
            print(f"[Skipping] Missing structural files in: {folder}")
            continue

        # FIX: Added "batch_analysis" into the path hierarchy
        day_output_dir = OUTPUTS / "batch_analysis" / dataset_id / day
        day_output_dir.mkdir(parents=True, exist_ok=True)

        sample_name = f"{dataset_id}_{day}"
        print(f"\nProcessing Timepoint -> {sample_name}")
        print(f"Saving outputs to    -> {day_output_dir}")

        # 1. Process segmentation file
        process_segmentation(folder, day_output_dir, sample_name)

        # 2. Process mass spectrometry files
        process_mass_spec(folder, day_output_dir, sample_name)

    print(f"\nProcessing complete! All figures saved inside: {OUTPUTS / 'batch_analysis' / dataset_id}")

if __name__ == "__main__":
    main()