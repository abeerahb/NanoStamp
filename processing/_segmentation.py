"""Mass spectrometry image segmentation

This module contains:
    * SegmentationInterface
    * MeanSegmentation

"""

import numpy as np
from typing import List
from abc import ABC, abstractmethod
from skimage import filters
from skimage.morphology import disk


class SegmentationInterface(ABC):
    """Interface for MSI segmentation."""

    @abstractmethod
    def segment(self, img: np.ndarray) -> np.ndarray:
        raise NotImplementedError



class MeanSegmentation(SegmentationInterface):
    """Segmentation using representative lipid peaks."""

    def __init__(
        self,
        mzs: np.ndarray,
        representative_peaks: List[float],
        mass_resolution: float,
        threshold_percentile: float = 50
    ) -> None:

        self.mzs = mzs
        self.peaks = representative_peaks
        self.res = mass_resolution
        self.threshold_percentile = threshold_percentile



    def segment(self, img: np.ndarray) -> np.ndarray:


        peak_images = []


        for peak in self.peaks:

            idx = (
                (self.mzs >= peak - self.res) &
                (self.mzs <= peak + self.res)
            )


            if not np.any(idx):
                raise ValueError(
                    f"No m/z values found around peak {peak}"
                )


            peak_img = img[:, :, idx].sum(axis=-1)


            # z-score lipid image
            peak_img = (
                peak_img - peak_img.mean()
            ) / peak_img.std()


            peak_images.append(
                peak_img
            )


        # combine lipid signal
        peaks_img = np.mean(
            peak_images,
            axis=0
        )


        # smooth
        smooth = filters.median(
            peaks_img,
            disk(2)
        )


        # 50% threshold
        threshold = np.percentile(
            smooth,
            self.threshold_percentile
        )


        segmentation = smooth > threshold


        print(
            "Segmentation threshold:",
            threshold
        )

        print(
            "Tumour pixels:",
            segmentation.sum(),
            "/",
            segmentation.size
        )


        return segmentation