"""Image overlays and simple per-channel pixel SNR."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tifffile
from matplotlib.patches import Patch
from skimage.filters import threshold_otsu

from .base import QCObject


class ImageQC(QCObject):
    def run(self, root, panel, images, cells, channels, report):
        names = {str(n).strip(): str(n).strip() for n in panel.name}
        names.update({str(c).strip(): str(n).strip() for c, n in zip(panel.channel, panel.name)})
        if "CD3E" in names:
            names.setdefault("CD3", names["CD3E"])
        channels = [names[x.strip()] for x in channels.split(",") if x.strip()] if channels else panel.name.iloc[:3].tolist()
        if not channels or len(channels) > 3:
            raise ValueError("Choose one to three --image-channels")
        rows = []
        for image in images.itertuples(index=False):
            filename = Path(str(image.image)).name
            img_path = root / self.args.images / filename
            mask_path = root / self.args.masks / filename
            if not mask_path.exists() and self.args.masks == "masks" and (root / "masks_deepcell" / filename).exists():
                mask_path = root / "masks_deepcell" / filename
            if not img_path.exists() or not mask_path.exists():
                raise ValueError(f"Missing image/mask: {img_path} or {mask_path}")
            stack, mask = tifffile.imread(img_path), tifffile.imread(mask_path)
            if stack.ndim == 2: stack = stack[None]
            elif stack.shape[-1] == len(panel): stack = np.moveaxis(stack, -1, 0)
            if stack.shape != (len(panel), int(image.height_px), int(image.width_px)) or mask.shape != stack.shape[1:]:
                raise ValueError(f"Image/mask dimensions or panel order differ for {image.roi_id}")
            edge = np.zeros(mask.shape, bool)
            edge[1:] |= mask[1:] != mask[:-1]; edge[:-1] |= mask[:-1] != mask[1:]
            edge[:, 1:] |= mask[:, 1:] != mask[:, :-1]; edge[:, :-1] |= mask[:, :-1] != mask[:, 1:]
            edge &= mask > 0
            rgb = np.zeros((*mask.shape, 3)); legend = []
            for i, marker in enumerate(channels):
                channel = panel.name.tolist().index(marker)
                values = stack[channel].astype(float); lo, hi = np.percentile(values, [1, 99.8])
                normalized = np.clip((values - lo) / (hi - lo), 0, 1) if hi > lo else values * 0
                color = np.eye(3)[i]; rgb += normalized[..., None] * color
                threshold = threshold_otsu(values)
                positive = values > threshold
                signal, noise = values[positive].mean() if positive.any() else np.nan, values[~positive].mean()
                rows.append({"roi_id": image.roi_id, "marker": marker, "signal": signal, "noise": noise,
                             "snr": signal / noise if noise > 0 else np.inf, "threshold": threshold})
                legend.append(Patch(color=color, label=marker))
            rgb[edge] = 1
            fig, ax = plt.subplots(figsize=(8, 8)); ax.imshow(np.clip(rgb, 0, 1)); ax.axis("off")
            ax.legend(handles=legend + [Patch(color="white", label="Cell outline")], loc="upper left", bbox_to_anchor=(1, 1))
            report.save_fig(fig, f"segmentation_{image.roi_id}", "Percentile-scaled marker channels with white mask outlines.")
        return rows

