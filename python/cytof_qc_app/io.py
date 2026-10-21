"""Read Steinbock tables and create the retained AnnData object."""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from .base import QCObject


class DataInput(QCObject):
    def run(self):
        root = Path(self.args.input).expanduser().resolve()
        panel = pd.read_csv(root / "panel.csv")
        panel = panel[panel.get("keep", 1).astype(str).str.lower().isin(["1", "true"])].copy() if "keep" in panel else panel.copy()
        panel[["name", "channel"]] = panel[["name", "channel"]].astype(str).apply(lambda c: c.str.strip())
        if panel.empty or panel.name.duplicated().any() or panel.channel.duplicated().any():
            raise ValueError("Retained panel names and channels must be unique")
        images = pd.read_csv(root / "images.csv")
        images["roi_id"] = images.image.map(lambda p: Path(str(p)).stem)
        mapping = pd.read_csv(self.args.samples).set_index("roi_id")
        missing = set(images.roi_id) - set(mapping.index)
        if missing:
            raise ValueError("Unmapped ROI IDs: " + ", ".join(sorted(missing)))
        excluded = {x.strip() for x in self.args.exclude.split(",") if x.strip()}
        if excluded - set(panel.name):
            raise ValueError("Unknown excluded markers: " + ", ".join(sorted(excluded - set(panel.name))))
        selected = panel.get("use_channel", pd.Series(True, index=panel.index)).astype(str).str.lower().isin(["1", "true"]).to_numpy()
        selected &= ~panel.name.isin(excluded).to_numpy()

        matrices, metadata = [], []
        for image in images.itertuples(index=False):
            intensity = pd.read_csv(root / self.args.intensities / f"{image.roi_id}.csv")
            props = pd.read_csv(root / self.args.regionprops / f"{image.roi_id}.csv")
            intensity.columns = intensity.columns.astype(str).str.strip()
            markers = panel.name.tolist() if set(panel.name) <= set(intensity) else panel.channel.tolist()
            if "Object" not in intensity or not {"Object", "area"} <= set(props) or not set(markers) <= set(intensity):
                raise ValueError(f"Invalid intensity/region-property tables for {image.roi_id}")
            if intensity.Object.duplicated().any() or props.Object.duplicated().any():
                raise ValueError(f"Duplicate Object IDs for {image.roi_id}")
            props = props.set_index("Object").loc[intensity.Object].reset_index()
            values = intensity[markers].to_numpy(dtype=np.float32)
            if not np.isfinite(values).all() or (values < 0).any():
                raise ValueError(f"Invalid intensities for {image.roi_id}")
            meta = props.copy()
            meta["roi_id"] = image.roi_id
            meta["cell_id"] = [f"{image.roi_id}::{obj}" for obj in intensity.Object]
            for col, value in mapping.loc[image.roi_id].items():
                meta[col] = value
            meta["width_px"], meta["height_px"], meta["image"] = image.width_px, image.height_px, image.image
            matrices.append(values)
            metadata.append(meta)
        cells = pd.concat(metadata, ignore_index=True)
        cells["area"] = cells.area.astype(float)
        cells["kept"] = cells.area >= self.args.min_area
        if cells.cell_id.duplicated().any():
            raise ValueError("Cell IDs must be unique")
        var = panel.set_index("name", drop=False)
        var["selected"] = selected
        counts = np.concatenate(matrices)
        adata = ad.AnnData(X=counts, obs=cells.set_index("cell_id", drop=False), var=var)
        adata.layers["counts"] = counts.copy()
        return root, panel, images, cells, selected, adata

