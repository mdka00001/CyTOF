"""Coordinate the small input, analysis, and report objects."""
import json
from pathlib import Path

import pandas as pd

from .analysis import Analysis
from .base import QCObject
from .io import DataInput
from .report import QCReport


class Pipeline(QCObject):
    def run(self):
        root, panel, images, cells, selected, all_data = DataInput(self.args).run()
        out = Path(self.args.output).expanduser().resolve()
        if out.exists() and any(out.iterdir()):
            raise ValueError("Output directory must be empty")
        for name in ("plots", "tables", "matrices", "coordinates", "objects"):
            (out / name).mkdir(parents=True, exist_ok=True)
        cells.to_csv(out / "tables/cells_all.csv", index=False)
        panel.to_csv(out / "tables/panel.csv", index=False)
        pd.DataFrame(all_data.layers["counts"], columns=panel.name).assign(cell_id=cells.cell_id).to_csv(out / "matrices/counts_all.csv.gz", index=False)
        transforms = list(Analysis(self.args).run(all_data, panel.name.to_numpy(), selected))
        for method, data, _, _ in transforms:
            pd.DataFrame(data.X, columns=panel.name).assign(cell_id=data.obs.index).to_csv(out / f"matrices/{method}.csv.gz", index=False)
            data.write_h5ad(out / f"objects/{method}.h5ad")
            for key, name in (("X_umap", "UMAP"), ("X_tsne", "TSNE")):
                if key in data.obsm:
                    coords = pd.DataFrame(data.obsm[key], columns=[name + "1", name + "2"], index=data.obs.index)
                    coords.insert(0, "cell_id", coords.index)
                    for col in ("roi_id", "sample_name", "patient_id"):
                        if col in data.obs: coords[col] = data.obs[col].to_numpy()
                    coords.to_csv(out / f"coordinates/{method}_{name}.csv", index=False)
        report_data = transforms[0][1]
        for method, data, _, _ in transforms:
            report_data.layers[method] = data.X.copy()
        QCReport(self.args).run(out, report_data, panel, selected, cells, images, transforms[0][2], transforms)
        report_data.write_h5ad(out / "objects/qc_data.h5ad")
        (out / "objects/config.json").write_text(json.dumps(vars(self.args), indent=2))
        (out / "STATUS").write_text("complete\n")
        print(f"Complete: {out}")

