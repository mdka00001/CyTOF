"""Small HTML/PDF report and marker plots."""
import html
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.backends.backend_pdf import PdfPages

from .base import QCObject
from .images import ImageQC


class QCReport(QCObject):
    def run(self, out, data, panel, selected, cells, images, snr, transforms):
        self.out, self.items = out, []
        with PdfPages(out / "report.pdf") as pdf:
            self.pdf = pdf
            self.page("Run overview", f"{len(panel)} markers; {len(cells)} input cells; {data.n_obs} retained cells; {len(images)} ROIs.",
                      lambda ax: ax.text(.05, .8, "SNR and QC flags are suggestions. Review marker localization and image quality.", wrap=True))
            metrics = pd.DataFrame([{"roi_id": r.roi_id,
                "cells_before": int((cells.roi_id == r.roi_id).sum()),
                "cells_after": int(((cells.roi_id == r.roi_id) & cells.kept).sum()),
                "coverage": cells.loc[cells.roi_id == r.roi_id, "area"].sum() / (r.width_px * r.height_px),
                "cells_per_mm2": int(((cells.roi_id == r.roi_id) & cells.kept).sum()) / (r.width_px * r.height_px * self.args.pixel_size ** 2 / 1e6)} for r in images.itertuples()])
            metrics["roi_label"] = [f"ROI {i:02d}" for i in range(1, len(metrics) + 1)]
            metrics.to_csv(out / "tables" / "roi_qc.csv", index=False)
            for col in ("coverage", "cells_before", "cells_after", "cells_per_mm2"):
                self.page(f"ROI {col}", f"Per-ROI {col}; ROI labels show the final numeric ID. Full IDs are in tables/roi_qc.csv.",
                          lambda ax, c=col: self.roi_plot(ax, metrics, c))
            for method, transformed, snr, group in transforms:
                snr.to_csv(out / "tables" / f"{method}_cell_snr.csv", index=False)
                self.page(f"{method} cell SNR", "Two-component Gaussian mixture fitted to this transformed matrix. SNR = (high-component mean - low-component mean) / low-component standard deviation. Values are computed on each transformed matrix after per-marker scaling for numerical stability. The score is scale-invariant, so arcsinh and its z-score should agree. Axes use log2 as in the R report; inspect the CSV for exact values.",
                          lambda ax, d=snr: self.snr_plot(ax, d))
                matrix = np.asarray(transformed.X)
                self.marker_plots(transformed, panel, selected, group, method)
                self.marker_distributions(matrix, panel, selected, group, method)
                for key, label in (("X_umap", "UMAP"), ("X_tsne", "t-SNE")):
                    if key not in transformed.obsm:
                        continue
                    self.page(f"{method} {label} by patient", "All retained cells; patient_id, falling back to sample_name.",
                              lambda ax, k=key, d=transformed, g=group: self.embedding(ax, d.obsm[k], g))
                    self.embedding_markers(transformed, panel, selected, method, key, label)
            pd.concat([result[2] for result in transforms], ignore_index=True).to_csv(out / "tables/cell_snr_marker_review.csv", index=False)
            if not self.args.skip_images:
                pixels = ImageQC(self.args).run(Path(self.args.input).expanduser().resolve(), panel, images, cells, self.args.image_channels, self)
                pixels = pd.DataFrame(pixels)
                pixels.to_csv(out / "tables/pixel_snr.csv", index=False)
                for suffix, subset in (("all", pixels), ("signal_filtered", pixels[pixels.signal > self.args.min_signal])):
                    summary = subset.groupby("marker", as_index=False).agg(signal=("signal", "mean"), snr=("snr", "mean"))
                    summary.to_csv(out / "tables" / f"pixel_snr_summary_{suffix}.csv", index=False)
                    self.page(f"Pixel SNR {suffix}", "Per-marker ROI means. Both axes use log2 as in the R report; the filtered view excludes ROI-marker pairs with weak positive signal.",
                              lambda ax, d=summary: self.snr_plot(ax, d.rename(columns={"signal": "signal_mean"}), "Signal mean (log2)"))
        content = ["<!doctype html><meta charset='utf-8'><title>CyTOF QC</title><h1>CyTOF QC</h1><p><a href='report.pdf'>PDF report</a></p>"]
        content += [f"<h2>{html.escape(t)}</h2><p>{html.escape(e)}</p><img style='width:100%' src='{p}'>" for t, e, p in self.items]
        (out / "report.html").write_text("\n".join(content), encoding="utf-8")
        return metrics

    @staticmethod
    def roi_plot(ax, metrics, column):
        ax.barh(metrics.roi_label, metrics[column], color="#8a8a8a", edgecolor="#222", linewidth=.4)
        ax.invert_yaxis()
        ax.set_ylabel("ROI suffix (full ID in table)")
        ax.set_xlabel(column.replace("_", " "))
        ax.grid(axis="x", color="#ddd", linewidth=.6)
        ax.set_axisbelow(True)

    @staticmethod
    def snr_plot(ax, data, signal_label="Positive component mean (log2 transformed units)"):
        if data.empty:
            ax.text(.5, .5, "No markers meet the signal filter", ha="center", va="center", transform=ax.transAxes)
            return
        good = data[(data.signal_mean > 0) & (data.snr > 0) & np.isfinite(data.signal_mean) & np.isfinite(data.snr)]
        if good.empty:
            ax.text(.5, .5, "No finite positive SNR values", ha="center", va="center", transform=ax.transAxes)
            return
        x, y = np.log2(good.signal_mean), np.log2(good.snr)
        ax.scatter(x, y, color="#286b9e", s=24)
        for xi, yi, marker in zip(x, y, good.marker):
            ax.annotate(marker, (xi, yi), xytext=(3, 3), textcoords="offset points", fontsize=7)
        ax.set_xlabel(signal_label)
        ax.set_ylabel("SNR (log2)")
        ax.grid(color="#ddd", linewidth=.6)

    def page(self, title, explanation, draw):
        fig, ax = plt.subplots(figsize=(12, 8)); draw(ax); fig.suptitle(title); fig.tight_layout()
        image = f"plots/{len(self.items)+1:03d}_{re.sub('[^A-Za-z0-9_-]', '_', title)}.png"
        fig.savefig(self.out / image, dpi=120, bbox_inches="tight"); self.pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
        self.items.append((title, explanation, image))

    def marker_plots(self, data, panel, selected, group, method):
        markers = panel.name.to_numpy()[selected]
        plot_data = data[:, markers].copy(); plot_data.obs["qc_group"] = pd.Categorical(group.astype(str))
        # Keep individual-cell heatmaps readable; grouped Scanpy plots use all cells.
        groups = plot_data.obs.qc_group
        rng = np.random.default_rng(self.args.seed)
        per_group = max(1, 500 // groups.nunique())
        sample = np.concatenate([rng.choice(np.flatnonzero(groups.to_numpy() == label),
                                            min(per_group, (groups == label).sum()), replace=False)
                                 for label in groups.cat.categories])
        for kind, plot in (("heatmap", sc.pl.heatmap), ("dotplot", sc.pl.dotplot), ("matrixplot", sc.pl.matrixplot)):
            current = plot_data[sample].copy() if kind == "heatmap" else plot_data
            with plt.rc_context({"figure.figsize": (max(12, len(markers) * .4), 7), "font.size": 6}):
                if kind == "heatmap":
                    plot(current, markers, groupby="qc_group", show=False, swap_axes=True, standard_scale="var", cmap="RdBu_r")
                else:
                    plot(current, markers, groupby="qc_group", show=False, standard_scale="var")
                fig = plt.gcf(); fig.suptitle(f"{method} marker {kind}")
                self.save_fig(fig, f"{method}_{kind}", f"Grouped Scanpy {kind} for selected markers.")

    def marker_distributions(self, matrix, panel, selected, group, method):
        markers, groups = panel.name.to_numpy()[selected], group.astype(str).to_numpy()
        labels = pd.unique(groups)
        for start in range(0, len(markers), 20):
            current = markers[start:start + 20]
            fig, axes = plt.subplots(5, 4, figsize=(20, 18), squeeze=False)
            for offset, marker in enumerate(current):
                ax = axes.flat[offset]
                index = np.flatnonzero(panel.name.to_numpy() == marker)[0]
                ax.boxplot([matrix[groups == label, index] for label in labels], showfliers=False)
                ax.set_title(marker, fontsize=9)
                ax.set_xticks(range(1, len(labels) + 1), labels, rotation=45, ha="right", fontsize=7)
                ax.grid(axis="y", color="#ddd", linewidth=.5)
            for ax in axes.flat[len(current):]: ax.set_axis_off()
            fig.suptitle(f"{method} marker distributions by patient/sample", fontsize=14)
            fig.tight_layout(rect=(0, 0, 1, .97))
            self.save_fig(fig, f"{method}_marker_distributions_{start // 20 + 1}",
                          "Grouped by patient_id, falling back to sample_name; at most 20 markers per page.")

    def embedding(self, ax, coords, group):
        for label in pd.unique(group):
            mask = np.asarray(group == label)
            ax.scatter(coords[mask, 0], coords[mask, 1], s=3, label=label)
        ax.legend(markerscale=3, fontsize=7)

    def embedding_markers(self, data, panel, selected, method, key, label):
        markers = panel.name.to_numpy()[selected]; coords = data.obsm[key]
        cols, rows = min(5, len(markers)), math.ceil(len(markers) / 5)
        fig, axes = plt.subplots(rows, cols, figsize=(cols * 3, rows * 2.5), squeeze=False)
        for i, marker in enumerate(markers):
            axes.flat[i].scatter(coords[:, 0], coords[:, 1], c=data[:, marker].X.toarray().ravel() if hasattr(data[:, marker].X, "toarray") else data[:, marker].X.ravel(), s=2, cmap="viridis")
            axes.flat[i].set_title(marker, fontsize=8); axes.flat[i].set_axis_off()
        for ax in axes.flat[len(markers):]: ax.set_axis_off()
        fig.suptitle(f"{method} {label} marker expression"); self.save_fig(fig, f"{method}_{label}_markers", "Selected marker expression with continuous viridis scales.")

    def save_fig(self, fig, name, explanation):
        filename = f"plots/{len(self.items)+1:03d}_{name}.png"
        fig.savefig(self.out / filename, dpi=120, bbox_inches="tight"); self.pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
        self.items.append((name, explanation, filename))

