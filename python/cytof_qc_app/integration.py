"""Integrate multiple AnnData objects with scVI and/or Harmony."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

from .base import QCObject


class IntegrationInput(QCObject):
    def run(self):
        import anndata as ad
        objects = [ad.read_h5ad(path) for path in self.args.inputs]
        markers = objects[0].var_names
        if any(not markers.isin(obj.var_names).all() for obj in objects[1:]):
            markers = markers.intersection(objects[0].var_names, sort=False)
            for obj in objects[1:]:
                markers = markers.intersection(obj.var_names, sort=False)
        if len(markers) < 2:
            raise ValueError("Input objects must share at least two markers")
        batches = ([name.strip() for name in self.args.batch_names.split(",")]
                   if self.args.batch_names else [Path(path).stem for path in self.args.inputs])
        if len(batches) != len(objects) or any(not name for name in batches):
            raise ValueError("--batch-names must provide one non-empty name per input object")
        if len(set(batches)) != len(batches):
            raise ValueError("Input file names must have unique stems")
        for obj, batch in zip(objects, batches):
            obj._inplace_subset_var(markers)
            obj.obs["batch"] = batch
        data = ad.concat(objects, join="inner", label=None, index_unique="::")
        data.obs_names.name = "obs_id"
        data.var = objects[0].var.loc[markers].copy()
        return data


class SCVIIntegration(QCObject):
    def run(self, data, out):
        try:
            import scvi
        except ImportError as error:
            raise ImportError("scvi-tools is required for --methods scvi; install it in this environment") from error
        if self.args.layer not in data.layers:
            raise ValueError(f"scVI layer '{self.args.layer}' is missing; inputs need a raw counts layer")
        import torch
        scvi.settings.seed = self.args.seed
        torch.set_float32_matmul_precision("high")
        scvi.model.SCVI.setup_anndata(data, layer=self.args.layer, batch_key="batch")
        model = scvi.model.SCVI(data, n_layers=2, n_latent=30, gene_likelihood="nb")
        model.train()
        data.obsm["X_scVI"] = model.get_latent_representation().astype(np.float32)
        model.save(out / "models/scvi", overwrite=True, save_anndata=False)
        return data.obsm["X_scVI"]


class HarmonyIntegration(QCObject):
    def run(self, data, out):
        import scanpy as sc
        try:
            import harmonypy
        except ImportError as error:
            raise ImportError("harmonypy is required for --methods harmony; install it in this environment") from error
        sc.pp.pca(data, n_comps=min(self.args.n_pcs, data.n_vars - 1, data.n_obs - 1))
        result = harmonypy.run_harmony(data.obsm["X_pca"], data.obs, "batch", random_state=self.args.seed)
        corrected = np.asarray(result.Z_corr)
        if corrected.shape == (data.obsm["X_pca"].shape[1], data.n_obs):
            corrected = corrected.T
        if corrected.shape[0] != data.n_obs:
            raise ValueError(f"Harmony returned unexpected embedding shape: {corrected.shape}")
        data.obsm["X_pca_harmony"] = corrected
        return data.obsm["X_pca_harmony"]


class IntegrationReport(QCObject):
    def run(self, data, methods, out):
        import matplotlib.pyplot as plt
        import scanpy as sc
        plots, tables = out / "plots", out / "tables"
        plots.mkdir(parents=True, exist_ok=True)
        tables.mkdir(parents=True, exist_ok=True)
        if self.args.cofactor <= 0:
            raise ValueError("--cofactor must be greater than zero")
        if "arcsinh_zscore" in data.layers:
            marker_layer = "arcsinh_zscore"
        elif "arcsinh" in data.layers:
            marker_layer = "arcsinh"
        else:
            expression = data.X.toarray() if hasattr(data.X, "toarray") else data.X
            data.layers["arcsinh"] = np.arcsinh(expression / self.args.cofactor).astype(np.float32)
            marker_layer = "arcsinh"
        group_cols = [c for c in self.args.group_cols.split(",") if c in data.obs]
        roi_col = self.args.roi_col if self.args.roi_col in data.obs else None
        group_cols = list(dict.fromkeys(["batch", *group_cols, *([roi_col] if roi_col else [])]))
        for method, rep in methods.items():
            neighbors = f"neighbors_{method}"
            sc.pp.neighbors(data, use_rep=rep, n_neighbors=min(self.args.neighbors, data.n_obs - 1), key_added=neighbors)
            sc.tl.umap(data, neighbors_key=neighbors, random_state=self.args.seed)
            data.obsm[f"X_umap_{method}"] = data.obsm["X_umap"].copy()
            sc.tl.tsne(data, use_rep=rep,
                       perplexity=min(self.args.perplexity, max(1, (data.n_obs - 1) // 3)),
                       random_state=self.args.seed)
            data.obsm[f"X_tsne_{method}"] = data.obsm["X_tsne"].copy()
            for color in group_cols:
                fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
                for ax, key, title in zip(axes, (f"X_umap_{method}", f"X_tsne_{method}"), ("UMAP", "t-SNE")):
                    coords = data.obsm[key]
                    categories = data.obs[color].astype("category")
                    for category in categories.cat.categories:
                        mask = categories == category
                        ax.scatter(coords[mask, 0], coords[mask, 1], s=3, alpha=.65, label=str(category), rasterized=True)
                    ax.set(title=f"{title} by {color}", xlabel=f"{title}1", ylabel=f"{title}2")
                    ax.legend(markerscale=3, fontsize=6, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
                fig.tight_layout()
                fig.savefig(plots / f"{method}_{color}_embeddings.png", dpi=180, bbox_inches="tight")
                plt.close(fig)
                coords = pd.DataFrame(data.obsm[f"X_umap_{method}"], index=data.obs_names, columns=["UMAP1", "UMAP2"])
                coords[["TSNE1", "TSNE2"]] = data.obsm[f"X_tsne_{method}"]
                coords.insert(0, "cell_id", coords.index)
                coords[color] = data.obs[color].to_numpy()
                coords.to_csv(tables / f"{method}_{color}_coordinates.csv", index=False)
            for group in group_cols:
                categories = data.obs[group].astype("category")
                per_group = max(1, 500 // categories.nunique())
                rng = np.random.default_rng(self.args.seed)
                sample = np.concatenate([
                    rng.choice(np.flatnonzero(categories.to_numpy() == label),
                               min(per_group, (categories == label).sum()), replace=False)
                    for label in categories.cat.categories
                ])
                for plot in (sc.pl.heatmap, sc.pl.dotplot, sc.pl.matrixplot):
                    name = plot.__name__.replace("plot", "plot")
                    current = data[sample].copy() if name == "heatmap" else data
                    options = {"standard_scale": "var", "cmap": "RdBu_r"}
                    if name == "heatmap":
                        options.update(swap_axes=True, figsize=(12, 30))
                    plot(current, var_names=data.var_names.tolist(), groupby=group,
                         layer=marker_layer, use_raw=False, show=False, **options)
                    fig = plt.gcf()
                    fig.savefig(plots / f"{method}_{group}_{name}.png", dpi=180, bbox_inches="tight")
                    plt.close(fig)
            score = self._snr(data.layers[marker_layer], data.var_names)
            score.to_csv(tables / f"{method}_snr.csv", index=False)
            fig, ax = plt.subplots(figsize=(max(6, len(score) * .22), 4))
            ax.bar(score.marker, score.snr, color="#3b82a0")
            ax.set(ylabel="SNR", xlabel="Marker", title=f"{method} marker SNR")
            ax.tick_params(axis="x", rotation=90, labelsize=7)
            fig.tight_layout()
            fig.savefig(plots / f"{method}_snr.png", dpi=180)
            plt.close(fig)

    def _snr(self, matrix, markers):
        if hasattr(matrix, "toarray"):
            matrix = matrix.toarray()
        x = np.asarray(matrix)
        rng = np.random.default_rng(self.args.seed)
        rows = []
        ids = rng.choice(len(x), min(len(x), self.args.snr_cells), replace=False)
        for j, marker in enumerate(markers):
            values = x[:, j]
            sd = values.std(ddof=1)
            z = ((values - values.mean()) / sd)[ids] if sd > 0 else np.zeros(len(ids))
            snr = np.nan
            if np.unique(z).size >= 3:
                try:
                    labels = GaussianMixture(2, random_state=self.args.seed, n_init=5).fit_predict(z[:, None])
                    means = [z[labels == i].mean() for i in (0, 1)]
                    high = labels == int(np.argmax(means))
                    noise = z[~high].std(ddof=1)
                    if noise > 1e-4:
                        snr = (z[high].mean() - z[~high].mean()) / noise
                except (ValueError, FloatingPointError):
                    pass
            rows.append({"marker": marker, "snr": snr})
        return pd.DataFrame(rows)


class IntegrationPipeline(QCObject):
    def run(self):
        out = Path(self.args.output).expanduser().resolve()
        if out.exists() and any(out.iterdir()):
            raise ValueError("Output directory must be empty")
        out.mkdir(parents=True, exist_ok=True)
        (out / "models").mkdir(exist_ok=True)
        data = IntegrationInput(self.args).run()
        methods = {}
        for method in self.args.methods.split(","):
            method = method.strip().lower()
            if method == "scvi":
                methods[method] = "X_scVI"
                SCVIIntegration(self.args).run(data, out)
            elif method == "harmony":
                methods[method] = "X_pca_harmony"
                HarmonyIntegration(self.args).run(data, out)
            else:
                raise ValueError(f"Unknown integration method: {method}")
        if not methods:
            raise ValueError("Select at least one method")
        IntegrationReport(self.args).run(data, methods, out)
        data.write_h5ad(out / "integrated.h5ad")
        print(f"Complete: {out}")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Integrate multiple AnnData objects")
    parser.add_argument("--inputs", nargs="+", required=True, help="Input .h5ad files; each file is one batch")
    parser.add_argument("--output", required=True)
    parser.add_argument("--methods", default="scvi,harmony", help="Comma-separated scvi and/or harmony")
    parser.add_argument("--batch-names", default="", help="Comma-separated batch names in input-file order")
    parser.add_argument("--group-cols", default="patient_id,sample_name")
    parser.add_argument("--roi-col", default="roi_id")
    parser.add_argument("--layer", default="counts", help="Raw count layer used by scVI")
    parser.add_argument("--cofactor", type=float, default=1, help="Cofactor for arcsinh marker plots and SNR")
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--perplexity", type=float, default=30)
    parser.add_argument("--snr-cells", type=int, default=30000)
    parser.add_argument("--seed", type=int, default=220225)
    args = parser.parse_args(argv)
    try:
        IntegrationPipeline(args).run()
    except Exception as error:
        parser.exit(1, f"ERROR: {error}\n")
