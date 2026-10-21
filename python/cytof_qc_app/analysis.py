"""Transformations, marker review and Scanpy embeddings."""
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.mixture import GaussianMixture

from .base import QCObject


class Analysis(QCObject):
    def run(self, adata, markers, selected):
        keep = adata.obs.kept.to_numpy()
        data = adata[keep].copy()
        raw = np.asarray(data.layers["counts"])
        group = data.obs.get("patient_id", data.obs["sample_name"]).astype(str)
        for method in filter(None, self.args.transforms.split(",")):
            matrix = np.arcsinh(raw / self.args.cofactor)
            if method.endswith("_zscore"):
                sd = matrix.std(axis=0, ddof=1)
                sd[~np.isfinite(sd) | (sd == 0)] = 1
                matrix = (matrix - matrix.mean(axis=0)) / sd
            elif method not in {"arcsinh", "rlog"}:
                raise ValueError(f"Unknown transform: {method}")
            if method.startswith("rlog"):
                raise ValueError("Python implementation does not approximate DESeq2 rlog")
            transformed = data.copy()
            transformed.X = matrix.astype(np.float32)
            transformed.layers[method] = transformed.X.copy()
            snr = self._snr(matrix, markers)
            snr["transform"] = method
            if not self.args.skip_embeddings:
                self._embeddings(transformed, selected)
            yield method, transformed, snr, group

    def _embeddings(self, data, selected):
        mask = selected & (np.var(data.X, axis=0) > 0)
        embedding = data[:, mask].copy()
        if embedding.n_vars < 2 or embedding.n_obs < 5:
            raise ValueError("UMAP/t-SNE need two variable markers and five cells")
        sc.pp.neighbors(embedding, n_neighbors=min(self.args.neighbors, embedding.n_obs - 1), n_pcs=0, use_rep="X")
        sc.tl.umap(embedding, random_state=self.args.seed)
        sc.tl.tsne(embedding, use_rep="X", perplexity=min(self.args.perplexity, max(1, (embedding.n_obs - 1) // 3)), random_state=self.args.seed)
        for key in ("X_umap", "X_tsne"):
            data.obsm[key] = embedding.obsm[key]

    def _snr(self, matrix, markers):
        rng = np.random.default_rng(self.args.seed)
        ids = rng.choice(len(matrix), min(len(matrix), self.args.snr_cells), replace=False)
        rows = []
        for j, marker in enumerate(markers):
            all_values = matrix[:, j]
            scale = all_values.std(ddof=1)
            if not np.isfinite(scale) or scale == 0:
                values = np.zeros(len(ids))
            else:
                # Scale first for numerical stability; z-scoring does not
                # change a signal/noise separation score.
                values = ((all_values - all_values.mean()) / scale)[ids]
            row = {"marker": marker, "signal_mean": np.nan, "background_mean": np.nan,
                   "noise_sd": np.nan, "snr": np.nan, "status": "ok"}
            if np.unique(values).size < 3:
                row["status"] = "insufficient_variation"
            else:
                try:
                    groups = GaussianMixture(2, random_state=self.args.seed, n_init=5).fit_predict(values[:, None])
                    means = [values[groups == i].mean() for i in (0, 1)]
                    positive = groups == int(np.argmax(means))
                    signal, background = values[positive], values[~positive]
                    noise_sd = background.std(ddof=1) if len(background) > 1 else np.nan
                    score = (signal.mean() - background.mean()) / noise_sd if noise_sd > 1e-4 else np.nan
                    row.update(signal_mean=signal.mean(), background_mean=background.mean(),
                               noise_sd=noise_sd, snr=score)
                    if not np.isfinite(score):
                        row["status"] = "near_zero_background_variation"
                except Exception:
                    row["status"] = "mixture_failed"
            row["review"] = ("Review: SNR unavailable" if not np.isfinite(row["snr"]) else
                             "Potential exclusion: weak separation" if row["snr"] < self.args.low_snr else
                             "High separation; inspect specificity" if row["snr"] >= self.args.high_snr else
                             "Intermediate separation")
            rows.append(row)
        return pd.DataFrame(rows)

