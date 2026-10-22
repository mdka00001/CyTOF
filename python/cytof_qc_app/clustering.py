"""Cluster an integrated AnnData object and report its marker patterns."""
import argparse
import json
from pathlib import Path

import numpy as np

from .base import QCObject


class ClusteringInput(QCObject):
    def run(self):
        import anndata as ad

        data = ad.read_h5ad(self.args.input)
        aliases = {"pca": "X_pca", "scvi": "X_scVI", "harmony": "X_pca_harmony"}
        requested = aliases.get(self.args.embedding.lower(), self.args.embedding)
        if requested == "auto":
            requested = next((key for key in ("X_scVI", "X_pca_harmony", "X_pca") if key in data.obsm), None)
        if requested not in data.obsm:
            raise ValueError(f"Embedding '{requested}' is missing from the input AnnData object")
        if data.n_obs < 3:
            raise ValueError("Clustering requires at least three cells")
        embedding = np.asarray(data.obsm[requested])
        if embedding.ndim != 2 or not np.isfinite(embedding).all():
            raise ValueError(f"Embedding '{requested}' must be a finite cell-by-dimension matrix")
        return data, requested


class ClusteringMethods(QCObject):
    def run(self, data, embedding):
        import scanpy as sc

        graph_key = "neighbors_clustering"
        sc.pp.neighbors(data, use_rep=embedding,
                        n_neighbors=min(self.args.neighbors, data.n_obs - 1), key_added=graph_key)
        sc.tl.umap(data, neighbors_key=graph_key, random_state=self.args.seed)
        outputs = []
        for method in (name.strip().lower() for name in self.args.methods.split(",")):
            if method == "leiden":
                key = "leiden"
                try:
                    sc.tl.leiden(data, neighbors_key=graph_key, resolution=self.args.leiden_resolution,
                                 random_state=self.args.seed, key_added=key)
                except ImportError as error:
                    raise ImportError("Leiden requires Scanpy's Leiden dependencies; install scanpy[leiden]") from error
            elif method == "phenograph":
                try:
                    options = {
                        "k": min(self.args.phenograph_k, data.n_obs - 1),
                        "directed": self.args.phenograph_directed,
                        "prune": self.args.phenograph_prune,
                        "min_cluster_size": self.args.phenograph_min_cluster_size,
                        "jaccard": self.args.phenograph_jaccard,
                        "primary_metric": self.args.phenograph_primary_metric,
                        "n_jobs": self.args.phenograph_n_jobs,
                        "q_tol": self.args.phenograph_q_tol,
                        "louvain_time_limit": self.args.phenograph_louvain_time_limit,
                        "nn_method": self.args.phenograph_nn_method,
                        "resolution_parameter": self.args.leiden_resolution if self.args.phenograph_resolution is None else self.args.phenograph_resolution,
                        "n_iterations": self.args.phenograph_n_iterations,
                        "use_weights": self.args.phenograph_use_weights,
                        "seed": self.args.seed if self.args.phenograph_seed is None else self.args.phenograph_seed,
                    }
                    if self.args.phenograph_partition_type:
                        import leidenalg
                        options["partition_type"] = getattr(leidenalg, self.args.phenograph_partition_type)
                    extra = json.loads(self.args.phenograph_options)
                    if not isinstance(extra, dict):
                        raise ValueError("--phenograph-options must be a JSON object")
                    options.update(extra)
                    communities, graph, score = sc.external.tl.phenograph(
                        data.obsm[embedding], clustering_algo="leiden", **options,
                    )
                except ImportError as error:
                    raise ImportError("PhenoGraph requires the optional package: pip install phenograph") from error
                except AttributeError as error:
                    raise ValueError(f"Unknown Leiden partition type: {self.args.phenograph_partition_type}") from error
                key = "phenograph_leiden"
                data.obs[key] = communities.astype(str)
                data.obsp["phenograph_jaccard_ig"] = graph.tocsr()
                data.uns["phenograph_jaccard_q"] = score
            else:
                raise ValueError(f"Unknown clustering method '{method}'; choose leiden and/or phenograph")
            outputs.append(key)
        if not outputs:
            raise ValueError("Select at least one method with --methods leiden,phenograph")
        return outputs


class ClusteringReport(QCObject):
    def run(self, data, embedding, clusters, out):
        import matplotlib.pyplot as plt
        import scanpy as sc

        plots = out / "plots"
        plots.mkdir(parents=True, exist_ok=True)
        for cluster in clusters:
            sc.pl.umap(data, color=cluster, show=False)
            plt.gcf().savefig(plots / f"{cluster}_umap.png", dpi=180, bbox_inches="tight")
            plt.close("all")
            layer = next((key for key in ("arcsinh_zscore", "arcsinh") if key in data.layers), None)
            plot_args = {"use_raw": False, "show": False, "standard_scale": "var"}
            if layer:
                plot_args["layer"] = layer
            sample = self._sample(data.obs[cluster], self.args.seed)
            for name, plot in (("heatmap", sc.pl.heatmap), ("dotplot", sc.pl.dotplot), ("matrixplot", sc.pl.matrixplot)):
                current = data[sample].copy() if name == "heatmap" else data
                options = dict(plot_args)
                if name == "heatmap":
                    options.update(swap_axes=True, figsize=(12, 30), cmap="RdBu_r")
                plot(current, var_names=data.var_names.tolist(), groupby=cluster, **options)
                plt.gcf().savefig(plots / f"{cluster}_{name}.png", dpi=180, bbox_inches="tight")
                plt.close("all")

    @staticmethod
    def _sample(labels, seed):
        groups = labels.astype("category")
        per_group = max(1, 500 // groups.nunique())
        rng = np.random.default_rng(seed)
        return np.concatenate([
            rng.choice(np.flatnonzero(groups.to_numpy() == label),
                       min(per_group, (groups == label).sum()), replace=False)
            for label in groups.cat.categories
        ])


class ClusteringPipeline(QCObject):
    def run(self):
        out = Path(self.args.output).expanduser().resolve()
        if out.exists() and any(out.iterdir()):
            raise ValueError("Output directory must be empty")
        out.mkdir(parents=True, exist_ok=True)
        data, embedding = ClusteringInput(self.args).run()
        clusters = ClusteringMethods(self.args).run(data, embedding)
        ClusteringReport(self.args).run(data, embedding, clusters, out)
        data.uns["clustering"] = {"embedding": embedding, "methods": clusters}
        data.write_h5ad(out / "clustered.h5ad")
        print(f"Complete: {out}")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Cluster and plot an AnnData object")
    parser.add_argument("--input", required=True, help="Integrated or PCA-embedded .h5ad object")
    parser.add_argument("--output", required=True, help="Empty output directory")
    parser.add_argument("--embedding", default="auto", help="auto, pca, scvi, harmony, or an obsm key")
    parser.add_argument("--methods", default="leiden,phenograph", help="Comma-separated leiden and/or phenograph")
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--phenograph-k", type=int, default=30)
    parser.add_argument("--phenograph-directed", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--phenograph-prune", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--phenograph-min-cluster-size", type=int, default=10)
    parser.add_argument("--phenograph-jaccard", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--phenograph-primary-metric", choices=("euclidean", "manhattan", "correlation", "cosine"), default="euclidean")
    parser.add_argument("--phenograph-n-jobs", type=int, default=-1)
    parser.add_argument("--phenograph-q-tol", type=float, default=0.001)
    parser.add_argument("--phenograph-louvain-time-limit", type=int, default=2000)
    parser.add_argument("--phenograph-nn-method", choices=("kdtree", "brute"), default="kdtree")
    parser.add_argument("--phenograph-partition-type", default="", help="leidenalg partition class name")
    parser.add_argument("--phenograph-resolution", type=float, default=None)
    parser.add_argument("--phenograph-n-iterations", type=int, default=-1)
    parser.add_argument("--phenograph-seed", type=int, default=None)
    parser.add_argument("--phenograph-use-weights", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--phenograph-options", default="{}", help="Additional Leiden options as a JSON object")
    parser.add_argument("--leiden-resolution", "--resolution", dest="leiden_resolution", type=float, default=1,
                        help="Leiden cluster granularity (higher values create more clusters)")
    parser.add_argument("--seed", type=int, default=220225)
    args = parser.parse_args(argv)
    try:
        ClusteringPipeline(args).run()
    except Exception as error:
        parser.exit(1, f"ERROR: {error}\n")
