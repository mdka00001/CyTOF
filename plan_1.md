This directory is a github repository that I am trying to implement.

I want to build a qc workflow for IMC CyTOF preprocessed data that I can incorporate with any CyTOF dataset.
My Input shall be the steinbock preprocessed CyTOF data directory. I have some scripts that I used for the postprocessing QC. I want you to read them and build me a modular application reproducing the QC steps in R that can 
return a html+pdf file with all the QC plots, the normalized count matrices as csv or any conventional file format, the umap and tsne coordinates in csv or any conventional file format.

The sample names will be passed as an input csv/txt file to map the ROIs to their name. 

The scripts are in /home/md-adnan-karim/Documents/git_repo/CyTOF/1.qc_script_1.R and /home/md-adnan-karim/Documents/git_repo/CyTOF/1.qc_script_2.R.

The QC workflow was described here in detail: https://bodenmillergroup.github.io/IMCDataAnalysis/image-and-cell-level-quality-control.html

The script covers most of the steps.


## TASK 1
I need you:

1. Implement a very easy, documented, simply modular executable R application that can be run from CLI.
2. The time and space complexiety should be minimal.
3. Include Rlog transformation as well along with arcsinh, arcsinh+zscore. (rlog, rlog+zscore).
4. Make a report writing explanation of all plots.
5. Give a breif breakdown of markers with high snr score and potential to exclude.
6. Return the dimentional reduction coordinates, objects, counts matrix.
7. THe report should be in HTML+PDF format.


## TASK 2
for this task, I want to make this script executable by adding the path to source so I dont need to specify the path every time. Would be best if this can be installed directly to source path and execution ready.


## TASK 3
In this task I want to use the script in /home/md-adnan-karim/Documents/git_repo/CyTOF/1.qc_script_3.R and build another sub-method  separate from the qc method to cytof-qc to construct the fcs files using spe normalized data. Also document it.

## TASK 4
I want to add a python version of the same QC and analysis everything as of above using the scanpy - anndata data structure.
This does not mean to replace the R implementation and keep both executable. Ask for further query regarding the implementation and execution approach.

## TASK 5
This task should implement methods and a separate pipeline execution arguement that takes multiple *.h5ad objects of choice and integrate them to a single object. The pipeline should run two batch correction methods defined by users: 1. scvi integration 2. harmony integration Run the integration by batch ( here each object shall be considered as a separate batch). Output shall return: 1. integrated h5ad file. 2. UMAP+TSNE plot for batch, groups (patient/sample). 3. Shall plot snr plot. 4. shall plot heatmap+dotplot+matrixplot by batch+group+roi. 5. shall store the scvi model. 6. shall store both scvi/harmony embeddings. Use the same abstract object and same code base skeleton for the implementation.

## TASK 6
This task should run scanpy based leiden clustering and scanpy phenograph extension for leiden clustering on anndata objects (*.h5ad). The latent embeddings can be the pca, scvi_embeddings, harmony_embeddings how you named and stored them in previous steps and run neighbours, umaps, clustering and plot them as umap, heatmap, dotplot, matrixplot.
If scvi/harmony embeddings are missing then return error. Make scanpy leiden and phenograph optional for user. Also keep the embedding choice optional for users.
return the new h5ad object.
Parameters:
1. phenograph:
    [31mSignature:[39m
sc.external.tl.phenograph(
    data: [33m'AnnData | np.ndarray | SpBase'[39m,
    clustering_algo: [33m"Literal['louvain', 'leiden'] | None"[39m = [33m'louvain'[39m,
    *,
    k: [33m'int'[39m = [32m30[39m,
    directed: [33m'bool'[39m = [38;5;28;01mFalse[39;00m,
    prune: [33m'bool'[39m = [38;5;28;01mFalse[39;00m,
    min_cluster_size: [33m'int'[39m = [32m10[39m,
    jaccard: [33m'bool'[39m = [38;5;28;01mTrue[39;00m,
    primary_metric: [33m"Literal['euclidean', 'manhattan', 'correlation', 'cosine']"[39m = [33m'euclidean'[39m,
    n_jobs: [33m'int'[39m = -[32m1[39m,
    q_tol: [33m'float'[39m = [32m0.001[39m,
    louvain_time_limit: [33m'int'[39m = [32m2000[39m,
    nn_method: [33m"Literal['kdtree', 'brute']"[39m = [33m'kdtree'[39m,
    partition_type: [33m'type[MutableVertexPartition] | None'[39m = [38;5;28;01mNone[39;00m,
    resolution_parameter: [33m'float'[39m = [32m1[39m,
    n_iterations: [33m'int'[39m = -[32m1[39m,
    use_weights: [33m'bool'[39m = [38;5;28;01mTrue[39;00m,
    seed: [33m'int | None'[39m = [38;5;28;01mNone[39;00m,
    copy: [33m'bool'[39m = [38;5;28;01mFalse[39;00m,
    **kargs: [33m'Any'[39m,
) -> [33m'tuple[np.ndarray | None, SpBase, float | None] | None'[39m
[31mDocstring:[39m
PhenoGraph clustering :cite:p:`Levine2015`.

**PhenoGraph** is a clustering method designed for high-dimensional single-cell
data. It works by creating a graph ("network") representing phenotypic similarities
between cells and then identifying communities in this graph. It supports both
Louvain_ and Leiden_ algorithms for community detection.

.. _Louvain: https://louvain-igraph.readthedocs.io/en/latest/

.. _Leiden: https://leidenalg.readthedocs.io/en/latest/reference.html

.. note::
   More information and bug reports `here
   <https://github.com/dpeerlab/PhenoGraph>`__.

Parameters
----------
data : 'AnnData | np.ndarray | SpBase'
    AnnData, or Array of data to cluster, or sparse matrix of k-nearest neighbor
    graph. If ndarray, n-by-d array of n cells in d dimensions. if sparse matrix,
    n-by-n adjacency matrix.
clustering_algo : "Literal['louvain', 'leiden'] | None", optional (default: 'louvain')
    Choose between `'Louvain'` or `'Leiden'` algorithm for clustering.
k : 'int', optional (default: 30)
    Number of nearest neighbors to use in first step of graph construction.
directed : 'bool', optional (default: False)
    Whether to use a symmetric (default) or asymmetric (`'directed'`) graph.
    The graph construction process produces a directed graph, which is symmetrized
    by one of two methods (see `prune` below).
prune : 'bool', optional (default: False)
    `prune=False`, symmetrize by taking the average between the graph and its
    transpose. `prune=True`, symmetrize by taking the product between the graph
    and its transpose.
min_cluster_size : 'int', optional (default: 10)
    Cells that end up in a cluster smaller than min_cluster_size are considered
    outliers and are assigned to -1 in the cluster labels.
jaccard : 'bool', optional (default: True)
    If `True`, use Jaccard metric between k-neighborhoods to build graph. If
    `False`, use a Gaussian kernel.
primary_metric : "Literal['euclidean', 'manhattan', 'correlation', 'cosine']", optional (default: 'euclidean')
    Distance metric to define nearest neighbors. Note that performance will be
    slower for correlation and cosine.
n_jobs : 'int', optional (default: -1)
    Nearest Neighbors and Jaccard coefficients will be computed in parallel using
    n_jobs. If 1 is given, no parallelism is used. If set to -1, all CPUs are used.
    For n_jobs below -1, `n_cpus + 1 + n_jobs` are used.
q_tol : 'float', optional (default: 0.001)
    Tolerance, i.e. precision, for monitoring modularity optimization.
louvain_time_limit : 'int', optional (default: 2000)
    Maximum number of seconds to run modularity optimization. If exceeded the best
    result so far is returned.
nn_method : "Literal['kdtree', 'brute']", optional (default: 'kdtree')
    Whether to use brute force or kdtree for nearest neighbor search.
    For very large high-dimensional data sets, brute force, with parallel
    computation, performs faster than kdtree.
partition_type : 'type[MutableVertexPartition] | None', optional (default: None)
    Defaults to :class:`~leidenalg.RBConfigurationVertexPartition`. For the
    available options, consult the documentation for
    :func:`~leidenalg.find_partition`.
resolution_parameter : 'float', optional (default: 1)
    A parameter value controlling the coarseness of the clustering in Leiden. Higher
    values lead to more clusters. Set to `None` if overriding `partition_type` to
    one that does not accept a `resolution_parameter`.
n_iterations : 'int', optional (default: -1)
    Number of iterations to run the Leiden algorithm. If the number of iterations is
    negative, the Leiden algorithm is run until an iteration in which there was no
    improvement.
use_weights : 'bool', optional (default: True)
    Use vertices in the Leiden computation.
seed : 'int | None', optional (default: None)
    Leiden initialization of the optimization.
copy : 'bool', optional (default: False)
    Return a copy or write to `adata`.
kargs : 'Any'
    Additional arguments passed to :func:`~leidenalg.find_partition` and the
    constructor of the `partition_type`.

Returns
-------
Depending on `copy`, returns or updates `adata` with the following fields:

**communities** - :class:`~numpy.ndarray` (:attr:`~anndata.AnnData.obs`, dtype `int`)
    integer array of community assignments for each row in data.

**graph** - :class:`~scipy.sparse.spmatrix` (:attr:`~anndata.AnnData.obsp`, dtype `float`)
    the graph that was used for clustering.

**Q** - `float` (:attr:`~anndata.AnnData.uns`, dtype `float`)
    the modularity score for communities on graph.

Example
-------
>>> from anndata import AnnData
>>> import scanpy as sc
>>> import scanpy.external as sce
>>> import numpy as np
>>> import pandas as pd

With annotated data as input:

>>> adata = sc.datasets.pbmc3k()
>>> sc.pp.normalize_per_cell(adata)

Then do PCA:

>>> sc.pp.pca(adata, n_comps=100)

Compute phenograph clusters:

**Louvain** community detection

>>> sce.tl.phenograph(adata, clustering_algo="louvain", k=30)

**Leiden** community detection

>>> sce.tl.phenograph(adata, clustering_algo="leiden", k=30)

Return only `Graph` object

>>> sce.tl.phenograph(adata, clustering_algo=None, k=30)

Now to show phenograph on tSNE (for example):

Compute tSNE:

>>> sc.tl.tsne(adata, random_state=7)

Plot phenograph clusters on tSNE:

>>> sc.pl.tsne(
...     adata,
...     color=["pheno_louvain", "pheno_leiden"],
...     s=100,
...     palette=sc.pl.palettes.vega_20_scanpy,
...     legend_fontsize=10,
... )

Cluster and cluster centroids for input Numpy ndarray

>>> df = np.random.rand(1000, 40)
>>> dframe = pd.DataFrame(df)
>>> dframe.index, dframe.columns = (
...     map(str, dframe.index),
...     map(str, dframe.columns),
... )
>>> adata = AnnData(dframe)
>>> sc.pp.pca(adata, n_comps=20)
>>> sce.tl.phenograph(adata, clustering_algo="leiden", k=50)
>>> sc.tl.tsne(adata, random_state=1)
>>> sc.pl.tsne(
...     adata,
...     color=["pheno_leiden"],
...     s=100,
...     palette=sc.pl.palettes.vega_20_scanpy,
...     legend_fontsize=10,
... )

2. sc.tl.leiden:
    [31mSignature:[39m
sc.tl.leiden(
    adata: [33m'AnnData'[39m,
    resolution: [33m'float'[39m = [32m1[39m,
    *,
    restrict_to: [33m'tuple[str, Sequence[str]] | None'[39m = [38;5;28;01mNone[39;00m,
    random_state: [33m'_LegacyRandom'[39m = [32m0[39m,
    key_added: [33m'str'[39m = [33m'leiden'[39m,
    adjacency: [33m'CSBase | None'[39m = [38;5;28;01mNone[39;00m,
    directed: [33m'bool | None'[39m = [38;5;28;01mNone[39;00m,
    use_weights: [33m'bool'[39m = [38;5;28;01mTrue[39;00m,
    n_iterations: [33m'int'[39m = -[32m1[39m,
    partition_type: [33m'type[MutableVertexPartition] | None'[39m = [38;5;28;01mNone[39;00m,
    neighbors_key: [33m'str | None'[39m = [38;5;28;01mNone[39;00m,
    obsp: [33m'str | None'[39m = [38;5;28;01mNone[39;00m,
    copy: [33m'bool'[39m = [38;5;28;01mFalse[39;00m,
    flavor: [33m"Literal['leidenalg', 'igraph'] | None"[39m = [38;5;28;01mNone[39;00m,
    **clustering_args,
) -> [33m'AnnData | None'[39m
[31mDocstring:[39m
Cluster cells into subgroups :cite:p:`Traag2019`.

Cluster cells using the Leiden algorithm :cite:p:`Traag2019`,
an improved version of the Louvain algorithm :cite:p:`Blondel2008`.
It was proposed for single-cell analysis by :cite:t:`Levine2015`.

This requires having run :func:`~scanpy.pp.neighbors` or
:func:`~scanpy.external.pp.bbknn` first.

.. array-support:: tl.leiden

Parameters
----------
adata : 'AnnData'
    The annotated data matrix.
resolution : 'float', optional (default: 1)
    A parameter value controlling the coarseness of the clustering.
    Higher values lead to more clusters.
    Set to `None` if overriding `partition_type`
    to one that doesn’t accept a `resolution_parameter`.
random_state : '_LegacyRandom', optional (default: 0)
    Change the initialization of the optimization.
restrict_to : 'tuple[str, Sequence[str]] | None', optional (default: None)
    Restrict the clustering to the categories within the key for sample
    annotation, tuple needs to contain `(obs_key, list_of_categories)`.
key_added : 'str', optional (default: 'leiden')
    `adata.obs` key under which to add the cluster labels.
adjacency : 'CSBase | None', optional (default: None)
    Sparse adjacency matrix of the graph, defaults to neighbors connectivities.
directed : 'bool | None', optional (default: None)
    Whether to treat the graph as directed or undirected.
use_weights : 'bool', optional (default: True)
    If `True`, edge weights from the graph are used in the computation
    (placing more emphasis on stronger edges).
n_iterations : 'int', optional (default: -1)
    How many iterations of the Leiden clustering algorithm to perform.
    Positive values above 2 define the total number of iterations to perform,
    -1 has the algorithm run until it reaches its optimal clustering.
    2 is faster and the default for underlying packages.
partition_type : 'type[MutableVertexPartition] | None', optional (default: None)
    Type of partition to use.
    Defaults to :class:`~leidenalg.RBConfigurationVertexPartition`.
    For the available options, consult the documentation for
    :func:`~leidenalg.find_partition`.
neighbors_key : 'str | None', optional (default: None)
    Use neighbors connectivities as adjacency.
    If not specified, leiden looks at .obsp['connectivities'] for connectivities
    (default storage place for pp.neighbors).
    If specified, leiden looks at
    .obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.
obsp : 'str | None', optional (default: None)
    Use .obsp[obsp] as adjacency. You can't specify both
    `obsp` and `neighbors_key` at the same time.
copy : 'bool', optional (default: False)
    Whether to copy `adata` or modify it inplace.
flavor : "Literal['leidenalg', 'igraph'] | None", optional (default: None)
    Which package's implementation to use.
**clustering_args
    Any further arguments to pass to :func:`~leidenalg.find_partition` (which in turn passes arguments to the `partition_type`)
    or :meth:`igraph.Graph.community_leiden` from `igraph`.

Returns
-------
Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

`adata.obs['leiden' | key_added]` : :class:`pandas.Series` (dtype ``category``)
    Array of dim (number of samples) that stores the subgroup id
    (``'0'``, ``'1'``, ...) for each cell.

`adata.uns['leiden' | key_added]['params']` : :class:`dict`
    A dict with the values for the parameters `resolution`, `random_state`,
    and `n_iterations`.

`adata.uns['leiden' | key_added]['modularity']` : :class:`float`
    The modularity score of the final clustering,
    as calculated by the `flavor`.
    Use :func:`scanpy.metrics.modularity`\ `(adata, mode='calculate' | 'update')`
    to calculate a score independent of `flavor`.
[31mFile:[39m      ~/miniconda3/envs/rapids_singlecell/lib/python3.13/site-packages/scanpy/tools/_leiden.py
[31mType:[39m      function
