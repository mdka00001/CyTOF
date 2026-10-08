# CyTOF QC

A modular R command-line workflow for **steinbock-preprocessed IMC data**. Produces explained HTML and PDF QC reports, raw and transformed matrices, marker review tables, UMAP/t-SNE coordinates, and reusable R objects. The original two scripts are preserved as reference material.

## Quick start

Use R >= 4.3 with a compatible Bioconductor release. Install dependencies once:

```sh
Rscript scripts/install.R
Rscript bin/cytof-qc.R \
  --input /path/to/steinbock \
  --samples /path/to/samples.csv \
  --output results/run_01 \
  --masks masks_deepcell \
  --image-channels CD3,CD20,DNA1
```

On the current workstation, use `/usr/bin/Rscript` for these commands: it has the required packages. The `Rscript` on the default PATH points to a separate Miniconda R installation.

Use an empty/new output directory. `Rscript bin/cytof-qc.R --help` lists all options and defaults. Paths are resolved relative to your working directory; the executable itself can be invoked from anywhere. Installation may need system development libraries for TIFF, FFTW and the Bioconductor dependencies; installation errors identify missing libraries. Neither Pandoc nor LaTeX is needed.

Dependencies: `mclust` for cell SNR; `uwot` and `Rtsne` for embeddings; `EBImage` and `tiff` for image QC; `DESeq2` for optional rlog. `--skip-images` and `--skip-embeddings` explicitly disable those modules and their dependency checks. Skipped modules are identified in the report.

## Install the `cytof-qc` command (TASK 2)

Install a standalone application copy and launcher once, then run `cytof-qc`
from any working directory. On this workstation:

```sh
bash scripts/install-cli.sh --rscript /usr/bin/Rscript \
  --shell-rc "$HOME/.bashrc" --shell-rc "$HOME/.profile"
export PATH="$HOME/.local/bin:$PATH"
cytof-qc --help
cytof-qc --input /path/to/steinbock --samples samples.csv \
  --output results/run_01 --masks masks_deepcell
```

The installer defaults to `~/.local/bin/cytof-qc` and copies the application to
`~/.local/lib/cytof-qc`. No sudo is required. `--shell-rc` adds a PATH entry to
an explicitly chosen startup file without duplicating it on reinstall. Bash
interactive sessions read `.bashrc`; login shells typically read `.profile`
(unless `.bash_profile` or `.bash_login` takes precedence). Use the appropriate
file for your shell, e.g. `--shell-rc "$HOME/.zshrc"` for zsh. Without
`--shell-rc`, no startup files are edited. The printed `export PATH=...` command
activates the installation in an already-open terminal; new terminals load it
from the configured startup file.

On another machine, omit `--rscript` to select the Rscript currently on PATH,
or specify its absolute path. The launcher remembers the selected executable,
so later Conda activation does not silently switch R libraries. R packages
are installed separately with that same executable and `scripts/install.R`.
The installer checks that the copied CLI starts, but does not install or
validate analysis packages; the QC command checks the packages it needs.

Custom location:

```sh
bash scripts/install-cli.sh --prefix /path/to/prefix --rscript /path/to/Rscript
export PATH="/path/to/prefix/bin:$PATH"
```

Run the installer again after changing branches or updating the repository to
refresh the installed copy. The installed command does not depend on this
checkout and continues to work if it is moved. Reinstallation replaces only
paths marked as managed by this installer; unrelated existing commands or
application directories are rejected. Input/output paths remain relative to
the directory where you run `cytof-qc`.

To uninstall the default installation, remove `~/.local/bin/cytof-qc` and
`~/.local/lib/cytof-qc`, then remove lines ending in `# cytof-qc PATH` from the
startup files you selected. This does not remove R packages or QC results.

Installer tests: `bash tests/install-cli.sh` (override the test R executable
with `CYTOF_TEST_RSCRIPT=/usr/bin/Rscript`).

## Input contract

```text
steinbock/
  panel.csv
  images.csv
  intensities/ROI_ID.csv
  regionprops/ROI_ID.csv
  img/ROI_ID.tiff
  masks/ROI_ID.tiff
```

- `panel.csv`: `channel,name`, optionally `keep` and `use_channel` (0/1 or TRUE/FALSE). `keep=0` rows are removed before channel matching. Retained names and channel IDs must be unique. TIFF pages follow retained panel order. `use_channel` controls heatmaps, distribution plots and embedding features; all measured markers remain exported and get SNR review and embedding overlays. If absent, all retained markers are selected. `--exclude CD3,CD20` additionally excludes named markers from those feature selections.
- `images.csv`: `image,width_px,height_px`, one row per ROI. ROI IDs are image basenames without their extension; no heuristic extraction of trailing numbers is used.
- Intensity tables: `Object` and columns matching either panel marker names or channel IDs. Values must be finite and nonnegative. These are typically mean intensities, even though the output calls the raw assay `counts`.
- Region-property tables: `Object,area`, plus any other consistent measurement columns such as `centroid-0,centroid-1`. Object IDs must be unique within an ROI and match the intensity table exactly. The join uses IDs, not row order. Metadata columns must have a consistent schema across ROIs.
- TIFF images: grayscale pages, one per retained marker. Supported TIFF storage is 8/16-bit integer or 32-bit floating-point; masks should use steinbock's 16-bit integer format. Masks: single grayscale page, background 0, positive integer cell IDs matching `Object`. TIFF reads preserve raw values. Mismatched image/mask dimensions or labels fail explicitly; there is no silent resizing.
- Override subdirectory names with `--intensities`, `--regionprops`, `--images`, `--masks`. Data from different segmentation runs must use matching intensity, region-property and mask directories.

Mapping CSV or tab-delimited TXT requires a header:

```csv
roi_id,sample_name,patient_id,batch
sample_001,Sample_A,Patient_A,batch_1
sample_002,Sample_A,Patient_A,batch_1
sample_003,Sample_B,Patient_B,batch_2
```

Every image ROI needs exactly one mapping row; multiple ROIs may share a sample. Optional metadata columns are retained. Duplicate and missing IDs fail; unused mapping rows emit a warning. See [examples/samples.csv](examples/samples.csv).

## Transformations and filtering

Default: `--transforms arcsinh,arcsinh_zscore`. Arcsinh is `asinh(x/cofactor)`, default cofactor 1. Z-scoring is per marker across **all retained cells**; constant markers become zero. Raw matrices are never modified. Each transform is processed and saved separately.

Cell areas below `--min-area 5` are excluded from transforms and embeddings. Set `--min-area 0` to retain all cells. Size distributions, coverage and cell SNR use pre-filter cells; density uses post-filter cells. `--pixel-size 1` is the pixel side length in micrometers; change it to match the acquisition. Density assumes square pixels.

For true DESeq2 rlog and rlog followed by marker-wise z-scoring:

```sh
Rscript bin/cytof-qc.R --input /path/to/integer-count-data \
  --samples samples.csv --output results/rlog \
  --transforms arcsinh,arcsinh_zscore,rlog,rlog_zscore
```

**Rlog requires nonnegative integer counts, not typical fractional mean intensities.** The program rejects fractional data instead of rounding it or substituting log1p. It uses DESeq2 `rlog(blind=TRUE, fitType="mean")`, with `poscounts` size factors. At least three cells, no all-zero cells, and a successful dispersion fit are required. Small marker panels or degenerate counts can fail DESeq2 estimation; such failures stop the run, without a replacement transformation. Rlog was designed for sequencing counts; integer-valued IMC measurements do not establish that its model assumptions hold.

Rlog is expensive across many cells. `--rlog-max-cells 1000` is a hard guard, not a subsampling command; explicitly increase it only when resources permit. The regularized log is not claimed to have linear complexity. Arcsinh is the practical default for large IMC datasets.

## QC and interpretation

| Output | What to inspect |
| --- | --- |
| Sampled cell heatmaps | Marker coexpression and staining specificity; clustering is restricted to the plotted subset. |
| ROI mean heatmaps | Mean **transformed** expression by ROI; outliers can be biological or technical. |
| RGB image previews with white mask outlines | Localization, merged/missed cells, debris and segmentation alignment. Up to three channels; display-only 99th-percentile clipping per image. |
| Pixel SNR, unfiltered and signal-filtered | Otsu signal/background ratio per image-marker; compare positive signal intensity as well as ratio. |
| Coverage, counts and density | Sparse ROIs, filtering impact, and potential segmentation problems. |
| Cell SNR | Two-component Gaussian mixture on arcsinh expression; raw positive/negative means define SNR. |
| Cell area distributions | Segmentation or biological size differences; filter threshold shown. |
| Sample expression boxplots | Staining shifts and population changes. Boxplots replace the reference scripts' ridgelines with a lightweight base-R distribution summary. |
| UMAP/t-SNE by metadata and all markers | Exploratory mixing/separation, batch patterns and marker localization. |

Each plot has an explanation in both reports. SNR suggestions use configurable `--low-snr 2`, `--high-snr 10`, and `--min-signal 2`. High ratios alone do not demonstrate antibody specificity: low signal or zero background can inflate them. Undefined/constant/failed-mixture markers are explicitly flagged. Infinite SNR remains in CSV, but is excluded from log scatterplots. Potential exclusions are suggestions for inspection, never automatic removal. Thresholds depend on signal units and are not universal biological cutoffs.

## Outputs

| Path | Contents |
| --- | --- |
| `report.html`, `plots/` | HTML with linked PNG figures; distribute them together. |
| `report.pdf` | Multipage PDF with all figures and explanation pages. |
| `matrices/counts_all.csv.gz` | All input cells, before filtering. |
| `matrices/counts_retained.csv.gz` | Raw retained cells. |
| `matrices/<transform>.csv.gz` | Transformed retained cells, one cell per row; first column `cell_id`. |
| `coordinates/<transform>_<UMAP\|TSNE>.csv` | `cell_id,roi_id,sample_name` and two coordinates; sampled cells only. |
| `tables/cells_all.csv` | Metadata, region properties, stable IDs and `kept` filter status. |
| `tables/roi_qc.csv` | Coverage, before/after counts and density, including empty post-filter ROIs. |
| `tables/*snr*.csv` | Pixel/cell SNR statistics and cell-based marker review. |
| `tables/<transform>_roi_means.csv` | Marker-by-ROI mean transformed expression. |
| `objects/qc_data.rds` | Named list of retained counts, cell metadata, panel, selected markers and QC tables. |
| `objects/<transform>.rds` | Marker-by-cell transformed matrix with dimnames. |
| `objects/<transform>_embeddings.rds` | Coordinates, cell IDs, features and effective neighbor/perplexity settings. |
| `config.txt`, `objects/config.rds`, `sessionInfo.txt` | Parameters and package versions. |
| `manifest.csv`, `STATUS` | File checksums/sizes and `complete` or failed-run state. |

IDs use `ROI_ID::Object`. Matrices and coordinates join through `cell_id`; R matrices have marker rows and cell columns, while matrix CSVs have cell rows and marker columns. The RDS data object is an ordinary R list (not SpatialExperiment). Spatial centroids are preserved as input metadata; neighborhood graphs are not imported or needed for these QC steps. Raw image objects are not duplicated on disk.

```r
qc <- readRDS("results/run_01/objects/qc_data.rds")
x <- readRDS("results/run_01/objects/arcsinh.rds")
umap <- read.csv("results/run_01/coordinates/arcsinh_UMAP.csv")
stopifnot(identical(colnames(x), qc$cells$cell_id))
```

## Resource use and reproducibility

For M markers and N cells, matrix operations and exports are O(MN), with O(MN) memory plus temporary copies; this is not an out-of-core assay engine. Images stream one ROI at a time, adding memory proportional to the largest ROI, rather than the whole image collection. CSV export transposes chunks of 1,000 cells. Arcsinh/z-score transforms are O(MN); DESeq2 rlog can be substantially more expensive.

Cell SNR uses at most `--snr-cells 30000`. Embeddings use a seeded subset of at most `--embedding-cells 30000`, PCA to at most 30 components, approximate UMAP neighbors and Barnes-Hut t-SNE. Coordinates are only returned for those cells; there is no fabricated projection for unsampled cells. Increase the cap to the retained cell count for full coordinates. Perplexity and neighbor counts adapt to small inputs; fewer than five cells or two variable selected markers require `--skip-embeddings`.

Plots use at most `--plot-cells 5000`; clustered cell heatmaps use `--heatmap-cells 500`. Heatmap clustering is quadratic in that cap, and ROI clustering is quadratic in the ROI count. Raising these limits increases costs. Embedding plot cells are the first bounded entries of the seeded embedding subset. Sampling is uniform, so very rare populations may be absent. `--seed` defaults to 220225; embeddings use one thread. Reproducibility also depends on package/platform versions recorded per run. Package versions are not locked.

## Development and checks

Modules are separated into CLI/configuration, IO, transformations, QC, streaming image processing, reporting and orchestration under `R/`. Run:

```sh
Rscript tests/run.R
```

Base-R tests cover joins, validation, zero-variance transforms, rlog guards, density/coverage, chunked matrix export, and HTML/PDF generation. Integration checks run when their dependencies are installed and otherwise print explicit skips. GitHub Actions installs analysis dependencies and runs the suite. No patient data is bundled.

## Sources

- [Bodenmiller IMC image and cell quality control workflow](https://bodenmillergroup.github.io/IMCDataAnalysis/image-and-cell-level-quality-control.html)
- [steinbock file formats](https://bodenmillergroup.github.io/steinbock/latest/file-types/) and [measurements](https://bodenmillergroup.github.io/steinbock/latest/cli/measurement/)
- [DESeq2 transformation methodology and assumptions](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
