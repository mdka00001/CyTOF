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

Use an empty/new output directory. `Rscript bin/cytof-qc.R --help` lists all options and defaults. Paths are resolved relative to your working directory; the executable itself can be invoked from anywhere. Installation may need system development libraries for TIFF, FFTW and the Bioconductor dependencies; installation errors identify missing libraries. Neither Pandoc nor LaTeX is needed. Install Poppler utilities (`sudo apt-get install poppler-utils` on Ubuntu) for `pdfunite`, which assembles report pages of different sizes.

Dependencies: `mclust` for cell SNR; `scater` with `uwot` and `Rtsne` for embeddings; `dittoSeq`, `SingleCellExperiment`, `patchwork`, `ggrastr`, `ggplot2`, `viridis` and `RColorBrewer` for the reference-script plots; `EBImage`, `tiff` and `cytomapper` for image QC; `DESeq2` for optional rlog. `--skip-images` and `--skip-embeddings` explicitly disable those modules and their dependency checks. Skipped modules are identified in the report.

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

- `panel.csv`: `channel,name`, optionally `keep` and `use_channel` (0/1 or TRUE/FALSE). `keep=0` rows are removed before channel matching. Leading/trailing whitespace is trimmed from panel names/channel IDs and intensity headers; ambiguous duplicate names after trimming are rejected. Retained names and channel IDs must be unique. TIFF pages follow retained panel order. `use_channel` controls heatmaps, distribution plots and embedding features; all measured markers remain exported and get SNR review and embedding overlays. If absent, all retained markers are selected. `--exclude CD3,CD20` additionally excludes named markers from those feature selections.
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
| Sampled cell heatmaps | `dittoHeatmap` with patient annotations and the reference blue-white-red scale (-4 to 4, saturating outside that range); clustering is restricted to the plotted subset. |
| ROI mean heatmaps | `dittoHeatmap` of mean **transformed** expression by ROI, viridis legend, and patient/ROI annotations; outliers can be biological or technical. |
| RGB image previews with white mask outlines | Localization, merged/missed cells, debris and segmentation alignment. Up to three channels, using the reference `cytomapper::plotPixels` approach with its channel color legend and display normalization. |
| Pixel SNR, unfiltered and signal-filtered | Otsu signal/background ratio per image-marker; compare positive signal intensity as well as ratio. |
| Coverage, counts and density | Sparse ROIs, filtering impact, and potential segmentation problems. |
| Cell SNR | Two-component Gaussian mixture on arcsinh expression; raw positive/negative means define SNR. |
| Cell area distributions | Segmentation or biological size differences; filter threshold shown. |
| Patient expression distributions | `multi_dittoPlot` ridgelines, as in the reference script, plus the requested grouped boxplots. Each plot type contains all selected markers on one page per transform; all retained cells are used. |
| UMAP/t-SNE by metadata and all markers | `scater::runUMAP/runTSNE` and `dittoDimPlot`, as in the reference scripts. Paired metadata panels share a page. All marker UMAP panels share one page; all marker t-SNE panels share another, with individual continuous viridis colorbars. |

Each plot has an explanation in both reports. Marker grids use four columns and 5-by-4-inch panels, as in the reference script, so a 42-marker grid is a single 20-by-44-inch PDF page. Zoom to read individual panels; HTML figures link to full-size images. Patient grouping uses `patient_id` when supplied, otherwise `sample_name`. SNR suggestions use configurable `--low-snr 2`, `--high-snr 10`, and `--min-signal 2`. High ratios alone do not demonstrate antibody specificity: low signal or zero background can inflate them. Undefined/constant/failed-mixture markers are explicitly flagged. Infinite SNR remains in CSV, but is excluded from log scatterplots. Potential exclusions are suggestions for inspection, never automatic removal. Thresholds depend on signal units and are not universal biological cutoffs.

## Outputs

| Path | Contents |
| --- | --- |
| `report.html`, `plots/` | HTML with linked PNG figures; distribute them together. |
| `report.pdf` | Multipage PDF with all figures and explanation pages. |
| `matrices/counts_all.csv.gz` | All input cells, before filtering. |
| `matrices/counts_retained.csv.gz` | Raw retained cells. |
| `matrices/<transform>.csv.gz` | Transformed retained cells, one cell per row; first column `cell_id`. |
| `coordinates/<transform>_<UMAP\|TSNE>.csv` | `cell_id,roi_id,sample_name,patient_id` and two coordinates; all retained cells by default. |
| `tables/cells_all.csv` | Metadata, region properties, stable IDs and `kept` filter status. |
| `tables/roi_qc.csv` | Coverage, before/after counts and density, including empty post-filter ROIs. |
| `tables/<transform>_embedding_patient_counts.csv` | Numbers of actually embedded/plotted cells per patient. |
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

## FCS export (separate subcommand)

TASK 3 adds FCS construction as a separate command. First complete a regular QC
run, then point the FCS command to its output directory:

```sh
cytof-qc fcs \
  --qc-dir results/run_01 \
  --output results/run_01/fcs \
  --base arcsinh \
  --scales minmax,zscore \
  --area-min 5 --area-max 250 \
  --split-by time_group
```

The source assay is `objects/arcsinh.rds` by default; use `--base rlog` after a
QC run that produced `objects/rlog.rds`. For each selected assay, the command
writes a combined FCS and, when the `time_group` column exists, one FCS per
group. If that column is absent, it reports that fact and writes the combined
file. For example, the dataset mapping used in the supplied workflow has no
`time_group`, so it produces combined files unless that metadata is added during
QC. Use `--split-by none` to disable group files explicitly.

The default area filter includes cells with area from 5 through 250, matching
the supplied script. Filtering occurs before export normalization. `minmax`
rescales each marker across included cells to [0,1]; constant markers become
zero. `zscore` standardizes each marker across included cells; constant markers
become zero. These are generated from the selected QC assay and do not alter
QC matrices. Choose `--area-min` and `--area-max` to suit the dataset, or use
`--area-min 0 --area-max <large-value>` to retain all sizes.

FCS events contain every marker, numeric region property and coordinate, plus
numeric channels for metadata. FCS channels must be numeric; each FCS has a
matching `_events.csv` sidecar that preserves the original cell IDs and text
metadata alongside the event row index. The output also contains an
`included_cells.csv` and a brief `README.txt`. Install the Bioconductor
`flowCore` package with `Rscript scripts/install.R` before exporting. See
`cytof-qc fcs --help` for all options.

## Resource use and reproducibility

For M markers and N cells, matrix operations and exports are O(MN), with O(MN) memory plus temporary copies; this is not an out-of-core assay engine. Images stream one ROI at a time, adding memory proportional to the largest ROI, rather than the whole image collection. CSV export transposes chunks of 1,000 cells. Arcsinh/z-score transforms are O(MN); DESeq2 rlog can be substantially more expensive.

Cell SNR uses at most `--snr-cells 30000`. The default `--embedding-cells 0` means all retained cells, following the original scripts. A positive limit explicitly requests a seeded subset. `scater::runUMAP` and `scater::runTSNE` receive the selected variable markers and transformed assay directly, using their standard preprocessing rather than a separate custom PCA. Perplexity and neighbor counts adapt to small inputs; fewer than five cells or two variable selected markers require `--skip-embeddings`.

All embedded cells are plotted. The earlier report took the first 5,000 cells from sorted embedding IDs, which could omit patients appearing later in the input; this truncation has been removed. `--plot-cells` remains accepted for compatibility but does not truncate marker distributions or embedding plots. Patient-count tables expose the population used in each embedding. If a positive embedding cap is set, sampling may still miss rare populations; the default avoids this.

Clustered cell heatmaps use `--heatmap-cells 500` to bound quadratic clustering memory, while using the original `dittoHeatmap` plotting approach. ROI clustering is also quadratic in ROI count. `--seed` defaults to 220225 and embeddings use one thread. Reproducibility also depends on package/platform versions recorded per run; versions are not locked.

## Development and checks

Modules are separated into CLI/configuration, IO, transformations, QC, streaming image processing, reporting and orchestration under `R/`. Run:

```sh
Rscript tests/run.R
```

Base-R tests cover joins, validation, zero-variance transforms, rlog guards, density/coverage, chunked matrix export, and HTML/PDF generation. Regression tests also verify that four input-ordered patients remain present in the embedding panels and marker panels have continuous color scales. Integration checks run when their dependencies are installed and otherwise print explicit skips. GitHub Actions installs analysis dependencies and runs the suite. No patient data is bundled.

## Sources

- [Bodenmiller IMC image and cell quality control workflow](https://bodenmillergroup.github.io/IMCDataAnalysis/image-and-cell-level-quality-control.html)
- [steinbock file formats](https://bodenmillergroup.github.io/steinbock/latest/file-types/) and [measurements](https://bodenmillergroup.github.io/steinbock/latest/cli/measurement/)
- [DESeq2 transformation methodology and assumptions](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
