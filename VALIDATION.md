# Validation

Validated in this workspace using `/usr/bin/Rscript` (R 4.5.1):

```sh
/usr/bin/Rscript tests/run.R
```

Passed:

- CLI parsing and invalid parameter errors.
- ROI/sample mapping, cell ID joins, shuffled region-property rows and unmapped-ROI rejection.
- Arcsinh/z-scoring, constant markers, rlog input and resource guards.
- Area filtering, coverage, physical cell density and SNR edge cases.
- Chunked compressed CSV round-trip and HTML/PDF generation, including a one-value heatmap.
- Pipeline integration with image and embedding modules explicitly disabled.
- UMAP and t-SNE pipeline integration for both default transformations.
- Raw integer TIFF label preservation, floating-point TIFF images, Otsu pixel SNR and segmentation previews.
- Actual DESeq2 rlog and rlog followed by z-scoring.

The synthetic fixture uses two ROIs and 60 cells; a separate integer fixture exercises rlog. Tests remove their temporary output directories. No real steinbock dataset was supplied, so large-dataset performance and biological interpretation have not been validated. GitHub Actions is configured but has not been run remotely.

The default PATH's Miniconda R 4.2.3 passed base-R tests but lacks analysis dependencies. Use `/usr/bin/Rscript` here. A temporary dependency installation was attempted under `/tmp/cytof-r-lib`; the repository does not depend on that library.

## Reference-script plotting regression

The revised report uses `dittoHeatmap`, `multi_dittoPlot`, `dittoDimPlot`,
`scater::runUMAP/runTSNE`, `cytomapper::plotPixels`, and patchwork grids.
The integration suite passed under `/usr/bin/Rscript`, including a regression
with four input-ordered patients: all 60 test cells and all four patient colors
are retained in embedding panels, and marker panels use continuous color scales.
All-marker distributions and embeddings use one grid page per plot type and
transform. Poppler `pdfunite` preserves each grid's page dimensions in the final PDF.

The revised pipeline was also run on the supplied synovium dataset: 42 markers,
27 ROIs, 50,655 input cells and 50,145 retained cells. Both transformed embeddings
contain all 50,145 retained cells and all four patients. The regenerated
`results_revised/report.pdf` has 104 pages (the prior report had 612). All-marker
boxplot/ridgeline and UMAP/t-SNE grids were checked in the PDF, image channel
colorbars and heatmap/embedding legends were visually inspected, and all output
manifest checksums verified. The original results directory was preserved.
