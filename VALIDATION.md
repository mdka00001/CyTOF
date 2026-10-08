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
