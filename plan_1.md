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
This task should implement methods and a separate pipeline execution arguement that takes multiple *.h5ad objects of choice and integrate them to a single object. The pipeline should run two batch correction methods defined by users:
    1. scvi integration
    2. harmony integration
Run the integration by batch ( here each object shall be considered as a separate batch). Output shall return:
    1. integrated h5ad file.
    2. UMAP+TSNE plot for batch, groups (patient/sample).
    3. Shall plot snr plot.
    4. shall plot heatmap+dotplot+matrixplot by batch+group+roi.
    5. shall store the scvi model.
    6. shall store both scvi/harmony embeddings.
Use the same abstract object and same code base skeleton for the implementation.