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