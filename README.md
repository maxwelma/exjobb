# exjobb

Github repo for master thesis project "Single-cell RNA-seq mapping of chicken leukocytes"

Scripts:
├──bash_scripts  
│	├──blast.sh  - blast all reads against a gene to find matches
│	├──(bwa.sh) 
│	├──count.sh  - run cellranger count
│	├──doublet_removal.sh  - activate conda and run python script
│	├──make_ref.sh – make celllranger reference  
│	└──setup_conda.sh – set up conda environment for doublet detection
├──python_scripts  
│	├──DoubletDetection  - python package files
│	├──doublet_removal_env.yml – environment for conda
│	└──doublet_removal.py – perform the doublet estimation and removal
└──R_scripts
	├──appendix.R  - produce appendix plots
	├──combined_data_sct.RDS  - processed Seurat object
	├──figures.R  - produce figures for rapport
	├──functions.R  - functions for analyses
	├──markers.R  - marker lists for figure plotting
	├──produce_excels.R  - produce lists of DE and GOs
	├──quality_control.R  - produce plots for quality control
	└──runfunctions.R – process data

Script order to produce processed data:
1.	make_ref.sh 
  a.	make cellranger reference from annotation and reference genome 
2.	count.sh
  a.	Run the counting algorithm on each sample using the generated reference, give sample id as input.
3.	doublet_removal.sh  
  a.	Run the python doublet_removal.scriptipt on each sample using sample id as input. Takes matrix files as input and gives edited matrix files out.
4.	Runfunctions.R
  a.	Give path to folders including matrix.mtx, features.tsv, and barcodes.tsv files. Call functions from functions.R to process data. Save processed data as .RDS object. 

