The testing dataset is a small mouse brain single-cell multiome dataset (221.2 MB) and can be downloaded from 10X Genomics website (https://www.10xgenomics.com/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-1-0-0). Key metrics of the dataset can be found on the website.

We will use the processed data for testing the pipeline. nformation about cell barcodes, features (gene names and locations, peak locations) can be found in the folder named "filtered_feature_bc_matrix". Specifically, 
1) snATAC fragments that passed quality contorl can be found in the file named e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz that comprises 4880 cells. IQuality fragments will be one of the input.
2) To study genes of interest, there are a total of 32245 different genes present in this dataset, which can be found in the file named "features.tsv.gz".
3) ATAC peak locations are stored in a file named "peaks.bed" that contains 144437 peaks and is provided on the same website also.

The first part of the pipeline will need 1) and 3) to calculate co-accessibility scores between any pair of peaks (enhancers), while the second part of the pipeline will take the output from the first part to quantify TEs of these enhancers.  

Another input file for the secod part of the pipeline is a text file of mouse TE consesnsus sequences that contains 1369 different TE families and their genomic coordinates. The TE annotation file can be downloaded from Dfam (https://www.dfam.org/releases/Dfam_3.8/annotations/mm10/) will be another input file.

