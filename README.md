# TE_CRE_network

1. Steamer
The core of the this pipeline is _steamer_ which was developed by the Welch lab at University of Michigan and is used to qnautification of TE counts per cell.

The example here is to use a mouse TE annotation database which can be downloaded from here: https://www.dfam.org/releases/Dfam_3.8/annotations/mm10/. 

The pipeline, _steamer_,  takes two input files: 1) fragmemt.tsv from single nucleus ATAC-seq and 2) annotated TE sequences in any text format. It outputs a cell-by-TE count matrix.

2. Cicero-monocle3: This R script is to calcuate co-accessibility scores betwween enhancers of a target gene. The example output from this script contains co-accessibility score between any two peaks and can be found in the example_data folder. 


