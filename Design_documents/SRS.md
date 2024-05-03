**1. Introduction:**

   *Purpose:* The objecitve of this tool is to streamline the quantification of TE involvement in significnatly correlated enhancers of a target gene.
   
**2. Overall Description:**
   
   *2.1 Scope:* This tool first takes one fragment text file (e.g., e18_mouse_brain_fresh_5k_atac_fragments.tsv) which can contain one sample's fragment file, and then calculate co-accessibility scores of enhancers of a target gene and also highly correalted enhancers of the target gene. Significant enhancers are defined if they are highly correlated with the expression of the target gene and also highly correlated with each other.

   For signigifncat enhaners, we can then quantify TEs that intersect with these enahcers. 
        So the second step is to use the significant enhancer file and a TE consesus annotaiton file as input. Note that users should select a species-specific TE consensus annotation file. For both input files, the following minimum information should be included: chromoson, start and end positons.
  
   *2.2 Product functions:* To determine transposable element (TE) invovlvement in the enhancers of a given gene of interest
   
   *2.3 User characteristics:* Anyone interested in characterizing cell type-specific TE-CREs. However, users should preprocess ATAC-seq data to ensure the fitness of the input fragment file prior to running the current tool.
   
**3. Requirements:**
   
   *3.1 Functional:* The required packages should be installed first and run first cicero_monocle3.R script. Also, another R function named "gene_enhancer_corr", shoudl be run to identify enahcers that are highly correlated with the target gene. Afterwards, steamer.py should be run to get a cell-by-TE matrix.
   
   *3.2 Useability:* Terminal to run the python module
   
