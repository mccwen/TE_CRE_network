**1. Introduction:**

   *Purpose:* The objecitve of this tool is to streamline the quantification of TE involvement in open chromatins, specifically at single-cell resolution.
   
**2. Overall Description:**
   
   *2.1 Scope:* This tool takes one fragment text file (e.g., e18_mouse_brain_fresh_5k_atac_fragments.tsv) which can contain one sample's or multiple samples' fragments and one TE consesus annotaiton file as input. Users should select a species-specific TE consensus annotation file. For both input files, the following minimum information should be included: chromoson, start and end positons
  
   *2.2 Product functions:* To determine transposable element (TE) invovlvement in the enhancer networks of interest at single-cell resolution
   
   *2.3 User characteristics:* Anyone interested in characterizing cell type-specific TE-CREs
   
**3. Requirements:**
   
   *3.1 Functional:* Run steamer.py to get a cell-by-TE matrix
   
   *3.2 Useability:* Terminal to run the python module
   
