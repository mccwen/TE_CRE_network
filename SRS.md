**1. Introduction:**

*1.1 Purpose:* The objecitve of this tool is to streamline the quantification of TE involvement in open chromatins.

<to determine transposable element (TE) invovlvement in the enhancer networks of interest at single-cell resolution>
   
<*1.2 Product Scope:* This tool is designed to quantify TE activity in target open chromatins using the single-nucleus assay for transposase-accessible chromatin with sequencing (snATAC-Seq).
   Although it can also be used for bulk-seq data, the tool development is primary with cell type specific TE-enhancer relationships in mind.>
   NB: Type of input files, number of input files, other parameters/ config
   
<*1.3 References:* Interested readers can check out the following publications on TEs or cell type-specific TEs:
   
      i) Storer J., et al. The Dfam community resource of transposable element families, sequence models, and genome annotations. Mobile DNA 12, 2 (2021). https://doi.org/10.1186/s13100-020-00230-y.
   
      ii) He, J., et al. Identifying transposable element expression dynamics and heterogeneity during development at the single-cell level with a processing pipeline scTE. Nat Commun 2021;12(1):1456. doi: 10.1038/s41467-021-21808-x.
   
      iii) Zu, S., et al. Single-cell analysis of chromatin accessibility in the adult mouse brain. Nature 2023; 624(7991):378-389. doi: 10.1038/s41586-023-06824-9
   
**2. Overall Description:**
   
   <*2.1 Product perspective:* To >
   *2.1 Scope: * This tool takes one fragment file which can contain one sample's or multiple samples' fragments and one TE consesus annotaiton file as input. Users can select the species-specific TE consensus annotation file. 
  < This tool will take fragments of target enhancers and a species-specific TE concesus, such as mm10 for mouse data from Dfam (see reference 1 above) as input, and output a cell-by-TE count matrix.
      i) System Interfaces ???
      ii) User interfaces: None
      iii) Hardware interfaces: CPU
      iv) Software interfaces: Python (current development) and R (for future add-on functions)
      v) Memory constraints: >= 84 GB>
   
   ~~(*2.2 Design constraints*
        Operations
        Site adaptation requirements)~~
   
   *2.3 Product functions:* To determine transposable element (TE) invovlvement in the enhancer networks of interest at single-cell resolution
   
   *2.4 User characteristics:* Anyone interested in characterizing cell-specific TE-CREs
   
   <2.5 Constraints, assumptions and dependencies: The tool uses the default parameter of bedtools (i.e., a minimum of 1 base pair overlap between a TE and a CRE seruqneces) to define TE-deriven  CRE. Nevertheless, users can change this as they see fit. To do so, please consult bedtools documentation (https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) for detailed instructions>
   
**3. Requirements:**
   
   *3.1 Functional:* Run steamer.py to get a cell-by-TE matricx
   <the input file should be in bed or tsv format and include chromosome, start, and end columns, ii) to run the tool using Python (version >= 3.4) packages, including pandas, scipy, numpy, and pybedtools, to construct functions for cell-by-TE count matrices>   
   
   *3.2 Useability:* Terminal to run the python module
   
   <*3.3 Interfaces:* Python modules or jupyternotebook>
