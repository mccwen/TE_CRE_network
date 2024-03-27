#This module was primarily developed by and taken from the Welch lab with some modifications

import pandas as pd
import numpy as np
from pybedtools import BedTool
import scipy.sparse
from fuc import pybed
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix
import os
import glob


def create_bed_for_TEs(filename): 
    """
    Takes in a file (ideally a file containing TEs) and converts it into a BED file.
    
    Args:
        filename: str, filepath to the input file (ideally .rph.hits)
        
    Returns:
        bf: pybed.BedFrame object (and creating the TE bed file)
    """
    
    # For example, col_names can be the mm10.nrph.hits file which is mouse TE annotations
    col_names = ["seq_name","ali-st","ali-en"]
    other_col = ['family_name','strand']
    #nrows set to 1000 for testing
    #Below we only take columns as specified in col_names in the input file and thus columns 0, 9, and 10 are extracted
    TE_df = pd.read_csv(filename,sep='\t',usecols=[0,9,10], names=col_names)
    TE_df=TE_df.loc[1:]
    
    # Similarly, we will need "family_name" and "strand" later; both are in columns 2 & 8 in the input file
    strand_name_df =pd.read_csv(filename,sep='\t',usecols=[2,8], names= other_col)
    strand_name_df=strand_name_df.loc[1:]
     
    #next we need to sort the strands out or else it will cause issues later on,i.e. the start strand cannot be greaer than the end
    TE_df["ali-st"] = TE_df["ali-st"].astype(int, errors='ignore')
    TE_df["ali-en"] = TE_df["ali-en"].astype(int, errors='ignore')
    for row in TE_df.itertuples():
        row_index,seq_name,start,end = row
        
        if start > end:
            TE_df.at[row_index,'ali-st'] = end
            TE_df.at[row_index,'ali-en'] = start
    #renaming columns to prepare for pybed.BedFrame.from_frame
    TE_df.rename(columns = {'seq_name':'Chromosome','ali-st':'Start','ali-en':'End'},inplace = True)
    TE_df = pd.concat([TE_df,strand_name_df],axis=1)
    bf = pybed.BedFrame.from_frame(meta = [],data = TE_df)
    bf=bf.sort()
    TE_bed_file = bf.to_file('TEs.bed')
    return bf


def create_bed_for_enhancers(filename):
    """
    Takes in a file (ideally a fragment file) and converts it into a BED file.
    
    Args:
        filename: str, filepath to the input file 
        
    Returns:
        bf: pybed.BedFrame object (and creating the Fragment bed file)
    """
    
    frag_df = pd.read_csv(filename)
    cell_barcode=frag_df["barcode"]
    
    bf = pybed.BedFrame.from_frame(meta=[], data=frag_df)
    bf=bf.sort()
    bf.to_file("Enhancer.bed")
    return bf


def intersection(TE_bed, frag_bed):
    """
    This function takes a TE annotaiton file and a fragment file (both in .bed format) and use bedtools 
    to intersect them for finding overlaps. 

    The output file is another bed file.

    """
    bed_intersect=TE_bed.intersect(frags_bed, wb=True, sorted=True)
    return bed_intersect


def make_cell_x_element_matrix(bed_intersect, cell_barcodes):
    # Initialize lists and dictionaries
    row_indices = []  # Row indices for sparse matrix
    col_indices = []  # Column indices for sparse matrix (corresponding to UniqueTE)
    fams_col_indices = []  # Column indices for sparse matrix (corresponding to TE_Fam)
    data = []  # Data values for sparse matrix
    unique_TEs_list = defaultdict(lambda: len(unique_TEs_list))  # Dictionary to map UniqueTE names to indices
    TEs_fam_dict = defaultdict(lambda: len(TEs_fam_dict))  # Dictionary to map TE_Fam names to indices
    barcode_dict = {barcode: i for i, barcode in enumerate(cell_barcodes)}  # Dictionary to map barcodes to indices
    
    # Iterate over bed_intersect
    for interval in bed_intersect:
        fam_name = interval.fields[3]
        barcode = interval.fields[8]
        TE_start, TE_end, chrom = interval.fields[1], interval.fields[2], interval.fields[0]
        TE_name_unique = f"{fam_name}({chrom}:{TE_start},{TE_end})"
        # Check if UniqueTE is already in unique_TEs_list, otherwise assign a new index
        UniqueTE_index = unique_TEs_list[TE_name_unique]
        # Check if TE_Fam is already in TEs_fam_dict, otherwise assign a new index
        TEFam_index = TEs_fam_dict[fam_name]
        # Get the barcode index from barcode_dict
        barcode_index = barcode_dict[barcode]
        # Append values to the respective lists
        row_indices.append(barcode_index)
        col_indices.append(UniqueTE_index)
        fams_col_indices.append(TEFam_index)
        data.append(1)
    # Create sparse matrices using coo_matrix
    UniqueTEs_sparse_matrix = coo_matrix((data, (row_indices, col_indices)))
    TE_Fams_sparse_matrix = coo_matrix((data, (row_indices, fams_col_indices)))
    # Convert sparse matrices to DataFrames and perform grouping
    UniqueTEs_sparse_df = pd.DataFrame({'barcode_index': UniqueTEs_sparse_matrix.row,
                                        'UniqueTE_index': UniqueTEs_sparse_matrix.col,
                                        'data': UniqueTEs_sparse_matrix.data})
    UniqueTEs_sparse_df = UniqueTEs_sparse_df.groupby(['barcode_index', 'UniqueTE_index'], as_index=False)['data'].sum()
    TE_Fams_sparse_df = pd.DataFrame({'barcode_index': TE_Fams_sparse_matrix.row,
                                      'FamTE_index': TE_Fams_sparse_matrix.col,
                                      'data': TE_Fams_sparse_matrix.data})
    TE_Fams_sparse_df = TE_Fams_sparse_df.groupby(['barcode_index', 'FamTE_index'], as_index=False)['data'].sum()
    # Return the DataFrames and dictionaries
    return (UniqueTEs_sparse_df, TE_Fams_sparse_df, dict(unique_TEs_list), dict(TEs_fam_dict), barcode_dict)

