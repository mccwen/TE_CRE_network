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
    #TE_df["ali-st"] = TE_df["ali-st"].astype(int, errors='ignore')
    #TE_df["ali-en"] = TE_df["ali-en"].astype(int, errors='ignore')
    for row in TE_df.itertuples():
        row_index,seq_name,start,end = row
        #print(row_index,start,end)
        if start > end:
            TE_df.at[row_index,'ali-st'] = end
            TE_df.at[row_index,'ali-en'] = start
    #renaming columns to prepare for pybed.BedFrame.from_frame
    TE_df.rename(columns = {'seq_name':'Chromosome','ali-st':'Start','ali-en':'End'},inplace = True)
    TE_df = pd.concat([TE_df,strand_name_df],axis=1)
    bf = pybed.BedFrame.from_frame(meta = [],data = TE_df)
    TE_bed_file = bf.to_file('TEs.bed')
    return bf


def create_bed_for_fragments(filename,quality_barcode_file=''):
    """
    Takes in a file (ideally a fragment file) and converts it into a BED file.
    
    Args:
        filename: str, filepath to the input file (ideally fragment file: .tsv)
        quality_barcode_file: str, filepath to the input file(file containing only the barcodes that passed QC)
        
    Returns:
        bf: pybed.BedFrame object (and creating the Fragment bed file)
    """
    
    col_names = ["Chromosome", "Start", "End", "barcode"]
    frag_df = pd.read_csv(filename, sep="\t", usecols=[0, 1, 2, 3], names=col_names)
    #frag_df = pd.read_csv(filename, sep=",", usecols=[0, 1, 2, 3], names=col_names)
    frag_df=frag_df.loc[1:]
    #frag_df = pd.read_csv(filename, sep=",", usecols=[0, 1, 2, 3])
    if quality_barcode_file == '' :
        bf = pybed.BedFrame.from_frame(meta=[], data=frag_df)
        bf.to_file("Frag.bed")
    else:
        quality_barcodes_df = pd.read_csv(quality_barcode_file, sep='\t',names=['barcode'])
        quality_barcodes = quality_barcodes_df['barcode'].tolist()
        frag_df = frag_df[frag_df['barcode'].isin(quality_barcodes)]
        bf = pybed.BedFrame.from_frame(meta=[], data=frag_df)
        bf.to_file("Frag.bed")
    return bf


def intersection(TE_bed, frag_bed):
    """
    This function takes a TE annotaiton file and a fragment file (both in .bed format) and use bedtools 
    to intersect them for finding overlaps. 

    The output file is another bed file.

    """
    bed_intersect=TE_bed.intersect(frags_bed, wb=True, sorted=True)


def make_cell_x_element_matrix(bed_interesect, TE_fams, cell_barcodes):
    """
    Takes in a BedTools object after being intersected and the lists of TEs and cell barcodes.
    It iterates through the BedTools object and creates a matrix based on the TE families and the cell barcodes.
    
    Args:
        bed_interesect: pyBedTools object, output from intersecting two BED files
        TE_fams: list, list of TE families
        cell_barcodes: list, list of cell barcodes
        
    Returns:
        sparse_matrix: scipy.sparse matrix, matrix with rows corresponding to cell barcodes and columns corresponding to TE families
    """
    counts = {}
    # To iterate each row in the intsersected file and add 1 count if there is one overlap 
    for interval in bed_interesect:
        fam_name = interval.name
        barcode = interval.fields[8]
        try:
            counts[(barcode, fam_name)] += 1
        except KeyError as e:
            counts[(barcode, fam_name)] = 1

    # use barcode as row index and TE family names as column names
    row_inds = []
    col_inds = []
    data = []
    for (barcode, fam_name), count in counts.items():
        row_inds.append(np.where(cell_barcodes == barcode)[0][0])
        col_inds.append(np.where(TE_fams == fam_name)[0][0])
        data.append(count)
    # Convert to a sparse matrix since the count matrix is full of 0s.
    sparse_matrix = csr_matrix((data, (row_inds, col_inds)), shape=(len(cell_barcodes), len(TE_fams)))
    return sparse_matrix
