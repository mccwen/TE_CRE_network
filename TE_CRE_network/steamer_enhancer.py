#This module was primarily developed by and taken from the Welch lab with some modifications
import pandas as pd
import numpy as np
import random
from pybedtools import BedTool
import scipy.sparse
from fuc import pybed
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix
from collections import defaultdict
import os
import glob
import csv
import gzip


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

def get_enhancers(file):    
    # Load the data into a DataFrame
    df = pd.read_csv(file, sep=',')
    # Split each pair of peaks into two sets of separate columns
    peaks_split = df['Peak1'].str.split('-', expand=True)
    peak1_chr = peaks_split[0].str.split(':', expand=True)[0]
    start_position = peaks_split[0].str.split(':', expand=True)[1]
    end_position = peaks_split[1]
    # peaks_split.columns = ['Peak1_start', 'Peak1_end']

    peaks_split2 = df['Peak2'].str.split('-', expand=True)
    peak2_chr = peaks_split2[0].str.split(':', expand=True)[0]
    start_position2 = peaks_split2[0].str.split(':', expand=True)[1]
    end_position2 = peaks_split2[1]

    # Vertically concatenate the two sets of columns
    new_df = pd.DataFrame({
        'chr': pd.concat([peak1_chr, peak2_chr], ignore_index=True),
        'start_position': pd.concat([start_position, start_position2], ignore_index=True),
        'end_position': pd.concat([end_position, end_position2], ignore_index=True)
    })
    new_df['start_position'] = new_df['start_position'].astype(int)
    new_df['end_position'] = new_df['end_position'].astype(int)
    return new_df

def get_nearby_enhancers(df,start_position,end_position):
    filtered_index = (df['start_position'] > start_position - 500000) & (df['end_position'] < end_position + 500000)
    enhancer = df[filtered_index]
    # Gnerate random barcode
    def generate_barcode():
        return ''.join(random.choices('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789', k=10))
    # Create a new DataFrame with the required columns
    bed_df = pd.DataFrame({
        'Chromosome': enhancer['chr'],
        'Start': enhancer['start_position'],
        'End': enhancer['end_position'],
        'Barcode': [generate_barcode() for _ in range(len(enhancer))]
    })
    barcode=bed_df["Barcode"]
    # Reorder columns to match the BED format
    bed_df = bed_df[['Chromosome', 'Start', 'End', 'Barcode']]
    barcode=bed_df["Barcode"]
    bed_df=pybed.BedFrame.from_frame(meta=[], data=bed_df)
    bed_df=bed_df.sort()
    bed_df.to_file('Frag.bed')    
    return bed_df, barcode

def get_nearby_enhancers(df,chrom, start_position,end_position):
    df_sub=df[df["chr"]==chrom]
    filtered_index = (df_sub['start_position'] > start_position - 500000) & (df_sub['end_position'] < end_position + 500000)
    enhancer = df_sub[filtered_index]
    # Gnerate random barcode
    def generate_barcode():
        return ''.join(random.choices('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789', k=10))
    # Create a new DataFrame with the required columns
    bed_df = pd.DataFrame({
        'Chromosome': enhancer['chr'],
        'Start': enhancer['start_position'],
        'End': enhancer['end_position'],
        'Barcode': [generate_barcode() for _ in range(len(enhancer))]
    })
    barcode=bed_df["Barcode"]
    # Reorder columns to match the BED format
    bed_df = bed_df[['Chromosome', 'Start', 'End', 'Barcode']]
    barcode=bed_df["Barcode"]
    bed_df=pybed.BedFrame.from_frame(meta=[], data=bed_df)
    bed_df=bed_df.sort()
    bed_df.to_file('Frag.bed')
    return bed_df, barcode

def intersection(TE_bed, frag_bed):
    """
    This function takes a TE annotaiton file and a fragment file (both in .bed format) and use bedtools 
    to intersect them for finding overlaps. 

    The output file is another bed file.

    """
    bed_intersect=TE_bed.intersect(frag_bed, wb=True, sorted=True)
    return bed_intersect
    #return TE_bed, frag_bed

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

def convert_cell_x_element_matrix_to_file(matrix, col, rows,filename='Cell_x_Element_Matrix.csv'):
    '''
    Converts a scipy sparse matrix into a CSV file.

    Args:
        matrix (csr_matrix): The input matrix to convert.
        col (list): A list of column labels for the output CSV file.
        rows (list): A list of row labels for the output CSV file.

    Returns:
        None
    '''
    df = pd.DataFrame.sparse.from_spmatrix(matrix, columns=col,index = rows)
    df.to_csv(filename)

def save_df_as_gz(df, filename):
    # Save the DataFrame as a gzipped file
    with gzip.open(filename, 'wt', encoding='utf-8') as gz_file:
        df.to_csv(gz_file, sep='\t', index=False, header=False)


def prepare_df(df):
    'Swaps the barcodes and TE column indexes to prepare for scanpy.'
    df[df.columns[0]], df[df.columns[1]] = df[df.columns[1]], df[df.columns[0]]
    return df


def convert_df_to_sparse(df):
    ''
    df = prepare_df(df)
    array = df.values
    rows = [val[0] for val in array]
    columns = [val[1] for val in array]
    values = [val[2] for val in array]
    csr_mat = csr_matrix((values, (rows, columns)))
    return csr_mat

def compress_tsv_file(file_path, output_dir,barcode_dict):
    '''
    Takes in a file_path,out_dir, and the barcodes to make a .gzip .tsv file
    '''
    # Create the output file path
    output_file_path = os.path.join(output_dir, f'{os.path.basename(file_path)}.gz')
    # Write the TSV file
    with open(file_path, 'w', newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        for barcode in barcode_dict.keys():
            writer.writerow([barcode])
    # Compress the TSV file to gzip format
    with open(file_path, 'rb') as f_in, gzip.open(output_file_path, 'wb') as f_out:
        f_out.writelines(f_in)
    # Remove the original TSV file
    os.remove(file_path)

def make_features_df(TE_dict):
    # Convert TEs_fams dictionary to a DataFrame
    fams_df = pd.DataFrame(list(TE_dict.items()), columns=['TE_Name', 'idx'])
    # Drop the 'idx' column
    fams_df.drop('idx', axis=1, inplace=True)
    # Move the 'pseudoID' column to the first position
    fams_df.insert(0, 'pseudoID', fams_df['TE_Name'])
    # Add the 'expression' column
    fams_df['expression'] = ['Gene Expression'] * len(TE_dict)
    return fams_df

def compress_sparse_matrix(matrix, file_path):
    # Save the sparse matrix in Matrix Market format to a temporary uncompressed file
    mmwrite(file_path,matrix)
    # Compress the file using gzip
    with open(file_path, 'rb') as f_in:
        with gzip.open(file_path + '.gz', 'wb') as f_out:
            f_out.writelines(f_in)
    
    # Remove the temporary uncompressed file
    os.remove(file_path)




