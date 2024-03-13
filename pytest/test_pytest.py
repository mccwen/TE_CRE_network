import pytest
import sys
 
# setting path
sys.path.append('../')
 
# importing
from .. import steamer
#import steamer as st
import pandas as pd
#from pybedtools import BedTool
from fuc import pybed


# Test the first function in steamer
def test_create_bed_for_TEs():
	test_input="test_data/test.tsv"
	result=st.create_bed_for_TEs(test_input)
	expected_result=pd.read_csv("test_data/test.bed", sep="\t")
	expected_bf = pybed.BedFrame.from_frame(meta=[], data= expected_result)
	#expected_bf = pybed.BedFrame.from_frame([],expected_result)
	assert result.to_string() == expected_bf.to_string(), "The tsv function result does not match the bed format."

# Test the second function in steamer
def test_create_bed_for_fragments():
	test_sample="test_data/test_sample.tsv"
	result=st.create_bed_for_fragments(test_sample)
	expected_result=pd.read_csv("test_data/test_sample.bed", sep="\t")
	expected_bf=pybed.BedFrame.from_frame(meta=[], data=expected_result)
	#expected_bf = pybed.BedFrame.from_frame([],expected_result)
	assert result.to_string()==expected_bf.to_string()


# Test the final function in steamer to ensure the sizes of TE names and of cell barcode match those in Bedtool object
# First, make up TE family names and cell barcodes to compare with:
TE_fam=pd.read_csv("test_data/test_TE_barcode.csv").loc[:, "TE"]
barcodes=pd.read_csv("test_data/test_TE_barcode.csv").loc[:, "barcode"]

def test_make_cell_x_element_matrix():
	intersect_db=pd.read_csv("test_data/test_BEDintersection_data.bed", sep="\t", header=None)
	intersect_db=intersect_db.rename(columns={0: "Chromosome", 1: "Start", 2: "End", 3: "TE_family", 4: "Strand", 5: "Chromosome_r", 
                         6: "Start_r", 7: "End_r", 8: "barcode"})

	TEs_from_intersetc=intersect_db.loc[:, "TE_family"]
	barcode_from_intersect=intersect_db.loc[:, "barcode"]
	assert TEs_from_intersetc.equals(other=TE_fam) & barcode_from_intersect.equals(other=barcodes)
	

# The tests below serves as egde cases to check if chromosome names in the input TE file are in the right range.
# First, we define correct chromosome names
mouse_chroms=["chr"+str(i) for i in range(1, 20)]
#chroms=map(range(1:19): lambda x: "chr"+str(x))
mouse_chroms.extend(["chrX", "chrY"])

@pytest.fixture
def define_TE_chromosome_names():
	chrom_names=pd.read_csv("test_data/test.tsv", sep="\t").iloc[:, 0]
	return chrom_names
	

@pytest.fixture
def define_sample_chromosome_names():
	#ch_list=[]
	sample_chrom_names=	pd.read_csv("test_data/test_sample.tsv", sep="\t").iloc[:, 0]
	return sample_chrom_names
	
	

# Next, we check if the chromosomes in the input files all match the correct chromosome names
def test_bed_for_TEs(define_TE_chromosome_names):
    #only taking the columns with the chromsome col names 	
    assert all(x in mouse_chroms for x in define_TE_chromosome_names)
	
# The same check for the input sample file
def test_bed_for_fragments(define_sample_chromosome_names):
	assert all(x in mouse_chroms for x in define_sample_chromosome_names), "Some chromosome names are wrong."
	


