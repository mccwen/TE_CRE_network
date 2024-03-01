import pytest
import steamer as st
import pandas as pd
#from pybedtools import BedTool
from fuc import pybed


# test the first function in steamer
def test_create_bed_for_TEs():
	test_input="test_data/test.tsv"
	#test_input="test_simple.tsv"
	result=st.create_bed_for_TEs(test_input)
	#expected_result=pd.read_csv("test_simple.bed", sep="\t")
	expected_result=pd.read_csv("test_data/test.bed", sep="\t")
	expected_bf = pybed.BedFrame.from_frame(meta=[], data= expected_result)
	#expected_bf = pybed.BedFrame.from_frame([],expected_result)
	assert result.to_string() == expected_bf.to_string(), "The tsv function result does not match the bed format."

# test the second function in steamer
def test_create_bed_for_fragments():
	test_sample="test_data/test_sample.tsv"
	result=st.create_bed_for_fragments(test_sample)
	expected_result=pd.read_csv("test_data/test_sample.bed", sep="\t")
	expected_bf=pybed.BedFrame.from_frame(meta=[], data=expected_result)
	#expected_bf = pybed.BedFrame.from_frame([],expected_result)
	assert result.to_string()==expected_bf.to_string()

# The test below serves as an egde case to check if chromosome names in the input TE file are in the right range.
# For example, below we provide the correct mouse chromosome names.

# create a test dataset
#TEs=pd.read_csv("mm10.nrph.hits.gz", sep="\t")
test_TEs= pd.DataFrame({"seq": ("chr1", "chr2", "chr6"), "start": (1390, 1400, 1560), "end": (1420, 1450, 1679)})
test_samples= pd.DataFrame({"Chromosome":('chr1', 'chrX', 'chrY100'), "start": (100, 200, 300), 
	"end":(120, 250, 310)})

chroms=["chr"+str(i) for i in range(1, 20)]
#chroms=map(range(1:19): lambda x: "chr"+str(x))
chroms.extend(["chrX", "chrY"])

@pytest.fixture
def define_TE_chromosome_names():
	chrom_names=test_TEs.iloc[:, 0]
	for name in chrom_names:
		if name in chroms:
			return name
	

@pytest.fixture
def define_sample_chromosome_names():
	sample_chrom_names=	test_samples.iloc[:, 0]
	for s_name in sample_chrom_names:
		if s_name in chroms:
			return s_name

# Next, we check if the chromosomes in the input files all match the correct chromosome names
def test_bed_for_TEs(define_TE_chromosome_names):
    #only taking the columns with the chromsome col names 	
	assert define_TE_chromosome_names in chroms



# The same check for the input sample file
def test_bed_for_fragments(define_sample_chromosome_names):
	assert define_sample_chromosome_names in chroms





