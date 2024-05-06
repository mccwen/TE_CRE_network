setwd("/Users/ming-chingwen/Documents/bioinformatics/umich/courses/Winter2024/bioinf576/ciceroEnhancer/")

library(monocle3)
library(cicero)
library(Matrix)
library(tibble)

indata <- readMM("filtered_feature_bc_matrix/matrix.mtx.gz")
# binarize the matrix
indata@x[indata@x > 0] <- 1

# format cell info
cellinfo <- read.table("filtered_feature_bc_matrix/barcodes.tsv.gz")
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"
# get features
features=read.table("filtered_feature_bc_matrix/features.tsv.gz", header=F,
                    row.names=1, sep="\t")
peaks_only <- rownames(features)[features$V3 =="Peaks"]

#Format peak info
peakinfo <- read.table("filtered_feature_bc_matrix/peaks.bed")
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
peakinfo$site_name <- sub("_", ":", peakinfo$site_name, fixed=T)
peakinfo$site_name <- sub("_", "-", peakinfo$site_name, fixed=T)
row.names(peakinfo) <- peakinfo$site_name

rownames(indata) <- row.names(features)
colnames(indata) <- row.names(cellinfo)

# create a function to get highly correlated enhancers of target gene
#gene_enhancer_corr <- function(features){
#  g_chr2=subset(features, (V3=="Gene Expression") &(V4=="chr2") & (V2=="Macrod2"))
#  p_chr2=subset(features, (V3=="Peaks") & (V4=="chr2"))
#  GEX=indata[rownames(g_chr2), ]
#  PEAKs=indata[rownames(p_chr2), ]
#  g_db=as.data.frame(as.matrix(GEX))
#  p_db=as.data.frame(as.matrix(PEAKs))
#  chr_peaks_T=t(p_db)
#  peaks_of_Macrod2=cbind(g_db, chr_peaks_T)
#  coefficients=cor(peaks_of_Macrod2[-1], peaks_of_Macrod2$V1)
#  coe_db=as.data.frame(coefficients)
#  colnames(coe_db) <- "coeff"
#  coe_db <- tibble::rownames_to_column(coe_db, "peak")
#  sig_coe=coe_db[coe_db$coeff >=0.15, ]
#  # drop NA rows
#  complete_dat<- sig_coe[complete.cases(sig_coe), ]
#  return(complete_dat)
#  #write.csv(complete_dat, "results/enhancers_corr_gt0.15_with_Macrod2.csv")
#}

res=gene_enhancer_corr(features, "chr2", "Macrod2")
write.csv(res, "results/enhancers_corr_gt0.15_with_Macrod2.csv", row.names = FALSE)


# only use rows that are peaks
new_indata <- indata[peaks_only, ]
input_cds = suppressWarnings(new_cell_data_set(new_indata,
                                               cell_metadata = cellinfo,
                                               gene_metadata = peakinfo))
# ensure no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) !=0,]
#create a cicero CDS
set.seed(2017)
input_cds <- monocle3::detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method="LSI")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP',
                 preprocess_method = "LSI")

plot_cells(input_cds)
umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# Here we use a mouse ATAC-seq data from 10X multiome as input, so we go with the mouse genom.
# Change to the respective genome if data is of a different species
data("mouse.mm9.genome")

# use only a small part of the genome for testing purposes
# Also, specify the chromosome where the target gene (e.g., chromosome 3 for Macrod) locates
sample_genome <- subset(mouse.mm9.genome, V1 == "chr2")

# run on 2 samples
conns <- run_cicero(cicero_cds, sample_genome, sample_num = 2)
conns_sub <- conns[conns$coaccess >=0.2, ]
# Select peaks with co-accessibility score above 0.2 and save the data
write.csv(conns_sub, "cexample_data/hr2_coaccess_gt0.2.csv", row.names=FALSE)

# Use all data of 100 samples for real
# conns <- run_cicero(cicero_cds, mouse.mm9.genome, sample_num = 100)

