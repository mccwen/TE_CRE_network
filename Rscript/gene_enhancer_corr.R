# create a function to get correlation between target gene and enhancers
gene_enhancer_corr <- function(features, chrom, gene_name){
  g_chr2=subset(features, (V3=="Gene Expression") &(V4==chrom) & (V2==gene_name))
  p_chr2=subset(features, (V3=="Peaks") & (V4==chrom))
  GEX=indata[rownames(g_chr2), ]
  PEAKs=indata[rownames(p_chr2), ]
  g_db=as.data.frame(as.matrix(GEX))
  p_db=as.data.frame(as.matrix(PEAKs))
  chr_peaks_T=t(p_db)
  peaks_of_Macrod2=cbind(g_db, chr_peaks_T)
  coefficients=cor(peaks_of_Macrod2[-1], peaks_of_Macrod2$V1)
  coe_db=as.data.frame(coefficients)
  colnames(coe_db) <- "coeff"
  coe_db <- tibble::rownames_to_column(coe_db, "peak")
  sig_coe=coe_db[coe_db$coeff >=0.15, ]
  # drop NA rows
  complete_dat<- sig_coe[complete.cases(sig_coe), ]
  return(complete_dat)
  #write.csv(complete_dat, "results/enhancers_corr_gt0.15_with_Macrod2.csv")
}
