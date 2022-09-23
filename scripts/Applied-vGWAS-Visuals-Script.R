###############################################################
############# Applied vGWAS Visualization Script #############
###############################################################


install.packages("CMplot")
library(CMplot)

##Read in results files
BFT_BG282_PH_results <- read.csv("~/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/Applied_vGWAS/BF_BG282_PH.csv")
DGLM_BG282_PH_results <- read.csv("~/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/Applied_vGWAS/Zm_PH_DGLM_Results.csv")
MLM_BG282_PH_results <- read.delim("~/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/Applied_vGWAS/zm_ph_GWAS_SNP55K_maize282_AGP3_20190419_Filter_+_Peiffer_hgt_blups_modified_stats.txt")

#We need to remove the first row of MLM results from TASSEL
MLM_BG282_PH_results <- MLM_BG282_PH_results[-c(1), ]

#This section is just a sanity check 
length(grep(TRUE, BFT_BG282_PH_results$snp == DGLM_BG282_PH_results$SNP))
length(grep(TRUE, BFT_BG282_PH_results$snp == MLM_BG282_PH_results$Marker))
length(grep(TRUE, DGLM_BG282_PH_results$SNP == MLM_BG282_PH_results$Marker))

as.factor(BFT_BG282_PH_results$snp)

##SNP and Chromosome needs to be factors
PH_all_results <- data.frame(BFT_BG282_PH_results[, 1:4], DGLM_BG282_PH_results$P.disp, MLM_BG282_PH_results$p) 
colnames(PH_all_results) <- c("SNP", "Chromosome", "Position", "BFT", "DGLM", "MLM")

#then, I need to switch the columns up 
PH_all_results <- PH_all_results[, c(1:3, 6, 4, 5)]

##Here, I will switch the wporking directory to my figures folder
setwd("~/Ph.D._Dissertation/Projects/Dissertation_research/C1/3_output/Figs/")

## I am going to create a Circular manhattan plot for all tests
CMplot(PH_all_results, type = "p", plot.type = "c", chr.labels = paste("Chr", c(1:10), sep = ""), r = 0.4, cir.legend = TRUE,
       outward = T, cir.legend.col = "black", cir.chr.h = 1.3, chr.den.col = "black", file = "jpg",
       memo = "", dpi = 300, file.output = TRUE, verbose = TRUE, width = 10, height = 10)


###Know, I will create QQPLOTS using qqman R package
qq(PH_all_results$DGLM, main = "DGLM")
qq(PH_all_results$BFT, main = "BFT")
qq(PH_all_results$MLM, main = "MLM")

##Here, I will perform fdr using Benjamini Hochberg Procedure 
PH_all_results$MLM.fdr <- p.adjust(PH_all_results$MLM, "fdr")
PH_all_results$BFT.fdr <- p.adjust(PH_all_results$BFT, "fdr")
PH_all_results$DGLM.fdr <- p.adjust(PH_all_results$DGLM, "fdr")

