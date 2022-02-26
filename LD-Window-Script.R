######################################################################
########################## LD Window Script #########################
######################################################################

#libraries
library(readr)
library(data.table)
library(dplyr)

#What does this script need to do?
#If snp is on different chromosome
#--False Positive 1
#If snp is on same chromosome and has r2 >= 0.1
#--True Positive
##How to calculate r^2
###Corr function in R
#-----But need to square it

#LD_rate_calc <- function(geno.numeric = NULL, qtn.info.dir.root = NULL, qtn.folder.name = NULL,  results.folder.root = NULL, results.folder.name = NULL) {
#  
#}
LD_rate_calc <- function(qtn.info.dir.root = NULL, qtn.folder.name = NULL, results.folder.root = NULL, results.folder.name = NULL, results = NULL, test = NULL, trait = NULL, species = NULL, size = NULL) {
  if(results == "NULL") {
    vGWAS_rate_results <- matrix(0, nrow = 101, ncol = 5)
    vGWAS_rate_results[101, 1] <- "Rates"
    colnames(vGWAS_rate_results) <- c("FILE", "FP_FDR5", "TP_FDR5", "FP_FDR10", "TP_FDR10")
    results.folder.output <- paste0(results.folder.root, results.folder.name, "/")
    print(results.folder.output)
    results.vec <- list.files(results.folder.output)
    for(file in 1:length(results.vec)) {
      print(file)
      result.file <- fread(paste0(results.folder.output, results.vec[file]))
      fdr.five.snps <- result.file %>%
        select(snp, chr, pos, p.value, fdr) %>%
        filter(fdr < 0.05 & p.value < fdr)
      
      fdr.ten.snps <- result.file %>%
        select(snp, chr, pos, p.value, fdr) %>%
        filter(fdr < 0.10 & p.value < fdr)
      vGWAS_rate_results[file, 1] <- results.vec[file]
      if(dim(fdr.five.snps)[1] > 0) {
        vGWAS_rate_results[file, 2] <- 1
      }
      if(dim(fdr.ten.snps)[1] > 0) {
        vGWAS_rate_results[file, 4] <- 1
      }
    }
    vGWAS_rate_results[101, 2] <- sum(as.numeric(vGWAS_rate_results[1:100, 2]))/100
    vGWAS_rate_results[101, 3] <- sum(as.numeric(vGWAS_rate_results[1:100, 3]))/100
    vGWAS_rate_results[101, 4] <- sum(as.numeric(vGWAS_rate_results[1:100, 4]))/100
    vGWAS_rate_results[101, 5] <- sum(as.numeric(vGWAS_rate_results[1:100, 5]))/100
    write.csv(vGWAS_rate_results, file = paste0(test, "_", species, "_", trait, "_", "detection.rates", ".csv"), row.names = F)
  }
  if(results == "vQTL") {
    vGWAS_rate_results <- matrix(NA, nrow = 101, ncol = 5)
    colnames(vGWAS_rate_results) <- c("FILE", "FP_FDR5", "TP_FDR5", "FP_FDR10", "TP_FDR10")
    simulated.qtn.folder <- paste0(qtn.info.dir.root, qtn.folder.name, "/")
    vQTN.info <- fread(file = paste0(simulated.qtn.folder, "var.QTN.genotypic.information.csv"))
    vQTN.info <- vQTN.info[, c(2, 4, 5)]
    simchr <- vQTN.info$chr
    simSNP <- vQTN.info$snp
    results.folder.output <- paste0(results.folder.root, results.folder.name, "/")
    results.vec <- list.files(results.folder.output)
    for(file in 1:length(results.vec)) {
      result.file <- fread(paste0(results.folder.output, results.vec[file]))
      fdr.five.snps <- result.file %>%
        select(snp, chr, pos, p.value, fdr) %>%
        filter(fdr < 0.05 & p.value < fdr)
      
      fdr.ten.snps <- result.file %>%
        select(snp, chr, pos, p.value, fdr) %>%
        filter(fdr < 0.10 & p.value < fdr)
      if(dim(fdr.five.snps[fdr.five.snps$chr != simchr, ])[1] > 0) {
        vGWAS_rate_results[file, 2] <- 1
      }
      if(dim(fdr.ten.snps[fdr.ten.snps$chr != simchr, ])[1] > 0) {
        vGWAS_rate_results[file, 4] <- 1
      }
      fdr.five.snps.samechr <- fdr.five.snps[fdr.five.snps$chr == simchr, ]
      fdr.ten.snps.samechr <- fdr.ten.snps[fdr.ten.snps$chr == simchr, ]
      snps.4.ld.fdr5 <- grep(TRUE, colnames(geno) %in% fdr.five.snps.samechr$snp)
      where.simulatedQTN.is.at.fdr.5 <- grep(TRUE, colnames(geno) %in% simSNP)
      geno.4.ld.fdr5 <- geno[, c(where.simulatedQTN.is.at.fdr5, snps.4.ld.fdr5)]
      genosnpnames.fdr5 <- colnames(geno.4.ld.fdr5)
      ldobject.fdr5 <- cor(geno.4.ld.fdr5)
      ldobject.fdr5 <- ldobject^2
      rownames(ldobject.fdr5) <- genosnpnames.fdr5
      colnames(ldobject.fdr5) <- genosnpnames.fdr5
      if(max(ldobject.fdr5[-c(1), 1]) >= 0.1) {
        vGWAS_rate_results[file, 4] <- 1
      }
    }
  }
#  simulated.qtn.folder <- paste0(qtn.info.dir.root, qtn.folder.name, "/")
#  vQTN.info <- fread(file = paste0(simulated.qtn.folder, "var.QTN.genotypic.information.csv"))
#  vQTN.info <- vQTN.info[, c(2, 4, 5)]
#  vGWAS_rate_results <- matrix(NA, nrow = 100, ncol = 5)
#  colnames(vGWAS_rate_results) <- c("file", "FP_FDR5", "TP_FDR5", "FP_FDR10", "TP_FDR10")
} 

BFT.results.full.root <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/vGWAS_results/BFT/FULL/"

LD_rate_calc(qtn.info.dir.root = NULL, qtn.folder.name = NULL, results.folder.root = BFT.results.full.root, results.folder.name = "BFT_sp.zm.null.full", results = "NULL", test = "BFT", species = "zm", trait = "NULL", size = "2815")





simulated.vQTNS.folder <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/zm_vQTL_simulations/"

results.files <- list.files(results.folder)

LD_rate_calc(qtn.info.dir.root = simulated.vQTNS.folder, qtn.folder.name = "zm_vmolike_MAF40_h263_vQTN50_FULL")

#test run with GxE
results.folder <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/vGWAS_results/BFT/FULL/"
results.folder.name <- paste0(results.folder, "BFT_sp.zm.GxE.FULL.mod")
results.files <- list.files(results.folder)
vGWAS_rate_results <- matrix(NA, nrow = 100, ncol = 5)
colnames(vGWAS_rate_results) <- c("file", "FP_FDR5", "TP_FDR5", "FP_FDR10", "TP_FDR10")
simulated.QTNs <- read.delim("~/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/zm_sp_results/sp_zm_GxE_FULL/Additive_Selected_QTNs.txt")
for(i in 1:nrow(simulated.QTNs)) {
  simQTN <- simulated.QTNs[i, ]
  simchr <- simQTN$chr
  simSNP <- simQTN$snp
  print(simchr)
  for(file in 1:length(results.files)) {
    result.file <- fread(paste0(results.folder, results.files[file]))
    vGWAS_rate_results[file, 1] <- results.files[file]
    fdr.five.snps <- result.file %>%
      select(snp, chr, pos, p.value, fdr) %>%
      filter(fdr < 0.05 & p.value < fdr)
  
    fdr.ten.snps <- result.file %>%
      select(snp, chr, pos, p.value, fdr) %>%
      filter(fdr < 0.10 & p.value < fdr)
#  simulated.QTNs <- read.delim("~/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/zm_sp_results/sp_zm_GxE_FULL/Additive_Selected_QTNs.txt")
#  for(i in 1:nrow(simulated.QTNs)) {
#    simQTN <- simulated.QTNs[i, ]
#    simchr <- simQTN$chr
#    simSNP <- simQTN$snp
#    print(simchr)
    if(dim(fdr.five.snps[fdr.five.snps$chr != simchr, ])[1] > 0) {
      vGWAS_rate_results[file, 2] <- 1
    }
    if(dim(fdr.ten.snps[fdr.ten.snps$chr != simchr, ])[1] > 0) {
      vGWAS_rate_results[file, 4] <- 1
    }
    write.csv(vGWAS_rate_results, file = paste(simSNP, ".csv", sep=""), sep = ",", row.names = FALSE)
  #  snpsonchr <- fdr.five.snps[fdr.five.snps$chr == simchr, ]
#  sigsnps <- grep(TRUE, colnames(zm.just.geno.full) %in% snpsonchr$snp)
#  geno.4.ld <- zm.just.geno.full[, c(sigsnps)]
#  jc <- cor(geno.4.ld)
#  creepers <- jc^2
  }
}


result.file <- fread(paste0(results.folder, results.files[1])) #iterate through loop
#fdr 
fdr.five.snps <- result.file %>%
  select(snp, chr, pos, p.value, fdr) %>%
  filter(fdr < 0.05 & p.value < fdr)

fdr.ten.snps <- vGWAS_results_raw %>%
  select(snp, chr, pos, p.value, fdr) %>%
  filter(fdr < 0.10 & p.value < fdr)

#Reading simulated QTNs
simulated.QTNs <- read.delim("~/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/zm_sp_results/sp_zm_GxE_FULL/Additive_Selected_QTNs.txt")
for(i in 1:nrow(simulated.QTNs)) {
  simQTN <- simulated.QTNs[i, ]
  simchr <- simQTN$chr
  print(simchr)
}

dim(fdr.five.snps[fdr.five.snps$chr != 6, ])[1] > 0

snpsonchr <- fdr.five.snps[fdr.five.snps$chr == 6, ]
sigsnps <- grep(TRUE, colnames(zm.just.geno.full) %in% snpsonchr$snp)
geno.4.ld <- zm.just.geno.full[, c(sigsnps)]
jc <- cor(geno.4.ld)
"S6_162917265" %in% colnames(jc)

################################


Rate.calc <- function(Result_dir = NULL, rf = NULL, Null_results = NULL, simvQTN = NULL, mQTN.wd = NULL, DGLM = NULL, Additive = NULL, trait.name = NULL, test = NULL) {
  if(Null_results == TRUE) {
    BP_five_rate_matrix <- matrix(NA, nrow = 100, ncol = 2)
    colnames(BP_five_rate_matrix) <- c("BP_5_FP", "BP_5_TP")
    BP_ten_rate_matrix <- matrix(NA, nrow = 100, ncol = 2)
    colnames(BP_ten_rate_matrix) <- c("BP_10_FP", "BP_10_TP")
    FDR_five_rate_matrix <- matrix(NA, nrow = 100, ncol = 2)
    colnames(FDR_five_rate_matrix) <- c("FDR_5_FP", "FDR_5_TP")
    FDR_ten_rate_matrix <- matrix(NA, nrow = 100, ncol = 2)
    colnames(FDR_ten_rate_matrix) <- c("FDR_10_FP", "FDR_10_TP")
    vGWAS_rate_results <- matrix(NA, nrow = 100, ncol = 9)
    dirresults <- paste0(Result_dir[rf], "/")
    Result_file_vec <- list.files(dirresults)
    for(i in 1:length(Result_file_vec)) {
      print(i)
      if(DGLM == TRUE) {
        vGWAS_results_raw <- read.csv(paste0(dirresults, Result_file_vec[i]))
        vGWAS_results_raw <- vGWAS_results_raw[, c(1, 2, 3, 7, 9)]
        colnames(vGWAS_results_raw) <- c("snp", "chr", "pos", "p.value", "fdr")
      } else {
        vGWAS_results_raw <- read.csv(paste0(dirresults, Result_file_vec[i]))
      }
      N <- nrow(vGWAS_results_raw)
      alpha_ten <- 0.10
      alpha_five <- 0.05
      Bon_five <- 0.05 / N
      Bon_ten <- 0.10 / N
      print("Finding number of significant SNPs")
      numb_sig_BPfive <- nrow(vGWAS_results_raw[vGWAS_results_raw$p.value < Bon_five, ])
      numb_sig_BPten <- nrow(vGWAS_results_raw[vGWAS_results_raw$p.value < Bon_ten, ])
      numb_fdrfive <- nrow(vGWAS_results_raw[vGWAS_results_raw$fdr < 0.05, ])
      numb_fdrten <- nrow(vGWAS_results_raw[vGWAS_results_raw$fdr < 0.10, ])
      if(numb_sig_BPfive >= 1) {
        BP_five_rate_matrix[i, 1] <- 1
        BP_five_rate_matrix[i, 2] <- 0
      } else {
        if(numb_sig_BPfive == 0)
          BP_five_rate_matrix[i, 1] <- 0
        BP_five_rate_matrix[i, 2] <- 0
      }
      if(numb_sig_BPten >= 1) {
        BP_ten_rate_matrix[i, 1] <- 1
        BP_ten_rate_matrix[i, 2] <- 0
      } else {
        if(numb_sig_BPten == 0)
          BP_ten_rate_matrix[i, 1] <- 0
        BP_ten_rate_matrix[i, 2] <- 0
      }
      if(numb_fdrfive >= 1) {
        FDR_five_rate_matrix[i, 1] <- 1
        FDR_five_rate_matrix[i, 2] <- 0
      } else {
        if(numb_fdrfive == 0)
          FDR_five_rate_matrix[i, 1] <- 0
        FDR_five_rate_matrix[i, 2] <- 0
      }
      if(numb_fdrten >= 1) {
        FDR_ten_rate_matrix[i, 1] <- 1
        FDR_ten_rate_matrix[i, 2] <- 0
      } else {
        if(numb_fdrten == 0)
          FDR_ten_rate_matrix[i, 1] <- 0
        FDR_ten_rate_matrix[i, 2] <- 0
      }
      vGWAS_rate_results[i, 1] <- i
      vGWAS_rate_results[i, 2:3] <- BP_five_rate_matrix[i, 1:2]
      vGWAS_rate_results[i, 4:5] <- BP_ten_rate_matrix[i, 1:2]
      vGWAS_rate_results[i, 6:7] <- FDR_five_rate_matrix[i, 1:2]
      vGWAS_rate_results[i, 8:9] <- FDR_ten_rate_matrix[i, 1:2]
      colnames(vGWAS_rate_results) <- c("Rep", "BP_5_FP", "BP_5_TP", "BP_10_FP", "BP_10_TP", "FDR_5_FP", "FDR_5_TP", "FDR_10_FP", "FDR_10_TP")
    }
    write.csv(vGWAS_rate_results, file = paste(test, trait.name, "_", ".csv", sep = ""), sep = ",", row.names = FALSE)
  } 
  print("Reading in Simulated Additive QTNs")
  if(Additive == TRUE) {
    simmQTN <- read.table(paste0(mQTN.wd, list.files(mQTN.wd)[1]), header = T)
  }
  dirresults <- paste0(Result_dir[rf], "/")
  Result_file_vec <- list.files(dirresults)
  for(h in 1:dim(simvQTN)) { 
    print(h)
    BP_five_rate_matrix <- matrix(NA, nrow = 100, ncol = 2)
    colnames(BP_five_rate_matrix) <- c("BP_5_FP", "BP_5_TP")
    BP_ten_rate_matrix <- matrix(NA, nrow = 100, ncol = 2)
    colnames(BP_ten_rate_matrix) <- c("BP_10_FP", "BP_10_TP")
    FDR_five_rate_matrix <- matrix(NA, nrow = 100, ncol = 2)
    colnames(FDR_five_rate_matrix) <- c("FDR_5_FP", "FDR_5_TP")
    FDR_ten_rate_matrix <- matrix(NA, nrow = 100, ncol = 2)
    colnames(FDR_ten_rate_matrix) <- c("FDR_10_FP", "FDR_10_TP")
    vGWAS_rate_results <- matrix(NA, nrow = 100, ncol = 9)
    for(i in 1:length(Result_file_vec)) {
      print(i)
      if(DGLM==TRUE) {
        vGWAS_results_raw <- read.csv(paste0(dirresults, Result_file_vec[i]))
        vGWAS_results_raw <- vGWAS_results_raw[, c(1, 2, 3, 7, 9)]
        colnames(vGWAS_results_raw) <- c("snp", "chr", "pos", "p.value", "fdr")
      } else {
        vGWAS_results_raw <- read.csv(paste0(dirresults, Result_file_vec[i]))
      }
      svQTNs <- grep(TRUE, vGWAS_results_raw$snp %in% simvQTN$snp)
      if(Additive == TRUE) {
        saQTNs <- grep(TRUE, vGWAS_results_raw$snp %in% simmQTN$snp)
        vGWAS_results_raw <- vGWAS_results_raw[-c(saQTNs, svQTNs), ]
      }
      vGWAS_results_raw <- vGWAS_results_raw[-c(svQTNs), ]
      single_simvQTN <- simvQTN[h, ]
      N <- nrow(vGWAS_results_raw)
      alpha_ten <- 0.10
      alpha_five <- 0.05
      Bon_five <- 0.05 / N
      Bon_ten <- 0.10 / N
      print("Finding number of significant SNPs")
      numb_sig_BPfive <- nrow(vGWAS_results_raw[vGWAS_results_raw$p.value < Bon_five, ])
      numb_sig_BPten <- nrow(vGWAS_results_raw[vGWAS_results_raw$p.value < Bon_ten, ])
      numb_fdrfive <- nrow(vGWAS_results_raw[vGWAS_results_raw$fdr < 0.05, ])
      numb_fdrten <- nrow(vGWAS_results_raw[vGWAS_results_raw$fdr < 0.10, ])
      
      #BF_significant snps
      bon_five_snps <-vGWAS_results_raw[vGWAS_results_raw$p.value < Bon_five, ]
      bon_ten_snps <- vGWAS_results_raw[vGWAS_results_raw$p.value < Bon_ten, ]
      
      #Snps not significant by BF
      bon_five_failed <-vGWAS_results_raw[vGWAS_results_raw$p.value > Bon_five, ]
      bon_ten_failed <- vGWAS_results_raw[vGWAS_results_raw$p.value > Bon_ten, ]
      
      print("Finding signficant snps of FDR")
      #fdr 
      fdr_five_snps <- vGWAS_results_raw %>%
        select(snp, chr, pos, p.value, fdr) %>%
        filter(fdr < 0.05 & p.value < fdr)
      
      fdr_ten_snps <- vGWAS_results_raw %>%
        select(snp, chr, pos, p.value, fdr) %>%
        filter(fdr < 0.10 & p.value < fdr)
      
      ###Snps not signficant by fdr
      fdr_five_failed <- vGWAS_results_raw %>%
        select(snp, chr, pos, p.value, fdr) %>%
        filter(fdr > 0.05)
      
      fdr_ten_failed <- vGWAS_results_raw %>%
        select(snp, chr, pos, p.value, fdr) %>%
        filter(fdr > 0.10)
      
      trueQTN <- vGWAS_results_raw[vGWAS_results_raw$chr == single_simvQTN$chr, ]
      trueQTN <- trueQTN[trueQTN$pos >= single_simvQTN$QTNwindow_lowerbound, ]
      trueQTN <- trueQTN[trueQTN$pos <= single_simvQTN$QTNwindow_upperbound, ]
      #This section is to see how many of my simulated QTNs are in the significant snps
      trueQTN_insig_fdr5 <- length(grep(TRUE, trueQTN$snp %in% fdr_five_snps$snp))
      trueQTN_insig_fdr10 <- length(grep(TRUE, trueQTN$snp %in% fdr_ten_snps$snp))
      trueQTN_insig_BP5 <- length(grep(TRUE, trueQTN$snp %in% bon_five_snps$snp))
      trueQTN_insig_BP10 <- length(grep(TRUE, trueQTN$snp %in% bon_ten_snps$snp))
      
      #This section is for true QTNs labelled as negative, even though they should be true positives
      trueQTN_notsig_fdr5 <- length(grep(TRUE, trueQTN$snp %in% fdr_five_failed$snp))
      trueQTN_notsig_fdr10 <- length(grep(TRUE, trueQTN$snp %in% fdr_ten_failed$snp))
      trueQTN_notsig_BP5 <- length(grep(TRUE, trueQTN$snp %in% bon_five_failed$snp))
      trueQTN_notsig_BP10 <- length(grep(TRUE, trueQTN$snp %in% bon_ten_failed$snp))
      if(length(grep(FALSE, bon_five_snps[, 1] %in% trueQTN$snp)) >= 1) {
        BP_five_rate_matrix[i, 1] <- 1
      } else {
        if(length(grep(FALSE, bon_five_snps[, 1] %in% trueQTN$snp)) == 0)
          BP_five_rate_matrix[i, 1] <- 0
      }
      if(length(grep(FALSE, bon_ten_snps[,1] %in% trueQTN$snp)) >=1 ) {
        BP_ten_rate_matrix[i, 1] <- 1
      } else {
        if(length(grep(FALSE, bon_ten_snps[, 1] %in% trueQTN$snp)) == 0)
          BP_ten_rate_matrix[i, 1] <- 0
      }
      if(length(grep(FALSE, fdr_five_snps[, 1] %in% trueQTN$snp)) >= 1) {
        FDR_five_rate_matrix[i, 1] <- 1
      } else {
        if(length(grep(FALSE, fdr_five_snps[, 1] %in% trueQTN$snp)) == 0)
          FDR_five_rate_matrix[i, 1] <- 0
      }
      if(length(grep(FALSE, fdr_ten_snps[, 1] %in% trueQTN$snp)) >= 1) {
        FDR_ten_rate_matrix[i, 1] <- 1
      } else {
        if(length(grep(FALSE, fdr_ten_snps[, 1] %in% trueQTN$snp)) ==0 )
          FDR_ten_rate_matrix[i, 1] <- 0
      }
      if(length(grep(TRUE, fdr_five_snps[, 1] %in% trueQTN$snp)) >=1 ) {
        FDR_five_rate_matrix[i, 2] <- 1
      } else {
        if(length(grep(TRUE, fdr_five_snps[, 1] %in% trueQTN$snp)) == 0)
          FDR_five_rate_matrix[i, 2] <- 0
      }
      if(length(grep(TRUE, fdr_ten_snps[, 1] %in% trueQTN$snp)) >= 1) {
        FDR_ten_rate_matrix[i, 2] <- 1
      } else {
        if(length(grep(TRUE, fdr_ten_snps[, 1] %in% trueQTN$snp)) == 0)
          FDR_ten_rate_matrix[i, 2] <- 0
      }
      if(length(grep(TRUE, bon_five_snps[, 1] %in% trueQTN$snp)) >= 1) {
        BP_five_rate_matrix[i, 2] <- 1
      } else {
        if(length(grep(TRUE, bon_five_snps[, 1] %in% trueQTN$snp)) == 0)
          BP_five_rate_matrix[i, 2] <- 0
      }
      if(length(grep(TRUE, bon_ten_snps[, 1] %in% trueQTN$snp)) >= 1) {
        BP_ten_rate_matrix[i, 2] <- 1
      } else {
        if(length(grep(TRUE, bon_ten_snps[, 1] %in% trueQTN$snp)) == 0)
          BP_ten_rate_matrix[i, 2] <- 0
      }
      vGWAS_rate_results[i, 1] <- Result_file_vec[i]
      vGWAS_rate_results[i, 2:3] <- BP_five_rate_matrix[i, 1:2]
      vGWAS_rate_results[i, 4:5] <- BP_ten_rate_matrix[i, 1:2]
      vGWAS_rate_results[i, 6:7] <- FDR_five_rate_matrix[i, 1:2]
      vGWAS_rate_results[i, 8:9] <- FDR_ten_rate_matrix[i, 1:2]
      colnames(vGWAS_rate_results) <- c("Rep", "BP_5_FP", "BP_5_TP", "BP_10_FP", "BP_10_TP", "FDR_5_FP", "FDR_5_TP", "FDR_10_FP", "FDR_10_TP")
    }
    write.csv(vGWAS_rate_results, file = paste(test, trait.name, "_", simvQTN[h, 2], ".csv", sep=""), sep = ",", row.names = FALSE)
    
  }
  return(vGWAS_rate_results)
}
