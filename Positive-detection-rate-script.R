############################################################################################
############################## Positive-Detection Rate Script ##############################
############################################################################################

#In this script, we will be using this script to process true and false positive detection
#rates

###As usual, load your libraries of interest###
library(dplyr)
library(data.table)

#Here, I will create paths to my working directories for each species where the BF and DGLM
#results live
wd$vGWAS_chapter_one_at <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_at/"
wd$at_vGWAS_results <- paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/")
wd$at_bf_results <- paste0(wd$at_vGWAS_results, "BF/")
wd$at_dglm_results <- paste0(wd$at_vGWAS_results, "DGLM/")

#This is for maize
wd$vGWAS_chapter_one_zm <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_zm/"
wd$zm_vGWAS_results <- paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/")
wd$zm_bf_results <- paste0(wd$zm_vGWAS_results, "BF/")
wd$zm_dglm_results <- paste0(wd$zm_vGWAS_results, "DGLM/")

###Need to remove selected vQTNs in this script
####In this section, I will upload the simulated vQTN information from maize and Arabidopsis
###A thaliana
wd$at_vQTL_simulations <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_at/at_vQTL_simulations/"

wd$at_voillike_womQTL_MAF10_500 <- paste0(wd$at_vQTL_simulations, "at_voillike_womQTL_MAF10_500/")
wd$at_voillike_womQTL_MAF40_500 <- paste0(wd$at_vQTL_simulations, "at_voillike_womQTL_MAF40_500/")
wd$at_voillike_womQTL_MAF10_full <- paste0(wd$at_vQTL_simulations, "at_voillike_womQTL_MAF10_full/")
wd$at_voillike_womQTL_MAF40_full <- paste0(wd$at_vQTL_simulations, "at_voillike_womQTL_MAF40_full/")
wd$at_vmolike_MAF40_500 <- paste0(wd$at_vQTL_simulations, "at_vmolike_MAF40_500/")
wd$at_vmolike_MAF40_FULL <- paste0(wd$at_vQTL_simulations, "at_vmolike_MAF40_full/")
wd$at_vmolike_womQTL_MAF40_500 <- paste0(wd$at_vQTL_simulations, "at_vmolike_womQTL_MAF40_500/")
wd$at_vmolike_womQTL_MAF40_FULL <- paste0(wd$at_vQTL_simulations, "at_vmolike_womQTL_MAF40_full/")

###maize
wd$zm_vQTL_simulations <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_zm/zm_vQTL_simulations/"
wd$zm_oillike_womQTL_MAF10_500 <- paste0(wd$zm_vQTL_simulations, "zm_voillike_womQTL_MAF10_500/")
wd$zm_oillike_womQTL_MAF10_full <- paste0(wd$zm_vQTL_simulations, "zm_voillike_womQTL_MAF10_full/")
wd$zm_oillike_womQTL_MAF40_500 <- paste0(wd$zm_vQTL_simulations, "zm_voillike_womQTL_MAF40_500/")
wd$zm_oillike_womQTL_MAF40_full <- paste0(wd$zm_vQTL_simulations, "zm_voillike_womQTL_MAF40_full/")
wd$zm_molike_MAF40_500 <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF40_500/")
wd$zm_molike_MAF40_full <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF40_full/")
wd$zm_molike_womQTL_MAF40_500 <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_womQTL_MAF40_500/")
wd$zm_molike_womQTL_MAF40_full <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_womQTL_MAF40_full/")


#We also need the selected additive QTN from simplePHENOTYPES
#A thaliana
wd$at_sp_simulations <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_at/at_sp_results/"

wd$sp_at_null <- paste0(wd$at_sp_simulations, "sp_at_null/")
wd$sp_at_null_500 <- paste0(wd$at_sp_simulations, "sp_at_null_500/")
wd$sp_at_null_full <- paste0(wd$at_sp_simulations, "sp_at_null_full/")
wd$sp_at_molike_MAF40_500 <- paste0(wd$at_sp_simulations, "sp_at_molike_wmQTL_500/")
wd$sp_at_molike_MAF40_full <- paste0(wd$at_sp_simulations, "sp_at_molike_wmQTL_full/")

#Maize
wd$zm_sp_simulations <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_zm/zm_sp_results/"

wd$sp_zm_null <- paste0(wd$Zm_SP_simulations, "sp_zm_null/")
wd$sp_zm_null_500 <- paste0(wd$zm_sp_simulations, "sp_zm_null_500/")
wd$sp_zm_molike_MAF40_500 <- paste0(wd$zm_sp_results, "sp_zm_molike_wmQTL_500/")
wd$sp_zm_molike_MAF40_full <- paste0(wd$zm_sp_results, "sp_zm_molike_wmQTL_full/")

###In this section below, I will upload all the simulated QTN information from both SP and
#my vQTLs

#A thaliana
##The oillike settings
at_var.QTN.genotypic.information_oillike_momQTL_MAF10_500 <- read.csv(paste0(wd$at_voillike_womQTL_MAF10_500, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_oillike_momQTL_MAF40_500 <- read.csv(paste0(wd$at_voillike_womQTL_MAF40_500, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_oillike_momQTL_MAF10_full <- read.csv(paste0(wd$at_voillike_womQTL_MAF10_full, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_oillike_momQTL_MAF40_full <- read.csv(paste0(wd$at_voillike_womQTL_MAF40_full, "var.QTN.genotypic.information.csv"))

##The Molike Settings
at_var.QTN.genotypic.information_molike_MAF40_500 <- read.csv(paste0(wd$At_vmolike_MAF40_500, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_molike_MAF40_full <- read.csv(paste0(wd$At_vmolike_MAF40_full, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_molike_womQTL_MAF40_500 <- read.csv(paste0(wd$At_vmolike_womQTL_MAF40_500, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_molike_womQTL_MAF40_full <- read.csv(paste0(wd$At_vmolike_womQTL_MAF40_full, "var.QTN.genotypic.information.csv"))

#Maize
##oillike settings
zm_var.QTN.genotypic.information_oillike_womQTL_MAF10_500 <- read.csv(paste0(wd$Zm_oillike_womQTL_MAF10_500, "var.QTN.genotypic.information.csv"))
zm_var.QTN.genotypic.information_oillike_womQTL_MAF10_full <- read.csv(paste0(wd$Zm_oillike_womQTL_MAF10_full, "var.QTN.genotypic.information.csv"))
zm_var.QTN.genotypic.information_oillike_womQTL_MAF40_500 <- read.csv(paste0(wd$Zm_oillike_womQTL_MAF40_500, "var.QTN.genotypic.information.csv"))
zm_var.QTN.genotypic.information_oillike_womQTL_MAF40_full <- read.csv(paste0(wd$Zm_oillike_womQTL_MAF40_full, "var.QTN.genotypic.information.csv"))

##molike settings
zm_var.QTN.genotypic.information_Molike_MAF40_500 <- read.csv(paste0(wd$Zm_molike_MAF40_500, "var.QTN.genotypic.information.csv"))
zm_var.QTN.genotypic.information_Molike_MAF40_full <- read.csv(paste0(wd$Zm_molike_MAF40_full, "var.QTN.genotypic.information.csv"))

#Molike_womQTL_MAF40
Zm_var.QTN.genotypic.information_Molike_womQTL_MAF40_sub500 <- read.csv(paste0(wd$Zm_molike_womQTL_MAF40_Sub500, "var.QTN.genotypic.information.csv"))
Zm_var.QTN.genotypic.information_Molike_womQTL_MAF40_FULL <- read.csv(paste0(wd$Zm_molike_womQTL_MAF40_FULL, "var.QTN.genotypic.information.csv"))


#Here, I have created a function that allows flexibility between A. thaliana and maize
##var.QTN.genotypic.information is the simulated vQTN information and window.size
##is how large your want your genomic window to be around your simulated vQTN. For example, the 
##window sizes used for Arabidopsis and maize was 100K and 250K, respectively. 
vQTN_window_maker <- function(var.QTN.genotypic.information = NULL, window.size = NULL) {
  var.QTN.genotypic.information <- var.QTN.genotypic.information[, 1:6]
  var.QTN.genotypic.information$QTNwindow_lowerbound <- var.QTN.genotypic.information[, 6] - window.size
  var.QTN.genotypic.information$QTNwindow_upperbound <- var.QTN.genotypic.information[, 6] + window.size
  return(var.QTN.genotypic.information)
}

#A thaliana
##oillike settings
at_vqgi_oillike_MAF10_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_oillike_momQTL_MAF10_500[, 1:5], window.size = 100000)
at_vqgi_oillike_MAF40_500  <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_oillike_momQTL_MAF40_500[, 1:5], window.size = 100000)
at_vqgi_oillike_MAF10_full  <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_oillike_momQTL_MAF10_FULL[, 1:5], window.size = 100000)
at_vqgi_oillike_MAF40_full  <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_oillike_momQTL_MAF40_FULL[, 1:5], window.size = 100000)

#molike settings
at_vqgi_molike_MAF40_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_MAF40_500[, 1:5], window.size = 100000)
at_vqgi_molike_MAF40_full <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_MAF40_full[, 1:5], window.size = 100000)
at_vqgi_molike_womQTL_MAF40_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_womQTL_MAF40_500[, 1:5], window.size = 100000)
at_vqgi_molike_womQTL_MAF40_full <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_womQTL_MAF40_full[, 1:5], window.size = 100000)

#maize
##oillike settings
zm_vqgi_oillike_womQTL_MAF10_500 <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_oillike_womQTL_MAF10_500[, 1:5], window.size = 250000)
zm_vqgi_oillike_womQTL_MAF10_full <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_oillike_womQTL_MAF10_full[, 1:5], window.size = 250000)
zm_vqgi_oillike_womQTL_MAF40_500 <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_oillike_womQTL_MAF40_500[, 1:5], window.size = 250000)
zm_vqgi_oillike_womQTL_MAF40_full <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_oillike_womQTL_MAF40_full[, 1:5], window.size = 250000)

#Molike settings
zm_vqgi_molike_MAF40_500 <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_Molike_MAF40_500[, 1:6], window.size = 250000)
zm_vqgi_molike_MAF40_full <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_Molike_MAF40_full[, 1:6], window.size = 250000)
Zm_vqgi_molike_womQTL_MAF40_500 <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_Molike_womQTL_MAF40_500[, 1:6], window.size = 250000)
Zm_vqgi_molike_womQTL_MAF40_full <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_Molike_womQTL_MAF40_full[, 1:6], window.size = 250000)

##Then, I will create shortcuts to result files
wd$at_bf_results <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/2_pipeline/0_load_data/out/vGWAS_chapter_one_At/vGWAS_Results/BF/"
wd$zm_bf_results <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/2_pipeline/0_load_data/out/vGWAS_chapter_one_Zm/vGWAS_Results/BF/"

wd$at_DGLM_results <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/2_pipeline/0_load_data/out/vGWAS_chapter_one_At/vGWAS_Results/DGLM/"
wd$zm_DGLM_results <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/2_pipeline/0_load_data/out/vGWAS_chapter_one_Zm/vGWAS_Results/DGLM/"

##Here, I will also create the switch to 3_output folder 
wd$output <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/3_output/"

#check to see if your folder has the 100 result files
at_bf_dirs <- list.dirs(wd$at_bf_results)
at_bf_dirs <- at_bf_dirs[c(3, 4, 6, 7, 9, 10, 12, 13, 15, 16)]

zm_bf_dirs <- list.dirs(wd$zm_bf_results)
zm_bf_dirs <- zm_bf_dirs[c(3, 4, 6, 7, 9, 10, 12, 13, 15, 16)]

at_DGLM_dirs <- list.dirs(wd$at_dglm_results)
at_DGLM_dirs <- at_DGLM_dirs[c(3, 4, 6, 7, 9, 10, 12, 13, 15, 16)]

zm_DGLM_dirs <- list.dirs(wd$zm_dglm_results)
zm_DGLM_dirs <- zm_DGLM_dirs[c(3, 4, 6, 7, 9, 10, 12, 13, 15, 16)]

#For arguments: Result file vec, QTN information both vQTN and simplePhenotypes, 
#output (make sure oil has multiple outputs)
#Here is what the arguments mean
##Results_dir: This refers to a predefined vector of results for a species and test combination
##rf: This refers to the number that in the results_dir vector
##Null results: This argument allows you to specify if you are processing a null setting or not
##simvQTN: This argument refers to the var.QTN information 
##mQTN.wd: This argument refers to the simplePHENOTYPE folder for additive effects. NULL if your setting does not have any additive effects. 
##DGLM: This argument specifies if you are processing DGLM (TRUE) or BFT (FALSE)
##Additive: TRUE for simulated additive QTNs or FALSE for no simulated additive QTNS
##trait.name: A character to what you want to label your output file as
##test: A character for BF_ or DGLM_

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
    vGWAS_rate_results <- matrix(NA, nrow=100, ncol = 9)
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

#To store all of my results, I will need to switch the working directory to my output folder
wd$output <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/3_output/"
setwd(wd$output)

################################ Detection Rates ############################################
#The first section will be Arabidopsis. For both Arabidopsis and maize, I will run Brown-Forsythe
#first then DGLM

#At
##BF
###Null
####500
Rate.calc(Result_dir = at_bf_dirs, rf = 9, simvQTN = NULL, mQTN.wd = NULL, DGLM = FALSE, Null_results = TRUE, trait.name = "at_null_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 10, simvQTN = NULL, mQTN.wd = NULL, DGLM = FALSE, Null_results = TRUE, trait.name = "at_null_full", test = "BF_", Additive = FALSE)

#Oillike MAF 10
##500 
rate.calc(Result_dir = at_bf_dirs, rf = 12, simvQTN = at_vqgi_oillike_MAF10_500, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "at_oillike_womQTL_MAF10_500", test = "BF_", Additive = FALSE)

##FULL
rate.calc(Result_dir = at_bf_dirs, rf = 13, simvQTN = at_vqgi_oillike_MAF10_full, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "at_oillike_womQTL_MAF10_full", test = "BF_", Additive = FALSE)

#Oillike MAF 40
#500
Rate.calc(Result_dir = at_bf_dirs, rf = 15, simvQTN = at_vqgi_oillike_MAF40_500, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "at_oillike_womQTL_MAF40_500", test = "BF_", Additive = FALSE)

##FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 16, simvQTN = at_vqgi_oillike_MAF40_full, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "at_oillike_womQTL_MAF40_full", test = "BF_", Additive = FALSE)

#Molike_MAF40
##500
Rate.calc(Result_dir = at_bf_dirs, rf = 3, simvQTN = at_vqgi_molike_MAF40_500, mQTN.wd = wd$SP_At_Molike_MAF40_500, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF40_500", test = "BF_", Additive = TRUE)

##FULL
Rate.calc(Result_dir = At_bf_dirs, rf = 4, simvQTN = at_vqgi_molike_MAF40_full, mQTN.wd = wd$SP_At_Molike_MAF40_FULL, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF40_full", test = "BF_", Additive = TRUE)

#Molike_womQTL_MAF40
##500
Rate.calc(Result_dir = at_bf_dirs, rf = 6, simvQTN = at_vqgi_molike_womQTL_MAF40_500, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_womQTL_MAF40_500", test = "BF_", Additive = FALSE)

##FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 7, simvQTN = at_vqgi_molike_womQTL_MAF40_full, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_womQTL_MAF40_full", test = "BF_", Additive = FALSE)

##DGLM Results
###Null
####500
Rate.calc(Result_dir = at_DGLM_dirs, rf = 9, simvQTN = NULL, mQTN.wd = NULL, DGLM = TRUE, Null_results = TRUE, trait.name = "at_null_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_DGLM_dirs, rf = 10, simvQTN = NULL, mQTN.wd = NULL, DGLM = TRUE, Null_results = TRUE, trait.name = "at_null_full", test = "DGLM_", Additive = FALSE)

###Oillike_mowQTL_MAF10
####500
Rate.calc(Result_dir = at_DGLM_dirs, rf = 12, simvQTN = at_vqgi_oillike_MAF10_500, mQTN.wd = NULL, DGLM = TRUE, Null_results = FALSE, trait.name = "at_oillike_womQTL_MAF10_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_DGLM_dirs, rf = 13, simvQTN = at_vqgi_oillike_MAF10_full, mQTN.wd = NULL, DGLM = TRUE, Null_results = FALSE, trait.name = "at_oillike_womQTL_MAF10_full", test = "DGLM_", Additive = FALSE)

###Oillike_mowQTL_MAF40
####Sub500
Rate.calc(Result_dir = at_DGLM_dirs, rf = 15, simvQTN = at_vqgi_oillike_MAF40_500, mQTN.wd = NULL, DGLM = TRUE, Null_results = FALSE, trait.name = "at_oillike_womQTL_MAF40_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_DGLM_dirs, rf = 16, simvQTN = at_vqgi_oillike_MAF40_full, mQTN.wd = NULL, DGLM = TRUE, Null_results = FALSE, trait.name = "at_oillike_womQTL_MAF10_full", test = "DGLM_", Additive = FALSE)

#Molike_MAF40
##Sub500
Rate.calc(Result_dir = at_DGLM_dirs, rf = 3, simvQTN = at_vqgi_molike_MAF40_500, mQTN.wd = wd$SP_At_Molike_MAF40_500, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF40_500", test = "DGLM_", Additive = TRUE)

##FULL
Rate.calc(Result_dir = at_DGLM_dirs, rf = 4, simvQTN = at_vqgi_molike_MAF40_full, mQTN.wd = wd$SP_At_Molike_MAF40_FULL, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF40_full", test = "DGLM_", Additive = TRUE)

#Molike_womQTL_MAF40
##500
Rate.calc(Result_dir = at_DGLM_dirs, rf = 6, simvQTN = at_vqgi_molike_womQTL_MAF40_500, mQTN.wd = NULL, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_womQTL_MAF40_500", test = "DGLM_", Additive = FALSE)

##FULL
Rate.calc(Result_dir = at_DGLM_dirs, rf = 7, simvQTN = at_vqgi_molike_womQTL_MAF40_full, mQTN.wd = NULL, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_womQTL_MAF40_full", test = "DGLM_", Additive = FALSE)

##### MAize #####
##Null
###500
Rate.calc(Result_dir = zm_bf_dirs, rf = 9, simvQTN = NULL, mQTN.wd = NULL, DGLM = FALSE, Null_results = TRUE, trait.name = "zm_null_500", test = "BF_", Additive = FALSE)

###FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 10, simvQTN = NULL, mQTN.wd = NULL, DGLM = FALSE, Null_results = TRUE, trait.name = "zm_null_full", test = "BF_", Additive = FALSE)

##Oillike_womQTL_MAF10
###500
Rate.calc(Result_dir = zm_bf_dirs, rf = 12, simvQTN = zm_vqgi_oillike_womQTL_MAF10_500, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_oillike_womQTL_MAF10_500", test = "BF_", Additive = FALSE)

###FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 13, simvQTN = zm_vqgi_oillikeFULL_womQTL, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_oillike_womQTL_MAF10_full", test = "BF_", Additive = FALSE)

#Oillike_womQTL_MAF40
##500 Diversity Panel
Rate.calc(Result_dir = zm_bf_dirs, rf = 15, simvQTN = zm_vqgi_oillike_womQTL_MAF40_500, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_oillike_womQTL_MAF40_500", test = "BF_", Additive = FALSE)

##FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 16, simvQTN = zm_vqgi_oillike_womQTL_MAF40_full, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_oillike_womQTL_MAF40_full", test = "BF_", Additive = FALSE)

#Molike MAF40 
##500
Rate.calc(Result_dir = zm_bf_dirs, rf = 3, simvQTN = zm_vqgi_molike_MAF40_500, mQTN.wd = wd$SP_Zm_Molike_MAF40_sub500, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_500", test = "BF_", Additive = TRUE)

##FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 4, simvQTN = Zm_vqgi_molike_MAF40_full, mQTN.wd = wd$SP_Zm_Molike_MAF40_FULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_full", test = "BF_", Additive = TRUE)

#Molike womQTL MAF40
##Sub500
Rate.calc(Result_dir = zm_bf_dirs, rf = 6, simvQTN = Zm_vqgi_molike_womQTL_MAF40_500, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_womQTL_MAF40_500", test = "BF_", Additive = FALSE)

##FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 7, simvQTN = zm_vqgi_molike_womQTL_MAF40_full, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_womQTL_MAF40_full", test = "BF_", Additive = FALSE)

#DGLM Results
##Null
###500
Rate.calc(Result_dir = zm_bf_dirs, rf = 9, simvQTN = NULL, mQTN.wd = NULL, DGLM = TRUE, Null_results = TRUE, trait.name = "zm_null_500", test = "DGLM_", Additive = FALSE)

###FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 10, simvQTN = NULL, mQTN.wd = NULL, DGLM = TRUE, Null_results = TRUE, trait.name = "zm_null_500", test = "DGLM_", Additive = FALSE)

##Oillike_womQTL_MAF10
###500
Rate.calc(Result_dir = zm_DGLM_dirs, rf = 12, simvQTN = zm_vqgi_oillike_womQTL_MAF10_500, mQTN.wd = NULL, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_oillike_womQTL_MAF10_500", test = "DGLM_", Additive = FALSE)

###FULL
Rate.calc(Result_dir = zm_DGLM_dirs, rf = 13, simvQTN = zm_vqgi_oillike_womQTL_MAF10_full, mQTN.wd = NULL, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_oillike_womQTL_MAF10_full", test = "DGLM_", Additive = FALSE)

##Oillike_womQTL_MAF40
###500
Rate.calc(Result_dir = zm_DGLM_dirs, rf = 15, simvQTN = zm_vqgi_oillike_womQTL_MAF40_500, mQTN.wd = NULL, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_oillike_womQTL_MAF40_500", test = "DGLM_", Additive = FALSE)

###FULL
Rate.calc(Result_dir = zm_DGLM_dirs, rf = 16, simvQTN = zm_vqgi_oillike_womQTL_MAF40_full, mQTN.wd = NULL, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_oillike_womQTL_MAF40_full", test = "DGLM_", Additive = FALSE)

#Molike MAF40
##500
Rate.calc(Result_dir = zm_DGLM_dirs, rf = 3, simvQTN = zm_vqgi_molike_MAF40_500, mQTN.wd = wd$SP_Zm_Molike_MAF40_sub500, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_Molike_MAF40_500", test = "DGLM_", Additive = TRUE)

##FULL
Rate.calc(Result_dir = zm_DGLM_dirs, rf = 4, simvQTN = zm_vqgi_molike_MAF40_full, mQTN.wd = wd$SP_Zm_Molike_MAF40_FULL, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_molike_MAF40_full", test = "DGLM_", Additive = TRUE)

#Molike womQTL MAF40
##500
Rate.calc(Result_dir = zm_DGLM_dirs, rf = 6, simvQTN = zm_vqgi_molike_womQTL_MAF40_500, mQTN.wd = NULL, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_molike_womQTL_MAF40_500", test = "DGLM_", Additive = FALSE)

##FULL
Rate.calc(Result_dir = zm_DGLM_dirs, rf = 7, simvQTN = zm_vqgi_molike_womQTL_MAF40_full, mQTN.wd = NULL, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_molike_womQTL_MAF40_full", test = "DGLM_", Additive = FALSE)
################################End########################################