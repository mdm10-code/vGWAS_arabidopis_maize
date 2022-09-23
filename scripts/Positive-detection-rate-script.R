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
wd$vGWAS_chapter_one_at_revisions <- paste0(wd$vGWAS_chapter_one_at, "Revisions/")
wd$at_vGWAS_results <- paste0(wd$vGWAS_chapter_one_at_revisions, "vGWAS_results/")
wd$at_bf_results <- paste0(wd$at_vGWAS_results, "BF/")
wd$at_dglm_results <- paste0(wd$at_vGWAS_results, "DGLM/")
wd$at_dglm_results_revisions <- paste0(wd$at_dglm_results, "Revisions/")

#This is for maize
wd$vGWAS_chapter_one_zm <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_maize_filtlike_/"
wd$zm_vGWAS_results <- paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/")
wd$zm_bf_results <- paste0(wd$zm_vGWAS_results, "BF/")
wd$zm_dglm_results <- paste0(wd$zm_vGWAS_results, "DGLM/")

#Maize Revisions
wd$vGWAS_chapter_one_zm <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/"
wd$zm_vGWAS_results <- paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/")
wd$zm_bft_results <- paste0(wd$zm_vGWAS_results, "BFT/")
wd$zm_dglm_results <- paste0(wd$vGWAS_chapter_one_zm_filtlike_at, "DGLM/")


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

wd$at_vpheno_molike_MAF40_h233_vQTN10_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h233_vQTN10_500/")
wd$at_vpheno_molike_MAF40_h233_vQTN10_full <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h233_vQTN10_full/")
wd$at_vpheno_molike_MAF40_h233_vQTN50_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h233_vQTN50_500/")
wd$at_vpheno_molike_MAF40_h233_vQTN50_full <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h233_vQTN50_full/")
wd$at_vpheno_molike_MAF40_h233_vQTN90_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h233_vQTN90_500/")
wd$at_vpheno_molike_MAF40_h233_vQTN90_full <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h233_vQTN90_full/")

#Molike with h2 = 0.63
wd$at_vpheno_molike_MAF40_h263_vQTN10_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h263_vQTN10_500/")
wd$at_vpheno_molike_MAF40_h263_vQTN10_full <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h263_vQTN10_full/")
wd$at_vpheno_molike_MAF40_h263_vQTN50_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h263_vQTN50_500/")
wd$at_vpheno_molike_MAF40_h263_vQTN50_full <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h263_vQTN50_full/")

##molike MAF10 h233
wd$at_vpheno_molike_MAF10_h233_vQTN10_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF10_h233_vQTN10_500/")
wd$at_vpheno_molike_MAF10_h233_vQTN10_full <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF10_h233_vQTN10_full/")
wd$at_vpheno_molike_MAF10_h233_vQTN50_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF10_h233_vQTN50_500/")
wd$at_vpheno_molike_MAF10_h233_vQTN50_full <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF10_h233_vQTN50_full/")
wd$at_vpheno_molike_MAF10_h233_vQTN90_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF10_h233_vQTN90_500/")
wd$at_vpheno_molike_MAF10_h233_vQTN90_full <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF10_h233_vQTN90_full/")

##molike MAF10 h263
wd$at_vpheno_molike_MAF10_h263_vQTN10_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF10_h263_vQTN10_500/")
wd$at_vpheno_molike_MAF10_h263_vQTN10_full <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF10_h263_vQTN10_full/")
wd$at_vpheno_molike_MAF10_h263_vQTN50_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF10_h263_vQTN50_500/")
wd$at_vpheno_molike_MAF10_h263_vQTN50_full <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF10_h263_vQTN50_full/")
wd$at_vpheno_molike_MAF10_h263_vQTN90_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF10_h263_vQTN90_500/")
wd$at_vpheno_molike_MAF10_h263_vQTN90_full <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF10_h263_vQTN90_full/")

##Molike MAF 40
###h233
wd$at_vpheno_molike_MAF40_h233_vQTN10_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h233_vQTN10_500/")
wd$at_vpheno_molike_MAF40_h233_vQTN10_full <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h233_vQTN10_full/")
wd$at_vpheno_molike_MAF40_h233_vQTN50_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h233_vQTN50_500/")
wd$at_vpheno_molike_MAF40_h233_vQTN50_full <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h233_vQTN50_full/")
wd$at_vpheno_molike_MAF40_h233_vQTN90_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h233_vQTN90_500/")
wd$at_vpheno_molike_MAF40_h233_vQTN90_full <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h233_vQTN90_full/")

##h263
wd$at_vpheno_molike_MAF40_h263_vQTN10_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h263_vQTN10_500/")
wd$at_vpheno_molike_MAF40_h263_vQTN10_full <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h263_vQTN10_full/")
wd$at_vpheno_molike_MAF40_h263_vQTN50_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h263_vQTN50_500/")
wd$at_vpheno_molike_MAF40_h263_vQTN50_full <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF40_h263_vQTN50_full/")


##molike MAF5 h263
wd$at_vpheno_molike_MAF5_h263_vQTN90_500 <- paste0(wd$at_vQTL_simulations, "at_vpheno_molike_MAF5_h263_vQTN90_500/")

###maize
wd$zm_vQTL_simulations <- paste0(wd$vGWAS_chapter_one_zm, "zm_vQTL_simulations/")

##These are for the new simulations
wd$zm_molike_MAF40_500 <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF40_500/")
wd$zm_molike_MAF40_full <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF40_full/")

###MAF10
####h233
wd$zm_molike_MAF10_h233_vQTN10_500 <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF10_h233_vQTN10_500/")
wd$zm_molike_MAF10_h233_vQTN10_FULL <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF10_h233_vQTN10_FULL/")
wd$zm_molike_MAF10_h233_vQTN50_500 <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF10_h233_vQTN50_500/")
wd$zm_molike_MAF10_h233_vQTN50_FULL <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF10_h233_vQTN50_FULL/")
wd$zm_molike_MAF10_h233_vQTN90_500 <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF10_h233_vQTN90_500/")
wd$zm_molike_MAF10_h233_vQTN90_FULL <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF10_h233_vQTN90_FULL/")

####h263
wd$zm_molike_MAF10_h263_vQTN10_500 <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF10_h263_vQTN10_500/")
wd$zm_molike_MAF10_h263_vQTN10_FULL <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF10_h263_vQTN10_FULL/")
wd$zm_molike_MAF10_h263_vQTN50_500 <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF10_h263_vQTN50_500/")
wd$zm_molike_MAF10_h263_vQTN50_FULL <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF10_h263_vQTN50_FULL/")
wd$zm_molike_MAF10_h263_vQTN90_500 <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF10_h263_vQTN90_500/")
wd$zm_molike_MAF10_h263_vQTN90_FULL <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF10_h263_vQTN90_FULL/")

##MAF40
###h233
wd$zm_molike_MAF40_h233_vQTN10_500 <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF40_h233_vQTN10_500_filtlike_at/")
wd$zm_molike_MAF40_h233_vQTN10_FULL <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF40_h233_vQTN10_FULL_filtlike_at/")
wd$zm_molike_MAF40_h233_vQTN50_500 <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF40_h233_vQTN50_500/")
wd$zm_molike_MAF40_h233_vQTN50_FULL <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF40_h233_vQTN50_FULL/")
wd$zm_molike_MAF40_h233_vQTN90_500 <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF40_h233_vQTN90_500/")
wd$zm_molike_MAF40_h233_vQTN90_FULL <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF40_h233_vQTN90_FULL/")

####h263
wd$zm_molike_MAF40_h263_vQTN10_500 <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF40_h263_vQTN10_500/")
wd$zm_molike_MAF40_h263_vQTN10_FULL <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF40_h263_vQTN10_FULL/")
wd$zm_molike_MAF40_h263_vQTN50_500 <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF40_h263_vQTN50_500/")
wd$zm_molike_MAF40_h263_vQTN50_FULL <- paste0(wd$zm_vQTL_simulations, "zm_vmolike_MAF40_h263_vQTN50_FULL/")


#We also need the selected additive QTN from simplePHENOTYPES
#A thaliana
wd$at_sp_simulations <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_at/at_sp_results/"

wd$sp_at_null <- paste0(wd$at_sp_simulations, "sp_at_null/")
wd$sp_at_null_500 <- paste0(wd$at_sp_simulations, "sp_at_null_500/")
wd$sp_at_null_full <- paste0(wd$at_sp_simulations, "sp_at_null_full/")
wd$sp_at_molike_MAF40_500 <- paste0(wd$at_sp_simulations, "sp_at_molike_wmQTL_500/")
wd$sp_at_molike_MAF40_full <- paste0(wd$at_sp_simulations, "sp_at_molike_wmQTL_full/")

#For epistasis

#500
wd$sp_at_epistasis_lowefs_500 <- paste0(wd$at_sp_results, "sp_at_epi_loweffectsize_500/")
wd$sp_at_epistasis_highefs_500 <- paste0(wd$at_sp_results, "sp_at_epi_higheffectsize_500/")

#FULL
wd$sp_at_epistasis_lowefs_FULL <- paste0(wd$at_sp_results, "sp_at_epi_loweffectsize_full/")
wd$sp_at_epistasis_highefs_FULL <- paste0(wd$at_sp_results, "sp_at_epi_higheffectsize_full/")

##GxE nocor
wd$sp_at_GxE_500 <- paste0(wd$at_sp_results, "sp_at_GxE_500_nocor/")
wd$sp_at_GxE_FULL <- paste0(wd$at_sp_results, "sp_at_GxE_FULL_nocor/")

##Epistasis lowh2 effect sizes
#500
wd$sp_at_epistasis_lowh2_lowefs_500 <- paste0(wd$at_sp_results, "sp_at_epi_loweffectsize_h230_500/")
wd$sp_at_epistasis_lowh2_highefs_500 <- paste0(wd$at_sp_results, "sp_at_epi_higheffectsize_h230_500/")

#FULL
wd$sp_at_epistasis_lowh2_lowefs_FULL <- paste0(wd$at_sp_results, "sp_at_epi_loweffectsize_h230_full/")
wd$sp_at_epistasis_lowh2_highefs_FULL <- paste0(wd$at_sp_results, "sp_at_epi_higheffectsize_h230_full/")


#Maize
wd$zm_sp_results <- paste0(wd$vGWAS_chapter_one_zm, "zm_sp_results/")

wd$sp_zm_null <- paste0(wd$zm_sp_results, "sp_zm_null_full/")
wd$sp_zm_null_500 <- paste0(wd$zm_sp_results, "sp_zm_null_500/")
wd$sp_zm_molike_MAF40_500 <- paste0(wd$zm_sp_results, "sp_zm_molike_wmQTL_500/")
wd$sp_zm_molike_MAF40_full <- paste0(wd$zm_sp_results, "sp_zm_molike_wmQTL_full/")

wd$sp_zm_molike_wmQTL_500_filtlike_at <- paste0(wd$zm_sp_results, "sp_zm_molike_wmQTL_500_filtlike_at/")
wd$sp_zm_molike_wmQTL_FULL_filtlike_at <- paste0(wd$zm_sp_results, "sp_zm_molike_wmQTL_FULL_filtlike_at/")

##NULL filtered like at
wd$sp_zm_null_500_filt_like_at <- paste0(wd$zm_sp_simulations, "sp_zm_null_500_filt_like_at/")
wd$sp_zm_null_full_filt_like_at <- paste0(wd$zm_sp_simulations, "sp_zm_null_full_filt_like_at/")


#Epistasis
#For epistasis
##500
wd$sp_zm_epistasis_lowh2_500 <- paste0(wd$zm_sp_results, "sp_zm_epi_higheffectsize_h230_500/")
wd$sp_zm_epistasis_highh2_500 <- paste0(wd$zm_sp_results, "sp_zm_epi_higheffectsize_500/")

#FULL
wd$sp_zm_epistasis_lowh2_FULL <- paste0(wd$zm_sp_results, "sp_zm_epi_higheffectsize_h230_full/")
wd$sp_zm_epistasis_highh2_FULL <- paste0(wd$zm_sp_results, "sp_zm_epi_higheffectsize_full/")

##gXe
wd$sp_zm_GxE_500 <- paste0(wd$zm_sp_results, "sp_zm_GxE_500/")
wd$sp_zm_GxE_FULL <- paste0(wd$zm_sp_results, "sp_zm_GxE_FULL/")


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

##The Molike MAF40 h233 vQTN 0.10
at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h233_vQTN10_500 <- read.csv(paste0(wd$at_vpheno_molike_MAF40_h233_vQTN10_500, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h233_vQTN10_full <- read.csv(paste0(wd$at_vpheno_molike_MAF40_h233_vQTN10_full, "var.QTN.genotypic.information.csv"))

###vQTN 0.5
at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h233_vQTN50_500 <- read.csv(paste0(wd$at_vpheno_molike_MAF40_h233_vQTN50_500, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h233_vQTN50_full <- read.csv(paste0(wd$at_vpheno_molike_MAF40_h233_vQTN50_full, "var.QTN.genotypic.information.csv"))

###vQTN 0.9
at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h233_vQTN90_500 <- read.csv(paste0(wd$at_vpheno_molike_MAF40_h233_vQTN90_500, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h233_vQTN90_full <- read.csv(paste0(wd$at_vpheno_molike_MAF40_h233_vQTN90_full, "var.QTN.genotypic.information.csv"))

##The molike MAF40 h263 settings
###vQTN = 0.1
at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h263_vQTN10_500 <- read.csv(paste0(wd$at_vpheno_molike_MAF40_h263_vQTN10_500, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h263_vQTN10_full <- read.csv(paste0(wd$at_vpheno_molike_MAF40_h263_vQTN10_full, "var.QTN.genotypic.information.csv"))

###vQTN = 0.5
at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h263_vQTN50_500 <- read.csv(paste0(wd$at_vpheno_molike_MAF40_h263_vQTN50_500, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h263_vQTN50_full <- read.csv(paste0(wd$at_vpheno_molike_MAF40_h263_vQTN50_full, "var.QTN.genotypic.information.csv"))

##The molike MAF10 h233 settings
###vQTN = 0.1
at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h233_vQTN10_500 <- read.csv(paste0(wd$at_vpheno_molike_MAF10_h233_vQTN10_500, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h233_vQTN10_full <- read.csv(paste0(wd$at_vpheno_molike_MAF10_h233_vQTN10_full, "var.QTN.genotypic.information.csv"))

###vQTN = 0.5
at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h233_vQTN50_500 <- read.csv(paste0(wd$at_vpheno_molike_MAF10_h233_vQTN50_500, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h233_vQTN50_full <- read.csv(paste0(wd$at_vpheno_molike_MAF10_h233_vQTN50_full, "var.QTN.genotypic.information.csv"))

###vQTN = 0.9
at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h233_vQTN90_500 <- read.csv(paste0(wd$at_vpheno_molike_MAF10_h233_vQTN90_500, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h233_vQTN90_full <- read.csv(paste0(wd$at_vpheno_molike_MAF10_h233_vQTN90_full, "var.QTN.genotypic.information.csv"))

##MAF10 h263
###vQTN = 0.1
at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h263_vQTN10_500 <- read.csv(paste0(wd$at_vpheno_molike_MAF10_h263_vQTN10_500, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h263_vQTN10_full <- read.csv(paste0(wd$at_vpheno_molike_MAF10_h263_vQTN10_full, "var.QTN.genotypic.information.csv"))

###vQTN = 0.5
at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h263_vQTN50_500 <- read.csv(paste0(wd$at_vpheno_molike_MAF10_h263_vQTN50_500, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h263_vQTN50_full <- read.csv(paste0(wd$at_vpheno_molike_MAF10_h263_vQTN50_full, "var.QTN.genotypic.information.csv"))

###vQTN = 0.9
at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h263_vQTN90_500 <- read.csv(paste0(wd$at_vpheno_molike_MAF10_h263_vQTN90_500, "var.QTN.genotypic.information.csv"))
at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h263_vQTN90_full <- read.csv(paste0(wd$at_vpheno_molike_MAF10_h263_vQTN90_full, "var.QTN.genotypic.information.csv"))

#Maize

##molike settings
zm_var.QTN.genotypic.information_Molike_MAF40_500 <- read.csv(paste0(wd$zm_molike_MAF40_500, "var.QTN.genotypic.information.csv"))
zm_var.QTN.genotypic.information_Molike_MAF40_full <- read.csv(paste0(wd$zm_molike_MAF40_full, "var.QTN.genotypic.information.csv"))

#######Rest of the molike settings
##The Molike MAF40 h233 vQTN 0.10
zm_var.QTN.genotypic.information_molike_MAF40_h233_vQTN10_500 <- read.csv(paste0(wd$zm_molike_MAF40_h233_vQTN10_500, "var.QTN.genotypic.information.csv"))
zm_var.QTN.genotypic.information_molike_MAF40_h233_vQTN10_full <- read.csv(paste0(wd$zm_molike_MAF40_h233_vQTN10_FULL, "var.QTN.genotypic.information.csv"))

###vQTN 0.5
zm_var.QTN.genotypic.information_molike_MAF40_h233_vQTN50_500 <- read.csv(paste0(wd$zm_molike_MAF40_h233_vQTN50_500, "var.QTN.genotypic.information.csv"))
zm_var.QTN.genotypic.information_molike_MAF40_h233_vQTN50_full <- read.csv(paste0(wd$zm_molike_MAF40_h233_vQTN50_FULL, "var.QTN.genotypic.information.csv"))

###vQTN 0.9
zm_var.QTN.genotypic.information_molike_MAF40_h233_vQTN90_500 <- read.csv(paste0(wd$zm_molike_MAF40_h233_vQTN90_500, "var.QTN.genotypic.information.csv"))
zm_var.QTN.genotypic.information_molike_MAF40_h233_vQTN90_full <- read.csv(paste0(wd$zm_molike_MAF40_h233_vQTN90_FULL, "var.QTN.genotypic.information.csv"))

##The molike MAF40 h263 settings
###vQTN = 0.1
zm_var.QTN.genotypic.information_molike_MAF40_h263_vQTN10_500 <- read.csv(paste0(wd$zm_molike_MAF40_h263_vQTN10_500, "var.QTN.genotypic.information.csv"))
zm_var.QTN.genotypic.information_molike_MAF40_h263_vQTN10_full <- read.csv(paste0(wd$zm_molike_MAF40_h263_vQTN10_FULL, "var.QTN.genotypic.information.csv"))

###vQTN = 0.5
zm_var.QTN.genotypic.information_molike_MAF40_h263_vQTN50_500 <- read.csv(paste0(wd$zm_molike_MAF40_h263_vQTN50_500, "var.QTN.genotypic.information.csv"))
zm_var.QTN.genotypic.information_molike_MAF40_h263_vQTN50_full <- read.csv(paste0(wd$zm_molike_MAF40_h263_vQTN50_FULL, "var.QTN.genotypic.information.csv"))

##The molike MAF10 h233 settings
###vQTN = 0.1
zm_var.QTN.genotypic.information_molike_MAF10_h233_vQTN10_500 <- read.csv(paste0(wd$zm_molike_MAF10_h233_vQTN10_500, "var.QTN.genotypic.information.csv"))
zm_var.QTN.genotypic.information_molike_MAF10_h233_vQTN10_full <- read.csv(paste0(wd$zm_molike_MAF10_h233_vQTN10_FULL, "var.QTN.genotypic.information.csv"))

###vQTN = 0.5
zm_var.QTN.genotypic.information_molike_MAF10_h233_vQTN50_500 <- read.csv(paste0(wd$zm_molike_MAF10_h233_vQTN50_500, "var.QTN.genotypic.information.csv"))
zm_var.QTN.genotypic.information_molike_MAF10_h233_vQTN50_full <- read.csv(paste0(wd$zm_molike_MAF10_h233_vQTN50_FULL, "var.QTN.genotypic.information.csv"))

###vQTN = 0.9
zm_var.QTN.genotypic.information_molike_MAF10_h233_vQTN90_500 <- read.csv(paste0(wd$zm_molike_MAF10_h233_vQTN90_500, "var.QTN.genotypic.information.csv"))
zm_var.QTN.genotypic.information_molike_MAF10_h233_vQTN90_full <- read.csv(paste0(wd$zm_molike_MAF10_h233_vQTN90_FULL, "var.QTN.genotypic.information.csv"))

##MAF10 h263
###vQTN = 0.1
zm_var.QTN.genotypic.information_molike_MAF10_h263_vQTN10_500 <- read.csv(paste0(wd$zm_molike_MAF10_h263_vQTN10_500, "var.QTN.genotypic.information.csv"))
zm_var.QTN.genotypic.information_molike_MAF10_h263_vQTN10_full <- read.csv(paste0(wd$zm_molike_MAF10_h263_vQTN10_FULL, "var.QTN.genotypic.information.csv"))

###vQTN = 0.5
zm_var.QTN.genotypic.information_molike_MAF10_h263_vQTN50_500 <- read.csv(paste0(wd$zm_molike_MAF10_h263_vQTN50_500, "var.QTN.genotypic.information.csv"))
zm_var.QTN.genotypic.information_molike_MAF10_h263_vQTN50_full <- read.csv(paste0(wd$zm_molike_MAF10_h263_vQTN50_FULL, "var.QTN.genotypic.information.csv"))

###vQTN = 0.9
zm_var.QTN.genotypic.information_molike_MAF10_h263_vQTN90_500 <- read.csv(paste0(wd$zm_molike_MAF10_h263_vQTN90_500, "var.QTN.genotypic.information.csv"))
zm_var.QTN.genotypic.information_molike_MAF10_h263_vQTN90_full <- read.csv(paste0(wd$zm_molike_MAF10_h263_vQTN90_FULL, "var.QTN.genotypic.information.csv"))


#Here, I have created a function that allows flexibility between A. thaliana and maize
##var.QTN.genotypic.information is the simulated vQTN information and window.size
##is how large your want your genomic window to be around your simulated vQTN. For example, the 
##window sizes used for Arabidopsis and maize was 100K and 250K, respectively. 
vQTN_window_maker <- function(var.QTN.genotypic.information = NULL, window.size = NULL) {
  var.QTN.genotypic.information <- var.QTN.genotypic.information[, 1:5]
  var.QTN.genotypic.information$QTNwindow_lowerbound <- var.QTN.genotypic.information$pos - window.size
  var.QTN.genotypic.information$QTNwindow_upperbound <- var.QTN.genotypic.information$pos + window.size
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

#molike h233 MAF40 vQTN 10
at_vqgi_molike_wmQTL_MAF40_h233_vQTN10_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h233_vQTN10_500[, 1:6], window.size = 100000)
at_vqgi_molike_wmQTL_MAF40_h233_vQTN10_full <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h233_vQTN10_full[, 1:6], window.size = 100000)

#molike h233 MAF40 vQTN 50
at_vqgi_molike_wmQTL_MAF40_h233_vQTN50_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h233_vQTN50_500[, 1:6], window.size = 100000)
at_vqgi_molike_wmQTL_MAF40_h233_vQTN50_full <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h233_vQTN50_full[, 1:6], window.size = 100000)

#molike h233 MAF40 vQTN 90
at_vqgi_molike_wmQTL_MAF40_h233_vQTN90_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h233_vQTN90_500[, 1:6], window.size = 100000)
at_vqgi_molike_wmQTL_MAF40_h233_vQTN90_full <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h233_vQTN90_full[, 1:6], window.size = 100000)

#molike h263 MAF40 vQTN10
at_vqgi_molike_wmQTL_MAF40_h263_vQTN10_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h263_vQTN10_500[, 1:6], window.size = 100000)
at_vqgi_molike_wmQTL_MAF40_h263_vQTN10_full <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h263_vQTN10_full[, 1:6], window.size = 100000)

#molike h263 MAF40 vQTN 50
at_vqgi_molike_wmQTL_MAF40_h263_vQTN50_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h263_vQTN50_500[, 1:6], window.size = 100000)
at_vqgi_molike_wmQTL_MAF40_h263_vQTN50_full <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF40_h263_vQTN50_full[, 1:6], window.size = 100000)

#molike h233 MAF10 vQTN 10
at_vqgi_molike_wmQTL_MAF10_h233_vQTN10_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h233_vQTN10_500[, 1:6], window.size = 100000)
at_vqgi_molike_wmQTL_MAF10_h233_vQTN10_full <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h233_vQTN10_full[, 1:6], window.size = 100000)

#molike h233 MAF10 vQTN 50
at_vqgi_molike_wmQTL_MAF10_h233_vQTN50_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h233_vQTN50_500[, 1:6], window.size = 100000)
at_vqgi_molike_wmQTL_MAF10_h233_vQTN50_full <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h233_vQTN50_full[, 1:6], window.size = 100000)

#molike h233 MAF10 vQTN90
at_vqgi_molike_wmQTL_MAF10_h233_vQTN90_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h233_vQTN90_500[, 1:6], window.size = 100000)
at_vqgi_molike_wmQTL_MAF10_h233_vQTN90_full <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h233_vQTN90_full[, 1:6], window.size = 100000)

##Molike h263 MAF10
#molike h263 MAF10 vQTN 10
at_vqgi_molike_wmQTL_MAF10_h263_vQTN10_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h263_vQTN10_500[, 1:6], window.size = 100000)
at_vqgi_molike_wmQTL_MAF10_h263_vQTN10_full <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h263_vQTN10_full[, 1:6], window.size = 100000)

#molike h263 MAF10 vQTN 50
at_vqgi_molike_wmQTL_MAF10_h263_vQTN50_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h263_vQTN50_500[, 1:6], window.size = 100000)
at_vqgi_molike_wmQTL_MAF10_h263_vQTN50_full <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h263_vQTN50_full[, 1:6], window.size = 100000)

#molike h263 MAF10 vQTN90
at_vqgi_molike_wmQTL_MAF10_h263_vQTN90_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h263_vQTN90_500[, 1:6], window.size = 100000)
at_vqgi_molike_wmQTL_MAF10_h263_vQTN90_full <- vQTN_window_maker(var.QTN.genotypic.information = at_var.QTN.genotypic.information_molike_wmQTL_MAF10_h263_vQTN90_full[, 1:6], window.size = 100000)

#For The epistasis settings, I am going to have to do some data wrangling
#First, I will need to manipulate the simplePHENOTYPES information to match
#my functions. I will need to read in that information first.

###Heritability
#500
at_epi.QTN.genotypic.information_lowefs_500 <- read.delim(paste0(wd$sp_at_epistasis_lowefs_500, "Epistatic_QTNs.txt"))
at_epi.QTN.genotypic.information_highefs_500 <- read.delim(paste0(wd$sp_at_epistasis_highefs_500, "Epistatic_QTNs.txt"))

at_epi_lowefs_500_forwindow <- at_epi.QTN.genotypic.information_lowefs_500[, 3:8]
at_epi_highefs_500_forwindow <- at_epi.QTN.genotypic.information_highefs_500[, 3:8]


at_eqgi_epistasis_lowefs_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_epi_lowefs_500_forwindow[, 1:6], window.size = 100000)
at_eqgi_epistasis_highefs_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_epi_highefs_500_forwindow[, 1:6], window.size = 100000)


#FULL
at_epi.QTN.genotypic.information_lowefs_full <- read.delim(paste0(wd$sp_at_epistasis_lowefs_FULL, "Epistatic_QTNs.txt"))
at_epi.QTN.genotypic.information_highefs_full <- read.delim(paste0(wd$sp_at_epistasis_highefs_FULL, "Epistatic_QTNs.txt"))

at_epi_lowefs_full_forwindow <- at_epi.QTN.genotypic.information_lowefs_full[, 3:8]
at_epi_highefs_full_forwindow <- at_epi.QTN.genotypic.information_highefs_full[, 3:8]


at_eqgi_epistasis_lowefs_FULL <- vQTN_window_maker(var.QTN.genotypic.information = at_epi_lowefs_full_forwindow[, 1:6], window.size = 100000)
at_eqgi_epistasis_highefs_FULL <- vQTN_window_maker(var.QTN.genotypic.information = at_epi_highefs_full_forwindow[, 1:6], window.size = 100000)

###This is for lowh2 effect sizes
at_epi.QTN.genotypic.information_h230_lowefs_500 <- read.delim(paste0(wd$sp_at_epistasis_lowh2_lowefs_500, "Epistatic_QTNs.txt"))
at_epi.QTN.genotypic.information_h230_highefs_500 <- read.delim(paste0(wd$sp_at_epistasis_lowh2_highefs_500, "Epistatic_QTNs.txt"))

at_epi_lowh2_lowefs_500_forwindow <- at_epi.QTN.genotypic.information_h230_lowefs_500[, 3:8]
at_epi_lowh2_highefs_500_forwindow <- at_epi.QTN.genotypic.information_h230_highefs_500[, 3:8]


at_eqgi_epistasis_lowh2_lowefs_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_epi_lowh2_lowefs_500_forwindow[, 1:6], window.size = 100000)
at_eqgi_epistasis_lowh2_highefs_500 <- vQTN_window_maker(var.QTN.genotypic.information = at_epi_lowh2_highefs_500_forwindow[, 1:6], window.size = 100000)

###This is for lowh2 effect sizes FULL
at_epi.QTN.genotypic.information_h230_lowefs_FULL <- read.delim(paste0(wd$sp_at_epistasis_lowh2_lowefs_FULL, "Epistatic_QTNs.txt"))
at_epi.QTN.genotypic.information_h230_highefs_FULL <- read.delim(paste0(wd$sp_at_epistasis_lowh2_highefs_FULL, "Epistatic_QTNs.txt"))

at_epi_lowh2_lowefs_FULL_forwindow <- at_epi.QTN.genotypic.information_h230_lowefs_FULL[, 3:8]
at_epi_lowh2_highefs_FULL_forwindow <- at_epi.QTN.genotypic.information_h230_highefs_FULL[, 3:8]

at_eqgi_epistasis_lowh2_lowefs_FULL <- vQTN_window_maker(var.QTN.genotypic.information = at_epi_lowh2_lowefs_FULL_forwindow[, 1:6], window.size = 100000)
at_eqgi_epistasis_lowh2_highefs_FULL <- vQTN_window_maker(var.QTN.genotypic.information = at_epi_lowh2_highefs_FULL_forwindow[, 1:6], window.size = 100000)



##I will also need to do the same for GxE
at_gxe.QTN.genotypic.information_GxE_500 <- read.delim(paste0(wd$sp_at_GxE_500, "Additive_QTNs.txt"))
at_gxe.QTN.genotypic.information_GxE_FULL <- read.delim(paste0(wd$sp_at_GxE_FULL, "Additive_QTNs.txt"))

at_gxe_diff_500_forwindow <- at_gxe.QTN.genotypic.information_GxE_500[c(1, 2, 4), 4:9]
at_gxe_diff_full_forwindow <- at_gxe.QTN.genotypic.information_GxE_FULL[c(1, 2, 4), 4:9]

##Maize
zm_gxe.QTN.genotypic.information_GxE_500 <- read.delim(paste0(wd$sp_zm_GxE_500, "Additive_Selected_QTNs.txt"))
zm_gxe.QTN.genotypic.information_GxE_FULL <- read.delim(paste0(wd$sp_zm_GxE_FULL, "Additive_Selected_QTNs.txt"))

zm_gxe_diff_500_forwindow <- zm_gxe.QTN.genotypic.information_GxE_500[c(1, 2, 4), 4:9]
zm_gxe_diff_full_forwindow <- zm_gxe.QTN.genotypic.information_GxE_FULL[c(1, 2, 4), 4:9]


#Now I repeat the same for maize
##
zm_epi.QTN.genotypic.information_lowh2_500 <- read.delim(paste0(wd$sp_zm_epistasis_lowh2_500, "Epistatic_Selected_QTNs.txt"))
zm_epi.QTN.genotypic.information_highh2_500 <- read.delim(paste0(wd$sp_zm_epistasis_highh2_500, "Epistatic_Selected_QTNs.txt"))

zm_epi.QTN.genotypic.information_lowh2_full <- read.delim(paste0(wd$sp_zm_epistasis_lowh2_FULL, "Epistatic_Selected_QTNs.txt"))
zm_epi.QTN.genotypic.information_highh2_full <- read.delim(paste0(wd$sp_zm_epistasis_highh2_FULL, "Epistatic_Selected_QTNs.txt"))

#500
zm_epi_lowh2_500_forwindow <- zm_epi.QTN.genotypic.information_lowh2_500[, 3:8]
zm_epi_highh2_500_forwindow <- zm_epi.QTN.genotypic.information_highh2_500[, 3:8]

#FULL
zm_epi_lowh2_full_forwindow <- zm_epi.QTN.genotypic.information_lowh2_full[, 3:8]
zm_epi_highh2_full_forwindow <- zm_epi.QTN.genotypic.information_highh2_full[, 3:8]

###Now, repeat with 1MB
#500
zm_eqgi_epistasis_lowh2_500_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_epi_lowh2_500_forwindow[, 1:6], window.size = 250000)
zm_eqgi_epistasis_highh2_500_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_epi_highh2_500_forwindow[, 1:6], window.size = 250000)

#FULL
zm_eqgi_epistasis_lowh2_FULL_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_epi_lowh2_full_forwindow[, 1:6], window.size = 250000)
zm_eqgi_epistasis_highh2_FULL_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_epi_highh2_full_forwindow[, 1:6], window.size = 250000)



##GxE
at_gxeqgi_diff_500_100K <- vQTN_window_maker(var.QTN.genotypic.information = at_gxe_diff_500_forwindow[, 1:6], window.size = 100000)
at_gxeqgi_diff_FULL_100K <- vQTN_window_maker(var.QTN.genotypic.information = at_gxe_diff_full_forwindow[, 1:6], window.size = 100000)

#Molike filtered like Arabidopsis
zm_vqgi_molike_MAF40_500 <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_Molike_MAF40_500[, 1:6], window.size = 250000)
zm_vqgi_molike_MAF40_full <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_Molike_MAF40_full[, 1:6], window.size = 250000)

##GxE

#Repeating above for maize for 1MB windows
zm_gxeqgi_diff_500_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_gxe_diff_500_forwindow[, 1:6], window.size = 250000)
zm_gxeqgi_diff_FULL_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_gxe_diff_full_forwindow[, 1:6], window.size = 250000)

#molike h233 MAF40 vQTN 10
zm_vqgi_molike_MAF40_h233_vQTN10_500_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF40_h233_vQTN10_500[, 1:6], window.size = 250000)
zm_vqgi_molike_MAF40_h233_vQTN10_full_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF40_h233_vQTN10_full[, 1:6], window.size = 250000)

#molike h233 MAF40 vQTN 50
zm_vqgi_molike_MAF40_h233_vQTN50_500_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF40_h233_vQTN50_500[, 1:6], window.size = 250000)
zm_vqgi_molike_MAF40_h233_vQTN50_full_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF40_h233_vQTN50_full[, 1:6], window.size = 250000)

#molike h233 MAF40 vQTN 90
zm_vqgi_molike_MAF40_h233_vQTN90_500_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF40_h233_vQTN90_500[, 1:6], window.size = 250000)
zm_vqgi_molike_MAF40_h233_vQTN90_full_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF40_h233_vQTN90_full[, 1:6], window.size = 250000)

#molike h263 MAF40 vQTN10
zm_vqgi_molike_MAF40_h263_vQTN10_500_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF40_h263_vQTN10_500[, 1:6], window.size = 250000)
zm_vqgi_molike_MAF40_h263_vQTN10_full_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF40_h263_vQTN10_full[, 1:6], window.size = 250000)

#molike h263 MAF40 vQTN 50
zm_vqgi_molike_MAF40_h263_vQTN50_500_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF40_h263_vQTN50_500[, 1:6], window.size = 250000)
zm_vqgi_molike_MAF40_h263_vQTN50_full_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF40_h263_vQTN50_full[, 1:6], window.size = 250000)

#molike h233 MAF10 vQTN 10
zm_vqgi_molike_MAF10_h233_vQTN10_500_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF10_h233_vQTN10_500[, 1:6], window.size = 250000)
zm_vqgi_molike_MAF10_h233_vQTN10_full_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF10_h233_vQTN10_full[, 1:6], window.size = 250000)

#molike h233 MAF10 vQTN 50
zm_vqgi_molike_MAF10_h233_vQTN50_500_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF10_h233_vQTN50_500[, 1:6], window.size = 250000)
zm_vqgi_molike_MAF10_h233_vQTN50_full_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF10_h233_vQTN50_full[, 1:6], window.size = 250000)

#molike h233 MAF10 vQTN90
zm_vqgi_molike_MAF10_h233_vQTN90_500 <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF10_h233_vQTN90_500[, 1:6], window.size = 250000)
zm_vqgi_molike_MAF10_h233_vQTN90_full <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF10_h233_vQTN90_full[, 1:6], window.size = 250000)

##Molike h263 MAF10
#molike h263 MAF10 vQTN 10
zm_vqgi_molike_MAF10_h263_vQTN10_500_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF10_h263_vQTN10_500[, 1:6], window.size = 250000)
zm_vqgi_molike_MAF10_h263_vQTN10_full_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF10_h263_vQTN10_full[, 1:6], window.size = 250000)

#molike h263 MAF10 vQTN 50
zm_vqgi_molike_MAF10_h263_vQTN50_500_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF10_h263_vQTN50_500[, 1:6], window.size = 250000)
zm_vqgi_molike_MAF10_h263_vQTN50_full_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF10_h263_vQTN50_full[, 1:6], window.size = 250000)

#molike h263 MAF10 vQTN90
zm_vqgi_molike_MAF10_h263_vQTN90_500_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF10_h263_vQTN90_500[, 1:6], window.size = 250000)
zm_vqgi_molike_MAF10_h263_vQTN90_full_250K <- vQTN_window_maker(var.QTN.genotypic.information = zm_var.QTN.genotypic.information_molike_MAF10_h263_vQTN90_full[, 1:6], window.size = 250000)


##Then, I will create shortcuts to result files
wd$at_bf_results <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_at/vGWAS_results/BF/Revisions"
wd$zm_bf_results <- paste0(wd$zm_vGWAS_results, "BFT/")

wd$at_DGLM_results <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_At/vGWAS_Results/DGLM/Revisions/"
wd$zm_DGLM_results <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_maize_filtlike_at/DGLM/"

#GAPIT results
wd$at_gapit_results <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_at/vGWAS_results/GAPIT/At/"
wd$zm_gapit_results <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_maize_filtlike_at/GAPIT/Zm/"

#
wd$at_DGLM_results_revisions <- paste0(wd$at_DGLM_results, "Revisions/")

##Here, I will also create the switch to 3_output folder 
wd$output <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/3_output/"

#check to see if your folder has the 100 result files
at_bf_dirs <- list.dirs(wd$at_bf_results)
at_bf_dirs <- at_bf_dirs[c(9, 10)]

zm_bf_dirs <- list.dirs(wd$zm_bf_results)
zm_bf_dirs <- zm_bf_dirs[c(3, 4, 6, 7, 9, 10, 12, 13, 15, 16)]

at_dglm_dirs <- list.dirs(wd$at_DGLM_results)

zm_dglm_dirs <- list.dirs(wd$zm_DGLM_results)


at_gapit_dirs <- list.dirs(wd$at_gapit_results)
zm_gapit_dirs <- list.dirs(wd$zm_gapit_results)

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
wd$output <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/3_output/"
setwd(wd$output)

################################ Detection Rates ############################################
#The first section will be Arabidopsis. For both Arabidopsis and maize, I will run Brown-Forsythe
#first then DGLM

#At
##BF
###Null
####500
Rate.calc(Result_dir = at_bf_dirs, rf = 57, simvQTN = NULL, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = TRUE, trait.name = "at_null_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 56, simvQTN = NULL, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = TRUE, trait.name = "at_null_full", test = "BF_", Additive = FALSE)

#Molike_MAF40
##500
Rate.calc(Result_dir = at_bf_dirs, rf = 37, simvQTN = at_vqgi_molike_MAF40_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF40_500", test = "BF_", Additive = F)

##FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 36, simvQTN = at_vqgi_molike_MAF40_full, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF40_full", test = "BF_", Additive = F)

###Molike MAF40 h233 vQTN 10
####500
Rate.calc(Result_dir = at_bf_dirs, rf = 39, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN10_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF40_h233_vQTN10_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 40, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN10_full, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF40_h233_vQTN10_FULL", test = "BF_", Additive = FALSE)

###Molike MAF40 h233 vQTN 50
####500
Rate.calc(Result_dir = at_bf_dirs, rf = 42, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN50_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF40_h233_vQTN50_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 43, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN50_full, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF40_h233_vQTN50_full", test = "BF_", Additive = FALSE)

###Molike MAF 40 vQTN 90
####500
Rate.calc(Result_dir = at_bf_dirs, rf = 45, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN90_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF40_h233_vQTN90_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 46, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN90_full, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF40_h233_vQTN90_full", test = "BF_", Additive = FALSE)

###Molike MAF 40 h263 vQTN 10
####500
Rate.calc(Result_dir = at_bf_dirs, rf = 48, simvQTN = at_vqgi_molike_wmQTL_MAF40_h263_vQTN10_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF40_h263_vQTN10_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 49, simvQTN = at_vqgi_molike_wmQTL_MAF40_h263_vQTN10_full, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF40_h263_vQTN10_full", test = "BF_", Additive = FALSE)

###Molike MAF40 h263 vQTN 50 
####500
Rate.calc(Result_dir = at_bf_dirs, rf = 51, simvQTN = at_vqgi_molike_wmQTL_MAF40_h263_vQTN50_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF40_h263_vQTN50_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 52, simvQTN = at_vqgi_molike_wmQTL_MAF40_h263_vQTN50_full, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF40_h263_vQTN50_full", test = "BF_", Additive = FALSE)

###Molike MAF10 h233 vQTN10
####500
####STILL NEED TO RUN
Rate.calc(Result_dir = at_bf_dirs, rf = 18, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN10_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF10_h233_vQTN10_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 19, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN10_full, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF10_h233_vQTN10_full", test = "BF_", Additive = FALSE)

###Molike MAF10 h233 vQTN50
####500
Rate.calc(Result_dir = at_bf_dirs, rf = 21, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN50_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF10_h233_vQTN50_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 22, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN50_full, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF10_h233_vQTN50_full", test = "BF_", Additive = FALSE)

###Molike MAF10 h233 vQTN90
####500
Rate.calc(Result_dir = at_bf_dirs, rf = 24, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN90_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF10_h233_vQTN90_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 25, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN90_full, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF10_h233_vQTN90_full", test = "BF_", Additive = FALSE)

####Molike MAF10 h263
###Molike MAF10 h263 vQTN10
####500
####STILL NEED TO RUN
Rate.calc(Result_dir = at_bf_dirs, rf = 27, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN10_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF10_h263_vQTN10_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 28, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN10_full, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF10_h263_vQTN10_full", test = "BF_", Additive = FALSE)

###Molike MAF10 h263 vQTN50
####500
Rate.calc(Result_dir = at_bf_dirs, rf = 30, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN50_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF10_h263_vQTN50_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 31, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN50_full, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF10_h263_vQTN50_full", test = "BF_", Additive = FALSE)

###Molike MAF10 h263 vQTN90
####500
Rate.calc(Result_dir = at_bf_dirs, rf = 33, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN90_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF10_h263_vQTN90_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 34, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN90_full, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_molike_MAF10_h263_vQTN90_full", test = "BF_", Additive = FALSE)

###Epistasis settings
#low effect size
####500
Rate.calc(Result_dir = at_bf_dirs, rf = 6, simvQTN = at_eqgi_epistasis_lowefs_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_epistasis_lowefs_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 7, simvQTN = at_eqgi_epistasis_lowefs_FULL, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_epistasis_lowefs_FULL", test = "BF_", Additive = FALSE)

#high effect size
####500
Rate.calc(Result_dir = at_bf_dirs, rf = 3, simvQTN = at_eqgi_epistasis_highefs_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_epistasis_highefs_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 4, simvQTN = at_eqgi_epistasis_highefs_FULL, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_epistasis_highefs_FULL", test = "BF_", Additive = FALSE)

#low h2 low effect size
####500
Rate.calc(Result_dir = at_bf_dirs, rf = 12, simvQTN = at_eqgi_epistasis_lowh2_lowefs_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_epistasis_lowh2_lowefs_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 13, simvQTN = at_eqgi_epistasis_lowh2_lowefs_FULL, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_epistasis_lowh2_lowefs_FULL", test = "BF_", Additive = FALSE)

#low h2 high effect size
####500
Rate.calc(Result_dir = at_bf_dirs, rf = 9, simvQTN = at_eqgi_epistasis_lowh2_highefs_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_epistasis_lowh2_highefs_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 10, simvQTN = at_eqgi_epistasis_lowh2_highefs_FULL, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_epistasis_lowh2_highefs_FULL", test = "BF_", Additive = FALSE)

#GxE settings
####500
Rate.calc(Result_dir = at_bf_dirs, rf = 15, simvQTN = at_gxeqgi_diff_500_100K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_gxe_diff_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_bf_dirs, rf = 16, simvQTN = at_gxeqgi_diff_FULL_100K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "at_gxe_diff_FULL", test = "BF_", Additive = FALSE)

#Maize
#GxE settings
####500
zm_gxeqgi_diff_500_250K <- zm_gxeqgi_diff_500_250K[, c(2, 1, 3:7)]

Rate.calc(Result_dir = zm_bf_dirs, rf = 3, simvQTN = zm_gxeqgi_diff_500_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_gxe_diff_500", test = "BF_", Additive = FALSE)

####FULL
zm_gxeqgi_diff_FULL_250K <- zm_gxeqgi_diff_FULL_250K[, c(2, 1, 3:7)]

Rate.calc(Result_dir = zm_bf_dirs, rf = 22, simvQTN = zm_gxeqgi_diff_FULL_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_gxe_diff_FULL", test = "BF_", Additive = FALSE)


#maize
##high effect size H2 80
####500 
zm_eqgi_epistasis_highh2_500_250K <- zm_eqgi_epistasis_highh2_500_250K[, c(2, 1, 3:7)]
Rate.calc(Result_dir = zm_bf_dirs, rf = 5, simvQTN = zm_eqgi_epistasis_highh2_500_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_epistasis_higheffectsize_500_250K", test = "BF_", Additive = FALSE)

####FULL
zm_eqgi_epistasis_highh2_FULL_250K <- zm_eqgi_epistasis_highh2_FULL_250K[, c(2, 1, 3:7)]

Rate.calc(Result_dir = zm_bf_dirs, rf = 24, simvQTN = zm_eqgi_epistasis_highh2_FULL_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_epistasis_highh280_FULL_250K", test = "BF_", Additive = FALSE)

#low h2 low effect size
####500
zm_eqgi_epistasis_lowh2_500_250K <- zm_eqgi_epistasis_lowh2_500_250K[, c(2, 1, 3:7)]

Rate.calc(Result_dir = zm_bf_dirs, rf = 6, simvQTN = zm_eqgi_epistasis_lowh2_500_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_epistasis_lowh2_500_250K", test = "BF_", Additive = FALSE)

####FULL
zm_eqgi_epistasis_lowh2_FULL_250K <- zm_eqgi_epistasis_lowh2_FULL_250K[, c(2, 1, 3:7)]

Rate.calc(Result_dir = zm_bf_dirs, rf = 25, simvQTN = zm_eqgi_epistasis_lowh2_FULL_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_epistasis_lowh2_FULL_250K", test = "BF_", Additive = FALSE)


###Molike MAF40 h233 vQTN 10
####500
Rate.calc(Result_dir = zm_bf_dirs, rf = 14, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN10_500_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_h233_vQTN10_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 32, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN10_full_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_h233_vQTN10_FULL", test = "BF_", Additive = FALSE)

###Molike MAF40 h233 vQTN 50
####500
Rate.calc(Result_dir = zm_bf_dirs, rf = 15, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN50_500_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_h233_vQTN50_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 33, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN50_full_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_h233_vQTN50_full", test = "BF_", Additive = FALSE)

###Molike MAF 40 vQTN 90
####500
Rate.calc(Result_dir = zm_bf_dirs, rf = 16, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN90_500_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_h233_vQTN90_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 34, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN90_full_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_h233_vQTN90_full", test = "BF_", Additive = FALSE)

###Molike MAF 40 h263 vQTN 10
####500
Rate.calc(Result_dir = zm_bf_dirs, rf = 17, simvQTN = zm_vqgi_molike_MAF40_h263_vQTN10_500_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_h263_vQTN10_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 37, simvQTN = zm_vqgi_molike_MAF40_h263_vQTN10_full_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_h263_vQTN10_full", test = "BF_", Additive = FALSE)

###Molike MAF40 h263 vQTN 50 
####500
Rate.calc(Result_dir = zm_bf_dirs, rf = 18, simvQTN = zm_vqgi_molike_MAF40_h263_vQTN50_500_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_h263_vQTN50_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 35, simvQTN = zm_vqgi_molike_MAF40_h263_vQTN50_full_250K, mQTN.wd = NULL, GAPIT = F,  DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_h263_vQTN50_full", test = "BF_", Additive = FALSE)

##Molike MAF40 h263 vQTN90
####500
Rate.calc(Result_dir = zm_bf_dirs, rf = 13, simvQTN = zm_vqgi_molike_MAF40_h263_vQTN50_500_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_h263_vQTN50_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 31, simvQTN = zm_vqgi_molike_MAF40_full, mQTN.wd = NULL, GAPIT = F,  DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_h263_vQTN90_full", test = "BF_", Additive = FALSE)


###Molike MAF10 h233 vQTN10
####500
####STILL NEED TO RUN
Rate.calc(Result_dir = zm_bf_dirs, rf = 7, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN10_500_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF10_h233_vQTN10_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 26, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN10_full_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF10_h233_vQTN10_full", test = "BF_", Additive = FALSE)

###Molike MAF10 h233 vQTN50
####500
Rate.calc(Result_dir = zm_bf_dirs, rf = 8, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN50_500_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF10_h233_vQTN50_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 27, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN50_full_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF10_h233_vQTN50_full", test = "BF_", Additive = FALSE)

###Molike MAF10 h233 vQTN90
####500
Rate.calc(Result_dir = zm_bf_dirs, rf = 9, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN90_500, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF10_h233_vQTN90_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 21, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN90_full, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF10_h233_vQTN90_full", test = "BF_", Additive = FALSE)

####Molike MAF10 h263
###Molike MAF10 h263 vQTN10
####500
####STILL NEED TO RUN
Rate.calc(Result_dir = zm_bf_dirs, rf = 10, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN10_500_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF10_h263_vQTN10_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 28, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN10_full, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF10_h263_vQTN10_full", test = "BF_", Additive = FALSE)

###Molike MAF10 h263 vQTN50
####500
Rate.calc(Result_dir = zm_bf_dirs, rf = 11, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN50_500_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF10_h263_vQTN50_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 29, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN50_full_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF10_h263_vQTN50_full", test = "BF_", Additive = FALSE)

###Molike MAF10 h263 vQTN90
####500
Rate.calc(Result_dir = zm_bf_dirs, rf = 12, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN90_500_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF10_h263_vQTN90_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 30, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN90_full_250K, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF10_h263_vQTN90_full", test = "BF_", Additive = FALSE)

###NULL
####500
Rate.calc(Result_dir = zm_bf_dirs, rf = 4, simvQTN = NULL, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = T, trait.name = "zm_NULL_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 23, simvQTN = NULL, mQTN.wd = NULL, GAPIT = F,  DGLM = FALSE, Null_results = T, trait.name = "zm_NULL_FULL", test = "BF_", Additive = FALSE)

##Molike MAF40
###500
Rate.calc(Result_dir = zm_bf_dirs, rf = 51, simvQTN = zm_vqgi_molike_MAF40_500_filtlike_at, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = F, trait.name = "zm_molike_MAF40_500", test = "BF_", Additive = FALSE)

###FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 52, simvQTN = zm_vqgi_molike_MAF40_full_filtlike_at, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = F, trait.name = "zm_molike_MAF40_FULL", test = "BF_", Additive = FALSE)


#####DGLM
####500
Rate.calc(Result_dir = at_dglm_dirs, rf = 15, simvQTN = at_gxeqgi_diff_500_100K, mQTN.wd = NULL, DGLM = TRUE, Null_results = FALSE, trait.name = "at_gxe_diff_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 16, simvQTN = at_gxeqgi_diff_FULL_100K, mQTN.wd = NULL, DGLM = TRUE, Null_results = FALSE, trait.name = "at_gxe_diff_FULL", test = "DGLM_", Additive = FALSE)

#maize
#With 1MB
Rate.calc(Result_dir = zm_dglm_dirs, rf = 15, simvQTN = zm_gxeqgi_diff_500_1MB, mQTN.wd = NULL, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_gxe_diff_500_1MB", test = "DGLM_", Additive = FALSE)

Rate.calc(Result_dir = zm_dglm_dirs, rf = 16, simvQTN = zm_gxeqgi_diff_FULL_1MB, mQTN.wd = NULL, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_gxe_diff_FULL_1MB", test = "DGLM_", Additive = FALSE)


##DGLM Results
###Null
####500
Rate.calc(Result_dir = at_dglm_dirs, rf = 60, simvQTN = NULL, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = TRUE, trait.name = "at_null_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 59, simvQTN = NULL, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = TRUE, trait.name = "at_null_full", test = "DGLM_", Additive = FALSE)

#Molike_MAF40
##Sub500
Rate.calc(Result_dir = at_dglm_dirs, rf = 56, simvQTN = at_vqgi_molike_MAF40_500, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF40_500", test = "DGLM_", Additive = F)

##FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 57, simvQTN = at_vqgi_molike_MAF40_full, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF40_full", test = "DGLM_", Additive = F)

#Epistasis
#low effect size
####
Rate.calc(Result_dir = at_dglm_dirs, rf = 6, simvQTN = at_eqgi_epistasis_lowefs_500, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_epistasis_lowefs_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 7, simvQTN = at_eqgi_epistasis_lowefs_FULL, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_epistasis_lowefs_FULL", test = "DGLM_", Additive = FALSE)

#high effect size
###500
Rate.calc(Result_dir = at_dglm_dirs, rf = 3, simvQTN = at_eqgi_epistasis_highefs_500, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_epistasis_highefs_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 4, simvQTN = at_eqgi_epistasis_highefs_FULL, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_epistasis_highefs_FULL", test = "DGLM_", Additive = FALSE)

#Molike MAF40 h263 vQTN10
##500
Rate.calc(Result_dir = at_dglm_dirs, rf = 45, simvQTN = at_vqgi_molike_wmQTL_MAF40_h263_vQTN10_500, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF40_h263_vQTN10_500", test = "DGLM_", Additive = FALSE)

##FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 46, simvQTN = at_vqgi_molike_wmQTL_MAF40_h263_vQTN10_full, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF40_h263_vQTN10_FULL", test = "DGLM_", Additive = FALSE)

#Molike MAF40 h263 vQTN50
##500
Rate.calc(Result_dir = at_dglm_dirs, rf = 48, simvQTN = at_vqgi_molike_wmQTL_MAF40_h263_vQTN50_500, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF40_h263_vQTN50_500", test = "DGLM_", Additive = FALSE)

##FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 49, simvQTN = at_vqgi_molike_wmQTL_MAF40_h263_vQTN50_full, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF40_h263_vQTN50_FULL", test = "DGLM_", Additive = FALSE)

############DGLM
#low h2 low effect size
####500
Rate.calc(Result_dir = at_dglm_dirs, rf = 12, simvQTN = at_eqgi_epistasis_lowh2_lowefs_500, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_epistasis_lowh2_lowefs_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 13, simvQTN = at_eqgi_epistasis_lowh2_lowefs_FULL, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_epistasis_lowh2_lowefs_FULL", test = "DGLM_", Additive = FALSE)

#low h2 high effect size
####500
Rate.calc(Result_dir = at_dglm_dirs, rf = 9, simvQTN = at_eqgi_epistasis_lowh2_highefs_500, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_epistasis_lowh2_highefs_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 10, simvQTN = at_eqgi_epistasis_lowh2_highefs_FULL, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_epistasis_lowh2_highefs_FULL", test = "DGLM_", Additive = FALSE)

###Molike MAF10 h233 vQTN10
####500
####STILL NEED TO RUN
Rate.calc(Result_dir = at_dglm_dirs, rf = 18, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN10_500, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "at_molike_MAF10_h233_vQTN10_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 19, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN10_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "at_molike_MAF10_h233_vQTN10_full", test = "DGLM_", Additive = FALSE)

###Molike MAF10 h233 vQTN50
####500
Rate.calc(Result_dir = at_dglm_dirs, rf = 21, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN50_500, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF10_h233_vQTN50_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 22, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN50_full, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF10_h233_vQTN50_full", test = "DGLM_", Additive = FALSE)

###Molike MAF10 h233 vQTN90
####500
Rate.calc(Result_dir = at_dglm_dirs, rf = 24, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN90_500, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF10_h233_vQTN90_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 25, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN90_full, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF10_h233_vQTN90_full", test = "DGLM_", Additive = FALSE)

####Molike MAF10 h263
###Molike MAF10 h263 vQTN10
####500
####STILL NEED TO RUN
Rate.calc(Result_dir = at_dglm_dirs, rf = 27, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN10_500, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "at_molike_MAF10_h263_vQTN10_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 28, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN10_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "at_molike_MAF10_h263_vQTN10_full", test = "DGLM_", Additive = FALSE)

###Molike MAF10 h263 vQTN50
####500
Rate.calc(Result_dir = at_dglm_dirs, rf = 30, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN50_500, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF10_h263_vQTN50_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 31, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN50_full, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF10_h263_vQTN50_full", test = "DGLM_", Additive = FALSE)

###Molike MAF10 h263 vQTN90
####500
Rate.calc(Result_dir = at_dglm_dirs, rf = 33, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN90_500, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF10_h263_vQTN90_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 34, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN90_full, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF10_h263_vQTN90_full", test = "DGLM_", Additive = FALSE)

#########MAF40
###Molike MAF40 h233 vQTN 10
####500
Rate.calc(Result_dir = at_dglm_dirs, rf = 36, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN10_500, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF40_h233_vQTN10_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 37, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN10_full, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "at_molike_MAF40_h233_vQTN10_FULL", test = "DGLM_", Additive = FALSE)

###Molike MAF40 h233 vQTN 50
####500
Rate.calc(Result_dir = at_dglm_dirs, rf = 39, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN50_500, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "at_molike_MAF40_h233_vQTN50_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 40, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN50_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "at_molike_MAF40_h233_vQTN50_full", test = "DGLM_", Additive = FALSE)

###Molike MAF 40 h233 vQTN 90
####500
Rate.calc(Result_dir = at_dglm_dirs, rf = 42, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN90_500, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "at_molike_MAF40_h233_vQTN90_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 43, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN90_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "at_molike_MAF40_h233_vQTN90_full", test = "DGLM_", Additive = FALSE)

###Molike MAF 40 h263 vQTN 10
####500
Rate.calc(Result_dir = at_dglm_dirs, rf = 45, simvQTN = at_vqgi_molike_wmQTL_MAF40_h263_vQTN10_500, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "at_molike_MAF40_h263_vQTN10_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 46, simvQTN = at_vqgi_molike_wmQTL_MAF40_h263_vQTN10_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "at_molike_MAF40_h263_vQTN10_full", test = "DGLM_", Additive = FALSE)

###Molike MAF40 h263 vQTN 50 
####500
Rate.calc(Result_dir = at_dglm_dirs, rf = 48, simvQTN = at_vqgi_molike_wmQTL_MAF40_h263_vQTN50_500, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "at_molike_MAF40_h263_vQTN50_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_dglm_dirs, rf = 49, simvQTN = at_vqgi_molike_wmQTL_MAF40_h263_vQTN50_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "at_molike_MAF40_h263_vQTN50_full", test = "DGLM_", Additive = FALSE)

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
Rate.calc(Result_dir = zm_bf_dirs, rf = 51, simvQTN = zm_vqgi_molike_MAF40_500_filtlike_at, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_500", test = "BF_", Additive = FALSE)

##FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 52, simvQTN = zm_vqgi_molike_MAF40_full_filtlike_at, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_full", test = "BF_", Additive = FALSE)

#Molike womQTL MAF40
##Sub500
Rate.calc(Result_dir = zm_bf_dirs, rf = 6, simvQTN = Zm_vqgi_molike_womQTL_MAF40_500, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_womQTL_MAF40_500", test = "BF_", Additive = FALSE)

##FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 7, simvQTN = zm_vqgi_molike_womQTL_MAF40_full, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_womQTL_MAF40_full", test = "BF_", Additive = FALSE)

#Molike filtered like Arabidopsis
##500
Rate.calc(Result_dir = zm_bf_dirs, rf = 48, simvQTN = zm_vqgi_molike_MAF40_500_filtlike_at, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_500_filtlike_at", test = "BF_", Additive = FALSE)

##FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 52, simvQTN = zm_vqgi_molike_MAF40_full_filtlike_at, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_MAF40_full_filtlike_at", test = "BF_", Additive = FALSE)

#Molike womQTL MAF40
##Sub500
Rate.calc(Result_dir = zm_bf_dirs, rf = 12, simvQTN = zm_vqgi_molike_womQTL_MAF40_500_filtlike_at, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_womQTL_MAF40_500_filtlike_at", test = "BF_", Additive = FALSE)

##FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 13, simvQTN = zm_vqgi_molike_womQTL_MAF40_full_filtlike_at, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_molike_womQTL_MAF40_full_filtlike_at", test = "BF_", Additive = FALSE)

###Epistasis settings
#low effect size
####500
Rate.calc(Result_dir = zm_bf_dirs, rf = 6, simvQTN = zm_eqgi_epistasis_lowefs_500_1MB, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_epistasis_lowefs_500", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 7, simvQTN = zm_eqgi_epistasis_lowefs_FULL_1MB, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_epistasis_lowefs_FULL_1MB", test = "BF_", Additive = FALSE)

#high effect size
###500
Rate.calc(Result_dir = zm_bf_dirs, rf = 3, simvQTN = zm_eqgi_epistasis_highefs_500_1MB, mQTN.wd = NULL, GAPIT = F, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_epistasis_highefs_500_1MB", test = "BF_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_bf_dirs, rf = 4, simvQTN = zm_eqgi_epistasis_highefs_FULL_1MB, mQTN.wd = NULL, DGLM = FALSE, Null_results = FALSE, trait.name = "zm_epistasis_highefs_FULL_1MB", test = "BF_", Additive = FALSE)


#DGLM Results
##Null
###500
Rate.calc(Result_dir = zm_dglm_dirs, rf = 54, simvQTN = NULL, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = TRUE, trait.name = "zm_null_500", test = "DGLM_", Additive = FALSE)

###FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 55, simvQTN = NULL, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = TRUE, trait.name = "zm_null_500", test = "DGLM_", Additive = FALSE)

#Molike MAF40
##500
Rate.calc(Result_dir = zm_dglm_dirs, rf = 36, simvQTN = zm_vqgi_molike_MAF40_500_filtlike_at, mQTN.wd = NULL, GAPIT = T, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_Molike_MAF40_500", test = "DGLM_", Additive = FALSE)

##FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 37, simvQTN = zm_vqgi_molike_MAF40_full_filtlike_at, mQTN.wd = NULL, GAPIT = T, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_molike_MAF40_full", test = "DGLM_", Additive = FALSE)

###Epistasis settings
#low effect size
####500
Rate.calc(Result_dir = zm_dglm_dirs, rf = 6, simvQTN = zm_eqgi_epistasis_lowefs_500_1MB, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_epistasis_lowefs_500_1MB", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 7, simvQTN = zm_eqgi_epistasis_lowefs_FULL_1MB, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_epistasis_lowefs_FULL_1MB", test = "DGLM_", Additive = FALSE)

#high effect sizes
###500
Rate.calc(Result_dir = zm_dglm_dirs, rf = 3, simvQTN = zm_eqgi_epistasis_highefs_500_1MB, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_epistasis_highefs_500_1MB", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 4, simvQTN = zm_eqgi_epistasis_highefs_FULL_1MB, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_epistasis_highefs_FULL_1MB", test = "DGLM_", Additive = FALSE)

#Epistasis different effect sizes and low heritability
#low h2 low effect size
####500
Rate.calc(Result_dir = zm_dglm_dirs, rf = 12, simvQTN = zm_eqgi_epistasis_lowh2_lowefs_500_1MB, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_epistasis_lowh2_lowefs_500_1MB", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 13, simvQTN = zm_eqgi_epistasis_lowh2_lowefs_FULL_1MB, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_epistasis_lowh2_lowefs_FULL_1MB", test = "DGLM_", Additive = FALSE)

#low h2 high effect size
####500
Rate.calc(Result_dir = zm_dglm_dirs, rf = 9, simvQTN = zm_eqgi_epistasis_lowh2_highefs_500_1MB, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_epistasis_lowh2_highefs_500_1MB", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 10, simvQTN = zm_eqgi_epistasis_lowh2_highefs_FULL_1MB, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_epistasis_lowh2_highefs_FULL_1MB", test = "DGLM_", Additive = FALSE)

###Molike MAF40 h233 vQTN 10
####500
Rate.calc(Result_dir = zm_dglm_dirs, rf = 39, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN10_500, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF40_h233_vQTN10_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 40, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN10_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF40_h233_vQTN10_FULL", test = "DGLM_", Additive = FALSE)

###Molike MAF40 h233 vQTN 50
####500
Rate.calc(Result_dir = zm_dglm_dirs, rf = 42, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN50_500, mQTN.wd = NULL, GAPIT = F, DGLM = TRUE, Null_results = FALSE, trait.name = "zm_molike_MAF40_h233_vQTN50_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 43, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN50_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF40_h233_vQTN50_full", test = "DGLM_", Additive = FALSE)

###Molike MAF 40 vQTN 90
####500
Rate.calc(Result_dir = zm_dglm_dirs, rf = 45, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN90_500, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF40_h233_vQTN90_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 46, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN90_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF40_h233_vQTN90_full", test = "DGLM_", Additive = FALSE)

###Molike MAF 40 h263 vQTN 10
####500
Rate.calc(Result_dir = zm_dglm_dirs, rf = 48, simvQTN = zm_vqgi_molike_MAF40_h263_vQTN10_500, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF40_h263_vQTN10_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 49, simvQTN = zm_vqgi_molike_MAF40_h263_vQTN10_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF40_h263_vQTN10_full", test = "DGLM_", Additive = FALSE)

###Molike MAF40 h263 vQTN 50 
####500
Rate.calc(Result_dir = zm_dglm_dirs, rf = 51, simvQTN = zm_vqgi_molike_MAF40_h263_vQTN50_500, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF40_h263_vQTN50_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 52, simvQTN = zm_vqgi_molike_MAF40_h263_vQTN50_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF40_h263_vQTN50_full", test = "DGLM_", Additive = FALSE)

###Molike MAF10 h233 vQTN10
####500
####STILL NEED TO RUN
Rate.calc(Result_dir = zm_dglm_dirs, rf = 18, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN10_500, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF10_h233_vQTN10_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 19, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN10_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF10_h233_vQTN10_full", test = "DGLM_", Additive = FALSE)

###Molike MAF10 h233 vQTN50
####500
Rate.calc(Result_dir = zm_dglm_dirs, rf = 21, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN50_500, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF10_h233_vQTN50_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 22, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN50_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF10_h233_vQTN50_full", test = "DGLM_", Additive = FALSE)

###Molike MAF10 h233 vQTN90
####500
Rate.calc(Result_dir = zm_dglm_dirs, rf = 24, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN90_500, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF10_h233_vQTN90_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 25, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN90_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF10_h233_vQTN90_full", test = "DGLM_", Additive = FALSE)

####Molike MAF10 h263
###Molike MAF10 h263 vQTN10
####500
####STILL NEED TO RUN
Rate.calc(Result_dir = zm_dglm_dirs, rf = 27, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN10_500, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF10_h263_vQTN10_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 28, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN10_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF10_h263_vQTN10_full", test = "DGLM_", Additive = FALSE)

###Molike MAF10 h263 vQTN50
####500
Rate.calc(Result_dir = zm_dglm_dirs, rf = 30, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN50_500, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF10_h263_vQTN50_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 31, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN50_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF10_h263_vQTN50_full", test = "DGLM_", Additive = FALSE)

###Molike MAF10 h263 vQTN90
####500
Rate.calc(Result_dir = zm_dglm_dirs, rf = 33, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN90_500, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF10_h263_vQTN90_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 34, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN90_full, mQTN.wd = NULL, GAPIT = F, DGLM = T, Null_results = FALSE, trait.name = "zm_molike_MAF10_h263_vQTN90_full", test = "DGLM_", Additive = FALSE)

###NULL
####500
Rate.calc(Result_dir = zm_dglm_dirs, rf = 54, simvQTN = NULL, mQTN.wd = NULL, DGLM = T, Null_results = T, trait.name = "zm_NULL_500", test = "DGLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_dglm_dirs, rf = 55, simvQTN = NULL, mQTN.wd = NULL, DGLM = T, Null_results = T, trait.name = "zm_NULL_FULL", test = "DGLM_", Additive = FALSE)


####This section is for GAPIT
#####Maize
######NULL
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 52, simvQTN = NULL, mQTN.wd = NULL, GAPIT = T, DGLM = FALSE, Null_results = T, trait.name = "zm_null_500", test = "GAPIT_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 53, simvQTN = NULL, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = T, trait.name = "zm_null_full", test = "GAPIT_", Additive = FALSE)

######Epistasis_higheffectsize_h280
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 3, simvQTN = zm_eqgi_epistasis_highefs_500_1MB, mQTN.wd = NULL, GAPIT = T, DGLM = FALSE, Null_results = F, trait.name = "zm_epistasis_h280_higheffectsize_500", test = "MLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 4, simvQTN = zm_eqgi_epistasis_highefs_FULL_1MB, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_epistasis_h280_higheffectsize_full", test = "MLM_", Additive = FALSE)

######Epistasis lowh2 higheffectsize
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 6, simvQTN = zm_eqgi_epistasis_lowh2_highefs_500_1MB, mQTN.wd = NULL, GAPIT = T, DGLM = FALSE, Null_results = F, trait.name = "zm_epistasis_lowh2_higheffectsizes_500", test = "GAPIT_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 7, simvQTN = zm_eqgi_epistasis_lowh2_highefs_FULL_1MB, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_epistasis_lowh2_higheffectsizes_full", test = "MLM_", Additive = FALSE)

#####Epistasis loweffectsize H2=0.80
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 9, simvQTN = zm_eqgi_epistasis_lowefs_500_1MB, mQTN.wd = NULL, GAPIT = T, DGLM = FALSE, Null_results = F, trait.name = "zm_epistasis_h280_loweffectsizes_500", test = "MLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 10, simvQTN = zm_eqgi_epistasis_lowefs_FULL_1MB, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_epistasis_h280_loweffectsizes_full", test = "MLM_", Additive = FALSE)

#######Epistasis loweffectsize H2=0.30
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 12, simvQTN = zm_eqgi_epistasis_lowh2_lowefs_500_1MB, mQTN.wd = NULL, GAPIT = T, DGLM = FALSE, Null_results = F, trait.name = "zm_epistasis_lowh2_loweffectsizes_500", test = "MLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 13, simvQTN = zm_eqgi_epistasis_lowh2_lowefs_FULL_1MB, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_epistasis_lowh2_loweffectsizes_full", test = "MLM_", Additive = FALSE)

###GxE differences
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 15, simvQTN = zm_gxeqgi_diff_500_1MB, mQTN.wd = NULL, GAPIT = T, DGLM = FALSE, Null_results = F, trait.name = "zm_gxe_diff_500", test = "MLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 16, simvQTN = zm_gxeqgi_diff_FULL_1MB, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_gxe_diff_full", test = "MLM_", Additive = FALSE)

###Molike MAF10 h233 vQTN10
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 18, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN10_500, mQTN.wd = NULL, GAPIT = T, DGLM = FALSE, Null_results = F, trait.name = "zm_molike_MAF10_h233_vQTN10_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 19, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN10_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF10_h233_vQTN10_full", test = "MLM_", Additive = F)

###Molike MAF10 h233 vQTN50
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 21, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN50_500, mQTN.wd = NULL, GAPIT = T, DGLM = FALSE, Null_results = F, trait.name = "zm_molike_MAF10_h233_vQTN50_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 22, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN50_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF10_h233_vQTN50_full", test = "MLM_", Additive = F)

###Molike MAF10 h233 vQTN90
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 24, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN90_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF10_h233_vQTN90_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 25, simvQTN = zm_vqgi_molike_MAF10_h233_vQTN90_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF10_h233_vQTN90_full", test = "MLM_", Additive = F)

###Molike MAF10 h263 vQTN10
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 27, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN10_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF10_h263_vQTN10_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 28, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN10_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF10_h263_vQTN10_full", test = "MLM_", Additive = F)

###Molike MAF10 h263 vQTN50
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 30, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN50_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF10_h263_vQTN50_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 31, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN50_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF10_h263_vQTN50_full", test = "MLM_", Additive = F)

###Molike MAF10 h263 vQTN90
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 32, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN90_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF10_h263_vQTN90_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 33, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN90_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF10_h263_vQTN90_full", test = "MLM_", Additive = F)

###Molike MAF40 
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 35, simvQTN = zm_vqgi_molike_MAF40_500_filtlike_at, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF40_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 36, simvQTN = zm_vqgi_molike_MAF40_full_filtlike_at, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF40_full", test = "MLM_", Additive = F)

###Molike MAF40 h233 vQTN10
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 38, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN10_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF40_h233_vQTN10_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 39, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN10_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF40_h233_vQTN10_full", test = "MLM_", Additive = F)

###Molike MAF40 h233 vQTN50
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 41, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN50_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF40_h233_vQTN50_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 42, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN50_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF40_h233_vQTN50_full", test = "MLM_", Additive = F)

###Molike MAF40 h233 vQTN90
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 44, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN90_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF40_h233_vQTN90_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 45, simvQTN = zm_vqgi_molike_MAF40_h233_vQTN90_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF40_h233_vQTN90_full", test = "MLM_", Additive = F)

###Molike MAF40 h263 vQTN10
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 47, simvQTN = zm_vqgi_molike_MAF40_h263_vQTN10_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF40_h263_vQTN10_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 48, simvQTN = zm_vqgi_molike_MAF40_h263_vQTN10_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF40_h263_vQTN10_full", test = "MLM_", Additive = F)

###Molike MAF40 h263 vQTN50
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 50, simvQTN = zm_vqgi_molike_MAF40_h263_vQTN50_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF40_h263_vQTN50_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 51, simvQTN = zm_vqgi_molike_MAF10_h263_vQTN50_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_molike_MAF40_h263_vQTN50_full", test = "MLM_", Additive = F)

####GxE 
####500
Rate.calc(Result_dir = zm_gapit_dirs, rf = 15, simvQTN = zm_gxeqgi_diff_500_1MB, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_gxe_diff_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = zm_gapit_dirs, rf = 16, simvQTN = zm_gxeqgi_diff_FULL_1MB, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "zm_gxe_diff_full", test = "MLM_", Additive = F)


#################Arabidopsis
######NULL
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 54, simvQTN = NULL, mQTN.wd = NULL, GAPIT = T, DGLM = FALSE, Null_results = T, trait.name = "at_null_500", test = "GAPIT_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 55, simvQTN = NULL, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = T, trait.name = "at_null_full", test = "GAPIT_", Additive = FALSE)

######Epistasis_higheffectsize_h280
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 3, simvQTN = at_eqgi_epistasis_highefs_500, mQTN.wd = NULL, GAPIT = T, DGLM = FALSE, Null_results = F, trait.name = "at_epistasis_h280_higheffectsize_500", test = "MLM_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 4, simvQTN = at_eqgi_epistasis_highefs_FULL, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_epistasis_h280_higheffectsize_full", test = "MLM_", Additive = FALSE)

######Epistasis lowh2 higheffectsize
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 6, simvQTN = at_eqgi_epistasis_lowh2_highefs_500, mQTN.wd = NULL, GAPIT = T, DGLM = FALSE, Null_results = F, trait.name = "at_epistasis_lowh2_higheffectsizes_500", test = "GAPIT_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 7, simvQTN = at_eqgi_epistasis_lowh2_highefs_FULL, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_epistasis_lowh2_higheffectsizes_full", test = "GAPIT_", Additive = FALSE)

#####Epistasis loweffectsize H2=0.80
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 9, simvQTN = at_eqgi_epistasis_lowefs_500, mQTN.wd = NULL, GAPIT = T, DGLM = FALSE, Null_results = F, trait.name = "at_epistasis_h280_loweffectsizes_500", test = "GAPIT_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 10, simvQTN = at_eqgi_epistasis_lowefs_FULL, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_epistasis_h280_loweffectsizes_full", test = "MLM_", Additive = FALSE)

#######Epistasis loweffectsize H2=0.30
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 12, simvQTN = at_eqgi_epistasis_lowh2_lowefs_500, mQTN.wd = NULL, GAPIT = T, DGLM = FALSE, Null_results = F, trait.name = "at_epistasis_lowh2_loweffectsizes_500", test = "GAPIT_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 13, simvQTN = at_eqgi_epistasis_lowh2_lowefs_FULL, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_epistasis_lowh2_loweffectsizes_full", test = "GAPIT_", Additive = FALSE)

###GxE differences
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 15, simvQTN = at_gxeqgi_diff_500_100K, mQTN.wd = NULL, GAPIT = T, DGLM = FALSE, Null_results = F, trait.name = "at_gxe_diff_500", test = "GAPIT_", Additive = FALSE)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 16, simvQTN = at_gxeqgi_diff_FULL_100K, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_gxe_diff_full", test = "MLM_", Additive = FALSE)

###Molike MAF10 h233 vQTN10
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 18, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN10_500, mQTN.wd = NULL, GAPIT = T, DGLM = FALSE, Null_results = F, trait.name = "at_molike_MAF10_h233_vQTN10_500", test = "GAPIT_", Additive = F)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 19, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN10_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF10_h233_vQTN10_full", test = "GAPIT_", Additive = F)

###Molike MAF10 h233 vQTN50
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 21, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN50_500, mQTN.wd = NULL, GAPIT = T, DGLM = FALSE, Null_results = F, trait.name = "at_molike_MAF10_h233_vQTN50_500", test = "GAPIT_", Additive = F)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 22, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN50_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF10_h233_vQTN50_full", test = "GAPIT_", Additive = F)

###Molike MAF10 h233 vQTN90
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 24, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN90_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF10_h233_vQTN90_500", test = "GAPIT_", Additive = F)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 25, simvQTN = at_vqgi_molike_wmQTL_MAF10_h233_vQTN90_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF10_h233_vQTN90_full", test = "GAPIT_", Additive = F)

###Molike MAF10 h263 vQTN10
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 27, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN10_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF10_h263_vQTN10_500", test = "GAPIT_", Additive = F)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 28, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN10_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF10_h263_vQTN10_full", test = "GAPIT_", Additive = F)

###Molike MAF10 h263 vQTN50
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 30, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN50_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF10_h263_vQTN50_500", test = "GAPIT_", Additive = F)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 31, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN50_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF10_h263_vQTN50_full", test = "GAPIT_", Additive = F)

###Molike MAF10 h263 vQTN90
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 33, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN90_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF10_h263_vQTN90_500", test = "GAPIT_", Additive = F)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 34, simvQTN = at_vqgi_molike_wmQTL_MAF10_h263_vQTN90_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF10_h263_vQTN90_full", test = "GAPIT_", Additive = F)

###Molike MAF40 
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 36, simvQTN = , mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF40_h233_vQTN10_500", test = "GAPIT_", Additive = F)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 37, simvQTN = , mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF40_h233_vQTN10_full", test = "GAPIT_", Additive = F)

###Molike MAF40 h233 vQTN10
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 39, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN10_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF40_h233_vQTN10_500", test = "GAPIT_", Additive = F)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 40, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN10_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF40_h233_vQTN10_full", test = "GAPIT_", Additive = F)

###Molike MAF40 h233 vQTN50
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 42, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN50_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF40_h233_vQTN50_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 43, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN50_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF40_h233_vQTN50_full", test = "MLM_", Additive = F)

###Molike MAF40 h233 vQTN90
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 45, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN90_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF40_h233_vQTN90_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 46, simvQTN = at_vqgi_molike_wmQTL_MAF40_h233_vQTN90_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF40_h233_vQTN90_full", test = "MLM_", Additive = F)

###Molike MAF40 h263 vQTN10
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 48, simvQTN = at_vqgi_molike_wmQTL_MAF40_h263_vQTN10_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF40_h263_vQTN10_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 49, simvQTN = at_vqgi_molike_wmQTL_MAF40_h263_vQTN10_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF40_h263_vQTN10_full", test = "MLM_", Additive = F)

###Molike MAF40 h263 vQTN50
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 51, simvQTN = at_vqgi_molike_wmQTL_MAF40_h263_vQTN50_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF40_h263_vQTN50_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 52, simvQTN = at_vqgi_molike_wmQTL_MAF40_h263_vQTN50_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF40_h263_vQTN50_full", test = "MLM_", Additive = F)

###Molike MAF40 h263 vQTN90
####500
Rate.calc(Result_dir = at_gapit_dirs, rf = 36, simvQTN = at_vqgi_molike_MAF40_500, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF40_h263_vQTN90_500", test = "MLM_", Additive = F)

####FULL
Rate.calc(Result_dir = at_gapit_dirs, rf = 37, simvQTN = at_vqgi_molike_MAF40_full, mQTN.wd = NULL, GAPIT = T, DGLM = F, Null_results = F, trait.name = "at_molike_MAF40_h263_vQTN90_full", test = "MLM_", Additive = F)


################################End########################################
