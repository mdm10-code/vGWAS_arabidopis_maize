#####Script 1: simplePHENOTYPES#####
#The purpose of this script is to simulate the additive genetic effects of our settings using 
#simplePHENOTYPES

#Remove everything from the environment and clear the console
rm(list())
cat("\014")

#If not done so, install the simplePHENOTYPES R Package
setRepositories(ind = 1:2)
devtools::install_github("samuelbfernandes/simplePHENOTYPES", build_vignettes = TRUE)

library(simplePHENOTYPES)

#And finally, we load our 
load("vGWAS_chapter_one.RData")

###Here, This is my way of organizing all my input/output directories based on Rich Pauloo's way of 
#organizing his data science projects. This is a preference. All of this is similar to setwd, but 
#eliminated the need of going back and forth between directories

#First, we will create an empty list for directories
wd <- list()

#This the directory for the genotypic dataset in maize 
wd$manual <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/0_data/Manual/"

#This directory is for the results (output)
wd$out <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/"

#In These next directories, I will create folders for all the simulations by species. Within each 
#species specific folder, I will create paths to the simplePHENOTYPES (SP) simulated results. First, 
#I will need to create the folders
dir.create(paste0(wd$out, "vGWAS_chapter_one_at/")) #Arabidopsis
dir.create(paste0(wd$out, "vGWAS_chapter_one_zm/")) #Maize

#These next two will be for the simplePHENOTYPES folders
dir.create(paste0(wd$out, "vGWAS_chapter_one_at/at_sp_results/")) 
dir.create(paste0(wd$out, "vGWAS_chapter_one_zm/zm_sp_results/"))

#Now, I will create working directory shortcuts
wd$vGWAS_chapter_one_at <- paste0(wd$out, "vGWAS_chapter_one_at/") #Arabidopsis
wd$at_sp_results <- paste0(wd$vGWAS_chapter_one_at, "at_sp_results/") 

###For this next section, We will simulate NULL settings and Molybdenum-like settings with a single 
#mQTL. In this first section, I will create the settings for A. thaliana

#For my NULL trait, I will go ahead and use my pruned GBS dataset and output it as a numeric file 
#for downstream analyzes I am just going to do this once
at_sp_null_full <- create_phenotypes(
  geno_path = paste0(wd$manual),
  prefix = "At1000G_LDpruned.bed",
  add_QTN_num = 0,
  h2 = 0,
  add_effect = 0,
  ntraits = 1,
  rep = 100,
  seed = 200,
  output_format = "wide",
  out_geno = "numeric",
  output_dir = "sp_at_null_full",
  model = "A",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$at_sp_results)
)

#Here, I will upload the numeric genotypic file created with the previous simulation
#setting. I will need to create that shortcut
wd$sp_at_null_full <- paste0(wd$at_sp_results, "sp_at_null_full/")
at_geno_numeric_full <- read.delim2(paste0(wd$sp_at_null, "out_geno_numeric.txt"), sep="\t", head=TRUE)

#For my 500 susbet, I will use R sample function on my At1000Geno
#!!!!!WARNING!!!!!! Only run this once!
at_sub500 <- sample(6:ncol(at_geno_numeric_full), 500, replace=FALSE)
at_geno_numeric_500 <- at_geno_numeric_full[, c(1:5, at_sub500)]

#Arabidopsis Null Setting 500 sample size diversity panel
at_sp_null_500 <- create_phenotypes(
  geno_obj = at_geno_numeric_500,
  add_QTN_num = 0,
  h2 = 0,
  add_effect = 0,
  ntraits = 1,
  rep = 100,
  seed = 200,
  output_format = "wide",
  output_dir = "sp_at_null_500",
  model = "A",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$at_sp_results)
)

###Molybdenum-like setting with MAF = 0.40, h2 = 0.63, and vQTN effect size = 0.9
#Sub500
at_sp_molike_wMQTL_500 <- create_phenotypes(
  geno_obj = at_geno_numeric_500,
  add_QTN_num = 1,
  h2 = 1.0,
  add_effect = 0.25,
  ntraits = 1,
  rep = 100,
  seed = 65903,
  output_format = "wide",
  output_dir = "sp_at_molike_wmQTL_500",
  constraints =  list(maf_above = 0.395, maf_below = 0.405, hets=NULL),
  model = "A",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$at_sp_results)
)

#Full dataset
at_sp_molike_wmQTL_full <- create_phenotypes(
  geno_obj = at_geno_numeric_full,
  add_QTN_num = 1,
  h2 = 1.0,
  add_effect = 0.25,
  ntraits = 1,
  rep = 100,
  seed = 748577,
  output_format = "wide",
  output_dir = "sp_at_molike_wmQTL_full",
  constraints =  list(maf_above = 0.395, maf_below = 0.405, hets=NULL),
  model = "A",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$at_sp_results)
)

#For these next settings, we will explore how vGWAS can detect Epistatic and GxE Interactions
#For the epistatic setting, we will create custome effect sizes for to be used for both
#both Arabidopsis and maize. These settings are loosely based on what was done in Forsberg and
#Carlborg 2017

#In these next few lines, I will be creating 36 different settings to see how smaller effect vQTLS
#effects detection rates

#Arabidopsis GxE setting
##500
sp_at_GxE_500_nocor <- create_phenotypes(
  geno_obj = at_geno_numeric_500,
  ntraits = 2,
  pleio_a = 1,
  trait_spec_a_QTN_num = c(1, 1),
  h2 = c(0.7, 0.7),
  add_effect = list( trait1 = c(0.2, 0.5), trait2 = c(0.9, 0.5)),
  rep = 100,
  output_dir = "sp_at_GxE_500_nocor",
  output_format = "wide",
  architecture = "partially",
  sim_method = "custom",
  constraints =  list(maf_above = 0.295, maf_below = 0.305, hets = NULL),
  to_r = TRUE,
  model = "A",
  home_dir = paste0(wd$at_sp_results)
)

##FULL
sp_at_GxE_FULL_nocor <- create_phenotypes(
  geno_obj = at_geno_numeric_full,
  ntraits = 2,
  pleio_a = 1,
  trait_spec_a_QTN_num = c(1, 1),
  h2 = c(0.7, 0.7),
  add_effect = list( trait1 = c(0.2, 0.5), trait2 = c(0.9, 0.5)),
  cor = cor_matrix,
  rep = 100,
  output_dir = "sp_at_GxE_FULL_nocor",
  output_format = "wide",
  architecture = "partially",
  sim_method = "custom",
  constraints =  list(maf_above = 0.295, maf_below = 0.305, hets = NULL),
  to_r = TRUE,
  model = "A",
  home_dir = paste0(wd$at_sp_results)
)

###Now, we need to wrangle the data. First, I will add a mean column to each row
#Arabidopsis
sp_at_500_GxE_diff <- (sp_at_GxE_500_mod_t1[, 2:101] - sp_at_GxE_500_mod_t2[, 2:101])
sp_at_500_GxE_diff <- data.frame(sp_at_GxE_500_mod_t1[, 1], sp_at_500_GxE_diff)

#FULL
sp_at_FULL_GxE_diff <- (sp_at_GxE_FULL_mod_t1[, 2:101] - sp_at_GxE_FULL_mod_t2[, 2:101])
sp_at_FULL_GxE_diff <- data.frame(sp_at_GxE_FULL_mod_t1[, 1], sp_at_FULL_GxE_diff)

#Now, I will use tidyverse to get a wide format
sp_at_GxE_500_mod <- sp_at_500_diff[, c(1, 4, 5)] %>%
  pivot_wider(names_from = Rep, values_from = diff)

sp_at_GxE_FULL_mod <- sp_at_FULL_GxE_diff[, c(1, 4, 5)] %>%
  pivot_wider(names_from = Rep, values_from = diff)

##This next section is for epistasis
####Heritability = 0.80
at_sp_epi_higheffectsize_500 <- create_phenotypes(
  geno_obj = at_geno_numeric_500,
  epi_QTN_num = 3,
  epi_effect = c(0.75, 0.75, 0.75),
  h2 = 0.80,
  sim_method = "custom",
  ntraits = 1,
  rep = 100,
  seed = 199822,
  output_format = "wide",
  output_dir = "sp_at_epi_higheffectsize_500",
  constraints =  list(maf_above = 0.095, maf_below = 0.105, hets = NULL),
  model = "E",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$at_sp_results)
)

##FULL diversity panel for High Heritability
at_sp_epi_higheffectsize_full <- create_phenotypes(
  geno_obj = at_geno_numeric_full,
  epi_QTN_num = 3,
  epi_effect = c(0.75, 0.75, 0.75),
  h2 = 0.80,
  sim_method = "custom",
  ntraits = 1,
  rep = 100,
  seed = 103633,
  output_format = "wide",
  output_dir = "sp_at_epi_higheffectsize_full",
  constraints =  list(maf_above = 0.095, maf_below = 0.105, hets = NULL),
  model = "E",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$at_sp_results)
)

##Epistasis Broad-Sense Heritability = 0.3
###Arabidopsis 500 Diversity Panel
at_sp_epi_higheffectsize_h230_500 <- create_phenotypes(
  geno_obj = at_geno_numeric_500,
  epi_QTN_num = 3,
  epi_effect = c(0.75, 0.75, 0.75),
  h2 = 0.30,
  sim_method = "custom",
  ntraits = 1,
  rep = 100,
  seed = 199822,
  output_format = "wide",
  output_dir = "sp_at_epi_higheffectsize_h230_500_II",
  constraints =  list(maf_above = 0.095, maf_below = 0.105, hets = NULL),
  model = "E",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$at_sp_results)
)

##FULL diversity panel with low Broad-Sense Heritability
at_sp_epi_higheffectsize_h230_full <- create_phenotypes(
  geno_obj = at_geno_numeric_full,
  epi_QTN_num = 3,
  epi_effect = c(0.75, 0.75, 0.75),
  h2 = 0.30,
  sim_method = "custom",
  ntraits = 1,
  rep = 100,
  seed = 103633,
  output_format = "wide",
  output_dir = "sp_at_epi_higheffectsize_h230_full",
  constraints =  list(maf_above = 0.095, maf_below = 0.105, hets = NULL),
  model = "E",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$at_sp_results)
)


####################################################################
#This next section is a request to rerun the maize settings with a genotypic file
#filtered the same way as Arabidopsis. To obtain this, I will go ahead
sp_zm_null_full_filt_like_at <- create_phenotypes(
  geno_path = paste0(wd$manual),
  prefix = "zm_ames_pruned.bed",
  add_QTN_num = 0,
  h2 = 0,
  add_effect = 0.00,
  ntraits = 1,
  rep = 100,
  seed = 276692,
  output_format = "wide",
  out_geno = "numeric", 
  output_dir = "sp_zm_null_full_filt_like_at",
  model = "A",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

#I created a numeric genotypic file that I will use for my later simulations. To upload it, I need
#to create the path to this numeric genotypic file
wd$sp_zm_null_full_filt_like_at <- paste0(wd$zm_sp_results, "sp_zm_null_full_filt_like_at/")
zm_geno_numeric_full_filt_like_at <- read.table(paste0(wd$sp_zm_null_full_filt_like_at, "out_geno_numeric.txt"), sep="", head=T)

#As with Arabidopsis, I will randomly random accessions using R's sample function. Just like
#Arabidopsis, ONLY RUN THIS ONCE!
#zm_sub500 <- sample(6:ncol(zm_geno_numeric_full), 500, replace=FALSE)
zm_geno_numeric_500_filt_like_at <- zm_geno_numeric_full_filt_like_at[, c(1:5, zm_sub500)]

#this is for the null trait
sp_zm_null_500_filt_like_at <- create_phenotypes(
  geno_obj = zm_geno_numeric_500_filt_like_at,
  add_QTN_num = 0,
  h2 = 0,
  add_effect = 0.00,
  ntraits = 1,
  rep = 100,
  seed = 276692,
  output_format = "wide",
  output_dir = "sp_zm_null_500_filt_like_at",
  model = "A",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

###This is to simulate a Molybdenum-like setting with a single mQTL
sp_zm_molike_wmQTL_500_filt <- create_phenotypes(
  geno_obj = zm_geno_numeric_500,
  add_QTN_num = 1,
  h2 = 1.00,
  add_effect = 0.25,
  ntraits = 1,
  rep = 100,
  seed=678630,
  output_format = "wide",
  output_dir = "sp_zm_molike_wmQTL_500",
  constraints =  list(maf_above = 0.395, maf_below = 0.405, hets=NULL),
  model = "A",
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

sp_zm_molike_wmQTL_FULL <- create_phenotypes(
  geno_obj = zm_geno_numeric_full,
  add_QTN_num = 1,
  h2 = 1.00,
  add_effect = 0.25,
  ntraits = 1,
  rep = 100,
  seed = 461803,
  output_format = "wide",
  output_dir = "sp_zm_molike_wmQTL_FULL",
  constraints =  list(maf_above = 0.395, maf_below = 0.405, hets=NULL),
  model = "A",
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)


#######################New Epistasis settings##########
#Epistasis with low effect sizes
#This is maize 500 Diversity Panel


#Now, we will simulate epistatic settings with a high heritability
zm_sp_epi_higheffectsize_500_filt_like_at <- create_phenotypes(
  geno_obj = zm_geno_numeric_500_filt_like_at,
  epi_QTN_num = 3,
  epi_effect = c(0.75, 0.75, 0.75),
  h2 = 0.80,
  sim_method = "custom",
  ntraits = 1,
  rep = 100,
  seed = 36421,
  output_format = "wide",
  output_dir = "sp_zm_epi_higheffectsize_500_filt_like_at",
  constraints =  list(maf_above = 0.095, maf_below = 0.105, hets = NULL),
  model = "E",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

#Full
zm_sp_epi_higheffectsize_full_filt_like_at <- create_phenotypes(
  geno_obj = zm_geno_numeric_full_filt_like_at,
  epi_QTN_num = 3,
  epi_effect = c(0.75, 0.75, 0.75),
  h2 = 0.80,
  sim_method = "custom",
  ntraits = 1,
  rep = 100,
  seed = 811334,
  output_format = "wide",
  output_dir = "sp_zm_epi_higheffectsize_full_filt_like_at",
  constraints =  list(maf_above = 0.095, maf_below = 0.105, hets = NULL),
  model = "E",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

#Now, we will simulate epistatic settings with a low heritability
zm_sp_epi_higheffectsize_h230_500_filt_like_at <- create_phenotypes(
  geno_obj = zm_geno_numeric_500_filt_like_at,
  epi_QTN_num = 3,
  epi_effect = c(0.75, 0.75, 0.75),
  h2 = 0.30,
  sim_method = "custom",
  ntraits = 1,
  rep = 100,
  seed = 36421,
  output_format = "wide",
  output_dir = "sp_zm_epi_higheffectsize_h230_500_filt_like_at",
  constraints =  list(maf_above = 0.095, maf_below = 0.105, hets = NULL),
  model = "E",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

#Full
zm_sp_epi_higheffectsize_h230_full_filt_like_at <- create_phenotypes(
  geno_obj = zm_geno_numeric_full_filt_like_at,
  epi_QTN_num = 3,
  epi_effect = c(0.75, 0.75, 0.75),
  h2 = 0.30,
  sim_method = "custom",
  ntraits = 1,
  rep = 100,
  seed = 811334,
  output_format = "wide",
  output_dir = "sp_zm_epi_higheffectsize_h230_full_filt_like_at",
  constraints =  list(maf_above = 0.095, maf_below = 0.105, hets = NULL),
  model = "E",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

###This is to simulate a Molybdenum-like setting with a single mQTL, MAF = 0.40, h2=0.63, vQTN = 0.9
sp_zm_molike_wmQTL_500_filtlike_at <- create_phenotypes(
  geno_obj = zm_geno_numeric_500_filt_like_at,
  add_QTN_num = 1,
  h2 = 1.00,
  add_effect = 0.25,
  ntraits = 1,
  rep = 100,
  seed = 678630,
  output_format = "wide",
  output_dir = "sp_zm_molike_wmQTL_500_filtlike_at",
  constraints =  list(maf_above = 0.395, maf_below = 0.405, hets=NULL),
  model = "A",
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

#FULL
sp_zm_molike_wmQTL_FULL_filtlike_at <- create_phenotypes(
  geno_obj = zm_geno_numeric_full_filt_like_at,
  add_QTN_num = 1,
  h2 = 1.00,
  add_effect = 0.25,
  ntraits = 1,
  rep = 100,
  seed = 461803,
  output_format = "wide",
  output_dir = "sp_zm_molike_wmQTL_FULL_filtlike_at",
  constraints =  list(maf_above = 0.395, maf_below = 0.405, hets=NULL),
  model = "A",
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

#GxE setting
###500
sp_zm_GxE_500_nocor <- create_phenotypes(
  geno_obj = zm_geno_numeric_500_filt_like_at,
  ntraits = 2,
  pleio_a = 1,
  trait_spec_a_QTN_num = c(1, 1),
  h2 = c(0.7, 0.7),
  add_effect = list( trait1 = c(0.2, 0.5), trait2 = c(0.9, 0.5)),
  rep = 100,
  output_dir = "sp_zm_GxE_500_nocor",
  output_format = "wide",
  architecture = "partially",
  sim_method = "custom",
  constraints =  list(maf_above = 0.295, maf_below = 0.305, hets = NULL),
  to_r = TRUE,
  model = "A",
  home_dir = paste0(wd$zm_sp_results)
)

##FULL
sp_zm_GxE_FULL_nocor <- create_phenotypes(
  geno_obj = zm_geno_numeric_full_filt_like_at,
  ntraits = 2,
  pleio_a = 1,
  trait_spec_a_QTN_num = c(1, 1),
  h2 = c(0.7, 0.7),
  add_effect = list( trait1 = c(0.2, 0.5), trait2 = c(0.9, 0.5)),
  rep = 100,
  output_dir = "sp_zm_GxE_FULL_nocor",
  output_format = "wide",
  architecture = "partially",
  sim_method = "custom",
  constraints =  list(maf_above = 0.295, maf_below = 0.305, hets = NULL),
  to_r = TRUE,
  model = "A",
  home_dir = paste0(wd$zm_sp_results)
)

###In this next section, I will need to modify maize GxE objects to get the differnece
sp_zm_500_GxE_diff <- (sp_zm_GxE_500_mod_t1[, 2:101] - sp_zm_GxE_500_mod_t2[, 2:101])
sp_zm_500_GxE_diff <- data.frame(sp_zm_GxE_500_mod_t1[, 1], sp_zm_500_GxE_diff)

#FULL
sp_zm_FULL_GxE_diff <- (sp_zm_GxE_FULL_mod_t1[, 2:101] - sp_zm_GxE_FULL_mod_t2[, 2:101])
sp_zm_FULL_GxE_diff <- data.frame(sp_zm_GxE_FULL_mod_t1[, 1], sp_zm_FULL_GxE_diff)

sp_zm_GxE_500_mod <- sp_zm_500_GxE_diff[, c(1, 4, 5)] %>%
  pivot_wider(names_from = Rep, values_from = diff)

sp_zm_GxE_FULL_mod <- sp_zm_FULL_GxE_diff[, c(1, 4, 5)] %>%
  pivot_wider(names_from = Rep, values_from = diff)

############################################################################
#This last section is to output a numeric dataset for the Applied vGWAS in Maize
create_phenotypes(
  geno_path = "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/0_data/Manual/",
  prefix = "SNP55K_maize282_AGP3_20190419_Filtered.hmp",
  add_QTN_num = 0,
  h2 = 0,
  add_effect = 0.00,
  ntraits = 1,
  rep = 100,
  output_format = "wide",
  out_geno = "numeric", 
  output_dir = "applied",
  model = "A",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

save.image("vGWAS_chapter_one.RData")
save.image("vGWAS_chapter_one_copy.RData")
