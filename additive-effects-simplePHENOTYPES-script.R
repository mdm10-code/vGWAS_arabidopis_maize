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

wd$vGWAS_chapter_one_zm <- paste0(wd$out, "vGWAS_chapter_one_zm/") #Maize
wd$zm_sp_results <- paste0(wd$vGWAS_chapter_one_zm, "zm_sp_results/")

#Here, I will load my maize genotypic file as Zm_geno_hapmap
zm_geno_hapmap_full <- read.delim2(paste0(wd$manual, "ames_2532_mind09_maf01.hmp.txt"), sep="")

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
  seed=748577,
  output_format = "wide",
  output_dir = "sp_at_molike_wmQTL_full",
  constraints =  list(maf_above = 0.395, maf_below = 0.405, hets=NULL),
  model = "A",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$at_sp_results)
)

##############################This next section is for maize########################################

sp_zm_null_full <- create_phenotypes(
  geno_obj = zm_geno_hapmap_full,
  add_QTN_num = 0,
  h2 = 0,
  add_effect = 0.00,
  ntraits = 1,
  rep = 100,
  seed = 276692,
  output_format = "wide",
  out_geno = "numeric", 
  output_dir = "sp_zm_null_full",
  model = "A",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

#I created a numeric genotypic file that I will use for my later simulations. To upload it, I need
#to create the path to this numeric genotypic file
wd$sp_zm_null <- paste0(wd$zm_sp_results, "sp_zm_null_full/")
zm_geno_numeric_full <- read.table(paste0(wd$sp_zm_null, "zm_geno_hapmap_full_numeric.txt"), sep="", head=T)

#As with Arabidopsis, I will randomly random accessions using R's sample function. Just like
#Arabidopsis, ONLY RUN THIS ONCE!
zm_sub500 <- sample(6:ncol(zm_geno_numeric_full), 500, replace=FALSE)
zm_geno_numeric_500 <- zm_geno_numeric_full[, c(1:5, zm_sub500)]

#this is for the null trait
sp_zm_null_500 <- create_phenotypes(
  geno_obj = zm_geno_numeric_500,
  add_QTN_num = 0,
  h2 = 0,
  add_effect = 0.00,
  ntraits = 1,
  rep = 100,
  seed=276692,
  output_format = "wide",
  output_dir = "sp_zm_null_500",
  model = "A",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

###This is to simulate a Molybdenum-like setting with a single mQTL
sp_zm_molike_wmQTL_500 <- create_phenotypes(
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
  seed=461803,
  output_format = "wide",
  output_dir = "sp_zm_molike_wmQTL_FULL",
  constraints =  list(maf_above = 0.395, maf_below = 0.405, hets=NULL),
  model = "A",
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

save.image("vGWAS_chapter_one.RData")
