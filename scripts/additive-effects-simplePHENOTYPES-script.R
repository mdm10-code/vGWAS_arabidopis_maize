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

wd$vGWAS_chapter_one_zm <- paste0(wd$out, "vGWAS_chapter_one_zm/") #Maize
wd$zm_sp_results <- paste0(wd$vGWAS_chapter_one_zm, "zm_sp_results/")

#Here, I will load my maize genotypic file as Zm_geno_hapmap
zm_geno_hapmap_full <- read.delim2(paste0(wd$manual, "ames_2532_mind09_maf01.hmp.txt"), sep = "")

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
at_geno_numeric_full <- read.delim2(paste0(wd$sp_at_null, "out_geno_numeric.txt"), sep = "\t", head = TRUE)

#For my 500 susbet, I will use R sample function on my At1000Geno
#!!!!!WARNING!!!!!! Only run this once!
at_sub500 <- sample(6:ncol(at_geno_numeric_full), 500, replace = FALSE)
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

###These next two settings are for molike settings. For the molike settings
##The MAF of the additive QTN is 0.40. Also, this will be the same QTN used for
##additive effects for ALL molike simulations
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
  constraints =  list(maf_above = 0.395, maf_below = 0.405, hets = NULL),
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
  constraints =  list(maf_above = 0.395, maf_below = 0.405, hets = NULL),
  model = "A",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$at_sp_results)
)

###########################GxE setting#######################

#Arabidopsis GxE setting
##500
sp_at_GxE_500_nocor <- create_phenotypes(
  geno_obj = at_geno_numeric_500,
  ntraits = 2,
  pleio_a = 1,
  trait_spec_a_QTN_num = c(1, 1),
  h2 = c(0.7, 0.7),
  add_effect = list(trait1 = c(0.2, 0.5), trait2 = c(0.9, 0.5)),
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
  add_effect = list(trait1 = c(0.2, 0.5), trait2 = c(0.9, 0.5)),
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
sp_at_GxE_500$mean <- (sp_at_GxE_500$Trait_1_H2_0.7 + sp_at_GxE_500$Trait_2_H2_0.7) / 2
sp_at_GxE_FULL$mean <- (sp_at_GxE_FULL$Trait_1_H2_0.7 + sp_at_GxE_FULL$Trait_2_H2_0.7) / 2

#Now, I will use tidyverse to get a wide format
sp_at_GxE_500_mod <- sp_at_GxE_500[, c(1, 4, 5)] %>%
  pivot_wider(names_from = Rep, values_from = mean)

sp_at_GxE_FULL_mod <- sp_at_GxE_FULL[, c(1, 4, 5)] %>%
  pivot_wider(names_from = Rep, values_from = mean)

####################################################################
##This next section is for the Epistasis settings

#This is Arabidopsis 500 Diversity Panel with high effect sizes and high heritability
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

##FULL diversity panel with high effect size and high heritability
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


#This is Arabidopsis 500 Diversity Panel with high effect sizes and low heritability
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

#FULL diversity panel for 
#full diversity panel


##FULL diversity panel with high effect size
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

########################Maize Simulations###########################
##Note, the following simulations were done for a revision. The
##paths might not be the same as above
#And finally, we load our 
load("vGWAS_chaptI_RII.RData")

###Here, This is my way of organizing all my input/output directories based on Rich Pauloo's way of 
#organizing his data science projects. This is a preference. All of this is similar to setwd, but 
#eliminated the need of going back and forth between directories

#First, we will create an empty list for directories
wd <- list()

#This the directory for the genotypic dataset in maize 
wd$RII <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/"

#I will create the directories for my project
##Data file
wd$exgeno <- paste0(wd$RII, "0_data/External/")
wd$out <- paste0(wd$RII, "2_pipeline/out/")

#In These next directories, I will create folders for all the simulations by species. Within each 
#species specific folder, I will create paths to the simplePHENOTYPES (SP) simulated results. First, 
#I will need to create the folders
dir.create(paste0(wd$out, "vGWAS_chapter_one_at/")) #Arabidopsis
dir.create(paste0(wd$out, "vGWAS_chapter_one_zm/")) #Maize

#These next two will be for the simplePHENOTYPES folders
dir.create(paste0(wd$out, "vGWAS_chapter_one_at/at_sp_results/")) 
dir.create(paste0(wd$out, "vGWAS_chapter_one_zm/zm_sp_results/"))

#Now, I will create working directory shortcuts
##Arabidopsis
wd$vGWAS_chapter_one_at <- paste0(wd$out, "vGWAS_chapter_one_at/") #Arabidopsis
wd$at_sp_results <- paste0(wd$vGWAS_chapter_one_at, "at_sp_results/") 

##Maize
wd$vGWAS_chapter_one_zm <- paste0(wd$out, "vGWAS_chapter_one_zm/") #Arabidopsis
wd$zm_sp_results <- paste0(wd$vGWAS_chapter_one_zm, "zm_sp_results/") 

##
###This is just for filtering purposes
create_phenotypes(
  geno_path = paste0(wd$exgeno, "Genotypic/"),
  prefix = "ames_2815_mac5_num",
  add_QTN_num = 0,
  h2 = 0,
  add_effect = 0.00,
  ntraits = 1,
  rep = 1,
  output_format = "wide",
  out_geno = "BED", 
  output_dir = "zm.for.bed",
  model = "A",
  vary_QTN = F,
  to_r = F,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

###Now I convert my the filtered bed file from PLINK into a numeric genotypic object
###I will also do the first maize setting as null


##I will create wd
wd$ingeno <- paste0(wd$RII, "0_data/Internal/")

sp.zm.null.full <- create_phenotypes(
  geno_path = paste0(wd$ingeno),
  prefix = "ames2815_r220_100_25.bed",
  add_QTN_num = 0,
  h2 = 0,
  add_effect = 0.00,
  ntraits = 1,
  rep = 100,
  seed = 9399,
  output_format = "wide",
  out_geno = "numeric", 
  output_dir = "sp_zm_null_full",
  model = "A",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

#To do the null subset 500, I will need to upload the out_geno
#from above
wd$sp_zm_null_full <- paste0(wd$zm_sp_results, "sp_zm_null_full/")
zm.geno.numeric.full <- read.delim("~/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/zm_sp_results/sp_zm_null_full/out_geno_numeric.txt")

#As with Arabidopsis, I will randomly random accessions using R's sample function. Just like
#Arabidopsis, ONLY RUN THIS ONCE!
zm.500 <- sample(6:ncol(zm.geno.numeric.full), 500, replace = FALSE)
zm.geno.numeric.500 <- zm.geno.numeric.full[, c(1:5, zm.500)]

zm.geno.numeric.500 <- as.data.frame(zm.geno.numeric.500)

#this is for the null trait subset 500
sp.zm.null.500.filt <- create_phenotypes(
  geno_obj = zm.geno.numeric.500,
  add_QTN_num = 0,
  h2 = 0,
  add_effect = 0.00,
  ntraits = 1,
  rep = 100,
  seed = 839350,
  output_format = "wide",
  output_dir = "sp_zm_null_500",
  model = "A",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

###This is to simulate a Molybdenum-like setting with a single mQTL
sp.zm.molike.500 <- create_phenotypes(
  geno_obj = zm.geno.numeric.500,
  add_QTN_num = 1,
  h2 = 1.00,
  add_effect = 0.25,
  ntraits = 1,
  rep = 100,
  seed = 136068,
  output_format = "wide",
  output_dir = "sp_zm_molike_wmQTL_500",
  constraints =  list(maf_above = 0.395, maf_below = 0.405, hets = NULL),
  model = "A",
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

sp.zm.molike.FULL <- create_phenotypes(
  geno_obj = zm.geno.numeric.full,
  add_QTN_num = 1,
  h2 = 1.00,
  add_effect = 0.25,
  ntraits = 1,
  rep = 100,
  seed = 413385,
  output_format = "wide",
  output_dir = "sp_zm_molike_FULL",
  constraints = list(maf_above = 0.395, maf_below = 0.405, hets=NULL),
  model = "A",
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

#Now, we will simulate epistatic settings with a high heritability
##Not working 2-21-2022
epi_effect_vect <- c(0.75, 0.75, 0.75)
zm.sp.epi.higheffectsize.500 <- create_phenotypes(
  geno_obj = zm.geno.numeric.500,
  epi_QTN_num = 3,
  epi_effect = list(c(0.75, 0.75, 0.75)),
  h2 = 0.80,
  sim_method = "custom",
  seed = 395223,
  ntraits = 1,
  rep = 100,
  output_format = "wide",
  output_dir = "sp_zm_epi_higheffectsize_500",
  constraints =  list(maf_above = 0.095, maf_below = 0.105, hets = NULL),
  model = "E",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

#Full
zm.sp.epi.higheffectsize.full.h280 <- create_phenotypes(
  geno_obj = zm.geno.numeric.full,
  epi_QTN_num = 3,
  epi_effect = list(c(0.75, 0.75, 0.75)),
  h2 = 0.80,
  sim_method = "custom",
  ntraits = 1,
  rep = 100,
  seed = 324620,
  output_format = "wide",
  output_dir = "sp_zm_epi_higheffectsize_full",
  constraints =  list(maf_above = 0.095, maf_below = 0.105, hets = NULL),
  model = "E",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

#Now, we will simulate epistatic settings with a low heritability
zm.sp.epi.higheffectsize.h230.500 <- create_phenotypes(
  geno_obj = zm.geno.numeric.500,
  epi_QTN_num = 3,
  epi_effect = list(c(0.75, 0.75, 0.75)),
  h2 = 0.30,
  sim_method = "custom",
  ntraits = 1,
  rep = 100,
  seed = 395223,
  output_format = "wide",
  output_dir = "sp_zm_epi_higheffectsize_h230_500",
  constraints =  list(maf_above = 0.095, maf_below = 0.105, hets = NULL),
  model = "E",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

#Full
zm.sp.epi.higheffectsize.h230.full <- create_phenotypes(
  geno_obj = zm.geno.numeric.full,
  epi_QTN_num = 3,
  epi_effect = list(c(0.75, 0.75, 0.75)),
  h2 = 0.30,
  sim_method = "custom",
  ntraits = 1,
  rep = 100,
  seed = 324620,
  output_format = "wide",
  output_dir = "sp_zm_epi_higheffectsize_h230_full",
  constraints =  list(maf_above = 0.095, maf_below = 0.105, hets = NULL),
  model = "E",
  vary_QTN = FALSE,
  to_r = TRUE,
  SNP_effect = "Add",
  home_dir = paste0(wd$zm_sp_results)
)

#GxE setting
###500
sp.zm.GxE.500 <- create_phenotypes(
  geno_obj = zm.geno.numeric.500,
  ntraits = 2,
  pleio_a = 1,
  trait_spec_a_QTN_num = c(1, 1),
  h2 = c(0.7, 0.7),
  add_effect = list(trait1 = c(0.2, 0.5), trait2 = c(0.9, 0.5)),
  rep = 100,
  seed = 513018,
  output_dir = "sp_zm_GxE_500",
  output_format = "wide",
  architecture = "partially",
  sim_method = "custom",
  constraints =  list(maf_above = 0.295, maf_below = 0.305, hets = NULL),
  to_r = TRUE,
  model = "A",
  home_dir = paste0(wd$zm_sp_results)
)

##FULL
sp.zm.GxE.FULL <- create_phenotypes(
  geno_obj = zm.geno.numeric.full,
  ntraits = 2,
  pleio_a = 1,
  trait_spec_a_QTN_num = c(1, 1),
  h2 = c(0.7, 0.7),
  add_effect = list(trait1 = c(0.2, 0.5), trait2 = c(0.9, 0.5)),
  rep = 100,
  seed = 7713,
  output_dir = "sp_zm_GxE_FULL",
  output_format = "wide",
  architecture = "partially",
  sim_method = "custom",
  constraints =  list(maf_above = 0.295, maf_below = 0.305, hets = NULL),
  to_r = TRUE,
  model = "A",
  home_dir = paste0(wd$zm_sp_results)
)

###In this next section, I will need to modify maize GxE objects to get the differnece
#500
sp.zm.GxE.500$diff <- sp.zm.GxE.500$Trait_1_H2_0.7 - sp.zm.GxE.500$Trait_2_H2_0.7

#FULL
sp.zm.GxE.FULL$diff <- sp.zm.GxE.FULL$Trait_1_H2_0.7 - sp.zm.GxE.FULL$Trait_2_H2_0.7

sp.zm.GxE.500.mod <- sp.zm.GxE.500[, c(1, 4, 5)] %>%
  pivot_wider(names_from = Rep, values_from = diff)

sp.zm.GxE.FULL.mod <- sp.zm.GxE.FULL[, c(1, 4, 5)] %>%
  pivot_wider(names_from = Rep, values_from = diff)

##################################################################
save.image("vGWAS_chaptI_RII.RData")
save.image("vGWAS_chaptI_RII_copy.RData")
save.image("vGWAS_chaptI_RII_copy2.RData")

#############################################################
#This last section is to output a numeric dataset
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
