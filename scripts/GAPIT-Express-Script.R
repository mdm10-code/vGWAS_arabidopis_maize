#########################vGWAS-GAPIT SCRIPT#################
load("~/vGWAS_chaptI_RII.RData")

setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/RII/2_pipeline/out")

###Install and read in GAPIT v4
source("http://zzlab.net/GAPIT/GAPIT.library.R") 
source("http://zzlab.net/GAPIT/gapit_functions.txt")

###This is for your GAPIT results, change the "" everytime for naming
#Your GAPIT results
#For Arabdipsosis
##500 diversity panel
dir.create(paste0(wd$vGWAS_chapter_one_at, "at_GAPIT_500/"))
wd$at_GAPIT_500 <- paste0(wd$vGWAS_chapter_one_at, "at_GAPIT_500/")

##FULL
dir.create(paste0(wd$vGWAS_chapter_one_at, "at_GAPIT_full/"))
wd$at_GAPIT_full <- paste0(wd$vGWAS_chapter_one_at, "at_GAPIT_full/")

#For maize
##500 diversity panel
dir.create(paste0(wd$vGWAS_chapter_one_zm, "zm_GAPIT_500/"))
wd$zm_GAPIT_500 <- paste0(wd$vGWAS_chapter_one_zm, "zm_GAPIT_500/")

##FULL
dir.create(paste0(wd$vGWAS_chapter_one_zm, "zm_GAPIT_full/"))
wd$zm_GAPIT_full <- paste0(wd$vGWAS_chapter_one_zm, "zm_GAPIT_full/")

#Although there are a lot of ways to obtain PCAs, I will use GAPIT.
##Because I am only needing PCAs, I will just run GAPIT with one
##replication

####For separate GM and GD
####A way to read in files is using separate arguments for map and genotypic objects
###For downstream analyses, I will need to transpose just genotypic objects
#Arabidopsis
at_just_geno_500 <- t(at_just_geno_500)
colnames(at_just_geno_500) <- at_map$snp

at_just_geno_full <- t(at_just_geno_full)
colnames(at_just_geno_full) <- at_map$snp

#maize
zm.just.geno.500 <- t(zm.just.geno.500)
colnames(zm.just.geno.500) <- zm.map$snp

zm.just.geno.full <- t(zm.just.geno.full)
colnames(zm.just.geno.full) <- zm.map$snp

###A thaliana
#To manage the files better, I will use the wds 
setwd(wd$at_GAPIT_500)

#500 Diversity Panel
myGAPIT_at_500 <- GAPIT(
  Y = data.frame(at_sp_null_500[, c(1, 2)]), 
  GD = data.frame(at_sp_null_500[, 1], at_just_geno_500),
  GM = at_map[, c(1, 3, 4)],
  PCA.total = 10, 
  SNP.effect = "Add"
)

#Full diversity panel
setwd(wd$at_GAPIT_full)

myGAPIT_at_full <- GAPIT(
  Y = data.frame(at_sp_null_full[, c(1, 2)]), 
  GD = data.frame(at_sp_null_full[, 1], at_just_geno_full),
  GM = at_map[, c(1, 3, 4)],
  PCA.total = 10, 
  SNP.effect = "Add"
)

###We can actually save our PCAs from the myGAPIT lists
at_covar_500 <- myGAPIT_at_500[["PCA"]]
at_covar_FULL <- myGAPIT_at_full[["PCA"]]

#maize
##500 diversity panel
setwd(wd$zm_GAPIT_500)
myGAPIT.zm.500 <- GAPIT(
  Y = data.frame(sp.zm.null.500.filt[, c(1, 2)]), 
  GD = data.frame(sp.zm.null.500.filt[, 1], zm.just.geno.500),
  GM = zm.map[, c(1, 3, 4)],
  PCA.total = 10, 
  SNP.effect = "Add"
) 

##full diversity panel
setwd(wd$zm_GAPIT_full)
myGAPIT.zm.full <- GAPIT(
  Y = data.frame(sp.zm.null.full[, c(1, 2)]), 
  GD = data.frame(sp.zm.null.full[, 1], zm.just.geno.full),
  GM = zm.map[, c(1, 3, 4)],
  PCA.total = 10, 
  SNP.effect = "Add"
)

zm_covar_500 <- myGAPIT.zm.500[["PCA"]]
zm_covar_500 <- zm_covar_500[, c(1:5)]

zm_covar_full <- myGAPIT.zm.full[["PCA"]]
zm_covar_full <- zm_covar_full[, c(1:5)]

#Now I extract the kinship matrices
myKI.zm.500 <- myGAPIT.zm.500[["KI"]]
myKI.zm.full <- myGAPIT.zm.full[["KI"]]

############
#This next section is to run all simulations using GAPIT
#The following section is for Arabidopsis

##NULL
###500
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/NULL/500/")
GAPIT(
  Y = data.frame(at_sp_null_500), 
  GD = data.frame(at_sp_null_500[, 1], at_just_geno_500),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_500,
  group.from = 500,
  group.to = 500,
  SNP.effect = "Add",
  model = "MLM"
)

###FULL
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/NULL/FULL/")
GAPIT(
  Y = data.frame(at_sp_null_full), 
  GD = data.frame(at_sp_null_full[, 1], at_just_geno_full),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_full,
  group.from = 1087,
  group.to = 1087,
  SNP.effect = "Add",
  model = "MLM"
)

##Epistasis lowh2 
###500
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/Epistasis_lowh2/500/")
GAPIT(
  Y = data.frame(at_sp_epi_lowh2_500), 
  GD = data.frame(at_sp_epi_lowh2_500[, 1], at_just_geno_500),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_500,
  group.from = 500,
  group.to = 500,
  SNP.effect = "Add",
  model = "MLM"
)

###FULL
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/Epistasis_lowh2/FULL/")
GAPIT(
  Y = data.frame(at_sp_epi_lowh2_full), 
  GD = data.frame(at_sp_epi_lowh2_full[, 1], at_just_geno_full),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_full,
  group.from = 1087,
  group.to = 1087,
  SNP.effect = "Add",
  model = "MLM"
)

##Molike
###500
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/Molike/500/")
GAPIT(
  Y = data.frame(at_vpheno_molike_MAF40_500), 
  GD = data.frame(at_vpheno_molike_MAF40_500[, 1], at_just_geno_500),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_500,
  group.from = 500,
  group.to = 500,
  SNP.effect = "Add",
  model = "MLM"
)

###FULL
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/Molike/FULL/")
GAPIT(
  Y = data.frame(at_vpheno_molike_MAF40_full), 
  GD = data.frame(at_vpheno_molike_MAF40_full[, 1], at_just_geno_full),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_full,
  group.from = 1087,
  group.to = 1087,
  SNP.effect = "Add",
  model = "MLM"
)



##Epistasis highh2
###500
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/Epistasis_highh2/500/")
GAPIT(
  Y = data.frame(at_sp_epi_highh2_500), 
  GD = data.frame(at_sp_epi_highh2_500[, 1], at_just_geno_500),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_500,
  group.from = 500,
  group.to = 500,
  SNP.effect = "Add",
  model = "MLM"
)

###FULL
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/Epistasis_highh2/FULL/")
GAPIT(
  Y = data.frame(at_sp_epi_highh2_full), 
  GD = data.frame(at_sp_epi_highh2_full[, 1], at_just_geno_full),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_full,
  group.from = 1087,
  group.to = 1087,
  SNP.effect = "Add",
  model = "MLM"
)

##Molike MAF40 h233 vQTN10
###500
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/Molike_MAF40_h233_vQTN10/500/")
GAPIT(
  Y = data.frame(at_vpheno_molike_MAF40_h233_vQTN10_500), 
  GD = data.frame(at_vpheno_molike_MAF40_h233_vQTN10_500[, 1], at_just_geno_500),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_500,
  group.from = 500,
  group.to = 500,
  SNP.effect = "Add",
  model = "MLM"
)

###FULL
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/Molike_MAF40_h233_vQTN10/FULL/")
GAPIT(
  Y = data.frame(at_vpheno_molike_MAF40_h233_vQTN10_full), 
  GD = data.frame(at_vpheno_molike_MAF40_h233_vQTN10_full[, 1], at_just_geno_full),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_full,
  group.from = 1087,
  group.to = 1087,
  SNP.effect = "Add",
  model = "MLM"
)

##Molike MAF40 h233 vQTN50
###500
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/Molike_MAF40_h233_vQTN50/500/")
GAPIT(
  Y = data.frame(at_vpheno_molike_MAF40_h233_vQTN50_500), 
  GD = data.frame(at_vpheno_molike_MAF40_h233_vQTN50_500[, 1], at_just_geno_500),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_500,
  group.from = 500,
  group.to = 500,
  SNP.effect = "Add",
  model = "MLM"
)

###FULL
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/Molike_MAF40_h233_vQTN50/FULL/")
GAPIT(
  Y = data.frame(at_vpheno_molike_MAF40_h233_vQTN50_full), 
  GD = data.frame(at_vpheno_molike_MAF40_h233_vQTN50_full[, 1], at_just_geno_full),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_full,
  group.from = 1087,
  group.to = 1087,
  SNP.effect = "Add",
  model = "MLM"
)

#
##Molike MAF40 h233 vQTN90
###500
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/Molike_MAF40_h233_vQTN90/500/")
GAPIT(
  Y = data.frame(at_vpheno_molike_MAF40_h233_vQTN90_500), 
  GD = data.frame(at_vpheno_molike_MAF40_h233_vQTN90_500[, 1], at_just_geno_500),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_500,
  KI = myK_500,
  group.from = 500,
  group.to = 500,
  SNP.effect = "Add",
  model = "MLM"
)

###FULL
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/Molike_MAF40_h233_vQTN90/FULL/")
GAPIT(
  Y = data.frame(at_vpheno_molike_MAF40_h233_vQTN90_full), 
  GD = data.frame(at_vpheno_molike_MAF40_h233_vQTN90_full[, 1], at_just_geno_full),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_full,
  KI = myK_at_FULL,
  group.from = 1087,
  group.to = 1087,
  SNP.effect = "Add",
  model = "MLM"
)

#Molike MAF40 vQTN10 h263
###500
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/Molike_MAF40_h263_vQTN10/500/")
GAPIT(
  Y = data.frame(at_vpheno_molike_MAF40_h263_vQTN10_500), 
  GD = data.frame(at_vpheno_molike_MAF40_h263_vQTN10_500[, 1], at_just_geno_500),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_500,
  KI = myK_500,
  group.from = 500,
  group.to = 500,
  SNP.effect = "Add",
  model = "MLM"
)

###FULL
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/Molike_MAF40_h263_vQTN10/FULL/")
GAPIT(
  Y = data.frame(at_vpheno_molike_MAF40_h263_vQTN10_full), 
  GD = data.frame(at_vpheno_molike_MAF40_h263_vQTN10_full[, 1], at_just_geno_full),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_full,
  KI = myK_at_FULL,
  group.from = 1087,
  group.to = 1087,
  SNP.effect = "Add",
  model = "MLM"
)

#Molike MAF40 vQTN50 h263
###500
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/Molike_MAF40_h263_vQTN50/500/")
GAPIT(
  Y = data.frame(at_vpheno_molike_MAF40_h263_vQTN50_500), 
  GD = data.frame(at_vpheno_molike_MAF40_h263_vQTN50_500[, 1], at_just_geno_500),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_500,
  KI = myK_500,
  group.from = 500,
  group.to = 500,
  SNP.effect = "Add",
  model = "MLM"
)

###FULL
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/Molike_MAF40_h263_vQTN50/FULL/")
GAPIT(
  Y = data.frame(at_vpheno_molike_MAF40_h263_vQTN50_full), 
  GD = data.frame(at_vpheno_molike_MAF40_h263_vQTN50_full[, 1], at_just_geno_full),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_full,
  KI = myK_at_FULL,
  group.from = 1087,
  group.to = 1087,
  SNP.effect = "Add",
  model = "MLM"
)

#Before I proceed, I am going to load the kinship matrices
#from previous runs
myK <- read.csv("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/NULL/500/GAPIT.Kin.VanRaden.csv", header = F, sep = ",")

#GxE Settings
###500
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/GxE/500/")
GAPIT(
  Y = data.frame(sp_at_GxE_500_mod), 
  GD = data.frame(sp_at_GxE_500_mod[, 1], at_just_geno_500),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_500,
  KI = myK_500,
  group.from = 500,
  group.to = 500,
  SNP.effect = "Add",
  model = "MLM"
)

myK_at_FULL <- read.csv("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/NULL/FULL/GAPIT.Kin.VanRaden.csv", header = F, sep = ",")
myK_at_FULL$V1 <- gsub("^", "X", myK_at_FULL$V1)

###FULL
setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/2_pipeline/out/GAPIT_Results/Arabidopsis/GxE/FULL/")
GAPIT(
  Y = data.frame(sp_at_GxE_FULL_mod), 
  GD = data.frame(sp_at_GxE_FULL_mod[, 1], at_just_geno_full),
  GM = at_map[, c(1, 3, 4)],
  CV = at_covar_full,
  K = myK_at_FULL,
  group.from = 1087,
  group.to = 1087,
  SNP.effect = "Add",
  model = "MLM"
)

############################# GAPIT Expressway #############################################
#For this section, to run GAPIT over a course of days, I am going to create a for loop that
#iterates through a setting list
GAPIT.expressway <- function(wd = NULL, trait_list = NULL, GD = NULL, GM = NULL, CV = NULL, K = NULL) {
  for (h in 1:length(trait_list)) {
    print(paste0("this is trait", h))
    new_folder <- dir.create(paste0(wd, names(trait_list)[[h]]))
    setwd(paste0(wd, names(trait_list)[[h]]))
    myY <- trait_list[[h]]
    GAPIT(
      Y = data.frame(myY), 
      GD = GD,
      GM = GM,
      CV = CV,
      K = K,
      group.from = dim(myY)[1],
      group.to = dim(myY)[1],
      SNP.effect = "Add",
      model = "MLM"
    )
  }
}

#Here, I will create my first setting list
zm.pheno.FULL.listI <- list(sp.zm.null.full,
                            sp.zm.GxE.FULL.mod,
                            zm.vpheno.molike.MAF10.h233.vQTN10.FULL,
                            zm.vpheno.molike.MAF10.h233.vQTN50.FULL,
                            zm.vpheno.molike.MAF10.h233.vQTN90.FULL,
                            zm.vpheno.molike.MAF10.h263.vQTN10.FULL)

#For each element in my list, I am going to give them the setting name. Make sure to double-check
#to see each element is identitical to each individual setting
names(zm.pheno.FULL.listI)[[1]] <- "sp.zm.null.full"
names(zm.pheno.FULL.listI)[[2]] <- "sp.zm.GxE.FULL.mod"
names(zm.pheno.FULL.listI)[[3]] <- "zm.vpheno.molike.MAF10.h233.vQTN10.FULL"
names(zm.pheno.FULL.listI)[[4]] <- "zm.vpheno.molike.MAF10.h233.vQTN50.FULL"
names(zm.pheno.FULL.listI)[[5]] <- "zm.vpheno.molike.MAF10.h233.vQTN90.FULL"
names(zm.pheno.FULL.listI)[[6]] <- "zm.vpheno.molike.MAF10.h263.vQTN10.FULL"

#Change directory to maize settings for full association panels
wd$zm_gapit_FULL <- "~/projects/dissertation/RII/2_pipeline/out/GAPIT/zm/FULL/"

#It also appears that the formating of the taxa is different between all
#phenotypic objects. To make sure that GAPIT doesn't stop running
#I will modify the taxa as seen in sp.zm.null.full
sp.zm.GxE.FULL.mod <- as.data.frame(sp.zm.GxE.FULL.mod)

rownames(sp.zm.GxE.FULL.mod) <- sp.zm.null.full$`<Trait>`
rownames(zm.vpheno.molike.MAF10.h233.vQTN10.FULL) <- sp.zm.null.full$`<Trait>`
rownames(zm.vpheno.molike.MAF10.h233.vQTN50.FULL) <- sp.zm.null.full$`<Trait>`
rownames(zm.vpheno.molike.MAF10.h233.vQTN90.FULL) <- sp.zm.null.full$`<Trait>`
rownames(zm.vpheno.molike.MAF10.h263.vQTN10.FULL) <- sp.zm.null.full$`<Trait>`

#
sp.zm.GxE.FULL.mod$`<Trait>` <- sp.zm.null.full$`<Trait>`
zm.vpheno.molike.MAF10.h233.vQTN10.FULL$taxa <- sp.zm.null.full$`<Trait>`
zm.vpheno.molike.MAF10.h233.vQTN50.FULL$taxa <- sp.zm.null.full$`<Trait>`
zm.vpheno.molike.MAF10.h233.vQTN90.FULL$taxa <- sp.zm.null.full$`<Trait>`
zm.vpheno.molike.MAF10.h263.vQTN10.FULL$taxa <- sp.zm.null.full$`<Trait>`
zm.vpheno.molike.MAF40.h263.vQTN10.FULL$taxa <- sp.zm.null.full$`<Trait>`


#I am going to extract the first three PCAs for maize
zm_covar_full <- myGAPIT.zm.full[["PCA"]]
zm_covar_full <- zm_covar_full[, c(1:4)]

#Trait list I
GAPIT.expressway(wd = wd$zm_gapit_FULL, trait_list = zm.pheno.FULL.listI, GD = data.frame(sp.zm.null.full[, 1], zm.just.geno.full), GM = zm.map[, c(1, 3, 4)], CV = zm_covar_full, K = myKI.zm.full)

#################### Memory Issues ###I had to create another setting list for which replications did 
#not run
zm.pheno.FULL.listI <- list(
  sp.zm.GxE.FULL.mod[, c(1, 81:101)],
  zm.vpheno.molike.MAF10.h233.vQTN50.FULL)


names(zm.pheno.FULL.listI)[[1]] <- "sp.zm.GxE.FULL.mod"
names(zm.pheno.FULL.listI)[[2]] <- "zm.vpheno.molike.MAF10.h233.vQTN50.FULL"

GAPIT.expressway(wd = wd$zm_gapit_FULL, trait_list = zm.pheno.FULL.listI, GD = data.frame(sp.zm.null.full[, 1], zm.just.geno.full), GM = zm.map[, c(1, 3, 4)], CV = zm_covar_full, K = myKI.zm.full)

#####
zm.pheno.FULL.gapit.listX <- list(zm.vpheno.molike.MAF40.h263.vQTN10.FULL)

names(zm.pheno.FULL.gapit.listX)[[1]] <- "zm.vpheno.molike.MAF40.h263.vQTN10.FULL"

GAPIT.expressway(wd = wd$zm_gapit_FULL, trait_list = zm.pheno.FULL.gapit.listX, GD = data.frame(sp.zm.null.full[, 1], zm.just.geno.full), GM = zm.map[, c(1, 3, 4)], CV = zm_covar_full, K = myKI.zm.full)

############################ Maize 500 section ##############################
#We have to make sure the taxa match for downstream GAPIT
sp.zm.GxE.500.mod$`<Trait>` <- sp.zm.null.500.filt$`<Trait>`
zm.vpheno.molike.MAF10.h233.vQTN10.500$taxa <- sp.zm.null.500.filt$`<Trait>`
zm.vpheno.molike.MAF10.h233.vQTN50.500$taxa <- sp.zm.null.500.filt$`<Trait>`
zm.vpheno.molike.MAF10.h233.vQTN90.500$taxa <- sp.zm.null.500.filt$`<Trait>`
zm.vpheno.molike.MAF10.h263.vQTN10.500$taxa <- sp.zm.null.500.filt$`<Trait>`
zm.vpheno.molike.MAF40.h263.vQTN10.500$taxa <- sp.zm.null.500.filt$`<Trait>`

#I will create the setting list
zm.pheno.500.listI <- list(
  sp.zm.GxE.500.mod,
  zm.vpheno.molike.MAF10.h233.vQTN10.500,
  zm.vpheno.molike.MAF10.h233.vQTN50.500,
  zm.vpheno.molike.MAF10.h233.vQTN90.500,
  zm.vpheno.molike.MAF10.h263.vQTN10.500,
  zm.vpheno.molike.MAF40.h263.vQTN10.500
)

#And give the elements names
names(zm.pheno.500.listI)[[1]] <- "sp.zm.GxE.500.mod"
names(zm.pheno.500.listI)[[2]] <- "zm.vpheno.molike.MAF10.h233.vQTN10.500"
names(zm.pheno.500.listI)[[3]] <- "zm.vpheno.molike.MAF10.h233.vQTN50.500"
names(zm.pheno.500.listI)[[4]] <- "zm.vpheno.molike.MAF10.h233.vQTN90.500"
names(zm.pheno.500.listI)[[5]] <- "zm.vpheno.molike.MAF10.h263.vQTN10.500"
names(zm.pheno.500.listI)[[6]] <- "zm.vpheno.molike.MAF40.h263.vQTN10.500"

wd$zm_gapit_500 <- "~/projects/dissertation/RII/2_pipeline/out/GAPIT/zm/500/"

GAPIT.expressway(wd = wd$zm_gapit_500, trait_list = zm.pheno.500.listI, GD = data.frame(sp.zm.null.500.filt[, 1], zm.just.geno.500), GM = zm.map[, c(1, 3, 4)], CV = zm_covar_500, K = myKI.zm.500)

##############
zm.sp.epi.higheffectsize.h230.500$`<Trait>` <- sp.zm.null.500.filt$`<Trait>`
zm.sp.epi.higheffectsize.500$taxa <- sp.zm.null.500.filt$`<Trait>`
zm.vpheno.molike.MAF10.h263.vQTN50.500$taxa <- sp.zm.null.500.filt$`<Trait>`
zm.vpheno.molike.MAF10.h263.vQTN90.500$taxa <- sp.zm.null.500.filt$`<Trait>`
zm.vpheno.molike.MAF40.h233.vQTN10.500$taxa <- sp.zm.null.500.filt$`<Trait>`
zm.vpheno.molike.MAF40.h233.vQTN50.500$taxa <- sp.zm.null.500.filt$`<Trait>`
zm.vpheno.molike.MAF40.h233.vQTN90.500$taxa <- sp.zm.null.500.filt$`<Trait>`
zm.vpheno.molike.MAF40.h263.vQTN50.500$taxa <- sp.zm.null.500.filt$`<Trait>`
zm.vpheno.molike.MAF40.500$taxa <- sp.zm.null.500.filt$`<Trait>`
zm.sp.epi.higheffectsize.h230.500$`<Trait>` <- sp.zm.null.500.filt$`<Trait>`
zm.sp.epi.higheffectsize.500$`<Trait>` <- sp.zm.null.500.filt$`<Trait>`

#This is the second list for the 500 maize diversity panel
zm.pheno.500.listII <- list(
  zm.sp.epi.higheffectsize.h230.500,
  zm.sp.epi.higheffectsize.500,
  zm.vpheno.molike.MAF10.h263.vQTN50.500,
  zm.vpheno.molike.MAF10.h263.vQTN90.500,
  zm.vpheno.molike.MAF40.h233.vQTN10.500,
  zm.vpheno.molike.MAF40.h233.vQTN50.500,
  zm.vpheno.molike.MAF40.h233.vQTN90.500,
  zm.vpheno.molike.MAF40.h263.vQTN50.500,
  zm.vpheno.molike.MAF40.500
)

#Names for each element
names(zm.pheno.500.listII)[[1]] <- "zm.sp.epi.higheffectsize.h230.500"
names(zm.pheno.500.listII)[[2]] <- "zm.sp.epi.higheffectsize.500"
names(zm.pheno.500.listII)[[3]] <- "zm.vpheno.molike.MAF10.h263.vQTN50.500"
names(zm.pheno.500.listII)[[4]] <- "zm.vpheno.molike.MAF10.h263.vQTN90.500"
names(zm.pheno.500.listII)[[5]] <- "zm.vpheno.molike.MAF40.h233.vQTN10.500"
names(zm.pheno.500.listII)[[6]] <- "zm.vpheno.molike.MAF40.h233.vQTN50.500"
names(zm.pheno.500.listII)[[7]] <- "zm.vpheno.molike.MAF40.h233.vQTN90.500"
names(zm.pheno.500.listII)[[8]] <- "zm.vpheno.molike.MAF40.h263.vQTN50.500"
names(zm.pheno.500.listII)[[9]] <- "zm.vpheno.molike.MAF40.500"
wd$zm_gapit_500 <- "~/projects/dissertation/RII/2_pipeline/out/GAPIT/zm/500/"

GAPIT.expressway(wd = wd$zm_gapit_500, trait_list = zm.pheno.500.listII, GD = data.frame(sp.zm.null.500.filt[, 1], zm.just.geno.500), GM = zm.map[, c(1, 3, 4)], CV = zm_covar_500, K = myKI.zm.500)


######################
save.image("~/vGWAS_chaptI_RII.RData")
