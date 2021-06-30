#########################vGWAS-GAPIT SCRIPT#################
setwd("~/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/2_pipeline/0_load_data/out")

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
zm_just_geno_500 <- t(zm_just_geno_500)
colnames(zm_just_geno_500) <- zm_map$snp

zm_just_geno_full <- t(zm_just_geno_full)
colnames(zm_just_geno_full) <- zm_map$snp

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
  Y=data.frame(at_sp_null_full[, c(1, 2)]), 
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
myGAPIT_zm_500 <- GAPIT(
  Y = data.frame(sp_zm_null_500[, c(1, 2)]), 
  GD = data.frame(sp_zm_null_500[, 1], zm_just_geno_500),
  GM = zm_map[, c(1, 3, 4)],
  PCA.total = 10, 
  SNP.effect = "Add"
) 

##full diversity panel
setwd(wd$zm_GAPIT_full)
myGAPIT_zm_full <- GAPIT(
  Y = data.frame(sp_zm_null_full[, c(1, 2)]), 
  GD = data.frame(sp_zm_null_full[, 1], zm_just_geno_full),
  GM = zm_map[, c(1, 3, 4)],
  PCA.total = 10, 
  SNP.effect = "Add"
)

zm_covar_500 <- myGAPIT_zm_500[["PCA"]]
zm_covar_full <- myGAPIT_zm_full[["PCA"]]