############################### DGLM Code #######################

#First, we will need to upload the required packages
install.packages("foreach", "parallel", "dglm", "data.table", "dplyr", "doParallel")

#Then, we will load the required libraries
library(foreach)
library(parallel)
library(dglm)
library(data.table)
library(dplyr)
library(doParallel)

#Establish Working directory roots
dir.create(paste0(wd$vGWAS_chapter_one_at), "vGWAS_results/DGLM/")
wd$at_dglm <- paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/DGLM/") #DGLM root

#null
dir.create(paste0(wd$at_dglm, "null/"))
wd$at_null_dglm <- paste0(wd$at_dglm, "null/")

##500 diversity panel
dir.create(paste0(wd$at_null_dglm, "500/"))
wd$at_null_dglm_500 <- paste0(wd$at_null_dglm, "500/")

##full diversity panel
dir.create(paste0(wd$at_null_dglm, "full/"))
wd$at_null_dglm_full <- paste0(wd$at_null_dglm, "full/")

#molike wmQTL
dir.create(paste0(wd$at_dglm, "molike_wmQTL/"))
wd$at_molike_wmQTL_dglm <- paste0(wd$at_dglm, "molike_wmQTL/")

##500 diversity panel
dir.create(paste0(wd$at_molike_wmQTL_dglm, "500/"))
wd$at_molike_wmQTL_dglm_500 <- paste0(wd$at_molike_wmQTL_dglm, "500/")

##full diversity panel
dir.create(paste0(wd$at_molike_wmQTL_dglm, "full/"))
wd$at_molike_wmQTL_dglm_full <- paste0(wd$at_molike_wmQTL_dglm, "full/")

#molike womQTL
dir.create(paste0(wd$at_dglm, "molike_womQTL/"))
wd$at_molike_womQTL_dglm <- paste0(wd$at_dglm, "molike_womQTL/")

##500 diversity panel
dir.create(paste0(wd$at_molike_womQTL_dglm, "500/"))
wd$at_molike_womQTL_dglm_500 <- paste0(wd$at_molike_womQTL_dglm, "500/")

##full diversity panel
dir.create(paste0(wd$at_molike_womQTL_dglm, "full/"))
wd$at_molike_womQTL_dglm_full <- paste0(wd$at_molike_womQTL_dglm, "full/")

#oillike womQTL MAF10
dir.create(paste0(wd$at_dglm, "oillike_womQTL_MAF10/"))
wd$at_oillike_womQTL_MAF10_dglm <- paste0(wd$at_dglm, "oillike_womQTL_MAF10/")

##500 diversity panel
dir.create(paste0(wd$at_oillike_womQTL_MAF10_dglm, "500/"))
wd$at_oillike_womQTL_MAF10_dglm_500 <- paste0(wd$at_oillike_womQTL_MAF10_dglm, "500/")

##full diversity panel
dir.create(paste0(wd$at_oillike_womQTL_MAF10_dglm, "full/"))
wd$at_oillike_womQTL_MAF10_dglm_full <- paste0(wd$at_oillike_womQTL_MAF10_dglm, "full/")

#oillike womQTL MAF40
dir.create(paste0(wd$at_dglm, "oillike_womQTL_MAF40/"))
wd$at_oillike_womQTL_MAF40_dglm <- paste0(wd$at_dglm, "oillike_womQTL_MAF40/")

##500 diversity panel
dir.create(paste0(wd$at_oillike_womQTL_MAF40_dglm, "500/"))
wd$at_oillike_womQTL_MAF40_dglm_500 <- paste0(wd$at_oillike_womQTL_MAF40_dglm, "500/")

##full diversity panel
dir.create(paste0(wd$at_oillike_womQTL_MAF40_dglm, "full/"))
wd$at_oillike_womQTL_MAF40_dglm_full <- paste0(wd$at_oillike_womQTL_MAF40_dglm, "full/")

###############This next section is for maize
dir.create(paste0(wd$vGWAS_chapter_one_zm), "vGWAS_results/DGLM/")
wd$zm_dglm <- paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/DGLM/") #DGLM root

#null
dir.create(paste0(wd$zm_dglm, "null/"))
wd$zm_null_dglm <- paste0(wd$zm_dglm, "null/")

##500 diversity panel
dir.create(paste0(wd$zm_null_dglm, "500/"))
wd$zm_null_dglm_500 <- paste0(wd$zm_null_dglm, "500/")

##full diversity panel
dir.create(paste0(wd$zm_null_dglm, "full/"))
wd$zm_null_dglm_full <- paste0(wd$zm_null_dglm, "full/")

#molike wmQTL
dir.create(paste0(wd$zm_dglm, "molike_wmQTL/"))
wd$zm_molike_wmQTL_dglm <- paste0(wd$zm_dglm, "molike_wmQTL/")

##500 diversity panel
dir.create(paste0(wd$zm_molike_wmQTL_dglm, "500/"))
wd$zm_molike_wmQTL_dglm_500 <- paste0(wd$zm_molike_wmQTL_dglm, "500/")

##full diversity panel
dir.create(paste0(wd$zm_molike_wmQTL_dglm, "full/"))
wd$zm_molike_wmQTL_dglm_full <- paste0(wd$zm_molike_wmQTL_dglm, "full/")

#molike without mQTL
dir.create(paste0(wd$zm_dglm, "molike_womQTL/"))
wd$zm_molike_womQTL_dglm <- paste0(wd$zm_dglm, "molike_womQTL/")

##500 diversity panel
dir.create(paste0(wd$zm_molike_womQTL_dglm, "500/"))
wd$zm_molike_womQTL_dglm_500 <- paste0(wd$zm_molike_womQTL_dglm, "500/")

##full diversity panel
dir.create(paste0(wd$zm_molike_womQTL_dglm, "full/"))
wd$zm_molike_womQTL_dglm_full <- paste0(wd$zm_molike_womQTL_dglm, "full/")

#oillike womQTL MAF10
dir.create(paste0(wd$zm_dglm, "oillike_womQTL_MAF10/"))
wd$zm_oillike_womQTL_MAF10_dglm <- paste0(wd$zm_dglm, "oillike_womQTL_MAF10/")

##500 diversity panel
dir.create(paste0(wd$zm_oillike_womQTL_MAF10_dglm, "500/"))
wd$zm_oillike_womQTL_MAF10_dglm_500 <- paste0(wd$zm_oillike_womQTL_MAF10_dglm, "500/")

##full diversity panel
dir.create(paste0(wd$zm_oillike_womQTL_MAF10_dglm, "full/"))
wd$zm_oillike_womQTL_MAF10_dglm_full <- paste0(wd$zm_oillike_womQTL_MAF10_dglm, "full/")

#oillike womQTL MAF40
dir.create(paste0(wd$zm_dglm, "oillike_womQTL_MAF40/"))
wd$zm_oillike_womQTL_MAF40_dglm <- paste0(wd$zm_dglm, "oillike_womQTL_MAF40/")

##500 diversity panel
dir.create(paste0(wd$zm_oillike_womQTL_MAF40_dglm, "500/"))
wd$zm_oillike_womQTL_MAF40_dglm_500 <- paste0(wd$zm_oillike_womQTL_MAF40_dglm, "500/")

##full diversity panel
dir.create(paste0(wd$zm_oillike_womQTL_MAF40_dglm, "full/"))
wd$zm_oillike_womQTL_MAF40_dglm_full <- paste0(wd$zm_oillike_womQTL_MAF40_dglm, "full/")

####To get DGLM to run faster, I will be utilizing the foreach R package in conjuction with 
#doParallel

#We will see how many cores we can use with detectCores()
detectCores()

#Then, we need to register the cluster
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

###Load PCA data for Arabidopsis and maize
####Arabidopsis
##500
wd$at_GAPIT_500 <- paste0(wd$vGWAS_chapter_one_at, "at_GAPIT_500/")
at_covar_500 <- read.csv(paste0(wd$at_GAPIT_500, "GAPIT.PCA.csv"))
at_covar_500 <- at_covar_500[, c(1:4)]

##full
wd$at_GAPIT_full <- paste0(wd$vGWAS_chapter_one_at, "at_GAPIT_full/")
at_covar_full <- read.csv(paste0(wd$at_GAPIT_full, "GAPIT.PCA.csv"))
at_covar_full <- at_covar_full[, c(1:4)]

###Maize
##500
wd$zm_GAPIT_500 <- paste0(wd$vGWAS_chapter_one_zm, "zm_GAPIT_500/")
zm_covar_500 <- read.csv(paste0(wd$zm_GAPIT_500, "GAPIT.PCA.csv"))
zm_covar_500 <- zm_covar_500[, c(1:5)]

##full
wd$zm_GAPIT_full <- paste0(wd$vGWAS_chapter_one_zm, "zm_GAPIT_full/")
zm_covar_full <- read.csv(paste0(wd$zm_GAPIT_full, "GAPIT.PCA.csv"))
zm_covar_full <- zm_covar_full[, c(1:5)]

##Now, I will define the DGLM function. To get the foreach loop to work, I will need to double 
#define functions: One to run DGLM and the other to run in foreach

DGLM.maxima <- function(Geno = NULL, Phenos = NULL, mQTN = NULL, covar = NULL, Map = NULL, trait.name = NULL) {
  my.pdglm <- function(cT = NULL, i = NULL, Phenos = NULL, Geno = NULL, mQTN = NULL, covar = NULL, Map = NULL) {
    print(paste("--------- Fitting DGLM model for SNP ", i, " out of ", dim(Geno)[2], "  ----------", sep = ""))
    temp_123456789 <<- data.frame(y = Phenos[, cT], snp = Geno[, i], mQTL = Geno[, c(mQTN)], covar[, -1])
    #you dont need to change the formula anymore, it includes however number of PCs you have
    model <-
      dglm(
        formula = as.formula(paste("y ~ snp + ", paste0(colnames(temp_123456789[, -1:-2]), collapse = "+"))),
        ~ snp, data = temp_123456789,
        family = gaussian(link = "identity")
      )
    Map_info <- Map[i, ]
    P.mean <- summary(model)$coef["snp", "Pr(>|t|)"]  # Extarct p values for mean part
    an <- anova(model)
    P.disp <- pchisq(q = an["Dispersion model", "Adj.Chisq"], df = an["Dispersion model", "DF"], lower.tail = FALSE)
    s.model <- summary(model$dispersion.fit)
    beta <- s.model$coef["snp", "Estimate"]  # Extarct cofficients
    se <- s.model$coef["snp", "Std. Error"]  # Extract standard errors
    out <- data.frame(Map_info = Map_info, Beta = beta, SE = se, P.mean = P.mean, P.disp = P.disp, stringsAsFactors = FALSE)  # Save all the extracted variables in data frame out
    rm(list = "temp_123456789", envir = .GlobalEnv)
    return(out)
    print(out)
  }
  print(paste0("Welcome to DGLM Maxima"))
  print(paste0("Starting foreach loops"))
  for(cT in 2:dim(Phenos)[2]) {
    print(cT)
    results <- foreach(i = 1:dim(Geno)[2], .combine = 'rbind', .packages = 'dglm') %dopar%
      my.pdglm(cT = cT, i = i, Phenos = Phenos, Geno = Geno, mQTN = mQTN, covar = covar, Map = Map)
    results$fdr.pmean <- p.adjust(results$P.mean, "fdr")
    results$fdr.pdisp <- p.adjust(results$P.disp, "fdr")
    results <- results[order(results$fdr.pdisp), ]
    fwrite(results, file = paste0("DGLM_", trait.name, colnames(Phenos[cT]), ".csv", sep=""), sep=",", row.names = FALSE) 
    rm(results)
  }
}

##########################################################################
############################### Analyzes #################################
##########################################################################

#Now, we can run our analyzes
###Arabidopsis
####NULL
#####500
setwd(wd$at_null_dglm_500)
DGLM.maxima(Geno = at_just_geno_500, Phenos = at_sp_null_500, mQTN = NULL, covar = at_covar_500, Map = at_map[, c(1, 3, 4)], trait.name = "at_null_500")

####NULL
#####FULL
setwd(wd$at_null_dglm_full)
DGLM.maxima(Geno = at_just_geno_full, Phenos = at_sp_null_full, mQTN = NULL, covar = at_covar_full, Map = at_map[, c(1, 3, 4)], trait.name = "at_null_full")

####Molike without mean QTL maf40
#####500
setwd(wd$at_molike_womQTL_dglm_500)
DGLM.maxima(Geno = at_just_geno_500, Phenos = at_vpheno_molike_womQTL_MAF40_500, mQTN = NULL, covar = at_covar_500, Map = at_map[, c(1, 3, 4)], trait.name = "at_molike_womQTL_500")

#####FULL
setwd(wd$at_molike_womQTL_dglm_full)
DGLM.maxima(Geno = at_just_geno_full, Phenos = at_vpheno_molike_womQTL_MAF40_full, mQTN = NULL, covar = at_covar_full, Map = at_map[, c(1, 3, 4)], trait.name = "at_molike_womQTL_full")

####Molike with mean QTL maf40
#####500
setwd(wd$at_molike_wmQTL_dglm_500)
DGLM.maxima(Geno = at_just_geno_500, Phenos = at_vpheno_molike_MAF40_500, mQTN = 29513, covar = at_covar_500, Map = at_map[, c(1, 3, 4)], trait.name = "at_molike_wmQTL_500")

#####full
setwd(wd$at_molike_wmQTL_dglm_full)
DGLM.maxima(Geno = at_just_geno_full, Phenos = at_vpheno_molike_MAF40_full, mQTN = 16692, covar = at_covar_full, Map = at_map[, c(1, 3, 4)], trait.name = "at_molike_wmQTL_full")

####Oillike_womQTL_MAF10
#####500
setwd(wd$at_oillike_womQTL_MAF10_dglm_500)
DGLM.maxima(Geno = at_just_geno_500, Phenos = at_vpheno_oil_500_womQTL_MAF10, covar = at_covar_500, Map = At_map[, c(1, 3, 4)], trait.name = "At_oillike_womQTL_MAF10_500")

#####FULL
setwd(wd$at_oillike_womQTL_MAF10_dglm_full)
DGLM.maxima(Geno = at_just_geno_full, Phenos = at_vpheno_oil_full_womQTL_MAF10, covar = at_covar_full, Map = At_map[, c(1, 3, 4)], trait.name = "At_oillike_womQTL_MAF10_full")

####Oillike_womQTL_MAF40
#####500
setwd(wd$at_oillike_womQTL_MAF40_dglm_500)
DGLM.maxima(Geno = at_just_geno_500, Phenos = at_vpheno_oil_500_womQTL_MAF40, covar = at_covar_500, Map = At_map[, c(1, 3, 4)], trait.name = "At_oillike_womQTL_MAF40_500")

#####FULL
setwd(wd$at_oillike_womQTL_MAF40_dglm_full)
DGLM.maxima(Geno = at_just_geno_full, Phenos = at_vpheno_oil_full_womQTL_MAF40, covar = at_covar_full, Map = At_map[, c(1, 3, 4)], trait.name = "At_oillike_womQTL_MAF40_full")

############ Maize ##############3
####NULL
#####500
setwd(wd$zm_null_dglm_500)
DGLM.maxima(Geno = zm_just_geno_500, Phenos = sp_zm_null_500, mQTN = NULL, covar = zm_covar_500, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_null_500")

####NULL
#####FULL
setwd(wd$zm_null_dglm_full)
DGLM.maxima(Geno = zm_just_geno_full, Phenos = sp_zm_null_full, mQTN = NULL, covar = zm_covar_full, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_null_full")

####Molike without mean QTL maf40
#####500
setwd(wd$zm_molike_womQTL_dglm_500)
DGLM.maxima(Geno = zm_just_geno_500, Phenos = zm_vpheno_molike_womQTL_MAF40_500, mQTN = NULL, covar = zm_covar_500, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_molike_womQTL_500")

#####FULL
setwd(wd$zm_molike_womQTL_dglm_full)
DGLM.maxima(Geno = zm_just_geno_full, Phenos = zm_vpheno_molike_womQTL_MAF40_full, mQTN = NULL, covar = zm_covar_full, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_molike_womQTL_full")

####Molike with mean QTL maf40
#####500
setwd(wd$zm_molike_wmQTL_dglm_500)
DGLM.maxima(Geno = zm_just_geno_500, Phenos = zm_vpheno_molike_MAF40_500, mQTN = 30995, covar = zm_covar_500, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_molike_wmQTL_500")

#####full
setwd(wd$zm_molike_wmQTL_dglm_full)
DGLM.maxima(Geno = zm_just_geno_full, Phenos = zm_vpheno_molike_MAF40_full, mQTN = 14328, covar = zm_covar_full, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_molike_wmQTL_full")

####Oillike_womQTL_MAF10
#####500
setwd(wd$zm_oillike_womQTL_MAF10_dglm_500)
DGLM.maxima(Geno = zm_just_geno_500, Phenos = zm_vpheno_oil_500_womQTL_MAF10, mQTN = NULL, covar = zm_covar_500, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_oillike_womQTL_MAF10_500")

#####FULL
setwd(wd$zm_oillike_womQTL_MAF10_dglm_full)
DGLM.maxima(Geno = zm_just_geno_full, Phenos = zm_vpheno_oil_full_womQTL_MAF10, mQTN = NULL, covar = zm_covar_full, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_oillike_womQTL_MAF10_full")

####Oillike_womQTL_MAF40
#####500
setwd(wd$zm_oillike_womQTL_MAF40_dglm_500)
DGLM.maxima(Geno = zm_just_geno_500, Phenos = zm_vpheno_oil_500_womQTL_MAF40, mQTN = NULL, covar = zm_covar_500, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_oillike_womQTL_MAF40_500")

#####FULL
setwd(wd$zm_oillike_womQTL_MAF40_dglm_full)
DGLM.maxima(Geno = zm_just_geno_full, Phenos = zm_vpheno_oil_full_womQTL_MAF40, mQTN = NULL, covar = zm_covar_full, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_oillike_womQTL_MAF40_full")