#############################################################################
#################### DGLM Script ########################
#########################################################

#As always, load the required libraries
library(dglm)
library(foreach)
library(parallel)
library(data.table)

#Analogous to the BFT and GAPIT scripts, I am going to 
#create a DGLM script that iterates through a list
DGLM.maxima2 <- function(wd = NULL, Geno = NULL, trait_list = NULL, covar = NULL, Map = NULL) {
  my.pdglm <- function(cT = NULL, i = NULL, Phenos = NULL, Geno = NULL, covar = NULL, Map = NULL) {
    print(paste("--------- Fitting DGLM model for SNP ", i, " out of ", dim(Geno)[2], "  ----------", sep = ""))
    temp_123456789 <<- data.frame(y = Phenos[, cT], snp = Geno[, i], covar[, -1])
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
  for (h in 1:length(trait_list)) {
    print(paste0("this is trait", h))
    new_folder <- dir.create(paste0(wd, names(trait_list)[[h]]))
    setwd(paste0(wd, names(trait_list)[[h]]))
    Phenos <- trait_list[[h]]
    for (cT in 2:dim(Phenos)[2]) {
      print(cT)
      results <- foreach(i = 1:dim(Geno)[2], .combine = 'rbind', .errorhandling = 'remove', .packages = 'dglm') %dopar%
        my.pdglm(cT = cT, i = i, Phenos = Phenos, Geno = Geno, covar = covar, Map = Map)
      results$fdr.pmean <- p.adjust(results$P.mean, "fdr")
      results$fdr.pdisp <- p.adjust(results$P.disp, "fdr")
      results <- results[order(results$fdr.pdisp), ]
      fwrite(results, file = paste0("DGLM_", names(trait_list)[[h]], colnames(Phenos[cT]), ".csv", sep = ""), sep=",", row.names = FALSE) 
      rm(results)
    }
  }
}

##From here, I will subset the FULL list based on th covariate object needed
setwd("~/Chapter_one/2_pipeline/out/DGLM_Results/Zm/FULLII")
wd$zm_DGLM_FULL <- "~/projects/dissertation/RII/2_pipeline/out/DGLM/FULL/"

#Trait I. Because all of these settings use the same mQTn,
#I am going to make a new covariate object with this mQTN as
#a covariate
grep("S8_6701572", colnames(zm.just.geno.full))
zm_covar_full_molike <- data.frame(zm_covar_full, zm.just.geno.full[, 61805])
DGLM.maxima2(wd = wd$zm_DGLM_FULL, Geno = zm.just.geno.full, trait_list = zm.pheno.FULL.listI, covar = zm_covar_full_molike, Map = zm.map[, c(1, 3, 4)])

setwd("~/Chapter_one/2_pipeline/out/DGLM_Results/Zm/FULLII")
wd$zm_DGLM_FULLII <- "~/Chapter_one/2_pipeline/out/DGLM_Results/Zm/FULLII/"

DGLM.maxima2(wd = wd$zm_DGLM_FULLII, Geno = t(zm_just_geno_FULL_filtlike_at), trait_list = zm_settings_FULL_list_molike, covar = zm_covar_FULL_filtlike_at_molike, Map = zm_map_filtlike_at[, c(1, 3, 4)])

###MAF 40 list
zm_settings_FULL_list_molike40 <- list(zm_vpheno_molike_MAF40_h233_vQTN10_FULL_filtlike_at, 
                                       zm_vpheno_molike_MAF40_h233_vQTN50_FULL_filtlike_at,
                                       zm_vpheno_molike_MAF40_h233_vQTN90_FULL_filtlike_at,
                                       zm_vpheno_molike_MAF40_h263_vQTN10_FULL_filtlike_at,
                                       zm_vpheno_molike_MAF40_h263_vQTN50_FULL_filtlike_at,
                                       zm_vpheno_molike_MAF40_FULL_filtlike_at)

#Then, I will give each element its own unique name that corresponds to its setting
names(zm_settings_FULL_list_molike40)[[1]] <- "zm_vpheno_molike_MAF40_h233_vQTN10_FULL_filtlike_at"
names(zm_settings_FULL_list_molike40)[[2]] <- "zm_vpheno_molike_MAF40_h233_vQTN50_FULL_filtlike_at"
names(zm_settings_FULL_list_molike40)[[3]] <- "zm_vpheno_molike_MAF40_h233_vQTN90_FULL_filtlike_at"
names(zm_settings_FULL_list_molike40)[[4]] <- "zm_vpheno_molike_MAF40_h263_vQTN10_FULL_filtlike_at"
names(zm_settings_FULL_list_molike40)[[5]] <- "zm_vpheno_molike_MAF40_h263_vQTN50_FULL_filtlike_at"
names(zm_settings_FULL_list_molike40)[[6]] <- "zm_vpheno_molike_MAF40_FULL_filtlike_at"

setwd("~/Chapter_one/2_pipeline/out/DGLM_Results/Zm/FULLII")
wd$zm_DGLM_FULLII <- "~/Chapter_one/2_pipeline/out/DGLM_Results/Zm/FULLII/"

DGLM.maxima2(wd = wd$zm_DGLM_FULLII, Geno = t(zm_just_geno_FULL_filtlike_at), trait_list = zm_settings_FULL_list_molike40, covar = zm_covar_FULL_filtlike_at_molike, Map = zm_map_filtlike_at[, c(1, 3, 4)])

##This next list is for GxE and Null. The Null setting below did not complete to running due to
#memory issues
zm.pheno.FULL.listVI <- list(sp.zm.null.full[, c(1, 30:101)],
                             sp.zm.GxE.FULL.mod)

names(zm.pheno.FULL.listVI)[[1]] <- "sp.zm.null.full"
names(zm.pheno.FULL.listVI)[[2]] <- "sp.zm.GxE.FULL.mod"

#check dimensions of CV before running
dim(zm_covar_full)
DGLM.maxima2(wd = wd$zm_DGLM_FULL, Geno = zm.just.geno.full, trait_list = zm.pheno.FULL.listVI, covar = zm_covar_full, Map = zm.map[, c(1, 3, 4)])

#another setting list
zm.pheno.FULL.listVIII <- list(zm.vpheno.molike.MAF40.h233.vQTN10.FULL,
                               zm.vpheno.molike.MAF40.h233.vQTN50.FULL,
                               zm.vpheno.molike.MAF40.h233.vQTN90.FULL,
                               zm.vpheno.molike.MAF40.h263.vQTN10.FULL)

names(zm.pheno.FULL.listVIII)[[1]] <- "zm.vpheno.molike.MAF40.h233.vQTN10.FULL"
names(zm.pheno.FULL.listVIII)[[2]] <- "zm.vpheno.molike.MAF40.h233.vQTN50.FULL"
names(zm.pheno.FULL.listVIII)[[3]] <- "zm.vpheno.molike.MAF40.h233.vQTN90.FULL"
names(zm.pheno.FULL.listVIII)[[4]] <- "zm.vpheno.molike.MAF40.h263.vQTN10.FULL"

#
DGLM.maxima2(wd = wd$zm_DGLM_FULL, Geno = zm.just.geno.full, trait_list = zm.pheno.FULL.listVIII, covar = zm_covar_full_molike, Map = zm.map[, c(1, 3, 4)])

#####This next section is for the 500 diversity panel in maize
zm.pheno.500.listI.dglm <- list(sp.zm.null.500.filt,
                                sp.zm.GxE.500.mod)

names(zm.pheno.500.listI.dglm)[[1]] <- "sp.zm.null.500.filt"
names(zm.pheno.500.listI.dglm)[[2]] <- "sp.zm.GxE.500.mod"

#Path for 500 dglm in maize
wd$zm_dglm_500 <- "~/projects/dissertation/RII/2_pipeline/out/DGLM/500/"

zm_covar_500 <- zm_covar_500[, c(1:4)]

DGLM.maxima2(wd = wd$zm_dglm_500, Geno = zm.just.geno.500, trait_list = zm.pheno.500.listI.dglm, covar = zm_covar_500, Map = zm.map[, c(1, 3, 4)])

##Rep19 is no good
zm.pheno.500.listI.dglm <- list(sp.zm.null.500.filt[, c(1, 46:101)],
                                sp.zm.GxE.500.mod)


names(zm.pheno.500.listI.dglm)[[1]] <- "sp.zm.null.500.filt"
names(zm.pheno.500.listI.dglm)[[2]] <- "sp.zm.GxE.500.mod"

DGLM.maxima2(wd = wd$zm_dglm_500, Geno = zm.just.geno.500, trait_list = zm.pheno.500.listI.dglm, covar = zm_covar_500, Map = zm.map[, c(1, 3, 4)])

###Rep 6 clunked out for GxE
zm.pheno.500.listI.dglm <- list(
  sp.zm.GxE.500.mod[, c(1, 8:101)])

names(zm.pheno.500.listI.dglm)[[1]] <- "sp.zm.GxE.500.mod"

DGLM.maxima2(wd = wd$zm_dglm_500, Geno = zm.just.geno.500, trait_list = zm.pheno.500.listI.dglm, covar = zm_covar_500, Map = zm.map[, c(1, 3, 4)])

###Rep 19 clunked out
zm.pheno.500.listI.dglm <- list(
  sp.zm.GxE.500.mod[, c(1, 20, 83:101)])

names(zm.pheno.500.listI.dglm)[[1]] <- "sp.zm.GxE.500.mod"

DGLM.maxima2(wd = wd$zm_dglm_500, Geno = zm.just.geno.500, trait_list = zm.pheno.500.listI.dglm, covar = zm_covar_500, Map = zm.map[, c(1, 3, 4)])

####Rerunning NULL for reps that clunked out
zm.pheno.500.listI.dglm <- list(
  sp.zm.null.500.filt[, c(1, 20, 27, 45)])

names(zm.pheno.500.listI.dglm)[[1]] <- "sp.zm.null.500.filt"

DGLM.maxima2(wd = wd$zm_dglm_500, Geno = zm.just.geno.500, trait_list = zm.pheno.500.listI.dglm, covar = zm_covar_500, Map = zm.map[, c(1, 3, 4)])

####
zm.pheno.500.listII.dglm <- list(zm.sp.epi.higheffectsize.h230.500[, c(1, 50:101)],
                                 zm.sp.epi.higheffectsize.500)

#names for this list
names(zm.pheno.500.listII.dglm)[[1]] <- "zm.sp.epi.higheffectsize.h230.500"
names(zm.pheno.500.listII.dglm)[[2]] <- "zm.sp.epi.higheffectsize.500"

zm_covar_500 <- zm_covar_500[, c(1:4)]
DGLM.maxima2(wd = wd$zm_dglm_500, Geno = zm.just.geno.500, trait_list = zm.pheno.500.listII.dglm, covar = zm_covar_500, Map = zm.map[, c(1, 3, 4)])

####Reps that did not run for the lowh2 hertiability setting
zm.pheno.500.listII.dglm <- list(zm.sp.epi.higheffectsize.h230.500[, c(1, 25)])

#name of this list
names(zm.pheno.500.listII.dglm)[[1]] <- "zm.sp.epi.higheffectsize.h230.500"

zm_covar_500 <- zm_covar_500[, c(1:4)]
DGLM.maxima2(wd = wd$zm_dglm_500, Geno = zm.just.geno.500, trait_list = zm.pheno.500.listII.dglm, covar = zm_covar_500, Map = zm.map[, c(1, 3, 4)])

##
zm.pheno.500.listI.dglm <- list(
  sp.zm.null.500.filt[, c(1, 26)])

names(zm.pheno.500.listI.dglm)[[1]] <- "sp.zm.null.500.filt"

DGLM.maxima2(wd = wd$zm_dglm_500, Geno = zm.just.geno.500, trait_list = zm.pheno.500.listI.dglm, covar = zm_covar_500, Map = zm.map[, c(1, 3, 4)])

###Redoign Molike
zm_covar_500 <- zm_covar_500[, c(1:4)]
grep("S2_19265645", colnames(zm.just.geno.500))
zm_covar_500_molike <- data.frame(zm_covar_500, zm.just.geno.500[, 22449])

zm.pheno.500.listIII.dglm <- list(zm.vpheno.molike.MAF10.h233.vQTN10.500,
                                  zm.vpheno.molike.MAF10.h233.vQTN50.500,
                                  zm.vpheno.molike.MAF10.h233.vQTN90.500,
                                  zm.vpheno.molike.MAF10.h263.vQTN10.500,
                                  zm.vpheno.molike.MAF10.h263.vQTN50.500,
                                  zm.vpheno.molike.MAF10.h263.vQTN90.500)

names(zm.pheno.500.listIII.dglm)[[1]] <- "zm.vpheno.molike.MAF10.h233.vQTN10.500"
names(zm.pheno.500.listIII.dglm)[[2]] <- "zm.vpheno.molike.MAF10.h233.vQTN50.500"
names(zm.pheno.500.listIII.dglm)[[3]] <- "zm.vpheno.molike.MAF10.h233.vQTN90.500"
names(zm.pheno.500.listIII.dglm)[[4]] <- "zm.vpheno.molike.MAF10.h263.vQTN10.500"
names(zm.pheno.500.listIII.dglm)[[5]] <- "zm.vpheno.molike.MAF10.h263.vQTN50.500"
names(zm.pheno.500.listIII.dglm)[[6]] <- "zm.vpheno.molike.MAF10.h263.vQTN90.500"


DGLM.maxima2(wd = wd$zm_dglm_500, Geno = zm.just.geno.500, trait_list = zm.pheno.500.listIII.dglm, covar = zm_covar_500_molike, Map = zm.map[, c(1, 3, 4)])
