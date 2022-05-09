###################################################################
######################## ROC Curve Script ########################
##################################################################

###Load workspace
load("~/vGWAS_chaptI_RII.RData")

##As always, upload required libraries
library(dplyr)
library(pROC)

###Now, I will need to create a function that random selects 10 replications for each setting
#Here is what it needs to do
##1.) Select 10 random numbers
##2.) Read in BFT and DGLM files that correspond to above (Note this is not necessarily that rep!)
##3.) Take out simulated vQTN
##3.) Apply ROC
##4.) Feed simulated QTN at some point
##As far as the arguments are concerned
#repsel: this is for whatever you want to randomly select replications or input a previous run's random numbers
#randvec: if repsel = F, this is a vector of random numbers from a previous run
#BFT.path.root: This is the path that contains all BFT results
#DGLM.path,root: This is the path that contains all DGLM folders
#BFT.fold.nam: This is the name of the specific results folder for BFT
#DGLM.fold.nam: This is the name of the specifiv results folder for DGLM
#vQTN.info: This is the simulated vQTN information with physical window boundaries
#trait.name: The setting name that you are analysing

generate_ROC <- function(repsel = NULL, randvec = NULL, BFT.path.root = NULL, DGLM.path.root = NULL, BFT.fold.nam = NULL, DGLM.fold.nam = NULL, vQTN.info = NULL, trait.name = NULL) {
  if (is.null(repsel) == F) { #This if statement is for either randomly selecting replications or inputting previous reps
    randreps <- randvec
  } else {
    randreps <- sample(1:100, size = 10, replace = F)
    print(randreps)
  }
  AUC.mat <- matrix(NA, nrow = 10, ncol = 5)
  colnames(AUC.mat) <- c("randnumber", "BFT.file", "DGLM.file", "BFT.AUC", "DGLM.AUC")
  print(randreps)
  sim.snp <- vQTN.info$snp
  print(sim.snp)
  BFT.res.fold <- paste0(BFT.path.root, BFT.fold.nam, "/")
  DGLM.res.fold <- paste0(DGLM.path.root, DGLM.fold.nam, "/")
  print(DGLM.res.fold)
  BFT.files <- list.files(BFT.res.fold)
  DGLM.files <- list.files(DGLM.res.fold)
  for (i in 1:length(randreps)) { #This for loop is required to iterate through the files
    print(randreps[i])
    randrep <- randreps[i]
    AUC.mat[i, 1] <- randrep
    print(paste0(BFT.res.fold, BFT.files[c(randrep)]))
    AUC.mat[i, 2] <- BFT.files[c(randrep)]
    BFT.results <- read.csv(paste0(BFT.res.fold, BFT.files[c(randrep)]), sep = ",")
    remQTN <- grep(TRUE, BFT.results$snp == sim.snp)
    BFT.results <- BFT.results[-c(remQTN), ]
    print(paste0(DGLM.res.fold, DGLM.files[c(randrep)]))
    AUC.mat[i, 3] <- DGLM.files[c(randrep)]
    DGLM.results <- read.csv(paste0(DGLM.res.fold, DGLM.files[c(randrep)]), sep = ",")
    
    remQTNII <- grep(TRUE, DGLM.results$Map_info.snp == sim.snp)
    DGLM.results <- DGLM.results[-c(remQTNII), ]
    
    BFT.results$class <- 0
    DGLM.results$class <- 0
    BFT.trueQTN <- BFT.results[BFT.results$chr == vQTN.info$chr, ]
    BFT.trueQTN <- BFT.trueQTN[BFT.trueQTN$pos >= as.numeric(vQTN.info$QTNwindow_lowerbound), ]
    BFT.trueQTN <- BFT.trueQTN[BFT.trueQTN$pos <= as.numeric(vQTN.info$QTNwindow_upperbound), ]
    BFT.trueQTN$class <- 1
    id <- grep(TRUE, BFT.results$snp %in% BFT.trueQTN$snp)
    BFT.results[c(id), ] <- BFT.trueQTN
    
    DGLM.trueQTN <- DGLM.results[DGLM.results$Map_info.chr == vQTN.info$chr, ]
    DGLM.trueQTN <- DGLM.trueQTN[DGLM.trueQTN$Map_info.pos >= as.numeric(vQTN.info$QTNwindow_lowerbound), ]
    DGLM.trueQTN <- DGLM.trueQTN[DGLM.trueQTN$Map_info.pos <= as.numeric(vQTN.info$QTNwindow_upperbound), ]
    print(dim(DGLM.trueQTN))
    
    DGLM.trueQTN$class <- 1
    id <- grep(TRUE, DGLM.results$Map_info.snp %in% DGLM.trueQTN$Map_info.snp)
    DGLM.results[c(id), ] <- DGLM.trueQTN
    
    BFT.glm.fit <- glm(BFT.results$class ~ BFT.results$p.value, family = binomial)
    BFT.ROC <- roc(BFT.results$class, BFT.glm.fit$fitted.values, plot = T, legacy.axes = T, percent = T, xlab = "False Positive Percentage", ylab = "True Postive Percentage")
    AUC.mat[i, 4] <- BFT.ROC$auc
    
    DGLM.glm.fit <- glm(DGLM.results$class ~ DGLM.results$P.disp, family = binomial)
    DGLM.ROC <- roc(DGLM.results$class, DGLM.glm.fit$fitted.values, plot = T, legacy.axes = T, percent = T, xlab = "False Positive Percentage", ylab = "True Postive Percentage")
    AUC.mat[i, 5] <- DGLM.ROC$auc
 #   tiff(file = paste0(randrep, "_", trait.name, ".tiff"),
#         width = 800, height = 400, res = 100)
#      par(pty = "s")
#      roc(BFT.results$class, BFT.glm.fit$fitted.values, plot = T, legacy.axes = T, percent = T, xlab = "False Positive %", ylab = "True Postive %", col = "#377eb8", lwd = 2, print.auc = T, print.auc.x = 45)
#      plot.roc(DGLM.results$class, DGLM.glm.fit$fitted.values, percent = T, col = "#4daf4a", lwd = 2, print.auc = T, add = T, print.auc.y = 35)
#    dev.off()
  }
  write.csv(AUC.mat, paste0(trait.name, ".", "AUC.mat.test.csv"), row.names = F)
}

####Rerun with GAPIT
generate_ROC <- function(repsel = NULL, randvec = NULL, BFT.path.root = NULL, DGLM.path.root = NULL, MLM.path.root = NULL, BFT.fold.nam = NULL, DGLM.fold.nam = NULL, MLM.fold.nam = NULL, vQTN.info = NULL, trait.name = NULL) {
  if (is.null(repsel) == F) { #This if statement is for either randomly selecting replications or inputting previous reps
    randreps <- randvec
  } else {
    randreps <- sample(1:100, size = 10, replace = F)
    print(randreps)
  }
#  AUC.mat <- matrix(NA, nrow = 10, ncol = 7)
#  colnames(AUC.mat) <- c("randnumber", "BFT.file", "DGLM.file", "MLM.file", "BFT.AUC", "DGLM.AUC", "MLM.AUC")
  print(randreps)
#  sim.snp <- vQTN.info$snp
#  print(sim.snp)
  BFT.res.fold <- paste0(BFT.path.root, BFT.fold.nam, "/")
  DGLM.res.fold <- paste0(DGLM.path.root, DGLM.fold.nam, "/")
  MLM.res.fold <- paste0(MLM.path.root, MLM.fold.nam, "/")
  BFT.files <- list.files(BFT.res.fold)
  DGLM.files <- list.files(DGLM.res.fold)
  MLM.files <- list.files(MLM.res.fold)
  for (h in 1:nrow(vQTN.info)){
    print(h)
    vQTN.infoh <- vQTN.info[h, ]
    sim.snp <- vQTN.infoh$snp
    print(sim.snp)
    print(dim(vQTN.infoh))
    AUC.mat <- matrix(NA, nrow = 10, ncol = 7)
    colnames(AUC.mat) <- c("randnumber", "BFT.file", "DGLM.file", "MLM.file", "BFT.AUC", "DGLM.AUC", "MLM.AUC")
    for (i in 1:length(randreps)) { #This for loop is required to iterate through the files
      print(randreps[i])
      randrep <- randreps[i]
      AUC.mat[i, 1] <- randrep
      print(paste0(BFT.res.fold, BFT.files[c(randrep)]))
      AUC.mat[i, 2] <- BFT.files[c(randrep)]
      BFT.results <- read.csv(paste0(BFT.res.fold, BFT.files[c(randrep)]), sep = ",")
      remQTN <- grep(TRUE, BFT.results$snp == sim.snp)
      BFT.results <- BFT.results[-c(remQTN), ]
      print(paste0(DGLM.res.fold, DGLM.files[c(randrep)]))
      print(paste0(MLM.res.fold, MLM.files[c(randrep)]))
    
      AUC.mat[i, 3] <- DGLM.files[c(randrep)]
      DGLM.results <- read.csv(paste0(DGLM.res.fold, DGLM.files[c(randrep)]), sep = ",")
    
      remQTNII <- grep(TRUE, DGLM.results$Map_info.snp == sim.snp)
      DGLM.results <- DGLM.results[-c(remQTNII), ]
    
      AUC.mat[i, 4] <- MLM.files[c(randrep)]
      MLM.results <- read.csv(paste0(MLM.res.fold, MLM.files[c(randrep)]), sep = ",")
    
      remQTNIII <- grep(TRUE, MLM.results$SNP == sim.snp)
      MLM.results <- MLM.results[-c(remQTNIII), ]
    
      BFT.results$class <- 0
      DGLM.results$class <- 0
      MLM.results$class <- 0
    
      BFT.trueQTN <- BFT.results[BFT.results$chr == vQTN.info$chr, ]
      BFT.trueQTN <- BFT.trueQTN[BFT.trueQTN$pos >= as.numeric(vQTN.info$QTNwindow_lowerbound), ]
      BFT.trueQTN <- BFT.trueQTN[BFT.trueQTN$pos <= as.numeric(vQTN.info$QTNwindow_upperbound), ]
      BFT.trueQTN$class <- 1
      id <- grep(TRUE, BFT.results$snp %in% BFT.trueQTN$snp)
      BFT.results[c(id), ] <- BFT.trueQTN
    
      DGLM.trueQTN <- DGLM.results[DGLM.results$Map_info.chr == vQTN.info$chr, ]
      DGLM.trueQTN <- DGLM.trueQTN[DGLM.trueQTN$Map_info.pos >= as.numeric(vQTN.info$QTNwindow_lowerbound), ]
      DGLM.trueQTN <- DGLM.trueQTN[DGLM.trueQTN$Map_info.pos <= as.numeric(vQTN.info$QTNwindow_upperbound), ]

      DGLM.trueQTN$class <- 1
      id <- grep(TRUE, DGLM.results$Map_info.snp %in% DGLM.trueQTN$Map_info.snp)
      DGLM.results[c(id), ] <- DGLM.trueQTN
    
      MLM.trueQTN <- MLM.results[MLM.results$Chromosome == vQTN.info$chr, ]
      MLM.trueQTN <- MLM.trueQTN[MLM.trueQTN$Position >= as.numeric(vQTN.info$QTNwindow_lowerbound), ]
      MLM.trueQTN <- MLM.trueQTN[MLM.trueQTN$Position <= as.numeric(vQTN.info$QTNwindow_upperbound), ]
    
      MLM.trueQTN$class <- 1
      id <- grep(TRUE, MLM.results$SNP %in% MLM.trueQTN$SNP)
      MLM.results[c(id), ] <- MLM.trueQTN
    
    
      BFT.glm.fit <- glm(BFT.results$class ~ BFT.results$p.value, family = binomial)
      BFT.ROC <- roc(BFT.results$class, BFT.glm.fit$fitted.values, plot = T, legacy.axes = T, percent = T, xlab = "False Positive Percentage", ylab = "True Postive Percentage")
      AUC.mat[i, 5] <- BFT.ROC$auc
    
      DGLM.glm.fit <- glm(DGLM.results$class ~ DGLM.results$P.disp, family = binomial)
      DGLM.ROC <- roc(DGLM.results$class, DGLM.glm.fit$fitted.values, plot = T, legacy.axes = T, percent = T, xlab = "False Positive Percentage", ylab = "True Postive Percentage")
      AUC.mat[i, 6] <- DGLM.ROC$auc
    
      MLM.glm.fit <- glm(MLM.results$class ~ MLM.results$P.value, family = binomial)
      MLM.ROC <- roc(MLM.results$class, MLM.glm.fit$fitted.values, plot = T, legacy.axes = T, percent = T, xlab = "False Positive Percentage", ylab = "True Postive Percentage")
      AUC.mat[i, 7] <- MLM.ROC$auc
    
      tiff(file = paste0(randrep, "_", trait.name, vQTN.infoh[, 2], ".tiff"),
          width = 800, height = 400, res = 100)
      par(pty = "s")
      roc(BFT.results$class, BFT.glm.fit$fitted.values, plot = T, legacy.axes = T, percent = T, xlab = "False Positive %", ylab = "True Postive %", col = "darkorange3", lwd = 2, print.auc = T, print.auc.x = 45)
      plot.roc(DGLM.results$class, DGLM.glm.fit$fitted.values, percent = T, col = "navyblue", lwd = 2, print.auc = T, add = T, print.auc.y = 35)
      plot.roc(MLM.results$class, MLM.glm.fit$fitted.values, percent = T, col = "magenta4", lwd = 2, print.auc = T, add = T, print.auc.y = 25)
      dev.off()
    }
    write.csv(AUC.mat, paste0(trait.name, ".", "updated", ".", vQTN.infoh[, 2], "AUC.mat.test.csv"), row.names = F)
    
}
 # write.csv(AUC.mat, paste0(trait.name, ".", "updated", ".", "AUC.mat.test.csv"), row.names = F)
}

#Example with Molike MAF = 0.4 vQTN = 0.90, h2 = 0.63
##500
BFT.path.root.500 <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/vGWAS_results/BFT/500/"
DGLM.path.root.500 <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/vGWAS_results/DGLM/500/"
MLM.path.root.500 <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/vGWAS_results/GAPIT/500/"

##FULL
BFT.path.root.FULL <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/vGWAS_results/BFT/FULL/"
DGLM.path.root.FULL <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/vGWAS_results/DGLM/FULL/"
MLM.path.root.FULL <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/vGWAS_results/GAPIT/FULL/"

list.files(BFT.path.root.FULL)
list.files(DGLM.path.root.FULL)
list.files(MLM.path.root.FULL)

##Now, we do the same for Arabidopsis


####Maize
######Molike MAF40 h263 vQTN90
###500

#####Redoing MAF40 vQT90 h263 for GAPIT


#######Molike MAF10 h233 vQTN10
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.500, DGLM.path.root = DGLM.path.root.500, MLM.path.root = MLM.path.root.500, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF10.h233.vQTN10.500", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF10.h233.vQTN10.500", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF10.h233.vQTN10.500", vQTN.info = zm_vqgi_molike_MAF10_h233_vQTN10_500_250K, trait.name = "zm.molike.MAF10.h233.vQTN10.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.FULL, DGLM.path.root = DGLM.path.root.FULL, MLM.path.root = MLM.path.root.FULL, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF10.h233.vQTN10.FULL", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF10.h233.vQTN10.FULL", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF10.h233.vQTN10.FULL", vQTN.info = zm_vqgi_molike_MAF10_h233_vQTN10_full_250K, trait.name = "zm.molike.MAF10.h233.vQTN10.FULL")

#######Molike MAF10 h233 vQTN50
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.500, DGLM.path.root = DGLM.path.root.500, MLM.path.root = MLM.path.root.500, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF10.h233.vQTN50.500", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF10.h233.vQTN50.500", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF10.h233.vQTN50.500", vQTN.info = zm_vqgi_molike_MAF10_h233_vQTN50_500_250K, trait.name = "zm.molike.MAF10.h233.vQTN50.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.FULL, DGLM.path.root = DGLM.path.root.FULL, MLM.path.root = MLM.path.root.FULL, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF10.h233.vQTN50.FULL", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF10.h233.vQTN50_FULL", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF10.h233.vQTN50.FULL", vQTN.info = zm_vqgi_molike_MAF10_h233_vQTN50_full_250K, trait.name = "zm.molike.MAF10.h233.vQTN50.FULL")

######Molike MAF10 h233 vQTN90
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.500, DGLM.path.root = DGLM.path.root.500, MLM.path.root = MLM.path.root.500, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF10.h233.vQTN90.500", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF10.h233.vQTN90.500", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF10.h233.vQTN90.500", vQTN.info = zm_vqgi_molike_MAF10_h233_vQTN50_500_250K, trait.name = "zm.molike.MAF10.h233.vQTN90.500")


##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.FULL, DGLM.path.root = DGLM.path.root.FULL, MLM.path.root = MLM.path.root.FULL, BFT.fold.nam = "BFT.zm.vpheno.molike.MAF10.h233.vQTN90.FULL", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF10.h233.vQTN90_FULL", MLM.fold.nam = "MLM_zm.vpheno.molike.MAF10.h233.vQTN90.FULL", vQTN.info = zm_vqgi_molike_MAF10_h233_vQTN50_full_250K, trait.name = "zm.molike.MAF10.h233.vQTN90.FULL")


######Molike MAF10 h263 vQTN10
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.500, DGLM.path.root = DGLM.path.root.500, MLM.path.root = MLM.path.root.500, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF10.h263.vQTN10.500", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF10.h263.vQTN10.500", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF10.h263.vQTN10.500", vQTN.info = zm_vqgi_molike_MAF10_h263_vQTN10_500_250K, trait.name = "zm.molike.MAF10.h263.vQTN10.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.FULL, DGLM.path.root = DGLM.path.root.FULL, MLM.path.root = MLM.path.root.FULL, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF10.h263.vQTN10.FULL", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF10.h263.vQTN10_FULL", MLM.fold.nam = "zm.vpheno.molike.MAF10.h263.vQTN10.FULL", vQTN.info = zm_vqgi_molike_MAF10_h263_vQTN10_full_250K, trait.name = "zm.molike.MAF10.h263.vQTN10.FULL")

######Molike MAF10 h263 vQTN50
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.500, DGLM.path.root = DGLM.path.root.500, MLM.path.root = MLM.path.root.500, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF10.h263.vQTN50.500", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF10.h263.vQTN50.500", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF10.h263.vQTN50.500", vQTN.info = zm_vqgi_molike_MAF10_h263_vQTN50_500_250K, trait.name = "zm.molike.MAF10.h263.vQTN50.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.FULL, DGLM.path.root = DGLM.path.root.FULL, MLM.path.root = MLM.path.root.FULL, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF10.h263.vQTN50.FULL", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF10.h263.vQTN50.FULL", MLM.fold.nam = "GAPIT_MLM.zm.vpheno.molike.MAF10.h263.vQTN50.FULL", vQTN.info = zm_vqgi_molike_MAF10_h263_vQTN50_full_250K, trait.name = "zm.molike.MAF10.h263.vQTN50.FULL")


######Molike MAF10 h263 vQTN90 ###Come back later
#GAPIT_zm.vpheno.molike.MAF10.h263.vQTN90
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.500, DGLM.path.root = DGLM.path.root.500, MLM.path.root = MLM.path.root.500, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF10.h263.vQTN90.500", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF10.h263.vQTN90.500", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF10.h263.vQTN90.500", vQTN.info = zm_vqgi_molike_MAF10_h263_vQTN90_500_250K, trait.name = "zm.molike.MAF10.h263.vQTN90.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.FULL, DGLM.path.root = DGLM.path.root.FULL, MLM.path.root = MLM.path.root.FULL, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF10.h263.vQTN90.FULL", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF10.h263.vQTN90.FULL", MLM.fold.nam = "GAPIT_zm.vpheno.molike.MAF10.h263.vQTN90.FULL", vQTN.info = zm_vqgi_molike_MAF10_h263_vQTN50_full_250K, trait.name = "zm.molike.MAF10.h263.vQTN90.FULL")

######Molike MAF40 h233 vQTN10 
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.500, DGLM.path.root = DGLM.path.root.500, MLM.path.root = MLM.path.root.500, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF40.h233.vQTN10.500", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF40.h233.vQTN10.500", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF40.h233.vQTN10.500", vQTN.info = zm_vqgi_molike_MAF40_h233_vQTN10_500_250K, trait.name = "zm.molike.MAF40.h233.vQTN10.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.FULL, DGLM.path.root = DGLM.path.root.FULL, MLM.path.root = MLM.path.root.FULL, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF40.h233.vQTN10.FULL", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF40.h233.vQTN10.FULL", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF40.h233.vQTN10.FULL", vQTN.info = zm_vqgi_molike_MAF40_h233_vQTN10_full_250K, trait.name = "zm.molike.MAF40.h233.vQTN10.FULL")

######Molike MAF40 h233 vQTN50 
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.500, DGLM.path.root = DGLM.path.root.500, MLM.path.root = MLM.path.root.500, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF40.h233.vQTN50.500", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF40.h233.vQTN50.500", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF40.h233.vQTN50.500", vQTN.info = zm_vqgi_molike_MAF40_h233_vQTN50_500_250K, trait.name = "zm.molike.MAF40.h233.vQTN50.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.FULL, DGLM.path.root = DGLM.path.root.FULL, MLM.path.root = MLM.path.root.FULL, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF40.h233.vQTN50.FULL", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF40.h233.vQTN50.FULL", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF40.h233.vQTN50.FULL", vQTN.info = zm_vqgi_molike_MAF40_h233_vQTN50_full_250K, trait.name = "zm.molike.MAF40.h233.vQTN50.FULL")

######Molike MAF40 h233 vQTN90 
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.500, DGLM.path.root = DGLM.path.root.500, MLM.path.root = MLM.path.root.500, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF40.h233.vQTN90.500", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF40.h233.vQTN90.500", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF40.h233.vQTN90.500", vQTN.info = zm_vqgi_molike_MAF40_h233_vQTN90_500_250K, trait.name = "zm.molike.MAF40.h233.vQTN90.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.FULL, DGLM.path.root = DGLM.path.root.FULL, MLM.path.root = MLM.path.root.FULL, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF40.h233.vQTN90.FULL", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF40.h233.vQTN90.FULL", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF40.h233.vQTN90.FULL", vQTN.info = zm_vqgi_molike_MAF40_h233_vQTN90_full_250K, trait.name = "zm.molike.MAF40.h233.vQTN90.FULL")

######Molike MAF40 h263 vQTN10 ####Come back later
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.500, DGLM.path.root = DGLM.path.root.500, MLM.path.root = MLM.path.root.500, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF40.h263.vQTN10.500", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF40.h263.vQTN10.500", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF40.h263.vQTN10.500", vQTN.info = zm_vqgi_molike_MAF40_h263_vQTN10_500_250K, trait.name = "zm.molike.MAF40.h263.vQTN10.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.FULL, DGLM.path.root = DGLM.path.root.FULL, MLM.path.root = MLM.path.root.FULL, BFT.fold.nam = "zm.vpheno.molike.MAF40.h263.vQTN10.FULL", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF40.h263.vQTN10.FULL", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF40.h263.vQTN10.FULL", vQTN.info = zm_vqgi_molike_MAF40_h263_vQTN10_full_250K, trait.name = "zm.molike.MAF40.h263.vQTN10.FULL")

######Molike MAF40 h263 vQTN50 
#GAPIT_MLM.zm.vpheno.molike.MAF40.h263.vQTN50.FULL
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.500, DGLM.path.root = DGLM.path.root.500, MLM.path.root = MLM.path.root.500, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF40.h263.vQTN50.500", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF40.h263.vQTN50.500", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF40.h263.vQTN50.500", vQTN.info = zm_vqgi_molike_MAF40_h263_vQTN50_500_250K, trait.name = "zm.molike.MAF40.h263.vQTN50.500")


##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.FULL, DGLM.path.root = DGLM.path.root.FULL, MLM.path.root = MLM.path.root.FULL, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF40.h263.vQTN50.FULL", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF40.h263.vQTN50.FULL", MLM.fold.nam = "GAPIT_MLM.zm.vpheno.molike.MAF40.h263.vQTN50.FULL", vQTN.info = zm_vqgi_molike_MAF40_h263_vQTN50_full_250K, trait.name = "zm.molike.MAF40.h263.vQTN50.FULL")


######Molike MAF40 h263 vQTN90 
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.500, DGLM.path.root = DGLM.path.root.500, MLM.path.root = MLM.path.root.500, BFT.fold.nam = "BFT_zm.vPheno.molike.MAF40.500", DGLM.fold.nam = "DGLM_zm.vpheno.molike.MAF40.500", MLM.fold.nam = "MLM.zm.vpheno.molike.MAF40.h263.vQTN90.500", vQTN.info = zm_vqgi_molike_MAF40_h263_vQTN50_500_250K, trait.name = "zm.molike.MAF40.h263.vQTN90.500")


##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.FULL, DGLM.path.root = DGLM.path.root.FULL, MLM.path.root = MLM.path.root.FULL, BFT.fold.nam = "BFT_zm.vpheno.molike.MAF40.full", DGLM.fold.nam = "DGLM_zm.vPheno.molike.MAF40.full", MLM.fold.nam = "GAPIT.MLM.zm.vPheno.molike.MAF40.full", vQTN.info = zm_vqgi_molike_MAF40_h263_vQTN50_full_250K, trait.name = "zm.molike.MAF40.h263.vQTN90.FULL")



#####Epistasis lowh2
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.500, DGLM.path.root = DGLM.path.root.500, MLM.path.root = MLM.path.root.500, BFT.fold.nam = "BFT_zm.sp.epi.higheffectsize.h230.500", DGLM.fold.nam = "DGLM_zm.sp.epi.higheffectsize.h230.500", MLM.fold.nam = "MLM.zm.epistasis.higheffectsize.h230.500", vQTN.info = zm_eqgi_epistasis_lowh2_500_250K, trait.name = "zm.epistasis.lowh2.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.FULL, DGLM.path.root = DGLM.path.root.FULL, MLM.path.root = MLM.path.root.FULL, BFT.fold.nam = "BFT_zm.sp.epi.higheffectsize.h230.full", DGLM.fold.nam = "DGLM_zm.sp.epi.higheffectsize.h230.full", MLM.fold.nam = "GAPIT.MLM.zm.sp.epi.higheffectsize.h230.full", vQTN.info = zm_eqgi_epistasis_lowh2_FULL_250K, trait.name = "zm.epistasis.lowh2.FULL")

#####Epistasis highh2
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.500, DGLM.path.root = DGLM.path.root.500, MLM.path.root = MLM.path.root.500, BFT.fold.nam = "BFT_zm.sp.epi.higheffectsize.500", DGLM.fold.nam = "DGLM_zm.sp.epi.higheffectsize.500", MLM.fold.nam = "MLM.zm.sp.epi.higheffectsize.h280.500", vQTN.info = zm_eqgi_epistasis_highh2_500_250K, trait.name = "zm.epistasis.highh2.500")

###
spaceout <- c(83, 86, 88, 78, 23, 62, 61, 42, 33, 94)
zm_eqgi_epistasis_highh2_500_spaceout_250K <- zm_eqgi_epistasis_highh2_500_250K[5:6, ]

generate_ROC(repsel = TRUE, randvec = spaceout, BFT.path.root = BFT.path.root.500, DGLM.path.root = DGLM.path.root.500, MLM.path.root = MLM.path.root.500, BFT.fold.nam = "BFT_zm.sp.epi.higheffectsize.500", DGLM.fold.nam = "DGLM_zm.sp.epi.higheffectsize.500", MLM.fold.nam = "MLM.zm.sp.epi.higheffectsize.h280.500", vQTN.info = zm_eqgi_epistasis_highh2_500_spaceout_250K, trait.name = "zm.epistasis.highh2.500")


##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.FULL, DGLM.path.root = DGLM.path.root.FULL, MLM.path.root = MLM.path.root.FULL, BFT.fold.nam = "BFT_zm.sp.epi.higheffectsize.full.h280", DGLM.fold.nam = "DGLM_zm.sp.epi.h280.FULL", MLM.fold.nam = "GAPIT.MLM.zm.sp.epi.higheffectsize.full.h280", vQTN.info = zm_eqgi_epistasis_highh2_FULL_250K, trait.name = "zm.epistasis.highh2.FULL")

#####GxE
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.500, DGLM.path.root = DGLM.path.root.500, MLM.path.root = MLM.path.root.500, BFT.fold.nam = "BFT_sp.zm.GxE.500.mod", DGLM.fold.nam = "DGLM_sp.zm.GxE.500.mod", MLM.fold.nam = "MLM.zm.sp.GxE.500", vQTN.info = zm_gxeqgi_diff_500_250K, trait.name = "zm.GxE.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = BFT.path.root.FULL, DGLM.path.root = DGLM.path.root.FULL, MLM.path.root = MLM.path.root.FULL, BFT.fold.nam = "BFT_sp.zm.GxE.FULL.mod", DGLM.fold.nam = "DGLM_sp.zm.GxE.FULL.mod", MLM.fold.nam = "MLM.zm.GxE.FULL", vQTN.info = zm_gxeqgi_diff_FULL_250K, trait.name = "zm.GxE.FULL")


#####################
#This next section is for Arabidopsis

  ##FULL
at.BFT.path.root <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_at/vGWAS_results/BF/Revisions/"
at.DGLM.path.root <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_at/vGWAS_results/DGLM/Revisions/"
at.MLM.path.root <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_at/vGWAS_results/GAPIT/At/"

list.files(at.BFT.path.root)
list.files(at.DGLM.path.root)
list.files(at.MLM.path.root)

#######Molike MAF10 h233 vQTN10
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF10_h233_vQTN10/500", DGLM.fold.nam = "Molike_MAF10_h233_vQTN10/500", MLM.fold.nam = "Molike_MAF10_h233_vQTN10/500", vQTN.info = at_vqgi_molike_wmQTL_MAF10_h233_vQTN10_500, trait.name = "at.molike.MAF10.h233.vQTN10.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF10_h233_vQTN10/FULL", DGLM.fold.nam = "Molike_MAF10_h233_vQTN10/FULL", MLM.fold.nam = "Molike_MAF10_h233_vQTN10/FULL", vQTN.info = at_vqgi_molike_wmQTL_MAF10_h233_vQTN10_full, trait.name = "at.molike.MAF10.h233.vQTN10.FULL")

#######Molike MAF10 h233 vQTN50
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF10_h233_vQTN50/500", DGLM.fold.nam = "Molike_MAF10_h233_vQTN50/500", MLM.fold.nam = "Molike_MAF10_h233_vQTN50/500", vQTN.info = at_vqgi_molike_wmQTL_MAF10_h233_vQTN50_500, trait.name = "at.molike.MAF10.h233.vQTN50.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF10_h233_vQTN50/FULL", DGLM.fold.nam = "Molike_MAF10_h233_vQTN50/FULL", MLM.fold.nam = "Molike_MAF10_h233_vQTN50/FULL", vQTN.info = at_vqgi_molike_wmQTL_MAF10_h233_vQTN50_full, trait.name = "at.molike.MAF10.h233.vQTN50.FULL")

#######Molike MAF10 h233 vQTN90
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF10_h233_vQTN90/500", DGLM.fold.nam = "Molike_MAF10_h233_vQTN90/500", MLM.fold.nam = "Molike_MAF10_h233_vQTN90/500", vQTN.info = at_vqgi_molike_wmQTL_MAF10_h233_vQTN90_500, trait.name = "at.molike.MAF10.h233.vQTN90.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF10_h233_vQTN90/FULL", DGLM.fold.nam = "Molike_MAF10_h233_vQTN90/FULL", MLM.fold.nam = "Molike_MAF10_h233_vQTN90/FULL", vQTN.info = at_vqgi_molike_wmQTL_MAF10_h233_vQTN90_full, trait.name = "at.molike.MAF10.h233.vQTN90.FULL")

#######Molike MAF10 h263 vQTN10
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF10_h263_vQTN10/500", DGLM.fold.nam = "Molike_MAF10_h263_vQTN10/500", MLM.fold.nam = "Molike_MAF10_h263_vQTN10/500", vQTN.info = at_vqgi_molike_wmQTL_MAF10_h263_vQTN10_500, trait.name = "at.molike.MAF10.h263.vQTN10.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF10_h263_vQTN10/FULL", DGLM.fold.nam = "Molike_MAF10_h263_vQTN10/FULL", MLM.fold.nam = "Molike_MAF10_h263_vQTN10/FULL", vQTN.info = at_vqgi_molike_wmQTL_MAF10_h263_vQTN10_full, trait.name = "at.molike.MAF10.h263.vQTN10.FULL")

#######Molike MAF10 h263 vQTN50
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF10_h263_vQTN50/500", DGLM.fold.nam = "Molike_MAF10_h263_vQTN50/500", MLM.fold.nam = "Molike_MAF10_h263_vQTN50/500", vQTN.info = at_vqgi_molike_wmQTL_MAF10_h263_vQTN50_500, trait.name = "at.molike.MAF10.h263.vQTN50.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF10_h263_vQTN50/FULL", DGLM.fold.nam = "Molike_MAF10_h263_vQTN50/FULL", MLM.fold.nam = "Molike_MAF10_h263_vQTN50/FULL", vQTN.info = at_vqgi_molike_wmQTL_MAF10_h263_vQTN50_full, trait.name = "at.molike.MAF10.h263.vQTN50.FULL")

#######Molike MAF10 h263 vQTN90
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF10_h263_vQTN90/500", DGLM.fold.nam = "Molike_MAF10_h263_vQTN90/500", MLM.fold.nam = "Molike_MAF10_h263_vQTN90/500", vQTN.info = at_vqgi_molike_wmQTL_MAF10_h263_vQTN90_500, trait.name = "at.molike.MAF10.h263.vQTN90.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF10_h263_vQTN90/FULL", DGLM.fold.nam = "Molike_MAF10_h263_vQTN90/FULL", MLM.fold.nam = "Molike_MAF10_h263_vQTN90/FULL", vQTN.info = at_vqgi_molike_wmQTL_MAF10_h263_vQTN90_full, trait.name = "at.molike.MAF10.h263.vQTN90.FULL")

#######Molike MAF40 h233 vQTN10
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF40_h233_vQTN10/500", DGLM.fold.nam = "Molike_MAF40_h233_vQTN10/500", MLM.fold.nam = "Molike_MAF40_h233_vQTN10/500", vQTN.info = at_vqgi_molike_wmQTL_MAF40_h233_vQTN10_500, trait.name = "at.molike.MAF40.h233.vQTN10.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF40_h233_vQTN10/FULL", DGLM.fold.nam = "Molike_MAF40_h233_vQTN10/FULL", MLM.fold.nam = "Molike_MAF40_h233_vQTN10/FULL", vQTN.info = at_vqgi_molike_wmQTL_MAF40_h233_vQTN10_full, trait.name = "at.molike.MAF40.h233.vQTN10.FULL")

#######Molike MAF40 h233 vQTN50
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF40_h233_vQTN50/500", DGLM.fold.nam = "Molike_MAF40_h233_vQTN50/500", MLM.fold.nam = "Molike_MAF40_h233_vQTN50/500", vQTN.info = at_vqgi_molike_wmQTL_MAF40_h233_vQTN50_500, trait.name = "at.molike.MAF40.h233.vQTN50.500")

##FULL
gene#######Molike MAF40 h233 vQTN90
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF40_h233_vQTN90/500", DGLM.fold.nam = "Molike_MAF40_h233_vQTN90/500", MLM.fold.nam = "Molike_MAF40_h233_vQTN90/500", vQTN.info = at_vqgi_molike_wmQTL_MAF40_h233_vQTN90_500, trait.name = "at.molike.MAF40.h233.vQTN90.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF40_h233_vQTN90/FULL", DGLM.fold.nam = "Molike_MAF40_h233_vQTN90/FULL", MLM.fold.nam = "Molike_MAF40_h233_vQTN90/FULL", vQTN.info = at_vqgi_molike_wmQTL_MAF40_h233_vQTN90_full, trait.name = "at.molike.MAF40.h233.vQTN90.FULL")
rate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF40_h233_vQTN50/FULL", DGLM.fold.nam = "Molike_MAF40_h233_vQTN50/FULL", MLM.fold.nam = "Molike_MAF40_h233_vQTN50/FULL", vQTN.info = at_vqgi_molike_wmQTL_MAF40_h233_vQTN50_full, trait.name = "at.molike.MAF40.h233.vQTN50.FULL")

#######Molike MAF40 h263 vQTN10
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF40_h263_vQTN10/500", DGLM.fold.nam = "Molike_MAF40_h263_vQTN10/500", MLM.fold.nam = "Molike_MAF40_h263_vQTN10/500", vQTN.info = at_vqgi_molike_wmQTL_MAF40_h263_vQTN10_500, trait.name = "at.molike.MAF40.h263.vQTN10.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF40_h263_vQTN10/FULL", DGLM.fold.nam = "Molike_MAF40_h263_vQTN10/FULL", MLM.fold.nam = "Molike_MAF40_h263_vQTN10/FULL", vQTN.info = at_vqgi_molike_wmQTL_MAF40_h263_vQTN10_full, trait.name = "at.molike.MAF40.h263.vQTN10.FULL")

#######Molike MAF40 h263 vQTN50
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF40_h263_vQTN50/500", DGLM.fold.nam = "Molike_MAF40_h263_vQTN50/500", MLM.fold.nam = "Molike_MAF40_h263_vQTN50/500", vQTN.info = at_vqgi_molike_wmQTL_MAF40_h263_vQTN50_500, trait.name = "at.molike.MAF40.h263.vQTN50.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF40_h263_vQTN50/FULL", DGLM.fold.nam = "Molike_MAF40_h263_vQTN50/FULL", MLM.fold.nam = "Molike_MAF40_h263_vQTN50/FULL", vQTN.info = at_vqgi_molike_wmQTL_MAF40_h263_vQTN50_full, trait.name = "at.molike.MAF40.h263.vQTN50.FULL")

######Molike MAF40 h263 vQTN90
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF40/Sub500", DGLM.fold.nam = "molike_wmQTL/500", MLM.fold.nam = "Molike_MAF40/500", vQTN.info = at_vqgi_molike_wmQTL_MAF40_h263_vQTN50_500, trait.name = "at.molike.MAF40.h263.vQTN90.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Molike_MAF40/FULL", DGLM.fold.nam = "molike_wmQTL/FULL", MLM.fold.nam = "Molike_MAF40/FULL", vQTN.info = at_vqgi_molike_wmQTL_MAF40_h263_vQTN50_full, trait.name = "at.molike.MAF40.h263.vQTN90.FULL")

#####Epistasis lowh2
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Epistasis_lowh2_higheffectsizes/500", DGLM.fold.nam = "Epistasis_lowh2_higheffectsizes/500", MLM.fold.nam = "Epistasis_higheffectsize_h230/500", vQTN.info = at_eqgi_epistasis_lowh2_highefs_500, trait.name = "at.epistasis.lowh2.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Epistasis_lowh2_higheffectsizes/FULL", DGLM.fold.nam = "Epistasis_lowh2_higheffectsizes/FULL", MLM.fold.nam = "Epistasis_higheffectsize_h230/FULL", vQTN.info = at_eqgi_epistasis_lowh2_highefs_FULL, trait.name = "at.epistasis.lowh2.FULL")

#####Epistasis highh2
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Epistasis_higheffectsize/500", DGLM.fold.nam = "Epistasis_higheffectsize/500", MLM.fold.nam = "Epistasis_higheffectsize/500", vQTN.info = at_eqgi_epistasis_lowh2_highefs_500, trait.name = "at.epistasis.highh2.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "Epistasis_higheffectsize/FULL", DGLM.fold.nam = "Epistasis_higheffectsize/FULL", MLM.fold.nam = "Epistasis_higheffectsize/FULL", vQTN.info = at_eqgi_epistasis_lowh2_highefs_FULL, trait.name = "at.epistasis.highh2.FULL")

#####GxE
##500
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "GxE_diff/500", DGLM.fold.nam = "GxE_diff/500", MLM.fold.nam = "GxE_diff/500", vQTN.info = at_gxeqgi_diff_500_100K, trait.name = "at.GxE.500")

##FULL
generate_ROC(repsel = NULL, randvec = NULL, BFT.path.root = at.BFT.path.root, DGLM.path.root = at.DGLM.path.root, MLM.path.root = at.MLM.path.root, BFT.fold.nam = "GxE_diff/FULL", DGLM.fold.nam = "GxE_diff/FULL", MLM.fold.nam = "GxE_diff/FULL", vQTN.info = at_gxeqgi_diff_FULL_100K, trait.name = "at.GxE.FULL")



##################################
save.image("~/vGWAS_chaptI_RII.RData")
