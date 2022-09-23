############################################
### Brown-Forsythe Test Express Script ###
############################################

load("vGWAS_chaptI_RII.RData")

library(vGWAS)
library(dplyr)

zm.pheno.FULL.listI <- list(sp.zm.null.full,
                            sp.zm.GxE.FULL.mod,
                            zm.vpheno.molike.MAF10.h233.vQTN10.FULL,
                            zm.vpheno.molike.MAF10.h233.vQTN50.FULL,
                            zm.vpheno.molike.MAF10.h233.vQTN90.FULL,
                            zm.vpheno.molike.MAF10.h263.vQTN10.FULL)

names(zm.pheno.FULL.listI)[[1]] <- "sp.zm.null.full"
names(zm.pheno.FULL.listI)[[2]] <- "sp.zm.GxE.FULL.mod"
names(zm.pheno.FULL.listI)[[3]] <- "zm.vpheno.molike.MAF10.h233.vQTN10.FULL"
names(zm.pheno.FULL.listI)[[4]] <- "zm.vpheno.molike.MAF10.h233.vQTN50.FULL"
names(zm.pheno.FULL.listI)[[5]] <- "zm.vpheno.molike.MAF10.h233.vQTN90.FULL"
names(zm.pheno.FULL.listI)[[6]] <- "zm.vpheno.molike.MAF10.h263.vQTN10.FULL"

#This is for list II
zm.pheno.FULL.listII <- list(zm.vpheno.molike.MAF10.h263.vQTN50.FULL,
                             zm.vpheno.molike.MAF10.h263.vQTN90.FULL,
                             zm.vpheno.molike.MAF40.h233.vQTN10.FULL,
                             zm.vpheno.molike.MAF40.h233.vQTN50.FULL,
                             zm.vpheno.molike.MAF40.h233.vQTN90.FULL,
                             zm.vpheno.molike.MAF40.h263.vQTN10.FULL)

names(zm.pheno.FULL.listII)[[1]] <- "zm.vpheno.molike.MAF10.h263.vQTN50.FULL"
names(zm.pheno.FULL.listII)[[2]] <- "zm.vpheno.molike.MAF10.h263.vQTN90.FULL"
names(zm.pheno.FULL.listII)[[3]] <- "zm.vpheno.molike.MAF40.h233.vQTN10.FULL"
names(zm.pheno.FULL.listII)[[4]] <- "zm.vpheno.molike.MAF40.h233.vQTN50.FULL"
names(zm.pheno.FULL.listII)[[5]] <- "zm.vpheno.molike.MAF40.h233.vQTN90.FULL"
names(zm.pheno.FULL.listII)[[6]] <- "zm.vpheno.molike.MAF40.h263.vQTN10.FULL"

setwd("C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/RII/2_pipeline/out/vGWAS_chapter_one_zm/BFT/FULL/")
wd$zm_bft_FULL <- "C:/Users/mdm10/Documents/Projects/Dissertation/Chapter_one/RII/2_pipeline/out/vGWAS_chapter_one_zm/BFT/FULL/"

BFT.deluxe2 <- function(wd = NULL, trait_list = NULL, Geno = NULL, Map = NULL) {
  for(h in 1:length(trait_list)) {
    print(paste0("this is trait", h))
    new_folder <- dir.create(paste0(wd, names(trait_list)[[h]]))
    setwd(paste0(wd, names(trait_list)[[h]]))
    vPheno <- trait_list[[h]]
    for (i in 2:ncol(vPheno)) { #This outer for loop is for the phenotypic object
      print(paste("--------- This is Phenotype ", i - 1, " out of ", dim(vPheno)[2] - 1, "  ----------", sep = ""))
      BF_Results <- matrix(NA, nrow = dim(Geno)[2], ncol = 1) #Create a column for the p-values
      colnames(BF_Results) <- c("p.value")
      BF_Results <- data.frame(Map, BF_Results) #Merge Map info (three columns) with BF Results
      for (j in 1:ncol(Geno)) { #This inner for loop iterates through the genotypic object
        vGWA1 <- brown.forsythe.test(vPheno[, i], as.factor(Geno[, j]), kruskal.test = FALSE)
        BF_Results[j, 4] <- as.numeric(vGWA1[2]) #populate the p-value column
        if(j %% 100 == 0) #This checks progress of BFT
          print(j)
      }
      BF_Results$fdr <- p.adjust(BF_Results$p.value, "fdr") #fdr adjustment
      BF_Results <- BF_Results[order(BF_Results$fdr), ] #order by fdr
      write.csv(BF_Results, file = paste("BF_", names(trait_list)[[h]], "_", colnames(vPheno[i]), ".csv", sep = ""), sep = ",", row.names = FALSE)
    }
  }
}

##Test run
BFT.deluxe2(wd = wd$zm_bft_FULL, trait_list = zm_bft_FULL_list, Geno = t(zm_just_geno_FULL_filtlike_at), Map = zm_map_filtlike_at[, c(1, 3, 4)])

##
BFT.deluxe2(wd = wd$zm_bft_FULL, trait_list = zm.pheno.FULL.listI, Geno = t(zm.just.geno.full), Map = zm.map[, c(1, 3, 4)])

##Trait list II
BFT.deluxe2(wd = wd$zm_bft_FULL, trait_list = zm.pheno.FULL.listII, Geno = zm.just.geno.full, Map = zm.map[, c(1, 3, 4)])

##Save workspaces
save.image("vGWAS_chaptI_RII.RData")
save.image("vGWAS_chaptI_RII_copy.RData")
