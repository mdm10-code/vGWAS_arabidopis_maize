#########################Brown-Forsythe Script######################

#First, load required libraries (install them if not done so already)
install.packages("car")

library(vGWAS)
library(car)
library(plyr)
library(readr)

###This is the code to run your Brown-Forysthe Test
###Because Brown-Forsythe is more robust to deviations from 
###normality compared to Levene's Test

###Here, I am going to try and make BF to efficiently run over all columns in simulated vQTL phenotypes
BFT.deluxe <- function(Geno=NULL, vPheno=NULL, Map=NULL, trait.name=NULL) {
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
    write.csv(BF_Results, file = paste("BF_", trait.name, "_", colnames(vPheno[i]), ".csv", sep = ""), sep = ",", row.names = FALSE)
  }
}

###For the next few lines, I will be running BFT.Deluxe for all of my variance traits.
###Preferably, you want to visualize the file hierarchy for all of your results
###To efficiently store my results, I will create result folders to store
###My Brown-Forsythe Test and DGLM results
###These next lines are for the Brown-Forsythe Test
dir.create(paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/")) #general results folder
dir.create(paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/")) #general results folder

#null results
dir.create(paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/null/")) #null results
dir.create(paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/null/500/")) #500 diversity panel
dir.create(paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/null/full/")) #full diversity panel

#Molike MAF40 wMQTL
dir.create(paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/molike_wmQTL/")) 
dir.create(paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/molike_wmQTL/500/")) #500 diversity panel
dir.create(paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/molike_wmQTL/full/")) #full diversity panel

#Molike MAF40 womQTL
dir.create(paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/molike_womQTL/")) 
dir.create(paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/molike_womQTL/500/")) #500 diversity panel
dir.create(paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/molike_womQTL/full/")) #full diversity panel

#Oillike MAF10 womQTL
dir.create(paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/oillike_womQTL_MAF10/")) 
dir.create(paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/oillike_womQTL_MAF10/500/")) #500 diversity panel
dir.create(paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/oillike_womQTL_MAF10/full/")) #full diversity panel

#Oillike MAF40 womQTL
dir.create(paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/oillike_womQTL_MAF40/")) 
dir.create(paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/oillike_womQTL_MAF40/500/")) #500 diversity panel
dir.create(paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/oillike_womQTL_MAF40/full/")) #full diversity panel

##For maize
#null results
dir.create(paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/null/")) #null results
dir.create(paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/null/500/")) #500 diversity panel
dir.create(paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/null/full/")) #full diversity panel

#Molike MAF40 wMQTL
dir.create(paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/molike_wmQTL/")) 
dir.create(paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/molike_wmQTL/500/")) #500 diversity panel
dir.create(paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/molike_wmQTL/full/")) #full diversity panel

#Molike MAF40 womQTL
dir.create(paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/molike_womQTL/")) 
dir.create(paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/molike_womQTL/500/")) #500 diversity panel
dir.create(paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/molike_womQTL/full/")) #full diversity panel

#Oillike MAF10 womQTL
dir.create(paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/oillike_womQTL_MAF10/")) 
dir.create(paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/oillike_womQTL_MAF10/500/")) #500 diversity panel
dir.create(paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/oillike_womQTL_MAF10/full/")) #full diversity panel

#Oillike MAF40 womQTL
dir.create(paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/oillike_womQTL_MAF40/")) 
dir.create(paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/oillike_womQTL_MAF40/500/")) #500 diversity panel
dir.create(paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/oillike_womQTL_MAF40/full/")) #full diversity panel

##Now, I can assign the main vGWAS separate wds
wd$at_vGWAS_results <- paste0(wd$vGWAS_chapter_one_at, "vGWAS_results/") #general results folder
wd$zm_vGWAS_results <- paste0(wd$vGWAS_chapter_one_zm, "vGWAS_results/") #general results folder

####Now, I can run the Brown-Forsythes test from the vGWAS R package. I will need the pure geno 
#objects for this to work

#For Arabidopsis
##NULL
###500 diversity panel
setwd(paste0(wd$at_vGWAS_results, "null/500/"))
BFT.deluxe(Geno = at_just_geno_500, vPheno = at_sp_null_500, Map = at_map[, c(1, 3, 4)], trait.name = "at_null_500")

###FULL 
setwd(paste0(wd$at_vGWAS_results, "null/full/"))
BFT.deluxe(Geno = at_just_geno_full, vPheno = at_sp_null_full, Map = at_map[, c(1, 3, 4)], trait.name = "at_null_full")

##Molike with mQTL MAF40
setwd(paste0(wd$at_vGWAS_results, "molike_wmQTL/500/"))
BFT.deluxe(Geno = at_just_geno_500, vPheno = at_vpheno_molike_MAF40_500, Map = at_map[, c(1, 3, 4)], trait.name = "at_molike_wmQTL_500")

###FULL 
setwd(paste0(wd$at_vGWAS_results, "molike_wmQTL/full/"))
BFT.deluxe(Geno = at_just_geno_full, vPheno = at_vpheno_molike_MAF40_full, Map = at_map[, c(1, 3, 4)], trait.name = "at_molike_wmQTL_full")

##Molike womQTL MAF40
###500
setwd(paste0(wd$at_vGWAS_results, "molike_womQTL/500/"))
BFT.deluxe(Geno = at_just_geno_500, vPheno = at_vpheno_molike_womQTL_MAF40_500, Map = at_map[, c(1, 3, 4)], trait.name = "At_Molike_womQTL_MAF40_500")

###FULL
setwd(paste0(wd$at_vGWAS_results, "molike_womQTL/full/"))
BFT.deluxe(Geno = at_just_geno_full, vPheno = at_vpheno_molike_womQTL_MAF40_full, Map = at_map[, c(1, 3, 4)], trait.name = "At_Molike_womQTL_MAF40_FULL")

##Oillike_womQTL_MAF10
###500
setwd(paste0(wd$at_vGWAS_results, "oillike_womQTL_MAF10/500/"))
BFT.deluxe(Geno = at_just_geno_500, vPheno = at_vpheno_oil_500_womQTL_MAF10, Map = at_map[, c(1, 3, 4)], trait.name = "at_oillike_womQTL_MAF10_500")

###FULL
setwd(paste0(wd$at_vGWAS_results, "oillike_womQTL_MAF10/full/"))
BFT.deluxe(Geno = at_just_geno_full, vPheno = at_vpheno_oil_full_womQTL_MAF10, Map = at_map[, c(1, 3, 4)], trait.name = "at_oillike_womQTL_MAF10_full")

##Oillike_womQTL_MAF40
###500
setwd(paste0(wd$at_vGWAS_results, "oillike_womQTL_MAF40/500/"))
BFT.deluxe(Geno = at_just_geno_500, vPheno = at_vpheno_oil_500_womQTL_MAF40, Map = at_map[, c(1, 3, 4)], trait.name = "at_oillike_womQTL_MAF40_500")

###FULL
setwd(paste0(wd$at_vGWAS_results, "oillike_womQTL_MAF10/full/"))
BFT.deluxe(Geno = at_just_geno_full, vPheno = at_vpheno_oil_full_womQTL_MAF40, Map = at_map[, c(1, 3, 4)], trait.name = "at_oillike_womQTL_MAF40_full")

#For Maize
##NULL
###500 diversity panel
setwd(paste0(wd$zm_vGWAS_results, "null/500/"))
BFT.deluxe(Geno = zm_just_geno_500, vPheno = sp_zm_null_500, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_null_500")

###FULL 
setwd(paste0(wd$zm_vGWAS_results, "null/full/"))
BFT.deluxe(Geno = zm_just_geno_full, vPheno = sp_zm_null_full, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_null_full")

##Molike with mQTL MAF40
###500 diversity panel
setwd(paste0(wd$zm_vGWAS_results, "molike_wmQTL/500/"))
BFT.deluxe(Geno = zm_just_geno_500, vPheno = zm_vpheno_molike_MAF40_500, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_molike_wmQTL_500")

###FULL 
setwd(paste0(wd$zm_vGWAS_results, "molike_wmQTL/full/"))
BFT.deluxe(Geno = zm_just_geno_full, vPheno = zm_vPheno_molike_MAF40_full, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_molike_wmQTL_full")

##Molike womQTL MAF40
###500
setwd(paste0(wd$zm_vGWAS_results, "molike_womQTL/500/"))
BFT.deluxe(Geno = zm_just_geno_500, vPheno = zm_vpheno_molike_womQTL_MAF40_500, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_molike_womQTL_MAF40_500")

###FULL
setwd(paste0(wd$zm_vGWAS_results, "molike_womQTL/full/"))
BFT.deluxe(Geno = zm_just_geno_full, vPheno = zm_vpheno_molike_womQTL_MAF40_full, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_molike_womQTL_MAF40_full")

##Oillike_womQTL_MAF10
###500
setwd(paste0(wd$zm_vGWAS_results, "oillike_womQTL_MAF10/500/"))
BFT.deluxe(Geno = zm_just_geno_500, vPheno = zm_vpheno_oil_500_womQTL_MAF10, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_oillike_womQTL_MAF10_500")

###FULL
setwd(paste0(wd$zm_vGWAS_results, "oillike_womQTL_MAF10/full/"))
BFT.deluxe(Geno = zm_just_geno_full, vPheno = zm_vpheno_oil_full_womQTL_MAF10, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_oillike_womQTL_MAF10_full")

##Oillike_womQTL_MAF40
###500
setwd(paste0(wd$zm_vGWAS_results, "oillike_womQTL_MAF40/500/"))
BFT.deluxe(Geno = zm_just_geno_500, vPheno = zm_vpheno_oil_500_womQTL_MAF40, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_oillike_womQTL_MAF40_500")

###FULL
setwd(paste0(wd$zm_vGWAS_results, "oillike_womQTL_MAF10/full/"))
BFT.deluxe(Geno = zm_just_geno_full, vPheno = zm_vpheno_oil_full_womQTL_MAF40, Map = zm_map[, c(1, 3, 4)], trait.name = "zm_oillike_womQTL_MAF40_full")