
library(data.table)
library(SNPRelate)
gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
library(tidyverse)
library(SNPRelate)
#---------- data summary -----------
library(here)
zm_bft_full_results <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/vGWAS_results/BFT/FULL/"
gds.folder <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/"
qtn.path <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/zm_vQTL_simulations/"

LD_window_tp_fp_rate_calc <- function(LD = NULL, folder = NULL, a = NULL, species = NULL, size = NULL, gds.folder = NULL, gds.file = NULL, qtn.path = NULL, qtn.fol.ind = NULL) {
  library(here)
  ld_threshold <- LD
  print(ld_threshold)
  resultsfolder <- dir(folder) #This is the specific results folder
  resultsfolder <- resultsfolder[a]
  print(resultsfolder)
  for (z in seq_along(folder)) {
    print(species)
    print(size)
    gds <- paste0(gds.folder, gds.file)
    genofile <- snpgdsOpen(gds)
    sims <- dir(qtn.path)[6]
    print(sims)
    qtn <- fread(
      paste0(qtn.path, sims, "/",  
             "var.QTN.genotypic.information.csv"), #simulated QTNs
      header = TRUE,
      sep = ",", #Change to reflect 
      data.table = F
    )[,1:5]
    print(qtn)
    BFT_results <- list.files(paste0(folder, resultsfolder, "/"))
    n_gwas <- length(BFT_results)
    gwas_BFT <- vector("list", n_gwas)
    detection <- tibble(
      species = species,
      size = size,
      BFT = logical(n_gwas))
    #print(detection)
    for (n in 1:n_gwas){
      #    rep <- gsub("X.*", n, BFT_results[n]) #?
      BFT <- fread( #change this to BFT
        paste0(folder, resultsfolder, "/", BFT_results[n]),
        header = TRUE,
        sep = ",",
        select = c("snp", "chr", "pos", "p.value"), # this needs to be changed to GAPIT variables
        data.table = F
      )
      BFT <- BFT %>% filter(!snp %in% qtn$snp) %>%
        arrange(p.value) %>% mutate(fdr =  p.adjust(p.value, "fdr")) %>%
        filter(fdr <= 0.05)
    gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
    }
  }
}

path = "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/R2/2_pipeline/out/vGWAS_chapter_one_zm/" #This is the path to the results folder
setwd(path) 
library(here)
folder <- dir(paste0(path, "/vGWAS_results/BFT/FULL")) #This is the specific results folder
ld_threshold <- 0.01 #LD threshold
pv <- 0.000000001
results <- list()
for (z in seq_along(folder)){
  species <- ifelse(grepl("zm", folder[z]), "maize", "ara")
  size <- ifelse(grepl("FULL",folder[z]), "FULL", "500")
  if (species == "maize") {
    gds <- ifelse(grepl("FULL",folder[z]),
                  "./zm.geno.numeric.full.gds",
                  "./data/maize/ames_2815_mac5_ld09.gds")
  } else {
    gds <- ifelse(grepl("500",folder[z]),
                  "./data/soybean/soy_500_mac5_ld09.gds",
                  "./data/soybean/soy_2815_mac5_ld09.gds")
  }
  genofile <- snpgdsOpen(gds)
  qtn <- fread(
    paste0(path, "zm_vQTL_simulations/",
           "zm_vmolike_MAF10_h233_vQTN90_FULL/",
           "var.QTN.genotypic.information.csv"), #simulated QTNs
    header = TRUE,
    sep = ",", #Change to reflect 
    data.table = F
  )[,1:5]
  
  BFT_results <- dir(paste0(path, "/vGWAS_results/BFT/FULL/", folder[z]))
#  GAPIT_results <- GAPIT_results[!grepl("log.txt",GAPIT_results)]
  
  n_gwas <- length(BFT_results)
  gwas_BFT <- vector("list", n_gwas)
  detection <- tibble(
    species = species,
    size = size,
    BFT = logical(n_gwas))
  
  for (n in 1:n_gwas){
#    rep <- gsub("X.*", n, BFT_results[n]) #?
    BFT <- fread( #change this to BFT
      paste0(path, "/vGWAS_results/BFT/FULL/", folder[z], "/", BFT_results[n]),
      header = TRUE,
      sep = ",",
      select = c("snp", "chr", "pos", "p.value"), # this needs to be changed to GAPIT variables
      data.table = F
    )
   BFT <- BFT %>% filter(!snp %in% qtn$snp) %>%
      arrange(p.value) %>% mutate(fdr =  p.adjust(p.value, "fdr")) %>%
     filter(fdr <= 0.05)
      snps <- read.gdsn(index.gdsn(genofile, "snp.id"))

      #true detection
      detected_BFT <- c()
      l <- 1
      for( i in qtn$snp) {
        #GAPIT
        if (i %in% BFT$snp) {
          detected_BFT[l] <- TRUE
        } else { ###I am at this point 3-3-2022
          set <-
            snps[which(snps %in% c(i, BFT$snp))] 

          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD^2
          colnames(ld) = rownames(ld) = snps[snps %in% c(i, BFT$snp)]
          detected_BFT[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        l <- l + 1
      }
      detection$BFT[n] <- all(detected_BFT)
    
    print(n)
  }
  results[[z]] <- detection
  gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
  cat("-----------done with ---------", z)
}
# fwrite(total_all, "total_all_count_p_fdr05.txt", sep = "\t", quote = F, row.names = F)
