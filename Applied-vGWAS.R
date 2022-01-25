######################################################
############## Applied vGWAS #########################
######################################################

#First, I will load the SNP55K array 
snp_55K_maize_filtered <- read.delim("~/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/vGWAS_chapter_one_zm/zm_sp_results/applied/out_geno_numeric.txt")

#Next, I will do some data manipulations
bg282_map <- snp_55K_maize_filtered[, c(1, 3, 4)]
bg282_geno <- snp_55K_maize_filtered[, -c(1:5)]
bg282_geno <- t(bg282_geno)

#Then, I will need to create a matrices of ones to get SNPs encoded as 0, 1, 2
just_ones_4_282 <- matrix(1, nrow = 279, ncol = 48880)

bg282_geno <- bg282_geno + just_ones_4_282

#We have to rename four accessions in new_df_4
new_df_4_tassel[1, 1] <- "X4226"
new_df_4_tassel[2, 1] <- "X4722"
new_df_4_tassel[3, 1] <- "X3316"
new_df_4_tassel[4, 1] <- "X3811"

#I will need to rearrange the genotypic data frame 
new_df_4_tassel <- new_df_4_tassel[c(3, 4, 1, 2, 5:282), ]

#I will use my BFT deluxe function for BFT
BFT.deluxe(Geno = bg282_geno, vPheno = zm_ph_GWAS_Residuals.for.PH._residuals[c(genphenmat), ], Map = bg282_map, trait.name = "BG282_")

##This section is to define DGLM function for Applied vGWAS
  my.pdglm <- function(cT = NULL, i = NULL, Phenos = NULL, Geno = NULL, Map = NULL) {
    print(paste("--------- Fitting DGLM model for SNP ", i, " out of ", dim(Geno)[2], "  ----------", sep = ""))
    temp_123456789 <<- data.frame(y = Phenos[, cT], snp = Geno[, i])
    #you dont need to change the formula anymore, it includes however number of PCs you have
    model <-
      dglm(
        formula = as.formula(paste("y ~ snp")),
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

#I need to match the dimensions of my phenotypic and genotypic objects
genphenmatch <- grep(TRUE, rownames(bg282_geno) %in% new_df_4_tassel$`<Trait>`)

#This is to create the DGLM result data file
TF <- matrix(NA, nrow = dim(bg282_geno)[2], ncol = 7)
row.names(TF) <- colnames(bg282_geno)
colnames(TF) <- c("Map_info.snp", "Map_info.chr", "Map_info.pos", "Beta", "SE", "P.mean", "P.disp")

##I will iterate through the genotypic objects using a for loop
for (i in 1:dim(bg282_geno)[2]) {
  try({
    outm <- my.pdglm(cT = 2, i = i, Phenos = zm_ph_GWAS_Residuals.for.PH._residuals[c(genphenmat), ], Geno = bg282_geno, Map = bg282_map)
    TF[i, ] <- as.numeric(outm)
    print(i)
  }, silent = TRUE)
}

# Now add map information and select required columns for Manhattan plot
row.names(TF) <- colnames(bg282_geno)
colnames(TF) <- c("Map_info.snp", "Map_info.chr", "Map_info.pos", "Beta", "SE", "P.mean", "P.disp")
TF2 <- merge(map, TF)
TF2 <- data.frame(merge(map, TF[, 3:4], by = "row.names", all.x = TRUE))  # add map info
colnames(TF2) <- c("marker", "chrom", "pos", "p.mean", "p.disp")

#I will define the path to where the out directory is located
wd$out <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/2_pipeline/out/"
setwd(wd$out)

write.csv(TF, "Zm_PH_DGLM_Results.csv", row.names = F, sep = ",")
