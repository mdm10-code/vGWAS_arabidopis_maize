#####Part II of vGWAS Pipeline: vQTL Simulation Code#####
#!!!!!Important note: If running vGWAS tests on NULL, skip running this file

#As always, load the required libraries for this script
library(dplyr)

#Now, I will break my Arabidopsis genotypic files into separate genotypic and map R objects. The 
#purpose of this is for downstream analyses. Because the map information will be the same, I will 
#create the map from only one numeric R object
at_map <- at_geno_numeric_500[, 1:5] #map information
at_just_geno_500 <- at_geno_numeric_500[, -c(1:5)] #Just genotypic information
at_just_geno_full <- at_geno_numeric_full[, -c(1:5)] #Just genotypic data from the full dataset

#This is for Arabidopsis
##500 diversity panel
at_justones_500 <- matrix(1, nrow = nrow(at_just_geno_500), ncol = 500)
at_just_geno_500 <- at_just_geno_500 + at_justones_500
at_geno_numeric_500_plusone <- data.frame(at_map, at_just_geno_500)

##FULL Diversity Panel
at_justones_full <- matrix(1, nrow = nrow(at_just_geno_full), ncol = 1087)
at_just_geno_full <- at_just_geno_full + at_justones_full
at_geno_numeric_full_plusone <- data.frame(at_map, at_just_geno_full)

#Now I will do the same for maize
zm_map <- zm_geno_numeric_500[, 1:5] #map information
zm_just_geno_500 <- zm_geno_numeric_500[, -c(1:5)] #Just genotypic information
zm_just_geno_full <- zm_geno_numeric_full[, -c(1:5)] #Just genotypic data from the full dataset

##500 diversity panel
zm_justone_500 <- matrix(1, nrow = nrow(zm_just_geno_500), ncol = 500)
Zm_just_geno_500 <- zm_just_geno_500 + zm_justone_500
zm_geno_numeric_500_plusone <- data.frame(zm_map, Zm_just_geno_500)

#For Full
zm_justones_full <- matrix(1, nrow = nrow(zm_geno_numeric_full), ncol = 2532)
Zm_just_geno_full <- zm_just_geno_full + zm_justones_full
zm_geno_numeric_full_plusone <- data.frame(zm_map, Zm_just_geno_full)

##############################vQTL simulations############################

###This is the function that I will be using to simulate
##variance genetic architecture
#wd is your working directory
#Geno is your genotypic object. This needs to be encoded as 0,1,2 or otherwise there will be NA values 
#Geno also needs the map information as well (the first five columns)
#SP_obj is your scaled simplePHENOTYPES R object. I will discuss this further below
#vQTN_sel specifices the options of whether you want to yu vQTNs to be randomly selected or have them 
##as a predefined vector vQTN_vec "defined"
#vQTN_vec is a vector of numbers that corresponds to the column of your Geno object to select vQTNs
#vQTN_num is the number of vQTNs that you want to simulate. This needs to match the length of vQTN_vec
#var_sim is how you assign effect sizes to your vQTNs. If var_sim=='geometric", this assumes a 
#geometric series. If "efvector" 
#randvec is your predefined vQTN effect sizes. Use this argument if you have vQTN_sel=defined
##OR if you have a single vQTN, use this in conjuection with var_sim=="geometric"
#vQTN_es is the effect size of your largest vQTN if var_sim=="geometric"
#h2 is your narrow-sense heritability. 
#purev is whether you are specificying if your trait is only controlled by variance effects. If true,
##h2 needs to be zero
#add_var tells the function whether you have an SP_obj or not
#reps is the number of replications that you want to simulate

vQTL <- function(wd=NULL,
  Dir=NULL, 
  Geno=NULL, 
  SP_obj=NULL,
  vQTN_sel=NULL,
  vQTN_vec=NULL,
  vQTN_num=NULL, 
  var_sim=NULL, 
  randvec=NULL, 
  vQTN_es=NULL, 
  h2=NULL,
  purev=NULL,
  add_var=NULL,
  reps=NULL) {
  setwd(wd)
  dir.create(Dir)
  setwd(Dir)
  if(vQTN_sel=="random"){ #randomly picks vQTNs
    seed.number <- sample(-1000000:1000000, 1)
    seed.number <- seed.number
    set.seed(seed.number)
    these.rows <- 1:nrow(Geno)
    n.vQTN <- vQTN_num
    these.markers.for.QTN <- sample(these.rows, n.vQTN)
    var.QTN.genotypic.information <- Geno[these.markers.for.QTN, ] #same markers as mQTNs and vQTNs
    write.csv(var.QTN.genotypic.information, file = "var.QTN.genotypic.information.csv" , sep = ",", row.names = FALSE)
  }
  if(vQTN_sel == "defined") { #predefined vQTNs
    var.QTN.genotypic.information <- Geno[c(vQTN_vec), ]
    var.QTN.genotypic.information <- data.frame(randvec, var.QTN.genotypic.information)
    write.csv(var.QTN.genotypic.information, file = "var.QTN.genotypic.information.csv", sep = ",", row.names = FALSE)
  }
  variance.effect.trait.object <- t(var.QTN.genotypic.information[, -c(1:6)])
  sigma.component <- as.data.frame(matrix(1, nrow = nrow(variance.effect.trait.object), ncol = 1))
  print(dim(sigma.component))
  if(var_sim == "efvector") { #a vector of custom effect sizes
    for(i in 1:vQTN_num) sigma.component <- sigma.component + ((randvec[i])*(variance.effect.trait.object[, i]))
  }
  if(var_sim == "geometric") { #geometric series
    for(i in 1:vQTN_num) sigma.component <- sigma.component + ((vQTN_es^i) * (variance.effect.trait.object[, i]))
  }
  rownames(sigma.component) <- rownames(variance.effect.trait.object)
  colnames(sigma.component) <- "Variance.effect"
  trait <- data.frame(matrix(NA, nrow = dim(Geno)[2]-5, ncol = reps))
  if(add_var == TRUE) {
    ksq <- ((var(SP_obj[, 2])/h2 - var(SP_obj[, 2]))/((median(sigma.component[, 1])^2)))
    k <- sqrt(ksq)
  } else {
    if(add_var == FALSE)
    k <- 1
  }
  for(i in 1:ncol(Geno) - 5) {
    print(paste0("simulating breeding value for individual ", i))
    trait[i, ] <- k * rnorm(reps, mean = 0, sd = sigma.component[i, ])
    if(purev == FALSE) {
      trait[i, ] <- SP_obj[i, 2:ncol(SP_obj)] + trait[i, ]
    } else {
      if(purev == TRUE)
        print("Simulating a pure variance genetic architecture!")
    }
    #markers selected to be QTNs
  }
  trait <- data.frame(taxa = colnames(Geno)[6:dim(Geno)[2]], trait) #colnames
  write.csv(trait, file = "trait.csv" , sep = ",", row.names = FALSE)
  return(trait)
}

####A. thaliana simulations####
##First I will create a directory in my A. thaliana results folder for vQTL simulations
list.dirs(wd$vGWAS_chapter_one_at)
dir.create(paste0(wd$vGWAS_chapter_one_at, "at_vQTL_simulations/"))
wd$at_vQTL_simulations <- paste0(wd$vGWAS_chapter_one_at, "at_vQTL_simulations/")

#This is for maize
list.dirs(wd$vGWAS_chapter_one_zm)
dir.create(paste0(wd$vGWAS_chapter_one_zm, "zm_vQTL_simulations/"))
wd$zm_vQTL_simulations <- paste0(wd$vGWAS_chapter_one_zm, "zm_vQTL_simulations/")


#In order for the simulations to work for both additive and
#epistatic effects, I will need to scale the additive simplePHENOTYPES
#objects. We will only need to do this for Molybdenum-like settings with a mQTL

##For Arabidopsis
###Molybdenum-like with mQTLs
####500 Diversity Panel
scaled_sp_at_molike_wMQTL_500 <- scale(at_sp_molike_wMQTL_500[, 2:101], center = TRUE, scale = TRUE)
scaled_sp_at_molike_wMQTL_500 <- data.frame(taxa = at_sp_molike_wMQTL_500[, 1], scaled_sp_at_molike_wMQTL_500)

####FULL Diversity Panel
scaled_sp_at_molike_wMQTL_full <- scale(at_sp_molike_wmQTL_full[, 2:101], center = TRUE, scale = TRUE)
scaled_sp_at_molike_wMQTL_full <- data.frame(taxa = at_sp_molike_wmQTL_full[, 1], scaled_sp_at_molike_wMQTL_full)

###Scaling for molike wMQTL
####500 diversity panel
scaled_sp_zm_molike_wMQTL_500 <- scale(sp_zm_molike_wmQTL_500[, 2:101], center = TRUE, scale = TRUE)
scaled_sp_Zm_molike_wMQTL_500 <- data.frame(taxa = sp_zm_molike_wmQTL_500[, 1], scaled_sp_zm_molike_wMQTL_500)

####FULL
scaled_sp_zm_molike_wMQTL_full <- scale(sp_zm_molike_wmQTL_full[, 2:101], center = TRUE, scale = TRUE)
scaled_sp_zm_molike_wMQTL_full <- data.frame(taxa = sp_zm_molike_wmQTL_full[, 1], scaled_sp_zm_molike_wMQTL_full)

#########In this section, I need to calculate an oil-like setting without mQTLs.

vQTN_oil_vector_es <- c(0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10)

#In this section, to get snps with MAF of 0.10 and 0.40, I will need to calculate the MAFs.
##Arabidopsis
###500 diversity panel
at_MAF500 <- apply(t(at_just_geno_500), 2, mean)
at_MAF500 <- matrix(at_MAF500, nrow = 1)
at_MAF500 <- apply(at_MAF500, 2, function(x) min(1-x/2, x/2))
at_MAF500 <- data.frame(at_map[, c(1, 3, 4)], at_MAF500)

###For the full diversity panel
at_MAFFULL <- apply(t(at_just_geno_full), 2, mean)
at_MAFFULL <- matrix(at_MAFFULL, nrow = 1)
at_MAFFULL <- apply(at_MAFFULL, 2, function(x) min(1-x/2, x/2))
at_MAFFULL <- data.frame(at_map[, c(1, 3, 4)], at_MAFFULL)

#This section is for maize
###500 diversity panl
zm_MAF500 <- apply(t(zm_just_geno_500), 2, mean)
zm_MAF500 <- matrix(zm_MAF500, nrow = 1)
zm_MAF500 <- apply(zm_MAF500, 2, function(x) min(1-x/2, x/2))
zm_MAF500 <- data.frame(zm_map[, c(1, 3, 4)], zm_MAF500)

###Full diversity panel
zm_MAFFULL <- apply(t(zm_just_geno_full), 2, mean)
zm_MAFFULL <- matrix(zm_MAFFULL, nrow = 1)
zm_MAFFULL <- apply(zm_MAFFULL, 2, function(x) min(1-x/2, x/2))
zm_MAFFULL <- data.frame(zm_map[, c(1, 3, 4)], zm_MAFFULL)

colnames(at_MAF500) <- c("SNP", "CHR", "POS", "MAF")
colnames(at_MAFFULL) <- c("SNP", "CHR", "POS", "MAF")
colnames(zm_MAF500) <- c("SNP", "CHR", "POS", "MAF")
colnames(zm_MAFFULL) <- c("SNP", "CHR", "POS", "MAF")

###In this section, I am going to select my snps between 0.095-0.105
#and 0.395 to 0.405
#Arabidopsis
#MAF 10
at_MAF500_MAF10 <- at_MAF500 %>%
  select("SNP", "CHR", "POS", "MAF") %>%
  filter(MAF < 0.105 & MAF > 0.095)

at_MAFFULL_MAF10 <- at_MAFFULL %>%
  select("SNP", "CHR", "POS", "MAF") %>%
  filter(MAF < 0.105 & MAF > 0.095)

#MAF 40
at_MAF500_MAF40 <- at_MAF500 %>%
  select("SNP", "CHR", "POS", "MAF") %>%
  filter(MAF < 0.405 & MAF > 0.395)

at_MAFFULL_MAF40 <- at_MAFFULL %>%
  select("SNP", "CHR", "POS", "MAF") %>%
  filter(MAF < 0.405 & MAF > 0.395)

##ZM
###MAF10
zm_MAF500_MAF10 <- zm_MAF500 %>%
  select("SNP", "CHR", "POS", "MAF") %>%
  filter(MAF < 0.105 & MAF > 0.095)

zm_MAFFULL_MAF10 <- zm_MAFFULL %>%
  select("SNP", "CHR", "POS", "MAF") %>%
  filter(MAF < 0.105 & MAF > 0.095)

#MAF 40
zm_MAF500_MAF40 <- zm_MAF500 %>%
  select("SNP", "CHR", "POS", "MAF") %>%
  filter(MAF < 0.405 & MAF > 0.395)

zm_MAFFULL_MAF40 <- zm_MAFFULL %>%
  select("SNP", "CHR", "POS", "MAF") %>%
  filter(MAF < 0.405 & MAF > 0.395)

###Now, extract out snps from the MAF objects for my simulations
##MAF 10
at_rows_for_snps500_MAF10 <- sample(nrow(at_MAF500_MAF10), 17, replace = FALSE)
at_rows_for_snpsfull_MAF10 <- sample(nrow(at_MAFFULL_MAF10), 17, replace = FALSE)

snps_forat500oil_MAF10 <- at_MAF500_MAF10[c(at_rows_for_snps500_MAF10), ]
snps_foratfulloil_MAF10 <- at_MAFFULL_MAF10[c(at_rows_for_snpsfull_MAF10), ]
snp_500_MAF10 <- as.character(snps_forat500oil_MAF10$SNP)
snp_full_MAF10 <- as.character(snps_foratfulloil_MAF10$SNP)

rows_for_snps_at500_MAF10 <- grep(TRUE, at_geno_numeric_500_plusone$snp %in% snp_500_MAF10)
rows_for_snps_atfull_MAF10 <- grep(TRUE, at_geno_numeric_full_plusone$snp %in% snp_full_MAF10)

#This randomly mixes the snps up
rows_for_snps_at500_MAF10 <- sample(rows_for_snps_at500_MAF10, 17, replace = FALSE)
#The rows used 23775 28028 11801 30263  4194 39371 13216 29324 20799 31382 31128 10433 29863 24972 
#39975 25037 23390

rows_for_snps_atfull_MAF10 <- sample(rows_for_snps_atfull_MAF10, 17, replace = FALSE)
#the rows used were 40816  5028   895 18996 32362 13703 17707 30706 10490 27512  8094
#9832 36827 26743 35424 40655 22123

##MAF 0.40
at_rows_for_snps500_MAF40 <- sample(nrow(at_MAF500_MAF40), 17, replace = FALSE)
at_rows_for_snpsfull_MAF40 <- sample(nrow(at_MAFFULL_MAF40), 17, replace = FALSE)

snps_forat500oil_MAF40 <- at_MAF500_MAF40[c(at_rows_for_snps500_MAF40), ]
snps_foratfulloil_MAF40 <- at_MAFFULL_MAF40[c(at_rows_for_snpsfull_MAF40), ]
snp_500_MAF40 <- as.character(snps_forat500oil_MAF40$SNP)
snp_full_MAF40 <- as.character(snps_foratfulloil_MAF40$SNP)

rows_for_snps_at500_MAF40 <- grep(TRUE, at_geno_numeric_500_plusone$snp %in% snp_500_MAF40)
rows_for_snps_atfull_MAF40 <- grep(TRUE, at_geno_numeric_full_plusone$snp %in% snp_full_MAF40)

#This mixes the snps up
rows_for_snps_at500_MAF40 <- sample(rows_for_snps_at500_MAF10, 17, replace = FALSE)
#snps that were used: 31899 26910  1687 10902 24395 13289 38494 27338 35654 32879 33451
#10091 32212 36030 678 26480 33913

rows_for_snps_atfull_MAF40 <- sample(rows_for_snps_atfull_MAF10, 17, replace = FALSE)
#snps that were used: 22343  7226 11839 27838  4653  4359 20400 37420  3315 17159 20662
#16218 24499 273 7601 17530 14378

###Zm
###MAF10
zm_rows_for_snps500_MAF10 <- sample(nrow(zm_MAF500_MAF10), 17, replace = FALSE)
zm_rows_for_snpsfull_MAF10 <- sample(nrow(zm_MAFFULL_MAF10), 17, replace = FALSE)

snps_forzm500oil_MAF10 <- zm_MAF500_MAF10[c(zm_rows_for_snps500_MAF10), ]
snps_forzmfulloil_MAF10 <- zm_MAFFULL_MAF10[c(zm_rows_for_snpsfull_MAF10), ]
zm_snp_500_MAF10 <- as.character(snps_forzm500oil_MAF10$SNP)
zm_snp_full_MAF10 <- as.character(snps_forZmfulloil_MAF10$SNP)

rows_for_snps_zm500_MAF10 <- grep(TRUE, zm_geno_numeric_500_plusone$snp %in% zm_snp_500_MAF10)
rows_for_snps_zmfull_MAF10 <- grep(TRUE, zm_geno_numeric_full_plusone$snp %in% zm_snp_full_MAF10)

#This mixes the snps up
rows_for_snps_Zm500_MAF10 <- sample(rows_for_snps_Zm500_MAF10, 17, replace = FALSE)
#SNPS used: 4361 12509 34266  9486  6168 20727 27876 22050 16574 25928 23811
#14332 28872  5307 16736 25245 36567

rows_for_snps_Zmfull_MAF10 <- sample(rows_for_snps_Zmfull_MAF10, 17, replace = FALSE)
#SNPs used: 36751  4291 37359  2526 13464 37984  6786 11368 15434 20933 37179
#36751  4291 37359  2526 13464 37984  6786 11368 15434 20933 37179

#MAF40
zm_rows_for_snps500_MAF40 <- sample(nrow(zm_MAF500_MAF40), 17, replace = FALSE)
zm_rows_for_snpsfull_MAF40 <- sample(nrow(zm_MAFFULL_MAF40), 17, replace = FALSE)

snps_forzm500oil_MAF40 <- zm_MAF500_MAF40[c(zm_rows_for_snps500_MAF40), ]
snps_forzmfulloil_MAF40 <- zm_MAFFULL_MAF40[c(zm_rows_for_snpsfull_MAF40), ]
zm_snp_500_MAF40 <- as.character(snps_forzm500oil_MAF40$SNP)
zm_snp_full_MAF40 <- as.character(SNPs_forZmFULLOil_MAF40$SNP)

rows_for_snps_zm500_MAF40 <- grep(TRUE, zm_geno_numeric_500_plusone$snp %in% zm_snp_500_MAF40)
rows_for_snps_ZmFULL_MAF40 <- grep(TRUE, zm_geno_numeric_full_plusone$snp %in% zm_snp_full_MAF40)

rows_for_snps_zm500_MAF40 <- sample(rows_for_snps_zm500_MAF40, 17, replace = FALSE)
#SNPS used: 8266  6206 30120 21679  3978 16616  2084 19595  9443 31223 17007
#13267  9126 33820 10522  3095 30922

rows_for_snps_zmfull_MAF40 <- sample(rows_for_snps_zmfull_MAF40, 17, replace = FALSE)
#SNPS used: 13009 36648 23006  8308  5980  2130 35853 15714 34305 22444 34541
#304 17263  5726 38418  6589 30227

#vQTL simulations
##Arabidopsis
##Molike_MAF40
###500
####Because I am controlling for MAF, I need to set aside one snp. These were the SNPS used
rows_for_snps_at500_molike <- 10091
rows_for_snps_atfull_molike <- 20662

at_vpheno_molike_MAF40_500 <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "At_vmolike_MAF40_500", 
  Geno = at_geno_numeric_500_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_500,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_at500_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
at_vpheno_molike_MAF40_full <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vmolike_MAF40_full", 
  Geno = at_geno_numeric_full_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_full,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_atfull_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

#Molybdenum-like womQTL
##Sub500
at_vpheno_molike_womQTL_MAF40_500 <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vmolike_womQTL_MAF40_500", 
  Geno = at_geno_numeric_500_plusone, 
  SP_obj = NULL,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_at500_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  h2 = 0,
  purev = TRUE,
  add_var = FALSE,
  reps = 100)

##FULL
at_vpheno_molike_womQTL_MAF40_full <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vmolike_womQTL_MAF40_full", 
  Geno = at_geno_numeric_full_plusone, 
  SP_obj = NULL,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_atfull_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  h2 = 0,
  purev = TRUE,
  add_var = FALSE,
  reps = 100)

#oillike setting without mQTLs MAF10
##500 diversity panel
at_vpheno_oil_womQTL_MAF10_500 <- vQTL(wd = wd$at_vQTL_simulations,
  Dir = "at_voillike_womQTL_MAF10_500", 
  Geno = at_geno_numeric_500_plusone, 
  SP_obj = NULL,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_at500_MAF10, 
  vQTN_num = 17, 
  var_sim = "efvector", 
  randvec = vQTN_oil_vector_es, 
  vQTN_es = NULL,
  purev = TRUE,
  add_var = FALSE,
  h2 = 0.0,
  reps = 100)

##full diversity panel
at_vpheno_oil_womQTL_MAF10_full <- vQTL(wd = wd$at_vQTL_simulations,
  Dir = "at_voillike_womQTL_MAF10_full", 
  Geno = at_geno_numeric_full_plusone, 
  SP_obj = NULL,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_atfull_MAF10, 
  vQTN_num = 17, 
  var_sim = "efvector", 
  randvec = vQTN_oil_vector_es, 
  vQTN_es = NULL,
  purev = TRUE,
  add_var = FALSE,
  h2 = 0.0,
  reps = 100)

#Oillike setting with vQTNs of MAF40
##500 diversity panel
at_vpheno_oil_womQTL_MAF40_500 <- vQTL(wd = wd$at_vQTL_simulations,
  Dir = "at_voillike_womQTL_MAF40_500", 
  Geno = at_geno_numeric_500_plusone, 
  SP_obj = NULL,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_at500_MAF40, 
  vQTN_num = 17, 
  var_sim = "efvector", 
  randvec = vQTN_oil_vector_es, 
  vQTN_es = NULL,
  purev = TRUE,
  add_var = FALSE,
  h2 = 0,
  reps = 100)

##full diversity panel
at_vpheno_oil_womQTL_MAF40_full <- vQTL(wd = wd$at_vQTL_simulations,
  Dir = "at_voillike_womQTL_MAF40_full", 
  Geno = at_geno_numeric_full_plusone, 
  SP_obj = NULL,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_atfull_MAF40, 
  vQTN_num = 17, 
  var_sim = "efvector", 
  randvec = vQTN_oil_vector_es, 
  vQTN_es = NULL,
  purev = TRUE,
  add_var = FALSE,
  h2 = 0,
  reps = 100)

##################Maize##############################################
#Molybdenum-like setting
##pick a snp that has MAF40
rows_for_snps_zm500_molike <- 6206
rows_for_snps_zmfull_molike <- 34305

#500 diversity panel
zm_vpheno_molike_MAF40_500 <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF40_500", 
  Geno = zm_geno_numeric_500_plusone, 
  SP_obj = scaled_sp_zm_molike_wMQTL_500,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_zm500_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  purev = FALSE,
  add_var = TRUE,
  h2 = 0.63,
  reps = 100)

#full diversity panel
zm_vPheno_molike_MAF40_full <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF40_full", 
  Geno = zm_geno_numeric_full_plusone, 
  SP_obj = scaled_sp_zm_molike_wMQTL_full,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_zmfull_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  purev = FALSE,
  add_var = TRUE,
  h2 = 0.63,
  reps = 100)

#Molybdenum-like without mQTLs
##500 diversity panel
zm_vpheno_molike_womQTL_MAF40_500 <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_womQTL_MAF40_500", 
  Geno = zm_geno_numeric_500_plusone, 
  SP_obj = NULL,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_zm500_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  purev = TRUE,
  add_var = FALSE,
  h2 = 0.0,
  reps = 100)

#full diversity panel
zm_vpheno_molike_womQTL_MAF40_full <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_womQTL_MAF40_full", 
  Geno = zm_geno_numeric_full_plusone, 
  SP_obj = NULL,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_zmfull_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  purev = TRUE,
  add_var = FALSE,
  h2 = 0.0,
  reps = 100)

#maize
##MAF10
###Sub500
zm_vpheno_oil_womQTL_MAF10_500 <- vQTL(wd = wd$zm_vQTL_simulations,
  Dir = "zm_voillike_womQTL_MAF10_500", 
  Geno = zm_geno_numeric_500_plusone, 
  SP_obj = NULL,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_zm500_MAF10, 
  vQTN_num = 17, 
  var_sim = "efvector", 
  randvec = vQTN_oil_vector_es, 
  vQTN_es = NULL,
  purev = TRUE,
  add_var = FALSE,
  h2 = 0,
  reps = 100)

###FULL
zm_vpheno_oil_womQTL_MAF10_full <- vQTL(wd = wd$zm_vQTL_simulations,
  Dir = "zm_voillike_womQTL_MAF10_full", 
  Geno = zm_geno_numeric_full_plusone, 
  SP_obj = NULL,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_zmfull_MAF10, 
  vQTN_num = 17, 
  var_sim = "efvector", 
  randvec = vQTN_oil_vector_es, 
  vQTN_es = NULL,
  purev = TRUE,
  add_var = FALSE,
  h2 = 0,
  reps = 100)

#oillike settings with MAF40
zm_vpheno_oil_womQTL_MAF40_500 <- vQTL(wd = wd$zm_vQTL_simulations,
  Dir = "zm_voillike_womQTL_MAF40_500", 
  Geno = zm_geno_numeric_500_plusone, 
  SP_obj = NULL,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_zm500_MAF40, 
  vQTN_num = 17, 
  var_sim = "efvector", 
  randvec = vQTN_oil_vector_es, 
  vQTN_es = NULL,
  purev = TRUE,
  add_var = FALSE,
  h2 = 0,
  reps = 100)

#full diversity panel
zm_vpheno_oil_womQTL_MAF40_full <- vQTL(wd = wd$zm_vQTL_simulations,
  Dir = "zm_voillike_womQTL_MAF40_full", 
  Geno = zm_geno_numeric_full_plusone, 
  SP_obj = NULL,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_zmfull_MAF40, 
  vQTN_num = 17, 
  var_sim = "efvector", 
  randvec = vQTN_oil_vector_es, 
  vQTN_es = NULL,
  purev = TRUE,
  add_var = FALSE,
  h2 = 0,
  reps = 100)