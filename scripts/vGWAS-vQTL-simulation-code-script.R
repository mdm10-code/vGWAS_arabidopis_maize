####################################################################################################
##################################### vQTL Simulation Script #######################################
####################################################################################################

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
#Maize filtered
zm.map <- zm.geno.numeric.full[, 1:5] #map information
zm.just.geno.500 <- zm.geno.numeric.500[, -c(1:5)] #Just genotypic information
zm.just.geno.full <- zm.geno.numeric.full[, -c(1:5)] #Just genotypic data from the full dataset

#For filtered maize
##500 diversity panel
zm.justone.500 <- matrix(1, nrow = nrow(zm.just.geno.500), ncol = 500)
zm.just.geno.500 <- zm.just.geno.500 + zm.justone.500
zm.geno.500.plusone <- data.frame(zm.map, zm.just.geno.500)

#For Full
zm.justone.FULL <- matrix(1, nrow = nrow(zm.just.geno.full), ncol = 2815)
zm.just.geno.full <- zm.just.geno.full + zm.justone.FULL
zm.geno.full.plusone <- data.frame(zm.map, zm.just.geno.full)

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
  if (vQTN_sel == "random") { #randomly picks vQTNs
    seed.number <- sample(-1000000:1000000, 1)
    seed.number <- seed.number
    set.seed(seed.number)
    these.rows <- 1:nrow(Geno)
    n.vQTN <- vQTN_num
    these.markers.for.QTN <- sample(these.rows, n.vQTN)
    var.QTN.genotypic.information <- Geno[these.markers.for.QTN, ] #same markers as mQTNs and vQTNs
    write.csv(var.QTN.genotypic.information, file = "var.QTN.genotypic.information.csv", sep = ",", row.names = FALSE)
  }
  if (vQTN_sel == "defined") { #predefined vQTNs
    var.QTN.genotypic.information <- Geno[c(vQTN_vec), ]
    var.QTN.genotypic.information <- data.frame(randvec, var.QTN.genotypic.information)
    write.csv(var.QTN.genotypic.information, file = "var.QTN.genotypic.information.csv", sep = ",", row.names = FALSE)
  }
  variance.effect.trait.object <- t(var.QTN.genotypic.information[, -c(1:6)])
  sigma.component <- as.data.frame(matrix(1, nrow = nrow(variance.effect.trait.object), ncol = 1))
  print(dim(sigma.component))
  if (var_sim == "efvector") { #a vector of custom effect sizes
    for(i in 1:vQTN_num) sigma.component <- sigma.component + ((randvec[i]) * (variance.effect.trait.object[, i]))
  }
  if (var_sim == "geometric") { #geometric series
    for(i in 1:vQTN_num) sigma.component <- sigma.component + ((vQTN_es^i) * (variance.effect.trait.object[, i]))
  }
  rownames(sigma.component) <- rownames(variance.effect.trait.object)
  colnames(sigma.component) <- "Variance.effect"
  trait <- data.frame(matrix(NA, nrow = dim(Geno)[2] - 5, ncol = reps))
  if (add_var == TRUE) {
    ksq <- ((var(SP_obj[, 2]) / h2 - var(SP_obj[, 2])) / ((median(sigma.component[, 1]) ^ 2)))
    k <- sqrt(ksq)
  } else {
    if (add_var == FALSE)
      k <- 1
  }
  for (i in 1:ncol(Geno) - 5) {
    print(paste0("simulating breeding value for individual ", i))
    trait[i, ] <- k * rnorm(reps, mean = 0, sd = sigma.component[i, ])
    if (purev == FALSE) {
      trait[i, ] <- SP_obj[i, 2:ncol(SP_obj)] + trait[i, ]
    } else {
      if (purev == TRUE)
        print("Simulating a pure variance genetic architecture!")
    }
    #markers selected to be QTNs
  }
  trait <- data.frame(taxa = colnames(Geno)[6:dim(Geno)[2]], trait) #colnames
  write.csv(trait, file = "trait.csv", sep = ",", row.names = FALSE)
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
###Scaling for molike wMQTL MAF 0.40 filtered like Arabidopsis for Maize
####500
scaled.sp.zm.molike.MAF40.500 <- scale(sp.zm.molike.500[, 2:101], center = TRUE, scale = TRUE)
scaled.sp.zm.molike.MAF40.500 <- data.frame(taxa = sp.zm.molike.500[, 1], scaled.sp.zm.molike.MAF40.500)

####FULL Diversity Panel
scaled.sp.zm.molike.MAF40.FULL <- scale(sp.zm.molike.FULL[, 2:101], center = TRUE, scale = TRUE)
scaled.sp.zm.molike.MAF40.FULL <- data.frame(taxa = sp.zm.molike.FULL[, 1], scaled.sp.zm.molike.MAF40.FULL)

#In this section, to get snps with MAF of 0.10 and 0.40, I will need to calculate the MAFs.
##Arabidopsis
###500 diversity panel
at_MAF500 <- apply(t(at_just_geno_500), 2, mean)
at_MAF500 <- matrix(at_MAF500, nrow = 1)
at_MAF500 <- apply(at_MAF500, 2, function(x) min(1 - x / 2, x / 2))
at_MAF500 <- data.frame(at_map[, c(1, 3, 4)], at_MAF500)

###For the full diversity panel
at_MAFFULL <- apply(t(at_just_geno_full), 2, mean)
at_MAFFULL <- matrix(at_MAFFULL, nrow = 1)
at_MAFFULL <- apply(at_MAFFULL, 2, function(x) min(1 - x / 2, x / 2))
at_MAFFULL <- data.frame(at_map[, c(1, 3, 4)], at_MAFFULL)

#This section is for maize
#Maize Filtered like Arabidopsis
###500 diversity panl
zm.MAF500 <- apply(t(zm.just.geno.500), 2, mean)
zm.MAF500 <- matrix(zm.MAF500, nrow = 1)
zm.MAF500 <- apply(zm.MAF500, 2, function(x) min(1 - x / 2, x / 2))
zm.MAF500 <- data.frame(zm.map[, c(1, 3, 4)], zm.MAF500)

###Full diversity panel
zm.MAFFULL <- apply(t(zm.just.geno.full), 2, mean)
zm.MAFFULL <- matrix(zm.MAFFULL, nrow = 1)
zm.MAFFULL <- apply(zm.MAFFULL, 2, function(x) min(1 - x / 2, x / 2))
zm.MAFFULL <- data.frame(zm.map[, c(1, 3, 4)], zm.MAFFULL)

#
colnames(at_MAF500) <- c("SNP", "CHR", "POS", "MAF")
colnames(at_MAFFULL) <- c("SNP", "CHR", "POS", "MAF")

colnames(zm.MAF500) <- c("SNP", "CHR", "POS", "MAF")
colnames(zm.MAFFULL) <- c("SNP", "CHR", "POS", "MAF")


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
#MAF 40 (for filtered maize)
zm.MAF500.MAF40 <- zm.MAF500 %>%
  select("SNP", "CHR", "POS", "MAF") %>%
  filter(MAF < 0.405 & MAF > 0.395)

zm.MAFFULL.MAF40 <- zm.MAFFULL %>%
  select("SNP", "CHR", "POS", "MAF") %>%
  filter(MAF < 0.405 & MAF > 0.395)

#MAF10 (for filtered maize)
zm.MAF500.MAF10 <- zm.MAF500 %>%
  select("SNP", "CHR", "POS", "MAF") %>%
  filter(MAF < 0.105 & MAF > 0.095)

zm.MAFFULL.MAF10 <- zm.MAFFULL %>%
  select("SNP", "CHR", "POS", "MAF") %>%
  filter(MAF < 0.105 & MAF > 0.095)

###Now, extract out snps from the MAF objects for my simulations

###Zm
###MAF10
zm.rows.for.snps500.MAF10 <- sample(nrow(zm.MAF500.MAF10), 1, replace = FALSE)
zm.rows.for.snpsfull.MAF10 <- sample(nrow(zm.MAFFULL.MAF10), 1, replace = FALSE)

rows.for.snps.zm500.MAF10 <- grep(TRUE, zm.geno.numeric.500$snp %in% zm.MAF500.MAF10[c(zm.rows.for.snps500.MAF10), 1])
rows.for.snps.zmfull.MAF10 <- grep(TRUE, zm.geno.numeric.full$snp %in% zm.MAFFULL.MAF10[c(zm.rows.for.snpsfull.MAF10), 1])

#MAF40
zm.rows.for.snps500.MAF40 <- sample(nrow(zm.MAF500.MAF40), 1, replace = FALSE)
zm.rows.for.snpsfull.MAF40 <- sample(nrow(zm.MAFFULL.MAF40), 1, replace = FALSE)

rows.for.snps.zm500.MAF40 <- grep(TRUE, zm.geno.numeric.500$snp %in% zm.MAF500.MAF40[c(zm.rows.for.snps500.MAF40), 1])
rows.for.snps.zmfull.MAF40 <- grep(TRUE, zm.geno.numeric.full$snp %in% zm.MAFFULL.MAF40[c(zm.rows.for.snpsfull.MAF40), 1])

#vQTL simulations
##Arabidopsis
##Molike_MAF40
###500
####Because I am controlling for MAF, I need to set aside one snp. These were the SNPS used
rows_for_snps_at500_molike <- sample(nrow(at_MAF500_MAF40), 1, replace = FALSE)

#Let's see what snp was randomly selected
at_MAF500_MAF40[80, ]

#now, let's see what where this snp is at in the general dataset
grep(TRUE, at_MAF500$SNP %in% "S1_26529682")

rows_for_snps_at500_molike <- 10091

#Now we do this for the full dataste
rows_for_snps_atFULL_molike <- sample(nrow(at_MAFFULL_MAF40), 1, replace = FALSE)

#Let's see what snp was randomly selected
at_MAFFULL_MAF40[c(rows_for_snps_atFULL_molike), ]

#now, let's see what where this snp is at in the general dataset
grep(TRUE, at_MAFFULL$SNP %in% "S3_10180078")

rows_for_snps_at500_molike <- 10091 #These numbers are from the grep function as seen above
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

#Molike setting h2 = 0.33, MAF=0.40, vQTN = 0.10
##500
at_vpheno_molike_MAF40_h233_vQTN10_500 <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF40_h233_vQTN10_500", 
  Geno = at_geno_numeric_500_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_500,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_at500_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
at_vpheno_molike_MAF40_h233_vQTN10_full <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF40_h233_vQTN10_full", 
  Geno = at_geno_numeric_full_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_full,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_atfull_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

#Molike setting h2 = 0.33, MAF=0.40, vQTN = 0.50
##500
at_vpheno_molike_MAF40_h233_vQTN50_500 <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF40_h233_vQTN50_500", 
  Geno = at_geno_numeric_500_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_500,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_at500_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
at_vpheno_molike_MAF40_h233_vQTN50_full <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF40_h233_vQTN50_full", 
  Geno = at_geno_numeric_full_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_full,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_atfull_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

#Molike setting h2 = 0.33, MAF = 0.40, vQTN = 0.90
##500
at_vpheno_molike_MAF40_h233_vQTN90_500 <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF40_h233_vQTN90_500", 
  Geno = at_geno_numeric_500_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_500,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_at500_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
at_vpheno_molike_MAF40_h233_vQTN90_full <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF40_h233_vQTN90_full", 
  Geno = at_geno_numeric_full_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_full,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_atfull_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

#Molike setting h2 = 0.63, MAF = 0.40, vQTN = 0.10
##500
at_vpheno_molike_MAF40_h263_vQTN10_500 <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF40_h263_vQTN10_500", 
  Geno = at_geno_numeric_500_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_500,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_at500_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
at_vpheno_molike_MAF40_h263_vQTN10_full <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF40_h263_vQTN10_full", 
  Geno = at_geno_numeric_full_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_full,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_atfull_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

#Molike setting h2 = 0.63, MAF = 0.40, vQTN = 0.50
##500
at_vpheno_molike_MAF40_h263_vQTN50_500 <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF40_h263_vQTN50_500", 
  Geno = at_geno_numeric_500_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_500,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_at500_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
at_vpheno_molike_MAF40_h263_vQTN50_full <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF40_h263_vQTN50_full", 
  Geno = at_geno_numeric_full_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_full,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_atfull_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

#Just like I did for a sinlge QTN for molike, I will need to select one
#snp with MAF of 0.10. The additive QTN from simplePHENOTYPES is S4_9181327
#We need to avoid a vQTN within this additive QTN
rows_for_snps_at500_molike_MAF10 <- sample(nrow(at_MAF500_MAF10), 1, replace = FALSE)

#Let's see what snp was randomly selected
at_MAF500_MAF10[1128, ]

#now, let's see what where this snp is at in the general dataset
grep(TRUE, at_MAF500$SNP %in% "S5_24850827")

rows_for_snps_at500_molike_MAF10 <- 40830

#Now we do this for the full dataste
rows_for_snps_atFULL_molike_MAF10 <- sample(nrow(at_MAFFULL_MAF10), 1, replace = FALSE)

#Let's see what snp was randomly selected
at_MAFFULL_MAF10[c(rows_for_snps_atFULL_molike_MAF10), ]

#now, let's see what where this snp is at in the general dataset
grep(TRUE, at_MAFFULL$SNP %in% "S1_16715058")

rows_for_snps_atfull_molike_MAF10 <- 5044

#Molike setting h2 = 0.33, MAF = 0.10, vQTN = 0.10 #Need one snp with allele of MAF 0.10
##500
at_vpheno_molike_MAF10_h233_vQTN10_500 <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF10_h233_vQTN10_500", 
  Geno = at_geno_numeric_500_plusone, 
  SP_obj = scaled_SP_At_molike_wmQTL_MAF10_root_sub500,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_at500_molike_MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
at_vpheno_molike_MAF40_h263_vQTN50_full <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF40_h263_vQTN50_full", 
  Geno = at_geno_numeric_full_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_full,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_atfull_molike, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

#Just like I did for a sinlge QTN for molike, I will need to select one
#snp with MAF of 0.10. The additive QTN from simplePHENOTYPES is S4_9181327
#We need to avoid a vQTN within this additive QTN

#Let's see what snp was randomly selected
at_MAFFULL_MAF10[c(rows_for_snps_atFULL_molike_MAF10), ]

#now, let's see what where this snp is at in the general dataset
grep(TRUE, at_MAFFULL$SNP %in% "S1_16715058")

rows_for_snps_atfull_molike_MAF10 <- 5044



#Molike MAF 40 for additive, but MAF 10 for vQTN
at_vpheno_molike_MAF10_h233_vQTN10_500 <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF10_h233_vQTN10_500", 
  Geno = at_geno_numeric_500_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_500,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_at500_molike_MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

##FULL
at_vpheno_molike_MAF10_h233_vQTN10_FULL <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF10_h233_vQTN10_FULL", 
  Geno = at_geno_numeric_full_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_full,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_atfull_molike_MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

#Molike MAF 40 for additive, but MAF 10 for vQTN and effect size 0.5
at_vpheno_molike_MAF10_h233_vQTN50_500 <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF10_h233_vQTN50_500", 
  Geno = at_geno_numeric_500_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_500,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_at500_molike_MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

##FULL
at_vpheno_molike_MAF10_h233_vQTN50_FULL <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF10_h233_vQTN50_FULL", 
  Geno = at_geno_numeric_full_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_full,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_atfull_molike_MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

#Molike MAF 40 for additive, but MAF 10 for vQTN and effect size 0.5
at_vpheno_molike_MAF10_h233_vQTN90_500 <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF10_h233_vQTN90_500", 
  Geno = at_geno_numeric_500_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_500,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_at500_molike_MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

##FULL
at_vpheno_molike_MAF10_h233_vQTN90_FULL <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF10_h233_vQTN90_FULL", 
  Geno = at_geno_numeric_full_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_full,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_atfull_molike_MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

#Now, This is with h2 0.63 and vQTN MAF of 0.10 with ef of 0.1
#Molike MAF 40 for additive, but MAF 10 for vQTN and effect size 0.5
at_vpheno_molike_MAF10_h263_vQTN10_500 <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF10_h263_vQTN10_500", 
  Geno = at_geno_numeric_500_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_500,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_at500_molike_MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

##FULL
at_vpheno_molike_MAF10_h263_vQTN10_FULL <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF10_h263_vQTN10_FULL", 
  Geno = at_geno_numeric_full_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_full,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_atfull_molike_MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

#Low MAF of vQTN, h2=0.63, ef of 0.5
at_vpheno_molike_MAF10_h263_vQTN50_500 <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF10_h263_vQTN50_500", 
  Geno = at_geno_numeric_500_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_500,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_at500_molike_MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

##FULL
at_vpheno_molike_MAF10_h263_vQTN50_FULL <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF10_h263_vQTN50_FULL", 
  Geno = at_geno_numeric_full_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_full,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_atfull_molike_MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

#h2 = 0.63, MAF 0f 0.10 for vQTN, ef = 0.9 for vQTN
at_vpheno_molike_MAF10_h263_vQTN90_500 <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF10_h263_vQTN90_500", 
  Geno = at_geno_numeric_500_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_500,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_at500_molike_MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

##FULL
at_vpheno_molike_MAF10_h263_vQTN90_FULL <- vQTL(
  wd = wd$at_vQTL_simulations,
  Dir = "at_vpheno_molike_MAF10_h263_vQTN90_FULL", 
  Geno = at_geno_numeric_full_plusone, 
  SP_obj = scaled_sp_at_molike_wMQTL_full,
  vQTN_sel = "defined",
  vQTN_vec = rows_for_snps_atfull_molike_MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###############################################Maize################################################
#Molybdenum-like setting

#500 diversity panel
zm.vpheno.molike.MAF40.500 <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF40_500", 
  Geno = zm.geno.500.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.500,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zm500.MAF40, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  purev = FALSE,
  add_var = TRUE,
  h2 = 0.63,
  reps = 100)

#full diversity panel
zm.vPheno.molike.MAF40.full <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF40_full", 
  Geno = zm.geno.full.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.FULL,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zmfull.MAF40, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  purev = FALSE,
  add_var = TRUE,
  h2 = 0.63,
  reps = 100)

#Molike MAF40 vQTN ef 0.1 h2 0.33 wmQTL
zm.vpheno.molike.MAF40.h233.vQTN10.500 <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF40_h233_vQTN10_500_filtlike_at", 
  Geno = zm.geno.500.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.500,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zm500.MAF40, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
zm.vpheno.molike.MAF40.h233.vQTN10.FULL <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF40_h233_vQTN10_FULL", 
  Geno = zm.geno.full.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.FULL,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zmfull.MAF40, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

#Molike MAF40 vQTN ef 0.9 h2 0.33 wmQTL
zm.vpheno.molike.MAF40.h233.vQTN90.500 <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF40_h233_vQTN90_500", 
  Geno = zm.geno.500.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.500,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zm500.MAF40, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
zm.vpheno.molike.MAF40.h233.vQTN90.FULL <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF40_h233_vQTN90_FULL", 
  Geno = zm.geno.full.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.FULL,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zmfull.MAF40, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

##Molike MAF40 h233 vQTN 50
zm.vpheno.molike.MAF40.h233.vQTN50.500 <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF40_h233_vQTN50_500", 
  Geno = zm.geno.500.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.500,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zm500.MAF40, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
zm.vpheno.molike.MAF40.h233.vQTN50.FULL <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF40_h233_vQTN50_FULL", 
  Geno = zm.geno.full.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.FULL,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zmfull.MAF40, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

#Molike h263 vQTN10
##500
zm.vpheno.molike.MAF40.h263.vQTN10.500 <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF40_h263_vQTN10_500", 
  Geno = zm.geno.500.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.500,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zm500.MAF40, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
zm.vpheno.molike.MAF40.h263.vQTN10.FULL <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF40_h263_vQTN10_FULL", 
  Geno = zm.geno.full.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.FULL,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zmfull.MAF40, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

##Molike MAF40 h263 vQTN50
##500
zm.vpheno.molike.MAF40.h263.vQTN50.500 <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF40_h263_vQTN50_500", 
  Geno = zm.geno.500.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.500,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zm500.MAF40, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
zm.vpheno.molike.MAF40.h263.vQTN50.FULL <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF40_h263_vQTN50_FULL", 
  Geno = zm.geno.full.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.FULL,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zmfull.MAF40, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

##Molike settings where vQTN has MAF of 10
##500
zm.vpheno.molike.MAF10.h233.vQTN10.500 <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF10_h233_vQTN10_500", 
  Geno = zm.geno.500.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.500,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zm500.MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
zm.vpheno.molike.MAF10.h233.vQTN10.FULL <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF10_h233_vQTN10_FULL", 
  Geno = zm.geno.full.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.FULL,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zmfull.MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

##Molike vQTN 50 MAF10
##500
zm.vpheno.molike.MAF10.h233.vQTN50.500 <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF10_h233_vQTN50_500", 
  Geno = zm.geno.500.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.500,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zm500.MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
zm.vpheno.molike.MAF10.h233.vQTN50.FULL <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF10_h233_vQTN50_FULL", 
  Geno = zm.geno.full.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.FULL,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zmfull.MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

#Molike vQTN 90 has MAF10
##500
zm.vpheno.molike.MAF10.h233.vQTN90.500 <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF10_h233_vQTN90_500", 
  Geno = zm.geno.500.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.500,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zm500.MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
zm.vpheno.molike.MAF10.h233.vQTN90.FULL <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF10_h233_vQTN90_FULL", 
  Geno = zm.geno.full.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.FULL,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zmfull.MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  h2 = 0.33,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###h263 MAF10
zm.vpheno.molike.MAF10.h263.vQTN10.500 <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF10_h263_vQTN10_500", 
  Geno = zm.geno.500.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.500,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zm500.MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
zm.vpheno.molike.MAF10.h263.vQTN10.FULL <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF10_h263_vQTN10_FULL", 
  Geno = zm.geno.full.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.FULL,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zmfull.MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.10, 
  vQTN_es = 0.10,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

##Molike vQTN 50 MAF10
##500
zm.vpheno.molike.MAF10.h263.vQTN50.500 <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF10_h263_vQTN50_500", 
  Geno = zm.geno.500.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.500,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zm500.MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
zm.vpheno.molike.MAF10.h263.vQTN50.FULL <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF10_h263_vQTN50_FULL", 
  Geno = zm.geno.full.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.FULL,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zmfull.MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.50, 
  vQTN_es = 0.50,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

#Molike vQTN 90 has MAF10
##500
zm.vpheno.molike.MAF10.h263.vQTN90.500<- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF10_h263_vQTN90_500", 
  Geno = zm.geno.500.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.500,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zm500.MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)

###FULL
zm.vpheno.molike.MAF10.h263.vQTN90.FULL <- vQTL(
  wd = wd$zm_vQTL_simulations,
  Dir = "zm_vmolike_MAF10_h263_vQTN90_FULL", 
  Geno = zm.geno.full.plusone, 
  SP_obj = scaled.sp.zm.molike.MAF40.FULL,
  vQTN_sel = "defined",
  vQTN_vec = rows.for.snps.zmfull.MAF10, 
  vQTN_num = 1, 
  var_sim = "geometric", 
  randvec = 0.90, 
  vQTN_es = 0.90,
  h2 = 0.63,
  purev = FALSE,
  add_var = TRUE,
  reps = 100)
