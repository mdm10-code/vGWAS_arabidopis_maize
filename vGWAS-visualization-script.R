#############################################################################################
################################### vGWAS Visualization Code ################################
#############################################################################################

### SPECIAL STATEMENT ###
#The result files that were used for the visualises were complied and reformatted in Excel, not R.

install.packages("ggplot2")
install.packages("ggthemes") #Install 

library("ggplot2")
library("tidyr")
library("ggthemes") # Load

###Read in rate files
wd <- list()
wd$three_output <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/C1/3_output/"
list.files(wd$three_output)

##For my graphic folder, I am going to create a setting working directory shortcut to 
##there
wd$three_output_figures <- paste0(wd$three_output, "Figs/")

###This is the Master File for all the rates
TnF_rates_master <- read.csv(paste0(wd$three_output, "TnF.Rates.Master.1.11.2021.csv"))

#
TnF_rates <- read.csv(paste0(wd$three_output, "Detection_rates_3_8_2022.csv"))

############################################################
######################### Images ########################
###########################################################

####GxE
at_500_GxE_rates <- TnF_rates_master %>%
  select(Test, Species, Trait, Sample.Size, SNP, SNP.Code.Name, Heritability, Effect.Size, MAF, Window.Size, BP_5_FP, BP_5_TP, BP_10_FP, BP_10_TP, FDR_5_FP, FDR_5_TP, FDR_10_FP, FDR_10_TP) %>%
  filter(Trait == "GxE" & Sample.Size == "500" & Species == "Arabidopsis")

##Arabidopsis N = 500
tiff(file = "True Positive Detection Rates-Arabidopsis 500 GxE.tiff",
     width = 641, height = 376, res = 100)
ggplot() +
  geom_bar(data = at_500_GxE_rates, aes(x = SNP.Code.Name, y = FDR_5_TP, fill = Test), stat = "identity", position = position_dodge()) +
  ylim(0, 1.00) + 
  labs(title = "Arabidpsis (N = 500) GxE Setting") +
  xlab("QTN") +
  ylab("% of Replications with a True Positive") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue", "MLM" = "magenta4"))
dev.off()

at_full_GxE_rates <- TnF_rates_master %>%
  select(Test, Species, Trait, Sample.Size, SNP, SNP.Code.Name, Heritability, Effect.Size, MAF, Window.Size, BP_5_FP, BP_5_TP, BP_10_FP, BP_10_TP, FDR_5_FP, FDR_5_TP, FDR_10_FP, FDR_10_TP) %>%
  filter(Trait == "GxE" & Sample.Size == "FULL" & Species == "Arabidopsis")

##Arabidopsis N = 1087
tiff(file = "True Positive Detection Rates-Arabidopsis FULL (N = 1087) GxE.tiff",
     width = 641, height = 376, res = 100)
ggplot() +
  geom_bar(data = at_full_GxE_rates, aes(x = SNP.Code.Name, y = FDR_5_TP, fill = Test), stat = "identity", position = position_dodge()) +
  ylim(0, 1.00) + 
  labs(title = "Arabidpsis (N = 1087) GxE Setting") +
  xlab("QTN") +
  ylab("% of Replications with a True Positive") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue", "MLM" = "magenta4"))
dev.off()

zm_500_GxE_rates <- TnF_rates %>%
  select(Test, Species, Trait, Sample.Size, SNP, SNP.Code.Name, Heritability, Effect.Size, MAF, Window.Size, BP_5_FP, BP_5_TP, BP_10_FP, BP_10_TP, FDR_5_FP, FDR_5_TP, FDR_10_FP, FDR_10_TP) %>%
  filter(Trait == "GxE" & Sample.Size == "500" & Species == "Maize")

for (i in 1:nrow(zm_500_GxE_rates)) {
  the.results <- exactci(x = zm_500_GxE_rates[i, 16] * 100, n = 100, conf.level = 0.95)
  zm_500_GxE_rates[i, 19] <- the.results$conf.int[1]
  zm_500_GxE_rates[i, 20] <- the.results$conf.int[2]
  rm(the.results)
}

colnames(zm_500_GxE_rates)[19:20] <- c("FDR5_lb", "FDR5_ub")

##Maize N = 500
tiff(file = "True Positive Detection Rates-Maize (N = 500) GxE.tiff",
     width = 641, height = 376, res = 100)
ggplot() +
  geom_bar(data = zm_500_GxE_rates, aes(x = SNP.Code.Name, y = FDR_5_TP, fill = Test), stat = "identity", position = position_dodge()) +
  ylim(0, 1.00) + 
  labs(title = "Maize (N = 500) GxE Setting") +
  xlab("QTN") +
  ylab("% of Replications with a True Positive") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue", "MLM" = "magenta4"))
dev.off()

###
tiff(file = "Maize_N500_True_Positives.tiff",
     width = 800, height = 400, res = 100)
ggplot(zm_500_GxE_rates, aes(x = SNP.Code.Name, y = FDR_5_TP, fill = Test)) +
  geom_bar(colour = "black",  stat = "identity", position = position_dodge(width = .8), size = .3, width = 0.8)+
  geom_errorbar(aes(x = SNP.Code.Name, ymin = FDR5_lb, ymax = FDR5_ub) , width = .2, position = position_dodge(.8))+
#  facet_grid(~ SNP.Code.Name, scales = "free_x", space = "free_x") +
  xlab("QTN") +
  labs(title = "Maize (N = 500) GxE Setting") +
  theme(plot.title = element_text(hjust = 0.50)) +
  ylab("% of Replications with a True Positive") +
  ylim(0, 1.00) + 
#  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") + 
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue"))+
  theme(panel.border = element_rect(fill = NA, linetype = "dashed")) +
  theme_gray(base_size = 12) + theme( panel.border = element_rect(fill = NA, colour = "black"))+
  theme(strip.background =  element_rect(fill = "white", colour = "white"))+
  theme(axis.text.x = element_text(colour = "black",family = "Times", size = 12), axis.ticks = element_blank(), axis.ticks.length = unit(0.1, "cm")) +
  theme(axis.title.y = element_text(vjust = 1.3,colour = "black",family = "Times", size = 12), axis.text.y = element_text(  colour = "black", family = "Times"), axis.ticks = element_blank())+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
dev.off()
###

zm_full_GxE_rates <- TnF_rates %>%
  select(Test, Species, Trait, Sample.Size, SNP, SNP.Code.Name, Heritability, Effect.Size, MAF, Window.Size, BP_5_FP, BP_5_TP, BP_10_FP, BP_10_TP, FDR_5_FP, FDR_5_TP, FDR_10_FP, FDR_10_TP) %>%
  filter(Trait == "GxE" & Sample.Size == "FULL" & Species == "Maize")

for (i in 1:nrow(zm_full_GxE_rates)) {
  the.results <- exactci(x = zm_full_GxE_rates[i, 16] * 100, n = 100, conf.level = 0.95)
  zm_full_GxE_rates[i, 19] <- the.results$conf.int[1]
  zm_full_GxE_rates[i, 20] <- the.results$conf.int[2]
  rm(the.results)
}

colnames(zm_full_GxE_rates)[19:20] <- c("FDR5_lb", "FDR5_ub")


##Maize N = 500
tiff(file = "True Positive Detection Rates-Maize FULL (N = 2532) GxE.tiff",
     width = 641, height = 376, res = 100)
ggplot() +
  geom_bar(data = zm_full_GxE_rates, aes(x = SNP.Code.Name, y = FDR_5_TP, fill = Test), stat = "identity", position = position_dodge()) +
  ylim(0, 1.00) + 
  labs(title = "Maize (N = 2532) GxE Setting") +
  xlab("QTN") +
  ylab("% of Replications with a True Positive") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue", "MLM" = "magenta4"))
dev.off()

###
tiff(file = "Maize_N2815_True_Positives.tiff",
     width = 800, height = 400, res = 100)
ggplot(zm_full_GxE_rates, aes(x = SNP.Code.Name, y = FDR_5_TP, fill = Test)) +
  geom_bar(colour = "black",  stat = "identity", position = position_dodge(width = .8), size = .3, width = 0.8)+
  geom_errorbar(aes(x = SNP.Code.Name, ymin = FDR5_lb, ymax = FDR5_ub) , width = .2, position = position_dodge(.8))+
  #  facet_grid(~ SNP.Code.Name, scales = "free_x", space = "free_x") +
  xlab("QTN") +
  labs(title = "Maize (N = 2815) GxE Setting") +
  theme(plot.title = element_text(hjust = 0.50)) +
  ylab("% of Replications with a True Positive") +
  ylim(0, 1.00) + 
  #  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") + 
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue"))+
  theme(panel.border = element_rect(fill = NA, linetype = "dashed")) +
  theme_gray(base_size = 12) + theme( panel.border = element_rect(fill = NA, colour = "black"))+
  theme(strip.background =  element_rect(fill = "white", colour = "white"))+
  theme(axis.text.x = element_text(colour = "black",family = "Times", size = 12), axis.ticks = element_blank(), axis.ticks.length = unit(0.1, "cm")) +
  theme(axis.title.y = element_text(vjust = 1.3,colour = "black",family = "Times", size = 12), axis.text.y = element_text(  colour = "black", family = "Times"), axis.ticks = element_blank())+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
dev.off()
###

####Epistasis
at_500_Epistasis_rates <- TnF_rates_master %>%
  select(Test, Species, Trait, Sample.Size, SNP, SNP.Code.Name, Heritability, Effect.Size, MAF, Window.Size, BP_5_FP, BP_5_TP, BP_10_FP, BP_10_TP, FDR_5_FP, FDR_5_TP, FDR_10_FP, FDR_10_TP) %>%
  filter(Trait == "Epistasis" & Sample.Size == "500" & Species == "Arabidopsis" & Effect.Size == 0.75)

at_500_Epistasis_rates$Heritability <- gsub(0.30, "H2 = 0.30", at_500_Epistasis_rates$Heritability)
at_500_Epistasis_rates$Heritability <- gsub(0.80, "H2 = 0.80", at_500_Epistasis_rates$Heritability)

##Arabidopsis N = 500
tiff(file = "True Positive Detection Rates-Arabidopsis (N = 500) Epistasis.tiff",
     width = 641, height = 376, res = 100)
ggplot() +
  geom_bar(data = at_500_Epistasis_rates, aes(x = SNP.Code.Name, y = FDR_5_TP, fill = Test), stat = "identity", position = position_dodge()) +
  ylim(0, 1.00) + 
  facet_grid(~ Heritability, scales = "free_x", space = "free_x") +
  labs(title = "Arabidopsis (N = 500) Epistasis Settings") +
  xlab("QTN") +
  ylab("% of Replications with a True Positive") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue", "MLM" = "magenta4"))
dev.off()

#######FULL
at_full_Epistasis_rates <- TnF_rates_master %>%
  select(Test, Species, Trait, Sample.Size, SNP, SNP.Code.Name, Heritability, Effect.Size, MAF, Window.Size, BP_5_FP, BP_5_TP, BP_10_FP, BP_10_TP, FDR_5_FP, FDR_5_TP, FDR_10_FP, FDR_10_TP) %>%
  filter(Trait == "Epistasis" & Sample.Size == "FULL" & Species == "Arabidopsis" & Effect.Size == 0.75)

at_full_Epistasis_rates$Heritability <- gsub(0.30, "H2 = 0.30", at_full_Epistasis_rates$Heritability)
at_full_Epistasis_rates$Heritability <- gsub(0.80, "H2 = 0.80", at_full_Epistasis_rates$Heritability)

##Arabidopsis N = 1087
tiff(file = "True Positive Detection Rates-Arabidopsis (N = 1087) Epistasis.tiff",
     width = 641, height = 376, res = 100)
ggplot() +
  geom_bar(data = at_full_Epistasis_rates, aes(x = SNP.Code.Name, y = FDR_5_TP, fill = Test), stat = "identity", position = position_dodge()) +
  ylim(0, 1.00) + 
  facet_grid(~ Heritability, scales = "free_x", space = "free_x") +
  labs(title = "Arabidopsis FULL (N = 1087) Epistasis Settings") +
  xlab("QTN") +
  ylab("% of Replications with a True Positive") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue", "MLM" = "magenta4"))
dev.off()

###Maize 
##500
zm_500_Epistasis_rates <- TnF_rates %>%
  select(Test, Species, Trait, Sample.Size, SNP, SNP.Code.Name, Heritability, Effect.Size, MAF, Window.Size, BP_5_FP, BP_5_TP, BP_10_FP, BP_10_TP, FDR_5_FP, FDR_5_TP, FDR_10_FP, FDR_10_TP) %>%
  filter(Trait == "Epistasis" & Sample.Size == "500" & Species == "Maize" & Effect.Size == 0.75)

zm_500_Epistasis_rates$Heritability <- gsub(0.30, "H2 = 0.30", zm_500_Epistasis_rates$Heritability)
zm_500_Epistasis_rates$Heritability <- gsub(0.80, "H2 = 0.80", zm_500_Epistasis_rates$Heritability)

for (i in 1:nrow(zm_500_Epistasis_rates)) {
  the.results <- exactci(x = zm_500_Epistasis_rates[i, 16] * 100, n = 100, conf.level = 0.95)
  zm_500_Epistasis_rates[i, 19] <- the.results$conf.int[1]
  zm_500_Epistasis_rates[i, 20] <- the.results$conf.int[2]
  rm(the.results)
}

colnames(zm_500_Epistasis_rates)[19:20] <- c("FDR5_lb", "FDR5_ub")


###
tiff(file = "Maize_N500_Epistasis_True_Positives.tiff",
     width = 800, height = 400, res = 100)
ggplot(zm_500_Epistasis_rates, aes(x = SNP.Code.Name, y = FDR_5_TP, fill = Test)) +
  geom_bar(colour = "black",  stat = "identity", position = position_dodge(width = .8), size = .3, width = 0.8)+
  geom_errorbar(aes(x = SNP.Code.Name, ymin = FDR5_lb, ymax = FDR5_ub) , width = .2, position = position_dodge(.8))+
  facet_grid(~ Heritability, scales = "free_x", space = "free_x") +
  xlab("QTN") +
  labs(title = "Maize (N = 500) Epistasis Setting") +
  theme(plot.title = element_text(hjust = 0.50)) +
  ylab("% of Replications with a True Positive") +
  ylim(0, 1.00) + 
  #  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") + 
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue"))+
  theme(panel.border = element_rect(fill = NA, linetype = "dashed")) +
  theme_gray(base_size = 12) + theme( panel.border = element_rect(fill = NA, colour = "black"))+
  theme(strip.background =  element_rect(fill = "white", colour = "white"))+
  theme(axis.text.x = element_text(colour = "black",family = "Times", size = 12), axis.ticks = element_blank(), axis.ticks.length = unit(0.1, "cm")) +
  theme(axis.title.y = element_text(vjust = 1.3,colour = "black",family = "Times", size = 12), axis.text.y = element_text(  colour = "black", family = "Times"), axis.ticks = element_blank())+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
dev.off()
##

##Maize N = 500
tiff(file = "True Positive Detection Rates-Maize (N = 500) Epistasis.tiff",
     width = 641, height = 376, res = 100)
ggplot() +
  geom_bar(data = zm_500_Epistasis_rates, aes(x = SNP.Code.Name, y = FDR_5_TP, fill = Test), stat = "identity", position = position_dodge()) +
  ylim(0, 1.00) + 
  facet_grid(~ Heritability, scales = "free_x", space = "free_x") +
  labs(title = "Maize (N = 500) Epistasis Settings") +
  xlab("QTN") +
  ylab("% of Replications with a True Positive") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue", "MLM" = "magenta4"))
dev.off()

###FULL
zm_full_Epistasis_rates <- TnF_rates %>%
  select(Test, Species, Trait, Sample.Size, SNP, SNP.Code.Name, Heritability, Effect.Size, MAF, Window.Size, BP_5_FP, BP_5_TP, BP_10_FP, BP_10_TP, FDR_5_FP, FDR_5_TP, FDR_10_FP, FDR_10_TP) %>%
  filter(Trait == "Epistasis" & Sample.Size == "FULL" & Species == "Maize" & Effect.Size == 0.75)

zm_full_Epistasis_rates$Heritability <- gsub(0.30, "H2 = 0.30", zm_full_Epistasis_rates$Heritability)
zm_full_Epistasis_rates$Heritability <- gsub(0.80, "H2 = 0.80", zm_full_Epistasis_rates$Heritability)

for (i in 1:nrow(zm_full_Epistasis_rates)) {
  the.results <- exactci(x = zm_full_Epistasis_rates[i, 16] * 100, n = 100, conf.level = 0.95)
  zm_full_Epistasis_rates[i, 19] <- the.results$conf.int[1]
  zm_full_Epistasis_rates[i, 20] <- the.results$conf.int[2]
  rm(the.results)
}

colnames(zm_full_Epistasis_rates)[19:20] <- c("FDR5_lb", "FDR5_ub")

######
tiff(file = "Maize_N2815_Epistasis_True_Positives.tiff",
     width = 800, height = 400, res = 100)
ggplot(zm_full_Epistasis_rates, aes(x = SNP.Code.Name, y = FDR_5_TP, fill = Test)) +
  geom_bar(colour = "black",  stat = "identity", position = position_dodge(width = .8), size = .3, width = 0.8)+
  geom_errorbar(aes(x = SNP.Code.Name, ymin = FDR5_lb, ymax = FDR5_ub) , width = .2, position = position_dodge(.8))+
  facet_grid(~ Heritability, scales = "free_x", space = "free_x") +
  xlab("QTN") +
  labs(title = "Maize (N = 2815) Epistasis Setting") +
  theme(plot.title = element_text(hjust = 0.50)) +
  ylab("% of Replications with a True Positive") +
  ylim(0, 1.00) + 
  #  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") + 
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue"))+
  theme(panel.border = element_rect(fill = NA, linetype = "dashed")) +
  theme_gray(base_size = 12) + theme( panel.border = element_rect(fill = NA, colour = "black"))+
  theme(strip.background =  element_rect(fill = "white", colour = "white"))+
  theme(axis.text.x = element_text(colour = "black",family = "Times", size = 12), axis.ticks = element_blank(), axis.ticks.length = unit(0.1, "cm")) +
  theme(axis.title.y = element_text(vjust = 1.3,colour = "black",family = "Times", size = 12), axis.text.y = element_text(  colour = "black", family = "Times"), axis.ticks = element_blank())+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
dev.off()
######

##Maize N = 2532
tiff(file = "True Positive Detection Rates-Maize (N = 2532) Epistasis.tiff",
     width = 641, height = 376, res = 100)
ggplot() +
  geom_bar(data = zm_full_Epistasis_rates, aes(x = SNP.Code.Name, y = FDR_5_TP, fill = Test), stat = "identity", position = position_dodge()) +
  ylim(0, 1.00) + 
  facet_grid(~ Heritability, scales = "free_x", space = "free_x") +
  labs(title = "Maize FULL (N = 2532) Epistasis Settings") +
  xlab("QTN") +
  ylab("% of Replications with a True Positive") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue", "MLM" = "magenta4"))
dev.off()

###This next section is for Molike
TnF_rates$MAF <- gsub(0.10, "MAF = 0.10", TnF_rates$MAF)
TnF_rates$MAF <- gsub(0.40, "MAF = 0.40", TnF_rates$MAF)

TnF_rates$Heritability <- gsub(0.33, "h2 = 0.33", TnF_rates$Heritability)
TnF_rates$Heritability <- gsub(0.63, "h2 = 0.63", TnF_rates$Heritability)

molike_zm <- TnF_rates%>%
  select(Test, Species, Trait, Sample.Size, SNP, SNP.Code.Name, Heritability, Effect.Size, MAF, Window.Size, BP_5_FP, BP_5_TP, BP_10_FP, BP_10_TP, FDR_5_FP, FDR_5_TP, FDR_10_FP, FDR_10_TP) %>%
  filter(Species == "Maize" & Trait == "Molike")

molike_zm$Sample.Size <- gsub("FULL", "N = 2815", molike_zm$Sample.Size)
molike_zm$Sample.Size <- gsub("500", "N = 500", molike_zm$Sample.Size)

molike_zm$Effect.Size <- as.factor(molike_zm$Effect.Size)
molike_zm$Sample.Size <- factor(molike_zm$Sample.Size, levels = c("N = 500", "N = 2815"))

for (i in 1:nrow(molike_zm)) {
  the.results <- exactci(x = molike_zm[i, 16] * 100, n = 100, conf.level = 0.95)
  molike_zm[i, 19] <- the.results$conf.int[1]
  molike_zm[i, 20] <- the.results$conf.int[2]
  rm(the.results)
}

colnames(molike_zm)[19:20] <- c("FDR5_lb", "FDR5_ub")


###
tiff(file = "Maize_N2815_Epistasis_True_Positives.tiff",
     width = 800, height = 400, res = 100)
ggplot(molike_zm, aes(x = Effect.Size, y = FDR_5_TP, fill = Test)) +
  geom_bar(colour = "black",  stat = "identity", position = position_dodge(width = .8), size = .3, width = 0.8)+
  geom_errorbar(aes(x = Effect.Size, ymin = FDR5_lb, ymax = FDR5_ub) , width = .2, position = position_dodge(.8))+
  facet_grid(MAF ~ Heritability + Sample.Size, scales = "free_x", space = "free_x") +
  xlab("QTN") +
  labs(title = "Maize Molike Setting") +
  theme(plot.title = element_text(hjust = 0.50)) +
  ylab("% of Replications with a True Positive") +
  ylim(0, 1.00) + 
  #  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") + 
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue"))+
  theme(panel.border = element_rect(fill = NA, linetype = "dashed")) +
  theme_gray(base_size = 12) + theme( panel.border = element_rect(fill = NA, colour = "black"))+
  theme(strip.background =  element_rect(fill = "white", colour = "white"))+
  theme(axis.text.x = element_text(colour = "black",family = "Times", size = 12), axis.ticks = element_blank(), axis.ticks.length = unit(0.1, "cm")) +
  theme(axis.title.y = element_text(vjust = 1.3,colour = "black",family = "Times", size = 12), axis.text.y = element_text(  colour = "black", family = "Times"), axis.ticks = element_blank())+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
dev.off()
###


tiff(file = "True Positive Detection Rates-Maize Molike.tiff",
     width = 641, height = 376, res = 100)
ggplot() +
  geom_bar(data = molike_zm, aes(x = Effect.Size, y = FDR_5_TP, fill = Test), stat = "identity", position = position_dodge()) +
  ylim(0, 1.00) + 
  facet_grid(MAF ~ Heritability + Sample.Size, scales = "free_x", space = "free_x") +
  labs(title = "Maize") +
  xlab("vQTN Effect Size") +
  ylab("% of Replications with a True Positive") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue", "MLM" = "magenta4"))
dev.off()

###Arabidopsis Molike
molike_at <- TnF_rates_master %>%
  select(Test, Species, Trait, Sample.Size, SNP, SNP.Code.Name, Heritability, Effect.Size, MAF, Window.Size, BP_5_FP, BP_5_TP, BP_10_FP, BP_10_TP, FDR_5_FP, FDR_5_TP, FDR_10_FP, FDR_10_TP) %>%
  filter(Species == "Arabidopsis" & Trait == "Molike")

molike_at$Sample.Size <- gsub("FULL", "N = 1087", molike_at$Sample.Size)
molike_at$Sample.Size <- gsub("500", "N = 500", molike_at$Sample.Size)

molike_at$Effect.Size <- as.factor(molike_at$Effect.Size)
molike_at$Sample.Size <- factor(molike_at$Sample.Size, levels = c("N = 500", "N = 1087"))

molike_at$Heritability <- gsub("H2 = 0.33", "h2 = 0.33", molike_at$Heritability)
molike_at$Heritability <- gsub("H2 = 0.63", "h2 = 0.63", molike_at$Heritability)

molike_at$MAF <- gsub(0.1, "MAF = 0.10", molike_at$MAF)
molike_at$MAF <- gsub(0.4, "MAF = 0.40", molike_at$MAF)

tiff(file = "True Positive Detection Rates-Arabidopsis Molike.tiff",
     width = 641, height = 376, res = 100)
ggplot() +
  geom_bar(data = molike_at, aes(x = Effect.Size, y = FDR_5_TP, fill = Test), stat = "identity", position = position_dodge()) +
  ylim(0, 1.00) + 
  facet_grid(MAF ~ Heritability + Sample.Size, scales = "free_x", space = "free_x") +
  labs(title = "Arabidopsis") +
  xlab("vQTN Effect Size") +
  ylab("% of Replications with a True Positive") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue", "MLM" = "magenta4"))
dev.off()

#################This next section is for NULL settings only
NULL_rates <- TnF_rates_master %>%
  select(Test, Species, Trait, Sample.Size, SNP, SNP.Code.Name, Heritability, Effect.Size, MAF, Window.Size, BP_5_FP, BP_5_TP, BP_10_FP, BP_10_TP, FDR_5_FP, FDR_5_TP, FDR_10_FP, FDR_10_TP) %>%
  filter(Trait == "Null")

NULL_rates$FDR5_lb <- NA
NULL_rates$FDR5_ub <- NA
  
nomlm <- grep("MLM", NULL_rates$Test)
NULL_rates <- NULL_rates[-c(nomlm), ] 



for (i in 1:nrow(NULL_rates)) {
  the.results <- exactci(x = NULL_rates[i, 15] * 100, n = 100, conf.level = 0.95)
  NULL_rates[i, 19] <- the.results$conf.int[1]
  NULL_rates[i, 20] <- the.results$conf.int[2]
  rm(the.results)
}


the.results <- exactci(x = 7, n = 100, conf.level = 0.95)

the.lower.bound <- the.results$conf.int[1]

the.upper.bound <- the.results$conf.int[2]

####NULL
NULL_rates$Sample.Size <- gsub("N = 500", "500", NULL_rates$Sample.Size)

NULL_rates[5:6, 4] <- "N = 1087"
NULL_rates[7:8, 4] <- "N = 2532"

grep("Arabidopsis", NULL_rates$Species)
grep("Maize", NULL_rates$Species)

Null_rates_at <- NULL_rates[c(grep("Arabidopsis", NULL_rates$Species)), ]

Null_rates_zm <- NULL_rates[c(grep("Maize", NULL_rates$Species)), ]

tiff(file = "False_Positive_Detection_Rates.tiff",
     width = 641, height = 376, res = 100)
ggplot(NULL_rates, aes(x = Sample.Size, y = FDR_5_FP, fill = Test)) +
  geom_bar(colour = "black",  stat = "identity", position = position_dodge(width = .8), size = .3, width = 0.8)+
  geom_errorbar(aes(x = Sample.Size, ymin = FDR5_lb, ymax = FDR5_ub) , width = .2, position = position_dodge(.8))+
  facet_grid(~ Species, scales = "free_x", space = "free_x") +
  xlab("Sample Size") +
  ylab("% of Replications with a False Positive") +
  ylim(0, 0.25) + 
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") + 
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue"))+
  theme(panel.border = element_rect(fill = NA, linetype = "dashed")) +
  theme_gray(base_size = 12) + theme( panel.border = element_rect(fill = NA, colour = "black"))+
  theme(strip.background =  element_rect(fill = "white", colour = "white"))+
  theme(axis.text.x = element_text(colour = "black",family = "Times", size = 12), axis.ticks = element_blank(), axis.ticks.length = unit(0.1, "cm")) +
  theme(axis.title.y = element_text(vjust = 1.3,colour = "black",family = "Times", size = 12), axis.text.y = element_text(  colour = "black", family = "Times"), axis.ticks = element_blank())+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
dev.off()

NULL_rates$Sample.Size <- factor(c("N = 1087", "N = 2532", "N = 500"),
                levels = c("N = 500", "N = 1087", "N = 2532"))