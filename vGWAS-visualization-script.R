#############################################################################################
################################### vGWAS Visualization Code ################################
#############################################################################################

### SPECIAL STATEMENT ###
#The result files that were used for the visualises were complied and reformatted in Excel, not R.

install.packages("ggplot2")
install.packages("ggthemes") # Install 

library("ggplot2")
library("tidyr")
library(gridExtra)
library("ggthemes") # Load

###Read in rate files
wd <- list()
wd$three_output <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_One/3_output/"
list.files(wd$three_output)

##For my graphic folder, I am going to create a setting working directory shortcut to 
##there
wd$three_output_figures <- paste0(wd$three_output, "Figs/")

#These two files are the traits that are not Oilllike
False_positives_wo_oil <- read.csv(paste0(wd$three_output_nonoil, "At_Zm_Oil_FALSE_positive_detection_rates.csv"), sep = ",")
True_positives_wo_oil <- read.csv(paste0(wd$three_output_nonoil, "At_Zm_Oil_True_positive_detection_rates.csv"), sep = ",")

False_positives_just_oil <- read.csv(paste0(wd$three_output, "False Positive Detection Rates-Just_Oil.csv"), sep = ",")
True_positives_just_oil <- read.csv(paste0(wd$three_output, "True Positive Detection Rate-Just_Oil.csv"), sep = ",")

############################################################
tiff(file = "False Positive Detection Rates-Just NULL.tiff",
     width = 600, height = 350, res = 100)
ggplot() +
  geom_bar(data=False_positive_wo_oil_longer, aes(x = Sample.size, y = FP.FDR.5, fill = Test), stat = "identity", position = position_dodge()) +
  facet_grid( ~Species) +
  xlab("Sample Size") +
  ylab("Proportion of Replications with a False Positive")
dev.off()

tiff(file="True Positive Detection Rates-Just Molike Setting.tiff",
     width = 641, height = 376, res = 100)
ggplot() +
  geom_bar(data = True_positive_detection_rates, aes(x = Sample.Size, y = FDR5, fill = Test), stat = "identity", position = position_dodge()) +
  facet_grid(Species ~ Trait) +
  labs(title = "Molybdenum-like Settings") +
  xlab("Sample Size") +
  ylab("Proportion of Replications with a True Positive") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual("Test", values = c("BFT" = "darkorange3", "DGLM" = "navyblue"))
dev.off()

tiff(file = "True Positive Detection Rates-Oil-like Settings.tiff",
     width = 750, height = 425, res = 100)
ggplot() +
  geom_bar(data = True_Positive_Detection_Rate_Just_Oil, aes(x = Effect.Size, y = FDR5, fill = Test), stat = "identity", position = position_dodge()) +
  facet_grid(Trait ~ Species + Sample.Size) +
  labs(title = "Oil-like Settings") +
  xlab("Effect Size") +
  ylab("Proportion of Replications with a True Positive") +
  theme(plot.title = element_text(hjust = 0.5)) 
dev.off()