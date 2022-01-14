# vGWAS_arabidopis_maize
This is the repository for "Simulation study evaluating the ability of two statistical approaches to identify variance quantitative trait loci Arabidopsis and maize" 

The following is an outline/order of the scripts used for this publication with a description of each file.

1.) additive-effects-simplePHENOTYPES-script.R : This script was used to produce the null settings, epistasis, GxE, and additive effects for the Molydenum-like settings with additive effects

2.) vGWAS-vQTL-simulation-code-script.R : This script was used to simulate the variance genetic architectures for the Molydenum-like settings

3.) vGWAS-Brown Forsythe Test Script.R : This script was used to run the Brown-Forsythe Test (BFT) on the different settings

4.) GAPIT-script.R : This script used GAPIT to run the Unified Linear Mixed Model (MLM)

5.) DGLM-script.R : This script was used to run DGLM on the different settings

6.) Positive-detection-rate-script.R : This script was used to process the false and true-positive detection rates of all the settings

7.) vGWAS-visualization-script.R : This script was used to visualize the false and true-positive detection rates
