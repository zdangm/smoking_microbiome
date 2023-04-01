

###gut_to_smoking_all.r
This script uses the "TwoSampleMR" (version 0.5.6) and "MRcML" (version 0.0.0.9) packages to perform a two-sample Mendelian Randomization analysis. The intent of these codes is to cyclically investigate the possible causal effect of each of the 211 microbial taxa on the five smoking phenotypes, using five MR functions: inverse variance weighted method, weighted median method, MR Egger regression, MR-PRESSO method, and cML-MA method. This script is implemented in R (version 4.1.2). This code is being shared to facilitate replication of our analysis, provide insight into the methods used, and build on our work. Feel free to use this code as a reference for your own Mendelian randomization analyses, and don't hesitate to contact me with any questions or comments.

###smoking_to_gut_all.r
This script uses exactly the same MR methods and it cyclically evaluates the potential causal effects of five smoking phenotypes on each of the 211 microbial taxa.

###gut_to_metabolites_to_smoking_MVMR.r
This script mainly employs the mr_mvivw and mr_mvegger functions in "MendelianRandomization" (version 0.7.0) packages to conduct a Multivariable Mendelian randomization (MVMR) analysis. The goal of this analysis is to uncover potential vertical pleiotropic mechanisms linking gut microbiota to smoking, which may arise from specific neurotransmitter-related or bacterial metabolites. This script is implemented in R (version 4.1.2).
