# Rab
R(A,B) Statistic Computation of mutational load between two population genomes

The program itself is able to take in a large-scale genomic .vcf file along with separate files outlining the sample IDs in population A and population B and the chromosome and position numbers of sites 1 and 2. The program then computes summations of the derived allelic frequencies in population A versus population B over sites 1 versus sites 2. 

This computation can then be applied to comparing genomic differences between two populations. For example, an application is seen via comparison of LoF (loss of function) mutational variants across genomes and genetic variants giving rise to missense protein coding sequences between populations. 

This current program draft has been completed in under a month and is still in the process of being finalized. Next steps will be to supplement the program with a possible machine learning model that will conduct bootstrapping sampling and logistic regression to supply a confidence interval to the point value returned. 

The formula for the statistic was taken from the scientific article "Supplementary Materials for 'Mountain gorilla genomes reveal the impact of long-term population decline and inbreeding'" by Yali Xue et. al 2015, page 17. Further can be read here: https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.aaa3952&file=xue.sm.pdf
