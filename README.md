# Rab Population Genetics Tool
Rab is a independently-built matrix-based Python genomic algorithm able to measure the ratio of derived allele frequencies in two populations through the parsing and extraction of genomic data all through user input via the command line.

The program itself is able to take in a large-scale genomic .vcf files along with separate files outlining the sample IDs in population A and population B and the chromosome + position numbers of sites 1 and 2. The program then computes summations of the derived allelic frequencies in population A versus population B over sites 1 versus sites 2 and plots jackknife and/or bootstrap distributions of the statistical point value based on the output specified by the user.

This computation can then be applied to comparing genomic differences between two populations. For example, an application is seen via comparison of LoF (loss of function) mutational variants across genomes and genetic variants giving rise to missense protein coding sequences between populations. The statistic itself is able to applied to understand human demographic history, evolutionary medicine, and the intersection between genomic architecture and environmental factors.

The formula for the statistic was taken from the scientific article "Supplementary Materials for 'Mountain gorilla genomes reveal the impact of long-term population decline and inbreeding'" by Yali Xue et. al 2015, page 17. Further can be read here: https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.aaa3952&file=xue.sm.pdf
