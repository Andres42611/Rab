# Rab: A Matrix-Based Genomic Ratio Analysis Program
'Rab' is a robust matrix-based genomic program designed for Python. Its primary function is to compute the ratio of derived allele frequencies between two populations by parsing and extracting data from genomic VCF files. This tool was developed by Andres del Castillo, an Undergraduate Researcher at the [Szpiech Lab](http://szpiech.com/index.html).

## Description
Rab is tailored to efficiently process large-scale genomic VCF files, supplemented with data delineating sample IDs for two populations (A & B) as well as chromosome and position numbers for two specific set of sites (1 & 2). The program performs comprehensive computations to derive allelic frequencies in population A relative to population B over the set of sites 1 and 2. Furthermore, Rab provides visualization capabilities, plotting jackknife and/or bootstrap distributions of the statistical point value as per user specifications.

This advanced computation is pivotal in contrasting genomic variances between two distinct populations. A practical application is evident in the comparative analysis of LoF (Loss of Function) mutational variants across genomes and genetic alterations leading to missense protein-coding sequences among populations. The inherent statistic offers insights into human demographic history, the expanding field of evolutionary medicine, and the intricate relationship between genomic structure and environmental influences.

The underlying formula driving the statistic has its origins in the scientific publication titled, "Supplementary Materials for 'Mountain gorilla genomes reveal the impact of long-term population decline and inbreeding'" by Yali Xue et al., 2015. For an in-depth understanding, the article can be accessed [here](https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.aaa3952&file=xue.sm.pdf).

Moreover, the program has been featured in the research publication "Genomic consequences of isolation and inbreeding in an island dingo population" by Ana V. Leon-Apodaca et al., 2023. For more insights into the real application of Rab program in research, the article can accessed [here](https://www.biorxiv.org/content/10.1101/2023.09.15.557950v1).


## Prerequisite Libraries:
* For computation of just the statistic, only argparse, scikit-allele, numpy, and pandas is needed:
`import os, sys, numpy, allel, argparse, pandas as pd, numpy.linalg as la`

* For visualization of bootstrapping and/or jackknife sampling, SciPy Stats and Matplotlib are needed in addition to the former bullet's modules:
`import matplotlib.pyplot as plt, random, scipy.stats as st; from scipy.stats import norm, t`
  

## Command Line Arguments:
Rab is integrated with a command-line interface, facilitating seamless interaction for its users. Below are the arguments to execute the program:
### Mandatory Inputs:

* -v, --vcf: Reference genomic VCF file (unzipped).
* -1, --site1: File detailing site 1 with chromosome in column 1 and position in column 2.
* -2, --site2: File detailing site 2 with chromosome in column 1 and position in column 2.
* -A, --popA: Subpopulation A file (sample ID in column 1, subpopulation ID in column 2).
* -B, --popB: Subpopulation B file (sample ID in column 1, subpopulation ID in column 2).

### Optional Parameters:

* -m, --mag: Logarithmic magnitude of R(A/B) v. R(B/A). Default: 'n'.
* -f, --per: Percentage removal during jackknife sampling. Default: 20%.
* -n, --iter: Number of sampling iterations. Default: 99.
* -c, --conf: Confidence interval for data distribution. Default: 95%.
* -o, --out: Output type. Options include 'jack', 'boot', 'pv', 'allg', 'pvj', 'pvb', 'all'; Default: 'all'. 

  #### Kinship Matrix Inputs:

  * -Akm, --Akinmatrix: Path to .king file (kinship matrix for popA).
  * -Bkm, --Bkinmatrix: Path to .king file (kinship matrix for popB).
  * -Akid, --Akinid: Path to .king file containing IDs (popA).
  * -Bkid, --Bkinid: Path to .king file containing IDs (popB).
