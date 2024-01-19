# Rab: Ratio of Mutational Load Between Subpopulations 'a' and 'b'
## Synopsis
'Rab' is a statistical genomic program designed with a primary function of computing the ratio of derived allele frequencies between two populations by parsing and extracting data from genomic VCF files, while also accounting kinship relationships using the Best Linear Unbiased Estimaton (BLUE) (McPeek et al. (2004)). This tool was developed by Andres del Castillo, an Undergraduate Researcher at the [Szpiech Lab](http://szpiech.com/index.html).

The underlying formula, named the 'RXY statistic', can be found in the supplementary materials for the research article "Mountain gorilla genomes reveal the impact of long-term population decline and inbreeding'" by Yali Xue et al., 2015, accessesible at [here](https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.aaa3952&file=xue.sm.pdf).

Use of the program can be found in the research publication "Genomic consequences of isolation and inbreeding in an island dingo population" by Ana V. Leon-Apodaca et al., 2023. For more insights into the real application of Rab program in research, the article can accessed [here](https://www.biorxiv.org/content/10.1101/2023.09.15.557950v1).


## Requirements:
* python 3.6.8
* scikit-allel 1.3.7
* argparse
* numpy
* scipy
* matplotlib

## Command Line Arguments:
### Mandatory Inputs:

* -v, --vcf: Reference genomic VCF file (unzipped).
* -1, --site1: File detailing site 1 with chromosome in column 1 and position in column 2 (tab-delimited).
* -2, --site2: File detailing site 2 with chromosome in column 1 and position in column 2 (tab-delimited).
* -A, --popA: Subpopulation A file (sample ID in column 1, subpopulation ID in column 2) (tab-delimited).
* -B, --popB: Subpopulation B file (sample ID in column 1, subpopulation ID in column 2) (tab-delimited).

### Optional Parameters:

* -f, --per: Percentage removal during jackknife sampling; Default: 20%.
* -n, --iter: Number of sampling iterations; Default: 99.
* -c, --conf: Confidence interval for data distribution; Default: 95%.
* -o, --out: Output type. Options include 'jack', 'boot', 'pv', 'allg', 'pvj', 'pvb', 'all'; Default: 'all'. 

  #### Kinship Parameters:

  * -Akm, --Akinmatrix: Path to .king file (kinship matrix for popA).
  * -Bkm, --Bkinmatrix: Path to .king file (kinship matrix for popB).
  * -Akid, --Akinid: Path to .king.id file containing IDs (popA).
  * -Bkid, --Bkinid: Path to .king.id file containing IDs (popB).

## Quick Example of Use:
After downloading the RAB-1.0.zip, ensuring the prerequsities are installed and up-to-date, and unzipping in the selected directory, if (for example) we wanted to have a output of the RXY statistic point value and a jackknife sampling plot with a 95% confidence interval of a 20% reduction over 1,000 iterations, while also taking into account kinship relationships, we use the command:<br><br>
```./Rab.py -v ./TESTVCF.vcf -1 ./SITES1.txt -2 ./SITES2.txt -A POPA.txt -B POPB.txt -Akm ./AKINSHIP.king -Bkm ./BKINSHIP.king -Akid ./AKINID.king.id -Bkid ./BKINID.king.id -f 20 -n 1000 -c 95 -o pvj```<br><br>

The specifics of these example input files can be found in the zip file under the folder 'Example'.
