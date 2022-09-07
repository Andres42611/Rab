#!/usr/bin/env python3

import sys
import numpy
import allel
import argparse
import fileinput

parser = argparse.ArgumentParser(description = 'Calculate the R(A,B) ratio')
parser.add_argument('-r', '--ref', type = str, metavar = '', required = True, help = 'reference genomic vcf file' )
parser.add_argument('-1', '--site1', type = str, metavar = '', required = True, help = 'sites 1 file with chr# in column 1 and pos # in column 2' )
parser.add_argument('-2', '--site2', type = str, metavar = '', required = True, help = 'sites 2 file with chr# in column 1 and pos # in column 2' )
parser.add_argument('-A', '--popA', type = str, metavar = '', required = True, help = 'subpopulation A file with sample ID in column 1 and subpop ID(A) in column 2' )
parser.add_argument('-B', '--popB', type = str, metavar = '', required = True, help = 'subpopulation B file with sample ID in column 1 and subpop ID(B) in column 2' )
args = parser.parse_args()


# STEP 1 - INITIALIZE REFERENCE DOCUMENT DATASET
def reffile(file):
    try:
        with open(file) as f:
            return f.name
    except OSError as error:
        raise SystemExit(error)
    
    except FileNotFoundError as e:
        raise SystemExit(e)

arg1 = reffile(args.ref)


#STEP 2 - INITIALIZE CHROMOSOME TO EXTRACT DATA FROM
def chrom_dict(file): #make reference dict with all chr # of all sites as keys and an array of POS within that chromosome as the values
    var = allel.read_vcf(file)['variants/CHROM']
    var_lst = [*set(var.tolist())]
    
    chrom_dict = {}
    for i in range(1, (1 + len(var_lst))):
        chrom = 'chr' + str(i)
        callset = ((allel.read_vcf(arg1, region = chrom))['variants/POS'])
        chrom_dict[chrom] = callset
    
    return chrom_dict #will be used as the total reference for chromosome # matching in making the ref dicts for popA and pop B as well as for af arrays

refchrom_dict = chrom_dict(arg1)

# STEP 3 - INITIALIZE SITE 1 AND SITE 2 SITES
def rfile_tolist(file, refdict, lst):
    try:
        if file == args.site1 or file == args.site2:
            with open(file) as infile:
                for line in infile:
                    pos_num = int(line.split()[1])
                    chrom_num = line.split()[0]
                    if chrom_num in refdict and pos_num in refdict[chrom_num]:
                        elem = pos_num
                        lst.append(elem)
                    else:
                        print("No valid input for chromosome number and/or position on: \n'{0}' ".format(line))
                        sys.exit()
        else:
            with open(file) as infile:
                for line in infile:
                    sample_ID = line.split()[0]
                    pop_ID = line.split()[1]
                    if sample_ID in allel.read_vcf(arg1)['samples'] and pop_ID in 'AB':
                        elem = sample_ID
                        lst.append(elem)
                    else:
                        print("No valid input for sample ID and/or pop ID on: \n'{0}' ".format(line))
                        sys.exit()
        return lst

    
    except OSError as error:
        raise SystemExit(error)
    
    except FileNotFoundError as e:
        raise SystemExit(e)
        
sites_1 = []
sites_1_chr = []

sites_2 = []
sites_2_chr = []

rfile_tolist(args.site1, refchrom_dict, sites_1)
rfile_tolist(args.site2, refchrom_dict, sites_2)

def sites_pchrom(file, chrom_set): #make site-specific list with all unqiue chr # within the sites directly from site files
    with open(file) as infile:
        for line in infile:
            elem = line.split()[0]
            chrom_set.append(elem)
    res = []
    [res.append(elem) for elem in chrom_set if elem not in res]
    
    return res #will be used as the unique callset for site-specific chromosome number

site1_chrlist = sites_pchrom(args.site1, sites_1_chr) #unique chromosome list of site 1 chromosomes
site2_chrlist = sites_pchrom(args.site2, sites_2_chr) #unique chromosome list of site 2 chromosomes


# STEP 4 - INITIALIZE POP A AND POP B
popA = []

popB = []

rfile_tolist(args.popA, refchrom_dict, popA)
rfile_tolist(args.popB, refchrom_dict, popB)


# STEP 5 - MAKE ALLELE FREQUENCY ARRAYS FOR POP A AND POP B
def af_arr(chr_list, ref_vcf, subpop):
    af_arr = []
    for elem in chr_list:
        garr = allel.GenotypeArray(allel.read_vcf(ref_vcf, region = str(elem), samples = subpop)['calldata/GT'])
        garr_af = (garr.count_alleles()/((len(garr[0,0]))*(len(subpop))))[ : , 1: ]
        af_arr.append(garr_af)
        
    return af_arr

#list of arrays with allele vlaues for chromosomes in site
popA_af_s1 = af_arr(site1_chrlist, arg1, popA)
popA_af_s2 = af_arr(site2_chrlist, arg1, popA)
popB_af_s1 = af_arr(site1_chrlist, arg1, popB)
popB_af_s2 = af_arr(site2_chrlist, arg1, popB)


# STEP 6 - MAKE SITE-SPECIFIC DICTS FOR POP A AND POP B AT ALL SITES
def arrgen_todict(af_arr, chrom_list, chrom_dict):
    splitArrs = ([individualArray] for individualArray in af_arr) #CREATES A GENERATOR FOR LARGE ARRAYS
    pop_lst = []
    pop_dict= {}
    
    for i in splitArrs:
        indv_arr = numpy.array(i)
        raw_dict = dict(enumerate(indv_arr.flatten(), 1))
        pop_lst.append(raw_dict)

    for j in range(len(chrom_list)):
        elem = chrom_list[j]
        key_lst = chrom_dict[elem]
        final_dict = dict(zip(key_lst, list((pop_lst[j]).values())))
        pop_dict[elem] = final_dict
    
    return pop_dict
    
#dict with all af values for all sites within the chromosomes called specifically in sites
Aaf_s1 = arrgen_todict(popA_af_s1, site1_chrlist, refchrom_dict) #e.g. dict of A af values for all sites in site 1 chroms
Aaf_s2 = arrgen_todict(popA_af_s2, site2_chrlist, refchrom_dict)
Baf_s1 = arrgen_todict(popB_af_s1, site1_chrlist, refchrom_dict)
Baf_s2 = arrgen_todict(popB_af_s2, site2_chrlist, refchrom_dict)


# STEP 7 - CALCULATE RATIO
def sum_af(prim_af_sx, sec_af_sx, sites):
    sum_af = 0
    for j in range(len(sites)):
        for i in prim_af_sx:
            if sites[j] in prim_af_sx[i]:
                sum_af = sum_af + (prim_af_sx[i][(sites[j])] * (1 - sec_af_sx[i][(sites[j])]))
            else:
                continue
    return sum_af

def R_AB(sum_X1, sum_X2, sum_Y1, sum_Y2):
    R_AB = (sum_X1/sum_X2)*(sum_Y2/sum_Y1)
    return 'The value of R_(A,B) is: ' + str(R_AB)

sum_A1 = sum_af(Aaf_s1, Baf_s1, sites_1)
sum_A2 = sum_af(Aaf_s2, Baf_s2, sites_2)

sum_B1 = sum_af(Baf_s1, Aaf_s1, sites_1)
sum_B2 = sum_af(Baf_s2, Aaf_s2, sites_2)

if __name__ == '__main__':
    print(R_AB(sum_A1, sum_A2, sum_B1, sum_B2))
