#!/usr/bin/env python3

import sys
import numpy
import matplotlib.pyplot as plt
import random
from scipy.stats import norm
import allel
import argparse
import fileinput
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category = RuntimeWarning)

parser = argparse.ArgumentParser(description = 'Calculate the R(A,B) ratio')
parser.add_argument('-v', '--vcf', type = str, metavar = '', required = True, help = 'reference genomic vcf file' )
parser.add_argument('-1', '--site1', type = str, metavar = '', required = True, help = 'sites 1 file with chr# in column 1 and pos # in column 2' )
parser.add_argument('-2', '--site2', type = str, metavar = '', required = True, help = 'sites 2 file with chr# in column 1 and pos # in column 2' )
parser.add_argument('-A', '--popA', type = str, metavar = '', required = True, help = 'subpopulation A file with sample ID in column 1 and subpop ID(A) in column 2' )
parser.add_argument('-B', '--popB', type = str, metavar = '', required = True, help = 'subpopulation B file with sample ID in column 1 and subpop ID(B) in column 2' )
parser.add_argument('-m', '--mag', type = str, nargs = '?', const = 'y', metavar = '', required = False, help = 'logarithmic magnitude of R(A/B) v. R(B/A)')
parser.add_argument('-f', '--per', type = int, nargs = '?', const = 20, metavar = '', required = False, help = 'percentage to be removed during jackknife sampling; default = 20')
parser.add_argument('-n', '--iter', type = int, nargs = '?', const = 9999, metavar = '', required = False, help = 'number of iterations for sampling; default = 9,999')
parser.add_argument('-o', '--out', type = str, nargs = '?', const = 'all', metavar = '', required = False,  help = ("output type (lowercase):" 
                                                                                                                    "'jack' - jackknife plot \n" 
                                                                                                                    "'boot' - bootstrap plot \n" 
                                                                                                                    "'pv' - R(A/B) point value \n" 
                                                                                                                    "'allg' - all plots \n" 
                                                                                                                    "'pvj' - point value and jackknife plot \n"
                                                                                                                    "'pvb' - point value and bootstrap plot \n"
                                                                                                                    "'all' - (default) all outputs "))
                                                                                                                

args = parser.parse_args()

if args.per is None:
    args.per = 20
if args.iter is None:
    args.iter = 9999
if args.out is None:
    args.out = 'all'


# STEP 1 - INITIALIZE REFERENCE DOCUMENT DATASET
def reffile(file):
    try:
        with open(file) as f:
            return f.name
    
    except OSError as error:
        raise SystemExit(error)
    
    except FileNotFoundError as e:
        raise SystemExit(e)

arg1 = reffile(args.vcf)



# STEP 2 - INITIALIZE CHROMOSOME TO EXTRACT DATA FROM
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



# STEP 3 - INITIALIZE SITES 1 AND SITES 2
def rfile_tolist(file, refdict, lst):
    try:
        if file == args.site1 or file == args.site2:
            with open(file) as infile:
                for line in infile:
                    try:
                        pos_num = int(line.split()[1])
                        chrom_num = line.split()[0]
                        if chrom_num in refdict and pos_num in refdict[chrom_num]:
                            elem = pos_num
                            lst.append(elem)
                        else:
                            print("No valid input for chromosome number and/or position on: \n'{0}' ".format(line))
                            sys.exit()
                   
                    except IndexError as i:
                        print("Check file format for file '{0}' matches specified input format".format(file.name))
                        raise SystemExit(i)
        else:
            with open(file) as infile:
                for line in infile:
                    try:
                        sample_ID = line.split()[0]
                        pop_ID = line.split()[1]
                        if sample_ID in allel.read_vcf(arg1)['samples'] and pop_ID in 'AB':
                            elem = sample_ID
                            lst.append(elem)
                    
                        else:
                            print("No valid input for sample ID and/or pop ID on: \n'{0}' ".format(line))
                            sys.exit()
                    except IndexError as i:
                        print("Check file format for file '{0}' matches specified input format".format(file.name))
                        raise SystemExit(i)
        
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



# STEP 7 - CALCULATE RATIO POINT-VALUE
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
    try:
        R_AB = (sum_X1/sum_X2)*(sum_Y2/sum_Y1)
        
        if pd.isna(R_AB) == False:
            return R_AB 
    
    except ZeroDivisionError:
        raise SystemExit(ZeroDivisionError)

sum_A1 = sum_af(Aaf_s1, Baf_s1, sites_1)
sum_A2 = sum_af(Aaf_s2, Baf_s2, sites_2)

sum_B1 = sum_af(Baf_s1, Aaf_s1, sites_1)
sum_B2 = sum_af(Baf_s2, Aaf_s2, sites_2)

pv = 'The value of R_(A,B) is: ' + str(round(R_AB(sum_A1, sum_A2, sum_B1, sum_B2), 4))

if args.mag == 'y':
    if __name__ == '__main__':
        R_AB = R_AB(sum_A1, sum_A2, sum_B1, sum_B2)
        R_BA = 1/((R_AB))
        mag_AB = round(numpy.log(R_AB), 4)
        mag_BA = round(numpy.log(R_BA), 4)
        print('The magnitude of R(A/B) to R(B/A) is {0} to {1}'.format(mag_AB, mag_BA))



# STEP 8 - DATA SAMPLING 
if args.out != 'pv':
  
    def stats(site1_lst, site2_lst, primpop1, secpop1, primpop2, secpop2, F, N, boot):
        lens1 = len(site1_lst)
        lens2 = len(site2_lst)
    
        y = round(lens1 * (int(F)/100))
        z = round(lens2 * (int(F)/100))
        max_Fval = int(((lens1 - 0.5)/(lens1)) * 100)
        min_Fval = int(100 * ((0.5)/(lens1)))
    
        if 1 <= y < lens1 and 1 <= z < lens2:
            idx1 = int(y)
            idx2 = int(z)
    
            Rab_est = []
        
            if not boot:
    
                for i in range(int(N)):
                    rand_lst1 = random.sample(site1_lst, (lens1 - idx1))
                    sum_A1 = sum_af(primpop1, secpop1, rand_lst1)
                    sum_B1 = sum_af(secpop1, primpop1, rand_lst1)
                
                    rand_lst2 = random.sample(site2_lst, (lens2 - idx2))
                    sum_A2 = sum_af(primpop2, secpop2, rand_lst2)
                    sum_B2 = sum_af(secpop2, primpop2, rand_lst2)
        
                    R_AB_est = R_AB(sum_A1, sum_A2, sum_B1, sum_B2)
        
                    if pd.isna(R_AB_est) == False:
                        Rab_est.append(R_AB(sum_A1, sum_A2, sum_B1, sum_B2))
                        
            else:
            
                for i in range(N):
                    rand_lst1 = random.choices(site1_lst, k = (lens1 - idx1))
                    rand_lst2 = random.choices(site2_lst, k = (lens2 - idx2))
                    
                    sum_A1 = sum_af(primpop1, secpop1, rand_lst1)
                    sum_B1 = sum_af(secpop1, primpop1, rand_lst1)
                    sum_A2 = sum_af(primpop2, secpop2, rand_lst2)
                    sum_B2 = sum_af(secpop2, primpop2, rand_lst2)

                    R_AB_est = R_AB(sum_A1, sum_A2, sum_B1, sum_B2)                    

                    if pd.isna(R_AB_est) == False:
                        Rab_est.append(R_AB(sum_A1, sum_A2, sum_B1, sum_B2))
                    
            return numpy.array(Rab_est)
    
        elif y < 1 or z < 1:
            print('F value of {0} is too small ... the value of F needs to be greater than {1} or a larger site dataset must be provided.'.format(F, min_Fval))
            sys.exit()
    
        else:
            print('F value of {0} is too large ... the value of F needs to be less than {1} or a larger site dataset must be provided.'.format(F, max_Fval))
            sys.exit()



# STEP 9 - PLOT DISTRIBUTIONS
if args.out == 'pvj' or args.out == 'jack' or args.out == 'all' or args.out == 'allg':
    jx = Rab_jvar = stats(sites_1, sites_2, Aaf_s1, Baf_s1, Aaf_s2, Baf_s2, args.per, args.iter, False)

    Jack_mean = numpy.mean(jx)
    Jack_fin_mean = round(Jack_mean, 4)

    Jack_sd = numpy.std(jx)
    Jack_fin_sd = round(Jack_sd, 4)

    jy = freq = norm.pdf(jx, Jack_mean, Jack_sd)

    plt.style.use('seaborn')
    plt.figure(figsize = (6, 6))
    plt.plot(jx, jy, color = 'green',
             linestyle = 'dashed')
    plt.fill_between(jx, jy, 0,
                     alpha = 0.2, color = 'green')

    plt.title('Jackknife Partial Estimates of R(A/B)')
    plt.xlabel('R(A/B) Value')
    plt.ylabel('Probability Density')

    jack = plt.show()


if args.out == 'pvb' or args.out == 'boot' or args.out == 'all' or args.out == 'allg':
    bx = Rab_bvar = stats(sites_1, sites_2, Aaf_s1, Baf_s1, Aaf_s2, Baf_s2, args.per, args.iter, True)

    Boot_mean = numpy.mean(bx)
    Boot_fin_mean = round(Boot_mean, 4)

    Boot_sd = numpy.std(bx)
    Boot_fin_sd = round(Boot_sd, 4)

    by = freq = norm.pdf(bx, Boot_mean, Boot_sd)

    plt.style.use('seaborn')
    plt.figure(figsize = (6, 6))
    plt.plot(bx, by, color = 'green',
             linestyle = 'dashed')
    plt.fill_between(bx, by, 0,
                     alpha = 0.2, color = 'green')

    plt.title('Bootstrapping Partial Estimates of R(A/B)')
    plt.xlabel('R(A/B) Value')
    plt.ylabel('Probability Density')

    boot = plt.show()



# STEP 10 - OUTPUT
if __name__ == '__main__':
    if args.out == 'pv':
        print(pv)
    elif args.out == 'jack':
        print('The mean R(A,B) value is {0} with a standard deviation of {1} via jackknife sampling.'.format(Jack_fin_mean, Jack_fin_sd))
    elif args.out == 'boot':
        print('The mean R(A,B) value is {0} with a standard deviation of {1} via bootstrap sampling.'.format(Boot_fin_mean, Boot_fin_sd))
    elif args.out == 'pvj':
        print(pv)
        print('The mean R(A,B) value is {0} with a standard deviation of {1} via jackknife sampling.'.format(Jack_fin_mean, Jack_fin_sd))
    elif args.out == 'pvb':
        print(pv)
        print('The mean R(A,B) value is {0} with a standard deviation of {1} via bootstrap sampling.'.format(Boot_fin_mean, Boot_fin_sd))
    elif args.out == 'all':
        print(pv)
        print('The mean R(A,B) value is {0} with a standard deviation of {1} via jackknife sampling.'.format(Jack_fin_mean, Jack_fin_sd))
        print('The mean R(A,B) value is {0} with a standard deviation of {1} via bootstrap sampling.'.format(Boot_fin_mean, Boot_fin_sd))
    else:
        print("Invalid output type")
        system.exit()
