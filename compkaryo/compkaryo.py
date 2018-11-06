#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 15 16:57:04 2018

@author: beccalove
"""

import allel
import argparse
import numpy as np
import sys

parser = argparse.ArgumentParser()
parser.add_argument("vcf", help="path to variant call format file \
                    containing the genotypes")
parser.add_argument("inversion", help="inversion to be classified",
                    choices=["2La","2Rj","2Rb","2Rc","2Rd","2Ru"])
parser.add_argument("-o", "--out", help="name of the results file")
args = parser.parse_args()

##read in the predictive SNPs for the inversion of interest

##extract the predictive SNPs from the supplied callset

##calculate the average genotype across the predictive SNPs

##return the average genotype and the number of sites across which it was calculated

inversionDict = {"2La" : ("2L","path/to/2La"),
                 "2Rj" : ("2R","path/to/2Rj"),
                 "2Rb" : ("2R","path/to/2Rb"),
                 "2Rc" : ("2R","path/to/2Rc"),
                 "2Rd" : ("2Rd","path/to/2Rd"),
                 "2Ru" : ("2Ru","path/to/2Ru")}


def import_data(callset_path):
    
    callset = allel.read_vcf(callset_path,
    fields = ['samples','calldata/GT','calldata/GQ','variants/CHROM',
              'variants/FILTER_PASS', 'variants/POS','variants/QUAL',
              'variants/REF','variants/ALT'])
    
    return callset

def import_inversion(inversion):
    
    path = inversionDict[inversion][1]
    
    concordant_pos_list = np.loadtxt(path)    
    
    return concordant_pos_list

def extract_vtbl_indices(concordant_pos_list, callset):
    
    chrom = inversionDict[args.inversion][0]
    
    indices = []

    for site in concordant_pos_list:
    
        where = np.where( (callset["variants/POS"] == site) &\
                         (callset["variants/CHROM"] == chrom))
    
        #print(where)
    
        if len(where[0]) > 0:
        
            #print(where[0][0])
        
            indices.append(where[0][0])
            
    return indices

def extract_biallelic_SNPs(callset, indices):
    
    bi_bool =\
    allel.GenotypeArray(callset["calldata/GT"][indices]).count_alleles().\
    is_biallelic()
    
    return bi_bool
    
def calculate_genotype_at_concordant_sites(callset, bi_bool, indices):
    
    genos = allel.GenotypeArray(callset["calldata/GT"][indices][bi_bool])
    
    alt_count = genos.to_n_alt()
    
    is_called = genos.is_called()
    
    av_gts = np.mean(np.ma.MaskedArray(
            alt_count, mask = ~is_called), axis=0).data
            
    total_sites = np.sum(is_called, axis=0)
    
    return av_gts, total_sites

def main():
    
    callset = import_data(args.vcf)
    
    target_list = import_inversion(args.inversion)
    
    indices = extract_vtbl_indices(target_list, callset)
    
    bi_bool = extract_biallelic_SNPs(callset, indices)
    
    av_gts, total_sites = calculate_genotype_at_concordant_sites(
            callset, bi_bool, indices)
    
    if not len(av_gts) == len(total_sites):
        
        raise ValueError("mean genotypes and total sites differ in length")
    
    if args.out is None:
        
        out = sys.stdout
        
    else:
        
        out = open(args.out, 'w')
        
    try:
        
        for i in range(len(av_gts)):
        
            out_record = (av_gts[i],"\t",total_sites[i],"\n")
                
            out.write(out_record)
            
    finally:
        
        if args.out is not None:
            
            out.close()
    
if __name__ == "__main":
    main()

