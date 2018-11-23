#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 15 16:57:04 2018

@author: beccalove
"""

import allel
import argparse
import numpy as np
import pkgutil
import sys


##read in the predictive SNPs for the inversion of interest

##extract the predictive SNPs from the supplied callset

##calculate the average genotype across the predictive SNPs

##return the average genotype and the # of sites across which it was calculated

inversionDict = {"2La" : ("2L","targets/targets_2La.txt"),
                 "2Rj" : ("2R","targets/targets_2Rj.txt"),
                 "2Rb" : ("2R","targets/targets_2Rb.txt"),
                 "2Rc" : ("2R","targets/targets_2Rc.txt"),
                 "2Rd" : ("2Rd","targets/targets_2Rd.txt"),
                 "2Ru" : ("2Ru","targets/targets_2Ru.txt"),
                 "test" : ("fooY","targets/test_targets_clean.txt")}

def parse_args(custom_args=None):
    
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="path to variant call format file \
                containing the genotypes")
    parser.add_argument("inversion", help="inversion to be classified",
                choices=["2La","2Rj","2Rb","2Rc","2Rd","2Ru","test"])
    parser.add_argument("-o", "--out",
                        help="name of the results file")
    parser.add_argument("-s", "--samples", help="samples to include", 
                        nargs='+')
    
    args = parser.parse_args(custom_args)

    return args

def import_data(callset_path):
    
    callset = allel.read_vcf(callset_path,
    fields = ['samples','calldata/GT','variants/CHROM',
              'variants/FILTER', 'variants/POS',
              'variants/REF','variants/ALT'])
    
    return callset

def import_inversion(inversion):
    
    path = inversionDict[inversion][1]
    
    targets_raw = pkgutil.get_data("compkaryo",path)

    targets = np.array([int(entry) 
    for entry 
    in targets_raw.decode().splitlines() 
    if not entry.startswith('#')])    
    
    return targets

def extract_vtbl_indices(targets, callset, chrom=None):
    
    if not chrom:
        chrom = inversionDict[args.inversion][0]
    
    indices = []

    for site in targets:
    
        where = np.where( (callset["variants/POS"] == site) &\
                         (callset["variants/CHROM"] == chrom))
    
        #print(where)
    
        if len(where[0]) > 0:
        
            #print(where[0][0])
        
            indices.append(where[0][0])
            
    return indices

def create_samples_bool(callset):
    
    if len(args.samples) == 1 and args.samples[0].endswith('.txt'):

        samples_file_handle = args.samples[0]
            
        samples = np.genfromtxt(samples_file_handle)
    
    else:
        
        samples = np.array(args.samples)
    
    samples_bool =\
    np.array([sample in callset["samples"] for sample in samples])
    
    return samples_bool
            
'''def extract_biallelic_SNPs(callset, indices):
    
    bi_bool =\
    allel.GenotypeArray(callset["calldata/GT"][indices]).count_alleles().\
    max_allele() <= 1
    
    return bi_bool'''
    
def calculate_genotype_at_concordant_sites(callset, indices, 
                                           samples_bool=None):
    
    genos = allel.GenotypeArray(callset["calldata/GT"][indices])
    
    if samples_bool is not None:
        
        genos = genos.subset(sel1 = samples_bool)
    
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
    
    samples_bool = None
    
    if args.samples:
        
        samples_bool = create_samples_bool(callset)
    
    #bi_bool = extract_biallelic_SNPs(callset, indices)
    
    av_gts, total_sites = calculate_genotype_at_concordant_sites(
            callset, indices, samples_bool)
    
    if not len(av_gts) == len(total_sites):
        
        raise ValueError("mean genotypes and total sites differ in length")
    
    if args.out is None:
        
        out = sys.stdout
        
    else:
        
        out = open(args.out, 'w')
    
    try:
                
        for i in range(len(av_gts)):
                            
            record = (str(av_gts[i]),str(total_sites[i]))
            
            out_record = "\t".join(record) + "\n"
                            
            out.write(out_record)
            
    finally:
        
        if args.out is not None:
            
            out.close()

if __name__ == "__main__":
    
    args = parse_args()
    main()