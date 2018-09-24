#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 15 16:57:04 2018

@author: beccalove
"""

import numpy as np
import allel
from collections import namedtuple

Inversion = namedtuple('Inversion',['arm','start','end',"marker_file_path"])

inversionDict = {}

inversionDict["2La"] = Inversion(b'2L',20524058,42165532,"foo")
inversionDict["2Rb"] = Inversion(b'2R',19023925,26758676,"foo")
inversionDict["2Rc"] = Inversion(b'2R',26750000,31473100,"foo")
inversionDict["2Ru"] = Inversion(b'2R',31473000,35505236,"foo")
inversionDict["2Rj"] = Inversion(b'2R',3262186,15750717,"foo")
inversionDict["2Rd"] = Inversion(b'2R',31480000,42600000,"foo")

def construct_filter(inversion):
    
    chrom = inversionDict[inversion].arm
    start = inversionDict[inversion].proximal_start
    end = inversionDict[inversion].distal_end
    
    flt = '((POS >= {start}) & (POS <= {end}) & (CHROM == {chrom}))'.format(
            start=start, end=end, chrom=chrom)
    
    return flt

def import_data(path_to_vtbl_or_h5):
    
    pass

def import_inversion(inversion):
    
    path = inversionDict[inversion].path
    
    with open(path) as f:
        
        lines = f.read().splitlines()
        
    try:
        
        float_list = [float(x) for x in lines]
        
    except ValueError:
        
        raise
        
    if not all(x.is_integer() for x in float_list):
        
        raise ValueError("Non-integers present in concordance positions")
        
    concordant_pos_list = [int(x) for x in float_list]
    
    return concordant_pos_list
        
def filter_objects(flt, vtbl, gts):
    
    if not vtbl.shape[0] == gts.shape[0]:
        
        raise ValueError("Variant table and genotypes must be same length")
        
    bool_flt = vtbl.eval(flt)
    
    vtbl_chunk = vtbl[bool_flt]
    
    gts_chunk = gts[bool_flt]
    
    if not vtbl_chunk.shape[0] == gts_chunk.shape[0]:
        
        raise ValueError("Something has gone wrong in filtering")
    
    return(bool_flt, vtbl_chunk, gts_chunk)

def extract_vtbl_indices(concordant_pos_list, vtbl):
    
    pos_indices_list = []
        
    vtbl_pos_list = vtbl["POS"]
    
    for pos in concordant_pos_list:
        
        get_index = (np.where(vtbl_pos_list == pos))[0]
        
        if len(get_index) == 1:
            
            pos_indices_list.append(get_index[0])
            
    return pos_indices_list
    
def calculate_genotype_at_concordant_sites(pos_indices_list, gts):
    
    pos_indices = np.array(pos_indices_list)
    
    alt_counts = gts[[pos_indices]].to_n_alt()
    
    av_gts = np.sum(alt_counts, axis = 0) / len(pos_indices)
    
    return av_gts, len(pos_indices)

def assign_karyos(avg_gts):
    
    karyos = []
    
    for sample in avg_gts:
        
        if not 0 <= sample <= 2:
            
            raise ValueError("value out of bounds, must be 0 >= value >= 2")
    
        if sample < 0.66:

            karyos.append(0)

        elif sample >= 0.66 and sample < 1.33:

            karyos.append(1)

        elif sample >= 1.33:

            karyos.append(2)
            
    return karyos