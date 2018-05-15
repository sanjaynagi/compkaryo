#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 15 17:13:26 2018

@author: beccalove
"""

import unittest
import compkaryo
import allel

class CalculateGenotypesTestCase(unittest.TestCase):
    
    def setUp(self):
        
        self.genotypes = allel.GenotypeArray([
              [[0,0],[1,1],[0,1],[0,1],[0,1],[0,1]],
              [[0,1],[2,2],[0,2],[0,2],[0,2],[0,2]],
              [[1,1],[0,0],[0,1],[0,1],[0,1],[0,1]],
              [[2,2],[0,1],[0,2],[1,2],[0,2],[1,2]],
              [[0,0],[1,1],[0,1],[0,1],[0,1],[0,1]],
              [[1,1],[0,0],[0,1],[0,1],[0,1],[0,1]],
              [[0,0],[1,1],[0,1],[0,1],[0,1],[0,1]],
              [[1,1],[1,1],[1,1],[1,1],[1,1],[1,1]],
              [[0,0],[1,2],[0,1],[0,2],[0,2],[0,2]],
              [[0,1],[0,1],[1,1],[0,0],[0,1],[0,0]]
              ], dtype='i1')
        
        self.name="2Rq"
        self.mock_inversion_dict = {}
        self.mock_inversion_dict[self.name] =\
        ingenos.Inversion(b'2R', 15,20,750,775)
        
    def test_complex_case(self):
        
        test_expression =\
        ingenos.construct_filter_expression(self.name,self.mock_inversion_dict,
                                            buffer=0)
        
        self.assertEqual(
                test_expression,
                '( ( (POS > 15) & (POS < 20) ) |\
        ( (POS > 750) & (POS < 775) ) )') 
    
    def test_simple_case(self):
        
        test_expression =\
        ingenos.construct_filter_expression(self.name,self.mock_inversion_dict,
                                            buffer=0,whole_inversion=True)
        
        self.assertEqual(test_expression,
                        '( (POS > 15) & (POS < 775) )')
        
    def test_negative_coordinates_fail(self):
        self.assertRaises(ValueError, ingenos.construct_filter_expression, 
                          self.name, self.mock_inversion_dict)

if __name__ == '__main__':
	unittest.main()
