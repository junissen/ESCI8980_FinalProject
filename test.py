#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Run command for Age Calculation function, with accompanying input functions. 

Created by Julia Nissen, 2017.

"""

import agecalc

inquiry_input = raw_input("Would you like to print as you go? [y/n] : ")

agecalc.Age_Calculation('72U_011517.xlsx', '72Th_011517.xlsx', '72U_wash_011517.xlsx', '72Th_wash_011517.xlsx', 
                '71U_chemblank.xlsx', '71Th_chemblank.xlsx', '71U_chemblank_wash.xlsx', '71Th_chemblank_wash.xlsx', 
                6E-7, 0.12040, 0.11730, 0.0026,
                'DIII-B', inquiry_input, 2017)



