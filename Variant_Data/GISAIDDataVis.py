#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 14:18:36 2022

@author: ryanweightman
"""


import pandas as pd

giSaidDF = pd.read_csv('hcov_north-america.tsv', sep='\t', header=0)
dfUS = giSaidDF[giSaidDF['country'] == 'USA']
dfNJ = giSaidDF[giSaidDF['division'] == 'New Jersey']
