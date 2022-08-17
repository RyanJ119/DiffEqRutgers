#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 09:06:29 2022

@author: ryanweightman
"""


import pandas as pd
import requests
import io
    
# Downloading the csv file from your GitHub account

url = "https://raw.githubusercontent.com/govex/COVID-19/master/data_tables/demographic_data/demographics_by_state_standardized.csv" # Make sure the url is the raw version of the file on GitHub"
download = requests.get(url).content

# Reading the downloaded content and turning it into a pandas dataframe

df = pd.read_csv(io.StringIO(download.decode('utf-8')))

# Printing out the first 5 rows of the dataframe

#print (df['State'] == 'NJ')

njdf  = df[df['State'] == 'NJ']

njdfCases = njdf[df['Category'] == 'Cases']
njVaccines = njdf[df['Category'] == 'Vaccines']
njdfDeaths= njdf[df['Category'] == 'Deaths']

#print (njdfCases)

