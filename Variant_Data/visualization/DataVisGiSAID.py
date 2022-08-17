#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 09:52:38 2022

@author: ryanweightman
"""
import pandas as pd
import requests
import io
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

# This file brings in GISAID data and visualizes it. 

giSaidDF = pd.read_csv('ALASKA_GISAID_DATA.csv', header=0)

####Create list of variants And sum them for bar chart 

listvariants = giSaidDF['virusType']
frequency = {}

for item in giSaidDF['virusType']:
   if item in frequency:
      # incrementing the counr
      frequency[item] += 1
   else:
      # initializing the count
      frequency[item] = 1
variants= list(frequency.keys() )  

sumvar = [None]*len(variants)
for i in range(len(variants)):
    
     sumvar[i] = len(giSaidDF.virusType[giSaidDF['virusType'] == variants[i]] )

####Create list of variants And sum them for bar chart 

# plt.bar(variants, sumvar, align='center', alpha=0.5)
# plt.xticks(variants)
# plt.ylabel('Cases Reorted')
# plt.title('Covid Cases by Variant')

# plt.show()
# plt.close()



####Sum by variant by date for sequences over time chart 
date = {}
for item in giSaidDF['date']:
   if item in date:
      # incrementing the counr
      date[item] += 1
   else:
      # initializing the count
      date[item] = 1
dates= list(date.keys() )  
dates.sort()
sumdates = [None]*len(dates)
test = giSaidDF[giSaidDF['virusType'] == variants[1]] 

#print(test.date[test.date == dates[2]].count())

####Sum by variant by date for sequences over time chart  



####Plot by date
for i in range(len(variants)):
    variantSpecificDF = giSaidDF[giSaidDF['virusType'] == variants[i]] 
    sumdates = [None]*len(dates)
    for j in range(len(dates)):
        if dates[j] in list(variantSpecificDF['date']):
            
            sumdates[j] = variantSpecificDF.date[variantSpecificDF.date == dates[j]].count()
        else:
            sumdates[j] =0
    plt.plot(dates,sumdates,label = variants[i])
    plt.xticks(dates[::20], rotation=70)
plt.legend()
plt.show()
plt.close()
    
        
    
####Plot by date
    
    
    
giSaidDFomicron = giSaidDF[giSaidDF['virusType'] =='omicron']
    





# listvariants = giSaidDF['pangolin_lineage']
# frequency = {}

# for item in giSaidDFomicron['pangolin_lineage']:
#    if item in frequency:
#       # incrementing the counr
#       frequency[item] += 1
#    else:
#       # initializing the count
#       frequency[item] = 1
# variants= list(frequency.keys() )  

# sumvar = [None]*len(variants)
# for i in range(len(variants)):
    
#      sumvar[i] = len(giSaidDFomicron.virusType[giSaidDFomicron['pangolin_lineage'] == variants[i]] )

# ####Create list of variants And sum them for bar chart 

# plt.bar(variants, sumvar, align='center', alpha=0.5)
# plt.xticks(variants)
# plt.ylabel('Cases Reported')
# plt.title('Covid Cases by Variant')

# plt.show()
# plt.close()


#organizing Omicron by subvariant



subvariant = []
for item in giSaidDFomicron['pangolin_lineage']:
    if item[ 0 : 4 ] == 'BA.1':
        subvariant = subvariant+['BA.1']
    elif item[ 0 : 4 ] == 'BA.2':
        subvariant = subvariant+['BA.2']
    elif item[ 0 : 4 ] == 'BA.3':
        subvariant = subvariant+['BA.3']
    elif item[ 0 : 4 ] == 'BA.4':
        subvariant = subvariant+['BA.4']
    elif item[ 0 : 4 ] == 'BA.5':
        subvariant = subvariant+['BA.5']
    else:
        subvariant = subvariant+['other']
giSaidDFomicron['subvariant'] = subvariant   



date = {}
for item in giSaidDFomicron['date']:
   if item in date:
      # incrementing the counr
      date[item] += 1
   else:
      # initializing the count
      date[item] = 1
dates= list(date.keys() ) 
dates.sort() 
sumdates = [None]*len(dates)
test = giSaidDFomicron[giSaidDFomicron['virusType'] == variants[1]] 



#       pangolin_lineage
subvaraintList = ['BA.1','BA.2','BA.3','BA.4','BA.5', 'other']

for i in range(len(subvaraintList)):
    variantSpecificDF = giSaidDFomicron[giSaidDFomicron['subvariant'] == subvaraintList[i]] 
    sumdates = [None]*len(dates)
    for j in range(len(dates)):
        if dates[j] in list(variantSpecificDF['date']):
            
            sumdates[j] = variantSpecificDF.date[variantSpecificDF.date == dates[j]].count()
        else:
            sumdates[j] =0
            
    plt.plot(dates,sumdates,label = subvaraintList[i])
    plt.xticks(dates[::20], rotation=70)
plt.legend()
plt.show()
plt.close()    
    