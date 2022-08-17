#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 09:46:29 2022

@author: ryanweightman
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 09:06:29 2022

@author: ryanweightman
"""


import pandas as pd
import requests
import io
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
    
# Downloading the csv file from your GitHub account

url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv" # Make sure the url is the raw version of the file on GitHub"
download = requests.get(url).content

# Reading the downloaded content and turning it into a pandas dataframe

df = pd.read_csv(io.StringIO(download.decode('utf-8')))

# Printing out the first 5 rows of the dataframe

#print (df['State'] == 'NJ')

njdf  = df[df['Province_State'] == 'Hawaii']

df2 = njdf.iloc[:, 11::].sum(axis=0)
njdf2 = pd.DataFrame(df2)

njdf2.reset_index(inplace=True)

njdf2.columns = ['date', 'cumulativeCases']

njdf2['dailyCases'] = njdf2['cumulativeCases'].diff().fillna(njdf2['cumulativeCases'])
for i in range(len(njdf2['dailyCases'])):
    if njdf2.dailyCases[i] < 0:
        njdf2.dailyCases[i] = 0
dates = list(njdf2['date'])     
#dates.sort()


# njdf2.plot(x ='date', y='dailyCases', kind = 'line')
# plt.plot(njdf2['date'],njdf2['dailyCases'])
# plt.xticks(njdf2['date'][::20], rotation=70)
# plt.legend()
# plt.show()
# plt.close() 



plt.plot(dates,njdf2['dailyCases'])
plt.xticks(njdf2['date'][::35], rotation=70)
plt.legend()
plt.show()
plt.close() 



