import pandas as pd
import requests
import io
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from math import nan, isnan


# Depending on state, certain dates must be removed, check line 140
#To change state, you must change state on line 28 and change csv on 29
##########JHU data pull and clean


# Downloading the csv file from your GitHub account

url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv" # Make sure the url is the raw version of the file on GitHub"
download = requests.get(url).content

# Reading the downloaded content and turning it into a pandas dataframe

df = pd.read_csv(io.StringIO(download.decode('utf-8')))

stateDF  = df[df['Province_State'] == 'New York']   

giSaidDF = pd.read_csv('NEW_YORK_GISAID_DATA.csv', header=0) #States should match! 




df2 = stateDF.iloc[:, 11::].sum(axis=0)
casesDF = pd.DataFrame(df2)
casesDF.reset_index(inplace=True)
casesDF.columns = ['date', 'cumulativeCases']

#Turn cumulative cases into daily cases  by using diff
casesDF['dailyCases'] = casesDF['cumulativeCases'].diff().fillna(casesDF['cumulativeCases'])
for i in range(len(casesDF['dailyCases'])):
    if casesDF.dailyCases[i] < 0:
        casesDF.dailyCases[i] = 0
dates = list(casesDF['date'])     

#Smooth daily cases using a 7 day rolling average
casesDF['Rolling_Avg']= casesDF['dailyCases'].rolling(7).mean()

##########JHU data pull and clean






##########giSaidDF data pull 



####Create list of variants 

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


####Create list of variants 



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
dates = [x for x in dates if str(x) != 'nan']
#print(dates)

dates.sort() #Comment for California
sumdates = [None]*len(dates)


variantsSummedDF = pd.DataFrame()  
variantsSummedDF['date'] = dates


####Plot variants by date
for i in range(len(variants)):
    variantSpecificDF = giSaidDF[giSaidDF['virusType'] == variants[i]] 
    sumdates = [None]*len(dates)
    for j in range(len(dates)):
        if dates[j] in list(variantSpecificDF['date']):
            
            sumdates[j] = variantSpecificDF.date[variantSpecificDF.date == dates[j]].count()
        else:
            sumdates[j] =0
    variantsSummedDF[variants[i]] = sumdates
    plt.plot(dates,sumdates,label = variants[i])
    plt.xticks(dates[::20], rotation=70)
plt.legend()
plt.show()
plt.close()
    
        
    
####Plot variants by date
    
    
            #### state speific date cleaning, follow comments
#dates.remove('2021')   #Uncomment for Missouri
dates.remove('2021-01') #Uncomment for New York and Florida
dates.remove('2021-02') #Uncomment for New York and Pennsylvania and Florida
#dates.remove('2021-03') #Uncomment for Pennsylvania and Florida
#dates.remove('2021-04') #Uncomment for Florida
#dates.remove('2021-05') #Uncomment for Pennsylvania and  Florida
#dates.remove('2021-06') #Uncomment for Pennsylvania and Florida
#dates.remove('2021-07') #Uncomment for Pennsylvania and Florida
#dates.remove('2021-08') #Uncomment for Florida
#dates.remove('2021-11') #Uncomment for Florida
#dates.remove('2021-12') #Uncomment for Florida
#dates.remove('2022-01') #Uncomment for Florida
#dates.remove('2022-03') #Uncomment for Florida
#dates.remove('2022-05') #Uncomment for Florida
#dates.remove('2021') #Uncomment for California
#dates.remove('2020') #Uncomment for Florida
#dates.remove('2020-12') #Uncomment for Florida
            #### state speific date cleaning, follow comments
            
#Fix dates into datetime objects for easy use, then expand variant DF to include all dates
for i in range(len(dates)):
    
    dates[i] = datetime.strptime(dates[i], '%Y-%m-%d')
    #dates[i] = dates[i].date


JHUcases = casesDF
 
JHUcases['date'] = JHUcases['date'].astype('datetime64[ns]')
variantsSummedDF['date'] = variantsSummedDF['date'].astype('datetime64[ns]')
#JHUcases['date'] = JHUcases['date'].dt.date
for j in variants:
    JHUcases[j] = 0

for i in JHUcases['date']:
    if i in dates:
        for j in variants:
            
            JHUcases.at[JHUcases[JHUcases['date'] == i].index[0], j]= variantsSummedDF.loc[variantsSummedDF[variantsSummedDF['date'] == i].index[0]][j]
           #print(JHUcases[JHUcases['date'] == i].index[0])
           #print(variantsSummedDF.loc[variantsSummedDF[variantsSummedDF['date'] == i].index[0]][j])


#Fix dates into datetime objects for easy use, then expand variant DF to include all dates







#Plot daily cses and rolling average 
    
plt.plot(JHUcases['date'],JHUcases['dailyCases'], label = 'dailyCases')
plt.plot(JHUcases['date'],JHUcases['Rolling_Avg'],label = 'Rolling_Avg')
plt.xticks(JHUcases['date'][::35], rotation=70)
plt.legend()
plt.show()
plt.close() 
#Plot daily cses and rolling average 


#Plot variants with daily cases and rolling average
plt.plot(JHUcases['date'],JHUcases['dailyCases']/50, label = 'dailyCases')
plt.plot(JHUcases['date'],JHUcases['Rolling_Avg']/50,label = 'Rolling_Avg')
for i in variants:
    plt.plot(JHUcases['date'],JHUcases[i],label = i)

plt.xticks(JHUcases['date'][::35], rotation=70)
plt.legend()
plt.show()
plt.close() 

#Plot variants with daily cases and rolling average


#Create normalized variant DF for percent plotting 

normalizedJHUCasesDF =  pd.DataFrame(JHUcases)
normalizedJHUCasesDF = normalizedJHUCasesDF.drop(columns=['cumulativeCases', 'date', 'dailyCases','Rolling_Avg' ], axis=1)
normalizedJHUCasesDF = normalizedJHUCasesDF.div(normalizedJHUCasesDF.sum(axis=1), axis=0) #normalize the rows

 
  
y = np.vstack(normalizedJHUCasesDF)
labels = normalizedJHUCasesDF.columns
fig, ax = plt.subplots()
ax.stackplot(JHUcases['date'].values, normalizedJHUCasesDF.T, labels = labels)  
plt.plot(JHUcases['date'],JHUcases['Rolling_Avg']/JHUcases['Rolling_Avg'].max(),label = 'Rolling_Avg')
ax.legend(loc ='upper left')
plt.xticks(JHUcases['date'][::35], rotation=70)

ax.set_title('Noramlized Variant Distribution')
plt.show()


#Create normalized variant DF for percent plotting 

listvariants = giSaidDF['pangolin_lineage']
frequency = {}


giSaidDFomicron = giSaidDF[giSaidDF['virusType'] =='omicron']
   
for item in giSaidDFomicron['pangolin_lineage']:
    if item in frequency:
      # incrementing the counr
      frequency[item] += 1
    else:
      # initializing the count
      frequency[item] = 1
variants= list(frequency.keys() )  

sumvar = [None]*len(variants)
for i in range(len(variants)):
    
      sumvar[i] = len(giSaidDFomicron.virusType[giSaidDFomicron['pangolin_lineage'] == variants[i]] )

####Create list of variants And sum them for bar chart 




#organizing Omicron by subvariant



subvariant = []
for item in giSaidDFomicron['pangolin_lineage']:
    item = str(item)
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
dates = [x for x in dates if str(x) != 'nan']
dates.sort() 
sumdates = [None]*len(dates)
 
subvariantsSummedDF = pd.DataFrame()  
subvariantsSummedDF['date'] = dates




#       pangolin_lineage
subvariantList = ['BA.1','BA.2','BA.3','BA.4','BA.5', 'other']

for i in range(len(subvariantList)):
    variantSpecificDF = giSaidDFomicron[giSaidDFomicron['subvariant'] == subvariantList[i]] 
    sumdates = [None]*len(dates)
    for j in range(len(dates)):
        if dates[j] in list(variantSpecificDF['date']):
            
            sumdates[j] = variantSpecificDF.date[variantSpecificDF.date == dates[j]].count()
        else:
            sumdates[j] =0
        subvariantsSummedDF[subvariantList[i]] = sumdates     
    plt.plot(dates,sumdates,label = subvariantList[i])
    plt.xticks(dates[::20], rotation=70)
plt.legend()
plt.show()
plt.close()   




for i in range(len(dates)):
    
    dates[i] = datetime.strptime(dates[i], '%Y-%m-%d')
    #dates[i] = dates[i].date


#JHUcases = casesDF
 
JHUcases['date'] = JHUcases['date'].astype('datetime64[ns]')
subvariantsSummedDF['date'] = subvariantsSummedDF['date'].astype('datetime64[ns]')
#JHUcases['date'] = JHUcases['date'].dt.date
for j in subvariantList:
    JHUcases[j] = 0

for i in JHUcases['date']:
    if i in dates:
        for j in subvariantList:
            
            JHUcases.at[JHUcases[JHUcases['date'] == i].index[0], j]= subvariantsSummedDF.loc[subvariantsSummedDF[subvariantsSummedDF['date'] == i].index[0]][j]
           #print(JHUcases[JHUcases['date'] == i].index[0])
           #print(variantsSummedDF.loc[variantsSummedDF[variantsSummedDF['date'] == i].index[0]][j])



del JHUcases['omicron']




#Create normalized variant DF for percent plotting 

normalizedJHUCasesDFOmicronSubVariants =  pd.DataFrame(JHUcases)
normalizedJHUCasesDFOmicronSubVariants = normalizedJHUCasesDFOmicronSubVariants.drop(columns=['cumulativeCases', 'date', 'dailyCases','Rolling_Avg' ], axis=1)
normalizedJHUCasesDFOmicronSubVariants = normalizedJHUCasesDFOmicronSubVariants.div(normalizedJHUCasesDFOmicronSubVariants.sum(axis=1), axis=0) #normalize the rows

 
  
y = np.vstack(normalizedJHUCasesDFOmicronSubVariants)
labels = normalizedJHUCasesDFOmicronSubVariants.columns
fig, ax = plt.subplots()
ax.stackplot(JHUcases['date'].values, normalizedJHUCasesDFOmicronSubVariants.T, labels = labels)  
plt.plot(JHUcases['date'],JHUcases['Rolling_Avg']/JHUcases['Rolling_Avg'].max(),color='black', label = 'Rolling_Avg')
ax.legend(loc ='upper left')
plt.xticks(JHUcases['date'][::35], rotation=70)

ax.set_title('Noramlized Variant Distribution')
plt.show()
