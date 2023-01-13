# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 09:55:58 2021

@author: geotw
"""

import numpy as np
import pandas as pd
from datetime import date, timedelta, datetime

file = r"D:\BarotseFloodplain\CATCHMAL\Analysis\Python\ISIMIP\Q\Zambezi_Q_orchidee_ewembi_rcp60_2041_2050_Data_v1-00.csv"
outfolder = r"D:\BarotseFloodplain\CATCHMAL\Analysis\Python\ISIMIP\Q"

##get the GCM name to be used in the output file name
name_a = file.split("\\")[-1]
name_b = name_a.split("_")[2]
name_c = name_a.split("_")[3]
name_d = name_a.split("_")[4]
name = str(name_b+"_"+name_c+"_"+name_d)

#open the file and extract the precipitation data
lines = pd.read_csv(file)
x = lines['Flow']
dates = pd.to_datetime(lines['Date'])
startdate = dates[0].strftime("%Y-%m-%d")
syear = dates[0].year
#create the dates and weeks using the ISO standard timestamps
dates = pd.date_range(startdate, periods=len(x), freq='D').tz_localize('GMT')
startweek = dates[244].strftime("%Y-%m-%d") 
weeks = list(pd.date_range(startweek, periods=53, freq='W').tz_localize('GMT'))
#convert to a string format to use as the index in the write out tables
#the first year is added to allow projection in excel
dict_date_x =  [week.strftime("%d/%m/%Y") for week in weeks]

##create the data frame with the main data
dfdata = pd.DataFrame()
dfdata['Dates'] = dates
dfdata['Flow'] = x
#take the years to create the main looping variable
years = pd.DatetimeIndex(dfdata['Dates']).year.unique()
## temporary fix 
years = years[1:]

##create outputs variables for store & write out
tab = []
prtab = np.zeros([53,len(years)])

#create the output file named
outname = str(outfolder+"\\RiverFlow_"+name+"_"+str(years.min())+"_"+str(years.max())+".csv")
flowname = str(outfolder+"\\Q_"+name+"_"+str(years.min())+"_"+str(years.max())+".csv")

for i, year in enumerate(years[:-1]):
    startyear = str(years[i])
    endyear = str(years[i+1])
    print(startyear)
    
    startdate = (startyear+"-09-01")
    enddate = (endyear+"-08-31")
    
    dfwet = dfdata.set_index(['Dates'])
    dfloc = dfwet[startdate:enddate]
    ## resample to weeks, then remove any days below a threshold
    dfloc = dfloc.resample('W', level=0).mean()
    dfweeks = pd.DataFrame(dfloc)
    # this line is moved to ensure that the full annual flow is taken into account, rather than capped values
    flow =  list(dfweeks.Flow)
    #cap the flow values to bankful
    dfweeks[dfweeks.Flow>500] = 500
    
    ## Identify groups of wet weeks - select largest and use as the wet season
    ##https://stackoverflow.com/questions/26911851/how-to-use-pandas-to-find-consecutive-same-data-in-time-series
    dfweeks['value_zero'] = (dfweeks.Flow.diff(1) != 0).astype('int').cumsum()
    dfweeks['Dates'] = dfweeks.index
    
    dfcountzero = pd.DataFrame({'BeginDate' : dfweeks.groupby('value_zero').Dates.first(), 
                  'EndDate' : dfweeks.groupby('value_zero').Dates.last(),
                  'Consecutive' : dfweeks.groupby('value_zero').size(), 
                  'No' : dfweeks.groupby('value_zero').Flow.first()}).reset_index(drop=True)
    
    startweek = dfcountzero.BeginDate[(dfcountzero['Consecutive'].argmax())].strftime("%Y-%m-%d")
    endweek = dfcountzero.EndDate[(dfcountzero['Consecutive'].argmax())].strftime("%Y-%m-%d")
    dfwetweeks = dfweeks[startweek:endweek]
    dfwetweeks = dfwetweeks[:-1]
    ##Output the rainy season weeks and cumulative amount 
    #cumpr =  (list(dfweeks.PR.cumsum()))
    #cumpr_perc = list(100*dfweeks.PR.cumsum()/dfweeks['PR'].sum())
    if (len(flow) == 53):
        for j in range(0,53):
            prtab[j,i] = flow[j]
    
        if (len(dfwetweeks) >0):#create and process outputs for main results
            start = dfwetweeks.Dates.min()
            end = dfwetweeks.Dates.max()
            start = start.tz_convert('GMT')
            startstr = start.strftime(format = '%d-%m')
            end = end.tz_convert('GMT')
            endstr = end.strftime(format = '%d-%m')
            length = abs((end - start).days)/7
            flowmax = max(flow)
        else:
            startstr = '25-12' 
            endstr = '01-01'
            length = 0
            flowmax = max(flow)
        tab.append([endyear,startstr,endstr,length,flowmax])
    else:
        startstr = '25-12' 
        endstr = '01-01'
        length = 0
        flowmax = max(flow)
        tab.append([endyear,startstr,endstr,length,flowmax])

#output the wet season length
outdf = pd.DataFrame(tab,columns=['Year','Start','End','Length','max'])
outdf.to_csv(outname,index=False)
##CUMULATIVE STEPS
#output the raw pr
outpr = pd.DataFrame(prtab,index = dict_date_x)
outpr.columns = [years[x] for x in outpr.columns]
outpr.to_csv(flowname,index=True)
