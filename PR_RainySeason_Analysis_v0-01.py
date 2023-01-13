# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 09:55:58 2021

@author: geotw
"""

import numpy as np
import pandas as pd
from datetime import date, timedelta, datetime

file = r"D:\BarotseFloodplain\CATCHMAL\Data\ClimateData\ISIMIP\PR Data\UpperZambezi_Average_pr_GFDL-ESM2M_historical_19810101_20051231_Data_v1-00.csv"
outfolder = r"D:\BarotseFloodplain\CATCHMAL\Analysis\Python\ISIMIP\PR"

##get the GCM name to be used in the output file name
name_a = file.split("\\")[-1]
name_b = name_a.split("_")[3]
name_c = name_a.split("_")[4]
name = str(name_b+"_"+name_c)

#open the file and extract the precipitation data
lines = pd.read_csv(file)
x = lines['PR']
dates = pd.to_datetime(lines['Date'])
startdate = dates[0].strftime("%Y-%m-%d")
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
dfdata['PR'] = x
#take the years to create the main looping variable
years = pd.DatetimeIndex(dfdata['Dates']).year.unique()

##create outputs variables for store & write out
tab = []
prtab = np.zeros([53,len(years)])
cumtab = np.zeros([53, len(years)])
cumperctab = np.zeros([53,len(years)])

#create the output file named
outname = str(outfolder+"\\WetSeasons_"+name+"_"+str(years.min())+"_"+str(years.max())+".csv")
prname = str(outfolder+"\\PR_"+name+"_"+str(years.min())+"_"+str(years.max())+".csv")
cumname = str(outfolder+"\\CumulativePR_"+name+"_"+str(years.min())+"_"+str(years.max())+".csv")
cumpername = str(outfolder+"\\CumulativePercentPR_"+name+"_"+str(years.min())+"_"+str(years.max())+".csv")


for i, year in enumerate(years[:-1]):
    startyear = str(years[i])
    endyear = str(years[i+1])
    print(startyear)
    
    startdate = (startyear+"-09-01")
    enddate = (endyear+"-08-31")
    
    dfwet = dfdata.set_index(['Dates'])
    dfloc = dfwet[startdate:enddate]
    ## resample to weeks, then remove any days below a threshold
    dfloc = dfloc.resample('W', level=0).sum()
    dfloc[dfloc.PR<5] = 0
    dfweeks = pd.DataFrame(dfloc)
    #dfweeks[dfweeks.PR<2] = 0
    
    ## Identify groups of wet weeks - select largest and use as the wet season
    ##https://stackoverflow.com/questions/26911851/how-to-use-pandas-to-find-consecutive-same-data-in-time-series
    dfweeks['value_zero'] = (dfweeks.PR.diff(1) == 0).astype('int').cumsum()
    dfweeks['Dates'] = dfweeks.index
    
    dfcountzero = pd.DataFrame({'BeginDate' : dfweeks.groupby('value_zero').Dates.first(), 
                  'EndDate' : dfweeks.groupby('value_zero').Dates.last(),
                  'Consecutive' : dfweeks.groupby('value_zero').size(), 
                  'No' : dfweeks.groupby('value_zero').PR.first()}).reset_index(drop=True)
    
    startweek = dfcountzero.BeginDate[(dfcountzero['Consecutive'].argmax())].strftime("%Y-%m-%d")
    endweek = dfcountzero.EndDate[(dfcountzero['Consecutive'].argmax())].strftime("%Y-%m-%d")
    dfwetweeks = dfweeks[startweek:endweek]
    dfwetweeks = dfwetweeks[:-1]
    ##Output the rainy season weeks and cumulative amount 
    pr =  list(dfweeks.PR)
    cumpr =  (list(dfweeks.PR.cumsum()))
    cumpr_perc = list(100*dfweeks.PR.cumsum()/dfweeks['PR'].sum())
    for j in range(0,53):
        prtab[j,i] = pr[j]
        cumtab[j,i] = cumpr[j]
        cumperctab[j,i] = cumpr_perc[j]

    #create and process outputs for main results
    start = dfwetweeks.Dates.min()
    end = dfwetweeks.Dates.max()
    start = start.tz_convert('GMT')
    startstr = start.strftime(format = '%d-%m')
    end = end.tz_convert('GMT')
    endstr = end.strftime(format = '%d-%m')
    length = abs((end - start).days)/7
    rainsum = dfwetweeks.PR.sum()
    tab.append([endyear,startstr,endstr,length,rainsum])

#output the wet season length
outdf = pd.DataFrame(tab,columns=['Year','Start','End','Length','sum'])
outdf.to_csv(outname,index=False)
##CUMULATIVE STEPS
#output the raw pr
outpr = pd.DataFrame(prtab,index = dict_date_x)
outpr.columns = [years[x] for x in outpr.columns]
outpr.to_csv(prname,index=True)

#output the cumulative pr
outcum = pd.DataFrame(cumtab,index = dict_date_x)
outcum.columns = [years[x] for x in outcum.columns]
outcum.to_csv(cumname,index=True)

#output the cumulative percentoutcum = pd.DataFrame(cumtab,index = dict_date_x)
outcumperc = pd.DataFrame(cumperctab,index = dict_date_x)
outcumperc.columns = [years[x] for x in outcumperc.columns]
outcumperc.to_csv(cumpername,index=True)