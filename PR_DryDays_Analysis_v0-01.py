# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 09:55:58 2021

@author: geotw
"""

import numpy as np
import pandas as pd
from datetime import date, timedelta, datetime

file = r"D:\BarotseFloodplain\CATCHMAL\Data\ClimateData\ISIMIP\PR Data\UpperZambezi_Average_pr_GFDL-ESM2M_rcp85_20410101_20501231_Data_v1-00.csv"
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
dates = pd.date_range(startdate, periods=9131, freq='D').tz_localize('GMT')

##create the data frame
dfdata = pd.DataFrame()
dfdata['Dates'] = dates
dfdata['PR'] = x

## resample to weeks, then remove any days below a threshold
dfweeks = dfdata.resample('W', on='Dates')['PR'].sum()
dfweeks = pd.DataFrame(dfweeks)
dfweeks[dfweeks.PR<2] = 0

## Find consecutive days of dry need
##https://stackoverflow.com/questions/26911851/how-to-use-pandas-to-find-consecutive-same-data-in-time-series
dfweeks['value_grp'] = (dfweeks.PR.diff(1) != 0).astype('int').cumsum()
dfweeks['Dates'] = dfweeks.index
dfcount = pd.DataFrame({'BeginDate' : dfweeks.groupby('value_grp').Dates.first(), 
              'EndDate' : dfweeks.groupby('value_grp').Dates.last(),
              'Consecutive' : dfweeks.groupby('value_grp').size(), 
              'No' : dfweeks.groupby('value_grp').PR.first()}).reset_index(drop=True)
dfdryweeks = dfcount.loc[dfcount['No'] == 0]
dfdryweeks['Year'] = pd.DatetimeIndex(dfdryweeks['BeginDate']).year
years = pd.DatetimeIndex(dfdryweeks['BeginDate']).year.unique()

##emptylist to write out
tab = []
outname = str(outfolder+"\\DryYears_"+name+"_"+str(years.min())+"_"+str(years.max())+".csv")
##to be looped by years
for i, year, in enumerate(years):
    include = dfdryweeks[dfdryweeks['BeginDate'].dt.year == year]
    start = include.BeginDate.min()
    end = include.EndDate.max()
    start = start.tz_convert('GMT')
    startstr = start.strftime(format = '%d-%m-%Y')
    end = end.tz_convert('GMT')
    endstr = end.strftime(format = '%d-%m-%Y')
    length = abs((end - start).days)/7
    tab.append([year,startstr,endstr,length])

outdf = pd.DataFrame(tab,columns=['Year','Start','End','Length'])
outdf.to_csv(outname,index=False)