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

#n = 6
##get the 
name_a = file.split("\\")[-1]
name = name_a.split("_")[3]

lines = pd.read_csv(file)
x = lines['PR']
dates = pd.to_datetime(lines['Date'])
startdate = dates[0].strftime("%Y-%m-%d")
dates = pd.date_range(startdate, periods=9131, freq='D').tz_localize('GMT')
#average_dates = dates.iloc[:-n]

#def moving_average(x, w):
#    return np.convolve(x, np.ones(w), 'valid') / w

#average_pr = moving_average(x,n+1)
#average_pr[average_pr<1] = 0

#df = pd.DataFrame()
#df['dates'] = average_dates
##df['dates'] = pd.to_datetime(df["dates"].dt.strftime('%Y-%m'))
#df['averagePR'] = average_pr

##weeks
##test for weeks
##df['weeks'] = pd.DatetimeIndex(dfweek['Dates']).week
dfweek = pd.DataFrame()
dfweek['Dates'] = dates
dfweek['PR'] = x
## try resample
dfweeks = dfweek.resample('W', on='Dates')['PR'].sum()
##from series to dataframe maybe integrate with the above 
dfweeks = pd.DataFrame(dfweeks)

## try groupby
##https://stackoverflow.com/questions/46562401/group-python-pandas-dataframe-per-weeks-starting-on-monday
#df1 = dfweek.groupby(pd.Grouper(freq='W', key='Dates'))['PR'].sum()
#df1[df1<2] = 0
##correct for small weeks
dfweeks[dfweeks.PR<2] = 0

dftest = dfweeks.sort_values(["Dates", "PR"])

## 1. testing mehtod for finding consecutive days in pandas
## see - https://stackoverflow.com/questions/59965641/finding-consecutive-days-in-the-pandas-dataframe

dbTest = r"D:\BarotseFloodplain\CATCHMAL\Analysis\Python\DBTest.csv"

dbtestlines = pd.read_csv(dbTest)
dbtestlines['ColB'] = pd.to_datetime(dbtestlines['ColB'])
##1
dbtestlines = dbtestlines.sort_values(["ColA", "ColB"])
dbtestlines["ColB_std"] = dbtestlines["ColB"] - pd.to_timedelta(range(dbtestlines.shape[0]), 'day')
s = dbtestlines.groupby(["ColA", "ColB_std"])["ColB_std"].count()
colA = s[ s >= 3 ].index.get_level_values(0).unique().values
dbtestlines[ dbtestlines["ColA"].isin(colA) ]
dbtestlines[ dbtestlines["ColA"].isin(colA) ]

##2 . second method for idetnfying the number of consecutive days @value
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

##emptylist to write to
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



