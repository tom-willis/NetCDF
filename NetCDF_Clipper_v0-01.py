# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 17:58:43 2018

@author: tomwillis
"""

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import ogr
import os
import copy
import time
from datetime import date, timedelta, datetime
import re
import math

start_time = time.time()

#key standard format variables
format_str = "%Y%m%d"
#nc_file_format = r'3B42_Daily.%Y%m%d.7.nc4.nc4'
#input file
trmm_dir = r"D:\BarotseFloodplain\CATCHMAL\Data\ClimateData\ISIMIP\PR Data"
#clip_domain = r"F:\TRMM\UpperZambezi_DrainageBasin_WGS_V0-01.shp"
clip_domain = r"D:\BarotseFloodplain\CATCHMAL\Data\SWAT\BFP_SWAT_Area_WGS_v2-00.shp"

# def yearlist(trmm_files): 
#     years = []
#     for i,j in enumerate(trmm_files):
#         date = j.split('.')[1]
#         file_date = datetime.datetime.strptime(date,format_str)
#         year = file_date.year
#         years.append(year)
#         years_unique = list(set(years))
#         list.sort(years_unique)
#         years = years_unique
#     return(years)

def days_list(nc_open): 
    daylist = []
    time_since = nc_open.variables['time'].units
    match = re.search(r'\d{4}-\d{1}-\d{1}', time_since)
    start = datetime.strptime(match.group(), '%Y-%m-%d').date()
    days = nc_open.variables['time'][:]
    #start = date(1850,1,1) #should be made call variable
    for i,day in enumerate(days):
        ##amend for ISIMIP - days provided in 32bit format,int required
        day = int(day)
        delta = timedelta(day)
        offset = start + delta
        daylist.append(offset)
    return(daylist)
    
def latlonlist(trmm_dir,trmm_files):
    infile = trmm_files[0]
    datain = Dataset(os.path.join(trmm_dir,infile))
    lat = datain.variables['lat'][:]
    ##sort descending
    lat[::-1].sort()
    lon = datain.variables['lon'][:]
    return(lat,lon)

# def file_years(start_year,end_year,trmm_dir,trmm_files):
#     datestart = str(start_year)+'1101'
#     dateend = str(end_year)+'1031'
#     start_date = datetime.datetime.strptime(datestart, format_str).date()
#     start_end = datetime.datetime.strptime(dateend, format_str).date()
#     delta_one_day = datetime.timedelta(days=1)
#     yearfiles = []
#     date = start_date
#     while date <= start_end:
#         file_data = date.strftime(nc_file_format)
#         datafile = os.path.join(trmm_dir,file_data)
#         if os.path.isfile(datafile):
#             yearfiles.append(datafile)
#         date += delta_one_day
#     return(yearfiles)

def clip_coords(clip_domain,lat,lon):
    xy_minmax = []
    shp = ogr.Open(clip_domain)
    lyr = shp.GetLayer()
    xy_minmax = lyr.GetExtent()
    xy_round = [math.ceil(y) if x == 1 else int(y) for x,y in enumerate(xy_minmax)]
    xmin = xy_round[0]
    xmax = xy_round[1]
    ymin = xy_round[2]
    ymax = xy_round[3]
    return(xmin,xmax,ymin,ymax)

def output_name(nc_file):
    ##name to use for writeoutfile - may need to change
    ##the number used from the string depending on the nc file name
    aname = str(nc_file.split('_')[0])
    bname = str(nc_file.split('_')[2])
    cname = str(nc_file.split('_')[3])
    dname = str(aname+"_"+bname+"_"+cname)
    return(dname)
    
def array_clip(clip_domain,lat,lon):
    shp = ogr.Open(clip_domain)
    lyr = shp.GetLayer()
    xmin,xmax,ymin,ymax = lyr.GetExtent()
    latbounds = [int(ymax),math.ceil(ymin)]##swapped added for southers
    lonbounds = [int(xmin),int(xmax)]
    latli = np.argmin( np.abs( lat - latbounds[0] ) )
    latui = np.argmin( np.abs( lat - latbounds[1] ) ) 
    # longitude lower and upper index
    lonli = np.argmin( np.abs( lon - lonbounds[0] ) )
    lonui = np.argmin( np.abs( lon - lonbounds[1] ) )
    return (latli,latui,lonli,lonui)
    
def TRMM_yearoutput(startyear,endyear,trmm_dir,nc_file,lat,lon,name):
    nc_open = Dataset(os.path.join(trmm_dir,nc_file))
    dayfiles = days_list(nc_open)
  ##yearfiles = file_years(start_year,end_year,trmm_dir,trmm_files)
    xmin,xmax,ymin,ymax = clip_coords(clip_domain,lat,lon)
    latmin,latmax,lonmin,lonmax = array_clip(clip_domain,lat,lon)
    if lonmin == 0: lonmax = lonmax + 1
    if latmin == 0: latmax = latmax + 1
    lat_clip = [y  for y in lat if (y > ymin) & (y < ymax)]
    lon_clip = [x  for x in lon if (x > xmin) & (x < xmax)]
    dict_of_precip={}
    for i,nc_file in enumerate(dayfiles):
        nc_date = dayfiles[i]
        precip = nc_open.variables['pr'][i][latmin:latmax, lonmin:lonmax]
        precip_rain = precip * 86400
        dict_of_precip[nc_date] = copy.deepcopy(precip_rain)
    
    tab = np.zeros([(len(lat_clip)*len(lon_clip)),len(dict_of_precip)+2])
    count = 0
    for i in range(0,len(lat_clip)):
        for j in range(0,len(lon_clip)):
            tab[count,0] = lat_clip[i]
            tab[count,1] = lon_clip[j]
            count = count + 1
            
    ##create a new tab to record the average values - use the current processes
    ##and calculate the averages when files are open, rather than reopen files
    av_tab = []
    ##take the dates from the dictionary to create a column variable
    keys = dict_of_precip.keys()
    for i,keyin in enumerate(keys):
        count = 0
        ##take all the data from individual tables in the key array to write out to tab
        precip_in = dict_of_precip[keyin]
        ##take the average value of the mask area 
        av_tab.append(precip_in.mean())
        for j,element in np.ndenumerate(precip_in):
            tab[count,i+2] = (element)
            count = count + 1

    dict_keys = list(dict_of_precip.keys())
    dict_date_x =  [x.strftime("%d/%m/%Y") for x in dict_keys]
    outdf = pd.DataFrame(data=tab)
    outdf.columns = [dict_date_x[x-2] if x > 1 else x for x in outdf.columns]
    outdf = outdf.rename(columns={0:'Latitude',1:'Longitude'})
    outdf.to_csv(os.path.join(trmm_dir,"UpperZambezi_Points_"+name+"_"+str(startyear)+"_"+str(endyear)+"_Data_v1-00.csv")
    ,sep=',',index=False)
    ###AVERAGE OUTPUT
    ## replace all mean values below 0.1 as zero
    av_tab = [0 if x<0.1 else x for x in av_tab]
    ###create datafram to output
    outav = pd.DataFrame(data = av_tab)
    outav['Date'] = dict_date_x
    outav = outav.rename(columns={0:'PR'})
    outav = outav[['Date','PR']]
    outav['Date'] = dict_date_x
    outav.to_csv(os.path.join(trmm_dir,"UpperZambezi_Average_"+name+"_"+str(startyear)+"_"+str(endyear)+"_Data_v1-00.csv")
    ,sep=',',index=False)
    return ()

def main():
    nc_files = [f for f in os.listdir(trmm_dir) if f.endswith('.nc4')]
    lat,lon = latlonlist(trmm_dir,nc_files)
    for i,nc_file in enumerate(nc_files):
        years_str_1 = nc_file.split('_')[-1]
        years_str_2 = years_str_1.split('.')[0]
        startyear = years_str_2.split('-')[0]
        endyear = years_str_2.split('-')[1]
        name = output_name(nc_file)
        print(nc_file)
        TRMM_yearoutput(startyear,endyear,trmm_dir,nc_file,lat,lon,name)
    print ("Processing Complete -  Time to Execute  %s:" % (time.time()
    - start_time))
    return ()

if __name__== '__main__':
    main()