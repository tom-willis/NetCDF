# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 17:17:44 2021

@author: geotw
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

##same as NetCDF_Clipper

#key standard format variables
format_str = "%Y%m%d"
#nc_file_format = r'3B42_Daily.%Y%m%d.7.nc4.nc4'
#input file
trmm_dir = r"D:\BarotseFloodplain\CATCHMAL\Data\ClimateData\ISIMIP\Q Data\orchidee_ipsl_rcp60"
output_dir = r"D:\BarotseFloodplain\CATCHMAL\Analysis\Python\ISIMIP\Q\orchidee_ipsl_rcp60"
#clip_domain = r"F:\TRMM\UpperZambezi_DrainageBasin_WGS_V0-01.shp"
#clip_domain = r"D:\BarotseFloodplain\CATCHMAL\Data\SWAT\BFP_SWAT_Area_WGS_v2-00.shp"
coords = [-13.48,23.06] #v1 
coords = [-14.5,23.5]

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
    bname = str(nc_file.split('_')[1])
    cname = str(nc_file.split('_')[3])
    dname = str(aname+"_"+bname+"_"+cname)
    return(dname)
    
def array_clip(coords,lat,lon):
    latli = np.argmin( np.abs( lat - coords[0] ) )
    # longitude lower and upper index
    lonli = np.argmin( np.abs( lon - coords[1] ) )
    return (latli,lonli)
    
def TRMM_yearoutput(startyear,endyear,trmm_dir,nc_file,lat,lon,name):
    nc_open = Dataset(os.path.join(trmm_dir,nc_file))
    dayfiles = days_list(nc_open)
  ##yearfiles = file_years(start_year,end_year,trmm_dir,trmm_files)
    #xmin,xmax,ymin,ymax = clip_coords(clip_domain,lat,lon)
    latmin,lonmin = array_clip(coords,lat,lon)
    tab = np.zeros([len(dayfiles)])
    for j,days in enumerate(dayfiles):
        tab[j] = nc_open.variables['dis'][j][latmin,lonmin]
    dict_date_x =  [x.strftime("%d/%m/%Y") for x in dayfiles]
    outdf = pd.DataFrame(data=tab,index = dict_date_x)
    outdf.columns=["Flow"]
    outdf.to_csv(os.path.join(output_dir,"Zambezi_Q_"+name+"_"+str(startyear)+"_"+str(endyear)+"_Data_v1-00.csv")
    ,sep=',',index=True, index_label = "Date")
    return ()

def main():
    nc_files = [f for f in os.listdir(trmm_dir) if f.endswith('.nc4')]
    lat,lon = latlonlist(trmm_dir,nc_files)
    for i,nc_file in enumerate(nc_files):
        years_str_1 = nc_file.split('_')[-1]
        years_str_2 = years_str_1.split('.')[0]
        endyear = years_str_2.split('-')[0]
        startyear = nc_file.split('_')[-2]
        name = output_name(nc_file)
        print(nc_file)
        TRMM_yearoutput(startyear,endyear,trmm_dir,nc_file,lat,lon,name)
    print ("Processing Complete -  Time to Execute  %s:" % (time.time()
    - start_time))
    return ()

if __name__== '__main__':
    main()