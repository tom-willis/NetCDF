# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 16:35:33 2023

Function to clip NetCDF files to a shapefile and return average values for area

Returns a CSV file with the daily average value for the climate variable

To change the climate variable that you want to look at, adjust line 92 to 
match the variable naming convention in the netcdf
(given on line 84 as netc_file.pr, where .pr refers to the Precipitation Rate)

This should be defined as a variable, but for the time being it is left as this

Requires 2 inputs - 
1) shapefile (Line 40) 
2) directory with netcdf files (specifically .nc4 files (Line41))

@author: geotw
"""

from netCDF4 import Dataset
import geopandas
import xarray
from shapely.geometry import mapping
import pandas as pd
import os
import re
import numpy as np
from datetime import timedelta, datetime

#global variable
filepath =  os.getcwd()

##https://gis.stackexchange.com/questions/354782/masking-netcdf-time-series-data-from-shapefile-using-python

##a variable that COULD be adjusted, depending on the CRS of the NC file
crs_code = "epsg:4326"

##INPUTS
netc_folder = r"C:\Users\geotw\OneDrive - University of Leeds\EnivroMAL\Data\BethMroz\Test"
catchshapefile = r"C:\Users\geotw\OneDrive - University of Leeds\EnivroMAL\Data\BethMroz\Test\BFP_SWAT_Area_WGS_v2-00.shp"

def days_list(netc_file): 
    daylist = []
    nc_open = Dataset(netc_file)
    time_since = nc_open.variables['time'].units
    match = re.search(r'\d{4}-\d{1}-\d{1}', time_since)
    start = datetime.strptime(match.group(), '%Y-%m-%d').date()
    days = nc_open.variables['time'][:]
    for i,day in enumerate(days):
        ##amend for ISIMIP - days provided in 32bit format,int required
        day = int(day)
        delta = timedelta(day)
        offset = start + delta
        daylist.append(offset)
    return(daylist)

def output_name(nc_file):
    ##name to use for outputs this is set to the ISIMIP naming convention
    aname = str(nc_file.split('_')[0])
    bname = str(nc_file.split('_')[2])
    cname = str(nc_file.split('_')[3])
    dname = str(aname+"_"+bname+"_"+cname)
    return(dname)

def file_name(nc_file):
    file = nc_file.split("\\")[-1]
    years_str_1 = nc_file.split('_')[-1]
    years_str_2 = years_str_1.split('.')[0]
    startyear = years_str_2.split('-')[0]
    endyear = years_str_2.split('-')[1]
    file_name = output_name(file)
    name = file_name+"_"+str(startyear)+"_"+str(endyear)
    return(name)

##get the netcdf files to process 
netc_files = [os.path.join(netc_folder,f) for f in os.listdir(netc_folder) if f.endswith('.nc4')]

#open the shapefile to get the geopandas set up
catchshape = geopandas.read_file(catchshapefile, crs=crs_code)
#get the name from the shapefile - used for naming output
shapename = catchshapefile.split("\\")[-1]
shapename = shapename.split(".")[0]
#set up the output lists
out_list = []
out_dates = []
for i,file in enumerate(netc_files):
    event_name = file_name(file)
    dates = days_list(file)
    out_dates.extend(dates)
    net_open = xarray.open_dataset(file)
    var_open = net_open.pr
    var_open.rio.set_spatial_dims(x_dim="lat", y_dim="lon", inplace=True)
    var_open.rio.write_crs(crs_code, inplace=True)
    for j in range(0,len(var_open.time)):
        clipped = var_open[j].rio.clip(catchshape.geometry.apply(mapping), drop=False)
        values = clipped.values
        maskedvalues = np.nanmean(values)*86400
        out_list.append(float(maskedvalues))
out_df =pd.DataFrame({'Date':out_dates,'PR':out_list})
out_df.to_csv(os.path.join(filepath,(str(shapename+"_"+event_name+".csv")))
              ,sep=',',index=False)

