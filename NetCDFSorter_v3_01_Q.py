# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 12:05:09 2019

@author: tomwillis
"""


from netCDF4 import Dataset
#from shapely import geometry
import numpy as np
import pandas as pd
import ogr,osr
import shapely
import shapely.geometry
import json
import os
import gdal
#import copy
import calendar
import time
import re
from datetime import date, timedelta, datetime

start_time = time.time()
filepath = r"D:\BarotseFloodplain\CATCHMAL\Analysis\Python\ISIMIP\PR Data"

#input files
netc_file = r"D:\BarotseFloodplain\CATCHMAL\Data\ClimateData\ISIMIP\jules-w1_ipsl-cm5a-lr_ewembi_historical_nosoc_co2_dis_global_daily_1981_1990.nc4"
evap_file = r"D:\BarotseFloodplain\CATCHMAL\Data\ClimateData\ISIMIP\clm45_hadgem2-es_ewembi_historical_2005soc_co2_evap_global_daily_1981_1990.nc4"
clip_domain = r"D:\BarotseFloodplain\CATCHMAL\Data\Shapefiles\RoughUpZamOutline_WGS_v0-01.shp"
#clip_domain = r"D:\BarotseFloodplain\CATCHMAL\Data\Shapefiles\BarotseFloodplain_WGS_v0-04.shp"
##filepath = r"D:\BarotseFloodplain\CATCHMAL\Data\ClimateData\ISIMIP\clm45_ipsl-cm5a-lr_ewembi_historical_2005soc_co2_dis_global_daily_1981_1990.nc4"

startyear = 1980
endyear = 1990

enddate = datetime(endyear + 1,1,1)

#Raster Parameters
driver = gdal.GetDriverByName('GTiff')
spatial_reference = 4326

def days_list(netc_file): 
    daylist = []
    nc_open = Dataset(netc_file)
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
    
def latlonlist(netc_file):
    datain = Dataset(netc_file)
    lat = datain.variables['lat'][:]
    ##sort descending
    lat[::-1].sort()
    lon = datain.variables['lon'][:]
    return(lat,lon)
    
def day_frame(daylist):
    yearlist = []
    monthlist = []
    daynumlist = []
    for i,day in enumerate(daylist):
        a = day.year
        b = day.month
        c = day.day
        yearlist.append(a)
        monthlist.append(b)
        daynumlist.append(c)
    dayframe = pd.DataFrame(daylist)
    dayframe = pd.DataFrame(pd.DatetimeIndex(daylist))
    dayframe = dayframe.rename(index =str,columns = {0:'date'})
    dayframe['frameref'] = range(0, len(dayframe))
    dayframe['year'] = yearlist
    dayframe['month'] = monthlist
    dayframe['day'] = daynumlist
    dayframe = dayframe.set_index(pd.DatetimeIndex(dayframe['date']))
    return (dayframe)
    
def yearrange(startyear,endyear,mdayframe):
    startdate = datetime.strptime((str(startyear)+"0101"),"%Y%m%d")
    enddate = datetime.strptime((str(endyear)+"1231"),"%Y%m%d")
    years = range(startdate.year, enddate.year+1,1)
    years = list(years)
    return(years)
    
def clipmask(clip_domain,lat,lon):
    tab = np.zeros([(len(lat)*len(lon)),2])
    count = 0
    for i in range(0,len(lat)):
        for j in range(0,len(lon)):
            tab[count,0] = lat[i]
            tab[count,1] = lon[j]
            count = count + 1
    shp = ogr.Open(clip_domain)
    lyr = shp.GetLayer(0)
    feat = lyr.GetFeature(0)
    first = feat.ExportToJson()
    first_a = json.loads(first)
    #shp_geom = shape(first_a['geometry'])
    poly_data = first_a["geometry"]["coordinates"][0]
    poly = shapely.geometry.Polygon(poly_data)
    mask = np.array([poly.contains(shapely.geometry.Point(y, x)) for x, y in tab])
    return(mask)

def runoff_calc(dateselect,nc_open,et_open,clip_domain,tablength):
    outtab = np.zeros([tablength,1])
    monthstart = min(dateselect.frameref)
    monthend = max(dateselect.frameref)+1 ##May need to add 1 here
    #latmin,latmax,lonmin,lonmax = array_clip(clip_domain,lat,lon)
    pr = nc_open.variables['dis'][:][monthstart:monthend]
   ## et = et_open.variables['et'][:][monthstart:monthend]
## loop through each variable
    outframe = np.zeros_like(pr[0])
    for i,frame in enumerate(pr):
        a = (frame>30)
        a = a.astype(int)
        outframe = outframe + a
    count = 0
    for j,element in np.ndenumerate(outframe):
        outtab[count,0] = (element)
        count = count + 1
    return(outtab)

def num_dim(max_dim,min_dim,cellsize):
    num_out = int((max_dim - min_dim)/cellsize)
    return(num_out)
 
def spatialref(sp_rf):
    sr = osr.SpatialReference()
    sr.ImportFromEPSG(sp_rf)
    sr = sr.ExportToWkt()
    return(sr)
    
def monthly_average(tab,months,tablength):
    #drop the lat/lon columns
    tabs = tab[:,2:]
    month_out = []
    #create output tab
    monthtab = tab[:,0:2]
    for i,month in enumerate(months):
        month_out.append(str((calendar.month_name[month])))
        tabslice = tabs[:, i::12]
        monthmeans = np.mean(tabslice,axis=1)
        monthmeans = np.reshape(monthmeans,(tablength,1))
        monthtab = np.append(monthtab,monthmeans,axis=1)
    outdf = pd.DataFrame(data=monthtab)
    outdf = outdf.dropna(how='any')
    outdf.columns = [month_out[x-2] if x > 1 else x for x in outdf.columns]
    outdf = outdf.rename(columns={0:'Latitude',1:'Longitude'})
    outdf.to_csv(os.path.join(filepath,"NETCDF_MonthlyAverageDischarge_Data_v1-00.csv") ,sep=',',index=False)
    print("Monthly Averages Computed")
    return (outdf)

def monthly_write(month_table,months,lat,lon):
    ##figure out the dimensions of the data for the gtiff
    sr = spatialref(spatial_reference)
    cellsize = lat[0] - lat[1]
    xmin = min(lon)
    xmax = max(lon)
    ymin = min(lat)
    ymax = max(lat)
    num_cols = num_dim(xmax,xmin,cellsize)+1
    num_rows = num_dim(ymax,ymin,cellsize)+1
    ##going to sort the table to write out in the correct way
    #month_sort = month_table.sort_values(by = ['Latitude','Longitude'],ascending=False)
    month_sort = month_table
    for i,month in enumerate(months):
        month_name = str((calendar.month_name[month]))
        ##create file names
        tiff_name = os.path.join(filepath,str("NETCDF_Q"+month_name+"_Average.tif"))
        # get the data, split to array
        month_data = np.array(month_sort[[month_name]])
        month_array = month_data.reshape((len(lat),len(lon)))
        new_dataset = driver.Create(tiff_name,
                                    num_cols,    # number of columns
                                    num_rows,    # number of rows
                                    1,              # number of bands
                                    gdal.GDT_Float32)  # datatype of the raster
        new_dataset.SetProjection(sr)
        new_dataset.SetGeoTransform((xmin,cellsize,0,ymax,0,-cellsize))
        new_band = new_dataset.GetRasterBand(1)
        new_band.SetNoDataValue(-9999)
        new_band.WriteArray(month_array)
        new_band.FlushCache()
        new_dataset = None
    print("Monthly Raster WriteOut Complete")
    return()

def main():
    daylist = days_list(netc_file)
    #daylist = [d for d in daylist if d < datetime.date(enddate)]
    dayframe = day_frame(daylist)   
    #years = list(dayframe.year.unique())
    #years = yearrange(startyear,endyear,dayframe)
    months = list(dayframe.month.unique())
    
    ##new method for ordering months and years - extract unique year month
    date_strings = [d.strftime('%Y/%m') for d in daylist]
    date_months = list(set(date_strings))
    date_months.sort()
    #years = [np.int64(monthyear.split('/')[0]) for monthyear in date_months]
    #months = list(dayframe.month.unique())
    
    nc_open = Dataset(netc_file)
    et_open = Dataset(evap_file)
    lat = nc_open.variables['lat'][:]
    lon = nc_open.variables['lon'][:]
# build the output table
    tab = np.zeros([(len(lat)*len(lon)),2])
    tablength = (len(lat)*len(lon))
    dates = []
    count = 0
    for i in range(0,len(lat)):
        for j in range(0,len(lon)):
            tab[count,0] = lat[i]
            tab[count,1] = lon[j]
            count = count + 1
    #now for each year and each month, calculate wet days
    for i,monthyear in enumerate(date_months):
        month = np.int64(monthyear.split('/')[1])
        year = np.int64(monthyear.split('/')[0])
        print(str((calendar.month_name[month]+" "+str(year))))
        dates.append(str((calendar.month_name[month]+" "+str(year))))
        dateselect = dayframe.loc[(dayframe['year'] == year) & (dayframe['month'] == month)]
        runoff_component = runoff_calc(dateselect,nc_open,et_open,clip_domain,tablength)
        tab = np.append(tab,runoff_component,axis=1)
    month_table = monthly_average(tab,months,tablength)
    monthly_write(month_table,months,lat,lon)
    outdf = pd.DataFrame(data=tab)
    outdf = outdf.dropna(how='any')
    outdf.columns = [dates[x-2] if x > 1 else x for x in outdf.columns]
    outdf = outdf.rename(columns={0:'Latitude',1:'Longitude'})
    outdf.to_csv(os.path.join(filepath,"NETCDF_MonthlyDischargeDays_Data_v1-00.csv") ,sep=',',index=False)
    print ("Processing Complete -  Time to Execute  %s:" % (time.time()
    - start_time))
    return ()

if __name__== '__main__':
    main()