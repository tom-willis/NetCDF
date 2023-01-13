# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 12:05:09 2019

@author: tomwillis
"""


from netCDF4 import Dataset
import netCDF4
#from shapely import geometry
import numpy as np
import pandas as pd
import ogr,osr
import shapely
import json
import os
import gdal
import time
import datetime
from datetime import date, timedelta

###LAT and LON SWITCHED IN ARRAY CLIP!

start_time = time.time()
filepath = r"D:\BarotseFloodplain\CATCHMAL\Data\SWAT\WeatherDataInputs"
format_str = "%Y-%m-%d"

#testing
#filepath = r"C:\Users\geotw\OneDrive - University of Leeds\Python"
#test = r"C:\Users\geotw\OneDrive - University of Leeds\Python\TestGrid.csv"
#month_table = pd.read_csv(test)

#input files
tempmax_file = r"D:\BarotseFloodplain\CATCHMAL\Data\ClimateData\ISIMIP\tasmax_day_IPSL-CM5A-LR_historical_r1i1p1_EWEMBI_19810101-19901231.nc4"
tempmin_file = r"D:\BarotseFloodplain\CATCHMAL\Data\ClimateData\ISIMIP\tasmin_day_IPSL-CM5A-LR_historical_r1i1p1_EWEMBI_19810101-19901231.nc4"
netc_file = r"D:\BarotseFloodplain\CATCHMAL\Data\ClimateData\ISIMIP\pr_day_IPSL-CM5A-LR_historical_r1i1p1_EWEMBI_19810101-19901231.nc4"
relhum_file = r"D:\BarotseFloodplain\CATCHMAL\Data\ClimateData\ISIMIP\hurs_day_IPSL-CM5A-LR_historical_r1i1p1_EWEMBI_19810101-19901231.nc4"
wind_file = r"D:\BarotseFloodplain\CATCHMAL\Data\ClimateData\ISIMIP\sfcWind_day_IPSL-CM5A-LR_historical_r1i1p1_EWEMBI_19810101-19901231.nc4"
shsol_file = r"D:\BarotseFloodplain\CATCHMAL\Data\ClimateData\ISIMIP\rsds_day_IPSL-CM5A-LR_historical_r1i1p1_EWEMBI_19810101-19901231.nc4"
lowsol_file = r"D:\BarotseFloodplain\CATCHMAL\Data\ClimateData\ISIMIP\rlds_day_IPSL-CM5A-LR_historical_r1i1p1_EWEMBI_19810101-19901231.nc4"

evap_file =r"D:\BarotseFloodplain\CATCHMAL\Data\ClimateData\ISIMIP\clm45_ipsl-cm5a-lr_ewembi_historical_2005soc_co2_evap_global_daily_1981_1990.nc4"

clip_domain = r"D:\BarotseFloodplain\CATCHMAL\Data\SWAT\BFP_SWAT_Area_500m__wgs_v2-00.shp"
basinshapes = r"D:\BarotseFloodplain\CATCHMAL\Data\SWAT\BFP_SWAT_Area_500m__wgs_v2-00.shp"

#Raster Parameters
driver = gdal.GetDriverByName('GTiff')
spatial_reference = 4326

def days_list(netc_file): 
    daylist = []
    nc_open = Dataset(netc_file)
    days = nc_open.variables['time'][:]
    start = datetime.datetime(1661,1,1,0,0,0) #should be made call variable
    for i,day in enumerate(days):
        delta = timedelta(day)
        offset = start + delta
        daylist.append(offset)
    return(daylist)
    
def days_list_greg(netc_file):
    ##converts the fucking stupid 360 day year to something humans use
    nc_open = netCDF4.Dataset(netc_file)
    time_var = nc_open.variables['time']
    time_convert = netCDF4.num2date(time_var[:], time_var.units, time_var.calendar)
    no_days =range(0,len(time_convert))
    df = pd.DataFrame(no_days)
    df['dates'] = (time_convert)
    return()


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

def array_clip(geom,lat,lon):
    xmin,xmax,ymin,ymax = geom.GetEnvelope()
    latbounds = [int(ymin),round(ymax)]##swapped added for southers
    lonbounds = [int(xmin),int(xmax)]
    latli = np.argmin( np.abs( lat - latbounds[0] ) )
    latui = np.argmin( np.abs( lat - latbounds[1] ) ) 
    # longitude lower and upper index
    lonli = np.argmin( np.abs( lon - lonbounds[0] ) )
    lonui = np.argmin( np.abs( lon - lonbounds[1] ) )
    return (latui,latli,lonli,lonui)

##function to collect the dimensions of the output raster
def num_dim(max_dim,min_dim,cellsize):
    num_out = int((max_dim - min_dim)/cellsize)
    return(num_out)

def featurecount(grid_domain):
    shp = ogr.Open(grid_domain)
    lyr = shp.GetLayer()
    grid_count = lyr.GetFeatureCount()
    return(grid_count)
 
def spatialref(sp_rf):
    sr = osr.SpatialReference()
    sr.ImportFromEPSG(sp_rf)
    sr = sr.ExportToWkt()
    return(sr)

def SWAT_input(start_year,no_days,nc_open,latmin,latmax,lonmin,lonmax,subbasinid):
    start_point = start_year
    tab = np.zeros(no_days+1)
    tab[0] = start_point
    for i in range(0,no_days):
        precip_frame = nc_open.variables['pr'][i,latmin:latmax, lonmin:lonmax]
        precip_frame = precip_frame
        day_mean = precip_frame.mean() * 86400
        tab[i+1] = day_mean
    np.savetxt(os.path.join(filepath,"p_subbasin"+str(subbasinid)+".txt"),tab,fmt='%10.5f',delimiter=',')
    return()

def SWAT_temp_input(start_year,no_days,tempmax_open,tempmin_open,latmin,latmax,lonmin,lonmax,subbasinid):
    start_point = start_year
    tab = np.zeros(shape=(no_days+1,2))
    tab[0,0] = start_point
    for i in range(0,no_days):
        tempmax_frame = tempmax_open.variables['tasmax'][i,latmin:latmax, lonmin:lonmax]
        tempmin_frame = tempmin_open.variables['tasmin'][i,latmin:latmax, lonmin:lonmax]
        day_max = tempmax_frame.mean() - 273.15
        day_min = tempmin_frame.mean() - 273.15
        tab[i+1,0] = day_max
        tab[i+1,1] = day_min
    np.savetxt(os.path.join(filepath,"t_subbasin"+str(subbasinid)+".txt"),tab,fmt='%10.5f',delimiter=',')
    return()

def SWAT_relhum_input(start_year,no_days,relhum_open,latmin,latmax,lonmin,lonmax,subbasinid):
    start_point = start_year
    tab = np.zeros(no_days+1)
    tab[0] = start_point
    for i in range(0,no_days):
        tempmax_frame = relhum_open.variables['hurs'][i,latmin:latmax, lonmin:lonmax]
        day_max = tempmax_frame.mean() / 100
        tab[i+1] = day_max
    np.savetxt(os.path.join(filepath,"r_subbasin"+str(subbasinid)+".txt"),tab,fmt='%10.5f',delimiter=',')
    return()

def SWAT_wind_input(start_year,no_days,wind_open,latmin,latmax,lonmin,lonmax,subbasinid):
    start_point = start_year
    tab = np.zeros(no_days+1)
    tab[0] = start_point
    for i in range(0,no_days):
        tempmax_frame = wind_open.variables['sfcWind'][i,latmin:latmax, lonmin:lonmax]
        day_max = tempmax_frame.mean()
        tab[i+1] = day_max
    np.savetxt(os.path.join(filepath,"w_subbasin"+str(subbasinid)+".txt"),tab,fmt='%10.5f',delimiter=',')
    return()

def SWAT_evap_input(start_year,no_days,shsol_open,lowsol_open,latmin,latmax,lonmin,lonmax,subbasinid):
    start_point = start_year
    tab = np.zeros(no_days+1)
    tab[0] = start_point
    for i in range(0,no_days):
        tempmax_frame = shsol_open.variables['rsds'][i,latmin:latmax, lonmin:lonmax]
        templon_frame = lowsol_open.variables['rlds'][i,latmin:latmax, lonmin:lonmax]
        shsol_max = tempmax_frame.mean() /106
        lowsol_max = templon_frame.mean() / 106
        day_max = shsol_max + lowsol_max
        tab[i+1] = day_max
    np.savetxt(os.path.join(filepath,"s_subbasin"+str(subbasinid)+".txt"),tab,fmt='%10.5f',delimiter=',')
    return()

def main():
    #daylist = days_list(netc_file)
    #dayframe = day_frame(daylist)
    #years = list(dayframe.year.unique())
    year = netc_file.split('-')[-2]
    year = int(year.split('_')[-1])
    basincount = featurecount(basinshapes)
    nc_open = Dataset(netc_file)
    tempmax_open = Dataset(tempmax_file)
    tempmin_open = Dataset(tempmin_file)
    relhum_open = Dataset(relhum_file)
    wind_open = Dataset(wind_file)
    shsol_open = Dataset(shsol_file)
    lowsol_open = Dataset(lowsol_file)
    lat = nc_open.variables['lat'][:]
    lon = nc_open.variables['lon'][:]
    time = nc_open.variables['time'][:]
    no_days = len(time)
    start_year = year
    basinshp = ogr.Open(basinshapes)
    lyr = basinshp.GetLayer()    
    for i in range(0,basincount):
        print("Subbasin "+str(i))
        #start_year = start_year
        feat = lyr.GetFeature(i)
        geom = feat.GetGeometryRef()
        latmin,latmax,lonmin,lonmax=array_clip(geom,lat,lon)
        subbasinid = feat.GetField('Subbasin')
        SWAT_input(start_year,no_days,nc_open,latmin,latmax,lonmin,lonmax,subbasinid)
        SWAT_temp_input(start_year,no_days,tempmax_open,tempmin_open,latmin,latmax,lonmin,lonmax,subbasinid)
        SWAT_relhum_input(start_year,no_days,relhum_open,latmin,latmax,lonmin,lonmax,subbasinid)
        SWAT_wind_input(start_year,no_days,wind_open,latmin,latmax,lonmin,lonmax,subbasinid)
        SWAT_evap_input(start_year,no_days,shsol_open,lowsol_open,latmin,latmax,lonmin,lonmax,subbasinid)

        #TRMM_yearoutput(start_year,dayframe,nc_open,lat,lon) 
    print ("Processing Complete -  Time to Execute  %s:" % (time.time()
    - start_time))


# def main():
#     daylist = days_list(netc_file)
#     dayframe = day_frame(daylist)
#     years = list(dayframe.year.unique())
#     months = list(dayframe.month.unique())
#     nc_open = Dataset(netc_file)
#     lat = nc_open.variables['lat'][:]
#     lon = nc_open.variables['lon'][:]
# # build the output table
#     tab = np.zeros([(len(lat)*len(lon)),2])
#     tablength = (len(lat)*len(lon))
#     date = []
#     count = 0
#     for i in range(0,len(lat)):
#         for j in range(0,len(lon)):
#             tab[count,0] = lat[i]
#             tab[count,1] = lon[j]
#             count = count + 1
#     #now for each year and each month, calculate wet days
#     for i,year in enumerate(years):
#         for j, month in enumerate(months):
#             print(str((calendar.month_name[month]+" "+str(year))))
#             date.append(str((calendar.month_name[month]+" "+str(year))))
#             dateselect = dayframe.loc[(dayframe['year'] == year) & (dayframe['month'] == month)]
#             runoff_component = runoff_calc(dateselect,nc_open,et_open,clip_domain,tablength)
#             tab = np.append(tab,runoff_component,axis=1)
#     month_table = monthly_average(tab,months,tablength)
#     mothly_write(month_table,month,lat,lon)
#     outdf = pd.DataFrame(data=tab)
#     outdf = outdf.dropna(how='any')
#     outdf.columns = [date[x-2] if x > 1 else x for x in outdf.columns]
#     outdf = outdf.rename(columns={0:'Latitude',1:'Longitude'})
#     outdf.to_csv(os.path.join(filepath,"NETCDF_MonthlyWetDays_Data_v1-00.csv") ,sep=',',index=False)
#     print ("Processing Complete -  Time to Execute  %s:" % (time.time()
#     - start_time))
#     return ()

if __name__== '__main__':
    main()