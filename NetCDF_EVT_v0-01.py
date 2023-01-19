# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 13:04:12 2023

@author: geotw
"""


from netCDF4 import Dataset
import numpy as np
import pandas as pd
import geopandas as gpd
import ogr
import json
import shapely
import shapely.geometry
import os
from datetime import timedelta, datetime

netc_file = r"C:\Users\geotw\OneDrive - University of Leeds\CASTOR\Data\General\UK_PETData_EA\hydro-pe_hadukgrid_pet_1km_day_20050101-20050131.nc"
catchshapefile = r"C:\Users\geotw\OneDrive - University of Leeds\CASTOR\Analysis\Model Version 0\Lowther\Shapefiles\LowtherCatchment_v0-01.shp"

def latlonlist(netc_file):
    datain = Dataset(netc_file)
    lat = datain.variables['y'][:]
    ##sort descending
    lat[::-1].sort()
    lon = datain.variables['x'][:]
    return(lat,lon)

def days_list(netc_file): 
    daylist = []
    nc_open = Dataset(netc_file)
    #time_since = nc_open.variables['time'].units
    #match = re.search(r'\d{4}-\d{1}-\d{1}', time_since)'
    #start = datetime.strptime(match.group(), '%Y-%m-%d').date()'
    days = nc_open.variables['time'][:]
    start = datetime(1800,1,1,0,0,0) #should be made call variable
    for i,day in enumerate(days):
        ##amend for ISIMIP - days provided in 32bit format,int required
        day = int(day)
        delta = timedelta(hours = day)
        offset = start + delta
        daylist.append(offset)
    return(daylist)

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

nc_open = Dataset(netc_file)
lat,lon = latlonlist(netc_file)
mask = clipmask(catchshapefile,lat,lon)

tablength = len(tab)
outtab = np.zeros([tablength,2])
a = nc_open.variables['pet'][1]

count = 0
for j,element in np.ndenumerate(a):
    outtab[count,0] = (element)
    outtab[count,1] = (mask[count])
    count = count + 1

##clipping https://gis.stackexchange.com/questions/433287/clip-netcdf-files-by-shapefile
    
import geopandas
import rasterio
import rasterio.features
import xarray

crs = nc_open.variables['crsOSGB'].towgs84
crs = [1,0,0,0,1,0,0]
crs_aff = rasterio.Affine(crs[0],crs[1],crs[2],crs[3],crs[4],crs[5])
sf = geopandas.read_file(catchshape)
shape_mask = rasterio.features.geometry_mask(sf.geometry,
                                             out_shape = (len(lat),len(lon)),
                                             transform = crs_aff)
shape_mask = xarray.DataArray(shape_mask , dims=("y", "x"))
nc_file = xarray.open_dataset(netc_file)
masked_netcdf_file = nc_file.where(shape_mask == True, drop=True)

##https://gis.stackexchange.com/questions/354782/masking-netcdf-time-series-data-from-shapefile-using-python

import geopandas
import rioxarray
import xarray
from shapely.geometry import mapping

net_open = xarray.open_dataset(netc_file)
pet_open = net_open.pet
pet_open.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
pet_open.rio.write_crs("epsg:27700", inplace=True)
catchshape = geopandas.read_file(catchshapefile, crs="epsg:27700")

clipped = pet_open.rio.clip(catchshape.geometry.apply(mapping), catchshape.crs, drop=False)
