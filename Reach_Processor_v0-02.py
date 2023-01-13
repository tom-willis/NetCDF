# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 16:01:13 2020

@author: geotw
"""

import numpy as np
import pandas as pd
import os
import argparse

#startyear = '2001'

#daterange = "1999_2009"

#wustype = '20plus'

#PPE = "PPE4"

infile = r"D:\BarotseFloodplain\CATCHMAL\Analysis\SWAT\BarotseFP_CATCHMAL\Scenarios\Default\TxtInOut\output.rch"
sub_default = r"D:\BarotseFloodplain\CATCHMAL\Analysis\SWAT\Output\Subbasin_Default.csv"
outfolder = r"D:\BarotseFloodplain\CATCHMAL\Analysis\SWAT\Output"

def main():
    ##open input file and get inputs
    parser = argparse.ArgumentParser()
    parser.add_argument("--file","-f",type=str,required =True)
    args = parser.parse_args()
    file = args.file
    File = open(file)
    params = []
    for line in File:
        params.append(line)
    ## sort inputs to be used
    PPE = params[0].rstrip()
    wustype = params[1].rstrip()
    daterange = params[2].rstrip()
    startyear = params[3].rstrip()
    print(wustype)
    print(daterange)
    
    name = str(PPE+"_"+daterange+"_"+wustype+"_Subcatchment.csv")
    watname = str(PPE+"_"+daterange+"_"+wustype+"_Totaloutflow.csv")
    outfile = os.path.join(outfolder,name)
    watout = os.path.join(outfolder,watname)
    
    lines = np.loadtxt(infile,skiprows=9, usecols=[1,5],unpack=False)
    a = pd.DataFrame(lines)
    a = a.rename(columns={0: "a", 1: "SUM_FLOW"})
    a = a.set_index("a")
    result = a.groupby(a.index).sum()
    sub_def = pd.read_csv(sub_default)
    sub_def = sub_def.set_index("Subbasin")
    
    b = sub_def.join(result)
    b.to_csv(outfile)
    
    watfile = r"C:\Users\geotw\Desktop\MetOfficeCSSP\Modelling\SWAT\UYR_QSWAT_v1-02\Scenarios\Default\TxtInOut\watout.dat"
    #watpd = pd.read_csv(watfile,skiprows=2,delimiter=r"\s+")
    wat = np.loadtxt(watfile,skiprows=6, usecols=[3],unpack=False)
    #define a min below which the river cannot go
    wat[wat<10]=10
    #watdf = pd.DataFrame(wat)
    dates = wat.size
    
    #CHANGE START DATE
    startdate = str('01-01-'+startyear) 
    
    rng = pd.date_range(startdate, periods=dates, freq='D')

    df = pd.DataFrame({'Date': rng, 'OutFlow(m3/s)': wat})
    df.to_csv(watout,index=False)

if __name__== '__main__':
    main()