# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 16:06:37 2021

@author: wangplll
"""

"""
Data variables:
    raw:
         temperatures           (time) float64 ...  degree_Celsius
         dewpoints              (time) float64 ...  degree_Celsius
         slp                    (time) float64 ...  hPa
         stnlp                  (time) float64 ...  hPa
         windspeeds             (time) float64 ...  meters per second
         winddirs               (time) float64 ...  degree
         total_cloud_cover      (time) float64 ...  1
         low_cloud_cover        (time) float64 ...  1
         mid_cloud_cover        (time) float64 ...  1
         high_cloud_cover       (time) float64 ...  1
         precip1_depth          (time) float64 ...  mm Depth of Precipitation Reported in 1 hour
         precip2_depth          (time) float64 ...  mm Depth of Precipitation Reported in 2 hour
         precip3_depth          (time) float64 ...
         precip6_depth          (time) float64 ...
         precip9_depth          (time) float64 ...
         precip12_depth         (time) float64 ...
         precip15_depth         (time) float64 ...
         precip18_depth         (time) float64 ...
         precip24_depth         (time) float64 ...  mm Depth of Precipitation Reported in 24 hour
         cloud_base             (time) float64 ...  meters Cloud base of lowest cloud layer
         wind_gust              (time) float64 ...  meters per second   Wind Gust Speed at mast height (~10m)
         past_sigwx1            (time) float64 ...  1  Reported past significant weather phenomena
         input_station_id       (time) object ...
         quality_control_flags  (time, test) float64 ...
         flagged_obs            (time, flagged) float64 ...
         porting_stats        (reporting_v, reporting_t, reporting_2) float64 ...
         
         longitude, latitude, elevation
        
    heatstress:
        temperatures                (time) float64 ...
        dewpoints                   (time) float64 ...
        windspeeds                  (time) float64 ...
        temperature_humidity_index  (time) float64 ...
        wet_bulb_globe_temperature  (time) float64 ...
        humidex                     (time) float64 ...
        apparent_temperature        (time) float64 ...
        heat_index                  (time) float64 ...
        input_station_id            (time) object ...
        
    humidity:
        temperatures               (time) float64 ...
        dewpoints                  (time) float64 ...
        slp                        (time) float64 ...
        vapor_pressure             (time) float64 ...
        saturation_vapor_pressure  (time) float64 ...
        wet_bulb_temperature       (time) float64 ...
        specific_humidity          (time) float64 ...
        relative_humidity          (time) float64 ...
        input_station_id           (time) object ...
        
    Total stations: 9278 all over the world
"""

# In[0] Impport python libs
import xarray as xr
import os
import pandas as pd
import geopandas as gpd
import numpy as np
import datetime

# In[1] Import dataset
path0 = r"/home/wangplll/paper/HadISD_data"
filelist0 = sorted(os.listdir(path0))


def addCoords(data, lon, lat, alt, ID): # Add coordinates
    res = xr.DataArray(data.values,
                    coords={
                        "time" : data.time.values,
                        "lon" : lon,
                        "lat" : lat,
                        "alt" : alt,
                        "id"  : ID,
                        },
                    dims= "time",)
    return res


def cal_value(clm_factor, lon, lat, alt, ID):
    """
    xarray also supports a notion of “virtual” or “derived” coordinates 
    for datetime components implemented by pandas, including “year”, 
    “month”, “day”, “hour”, “minute”, “second”, “dayofyear”, “week”, 
    “dayofweek”, “weekday” and “quarter”, "season": 
        e.g., t_1973_2020["time.day"], t_1973_2020["time"].dt.day
    identify (labeled) per data belongs to selected datetime components(eg, "season")
    """
    # Select 1973-01-01 ->2020-12-31 (index by dimension coordinate labels)
    clm_factor_1973_2020 = clm_factor.sel(time=slice("1973-01-01", "2020-12-31"))
    # t_1973_01 = t.sel(time="1973-01") # select dataset of 1973-01
    # t.sel(time = datetime.time(12))   # select dataset of 12:00:00 of every day
        
        
    ## Extract day max, min and mean data from hour data within one day
    """
     Most DateOffsets have associated frequencies strings, or offset aliases, that 
     can be passed into freq keyword arguments. The available date offsets and 
     associated frequency strings can be found below:
         https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#dateoffset-objects
    eg:Day       'D'              one absolute day
       Hour      'H'              one hour
       Minute    'T' or 'min'     one minute
       Generally, these frequency strings are used with "resample" methods
       e.g.: temp = t_1973_2020.resample(time="1D").mean()
     """
    clm_factor_day_mean = clm_factor_1973_2020.resample(time="1D").mean()
    clm_factor_day_max = clm_factor_1973_2020.resample(time="1D").max()
    clm_factor_day_min = clm_factor_1973_2020.resample(time="1D").min()
        
    ## Add coords
    clm_factor_mean = addCoords(clm_factor_day_mean, lon, lat, alt, ID)
    clm_factor_max = addCoords(clm_factor_day_max, lon, lat, alt, ID)
    clm_factor_min = addCoords(clm_factor_day_min, lon, lat, alt, ID)
    
    return clm_factor_mean, clm_factor_max, clm_factor_min

def file_concat0(clm_factor_mean, clm_factor_max, clm_factor_min):
    clm_factor_mean0 = clm_factor_mean
    clm_factor_max0 = clm_factor_max
    clm_factor_min0 = clm_factor_min
    return clm_factor_mean0, clm_factor_max0, clm_factor_min0
    
def file_concat(clm_factor_mean0, clm_factor_max0, clm_factor_min0, clm_factor_mean, clm_factor_max, clm_factor_min):
    # Concatinate DataArray 
    clm_factor_mean0 = xr.concat([clm_factor_mean0, clm_factor_mean], dim = "station")
    clm_factor_max0 = xr.concat([clm_factor_max0, clm_factor_max], dim = "station")
    clm_factor_min0 = xr.concat([clm_factor_min0, clm_factor_min], dim = "station")
    
    return clm_factor_mean0, clm_factor_max0, clm_factor_min0


## Extract the requires data
ifile = 0
for filename0 in filelist0:
    path1 = path0 + "/" + filename0 + "/" + "decompress"
    filelist1 = sorted(os.listdir(path1))
    
    for filename1 in filelist1:
        path2 = path1 + "/" + filename1
        filename2 = os.listdir(path2)
        path3 = path2 + "/" + filename2[0]
        filename_in = filename0 + "_" + filename2[0][:-3]
        file = xr.open_dataset(path3)
        
        ## lat lon alt id
        lon = file.longitude.values[0]
        lat = file.latitude.values[0]
        alt = file.elevation.values[0]
        ID = file["input_station_id"].values[0]
        
        if int(str(file.time[-1].values)[:4])<2021:
            continue

        
        # Get specific climate factor (1973-2020) 
        if filename0[18:-4] == "heat_stress": # heat stress 
            wbgt = file["wet_bulb_globe_temperature"]
            at = file["apparent_temperature"]
            hi = file["heat_index"]
            
            # wet_bulb_globe_temperature
            [wbgtmean0, wbgtmax0, wbgtmin0] = cal_value(wbgt, lon, lat, alt, ID)
            # apparent_temperature
            [atmean0, atmax0, atmin0] = cal_value(at, lon, lat, alt, ID)
            # heat_index
            [himean0, himax0, himin0] = cal_value(hi, lon, lat, alt, ID)
            
            if ifile == 92:
                [wbgtmean, wbgtmax, wbgtmin] = file_concat0(wbgtmean0, wbgtmax0, wbgtmin0)
                [atmean, atmax, atmin]       = file_concat0(atmean0, atmax0, atmin0)
                [himean, himax, himin]       = file_concat0(himean0, himax0, himin0) 
            else:
                [wbgtmean, wbgtmax, wbgtmin] = file_concat(wbgtmean, wbgtmax, wbgtmin, wbgtmean0, wbgtmax0, wbgtmin0)
                [atmean, atmax, atmin]       = file_concat(atmean, atmax, atmin, atmean0, atmax0, atmin0)
                [himean, himax, himin]       = file_concat(himean, himax, himin, himean0, himax0, himin0)
            
           
        elif filename0[18:-4] == "humidity":  # humidity
            p = file["vapor_pressure"]
            wbt = file["wet_bulb_temperature"]
            rh = file["relative_humidity"]
            
            # vapor_pressure
            [pmean0, pmax0, pmin0] = cal_value(p, lon, lat, alt, ID)
            # wet_bulb_temperature
            [wbtmean0, wbtmax0, wbtmin0] = cal_value(wbt, lon, lat, alt, ID)
            # relative_humidity
            [rhmean0, rhmax0, rhmin0] = cal_value(rh, lon, lat, alt, ID)
            
            if ifile == 184:
                [pmean, pmax, pmin]       = file_concat0(pmean0, pmax0, pmin0)
                [wbtmean, wbtmax, wbtmin] = file_concat0(wbtmean0, wbtmax0, wbtmin0)
                [rhmean, rhmax, rhmin]    = file_concat0(rhmean0, rhmax0, rhmin0)
            else:
                [pmean, pmax, pmin]       = file_concat(pmean, pmax, pmin, pmean0, pmax0, pmin0)
                [wbtmean, wbtmax, wbtmin] = file_concat(wbtmean, wbtmax, wbtmin, wbtmean0, wbtmax0, wbtmin0)
                [rhmean, rhmax, rhmin]    = file_concat(rhmean, rhmax, rhmin, rhmean0, rhmax0, rhmin0)
            
            
        else:  # raw data
            t = file["temperatures"]
            windspeeds = file["windspeeds"]
            dewpoints = file["dewpoints"]
            stnlp = file["stnlp"]
            tcc = file["total_cloud_cover"]
            lcc = file["low_cloud_cover"]
                
            # temperature
            [tmean0, tmax0, tmin0] = cal_value(t, lon, lat, alt, ID)
            # dewpoints
            [dpmean0, dpmax0, dpmin0] = cal_value(dewpoints, lon, lat, alt, ID)
            # station surface air pressure
            [stnlpmean0, stnlpmax0, stnlpmin0] = cal_value(stnlp, lon, lat, alt, ID)
            # Total cloud cover
            [tccmean0, tccmax0, tccmin0] = cal_value(tcc, lon, lat, alt, ID)
            # Low cloud cover
            [lccmean0, lccmax0, lccmin0] = cal_value(lcc, lon, lat, alt, ID)
            # windspeeds
            [winmean0, winmax0, winmin0] = cal_value(windspeeds, lon, lat, alt, ID)
            
            if ifile == 0:
                [tmean, tmax, tmin]             = file_concat0(tmean0, tmax0, tmin0)
                [dpmean, dpmax, dpmin]          = file_concat0(dpmean0, dpmax0, dpmin0)
                [stnlpmean, stnlpmax, stnlpmin] = file_concat0(stnlpmean0, stnlpmax0, stnlpmin0)
                [tccmean, tccmax, tccmin]       = file_concat0(tccmean0, tccmax0, tccmin0)
                [lccmean, lccmax, lccmin]       = file_concat0(lccmean0, lccmax0, lccmin0)
                [winmean, winmax, winmin]       = file_concat0(winmean0, winmax0, winmin0)
            else:
                [tmean, tmax, tmin]             = file_concat(tmean, tmax, tmin, tmean0, tmax0, tmin0)
                [dpmean, dpmax, dpmin]          = file_concat(dpmean, dpmax, dpmin, dpmean0, dpmax0, dpmin0)
                [stnlpmean, stnlpmax, stnlpmin] = file_concat(stnlpmean, stnlpmax, stnlpmin, stnlpmean0, stnlpmax0, stnlpmin0)
                [tccmean, tccmax, tccmin]       = file_concat(tccmean, tccmax, tccmin, tccmean0, tccmax0, tccmin0)
                [lccmean, lccmax, lccmin]       = file_concat(lccmean, lccmax, lccmin, lccmean0, lccmax0, lccmin0)
                [winmean, winmax, winmin]       = file_concat(winmean, winmax, winmin, winmean0, winmax0, winmin0)
                    
        ifile = ifile + 1
        print(filename_in)
 
    print(filename0)    

# Output as .nc file
outpath = r"/home/wangplll/paper/Global_HPT/data2/"

# tmean.to_netcdf(outpath + "tmean.nc")
# tmax.to_netcdf(outpath + "tmax.nc")
# tmin.to_netcdf(outpath + "tmin.nc")

# dpmean.to_netcdf(outpath + "dpmean.nc")
# dpmax.to_netcdf(outpath + "dpmax.nc")
# dpmin.to_netcdf(outpath + "dpmin.nc")

# stnlpmean.to_netcdf(outpath + "stnlpmean.nc")
# stnlpmax.to_netcdf(outpath + "stnlpmax.nc")
# stnlpmin.to_netcdf(outpath + "stnlpmin.nc")

# tccmean.to_netcdf(outpath + "tccmean.nc")
# tccmax.to_netcdf(outpath + "tccmax.nc")
# tccmin.to_netcdf(outpath + "tccmin.nc")

# lccmean.to_netcdf(outpath + "lccmean.nc")
# lccmax.to_netcdf(outpath + "lccmax.nc")
# lccmin.to_netcdf(outpath + "lccmin.nc")

# winmean.to_netcdf(outpath + "winmean.nc")
# winmax.to_netcdf(outpath + "winmax.nc")
# winmin.to_netcdf(outpath + "winmin.nc")

# wbgtmean.to_netcdf(outpath + "wbgtmean.nc")
# wbgtmax.to_netcdf(outpath + "wbgtmax.nc")
# wbgtmin.to_netcdf(outpath + "wbgtmin.nc")

# atmean.to_netcdf(outpath + "atmean.nc")
# atmax.to_netcdf(outpath + "atmax.nc")
# atmin.to_netcdf(outpath + "atmin.nc")

# himean.to_netcdf(outpath + "himean.nc")
# himax.to_netcdf(outpath + "himax.nc")
# himin.to_netcdf(outpath + "himin.nc")

# pmean.to_netcdf(outpath + "pmean.nc")
# pmax.to_netcdf(outpath + "pmax.nc")
# pmin.to_netcdf(outpath + "pmin.nc")

# wbtmean.to_netcdf(outpath + "wbtmean.nc")
# wbtmax.to_netcdf(outpath + "wbtmax.nc")
# wbtmin.to_netcdf(outpath + "wbtmin.nc")

# rhmean.to_netcdf(outpath + "rhmean.nc")
# rhmax.to_netcdf(outpath + "rhmax.nc")
# rhmin.to_netcdf(outpath + "rhmin.nc")

# Notice that: something wrong with stnlpmean
