# -*- coding: utf-8 -*-
"""
Created on Fri May  7 20:15:08 2021

@author: wangplll
"""
# In[0. Import python libs]
import xarray as xr
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.io.shapereader import Reader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy import signal
from scipy import stats
import calendar

plt.rcParams['savefig.dpi'] = 300 #Í¼Æ¬ÏñË
plt.rcParams['figure.dpi'] = 300  

# In[1. Import dataset]
path0 = r"/home/wangplll/paper/Global_HPT/data3/tmean.nc"  # Tmean
path1 = r"/home/wangplll/paper/Global_HPT/data3/atmean.nc" # Apparent temperature
path2 = r"/home/wangplll/paper/Global_HPT/data3/wbgtmean.nc" # Wet_bulb_globe_temperature 
path3 = r"/home/wangplll/paper/Global_HPT/data3/rhmax.nc" # Relative humidity 
path4 = r"/home/wangplll/paper/Global_HPT/data3/himean.nc" # Heat index
path5 = r"/home/wangplll/paper/Global_HPT/data3/winmean.nc" # Wind
path6 = r"/home/wangplll/paper/Global_HPT/data3/pmean.nc" # Preciptation
path7 = r"/home/wangplll/paper/Global_HPT/data3/tmax.nc"  # Tmean

path_us = r"/home/wangplll/paper/Global_HPT/station/urban_station_retained.xlsx" # Urban station
path_rs = r"/home/wangplll/paper/Global_HPT/station/rural_station_retained.xlsx"

tmean = xr.open_dataset(path1)
urban_station = pd.read_excel(path_us)
rural_station = pd.read_excel(path_rs)

# Try Northren Hemisphere
urban_station = urban_station.loc[urban_station["Lat"]>0]
rural_station = rural_station.loc[rural_station["Lat"]>0]


# In[2. Preprocessing]
tmean = tmean.rename({"__xarray_dataarray_variable__":"tmean"})
lon = tmean["lon"].values
lat = tmean["lat"].values
alt = tmean["alt"].values

ID = np.arange(len(lon))
# add new coordinate to dimention---'station'
tmean = tmean.assign_coords(ID=("station", ID))

if False:  # Relative humidity---->True else False
    tmean1 = tmean.where(tmean>0, other=0)
    tmean2 = tmean1.where(tmean<100, other=100)
    tmean = tmean2 + tmean - tmean
else:
    tmean = tmean.where(tmean>-100)
    tmean = tmean.where(tmean<100)

## Processing missing value
# Count the days of exsiting values within a year (is not NaN value)
Year_count = tmean.resample(time="1Y", skipna=True).count(dim="time")
# Count the missing values =<10% within a year, i.e., >=330 days
Year_count_gt330 = Year_count.where(Year_count>=330).count(dim="time")
# Select the maintained stations whose missing year <=18 year 
Year_count_gt30 = Year_count_gt330.where(Year_count_gt330>=30, drop=True)
tmean_sel = tmean.where(tmean.ID.isin(Year_count_gt30.ID.values), drop=True)

# Calculate Anomalies
tmean_73_03 = tmean_sel.sel(time=slice("1973-01-01", "2003-12-31"))
tmean_73_03_mean = tmean_73_03.groupby("time.dayofyear").mean(dim="time")

for yr in np.arange(1973, 2021):
    year_days = 366 if calendar.isleap(yr) else 365
    if year_days==365:
        tmean_yr = tmean_73_03_mean.drop_isel({"dayofyear":59})
    else:
        tmean_yr = tmean_73_03_mean
    
    if yr == 1973:
        tmean_concat = tmean_yr
    else:
        tmean_concat = xr.concat([tmean_concat, tmean_yr], dim="dayofyear")

tmean_concat = tmean_concat.rename({"dayofyear":"time"})
tmean_concat = tmean_concat.transpose("station", "time")
tmean_concat["time"] = tmean_sel["time"].values
tmean_concat = tmean_concat.assign_coords(ID=("station", tmean_sel.ID.values))
tmean_concat = tmean_concat.assign_coords(lon=("station", tmean_sel.lon.values))
tmean_concat = tmean_concat.assign_coords(lat=("station", tmean_sel.lat.values))
tmean_concat = tmean_concat.assign_coords(alt=("station", tmean_sel.alt.values))
tmean_concat = tmean_concat.assign_coords(id=("station", tmean_sel.id.values))

anomalies = tmean_sel - tmean_concat
tmean_sel = anomalies

# In[3. Define functions]
def linear_trend(x, y):
    # pf = np.polyfit(x, y, 1)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    
    return xr.DataArray(slope) #, xr.DataArray(intercept)


# define a function to plot filled contour map over China
def scatter_map(lon, lat, values, levels=np.arange(-1, 1.1, 0.1), cmap='RdBu_r', 
             size=5, clabel="", title="Plot", grid=False): 
    ax = plt.axes(projection=ccrs.PlateCarree()) # projection=ccrs.Mollweide()
    
    
    hc =  plt.scatter(lon, lat, s=size, c=values,  edgecolors='None', cmap=cmap)
    #edgecolors=[0.25, 0.25, 0.25]
    # colorbar
    cb = plt.colorbar(hc, orientation="vertical", fraction=0.023, pad=0.02, extend="both")
    cb.set_label(label=clabel)
    plt.clim(min(levels), max(levels))

    
    # Countries boundary    
    cn_boarder = cfeat.ShapelyFeature(Reader(r'/home/wangplll/paper/Global_HPT/spt_data/countyies_bou.shp').geometries(), ccrs.PlateCarree())
    ax.add_feature(cn_boarder, linewidth=0.25, edgecolor='k', facecolor='none')   
    
   
    # set major and minor ticks
    ax.set_xticks(np.arange(-180, 180, 60), crs=ccrs.PlateCarree())
    ax.xaxis.set_minor_locator(plt.MultipleLocator(5))
    ax.set_yticks(np.arange(-90, 90, 30), crs=ccrs.PlateCarree())
    ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    
    # show grid lines
    if grid:
        ax.grid(which='major', linestyle='--', linewidth=0.5, color=[0.25, 0.25, 0.25])
        
    
    # set axis limits
    ax.set_ylim(-90, 90)
    ax.set_xlim(-180, 180)    # whole world
    
    # set axis labels
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    
    # set figure title
    ax.set_title(title, weight='bold')
    
    # return
    return hc, ax # cb,

# In[4. Calculate annual and seasonal mean]
Monthly = tmean_sel.resample(time="1M", skipna=True).mean()   # Monthly mean
Annual = Monthly.groupby("time.year").mean() # Annual mean

annual_trend = xr.apply_ufunc(linear_trend, Annual.year, Annual['tmean'], \
                        vectorize=True, input_core_dims=[['year'], ['year']], \
                            output_core_dims=[[] for _ in range(1)])  
scatter_map(tmean_sel.lon.values,tmean_sel.lat.values,annual_trend*10, levels = np.arange(-1, 1, 0.1), cmap='RdYlBu_r',
         size=3,clabel="Temperature (¡ãC)   ", title="Trend of Annual mean AT", grid=False)
 
# Seasonal trend
for ssn in ['MAM', 'JJA', 'SON', 'DJF']:
    Seasonal_monthes = Monthly.sel(time=(Monthly.time.dt.season==ssn))
    Seasonal = Seasonal_monthes.groupby('time.year').mean()
    # Trend of season
    seasonal_trend = xr.apply_ufunc(linear_trend, Seasonal.year, Seasonal["tmean"], \
                                    vectorize=True, input_core_dims=[['year'], ['year']], \
                                    output_core_dims=[[] for _ in range(1)])
    
    fig = plt.figure(figsize=(5.8, 4.1))  ## A6=4.1x5.8, A5=5.8x8.3, A4=8.3x11.7
    scatter_map(tmean_sel.lon.values,tmean_sel.lat.values,seasonal_trend*10, levels = np.arange(-1, 1, 0.1), cmap='RdYlBu_r',
         size=3,clabel="Temperature (¡ãC)", title="Trend of {} AT".format(ssn), grid=False)


# In[5. Long term trend]
tmean_u = tmean_sel.where(tmean_sel.ID.isin(urban_station["ID"]), drop=True)
tmean_r = tmean_sel.where(tmean_sel.ID.isin(rural_station["ID"]), drop=True)

tmean_u_mt = tmean_u.resample(time="1M", skipna=True).mean() # Monthly mean
tmean_u_yr = tmean_u_mt.groupby("time.year").mean() # Annual mean
# Regional mean
tmean_u_all = tmean_u_yr.mean(dim="station")

tmean_r_mt = tmean_r.resample(time="1M", skipna=True).mean() # Monthly mean
tmean_r_yr = tmean_r_mt.groupby("time.year").mean() # Annual mean
tmean_r_all = tmean_r_yr.mean(dim="station")

plt.plot(tmean_u_all.year.values, tmean_u_all["tmean"].values)
plt.plot(tmean_r_all.year.values, tmean_r_all["tmean"].values)

# Seasonal mean
for ssn in ['MAM', 'JJA', 'SON', 'DJF']:
    tmean_u_monthes = tmean_u_mt.sel(time=(tmean_u_mt.time.dt.season==ssn))
    tmean_u_seasonal = tmean_u_monthes.groupby('time.year').mean()
    tmean_u_seasonal_mean = tmean_u_seasonal.mean(dim="station")
    
    tmean_r_monthes = tmean_r_mt.sel(time=(tmean_r_mt.time.dt.season==ssn))
    tmean_r_seasonal = tmean_r_monthes.groupby('time.year').mean()
    tmean_r_seasonal_mean = tmean_r_seasonal.mean(dim="station")
    
    xx = np.arange(1973,2021)
    yy_u = tmean_u_seasonal_mean['tmean'].values
    slope, intercept, r_value, p_value, std_err = stats.linregress(xx, yy_u)
    
    yy_r = tmean_r_seasonal_mean['tmean'].values
    slope_r, intercept_r, r_value_r, p_value_r, std_err_r = stats.linregress(xx, yy_r)
    
    fig = plt.figure(figsize=(5.8, 4.1))  ## A6=4.1x5.8, A5=5.8x8.3, A4=8.3x11.7
    plt.plot(tmean_u_seasonal_mean.year.values, tmean_u_seasonal_mean["tmean"].values)
    plt.plot(xx, xx*slope + intercept, linestyle='--', linewidth=1)
    plt.plot(xx, xx*slope_r + intercept_r, linestyle='--', linewidth=1)
    
    plt.plot(tmean_r_seasonal_mean.year.values, tmean_r_seasonal_mean["tmean"].values)
    plt.title("Trend of {} mean AT".format(ssn))
    plt.show()
    
    print(slope, slope_r)
    
    
    