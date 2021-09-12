# -*- coding: utf-8 -*-
"""
Created on Sat Sep  4 15:31:28 2021

@author: wangplll
"""

# In[0. Import python libs]
import data_preprocessing as dpp
import data_processing as dp
import geoplot as gp
import xarray as xr

# In[1. Data preprocessing]
# Data preprocess and save in disk
""" 
#Dataset name
dataname = "tmean"
filepath = r"F:\data\data3" + "\\" + dataname + ".nc"
# Import dataset
data = xr.open_dataset(filepath)
# Change varible name
data = dpp.change_varname(data, dataname)

# Add ID
ID = np.arange(len(data.lon.values))
data = dpp.add_coords(data, ID)

# Processing missing values
data_sel = dpp.process_missingval(data)

# Save data
outputpath = r"F:\data\data3\sel_data"
outname = dataname + "_sel.nc"
dpp.save_data(data_sel, outputpath, outname)
"""

# In[2. Data processing]
# Import data
filepath1 = r"F:\data\data3\sel_data\tmean_sel.nc"
tmean_sel = xr.open_dataset(filepath1)
ds_ = tmean_sel
ds_clim= ds_.sel(time=slice('1973-01-01','2003-12-31'))
ds_clim = ds_clim.sel(time = (ds_clim['time.month']>=6)&(ds_clim['time.month']<=8))  # JJA
ds_clim = ds_clim.where(ds_clim>=0.1)

ds_JJA = ds_.sel(time = (ds_['time.month']>=6)&(ds_['time.month']<=8))

pre_trend_concat = dp.sel_pct_trend(ds_clim, ds_JJA)


# In[3. Plot]
## All station
pre_trend_all = pre_trend_concat*10
gp.box_pot(pre_trend_all, title="All Station")    