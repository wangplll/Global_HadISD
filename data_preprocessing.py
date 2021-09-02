## Global climate stations data preprocessing


# In[0. Import python libs]
import xarray as xr
import numpy as np


# In[1. Data preprocessing]
## Change data variable name
def change_varname(data, varname):
    data = data.rename({"__xarray_dataarray_variable__":varname})
    return data

## Add ID 
def add_coords(data, *coords):
    data = data.assign_coords({
        "ID":("station", coords[0])})  # Only add 1 coord to dim--"station"
    return data

## Process missing data
def process_missingval(data):
    # Count the days of exsiting values within a year (is not NaN value)
    Year_count = data.resample(time="1Y", skipna=True).count(dim="time")
    # Count the missing values =<10% within a year, i.e., >=330 days
    Year_count_gt330 = Year_count.where(Year_count>=330).count(dim="time")
    # Select the maintained stations whose missing year <=18 year 
    Year_count_gt30 = Year_count_gt330.where(Year_count_gt330>=30, drop=True)
    data_sel = data.where(data.ID.isin(Year_count_gt30.ID.values), drop=True)
    return data_sel


# In[main]
if __name__ == '__main__':
    # Define dataset name
    dataname = "tmean"
    filepath = r"F:\data\data3" + "\\" + dataname + ".nc"
    # Import dataset
    data = xr.open_dataset(filepath)
    # Change varible name
    data = change_varname(data, dataname)
    
    # Add ID
    ID = np.arange(len(data.lon.values))
    data = add_coords(data, ID)
    
    # Processing missing values
    data_sel = process_missingval(data)