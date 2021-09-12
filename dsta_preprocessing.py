import xarray as xr
import numpy as np
import os


path0 = r"/home/wangplll/paper/Global_HPT/data3/"
filelist = os.listdir(path0)

for fn in filelist:
    path1 = path0 + fn
    print (path1)
    data = xr.open_dataset(path1)
    
    data = data.rename({"__xarray_dataarray_variable__":fn[:-3]})
    lon = data["lon"].values
    lat = data["lat"].values
    alt = data["alt"].values
    
    ID = np.arange(len(lon))
    # add new coordinate to dimention---'station'
    data = data.assign_coords(ID=("station", ID))
    
    data.to_netcdf(r"/home/wangplll/paper/Global_HPT/data3_preprocessing/{}".format(fn))


