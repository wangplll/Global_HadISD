import xarray as xr
import pandas as pd
import numpy as np

path0 = r"F:\data\data3\tmean.nc"

tmean = xr.open_dataset(path0)

print(tmean)