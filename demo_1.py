#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  4 10:25:22 2021

@author: wangplll
"""
import pandas as pd
import xarray as xr

path0 = r"/home/wangplll/paper/Global_HPT/data3/tmean.nc"
tmean = xr.open_dataset(path0)
tmean = tmean.rename({'__xarray_dataarray_variable__': 'tmean'})

pd.date_range('2021-02-28', periods=100)