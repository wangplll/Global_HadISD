#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  3 10:49:44 2021

@author: wangplll
"""

import xarray as xr
import pandas as pd
import numpy as np
import os

path0 = r"/home/wangplll/paper/Global_HPT/data"
path1 = r"/home/wangplll/paper/Global_HPT/data1"
path2 = r"/home/wangplll/paper/Global_HPT/data2"
path3 = r"/home/wangplll/paper/Global_HPT/data3/"

filelist0 = sorted(os.listdir(path0))
for fl0 in filelist0:
    pt0 = path0 + "/" +fl0
    file0 = xr.open_dataset(pt0)
    if len(file0.station.values)>4988:
        file0 = file0.sel(station=slice(0, 4988))
        
    filelist1 = sorted(os.listdir(path1))
    for fl1 in filelist1:
        pt1 = path1 + "/" + fl1
        file1 = xr.open_dataset(pt1)
        if len(file1.station.values)>424:
            file1 = file1.sel(station=slice(0, 424))
        if fl0==fl1:
            file = xr.concat([file0,file1], dim="station")
        
        filelist2 = sorted(os.listdir(path2))
        for fl2 in filelist2:
            pt2 = path2 + "/" + fl2
            file2 = xr.open_dataset(pt2)
            if fl1==fl2 and fl0==fl2:
                file_ = xr.concat([file, file2], dim="station")
                file_.to_netcdf(path3+fl2)
                print(file_)
           
            
        
    