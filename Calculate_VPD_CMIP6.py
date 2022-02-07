# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 10:52:00 2021

@author: wangplll
"""
#%% Import libs
import xarray as xr
import os

#%% Import dataset
fp = r"F:\data\CMIP6\VPD_recalculate\Add_coords_new\new"
model_path = os.listdir(fp)

for mp in model_path:
    model_dir = fp +"\\" + mp
    run_path = os.listdir(model_dir)
    for rp in run_path:
        run_dir = model_dir + "\\" + rp
        data_path = os.listdir(run_dir)
        for dp in range(int(len(data_path)/2)):
            satvpr_dir = run_dir + "\\" + data_path[dp+int(len(data_path)/2)]
            satvpr = xr.open_dataset(satvpr_dir)
            
            actvpr_dir = run_dir + "\\" + data_path[dp]
            actvpr = xr.open_dataset(actvpr_dir)
            
            if (data_path[dp+int(len(data_path)/2)].split("_")[-1]==data_path[dp].split("_")[-1]):
                vpd = satvpr["satvpr"] - actvpr["actvpr"]
                vpd.to_netcdf(r"F:\data\CMIP6\VPD_recalculate\recalculate_vpd\new" + "\\" + mp  + "_" + rp +"_vpd_" + data_path[dp].split("_")[-1])
            else:
                print("Check the dataset Please!!!")
          
# piControl
fp = r"F:\data\CMIP6\VPD_recalculate\Add_coords_new\piControl"
model_path = os.listdir(fp)

for mp in model_path:
    data_dir = fp + "\\" + mp
    data_list = os.listdir(data_dir)
    for dl in range(int(len(data_list)/2)):
        satvpr_dir = data_dir + "\\" + data_list[dl+int(len(data_list)/2)]
        satvpr = xr.open_dataset(satvpr_dir)
            
        actvpr_dir = data_dir + "\\" + data_list[dl]
        actvpr = xr.open_dataset(actvpr_dir)
        
        vpd = satvpr["satvpr"] - actvpr["actvpr"]
        vpd.to_netcdf(r"F:\data\CMIP6\VPD_recalculate\recalculate_vpd\piControl" + "\\" + mp  +"_vpd_" + data_list[dl].split("_")[-1])
       