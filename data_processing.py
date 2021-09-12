# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 17:47:16 2021

@author: wangplll
"""

# In[0. Import python libs]
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
import pymannkendall as mk
import geo_plot as gp


# Calculate linear trend and its significance
def linear_trend(x, y):
    """x: x-axis data, timing data(1973-2020)
       y: actual dataset 
    """
    mask = (~(np.isnan(x))) & (~(np.isnan(y)))
    
    # calculate the trend with >20% valid data
    if len(x[mask]) <= 0.2 * len(x):
        return np.nan, np.nan 
    else:
        slope, intercept, r_value, p_value, std_err = stats.linregress(x[mask], y[mask])
        result = mk.hamed_rao_modification_test(y[mask])  # mk test
        p_value1 = result.p
        return slope, p_value1

def pre_percentile(pre, pmin, pmax, p):
    # Search for all days during  pmin<pre<pmax
    """pre: present data
       pmin: minimum percentile
       pmax: maximum percentile
       p: judging pmin equals 0 or not
       """
    if p == 0:
        idx_pre = (pre>=pmin) & (pre<=pmax)  # threshold value
    else:
        idx_pre = (pre>pmin) & (pre<=pmax)
    return sum(idx_pre)


def sel_pct_trend(ds_ref, ds):  # trend of select percentile data 
    """pmin: min percentile,
       pmax: max percentile,
       ds_ref: reference data,
       ds: original data """
    
    for p in range(0, 10):
        pmin = 0.1*p
        pmax = pmin + 0.1
        print(pmin, pmax)
        # Caluculate min/max percentile of reference years
        ds_ref_pmin = ds_ref.quantile(pmin, dim='time') 
        ds_ref_pmax = ds_ref.quantile(pmax, dim='time')
        
        pre = xr.apply_ufunc(pre_percentile,
                                     ds.groupby("time.year"), 
                                     ds_ref_pmin,
                                     ds_ref_pmax,
                                     p, 
                                     vectorize=True, 
                                     input_core_dims=[['time'], [],[],[]], 
                                     output_core_dims=[[] for _ in range(1)],)
        
        pre_trend = xr.apply_ufunc(linear_trend, 
                                   pre.year, 
                                   pre['tmean'],
                                   vectorize=True, 
                                   input_core_dims=[['year'], ['year']],
                                   output_core_dims= [[] for _ in range(2)], )
        
        # Scatter
        pre_annual_trend = 10 * pre_trend[0]
        plt.figure(figsize=(5.8, 4.1))  ## A6=4.1x5.8, A5=5.8x8.3, A4=8.3x11.7
        gp.scatter_map(pre_annual_trend.lon, pre_annual_trend.lat,
                    pre_annual_trend, 
                    size= 3,
                    levels=np.linspace(-0.3, 0.3, 100), cmap='RdYlBu_r',
                    clabel=" (events per decade)",
                    title="Precipitation Frequency Trend" + " | P"  + str(int(100*pmin)) + '-P' + str(int(100*pmax)))
        plt.show()
        
        
        Pct = p + 1
        pre_trend = pre_trend[0].assign_coords(Pct=Pct)
        if p == 0:
            pre_trend_concat = pre_trend
        else:
            pre_trend_concat = xr.concat([pre_trend_concat, pre_trend], dim="Pct")

    return pre_trend_concat



   
    
    