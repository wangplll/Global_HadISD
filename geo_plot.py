# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 14:26:09 2021

@author: wangplll
"""

# In[0. Python libs]
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.io.shapereader import Reader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy import stats
import pandas as pd
import pymannkendall as mk

plt.rcParams['savefig.dpi'] = 300 #图片像素
plt.rcParams['figure.dpi'] = 300  #分辨率


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
    cn_boarder = cfeat.ShapelyFeature(Reader(r'F:\data\geodata\region\region.shp').geometries(), ccrs.PlateCarree())
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

def box_pot(data, title="All Station"):
    plt.figure() 
    ax = plt.subplot() 
    ax.boxplot([data.sel(Pct=1).values, data.sel(Pct=2).values, data.sel(Pct=3).values,
                    data.sel(Pct=4).values, data.sel(Pct=5).values, data.sel(Pct=6).values,
                    data.sel(Pct=7).values, data.sel(Pct=8).values, data.sel(Pct=9).values, 
                    data.sel(Pct=10).values],
                   sym = 'o',            # 异常点形状
                   vert = True,          # 是否垂直
                   whis=1.5,             # IQR
                   patch_artist = True,  # 上下四分位框是否填充
                   boxprops = {'color':'black','facecolor':'pink'},
                   meanline = False, showmeans = True,  # 是否有均值线及其形状
                   showbox = True,  # 是否显示箱线 
                   notch = False)    # 中间箱体是否缺口
    
    ax.set_ylabel('Trend (events/decade)')
    ax.set_xlabel('Categories')
    ax.set_title(title)
    ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9,10])
    ax.set_xticklabels(['<P10','P10-P20','P20-P30','P30-P40','P40-P50','P50-P60','P60-P70', \
                        'P70-P80','P80-P90','>P90'], fontsize=7.5)
    plt.show()


    
    
    
    
    
    
    
    