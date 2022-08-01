#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 09:36:44 2022

@author: harriswr
"""

import urllib
from io import StringIO
import csv
import numpy as np
import scipy.interpolate as spinterp
import matplotlib.pyplot as plt
from matplotlib import cm
import time
import pickle as pickle
#from mpl_toolkits.basemap import basemap

#controls
download_URL = False

#define min and max long/lat
minlat = 35
maxlat = 40
minlon = -65
maxlon = -60

#read data from noaa
#https://coastwatch.pfeg.noaa.gov/erddap/griddap/usgsCeSrtm30v6.csv?topo%5B(34):1:(42)%5D%5B(-58):1:(-68)%5D
#response = urllib.request.urlopen('http://coastwatch.pfeg.noaa.gov/erddap/griddap/usgsCeSrtm30v6.csv?topo[('+str(maxlat)+'):1:('+str(minlat)+')][('+str(minlon)+'):1:('+str(maxlon)+')]')
#response = urllib.request.urlopen('https://coastwatch.pfeg.noaa.gov/erddap/griddap/usgsCeSrtm30v6.csv?topo%5B(34):1:(42)%5D%5B(-58):1:(-68)%5D')
#data = StringIO.StringIO(response.read())

import pandas as pd

if download_URL:
    url = f'https://coastwatch.pfeg.noaa.gov/erddap/griddap/usgsCeSrtm30v6.csv?topo%5B({minlat}):50:({maxlat})%5D%5B({maxlon}):50:({minlon})%5D'
    start = time.time()
    print("...Reading data from NOAA.gov")
    df = pd.read_csv(url, low_memory=False)
    end = time.time()
    total_time = end-start
    print(f"...Read NOAA data in {total_time} seconds")
else:
    noaa_root = '/home/harriswr/projects/acoustics/cpies/noaa/'
    fn = noaa_root+'lat_35to40_long_-65to-60.csv'
    print(f'...Reading topology from {fn}')
    start = time.time()
    df = pd.read_csv(fn, low_memory=False)
    end = time.time()
    total_time = end-start
    print(f"...Read CSV in {total_time} seconds")
    
    
topos = df[1:].to_numpy(dtype='float64')
lats = topos[:,0]
longs = topos[:,1]
depth = topos[:,2]
resolution = 0.00833333

nx = (max(longs)-min(longs))//resolution
ny = (max(lats)-min(lats))//resolution

# Build 2 grids: One with lats and the other with lons
grid_x, grid_y = np.meshgrid(np.linspace(min(longs), max(longs), num=int(nx)), \
                             np.linspace(min(lats), max(lats), num=int(ny)))
    
# Interpolate topo into a grid (x by y dimesions)
grid_z = spinterp.griddata((longs,lats),depth,(grid_x,grid_y),method='linear')


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(grid_x, grid_y, grid_z,\
                       linewidth=0, antialiased=False)
plt.show()