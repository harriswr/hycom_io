#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 13:47:13 2022

@author: harriswr
"""
import os
import sys
sys.path.append(r'/home/harriswr/projects/acoustics/cpies/scripts/')

import hycom_io as hi
import numpy as np

start_dir = '/home/harriswr/projects/acoustics/cpies/hycom_data/NESM_tile_S_T_and_SSH/'
image_dir = '/home/harriswr/projects/acoustics/cpies/images/ssp_40_-60/'
try:
    os.mkdir(image_dir)
except:
    print('image_dir already exists.')

plot_variable = 'sea_water_temperature'
depth = 0

initial = 0

for adir in os.listdir(start_dir):#[:20]:
    if adir[0] == '2':
        time = 0
        while time < 24:
            print(f'...Working on {adir} at time {time}')
            uv_name = f'GLBy0.04_190_{adir}_t{str(time).zfill(4)}_uv3z-nesm.nc'
            ts_name = f'GLBy0.04_190_{adir}_t{str(time).zfill(4)}_ts3z-nesm.nc'
            ssh_name = f'../ssh/GLBy0.04s_190_{adir}_t{str(time).zfill(4)}_ssh-nesm.nc'
            data = hi.hycom_cubes(root_dir=f'{start_dir}/{adir}/',\
                                  uv_fname = uv_name,\
                                  st_fname = ts_name,\
                                  ssh_fname = ssh_name,\
                                  date = adir,\
                                  time = time, verbose=False)
            if initial == 0:
                var_min, var_max = data.get_var_lims(data.get_var(plot_variable))
                sal_min, sal_max = data.get_var_lims(data.get_var('sea_water_salinity'))
                t_min, t_max = data.get_var_lims(data.get_var('sea_water_temperature'))
            #data.plot_depth_slice(plot_variable, slice_depth=depth, save_path=image_dir,\
            #                      bar_range=(var_min, var_max))
            data.make_depth_profile(lat=40, long=-60, save_path=image_dir, sal_lims=(30., sal_max),\
                                    t_lims=(t_min, t_max))
            time += 1
            initial += 1
