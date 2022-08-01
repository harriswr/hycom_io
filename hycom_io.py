#!/bin/python

from tkinter import Y
from unittest import skip
import netCDF4 as nc
import numpy as np
import copy
import matplotlib.pyplot as plt
from matplotlib import cm, colors
#from mpl_toolkits.basemap import Basemap
import time
import math
import pandas as pd
import scipy.interpolate as spinterp 

import cartopy.crs as ccrs
import string
import iris
from iris.util import squeeze as isqueeze
import iris.plot as iplt
import iris.quickplot as qplt
import iris.analysis.cartography as iac
from mpl_toolkits.axes_grid1 import make_axes_locatable
    
class plotBathy():
    
    def __init__(self, csv_fname=None, long_lims=(-180,180), lat_lims=(-90,90)):
        if csv_fname != None:
            print(f'...Reading topology from {csv_fname}')
            start = time.time()
            self.df = pd.read_csv(csv_fname, low_memory=False)
            end = time.time()
            total_time = end-start
            print(f"...Read CSV in {total_time:.2f} seconds")
            
            self.long_lims = long_lims
            self.lat_lims = lat_lims
            
            print(f'...Generating bathymetry surface.')
            start = time.time()
            self.gen_topography_data()  
            end = time.time()
            total_time = end-start
            print(f"...Generated bathymetry surface in {total_time:.2f} seconds")
            
        else:
            print('...Error: need to download the bathymetry file from:/n\
                  https://coastwatch.pfeg.noaa.gov/erddap/griddap/usgsCeSrtm30v6.html')
    
    def gen_topography_data(self):
        topos = self.df[1:].to_numpy(dtype='float64')
        lats = topos[:,0]
        longs = topos[:,1]
        depth = topos[:,2]
        resolution = 0.00833333

        nx = (max(longs)-min(longs))//resolution
        ny = (max(lats)-min(lats))//resolution

        # Build 2 grids: One with lats and the other with lons
        self.bathy_x, self.bathy_y = np.meshgrid(np.linspace(min(longs), max(longs), num=int(nx)), \
                                     np.linspace(min(lats), max(lats), num=int(ny)))
            
        # Interpolate topo into a grid (x by y dimesions)
        self.bathy_z = spinterp.griddata((longs,lats),depth,(self.bathy_x,self.bathy_y),method='linear')

    def plot_bathy_3d(self, show=True):
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        surf = ax.plot_surface(self.bathy_x, self.bathy_y, self.bathy_z,\
                               linewidth=0, antialiased=False)
        
        
class hycom_cubes():
    
    def __init__(self, root_dir='/home/harriswr/projects/acoustics/cpies/hycom_data/NESM_tile_S_T_and_SSH/2019052012/',\
                       uv_fname='GLBy0.04_190_2019052012_t0000_uv3z-nesm.nc',\
                       st_fname='GLBy0.04_190_2019052012_t0000_ts3z-nesm.nc',\
                       ssh_fname='../ssh/GLBy0.04s_190_2019052012_t0000_ssh-nesm.nc',\
                       date = '2019052112',\
                       time = 0,
                       verbose = True):
        
        self.verbose = verbose
        self.subsets = []
        self.range_slice = []
        self.time = time
        self.date = date
        self.cubes  = iris.load(root_dir+uv_fname) + iris.load(root_dir+st_fname) + iris.load(root_dir+ssh_fname)
    
        #get rid of Tau element since we aren't using it
        i = 0
        while i < len(self.cubes):
            if self.cubes[i].name() == 'Tau':
                self.cubes.pop(i)
            else:
                self.cubes[i] = isqueeze(self.cubes[i])
                self.cubes[i].rename(self.cubes[i].name().lower().replace(' ','_').replace('-','_'))
                i += 1
                
        self.list_var()
        print('\n')
        self.get_lat_long_lims(self.cubes[0])
                
    def list_var(self):
        '''list the variables loaded from the HYCOM output'''
        
        self.variables = []
        
        if self.verbose == True:
            print('\nVariables in dataset:\n')
        for _ in self.cubes:
            if self.verbose == True:
                print(f"\t{_.name()}")
            self.variables.append(_.name())
    
    def get_var(self, var_name, cubeList = None):
        '''returns a cube of the specificed variable'''
        if cubeList == None:
            cubeList = self.cubes

        if not var_name in self.variables:
            print('ERROR: Invalid variable name.')
            self.list_var()
        else:
            var_cube = cubeList.extract(str(var_name))
            return var_cube[0]
        
    def get_lat_long_lims(self, var_cube):
        '''get the horizontal limits of the variable cube'''
        
        lat = self.get_var_array(var_cube, 'latitude')
        long = self.get_var_array(var_cube, 'longitude')        
        
        self.long_lims = (long.min(), long.max())
        self.lat_lims = (lat.min(), lat.max())
        
        if self.verbose == True:
            print(f'Longitude:\n\tMin: {self.long_lims[0]:.2f}\n\tMax: {self.long_lims[1]:.2f}\nLatitude:\n\tMin: {self.lat_lims[0]:.2f}\n\tMax: {self.lat_lims[1]:.2f}')
    
    def get_var_lims(self, var_cube, depth=None):
        '''return the minimum and maximum values of the data values contained within a cube'''
        
        a, b = var_cube.data.min(), var_cube.data.max()
        if depth == None:
            return a, b
        else:
            sliced_cube = var_cube.extract(iris.Constraint(depth=depth))
            a,b = sliced_cube.data.min(), sliced_cube.data.max()
            return a,b
        
    def get_subset(self, var_cube, lat = None, long = None):
        '''Pull depth profile at a set point in (long, lat)'''
        
        if True:
            if len(lat) == 1:
                lat = [lat+((-1)**n)*1e-4 for n in range(2)]
                long = [long+((-1)**n)*1e-4 for n in range(2)]
                
                long.sort()
                lat.sort()
                
            sliced_cube = var_cube.intersection(longitude=long,latitude=lat)
            self.subsets.append(sliced_cube)
        if False:
            print('Must specify lat and long')
    
    def get_var_array(self, var_cube, var_name):
        '''return an array of the datapoints as a numpy ndarray'''
        
        return var_cube.coords(var_name)[0].points
    
    def plot_depth_slice(self, var_cube, slice_depth=0, save_path=None, bar_range=None):
        '''plot variable at a set depth'''
        
        sliced_cube = var_cube.extract(iris.Constraint(depth=slice_depth))
        x = sliced_cube.slices(['latitude','longitude']).next()
            
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()
        var_name = var_cube.name()
        var_name_plot = string.capwords(var_name.replace('_',' '))
        plt.title(f'{var_name_plot} at {slice_depth} {str(var_cube[0].coord("depth").units)}')
        date_fmt = f'{self.date[4:6]}/{self.date[6:8]}/{self.date[0:4]}'
        iplt.citation(f'{date_fmt}: {str(self.time*100).zfill(4)}') 
        ax.gridlines(draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
            
        _ = iplt.contourf(x, 25, cmap=cm.jet, axes=ax)
            
            
        #colorbar formatting
        if bar_range == None:
            bar = plt.colorbar(orientation="horizontal")
        else:
            v = np.linspace(bar_range[0], bar_range[1], num=11)
            bar = plt.colorbar(orientation="horizontal", ticks=v)
            
        if str(var_cube[0].units) != 'unknown':
            bar.ax.set_xlabel(var_cube[0].units)
        elif 'sal' in var_name:
            bar.ax.set_xlabel('psu')
            
        if save_path != None:
            if save_path[-1] == '/':
                save_path = save_path[:-1]
            save_name = f'{save_path}/{var_name}_{slice_depth}_{self.date}_t{str(self.time).zfill(4)}.png'
            plt.savefig(save_name)
        
        plt.close()
            
    def plot_vert_slice(self, var_cube, x_axis='latitude', slice_val=None, save_path=None, bar_range=None):
        '''Plot a range-depth slice along one of the main horizontal axes (longitude or latitude)'''
        if slice_val == None:
            if x_axis == 'latitude':
                slice_val = self.lat_lims[1]
            else:
                slice_val = self.long_lims[1]
            
        if x_axis == 'longitude':
            sliced_cube = var_cube.extract(iris.Constraint(longitude=slice_val))
            x = sliced_cube.slices(['latitude','depth']).next()
            iplt.contourf(x, coords=['latitude', 'depth'],cmap='RdBu_r')
        else:
            sliced_cube = var_cube.extract(iris.Constraint(latitude=slice_val))
            x = sliced_cube.slices(['longitude','depth']).next()
            iplt.contourf(x, coords=['longitude', 'depth'],cmap=cm.jet) 
            
        var_name = var_cube.name()
        var_name_plot = string.capwords(var_name.replace('_',' '))
        date_fmt = f'{self.date[4:6]}/{self.date[6:8]}/{self.date[0:4]}'
        iplt.citation(f'{date_fmt}: {str(self.time*100).zfill(4)}') 
        plt.title(f'Slice of {var_name_plot} at {string.capwords(x_axis)} of {slice_val}')
        plt.xlabel(string.capwords(x_axis))
        plt.ylabel("Depth (m)")
        
        #iplt.show()
            
            
        #colorbar formatting
        if bar_range == None:
            bar = plt.colorbar(orientation="horizontal")
        else:
            v = np.linspace(bar_range[0], bar_range[1], num=11)
            bar = plt.colorbar(orientation="horizontal", ticks=v)
            
        if str(var_cube[0].units) != 'unknown':
            bar.ax.set_xlabel(var_cube[0].units)
        elif 'sal' in var_name:
            bar.ax.set_xlabel('psu')
                
        if save_path != None:
            if save_path[-1] == '/':
                save_path = save_path[:-1]
            save_name = f'{save_path}/{var_name}_{x_axis}_{slice_val}_{self.date}_t{str(self.time).zfill(4)}.png'
            plt.savefig(save_name)
                
        plt.close()
    
    def calc_ssp(self, salinity, temperature, depth):
        '''takes arrays of salinity temperature and depth and calculates the resulting SSP at each depth
            
            formula taken from COA page 3.
        '''
        
        return 1449.2 + 4.6*temperature - 0.055*temperature**2 + 0.00029*temperature**3 +\
            (1.34-0.01*temperature)*(salinity-35) + 0.016*depth
     
    def make_depth_profile(self, clear_subsets = True, lat=None, long=None, save_path=None, sal_lims=None, t_lims=None):
        '''Create 2D profile of salinity and temperature on one plot and sound speed on a second set of axes'''
        lat_label = str(lat)
        long_label = str(long)
        
        if lat is not tuple:
            lat = [lat+((-1)**n)*1e-4 for n in range(2)]
            long = [long+((-1)**n)*1e-4 for n in range(2)]
            
        lat.sort()
        long.sort()
        
        if clear_subsets:
            self.subsets = []
            
        self.get_subset(self.get_var('sea_water_salinity'), lat=lat, long=long)
        self.get_subset(self.get_var('sea_water_temperature'), lat=lat, long=long)
        
        
        depth = self.subsets[0].coord('depth').points
        salinity = self.subsets[0].data
        temp = self.subsets[1].data
        
        depth = depth.squeeze()
        salinity = salinity.squeeze()
        temp = temp.squeeze()
     
        
        ssp = self.calc_ssp(salinity, temp, depth)
        date_fmt = f'{self.date[4:6]}/{self.date[6:8]}/{self.date[0:4]}'
        
        
        fig, ax = plt.subplots(1,2)
        ax1 = ax[0]
        ax3 = ax[1]
        l1 = ax1.plot(salinity, depth, color='red', label='Salinity')
        ax1.set_xlabel('Salinity (psu)')
        ax1.set_xlim(sal_lims)
        ax1.set_ylabel('Depth (m)')
        ax1.invert_yaxis()
        ax2 = ax1.twiny()
        ax2.set_xlim(t_lims)
        l2 = ax2.plot(temp, depth, color = 'blue', label='Temp.')
        ax2.set_xlabel('Temperature (degC)')
        l3 = ax3.plot(ssp, depth, color = 'orange', label='Sound Speed')
        ax3.set_xlabel('Sound Speed (m/s)')
        ax3.set_ylabel('Depth (m)')
        ax3.invert_yaxis()
        
        lns = l1+l2
        labs = [l.get_label() for l in lns]
        ax1.legend(lns, labs, loc=0)

        plt.suptitle(f'Depth Profiles at ({lat_label}, {long_label})\n{date_fmt}: {str(self.time*100).zfill(4)}')
        #plt.figtext(0.7, 0.15, f'{date_fmt}: {str(self.time*100).zfill(4)}')
        fig.tight_layout()
        #plt.show()
        
        if save_path != None:
            if save_path[-1] == '/':
                save_path = save_path[:-1]
            lat_str = str(lat).replace(' ','_')
            long_str = str(long).replace(' ','_')
            save_name = f'{save_path}/ssp_{lat_str}_{long_str}_{self.date}_t{str(self.time).zfill(4)}.png'
            plt.savefig(save_name)
            
        plt.close()
    
    def get_sliced_cube(self, p1=(None,None), p2=(None, None), overwrite=False):
        ''' for each cube variable defined in the current object, define an interpolation slice between
            p1 and p2 (lat, long pairs) to the maximum depth. Then interpolate each of the variables onto that slice
            
            should probably only do this operation on a deepcopy of first cube so that can compare.'''
        if (p1 == (None, None)) or (p2 == (None,None)) or (len(p1)!=2 and len(p2)!=2):
            print('must specify latitude and longitude')
            return
        
        if self.range_slice != [] and overwrite == False:
            print('Overwriting existing range slice. Please hit overwrite.')
            return

        def lat(long, p1=p1, p2=p2):
            try:
                return ((p2[0]-p1[0])/(p2[1]-p1[1]))*(long - p1[1])+p1[0]
            except:
                print('invalid inputs for linear determination of latitude')
        #assume that p1 = a,b and p2 = c,d

        b, a = p1
        d, c = p2 

        dist = np.sqrt((d-b)**2 + (c-a)**2) * 100.
        theta = math.atan2(d-b, c-a)
        
        type(dist)
        type(theta)

        #hycom is 1/25th degree (roughly 4km)
        n = math.ceil(dist / 4.)
        dx = dist/n * np.cos(theta) / 100.
        dy = dist/n * np.sin(theta) / 100.

        #print(a,c,n,dy)

        dist_inc = np.linspace(0, dist, num=n)
        longs = [a + dx*i for i in range(n)]
        lats = [lat(along) for along in longs]

        #print(longs, lats)
        
        cube_slice_list = []
        for ai, acube in enumerate(self.cubes):
            temp_cube = acube
            
            #sample_pts = [('longitude', longs), ('latitude', lats)]
            #for ai, acube in enumerate(self.cubes):
            #    self.cubes[ai] = self.cubes[ai].interpolate(sample_pts, iris.analysis.Linear())
            interpolator = iris.analysis.Linear().interpolator(temp_cube, ['latitude','longitude'])
            slice_cubes = []
            for aj,(lat,lon) in enumerate(zip(lats,longs)):
                result = interpolator([lat,lon])
                coord = iris.coords.DimCoord(dist_inc[aj], var_name='Range', units='km')
                result.add_aux_coord(coord)
                slice_cubes.append(result)

            slice_cubes = iris.cube.CubeList(slice_cubes)
            cube_slice_list.append(slice_cubes.merge_cube())
        self.range_slice = iris.cube.CubeList(cube_slice_list)

    def plot_range_slice_2D(self, var_names, p1=(None,None), p2=(None,None), bar_range = None, save_path=None):
        '''Uses the get_sliced_cube method to develop range_depth and salinity contours'''
        if type(var_names) != list:
            n = 1
        else:
            n = len(var_names)
        
        self.get_sliced_cube(p1 = p1, p2 = p2, overwrite=True)

        fig, ax = plt.subplots(n, 1, figsize=(12,8))
        for i in range(n):
            var = var_names[i]
            if var == "ssp":
                temp = self.get_var('sea_water_temperature', cubeList = self.range_slice)
                sal = self.get_var('sea_water_salinity', cubeList = self.range_slice)

                temp_arr = temp.data
                sal_arr = sal.data
                depth_arr = temp.coord('depth').points
                range_arr = temp.coord('Range').points

                contour_var = self.calc_ssp(sal_arr, temp_arr, depth_arr)

                #iplt.contourf(ssp, coords=['latitude', 'depth'],cmap='RdBu_r')
            else:
                var_cube = self.get_var(var, cubeList=self.range_slice)
                contour_var = var_cube.data
                range_arr = var_cube.coord('Range').points
                depth_arr = var_cube.coord('depth').points

            if i == n-1:
                ax[i].set_xlabel('Range (km)')
            ax[i].set_ylabel('Depth (m)')
            ax[i].set_ylim(0, 6000)
            Range, Depth = np.meshgrid(range_arr, depth_arr)

            if contour_var.shape != Range.shape:
                contour_var = contour_var.T

            if 'sal' in var:
                unit_title = ' (psu)'
            elif 'temp' in var:
                unit_title = ' (degC)'

            if var == 'ssp':
                ax[i].text(0.7, 0.1, "Sound Speed (m/s)", transform=ax[i].transAxes)
            else:
                title_text = string.capwords(var.replace('_',' ')) + unit_title
                ax[i].text(0.7, 0.1, title_text, transform=ax[i].transAxes)

            l2 = ax[i].contourf(Range, Depth, contour_var, cmap=cm.jet)

            if bar_range == None:
                bar = plt.colorbar(l2, orientation="vertical", ax=ax[i], pad=0.01)
            else:
                v = np.linspace(bar_range[0], bar_range[1], num=11)
                bar = plt.colorbar(l2, orientation="vertical", ticks=v,ax=ax[i],pad = 0.01)

            ax[i].invert_yaxis()
            #fig.colorbar(cntr1, ax=ax1)
            point_label = f"({p1[0]}, {p1[1]}) to ({p2[0]}, {p2[1]})"
            date_fmt = f'{self.date[4:6]}/{self.date[6:8]}/{self.date[0:4]}'
            plt.suptitle(f'Depth Profiles from {point_label}\n{date_fmt}: {str(self.time*100).zfill(4)}', fontsize=20 )

            if str(var_cube[0].units) != 'unknown':
                bar.ax.set_xlabel(var_cube.units)
            elif 'sal' in var:
                bar.ax.set_xlabel('psu')
            elif 'ssp' in var:
                bar.ax.set_xlabel('m/s')

        plt.subplots_adjust(top = 0.9, bottom=0.08, hspace=0.3, wspace=0.01)
        
        if save_path != None:
            if save_path[-1] == '/':
                save_path = save_path[:-1]
            lat_str = str(lat).replace(' ','_')
            long_str = str(long).replace(' ','_')
            save_name = f'{save_path}/slice_{p1[0]},{p1[1]}_to_{p2[0]},{p2[1]}_{self.date}_t{str(self.time).zfill(4)}.png'
            plt.savefig(save_name)
            plt.close()
        else:
            plt.show()
