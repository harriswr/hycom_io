#!/bin/python

exec(open('./hycom_io.py').read())

####
x = hycom_cubes()
x.cubes
p1 = (40,-62)
p2 = (42, -60)
x.get_sliced_cube(p1=p1,p2=p2)
x.plot_range_slice_2D(['sea_water_temperature', 'sea_water_salinity', 'ssp'], p1=p1, p2=p2)
#x.plot_vert_slice('sea_water_temperature', x_axis='longitude', slice_val=-65, save_path='/home/harriswr/projects/acoustics/cpies/images/')
#x.plot_depth_slice('sea_water_temperature', 3000, save_path='/home/harriswr/projects/acoustics/cpies/images/')
