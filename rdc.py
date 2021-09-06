#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for River Discharge (and Runoff) (RDC) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 30                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2021-06-15 10:41:02 +0100 (Tue, 15 Jun #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import netCDF4 as ncdf
import numpy as np
import iris

import utils # RJHD utilities
import settings


data_loc = "{}/{}/data/RDC/".format(settings.ROOTLOC, settings.YEAR)
reanalysis_loc = "{}/{}/data/RNL/".format(settings.ROOTLOC, settings.YEAR)
image_loc = "{}/{}/images/".format(settings.ROOTLOC, settings.YEAR)

COORD_DICT = {"lats" : "latitude", "lons" : "longitude"}

#************************************************************************
def fix_coords(cube):

    for coord in ["latitude", "longitude"]:

        cube.coord(coord).guess_bounds()
        cube.coord(coord).units = "degrees"  
#        cube.coord(coord).standard_name = COORD_DICT[coord]

    return cube # fix_coords



#************************************************************************
# Discharge Map

#cube_list = iris.load(data_loc + "watflw_1981-2010_{}.nc".format(settings.YEAR))
# cube = fix_coords(cube_list[0])[0]
# cube.data = np.ma.masked_where(cube.data == 0, cube.data)
# mask = cube.data.mask

# 2020 build via netcdf as for some reason resetting the plot extent wiped the map
ncfile = ncdf.Dataset(data_loc + "watflw_1981-2010_{}.nc".format(settings.YEAR))

var=ncfile.variables["watflw"][:] # this is a masked array
nlons = ncfile.variables["lon"][:]
nlats = ncfile.variables["lat"][:]

cube = utils.make_iris_cube_2d(var[0], nlats, nlons, "RDC", "m3/s")

bounds = [-10000, -1000, -500, -100, -50, -25, 25, 50, 100, 500, 1000, 10000]

utils.plot_smooth_map_iris(image_loc + "RDC_discharge_{}_jra55".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (m"+r'$^{3}$'+" s"+r'$^{-1}$'+")")
utils.plot_smooth_map_iris(image_loc + "p2.1_RDC_discharge_{}_jra55".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (m"+r'$^{3}$'+" s"+r'$^{-1}$'+")", figtext="(p) River Discharge")


#************************************************************************
# Runoff Map

# cube_list = iris.load(data_loc + "runoffgw_1981-2010_{}.nc".format(settings.YEAR))

# cube = fix_coords(cube_list[0])
# cube.data = np.ma.array(cube.data)
# cube.data.mask = mask

# 2020 build via netcdf as for some reason resetting the plot extent wiped the map
ncfile = ncdf.Dataset(data_loc + "runoffgw_1981-2010_{}.nc".format(settings.YEAR))

var=ncfile.variables["runoffgw"][:] # this is a masked array
nlons = ncfile.variables["lon"][:]
nlats = ncfile.variables["lat"][:]

cube = utils.make_iris_cube_2d(var[0], nlats, nlons, "RDC", "mm/yr")
bounds = [-10000, -500, -250, -100, -50, -25, 25, 50, 100, 250, 500, 10000]

utils.plot_smooth_map_iris(image_loc + "RDC_runoff_{}_jra55".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (mm yr"+r'$^{-1}$'+")")
utils.plot_smooth_map_iris(image_loc + "p2.1_RDC_runoff_{}_jra55".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (mm yr"+r'$^{-1}$'+")", figtext="(o) Runoff")


#************************************************************************
#                                 END
#************************************************************************
