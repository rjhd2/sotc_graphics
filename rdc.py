#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for River Discharge (and Runoff) (RDC) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev::                                          $:  Revision of last commit
# $Author::                                       $:  Author of last commit
# $Date::                                         $:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import matplotlib.cm as mpl_cm
import matplotlib as mpl

import iris

import utils # RJHD utilities
import settings


data_loc = "/data/local/rdunn/SotC/{}/data/RDC/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

COORD_DICT = {"lat" : "latitude", "lon" : "longitude"}

def fix_coords(cube):

    for coord in ["lat", "lon"]:

        cube.coord(coord).guess_bounds()
        cube.coord(coord).units = "degrees"  
        cube.coord(coord).standard_name = COORD_DICT[coord]

    return cube # fix_coords



#************************************************************************
# Discharge Map

cube_list = iris.load(data_loc + "discharge0.5.anomaly{}.SOC15.JRA55_GPCCLW90.1958-{}.nc".format(settings.YEAR, settings.YEAR))

cube = fix_coords(cube_list[0])

#cube.data = np.ma.masked_where(np.logical_and(cube.data < 20, cube.data > -20), cube.data)
mask = cube.data.mask

bounds = [-10000, -1000, -500, -100, -50, -25, 25, 50, 100, 500, 1000, 10000]

utils.plot_smooth_map_iris(image_loc + "RDC_discharge_{}_jra55".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1958-2015 (m"+r'$^{3}$'+" s"+r'$^{-1}$'+")")
utils.plot_smooth_map_iris(image_loc + "p2.1_RDC_discharge_{}_jra55".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1958-2015 (m"+r'$^{3}$'+" s"+r'$^{-1}$'+")", figtext = "(k) River Discharge")

#************************************************************************
# Runoff Map

cube_list = iris.load(data_loc + "runoff.anomaly{}.SOC15.JRA55_GPCCLW90.1958-{}.nc".format(settings.YEAR, settings.YEAR))

cube = fix_coords(cube_list[0])
cube.data = np.ma.array(cube.data)
cube.data.mask= mask

bounds = [-10000, -500, -250, -100, -50, -25, 25, 50, 100, 250, 500, 10000]

utils.plot_smooth_map_iris(image_loc + "RDC_runoff_{}_jra55".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1958-2015 (mm yr"+r'$^{-1}$'+")")
utils.plot_smooth_map_iris(image_loc + "p2.1_RDC_runoff_{}_jra55".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1958-2015 (mm yr"+r'$^{-1}$'+")", figtext = "(j) Runoff")


#************************************************************************
#                                 END
#************************************************************************
