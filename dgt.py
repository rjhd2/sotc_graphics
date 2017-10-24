#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for Drought (DGT) section.
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
import calendar
import sys
import datetime as dt

import utils # RJHD utilities
import settings

data_loc = "/data/local/rdunn/SotC/{}/data/DGT/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)


#************************************************************************
# Anomaly Map

#************************************************************************
# Maps for selected years


cube_list = iris.load(data_loc + "BAMS.scPDSI.cru.3.25.2017.GLOBAL.annual.mean.2015.2016.nc")

cube = cube_list[0]
cube.coord('latitude').guess_bounds()
cube.coord('longitude').guess_bounds()

bounds = [-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100]

# select the year - plot 2016 (cube[1])

utils.plot_smooth_map_iris(image_loc + "DGT_{}_rel_1901-{}".format(settings.YEAR, settings.YEAR), cube[1], settings.COLOURMAP_DICT["hydrological"], bounds, "Categories relative to 1901-{} (self-calibrating PSDI)".format(settings.YEAR), cb_extra = ["Dry", "Wet"], contour = True)
utils.plot_smooth_map_iris(image_loc + "p2.1_DGT_{}_rel_1901-{}".format(settings.YEAR, settings.YEAR), cube[1], settings.COLOURMAP_DICT["hydrological"], bounds, "Categories relative to 1901-{} (self-calibrating PSDI)".format(settings.YEAR), figtext = "(q) Drought (self-calibrating PSDI)", cb_extra = ["Dry", "Wet"], contour = True, save_netcdf_filename = "{}DGT_for_NOAA_{}.nc".format(data_loc, dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y")))

# plot 2015
utils.plot_smooth_map_iris(image_loc + "DGT_{}_rel_1901-{}".format(int(settings.YEAR)-1, settings.YEAR), cube[0], settings.COLOURMAP_DICT["hydrological"], bounds, "Categories relative to 1901-{} (self-calibrating PSDI)".format(settings.YEAR), cb_extra = ["Dry", "Wet"], contour = True)


# from 2015 report 

sys.exit()

# select the last year
cube_list = []
YEARS = [1982, 1997]
for year in YEARS:
    
    loc, = np.where(cube.coord("time").points == year)
    cube_list += [cube[loc][0]]

utils.plot_smooth_map_iris_multipanel(image_loc + "DGR_{}_{}_{}".format(settings.YEAR, YEARS[0], YEARS[1]), cube_list, settings.COLOURMAP_DICT["hydrological"], bounds, "Categories relative to 1901-2015 (self-calibrating PSDI)", shape = (2, 1), title = [str(y) for y in YEARS], figtext = ["(a)","(b)"])

# select the last year
cube_list = []
YEARS = [1985, 1987]
for year in YEARS:
    
    loc, = np.where(cube.coord("time").points == year)
    cube_list += [cube[loc][0]]

utils.plot_smooth_map_iris_multipanel(image_loc + "DGR_{}_{}_{}".format(settings.YEAR, YEARS[0], YEARS[1]), cube_list, settings.COLOURMAP_DICT["hydrological"], bounds, "Categories relative to 1901-2015 (self-calibrating PSDI)", shape = (2, 1), title = [str(y) for y in YEARS], figtext = ["(a)","(b)"])



#************************************************************************
#                                 END
#************************************************************************

