#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for River Discharge (and Runoff) (RDC) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 28                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2020-04-09 11:37:08 +0100 (Thu, 09 Apr #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
from __future__ import absolute_import
from __future__ import print_function

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

    for coord in ["lats", "lons"]:

        cube.coord(coord).guess_bounds()
        cube.coord(coord).units = "degrees"  
        cube.coord(coord).standard_name = COORD_DICT[coord]

    return cube # fix_coords



#************************************************************************
# Discharge Map

cube_list = iris.load(data_loc + "discharge.{}.nc".format(settings.YEAR))

cube = fix_coords(cube_list[0])

#cube.data = np.ma.masked_where(np.logical_and(cube.data < 20, cube.data > -20), cube.data)
mask = cube.data.mask

bounds = [-10000, -1000, -500, -100, -50, -25, 25, 50, 100, 500, 1000, 10000]

utils.plot_smooth_map_iris(image_loc + "RDC_discharge_{}_jra55".format(settings.YEAR), cube[0], settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1958-{} (m".format(int(settings.YEAR)-1)+r'$^{3}$'+" s"+r'$^{-1}$'+")")
utils.plot_smooth_map_iris(image_loc + "p2.1_RDC_discharge_{}_jra55".format(settings.YEAR), cube[0], settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1958-{} (m".format(int(settings.YEAR)-1)+r'$^{3}$'+" s"+r'$^{-1}$'+")", figtext="(o) River Discharge")

#************************************************************************
# Runoff Map

cube_list = iris.load(data_loc + "runoff.{}.nc".format(settings.YEAR))

cube = fix_coords(cube_list[0])
cube.data = np.ma.array(cube.data)
cube.data.mask = mask

bounds = [-10000, -500, -250, -100, -50, -25, 25, 50, 100, 250, 500, 10000]

utils.plot_smooth_map_iris(image_loc + "RDC_runoff_{}_jra55".format(settings.YEAR), cube[0], settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1958-{} (mm yr".format(int(settings.YEAR)-1)+r'$^{-1}$'+")")
utils.plot_smooth_map_iris(image_loc + "p2.1_RDC_runoff_{}_jra55".format(settings.YEAR), cube[0], settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1958-{} (mm yr".format(int(settings.YEAR)-1)+r'$^{-1}$'+")", figtext="(p) Runoff")


#************************************************************************
#                                 END
#************************************************************************
