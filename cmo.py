#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for Carbon Monoxide (CMO) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 26                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2019-04-17 15:34:18 +0100 (Wed, 17 Apr #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import iris

import utils # RJHD utilities
import settings


data_loc = "{}/{}/data/CMO/".format(settings.ROOTLOC, settings.YEAR)
reanalysis_loc = "{}/{}/data/RNL/".format(settings.ROOTLOC, settings.YEAR)
image_loc = "{}/{}/images/".format(settings.ROOTLOC, settings.YEAR)

LW = 3

#************************************************************************
def read_ts(filename):
    '''
    Read the timeseries and return a Timeseries object

    :param str filename: file to read

    :returns: Timeseries object
    '''

    indata = np.genfromtxt(filename, dtype=(float), skip_header=1)

    year = indata[:, 0]
    month = indata[:, ]

    times = year + (month-1)/12.

    return utils.Timeseries("CO burden", times, indata[:, 2]) # read_ts

#************************************************************************
# Timeseries 

print("missing linear trend")
# cmo_monthly = read_ts(data_loc + "co_global_tg_mm_cams.txt")

# minor_tick_interval = 1
# minorLocator = MultipleLocator(minor_tick_interval)
# fig = plt.figure(figsize=(10, 8))
# ax = plt.axes([0.13, 0.07, 0.75, 0.86])

# plt.plot(cmo_monthly.times, cmo_monthly.data, 'k', ls='-', lw=LW)

# ax.set_xlim([None, int(settings.YEAR)+2])

# ax.xaxis.set_minor_locator(minorLocator)
# utils.thicken_panel_border(ax)
# for tick in ax.yaxis.get_major_ticks():
#     tick.label.set_fontsize(settings.FONTSIZE)
# for tick in ax.xaxis.get_major_ticks():
#     tick.label.set_fontsize(settings.FONTSIZE)

# fig.text(0.03, 0.5, "Tg", va='center', rotation='vertical', fontsize=settings.FONTSIZE)

# plt.savefig(image_loc + "CMO_ts{}".format(settings.OUTFMT))
# plt.close()

#************************************************************************
# Global Map
seasonal_list = []

cube_list = iris.load(data_loc + "TCCO_ANO_YEAR_JUL_mean_{}.nc".format(settings.YEAR))
names = np.array([cube.var_name for cube in cube_list])

selected_cube, = np.where(names == "tcco_ano")

cube = cube_list[selected_cube[0]]
cube.coord('latitude').guess_bounds()
cube.coord('longitude').guess_bounds()

bounds = [-100, -20, -15, -10, -5, 0, 5, 10, 15, 20, 100]

utils.plot_smooth_map_iris(image_loc + "CMO_{}_anoms".format(settings.YEAR), cube, settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-{} (%)".format(settings.YEAR[2:]))
utils.plot_smooth_map_iris(image_loc + "p2.1_CMO_{}_anoms".format(settings.YEAR), cube, settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-{} (%)".format(settings.YEAR[2:]), figtext="(ac) Carbon Monoxide")

# # Global Map - Jan-Jun

# selected_cube, = np.where(names == "tcco_ano_jan_jun")

# cube = cube_list[selected_cube[0]]
# cube.coord('latitude').guess_bounds()
# cube.coord('longitude').guess_bounds()

# seasonal_list += [cube]

# bounds = [-100, -20, -15, -10, -5, 0, 5, 10, 15, 20, 100]

# utils.plot_smooth_map_iris(image_loc + "CMO_{}_Jan-Jun_anoms".format(settings.YEAR), cube, settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-{} (%)".format(settings.YEAR[2:]), title="January - June {}".format(settings.YEAR))

# Global Map - Jul-Dec

selected_cube, = np.where(names == "tcco_ano_jul_sep")

cube = cube_list[selected_cube[0]]
cube.coord('latitude').guess_bounds()
cube.coord('longitude').guess_bounds()

seasonal_list += [cube]

bounds = [-100, -20, -15, -10, -5, 0, 5, 10, 15, 20, 100]

utils.plot_smooth_map_iris(image_loc + "CMO_{}_Jul_Sep_anoms".format(settings.YEAR), cube, settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-{} (%)".format(settings.YEAR[2:]), title="July - Sept {}".format(settings.YEAR))

#utils.plot_smooth_map_iris_multipanel(image_loc + "CMO_{}_season_anoms".format(settings.YEAR), seasonal_list, settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-{} (%)".format(settings.YEAR[2:]), shape=(2,1), title=["January - June {}".format(settings.YEAR), "July - December {}".format(settings.YEAR)], figtext=["(a)","(b)"])
#************************************************************************
# Trend Map

# cube_list = iris.load(data_loc + "SOTC_2015_CO_Trends_Map_Data.nc", "relative total column carbon monoxide linear trend 2003 2015 ") # final space in name necessary

# cube = cube_list[0]
# cube.coord('latitude').guess_bounds()
# cube.coord('longitude').guess_bounds()

# bounds = [-100, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 100]
# print("add zero line")
# utils.plot_smooth_map_iris(image_loc + "CMO_{}_trend".format(settings.YEAR), cube, settings.COLOURMAP_DICT["composition"], bounds, "Trend over 2003-15 (% yr"+r'$^{-1}$'+")")


#************************************************************************
#                                 END
#************************************************************************
