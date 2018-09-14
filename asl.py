#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for Aerosols (ASL) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 22                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2018-04-06 15:34:21 +0100 (Fri, 06 Apr #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import matplotlib.cm as mpl_cm
import matplotlib as mpl

import iris
import iris.quickplot as qplt
import cartopy

import calendar

import utils # RJHD utilities
import settings

data_loc = "/data/local/rdunn/SotC/{}/data/ASL/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

LW = 3
LEGEND_LOC = "upper left"

#************************************************************************
def read_data(filename, name):
    """
    Read data file and returns Timeseries
    
    :param str filename: file to process
    :param str name: name to give Timeseries

    :returns: Timeseries object
    """

    indata = np.genfromtxt(filename, dtype = (str))

    date = indata[:,0]
    data = indata[:,1].astype(float)

    month = [d[5:7] for d in date]
    year = [d[:4] for d in date]

    times = np.array(year).astype(int) + (np.array(month).astype(int) - 1)/12.

    return utils.Timeseries(name, times, data) # read_data


#************************************************************************
# Timeseries plot

monthly = read_data(data_loc + "monthmean", "AOD monthly")
annual = read_data(data_loc + "yearmean", "AOD annual")

minor_tick_interval = 1
minorLocator = MultipleLocator(minor_tick_interval)
COLOURS = settings.COLOURS["composition"]

fig = plt.figure(figsize = (10,5))
ax = plt.axes([0.13, 0.07, 0.75, 0.86])

plt.plot(monthly.times, monthly.data,  COLOURS[monthly.name], ls = '-', label = monthly.name, lw = LW)
plt.plot(annual.times, annual.data,  COLOURS[annual.name], ls = '-', label = annual.name, lw = LW)

ax.legend(loc = LEGEND_LOC, ncol = 1, frameon = False, prop = {'size':settings.LEGEND_FONTSIZE}, labelspacing = 0.1, columnspacing = 0.5)

# prettify
ax.xaxis.set_minor_locator(minorLocator)
utils.thicken_panel_border(ax)

fig.text(0.03, 0.5, "AOD", va='center', rotation='vertical', fontsize = settings.FONTSIZE)

plt.xlim([2002.5, int(settings.YEAR)+1.5])
plt.ylim([0.09, 0.23])
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(settings.FONTSIZE)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(settings.FONTSIZE)

plt.savefig(image_loc + "ASL_ts{}".format(settings.OUTFMT))
plt.close()

#************************************************************************
# Global maps

bounds = np.array([-10, -0.14, -0.10, -0.06, -0.04, -0.02, 0.02, 0.04, 0.06, 0.10, 0.14, 10])

# Biomass Burning
BB = iris.load(data_loc + "diffbb{}.nc".format(settings.YEAR))

utils.plot_smooth_map_iris(image_loc + "ASL_BB_anomalies_{}".format(settings.YEAR), BB[0][0], settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-20{} (AOD)".format(int(settings.YEAR[2:])-1))
utils.plot_smooth_map_iris(image_loc + "p2.1_ASL_BB_anomalies_{}".format(settings.YEAR), BB[0][0], settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-20{} (AOD)".format(int(settings.YEAR[2:])-1), figtext = "(ab) Biomass Burning Aerosol")

# Dust
dust = iris.load(data_loc + "diffdust{}.nc".format(settings.YEAR))

utils.plot_smooth_map_iris(image_loc + "ASL_dust_anomalies_{}".format(settings.YEAR), dust[0][0], settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-20{} (AOD)".format(int(settings.YEAR[2:])-1))
utils.plot_smooth_map_iris(image_loc + "p2.1_ASL_dust_anomalies_{}".format(settings.YEAR), dust[0][0], settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-20{} (AOD)".format(int(settings.YEAR[2:])-1), figtext = "(aa) Dust Aerosol")

# total
total = iris.load(data_loc + "difftotal{}.nc".format(settings.YEAR))

utils.plot_smooth_map_iris(image_loc + "ASL_total_anomalies_{}".format(settings.YEAR), total[0][0], settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-{} (AOD)".format(settings.YEAR[2:]))
utils.plot_smooth_map_iris(image_loc + "p2.1_ASL_total_anomalies_{}".format(settings.YEAR), total[0][0], settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-20{} (AOD)".format(int(settings.YEAR[2:])-1), figtext = "(z) Total Aerosol")


#************************************************************************
# Total, Trend and Event map

# as one colorbar per map, don't use the utils Iris routine

fig = plt.figure(figsize = (7,10))

#total = iris.load(data_loc + "Aerosol_Fig2a_Total_Averages.nc")
#trend = iris.load(data_loc + "Aerosol_Fig2b_Total_Trends.nc")
#event = iris.load(data_loc + "Aerosol_Fig2c_ExtremeDayCounts_{}.nc".format(settings.YEAR))

total = iris.load(data_loc + "total_meanfinal_{}.nc".format(settings.YEAR))
trend = iris.load(data_loc + "sig.nc")

bounds_a = np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0])
bounds_b = np.array([-10, -0.050, -0.010, -0.005, -0.002, -0.001, 0.001, 0.002, 0.005, 0.010, 0.050, 10])

all_cubes = [total, trend]
all_bounds = [bounds_a, bounds_b]

# do colourmaps by hand
cmap = [plt.cm.YlOrBr,settings.COLOURMAP_DICT["composition"]]
PLOTLABELS = ["(a)", "(b)"]
cb_label = ["AOD","AOD yr"+r'$^{-1}$']

# spin through axes
for a in range(2):  

    norm = mpl.cm.colors.BoundaryNorm(all_bounds[a], cmap[a].N)
    
    ax = plt.subplot(2, 1, a+1, projection=cartopy.crs.Robinson())
    
    ax.gridlines() #draw_labels=True)
    ax.add_feature(cartopy.feature.LAND, zorder = 0, facecolor = "0.9", edgecolor = "k")
    ax.coastlines()
        
    ext = ax.get_extent() # save the original extent

    mesh = iris.plot.pcolormesh(all_cubes[a][0][0], cmap = cmap[a], norm = norm, axes = ax)
 
    ax.set_extent(ext, ax.projection) # fix the extent change from colormesh
    ax.text(0.0, 1.0, PLOTLABELS[a], fontsize = settings.FONTSIZE * 0.8, transform=ax.transAxes)

    # sort the colourbar
    cb = fig.colorbar(mesh, ax = ax, orientation = 'horizontal', \
                          ticks = all_bounds[a][1:-1], label = cb_label[a], drawedges=True, pad = 0.05, fraction = 0.07, aspect = 30)
    cb.set_ticklabels(["{:g}".format(b) for b in all_bounds[a][1:-1]])
    cb.outline.set_linewidth(2)
    cb.dividers.set_color('k')
    cb.dividers.set_linewidth(2)

fig.subplots_adjust(bottom=0.05, top=0.95, left=0.04, right=0.95, wspace=0.02)

plt.savefig(image_loc + "ASL_Trends{}".format(settings.OUTFMT))

plt.close()
#************************************************************************
#                                 END
#************************************************************************
