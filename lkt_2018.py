#!/usr/local/sci/python
# python3
from __future__ import absolute_import
from __future__ import print_function
#************************************************************************
#
#  Plot figures and output numbers for lake temperature section.
#       For BAMS SotC 2017
#
#************************************************************************
#                    SVN Info
# $Rev:: 24                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2019-03-11 12:46:29 +0000 (Mon, 11 Mar #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************

import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator

import iris
import cartopy

import utils # RJHD utilities
import settings


data_loc = "{}/{}/data/LKT/".format(settings.ROOTLOC, settings.YEAR)
image_loc = "{}/{}/images/".format(settings.ROOTLOC, settings.YEAR)

CLIM_PERIOD = "8110"

FONTSIZE = 20

LAKE_COLOURS = {"Global" : "k", "Northern Hemisphere" : "r", "Southern Hemisphere" : "b", "Tropics" : "g"}

#************************************************************************
def read_ts(filename):
    '''
    Read timeseries
    '''

    indata = np.genfromtxt(filename, delimiter=',', dtype=(float), skip_header=1)

    years = indata[:, 0].astype(float)

    euro = utils.Timeseries("Lake", years, indata[:, 1])
    africa = utils.Timeseries("Lake", years, indata[:, 2])
    tibet = utils.Timeseries("Lake", years, indata[:, 3])
    canada = utils.Timeseries("Lake", years, indata[:, 4])

    return euro, africa, tibet, canada # read_ts

#************************************************************************
def read_lakes(filename):

    MDI = -99.9

    indata = np.genfromtxt(filename, delimiter=',', skip_header=1, dtype=(str))

    locs = np.where(indata == "NA")
    indata[locs] = MDI
    indata = indata.astype(float)

    lats = indata[:, 0]
    lons = indata[:, 1]

    data = indata[:, 2]

    data = np.ma.masked_where(data <= MDI, data)

    lons = np.ma.array(lons, mask=data.mask)
    lats = np.ma.array(lats, mask=data.mask)

    return lats, lons, data # read_lakes

#************************************************************************
def plot_lakes(ax, lons, lats, values, cmap, norm, cb_label, bounds, figtext):

    ax.gridlines() #draw_labels=True)
    ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
    ax.coastlines()

    ext = ax.get_extent() # save the original extent

    scatter = plt.scatter(lons, lats, c=values, cmap=cmap, norm=norm, s=25, \
                        transform=cartopy.crs.Geodetic(), edgecolor='0.5', linewidth=0.5)

    # colorbar
    cb = plt.colorbar(scatter, orientation='horizontal', ticks=bounds[1:-1], label=cb_label, \
                        drawedges=True, pad=0.05, fraction=0.05, aspect=30)

    # prettify
    cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
    cb.outline.set_color('k')
    cb.outline.set_linewidth(2)
    cb.dividers.set_color('k')
    cb.dividers.set_linewidth(2)

    ax.set_extent(ext, ax.projection)
    ax.text(0.03, 0.95, figtext, transform=ax.transAxes)

    return # plot_lakes

#************************************************************************
def read_hovmuller_style(filename):

    MDI = -99.9

    indata = np.genfromtxt(filename, delimiter=',', dtype=(str), skip_header=1)

    locs = np.ma.where(indata == "NA")
    indata[locs] = MDI
    locs = np.ma.where(indata == "NaN")
    indata[locs] = MDI

    data = indata[:, 2:].astype(float)
    data = np.ma.masked_where(data == MDI, data)
#    data = np.ma.masked_where(np.abs(data) < 0.25, data)
    data.fill_value = MDI

    data = np.flipud(data)

    return data # read_hovmuller_style




#***************
# Figure 1

euro, africa, tibet, canada = read_ts(data_loc + "Fig1_data_LSWT.csv")
euro_fit, africa_fit, tibet_fit, canada_fit = read_ts(data_loc + "Fig1_lines_LSWT.csv")

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(10, 13), sharex=True)

#***************
# the timeseries
LEGEND_LOC = ""

utils.plot_ts_panel(ax1, [euro], "-", "temperature", loc=LEGEND_LOC)
ax1.plot(euro_fit.times, euro_fit.data, c=settings.COLOURS["temperature"][euro_fit.name], lw=2, ls="--")
ax1.text(1995, 0.75, "Europe, 127 lakes", fontsize=settings.FONTSIZE)
utils.plot_ts_panel(ax2, [africa], "-", "temperature", loc=LEGEND_LOC)
ax2.plot(africa_fit.times, africa_fit.data, c=settings.COLOURS["temperature"][africa_fit.name], lw=2, ls="--")
ax2.text(1995, 0.75, "Africa, 68 lakes", fontsize=settings.FONTSIZE)
utils.plot_ts_panel(ax3, [tibet], "-", "temperature", loc=LEGEND_LOC)
ax3.plot(tibet_fit.times, tibet_fit.data, c=settings.COLOURS["temperature"][tibet_fit.name], lw=2, ls="--")
ax3.text(1995, 0.75, "Tibetan Plateau, 106 lakes", fontsize=settings.FONTSIZE)
utils.plot_ts_panel(ax4, [canada], "-", "temperature", loc=LEGEND_LOC)
ax4.plot(canada_fit.times, canada_fit.data, c=settings.COLOURS["temperature"][canada_fit.name], lw=2, ls="--")
ax4.text(1995, 0.75, "Canada, 245 lakes", fontsize=settings.FONTSIZE)


# prettify
for ax in [ax1, ax2, ax3, ax4]:
    ax.axhline(0, c='0.5', ls='--')
    utils.thicken_panel_border(ax)
    ax.set_ylim([-1, 1])
    ax.set_xlim([euro.times[0]-1, int(settings.YEAR)+1])
    ax.yaxis.set_ticks_position('left')
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)
for tick in ax4.xaxis.get_major_ticks():
    tick.label.set_fontsize(settings.FONTSIZE)

fig.text(0.01, 0.65, "Anomaly from 1996-2016 ("+r'$^\circ$'+"C)", fontsize=settings.FONTSIZE, rotation="vertical")
fig.subplots_adjust(right=0.95, top=0.95, hspace=0.001)

plt.savefig(image_loc+"LKT_ts{}".format(settings.OUTFMT))

plt.close()


#***************
# Anomaly Scatter map
anomalies = read_lakes(data_loc + "PlateX_data_LSWT.csv")

bounds = [-8, -2, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2, 8]
bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]

lons = np.arange(-90, 120, 30)
lats = np.arange(-180, 210, 30)
dummy = np.ma.zeros((len(lats), len(lons)))
dummy.mask = np.ones(dummy.shape)

cube = utils.make_iris_cube_2d(dummy, lats, lons, "blank", "m")

utils.plot_smooth_map_iris(image_loc + "LKT_anomaly", cube, settings.COLOURMAP_DICT["temperature"], \
                               bounds, "Anomalies from 1996-2016 ("+r"$^{\circ}$"+"C)", \
                               scatter=[anomalies[1], anomalies[0], anomalies[2]], figtext="", title="")

utils.plot_smooth_map_iris(image_loc + "p2.1_LKT_anomaly", cube, settings.COLOURMAP_DICT["temperature"], \
                               bounds, "Anomalies from 1996-2016 ("+r"$^{\circ}$"+"C)", \
                               scatter=[anomalies[1], anomalies[0], anomalies[2]], \
                               figtext="(b) Lake Temperatures", title="")


#***************
# Insets Scatter map

fig = plt.figure(figsize=(8, 7))
plt.clf()

anomalies = read_lakes(data_loc + "Fig2_data_LSWT.csv")

bounds = [-8, -2.5, -2, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2, 2.5, 8]
bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
cmap = settings.COLOURMAP_DICT["temperature"]
norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)
this_cmap = copy.copy(cmap)

first_cube = iris.load(data_loc + "amaps_1st_quarter_2018_250km.nc")[0]
third_cube = iris.load(data_loc + "amaps_3rd_quarter_2018_250km.nc")[0]
annual_cube = iris.load(data_loc + "amaps_annual_2018_250km.nc")[0]

# make axes by hand
axes = ([0.01, 0.55, 0.59, 0.41], [0.565, 0.45, 0.47, 0.50], [0.01, 0.13, 0.59, 0.41], [0.61, 0.07, 0.38, 0.41],[0.1, 0.08, 0.8, 0.03])

# Europe
ax = plt.axes(axes[0], projection=cartopy.crs.PlateCarree())

ax.gridlines() #draw_labels=True)
ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
ax.coastlines(resolution="50m")
#ax.add_feature(cartopy.feature.BORDERS.with_scale('110m'), linewidth=.5)
ax.set_extent([-25, 40, 34, 72], cartopy.crs.PlateCarree())

mesh = iris.plot.pcolormesh(third_cube, cmap=this_cmap, norm=norm, axes=ax)
plt.scatter(anomalies[1], anomalies[0], c=anomalies[2], cmap=this_cmap, norm=norm, s=25, \
                        transform=cartopy.crs.Geodetic(), edgecolor='0.2', linewidth=0.5, zorder=10)

ax.text(0.05, 1.05, "(a) Europe", fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)
utils.thicken_panel_border(ax)

# Africa
ax = plt.axes(axes[1], projection=cartopy.crs.PlateCarree())

ax.gridlines() #draw_labels=True)
ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
ax.coastlines(resolution="50m")
#ax.add_feature(cartopy.feature.BORDERS.with_scale('110m'), linewidth=.5)
ax.set_extent([-19, 43, -40, 33], cartopy.crs.PlateCarree())

lat_constraint = utils.latConstraint([25, 90]) 
nh_cube = third_cube.extract(lat_constraint)
lat_constraint = utils.latConstraint([-90, -25]) 
sh_cube = first_cube.extract(lat_constraint)
lat_constraint = utils.latConstraint([-25, 25]) 
trop_cube = annual_cube.extract(lat_constraint)

mesh = iris.plot.pcolormesh(nh_cube, cmap=this_cmap, norm=norm, axes=ax)
mesh = iris.plot.pcolormesh(trop_cube, cmap=this_cmap, norm=norm, axes=ax)
mesh = iris.plot.pcolormesh(sh_cube, cmap=this_cmap, norm=norm, axes=ax)

plt.scatter(anomalies[1], anomalies[0], c=anomalies[2], cmap=this_cmap, norm=norm, s=25, \
                        transform=cartopy.crs.Geodetic(), edgecolor='0.2', linewidth=0.5, zorder=10)

ax.text(0.05, 1.05, "(b) Africa", fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)
utils.thicken_panel_border(ax)

# Canada
ax = plt.axes(axes[2], projection=cartopy.crs.PlateCarree())

ax.gridlines() #draw_labels=True)
ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
ax.coastlines(resolution="50m")
#ax.add_feature(cartopy.feature.BORDERS.with_scale('110m'), linewidth=.5)
ax.set_extent([-140, -55, 42, 82], cartopy.crs.PlateCarree())

mesh = iris.plot.pcolormesh(third_cube, cmap=this_cmap, norm=norm, axes=ax)
plt.scatter(anomalies[1], anomalies[0], c=anomalies[2], cmap=this_cmap, norm=norm, s=25, \
                        transform=cartopy.crs.Geodetic(), edgecolor='0.2', linewidth=0.5, zorder=10)

ax.text(0.05, 1.05, "(c) Canada", fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)
utils.thicken_panel_border(ax)

# Tibet
ax = plt.axes(axes[3], projection=cartopy.crs.PlateCarree())

ax.gridlines() #draw_labels=True)
ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
ax.coastlines(resolution="50m")
ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linewidth=.5)
ax.set_extent([78, 102, 28, 39], cartopy.crs.PlateCarree())

mesh = iris.plot.pcolormesh(third_cube, cmap=this_cmap, norm=norm, axes=ax)
plt.scatter(anomalies[1], anomalies[0], c=anomalies[2], cmap=this_cmap, norm=norm, s=25, \
                        transform=cartopy.crs.Geodetic(), edgecolor='0.2', linewidth=0.5, zorder=10)

ax.text(0.05, 1.05, "(d) Tibetan Plateau", fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)
utils.thicken_panel_border(ax)


# colourbar
cb = plt.colorbar(mesh, cax=plt.axes(axes[4]), orientation='horizontal', ticks=bounds[1:-1], \
                    label="Anomalies from 1996-2016 ("+r"$^{\circ}$"+"C)", drawedges=True)

# prettify
cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
cb.outline.set_linewidth(2)
cb.dividers.set_color('k')
cb.dividers.set_linewidth(2)

plt.savefig(image_loc + "LKT_Regions_scatter_map{}".format(settings.OUTFMT))
plt.close()


#************************************************************************
#   END
#************************************************************************
