#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for lake temperature section.
#       For BAMS SotC 2017
#
#************************************************************************
#                    SVN Info
# $Rev:: 30                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2021-06-15 10:41:02 +0100 (Tue, 15 Jun #$:  Date of last commit
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

    indata = np.genfromtxt(filename, delimiter=',', dtype=(float))

    years = indata[:, 0].astype(float)

    anomalies = utils.Timeseries("Anomalies", years, indata[:, 1])
    uncertainties = utils.Timeseries("Uncertainties", years, indata[:, 2])

    return anomalies, uncertainties # read_ts

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
                        transform=cartopy.crs.Geodetic(), edgecolor='0.5', linewidth='0.5')

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

timeseries, uncertainties = read_ts(data_loc + "fig1a.csv")

fig = plt.figure(figsize=(8, 9))
plt.clf()

# make axes by hand
ax1 = plt.axes([0.1, 0.8, 0.75, 0.15])
ax2 = plt.axes([0.1, 0.25, 0.75, 0.5], sharex=ax1)
cbar_ax2 = plt.axes([0.88, 0.05, 0.03, 0.7])
ax3 = plt.axes([0.1, 0.05, 0.75, 0.15], sharex=ax1)

#***************
# the timeseries

ax1.plot(timeseries.times, timeseries.data, c="0.5", ls="-", lw=2)

ax1.fill_between(timeseries.times, timeseries.data+uncertainties.data, timeseries.data-uncertainties.data, color="0.75")

slope, upper, lower = utils.median_pairwise_slopes(timeseries.times, timeseries.data, -99.9, sigma=1.0)

plotx, ploty = utils.mpw_plot_points(slope, timeseries.times, timeseries.data)

ax1.plot(plotx, ploty, c="k", ls="-", lw=1)

# prettify
ax1.axhline(0, c='0.5', ls='--')
utils.thicken_panel_border(ax1)
ax1.set_ylim([-1, 1])
ax1.set_xlim([timeseries.times[0]-1, int(settings.YEAR)+1])
ax1.set_ylabel("Anomaly ("+r"$^{\circ}$"+"C)")
ax1.text(-0.1, 0.95, "(a)", transform=ax1.transAxes)
ax1.set_title("Global average lake temperature anomalies")
ax1.yaxis.set_ticks_position('left')

minorLocator = MultipleLocator(1)

#***************
# the first hovmuller

# adjust years for plotting
years = np.copy(timeseries.times)
years = np.append(years, years[-1] + 1)
years = years-0.5

# read in data
data = read_hovmuller_style(data_loc + "fig1b.csv")

# set up color ranges
bounds = [-8, -1.0, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1.0, 8]
cmap = settings.COLOURMAP_DICT["temperature"]
norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)
hov_cmap = copy.deepcopy(cmap)
hov_cmap.set_bad("0.5", 1.)

mesh = ax2.pcolormesh(years, np.arange(data.shape[0]+1), data, cmap=hov_cmap, norm=norm)

cb = plt.colorbar(mesh, cax=cbar_ax2, orientation='vertical', ticks=bounds[1:-1], \
                    label="Anomaly ("+r"$^{\circ}$"+"C)", drawedges=True)

# horizontal lines to delineate continents
ax2.axhline(689-206, c="0.5", ls="--")
ax2.axhline(689-406, c="0.5", ls="--")
ax2.axhline(689-516, c="0.5", ls="--")
ax2.axhline(689-585, c="0.5", ls="--")
ax2.axhline(689-678, c="0.5", ls="--")


# prettify
cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
cb.outline.set_linewidth(2)
cb.dividers.set_color('k')
cb.dividers.set_linewidth(2)
plt.xlim([1994, int(settings.YEAR)+1])
ax2.set_ylim([0, data.shape[0] + 1])
ax2.yaxis.set_ticks_position('none')

ax2.text(-0.1, 0.95, "(b)", transform=ax2.transAxes)
ax2.set_yticklabels("")
ax2.set_ylabel("Lakes ordered from North (top) to \n South (bottom), per continent")
ax2.set_title("Satellite-derived lake temperature anomalies")

utils.thicken_panel_border(ax2)

#***************
# the second hovmuller

# read in data
data = read_hovmuller_style(data_loc + "fig1c.csv")

mesh = ax3.pcolormesh(years, np.arange(data.shape[0]+1), data, cmap=hov_cmap, norm=norm)

ax3.text(-0.1, 0.95, "(c)", transform=ax3.transAxes)
ax3.set_yticklabels("")
ax3.set_title("In-situ lake temperature anomalies")
ax3.set_ylim([0, data.shape[0] + 1])
ax3.yaxis.set_ticks_position('none')

# horizontal lines to delineate continents
ax3.axhline(35-14, c="0.5", ls="--")
ax3.axhline(35-16, c="0.5", ls="--")
ax3.axhline(35-33, c="0.5", ls="--")

utils.thicken_panel_border(ax3)
ax3.set_ylabel("As (b)")
#
plt.savefig(image_loc + "LKT_ts_hovmuller{}".format(settings.OUTFMT))
plt.close()
plt.clf()

#***************
# Figure 2

anomalies = read_lakes(data_loc + "fig2_lake.csv")

#***************
# Anomaly Scatter map

bounds = [-8, -2, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2, 8]

lons = np.arange(-90, 120, 30)
lats = np.arange(-180, 210, 30)
dummy = np.ma.zeros((len(lats), len(lons)))
dummy.mask = np.ones(dummy.shape)

cube = utils.make_iris_cube_2d(dummy, lats, lons, "blank", "m")

utils.plot_smooth_map_iris(image_loc + "LKT_anomaly", cube, settings.COLOURMAP_DICT["temperature"], \
                               bounds, "Anomaly ("+r"$^{\circ}$"+"C)", \
                               scatter=[anomalies[1], anomalies[0], anomalies[2]], figtext="", title="")

utils.plot_smooth_map_iris(image_loc + "p2.1_LKT_anomaly", cube, settings.COLOURMAP_DICT["temperature"], \
                               bounds, "Anomaly ("+r"$^{\circ}$"+"C)", \
                               scatter=[anomalies[1], anomalies[0], anomalies[2]], \
                               figtext="(b) Lake Temperatures", title="")

#***************
# Insets Scatter map

fig = plt.figure(figsize=(6, 8.5))
plt.clf()


bounds = [-8, -2, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2, 8]
cmap = settings.COLOURMAP_DICT["temperature"]
norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)
this_cmap = copy.copy(cmap)

cube = iris.load(data_loc + "fig2_giss.nc")[0]

# make axes by hand
axes = ([0.1, 0.55, 0.8, 0.41], [0.1, 0.13, 0.8, 0.41], [0.1, 0.07, 0.8, 0.03])

# USA
ax = plt.axes(axes[0], projection=cartopy.crs.LambertConformal())

ax.gridlines() #draw_labels=True)
ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
ax.coastlines(resolution="50m")
ax.set_extent([-140, -50, 20, 70], cartopy.crs.PlateCarree())

mesh = iris.plot.pcolormesh(cube, cmap=this_cmap, norm=norm, axes=ax)
plt.scatter(anomalies[1], anomalies[0], c=anomalies[2], cmap=this_cmap, norm=norm, s=25, \
                        transform=cartopy.crs.Geodetic(), edgecolor='0.2', linewidth='0.5')

ax.text(-0.05, 1.05, "(a)", fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)
utils.thicken_panel_border(ax)

# Europe
ax = plt.axes(axes[1], projection=cartopy.crs.LambertConformal(central_longitude=10.0))

ax.gridlines() #draw_labels=True)
ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
ax.coastlines(resolution="50m")
ax.set_extent([-20, 70, 30, 70], cartopy.crs.PlateCarree())

mesh = iris.plot.pcolormesh(cube, cmap=this_cmap, norm=norm, axes=ax)
plt.scatter(anomalies[1], anomalies[0], c=anomalies[2], cmap=this_cmap, norm=norm, s=25, \
                        transform=cartopy.crs.Geodetic(), edgecolor='0.2', linewidth='0.5')

ax.text(-0.05, 1.05, "(b)", fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)
utils.thicken_panel_border(ax)

# colourbar
cb = plt.colorbar(mesh, cax=plt.axes(axes[2]), orientation='horizontal', ticks=bounds[1:-1], \
                    label="Anomaly ("+r"$^{\circ}$"+"C)", drawedges=True)

# prettify
cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
cb.outline.set_linewidth(2)
cb.dividers.set_color('k')
cb.dividers.set_linewidth(2)

plt.savefig(image_loc + "LKT_USA_EU_scatter_map{}".format(settings.OUTFMT))
plt.close()


#************************************************************************
#   END
#************************************************************************
