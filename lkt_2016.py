#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for lake temperature section.
#       For BAMS SotC 2015 
#
#************************************************************************
#                    SVN Info
# $Rev:: 26                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2019-04-17 15:34:18 +0100 (Wed, 17 Apr #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as mpl_cm
import matplotlib as mpl
import cartopy
from matplotlib.ticker import MultipleLocator

import utils # RJHD utilities
import settings

import copy

data_loc = "{}/{}/data/LKT/".format(settings.ROOTLOC, settings.YEAR)
image_loc = "{}/{}/images/".format(settings.ROOTLOC, settings.YEAR)

CLIM_PERIOD="8110"

FONTSIZE = 20

LAKE_COLOURS = {"Global" : "k", "Northern Hemisphere" : "r", "Southern Hemisphere" : "b", "Tropics" : "g"}

#************************************************************************
def read_ts(filename):
    '''
    Read timeseries
    '''

    MDI = -99.9

    indata = np.genfromtxt(filename, delimiter = ',', skip_header = 1, dtype = (str))
   
    years = indata[:,0].astype(float)

    locs = np.ma.where(indata == "NA")
    indata[locs] = MDI
    indata = indata.astype(float)
    indata = np.ma.masked_where(indata == MDI, indata)

    glob = utils.Timeseries("Global", years, indata[:,1])
    north = utils.Timeseries("Northern Hemisphere", years, indata[:,2])
    south = utils.Timeseries("Southern Hemisphere", years, indata[:,3])
    tropics = utils.Timeseries("Tropics", years, indata[:,4])
    
    return glob, north, south, tropics # read_ts


#************************************************************************
def read_lakes(filename):

    MDI = -99.9

    indata = np.genfromtxt(filename, delimiter = ',', skip_header = 1, dtype = (str))
    
    locs = np.where(indata == "NA")
    indata[locs] = MDI
    indata = indata.astype(float)

    lats = indata[:,1]
    lons = indata[:,2]
    
    data = indata[:,3]

    data = np.ma.masked_where(data <= MDI, data)
    
    lons = np.ma.array(lons, mask = data.mask)
    lats = np.ma.array(lats, mask = data.mask)


    return lats, lons, data # read_lakes


#************************************************************************
def plot_lakes(ax, lons, lats, values, cmap, norm, cb_label, bounds, figtext):

    ax.gridlines() #draw_labels=True)
    ax.add_feature(cartopy.feature.LAND, zorder = 0, facecolor = "0.9", edgecolor = "k")
    ax.coastlines()
    
    ext = ax.get_extent() # save the original extent
    
    scatter = plt.scatter(lons, lats, c = values, cmap = cmap, norm = norm, s=25, \
                        transform = cartopy.crs.Geodetic(), edgecolor = '0.5', linewidth='0.5')

    # colorbar
    cb=plt.colorbar(scatter,  orientation = 'horizontal', ticks = bounds[1:-1], label = cb_label, drawedges=True, pad = 0.05, fraction = 0.05, aspect = 30)
    
    # prettify
    cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
    cb.outline.set_color('k')
    cb.outline.set_linewidth(2)
    cb.dividers.set_color('k')
    cb.dividers.set_linewidth(2)
    
    ax.set_extent(ext, ax.projection) 
    ax.text(0.03, 0.95, figtext, transform = ax.transAxes)

    return # plot_lakes

#************************************************************************
def read_hovmuller_style(filename):

    MDI = -99.9

    indata = np.genfromtxt(filename, delimiter = ',', dtype = (str))
    
    years = indata[0,1:].astype(int)
    locations = indata[1:,0]

    locs = np.ma.where(indata == "NA")
    indata[locs] = MDI

    data = indata[1:,1:].astype(float)
    data = np.ma.masked_where(data == MDI, data)
#    data = np.ma.masked_where(np.abs(data) < 0.25, data)
    data.fill_value = MDI

    return locations, years, data # read_hovmuller_style

#***************
# Figure 1

locations, years, data = read_hovmuller_style(data_loc + "fig1.csv")

# adjust years for plotting
years = np.append(years, years[-1] + 1)
years = years-0.5


# set up color ranges
bounds = [-8, -1.5, -1.0, -0.5, -0.25, 0.25, 0.5, 1.0, 1.5, 8]
cmap = settings.COLOURMAP_DICT["temperature"]
norm=mpl.cm.colors.BoundaryNorm(bounds,cmap.N)

fig = plt.figure(figsize = (8,11))
plt.clf()

# make axes by hand
ax = plt.axes([0.15,0.05,0.8,0.88])

this_cmap = copy.copy(cmap)
this_cmap.set_bad("0.5", 1.)
mesh = ax.pcolormesh(years, np.arange(len(locations)+1), data, cmap = this_cmap, norm = norm)

cb=plt.colorbar(mesh,  orientation = 'horizontal', ticks = bounds[1:-1], label = "Lake Temperature Anomaly\n("+r"$^{\circ}$"+"C)", drawedges=True, pad = 0.05, fraction = 0.05, aspect = 20)
    
# prettify
cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
cb.outline.set_color('k')
cb.outline.set_linewidth(2)
cb.dividers.set_color('k')
cb.dividers.set_linewidth(2)

#
plt.xlim([1979.5,2016.5])
ax.axhline(12, c = 'k', ls = '-', lw = 2, zorder = 10)
ax.axhline(14, c = 'k', ls = '-', lw = 2, zorder = 10)
ax.axhline(31, c = 'k', ls = '-', lw = 2, zorder = 10)
ax.axhline(35, c = 'k', ls = '-', lw = 2, zorder = 10)
ax.axhline(39, c = 'k', ls = '-', lw = 2, zorder = 10)
ax.axhline(48, c = 'k', ls = '-', lw = 2, zorder = 10)
plt.ylim([0,48.1])

ax.set_yticks([6,13,22,33,37,42])
ax.set_yticklabels(["NZ", "Asia", "Europe", "Great Lakes", "US", "Canada"])

minorLocator = MultipleLocator(1)
ax.xaxis.set_minor_locator(minorLocator)

utils.thicken_panel_border(ax)

plt.savefig(image_loc + "LKT_hovmuller{}".format(settings.OUTFMT))
plt.close()

#***************
# Figure 2

timeseries = read_ts(data_loc + "fig2a.csv")

trends = read_lakes(data_loc + "fig2c.csv")

anomalies = read_lakes(data_loc + "fig2b.csv")

# borrow from slp.py
fig = plt.figure(figsize = (8,9))
plt.clf()

# make axes by hand
axes = ([0.15,0.62,0.8,0.34],[0.05,0.07,0.9,0.5])#,[0.05,0.05,0.8,0.3])
 
# the timeseries
ax = plt.axes(axes[0])
for ts in timeseries:
    ax.plot(ts.times, ts.data, c = LAKE_COLOURS[ts.name], ls = "-", label = ts.name, lw = 2)

ax.legend(loc = "lower center", ncol = 2, frameon = False, prop = {'size':settings.LEGEND_FONTSIZE*0.8}, labelspacing = 0.1, columnspacing = 0.5)

# prettify
ax.axhline(0, c = '0.5', ls = '--')
utils.thicken_panel_border(ax)
ax.set_ylim([-1,1])
ax.set_xlim([None,2017])
ax.set_ylabel("Lake Temperature Anomaly\n("+r"$^{\circ}$"+"C)")
ax.text(-0.1, 0.95, "(a)", transform = ax.transAxes)

minorLocator = MultipleLocator(1)
ax.xaxis.set_minor_locator(minorLocator)

# set up color ranges
bounds = [-8, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 8]
cmap = settings.COLOURMAP_DICT["temperature"]
norm=mpl.cm.colors.BoundaryNorm(bounds,cmap.N)

# the trends
ax = plt.axes(axes[1], projection=cartopy.crs.Robinson())     
plot_lakes(ax, trends[1], trends[0], trends[2]*10, cmap, norm, "Lake Temperature Trend ("+r"$^{\circ}$"+"C decade"+r"$^{-1}$"+")", bounds, "(b)")

# the anomalies --> moved to plate 2.1 - see below
# set up color ranges
#bounds = [-8, -2, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2, 8]
#norm=mpl.cm.colors.BoundaryNorm(bounds,cmap.N)
#ax = plt.axes(axes[2], projection=cartopy.crs.Robinson())     
#plot_lakes(ax, anomalies[1], anomalies[0], anomalies[2], cmap, norm, "Lake Temperature Anomaly\n("+r"$^{\circ}$"+"C)", bounds, "(c)")

plt.savefig(image_loc + "LKT_ts_maps{}".format(settings.OUTFMT))
plt.close()



#***************
# Anomaly Scatter map

bounds = [-8, -2, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2, 8]

lons = np.arange(-90,120,30)
lats = np.arange(-180,210,30)
dummy = np.ma.zeros((len(lats), len(lons)))
dummy.mask = np.ones(dummy.shape)

cube = utils.make_iris_cube_2d(dummy, lats, lons, "blank", "m")


utils.plot_smooth_map_iris(image_loc + "LKT_anomaly", cube, settings.COLOURMAP_DICT["temperature"], bounds, "Anomaly ("+r"$^{\circ}$"+"C)", scatter = [anomalies[1], anomalies[0], anomalies[2]], figtext = "", title = "")
utils.plot_smooth_map_iris(image_loc + "p2.1_LKT_anomaly", cube, settings.COLOURMAP_DICT["temperature"], bounds, "Anomaly ("+r"$^{\circ}$"+"C)", scatter = [anomalies[1], anomalies[0], anomalies[2]], figtext = "(b) Lake Temperatures", title = "")


#************************************************************************
#   END
#************************************************************************
