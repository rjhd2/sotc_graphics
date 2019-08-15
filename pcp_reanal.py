#!/bin/python
#************************************************************************
#
#  Plot figures and output numbers for Sidebar on Reanalyses and Precip.
#       For BAMS SotC 2018
#
#************************************************************************
#                    SVN Info
# $Rev:: 23                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2018-06-05 17:55:11 +0100 (Tue, 05 Jun #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import calendar
import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl

import iris
import cartopy
import cartopy.feature as cfeature


import utils # RJHD utilities
import settings


data_loc = "{}/{}/data/SIDE/".format(settings.ROOTLOC, settings.YEAR)
image_loc = "{}/{}/images/".format(settings.ROOTLOC, settings.YEAR)

LEGEND_LOC = 'lower right'
LW = 3

land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                            edgecolor='face',
                                            facecolor=cfeature.COLORS['land'])

# get data 
cube_list = iris.load(data_loc + "g4.accumulate.GPM_3IMERGHHL_05_precipitationCal.20180721-20180726.80W_36N_71W_44N.nc")
cube = cube_list[0]

# set up colour maps
bounds = [0, 10, 20, 40, 60, 80, 100, 150, 200, 250, 300]
cmap = settings.COLOURMAP_DICT["precip_sequential"]
norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)

# set up figure
fig = plt.figure(figsize=(8, 7))
ax = plt.axes([0.1, 0.1, 0.8, 0.8], projection=cartopy.crs.PlateCarree())
ax.set_extent([-79, -74, 38, 42], cartopy.crs.PlateCarree())
               
# extra fratures 
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')
ax.add_feature(states_provinces, edgecolor='gray', lw=1)
ax.coastlines(resolution="10m", linewidth=1, edgecolor="k")
ax.add_feature(land_10m, zorder=0, facecolor="0.9", edgecolor="k")

# add other features
gl = ax.gridlines(draw_labels=True)
gl.xlabel_style = {'size': settings.FONTSIZE*0.7}
gl.ylabel_style = {'size': settings.FONTSIZE*0.7}
ax.add_feature(cartopy.feature.BORDERS, zorder=0, facecolor="0.9", edgecolor="k")
ext = ax.get_extent()

# plot data
mesh = iris.plot.pcolormesh(cube, cmap=cmap, norm=norm)

# metadata
plt.scatter([-76.66], [39.16], c="k", s=75, transform=cartopy.crs.Geodetic(), edgecolor='0.5', linewidth=0.5)
plt.text(-76.7, 39.2, "BWI", transform=cartopy.crs.Geodetic(), fontsize=settings.FONTSIZE, color="k", ha="right")

locs, = np.argwhere(cube.data == np.max(cube.data))
lat = cube.coord("latitude").points[locs[0]]
lon = cube.coord("longitude").points[locs[1]]
plt.text(lon, lat, "x", transform=cartopy.crs.Geodetic(), fontsize=settings.FONTSIZE, color="0.9", ha="center", va="center")
plt.text(lon+0.1, lat, "Max = {:5.1f}mm".format(np.max(cube.data)), transform=cartopy.crs.Geodetic(), fontsize=settings.FONTSIZE*0.8, color="0.1", ha="left")


# and colourbar
cb = plt.colorbar(mesh, orientation='horizontal', pad=0.06, fraction=0.05, \
                        aspect=30, ticks=bounds[1:-1], drawedges=True)

# prettify colourbar
cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
cb.ax.tick_params(axis='x', labelsize=settings.FONTSIZE*0.6, direction='in')
cb.set_label(label="Accumulation 21-26 July 2018 (mm)", fontsize=settings.FONTSIZE*0.6)
cb.outline.set_linewidth(2)
cb.dividers.set_color('k')
cb.dividers.set_linewidth(2)

# ensure correct extent
ax.set_extent(ext, ax.projection) # fix the extent change from colormesh

# save
plt.savefig(image_loc + "SB2.1_{}_prcp_rnl_USA".format(settings.YEAR) + settings.OUTFMT)
plt.close()


