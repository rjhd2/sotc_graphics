#!/usr/bin/env python
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


DATALOC = "{}/{}/data/LKT/".format(settings.ROOTLOC, settings.YEAR)

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


#************************************************************************
def run_all_plots():


    if True:
        #***************
        # Figure 1

        euro, africa, tibet, canada = read_ts(DATALOC + "BAMS2021Fig1_data_LSWT.csv")
#        euro_fit, africa_fit, tibet_fit, canada_fit = read_ts(DATALOC + "Fig1_lines_LSWT.csv")

        euro_fit = utils.Timeseries("Lake", [1995, 2020], [-0.5028, (2020-1994)*0.0396 - 0.5028])
        africa_fit = utils.Timeseries("Lake", [1995, 2020], [-0.0157, (2020-1994)*0.0067 -0.0157])
        tibet_fit = utils.Timeseries("Lake", [1995, 2020], [0.0823, (2020-1994)*0.0174 + 0.0823])
        canada_fit = utils.Timeseries("Lake", [1995, 2020], [-0.1909, (2020-1994)*0.0177 - 0.1909])

        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(8, 8), sharex=True)

        #***************
        # the timeseries
        LEGEND_LOC = ""

        utils.plot_ts_panel(ax1, [euro], "-", "temperature", loc=LEGEND_LOC)
        ax1.plot(euro_fit.times, euro_fit.data, c=settings.COLOURS["temperature"][euro_fit.name], lw=2, ls="--")
        ax1.text(1995, 0.8, "(a) Europe, 127 lakes", fontsize=settings.FONTSIZE)
        utils.plot_ts_panel(ax2, [africa], "-", "temperature", loc=LEGEND_LOC)
        ax2.plot(africa_fit.times, africa_fit.data, c=settings.COLOURS["temperature"][africa_fit.name], lw=2, ls="--")
        ax2.text(1995, 0.8, "(b) Africa, 70 lakes", fontsize=settings.FONTSIZE)
        utils.plot_ts_panel(ax3, [tibet], "-", "temperature", loc=LEGEND_LOC)
        ax3.plot(tibet_fit.times, tibet_fit.data, c=settings.COLOURS["temperature"][tibet_fit.name], lw=2, ls="--")
        ax3.text(1995, 0.8, "(c) Tibetan Plateau, 104 lakes", fontsize=settings.FONTSIZE)
        utils.plot_ts_panel(ax4, [canada], "-", "temperature", loc=LEGEND_LOC)
        ax4.plot(canada_fit.times, canada_fit.data, c=settings.COLOURS["temperature"][canada_fit.name], lw=2, ls="--")
        ax4.text(1995, 0.8, "(d) Canada, 246 lakes", fontsize=settings.FONTSIZE)


        # prettify
        for ax in [ax1, ax2, ax3, ax4]:
            ax.axhline(0, c='0.5', ls='--')
            utils.thicken_panel_border(ax)
            ax.set_ylim([-1.29, 1.29])
            ax.set_xlim([euro.times[0]-2, int(settings.YEAR)+2])
            ax.yaxis.set_ticks_position('left')
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)
        for tick in ax4.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)

        fig.text(0.01, 0.65, "Anomaly from 1996-2016 ("+r'$^\circ$'+"C)", fontsize=settings.FONTSIZE, rotation="vertical")
        fig.subplots_adjust(bottom=0.03, right=0.96, top=0.99, hspace=0.001)

        plt.savefig(settings.IMAGELOC+"LKT_ts{}".format(settings.OUTFMT))

        plt.close()


    #***************
    # Anomaly Scatter map
    if True:
        anomalies = read_lakes(DATALOC + "PlateX_data_LSWT.csv")

        bounds = [-8, -2, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2, 8]
        bounds = [-8, -1.5, -1.0, -0.5, -0.25, 0, 0.25, 0.5, 1.0, 1.5, 8]
#        bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]

        lons = np.arange(-90, 120, 30)
        lats = np.arange(-180, 210, 30)
        dummy = np.ma.zeros((len(lats), len(lons)))
        dummy.mask = np.ones(dummy.shape)

        cube = utils.make_iris_cube_2d(dummy, lats, lons, "blank", "m")

        utils.plot_smooth_map_iris(settings.IMAGELOC + "LKT_anomaly", cube, settings.COLOURMAP_DICT["temperature"], \
                                       bounds, "Anomalies from 1996-2016 ("+r"$^{\circ}$"+"C)", \
                                       scatter=[anomalies[1], anomalies[0], anomalies[2]], figtext="", title="")

        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_LKT_anomaly", cube, settings.COLOURMAP_DICT["temperature"], \
                                       bounds, "Anomalies from 1996-2016 ("+r"$^{\circ}$"+"C)", \
                                       scatter=[anomalies[1], anomalies[0], anomalies[2]], \
                                       figtext="(b) Lake Temperatures", title="")


    #***************
    # Insets Scatter map
    if False:

        fig = plt.figure(figsize=(8, 7))
        plt.clf()

        anomalies = read_lakes(DATALOC + "Fig2_data_LSWT.csv")

        bounds = [-8, -2, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2, 8]
#        bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
        cmap = settings.COLOURMAP_DICT["temperature"]
        norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)
        this_cmap = copy.copy(cmap)

        # first_cube = iris.load(DATALOC + "amaps_1st_quarter_2018_250km.nc")[0]
        # third_cube = iris.load(DATALOC + "amaps_3rd_quarter_2018_250km.nc")[0]
        # annual_cube = iris.load(DATALOC + "amaps_annual_2018_250km.nc")[0]

        # cube = iris.load(DATALOC + "lswt_anom_1979_2019.nc")[0]
        # if settings.OUTFMT in [".eps", ".pdf"]:
        #     if cube.coord("latitude").points.shape[0] > 180 or cube.coord("longitude").points.shape[0] > 360:
        #         regrid_size = 1.0
        #         print("Regridding cube for {} output to {} degree resolution".format(settings.OUTFMT, regrid_size))
        #         print("Old Shape {}".format(cube.data.shape))
        #         plot_cube = utils.regrid_cube(cube, regrid_size, regrid_size)
        #         print("New Shape {}".format(plot_cube.data.shape))
        #     else:
        #         plot_cube = copy.deepcopy(cube)
        # else:
        #     plot_cube = copy.deepcopy(cube)
        

        # make axes by hand
        axes = ([0.01, 0.55, 0.59, 0.41], [0.565, 0.45, 0.47, 0.50], [0.01, 0.13, 0.59, 0.41], [0.61, 0.07, 0.38, 0.41],[0.1, 0.1, 0.8, 0.03])

        # Europe
        ax = plt.axes(axes[0], projection=cartopy.crs.PlateCarree())

        ax.gridlines() #draw_labels=True)
        ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
        ax.coastlines(resolution="50m")
        #ax.add_feature(cartopy.feature.BORDERS.with_scale('110m'), linewidth=.5)
        ax.set_extent([-25, 40, 34, 72], cartopy.crs.PlateCarree())

        # mesh = iris.plot.pcolormesh(plot_cube, cmap=this_cmap, norm=norm, axes=ax)
        plt.scatter(anomalies[1], anomalies[0], c=anomalies[2], cmap=this_cmap, norm=norm, s=25, \
                                transform=cartopy.crs.Geodetic(), edgecolor='0.1', linewidth=0.5, zorder=10)

        ax.text(0.05, 1.05, "(a) Europe", fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)
        utils.thicken_panel_border(ax)

        # Africa
        ax = plt.axes(axes[1], projection=cartopy.crs.PlateCarree())

        ax.gridlines() #draw_labels=True)
        ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
        ax.coastlines(resolution="50m")
        #ax.add_feature(cartopy.feature.BORDERS.with_scale('110m'), linewidth=.5)
        ax.set_extent([-19, 43, -40, 33], cartopy.crs.PlateCarree())

        # lat_constraint = utils.latConstraint([25, 90]) 
        # nh_cube = third_cube.extract(lat_constraint)
        # lat_constraint = utils.latConstraint([-90, -25]) 
        # sh_cube = first_cube.extract(lat_constraint)
        # lat_constraint = utils.latConstraint([-25, 25]) 
        # trop_cube = annual_cube.extract(lat_constraint)

        # mesh = iris.plot.pcolormesh(nh_cube, cmap=this_cmap, norm=norm, axes=ax)
        # mesh = iris.plot.pcolormesh(trop_cube, cmap=this_cmap, norm=norm, axes=ax)
        # mesh = iris.plot.pcolormesh(sh_cube, cmap=this_cmap, norm=norm, axes=ax)
        # mesh = iris.plot.pcolormesh(plot_cube, cmap=this_cmap, norm=norm, axes=ax)

        plt.scatter(anomalies[1], anomalies[0], c=anomalies[2], cmap=this_cmap, norm=norm, s=25, \
                                transform=cartopy.crs.Geodetic(), edgecolor='0.1', linewidth=0.5, zorder=10)

        ax.text(0.05, 1.05, "(b) Africa", fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)
        utils.thicken_panel_border(ax)

        # Canada
        ax = plt.axes(axes[2], projection=cartopy.crs.PlateCarree())

        ax.gridlines() #draw_labels=True)
        ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
        ax.coastlines(resolution="50m")
        #ax.add_feature(cartopy.feature.BORDERS.with_scale('110m'), linewidth=.5)
        ax.set_extent([-140, -55, 42, 82], cartopy.crs.PlateCarree())

        # mesh = iris.plot.pcolormesh(plot_cube, cmap=this_cmap, norm=norm, axes=ax)
        plt.scatter(anomalies[1], anomalies[0], c=anomalies[2], cmap=this_cmap, norm=norm, s=25, \
                                transform=cartopy.crs.Geodetic(), edgecolor='0.1', linewidth=0.5, zorder=10)

        ax.text(0.05, 1.05, "(c) Canada", fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)
        utils.thicken_panel_border(ax)

        # Tibet
        ax = plt.axes(axes[3], projection=cartopy.crs.PlateCarree())

        ax.gridlines() #draw_labels=True)
        ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
        ax.coastlines(resolution="50m")
        ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linewidth=.5)
        ax.set_extent([78, 102, 28, 39], cartopy.crs.PlateCarree())

        # mesh = iris.plot.pcolormesh(plot_cube, cmap=this_cmap, norm=norm, axes=ax)
        plt.scatter(anomalies[1], anomalies[0], c=anomalies[2], cmap=this_cmap, norm=norm, s=25, \
                                transform=cartopy.crs.Geodetic(), edgecolor='0.1', linewidth=0.5, zorder=10)

        ax.text(0.05, 1.05, "(d) Tibetan Plateau", fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)
        utils.thicken_panel_border(ax)


        # colourbar
        # cb = plt.colorbar(mesh, cax=plt.axes(axes[4]), orientation='horizontal', ticks=bounds[1:-1], drawedges=True)
        cb = plt.colorbar(cax=plt.axes(axes[4]), orientation='horizontal', ticks=bounds[1:-1], drawedges=True)

        # prettify
        cb.ax.tick_params(axis='x', labelsize=settings.FONTSIZE, direction='in', size=0)
        cb.set_label(label="Anomalies from 1996-2016 ("+r"$^{\circ}$"+"C)", fontsize=settings.FONTSIZE)
        cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
        cb.outline.set_linewidth(2)
        cb.dividers.set_color('k')
        cb.dividers.set_linewidth(2)

        plt.savefig(settings.IMAGELOC + "LKT_Regions_scatter_map{}".format(settings.OUTFMT))
        plt.close()

    #***************
    # Insets Scatter map - individual panels.
    if True:

        fig = plt.figure(figsize=(8, 12))
        plt.clf()

        anomalies = read_lakes(DATALOC + "Fig2_data_LSWT.csv")

        bounds = [-8, -2, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2, 8]
#        bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
        cmap = settings.COLOURMAP_DICT["temperature"]
        norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)
        this_cmap = copy.copy(cmap)    

        # make axes by hand
        axes = ([0.05, 0.64, 0.9, 0.36], [0.05, 0.355, 0.9, 0.29], [0.05, 0.075, 0.9, 0.28], [0.05, 0.045, 0.9, 0.03])

        # Europe
        ax = plt.axes(axes[0], projection=cartopy.crs.PlateCarree())

        ax.gridlines() #draw_labels=True)
        ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
        ax.coastlines(resolution="50m")
        ax.set_extent([-25, 40, 34, 72], cartopy.crs.PlateCarree())

        plt.scatter(anomalies[1], anomalies[0], c=anomalies[2], cmap=this_cmap, norm=norm, s=70, \
                                transform=cartopy.crs.Geodetic(), edgecolor='0.1', linewidth=0.5, zorder=10)

        ax.text(0.05, 0.9, "(e) Europe", fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)
        utils.thicken_panel_border(ax)

        # # Africa
        # ax = plt.axes(axes[1], projection=cartopy.crs.PlateCarree())

        # ax.gridlines() #draw_labels=True)
        # ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
        # ax.coastlines(resolution="50m")
        # ax.set_extent([-19, 43, -40, 33], cartopy.crs.PlateCarree())

        # plt.scatter(anomalies[1], anomalies[0], c=anomalies[2], cmap=this_cmap, norm=norm, s=25, \
        #                         transform=cartopy.crs.Geodetic(), edgecolor='0.1', linewidth=0.5, zorder=10)

        # ax.text(0.05, 1.05, "(b) Africa", fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)
        # utils.thicken_panel_border(ax)

        # Canada
        ax = plt.axes(axes[1], projection=cartopy.crs.PlateCarree())

        ax.gridlines() #draw_labels=True)
        ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
        ax.coastlines(resolution="50m")
        ax.set_extent([-140, -55, 42, 82], cartopy.crs.PlateCarree())

        plt.scatter(anomalies[1], anomalies[0], c=anomalies[2], cmap=this_cmap, norm=norm, s=70, \
                                transform=cartopy.crs.Geodetic(), edgecolor='0.1', linewidth=0.5, zorder=10)

        ax.text(0.05, 0.9, "(f) Canada", fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)
        utils.thicken_panel_border(ax)

        # Tibet
        ax = plt.axes(axes[2], projection=cartopy.crs.PlateCarree())

        ax.gridlines() #draw_labels=True)
        ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
        ax.coastlines(resolution="50m")
        ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linewidth=.5)
        ax.set_extent([78, 102, 28, 39], cartopy.crs.PlateCarree())

        plt.scatter(anomalies[1], anomalies[0], c=anomalies[2], cmap=this_cmap, norm=norm, s=70, \
                                transform=cartopy.crs.Geodetic(), edgecolor='0.1', linewidth=0.5, zorder=10)

        ax.text(0.05, 0.9, "(g) Tibetan Plateau", fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)
        utils.thicken_panel_border(ax)


        # colourbar
        # cb = plt.colorbar(mesh, cax=plt.axes(axes[4]), orientation='horizontal', ticks=bounds[1:-1], drawedges=True)
        cb = plt.colorbar(cax=plt.axes(axes[3]), orientation='horizontal', ticks=bounds[1:-1], drawedges=True)

        # prettify
        cb.ax.tick_params(axis='x', labelsize=settings.FONTSIZE, direction='in', size=0)
        cb.set_label(label="Anomalies from 1996-2016 ("+r"$^{\circ}$"+"C)", fontsize=settings.FONTSIZE)
        cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
        cb.outline.set_linewidth(2)
        cb.dividers.set_color('k')
        cb.dividers.set_linewidth(2)

        plt.savefig(settings.IMAGELOC + "LKT_Regions_scatter_map{}".format(settings.OUTFMT))
        plt.close()


#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
