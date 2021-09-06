#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for lake water level section.
#       For BAMS SotC 2019
#
#************************************************************************
#                    SVN Info
# $Rev:: 24                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2019-03-11 12:46:29 +0000 (Mon, 11 Mar #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import os
import datetime as dt
import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
import matplotlib.dates as mdates

import iris
import cartopy

import utils # RJHD utilities
import settings


DATALOC = "{}/{}/data/LWL/".format(settings.ROOTLOC, settings.YEAR)

#************************************************************************
def read_lakes(filename):

    MDI = -99.9

    indata = np.genfromtxt(filename, delimiter=',', skip_header=1, dtype=(str), encoding="latin-1")

    lats = indata[:, 2].astype(float)
    lons = indata[:, 3].astype(float)

    data = indata[:, 4].astype(float)

    return lats, lons, data # read_lakes

#************************************************************************
def read_hovmuller_style(filename):

    with open(filename, "r") as infile:
        
        times = np.arange(1993, int(settings.YEAR)+1)
        data = []
        names = []

        these_times = []
        these_data = []
        name = ""
        for line in infile:
            sline = line.split(",")
            if sline[0] == "lakename":
                continue

            else:
                if name != "" and name != sline[0]:
                    locs = np.in1d(times, np.array(these_times))
                    
                    matched_data = np.ones(times.shape[0])*-9999
                    matched_data[locs] = these_data

                    data += [matched_data]
                    names += [name]
                    
                    these_times = []
                    these_data = []

                # else just read in
                name = sline[0]
                these_times += [int(sline[1])]
                these_data += [float(sline[4])]

        # and last one
        locs = np.in1d(times, np.array(these_times))
                    
        matched_data = np.ones(times.shape[0])*-9999
        matched_data[locs] = these_data
        
        data += [matched_data]
        names += [name]

        # combine and mask
        data = np.array(data)
        data = np.ma.masked_where(data == -9999, data)
        names = np.array(names)

    return names, times, data # read_continuous_csv
                


#************************************************************************
def read_continuous_csv(filename):

    all_data = []

    with open(filename, "r") as infile:

        time = []
        data = []
        name = ""
        for line in infile:
            sline = line.split(",")
            if sline[0] == "date":
                continue
            else:
                if name != "" and name != sline[1]:
                    all_data += [utils.Timeseries(name, time, data)]
                    time = []
                    data = []
                # else just read in
                name = sline[1]
                time += [dt.datetime.strptime(sline[0], "%m/%d/%Y")]
                data += [float(sline[2])]
        # and last one
        all_data += [utils.Timeseries(name, time, data)]        
        
    return all_data # read_continuous_csv

#************************************************************************
def plot_smooth_map_iris(outname, cube, cmap, bounds, cb_label, scatter=[],\
                             figtext="", title=""):
    '''
    Standard scatter map

    :param str outname: output filename root
    :param array cube: cube to plot
    :param obj cmap: colourmap to use
    :param array bounds: bounds for discrete colormap
    :param str cb_label: label for colorbar
    :param str figtext: add text to figure label (top left)
    :param str title: add title
    '''

    norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)

    fig = plt.figure(figsize=(8, 5.5))
    plt.clf()
    ax = plt.axes([0.01, 0.12, 0.98, 0.88], projection=cartopy.crs.Robinson())

    ax.gridlines() #draw_labels=True)
    ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
    ax.coastlines()

    ext = ax.get_extent() # save the original extent

    if settings.OUTFMT in [".eps", ".pdf"]:
        if cube.coord("latitude").points.shape[0] > 180 or cube.coord("longitude").points.shape[0] > 360:
            regrid_size = 1.0
            print("Regridding cube for {} output to {} degree resolution".format(settings.OUTFMT, regrid_size))
            print("Old Shape {}".format(cube.data.shape))
            plot_cube = regrid_cube(cube, regrid_size, regrid_size)
            print("New Shape {}".format(plot_cube.data.shape))
        else:
            plot_cube = copy.deepcopy(cube)
    else:
        plot_cube = copy.deepcopy(cube)

    mesh = iris.plot.pcolormesh(plot_cube, cmap=cmap, norm=norm)

    lons, lats, data = scatter
    under = np.where(data <= 0)
    plt.scatter(lons[under], lats[under], c=data[under], cmap=cmap, norm=norm, s=40, marker="v",\
                transform=cartopy.crs.Geodetic(), edgecolor='0.2', linewidth=0.5, zorder=10)

    over = np.where(data > 0)
    plt.scatter(lons[over], lats[over], c=data[over], cmap=cmap, norm=norm, s=40, marker="^",\
                transform=cartopy.crs.Geodetic(), edgecolor='0.2', linewidth=0.5, zorder=10)


    cb = plt.colorbar(mesh, orientation='horizontal', pad=0.05, fraction=0.05, \
                        aspect=30, ticks=bounds[1:-1], label=cb_label, drawedges=True)

    # thicken border of colorbar and the dividers
    # http://stackoverflow.com/questions/14477696/customizing-colorbar-border-color-on-matplotlib
    cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
    cb.ax.tick_params(axis='x', labelsize=settings.FONTSIZE, direction='in', size=0)

    cb.set_label(label=cb_label, fontsize=settings.FONTSIZE)

    cb.outline.set_linewidth(2)
    cb.dividers.set_color('k')
    cb.dividers.set_linewidth(2)

    ax.set_extent(ext, ax.projection) # fix the extent change from colormesh

    plt.title(title, fontsize=settings.FONTSIZE)
    fig.text(0.01, 0.95, figtext, fontsize=settings.FONTSIZE)

    plt.savefig(outname + settings.OUTFMT)
    plt.close()

    return # plot_smooth_map_iris


#************************************************************************
def run_all_plots():

    if True:
    #***************
        # Anomaly Scatter map
        anomalies = read_lakes(DATALOC + "DATA_StateOfTheClimate{}_WaterLevelsAnom_v1.csv".format(settings.YEAR))

        bounds = [-2, -1, -0.75, -0.50, -0.25, 0, 0.25, 0.50, 0.75, 1, 2]
        bounds = [-5, -2, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2, 5]

        lons = np.arange(-90, 120, 30)
        lats = np.arange(-180, 210, 30)
        dummy = np.ma.zeros((len(lats), len(lons)))
        dummy.mask = np.ones(dummy.shape)

        cube = utils.make_iris_cube_2d(dummy, lats, lons, "blank", "m")

        # use local version to set different markers
        plot_smooth_map_iris(settings.IMAGELOC + "LWL_anomaly", cube, settings.COLOURMAP_DICT["hydrological"], \
                                       bounds, "Anomalies from 1993-2001 (m)", \
                                       scatter=[anomalies[1], anomalies[0], anomalies[2]], figtext="", title="")
        plot_smooth_map_iris(settings.IMAGELOC + "p2.1_LWL_{}_anomaly".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], \
                                       bounds, "Anomalies from 1993-2001 (m)", \
                                       scatter=[anomalies[1], anomalies[0], anomalies[2]], figtext="(m) Lake Water Level", title="")


    if True:
        #***************
        # multi-panel time series
        indata = read_continuous_csv(os.path.join(DATALOC, "DATA_StateOfTheClimate{}_WaterLevelsTimeSeriesPart2_v2.csv".format(settings.YEAR)))
        names = np.array([d.name for d in indata])

#        lake_order = ["Huron", "Superior", "Balkhash", "Tanganyika", "Caspian Sea", "Large Aral Sea", "Urmia", "Rukwa"]
        lake_order = ["Volta", "Victoria", "Tanganyika", "Huron/Michigan", "Superior", "Kara-Bogaz", "Caspian", "Kariba", "Urmia", "Aral"]

        # match HUM plot size

        fig, axes = plt.subplots(nrows=10, ncols=1, figsize=(8, 14), sharex=True)

        counter = 0
        for ax in axes:

            lake, = np.where(names == lake_order[counter])[0]
            lake = indata[lake]

            utils.thicken_panel_border(ax)
            ax.plot(lake.times, lake.data, "k-")

            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)

            ax.xaxis.set_major_locator(mdates.YearLocator(10))
            ax.xaxis.set_minor_locator(mdates.YearLocator(1))
            myFmt = mdates.DateFormatter('%Y')
            ax.xaxis.set_major_formatter(myFmt)

            ax.text(0.05, 0.8, lake.name, fontsize=settings.FONTSIZE, transform=ax.transAxes)

            ylim = ax.get_ylim()
            ax.set_ylim([ylim[0], ylim[1]*1.2])

            counter += 1

        fig.text(0.01, 0.5, "Lake Level (m)", rotation = "vertical", va="center", fontsize=settings.FONTSIZE)
        fig.subplots_adjust(left=0.1, right=0.99, bottom=0.04, top=0.99, hspace=0.001)
        plt.savefig(settings.IMAGELOC + "LWL_ts_{}{}".format(settings.YEAR, settings.OUTFMT))

    if False:
        #***************
        # Lake Cutouts

        EXTENTS = {"GL" : (86, 45, [-93, -78, 41, 49]), "RU" : (62, 42, [43, 81, 35, 48]), "AF" : (31, -6, [28.5, 33, -10, -3])}

        for region in ["GL", "RU", "AF"]:
            if region == "AF":
                fig = plt.figure(figsize=(6, 8))
            elif region == "RU":
                fig = plt.figure(figsize=(18, 8))
            elif region == "GL":
                fig = plt.figure(figsize=(12, 8))

            plt.clf()
            ax = plt.axes([0.05, 0.05, 0.9, 0.9], projection=cartopy.crs.PlateCarree(central_longitude=EXTENTS[region][0]))

            ax.set_extent(EXTENTS[region][2], cartopy.crs.PlateCarree())

            ax.add_feature(cartopy.feature.LAND.with_scale('10m'), zorder=0, facecolor="0.75", edgecolor="k")
            ax.add_feature(cartopy.feature.OCEAN.with_scale('10m'), zorder=0, facecolor="#e0ffff", edgecolor="k")
            ax.add_feature(cartopy.feature.LAKES.with_scale('10m'), zorder=0, facecolor="#e0ffff", edgecolor="k")
            ax.add_feature(cartopy.feature.RIVERS.with_scale('10m'), zorder=0, edgecolor="#e0ffff")
            ax.coastlines(resolution="10m", linewidth=0.5)

            # add other features
            gl = ax.gridlines(draw_labels=True)
            gl.xlabel_style = {'size': settings.FONTSIZE}
            gl.ylabel_style = {'size': settings.FONTSIZE}

            # add labels
            if region == "AF":
                ax.text(30, -5, "Tanganyika", fontsize=settings.FONTSIZE, transform=cartopy.crs.Geodetic())
                ax.text(32, -7.5, "Rukwa", fontsize=settings.FONTSIZE, transform=cartopy.crs.Geodetic())
            elif region == "RU":
                ax.text(46, 38, "Urmia", fontsize=settings.FONTSIZE, transform=cartopy.crs.Geodetic())
                ax.text(52, 44, "Caspian Sea", fontsize=settings.FONTSIZE, transform=cartopy.crs.Geodetic())
                ax.text(60.5, 45, "Large Aral Sea", fontsize=settings.FONTSIZE, transform=cartopy.crs.Geodetic())
                ax.text(75, 45.5, "Balkhash", fontsize=settings.FONTSIZE, transform=cartopy.crs.Geodetic())
            elif region == "GL":
                ax.text(-91, 46.2, "Superior", fontsize=settings.FONTSIZE, transform=cartopy.crs.Geodetic())
                ax.text(-89, 42, "Michigan", fontsize=settings.FONTSIZE, transform=cartopy.crs.Geodetic())
                ax.text(-84.5, 44.5, "Huron", fontsize=settings.FONTSIZE, transform=cartopy.crs.Geodetic())

            plt.savefig(settings.IMAGELOC + "LWL_map_{}".format(region) + settings.OUTFMT)
            plt.close()

    # hovmuller style
    if True:
        names, times, data = read_hovmuller_style(DATALOC + "DATA_StateOfTheClimate2020_WaterLevelsTimeSeries_v1.csv")

        # need to sort the data by the order of the current year anomalies
        this_year = data[:, -1]
        sort_order = np.argsort(this_year)

        data = data[sort_order, :]
        names = names[sort_order]

        # need n+1 times for pcolormesh
        times=np.append(times, int(settings.YEAR)+1)
        times = times-0.5 # centre on year      


        # set up color ranges
        bounds = [-8, -1.5, -1.0, -0.5, -0.25, 0.25, 0.5, 1.0, 1.5, 8]
        bounds = [-8, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 8]
        bounds = [-8, -2, -1.0, 0, 1, 2, 8]
        cmap = settings.COLOURMAP_DICT["hydrological"]
        norm=mpl.cm.colors.BoundaryNorm(bounds,cmap.N)
        
        fig = plt.figure(figsize = (8,14))
        plt.clf()

        # make axes by hand
        ax = plt.axes([0.17,0.04,0.81,0.93])

        this_cmap = copy.copy(cmap)
        this_cmap.set_bad("0.5", 1.)
        mesh = ax.pcolormesh(times, np.arange(data.shape[0]+1)-0.5, data, cmap=this_cmap, norm=norm)

        cb = plt.colorbar(mesh,  orientation='horizontal', ticks=bounds[1:-1],  drawedges=True, pad=0.02, fraction=0.03, aspect=20)

        # prettify
        cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
        cb.ax.tick_params(axis='x', labelsize=settings.FONTSIZE, direction='in', size=0)
        
        cb.set_label(label="Lake Water Level Anomaly (m)", fontsize=settings.FONTSIZE)
        cb.outline.set_linewidth(2)
        cb.dividers.set_color('k')
        cb.dividers.set_linewidth(2)

        minorLocator = MultipleLocator(1)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.set_xlim([1991.5, int(settings.YEAR)+1.5])
        ax.set_ylim([-5, 255])
        ax.set_ylabel("Lake", fontsize=settings.FONTSIZE, labelpad=-10)
        
        utils.thicken_panel_border(ax)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
#        for tick in ax.yaxis.get_major_ticks():
#            tick.set_visible(False)

        # set nice labels
        to_label = ["caspian","Aral_Sea_1","kara_bogaz_gol","Urmia_1","kariba", "tanganika","volta","superior","huron","michigan","victoria"]

        locs = np.in1d(names, to_label)

        ticks = np.arange(data.shape[0])
        ticks = ticks[locs]

        # process the tick names
        tick_names = []
        for rawname in names[ticks]:
            rawname = rawname.split("_")
            if rawname[-1] == "1":
                rawname = rawname[:-1]
            rawname = [rn.capitalize() for rn in rawname]
            tick_names += [" ".join(rawname)]

        ax.set_yticks(ticks)
        ax.set_yticklabels(tick_names)
        
        plt.savefig(settings.IMAGELOC + "LWL_hovmuller{}".format(settings.OUTFMT))

        plt.close()


    return # run_all_plots
#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
