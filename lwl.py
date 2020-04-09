#!/usr/bin/env python
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

    lats = indata[:, 1].astype(float)
    lons = indata[:, 2].astype(float)

    data = indata[:, 4].astype(float)

    return lats, lons, data # read_lakes

#************************************************************************
def read_continuous_csv(filename):

    NAMES = {"lake0082" : "Rukwa", "lake0115" : "Urmia", "lake0270" : "Caspian Sea", "lake0278" : "Balkhash", "lake0315" : "Tanganyika", "lake0336" : "Huron", "lake0337" : "Superior", "lake2327" : "Large Aral Sea"}
    
    all_data = []

    with open(filename, "r") as infile:

        time = []
        data = []
        name = ""
        for line in infile:
            sline = line.split(",")
            if sline[0] == "lakeID":
                continue
            else:
                if name != "" and name != sline[0]:
                    all_data += [utils.Timeseries(NAMES[name], time, data)]
                    time = []
                    data = []
                # else just read in
                name = sline[0]
                time += [dt.datetime.strptime(sline[1], "%Y-%m-%d")]
                data += [float(sline[2])]
        # and last one
        all_data += [utils.Timeseries(NAMES[name], time, data)]        
        
    return all_data # read_continuous_csv

#************************************************************************
def run_all_plots():

    #***************
    # Anomaly Scatter map
    anomalies = read_lakes(DATALOC + "DATA_StateOfTheClimate_WaterLevelsAnom_v1.csv")

    bounds = [-2, -1, -0.75, -0.50, -0.25, 0, 0.25, 0.50, 0.75, 1, 2]

    lons = np.arange(-90, 120, 30)
    lats = np.arange(-180, 210, 30)
    dummy = np.ma.zeros((len(lats), len(lons)))
    dummy.mask = np.ones(dummy.shape)

    cube = utils.make_iris_cube_2d(dummy, lats, lons, "blank", "m")

    utils.plot_smooth_map_iris(settings.IMAGELOC + "LWL_anomaly", cube, settings.COLOURMAP_DICT["hydrological"], \
                                   bounds, "Anomalies from 1992-2002 (m)", \
                                   scatter=[anomalies[1], anomalies[0], anomalies[2]], figtext="", title="")
    utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_LWL_{}_anomaly".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], \
                                   bounds, "Anomalies from 1992-2002 (m)", \
                                   scatter=[anomalies[1], anomalies[0], anomalies[2]], figtext="(m) Lake Water Level", title="")



    #***************
    # multi-panel time series
    indata = read_continuous_csv(os.path.join(DATALOC, "DATA_StateOfTheLakes_WaterLevelAnomTimeSeries_sub_v1.csv"))
    names = np.array([d.name for d in indata])

    lake_order = ["Huron", "Superior", "Balkhash", "Tanganyika", "Caspian Sea", "Large Aral Sea", "Urmia", "Rukwa"]

    # match HUM plot size

    fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(16, 8), sharex=True)

    counter = 0
    for row in axes:
        for ax in row:

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

            ax.text(0.1, 0.9, lake.name, fontsize=settings.FONTSIZE, transform=ax.transAxes)

            ylim = ax.get_ylim()
            ax.set_ylim([ylim[0], ylim[1]*1.2])
        
            counter += 1

    fig.text(0.01, 0.5, "Lake Level (m)", rotation = "vertical", va="center", fontsize=settings.FONTSIZE)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.05, top=0.99, hspace=0.001)
    plt.savefig(settings.IMAGELOC + "LWL_ts_{}{}".format(settings.YEAR, settings.OUTFMT))

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

    return # run_all_plots
#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
