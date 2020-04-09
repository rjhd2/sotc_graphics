#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for Biomass Burning (BOB) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 27                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2019-08-15 16:09:25 +0100 (Thu, 15 Aug #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
# python3
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import iris

import utils # RJHD utilities
import settings

DATALOC = "{}/{}/data/BOB/".format(settings.ROOTLOC, settings.YEAR)

LEGEND_LOC = 'upper right'
GFAS_FACTOR = 1.14

#************************************************************************
def read_data(filename, mdi):
    '''
    Read text format data into data array and return a cube
    Not used in SotC2016

    :param str filename: file to read

    :returns: cube
    '''    

    data = np.genfromtxt(filename, dtype=(float))

    nlats = data.shape[0]
    nlons = data.shape[1]

    latitudes = np.arange(90, -90, -180./nlats)
    longitudes = np.arange(0, 360, 360./nlons)

    data = np.ma.masked_where(data == mdi, data)

    cube = utils.make_iris_cube_2d(data, latitudes, longitudes, "BOB", "g C/m2/yr")

    return cube # read_data

#************************************************************************
def read_csv(filenameroot, name, make_annual=False):
    """
    Read user supplied CSV for LST into Timeseries object
    """

    indata = np.genfromtxt("{}_{}.dat".format(filenameroot, name), delimiter=' ', dtype=(str), skip_header=1)
   
    years = np.array([i[:4] for i in indata[:, 0]]).astype(int)
    months = np.array([i[-2:] for i in indata[:, 0]]).astype(int)

    times = years + (months - 1)/12.

    if name == "GFAS1p4":
        name = "GFASv1.4"
    elif name == "GFAS1p0":
        name = "GFASv1.0"

    if make_annual:

        monthly = indata[:, 1].astype(float)

        annual = np.sum(monthly.reshape(-1, 12), axis=1)/1000.

        years = (years.reshape(-1, 12)[:, 0])

        if years[-1] == int(settings.YEAR) + 1:
            print("skipping final year {} ".format(years[-1]))
            timeseries = utils.Timeseries(name, years[:-1], annual[:-1])
        else:
            timeseries = utils.Timeseries(name, years, annual)
    else:

        timeseries = utils.Timeseries(name, times, indata[:, 1].astype(float))

    return timeseries # read_csv

#************************************************************************
def read_gfas_csv(filename, name, make_annual=False):
    """
    Read user supplied CSV for LST into Timeseries object
    """

    indata = np.genfromtxt(filename, dtype=(str), skip_header=2)
   
    years = indata[:, 0].astype(int)
    months = indata[:, 1].astype(int)

    times = years + (months - 1)/12.

    if name == "global":
        column = 2
    elif name == "NAme":
        column = 3

    print("remove 1.14 factor after 2017 run")
    if make_annual:

        monthly = 1.14*indata[:, column].astype(float)

        annual = np.sum(monthly.reshape(-1, 12), axis=1)/1000.

        years = (years.reshape(-1, 12)[:, 0])

        if years[-1] == int(settings.YEAR) + 1:
            print("skipping final year {} ".format(years[-1]))
            timeseries = utils.Timeseries("GFASv1.4", years[:-1], annual[:-1]/1.e9)
        else:
            timeseries = utils.Timeseries("GFASv1.4", years, annual/1.e9)
    else:

        timeseries = utils.Timeseries("GFASv1.4", times, GFAS_FACTOR*indata[:, column].astype(float)/1.e9)

    return timeseries # read_gfas_csv

#************************************************************************
def read_gfed_csv(filename, name, make_annual=False):
    """
    Read user supplied CSV for LST into Timeseries object
    """

    indata = np.genfromtxt(filename, dtype=(str))
   
    years = indata[:, 0].astype(float)
    months = indata[:, 1].astype(float)

    times = years + (months - 1)/12.

    if name == "global":
        column = 2
    elif name == "NAme":
        column = 3

    if make_annual:

        monthly = indata[:, column].astype(float)

        annual = np.sum(monthly.reshape(-1, 12), axis=1)/1000.

        years = (years.reshape(-1, 12)[:, 0])

        if years[-1] == int(settings.YEAR) + 1:
            print("skipping final year {} ".format(years[-1]))
            timeseries = utils.Timeseries("GFED4s", years[:-1], annual[:-1])
        else:
            timeseries = utils.Timeseries("GFED4s", years, annual)
    else:

        timeseries = utils.Timeseries("GFED4s", times, indata[:, column].astype(float))

    return timeseries # read_gfed_csv

#************************************************************************
def run_all_plots():


    #************************************************************************
    # Timeseries - 2016 version
    if False:
        fig, (ax1, ax2) = plt.subplots(2, figsize = (8, 6.5), sharex=True)

        gfed = read_csv(DATALOC + "timeseries_NAme", "GFED3.1")
        gfas3 = read_csv(DATALOC + "timeseries_NAme", "GFAS1p3")
        gfas0 = read_csv(DATALOC + "timeseries_NAme", "GFAS1p0")

        utils.plot_ts_panel(ax1, [gfed, gfas3, gfas0], "-", "land_surface", loc = LEGEND_LOC)

        ax1.set_ylabel("Tg(C) per month", fontsize = settings.FONTSIZE)
        ax1.set_ylim([0,100])

        gfed = read_csv(DATALOC + "timeseries_TAsi", "GFED3.1")
        gfas3 = read_csv(DATALOC + "timeseries_TAsi", "GFAS1p3")
        gfas0 = read_csv(DATALOC + "timeseries_TAsi", "GFAS1p0")

        utils.plot_ts_panel(ax2, [gfed, gfas3, gfas0], "-", "land_surface", loc = LEGEND_LOC)

        ax2.set_ylabel("Tg(C) per month", fontsize = settings.FONTSIZE)
        ax2.set_ylim([0,450])

        # sort formatting
        plt.xlim([1997,2017])

        for tick in ax2.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
        for ax in [ax1, ax2]:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE) 

        # sort labelling
        ax1.text(0.02, 0.9, "(a) North America", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE)
        ax2.text(0.02, 0.9, "(b) Tropical Asia", transform = ax2.transAxes, fontsize = settings.LABEL_FONTSIZE)

        fig.subplots_adjust(right = 0.95, top = 0.95, hspace = 0.001)

        plt.savefig(settings.IMAGELOC+"BOB_ts{}".format(settings.OUTFMT))

        plt.close()

    #************************************************************************
    # Timeseries - 2017, 2018 version
    if True:
        fig = plt.figure(figsize=(8, 5))
        plt.clf()
        ax = plt.axes([0.11, 0.09, 0.88, 0.89])

        gfed = read_gfed_csv(DATALOC + "data4Johannes.txt", "global")
        gfas = read_csv(DATALOC + "timeseries_glob", "GFAS1p4")

        utils.plot_ts_panel(ax, [gfed, gfas], "-", "land_surface", loc=LEGEND_LOC)

        ax.set_ylabel("Tg(C) month"+r'$^{-1}$', fontsize=settings.FONTSIZE)
        ax.set_ylim([0, 840])

        # sort formatting
        plt.xlim([1997, int(settings.YEAR) + 2])

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 

        # sort labelling
        ax.text(0.02, 0.9, "Global", transform=ax.transAxes, fontsize=settings.LABEL_FONTSIZE)

        plt.savefig(settings.IMAGELOC+"BOB_ts{}".format(settings.OUTFMT))
        plt.close()

    #************************************************************************
    # Timeseries - 2019 version
    if True:
        majorLocator = MultipleLocator(5)
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize = (8, 10), sharex=True)

        data = read_csv(DATALOC + "timeseries_TAsi", "GFAS1p4")
        utils.plot_ts_panel(ax1, [data], "-", "land_surface", loc = "")
        ax1.text(0.02, 0.9, "(a) Tropical Asia", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE)
        ax1.set_ylim([0,160])

        data = read_csv(DATALOC + "timeseries_Arct", "GFAS1p4")
        utils.plot_ts_panel(ax2, [data], "-", "land_surface", loc = "")
        ax2.text(0.02, 0.9, "(b) Arctic", transform = ax2.transAxes, fontsize = settings.LABEL_FONTSIZE)
        ax2.set_ylim([0,12])
 
        data = read_csv(DATALOC + "timeseries_SEAu", "GFAS1p4")
        utils.plot_ts_panel(ax3, [data], "-", "land_surface", loc = "")
        ax3.text(0.02, 0.9, "(c) NSW & Victoria", transform = ax3.transAxes, fontsize = settings.LABEL_FONTSIZE)
        ax3.set_ylim([0,16])

        data = read_csv(DATALOC + "timeseries_SAme", "GFAS1p4")
        utils.plot_ts_panel(ax4, [data], "-", "land_surface", loc = "")
        ax4.text(0.02, 0.9, "(d) Southern America", transform = ax4.transAxes, fontsize = settings.LABEL_FONTSIZE)
        ax4.set_ylim([0,160])

        # sort formatting
        plt.xlim([2003, int(settings.YEAR)+2])
        fig.text(0.03, 0.5, "Tg(C) per month", fontsize = settings.FONTSIZE, rotation="vertical")
        ax4.xaxis.set_major_locator(majorLocator)
        for tick in ax4.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
        for ax in [ax1, ax2, ax3, ax4]:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE) 

        # sort labelling

        fig.subplots_adjust(right = 0.96, top = 0.98, bottom=0.05, hspace = 0.001)

        plt.savefig(settings.IMAGELOC+"BOB_regional_ts{}".format(settings.OUTFMT))

        plt.close()

    #************************************************************************
    # Global Anomalies
    if True:
    #    cube = read_data(DATALOC + "anomaly{}from2003-2015_smooth.txt".format(settings.YEAR), 2.1021584473146504e-05)
    #    cube_list = iris.load(DATALOC + "anomaly{}from2003-{}_coarse20.grb".format(settings.YEAR, int(settings.YEAR)-1))
        cube_list = iris.load(DATALOC + "anomaly{}from2003-{}_coarse20.grb".format(settings.YEAR, 2010))

        # 2016 - these grib files have 6 identical fields (email from J Kaiser 10/3/2017)
        cube = cube_list[0]
        # convert from kg(C) m-2 s-1 to g(C) m-2 a-1
        print("remove factor 1.14 after 2018 run")
        cube.data = cube.data * 1000 * 60 * 60 * 24 * 365. * GFAS_FACTOR


        bounds = [-1000, -160, -80, -40, -10, -5, 5, 10, 40, 80, 160, 1000]
        bounds = [-1000, -100, -40, -10, -5, -1, 1, 5, 10, 40, 100, 1000]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "BOB_anomalies{}".format(settings.YEAR), cube, settings.COLOURMAP_DICT["land_surface"], bounds, "Anomalies from 2003-{} (g C m".format(2010)+r'$^{-2}$'+" yr"+r'$^{-1}$'+")")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_BOB_anomalies{}".format(settings.YEAR), cube, settings.COLOURMAP_DICT["land_surface"], bounds, "Anomalies from 2003-{} (g C m".format(2010)+r'$^{-2}$'+" yr"+r'$^{-1}$'+")", figtext="(af) Carbon Emissions from Biomass Burning")

    #************************************************************************
    # Global actuals
    if True:
    #    cube = read_data(DATALOC + "year2016_smooth.txt", 0)
        cube_list = iris.load(DATALOC + "year{}_coarse20.grb".format(settings.YEAR))
        cube = cube_list[0]
        # convert from kg(C) m-2 s-1 to g(C) m-2 a-1
        print("remove factor 1.14 after 2018 run")
        cube.data = cube.data * 1000 * 60 * 60 * 24 * 365. * GFAS_FACTOR

        bounds = [0, 5, 10, 40, 80, 120, 160, 200, 240, 500]
        bounds = [0, 1, 5, 10, 40, 80, 120, 160, 200, 500]

        # adjust to mask out <1 gC/m2/yr (Feb 2018)
        cmap = plt.cm.YlOrBr
        cmaplist = [cmap(i) for i in range(cmap.N)]
        for i in range(30):
            cmaplist[i] = (0.75, 0.75, 0.75, 1.0)
        cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

        cube.data = np.ma.masked_where(cube.data < 1, cube.data)

        print("using hard coded sequential colourmap")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "BOB_actuals{}".format(settings.YEAR), cube, cmap, bounds, "g C m"+r'$^{-2}$'+" yr"+r'$^{-1}$'+"")

    #************************************************************************
    # Global climatology
    if False:
        cube_list = iris.load(DATALOC + "clim2003-2015_smooth.grb")
        cube = cube_list[0]
        # convert from kg(C) m-2 s-1 to g(C) m-2 a-1
        cube.data = cube.data * 1000 * 60 * 60 * 24 * 365.

        bounds = [0,5,10,40,80,120,160,200,240,500]

        print("using hard coded sequential colourmap")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "BOB_climatology{}".format(settings.YEAR), cube, plt.cm.YlOrBr, bounds, "g C m"+r'$^{-2}$'+" yr"+r'$^{-1}$'+"")



    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
