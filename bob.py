#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for Biomass Burning (BOB) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev::                                          $:  Revision of last commit
# $Author::                                       $:  Author of last commit
# $Date::                                         $:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.cm as mpl_cm
import matplotlib as mpl

import iris
import datetime as dt
import struct

import utils # RJHD utilities
import settings

data_loc = "/data/local/rdunn/SotC/{}/data/BOB/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

LEGEND_LOC = 'upper right'

#************************************************************************
def read_data(filename, mdi):
    '''
    Read text format data into data array and return a cube
    Not used in SotC2016

    :param str filename: file to read

    :returns: cube
    '''    

    data = np.genfromtxt(filename, dtype = (float))

    nlats = data.shape[0]
    nlons = data.shape[1]

    latitudes = np.arange(90,-90, -180./nlats)
    longitudes = np.arange(0,360, 360./nlons)

    data = np.ma.masked_where(data == mdi, data)

    cube = utils.make_iris_cube_2d(data, latitudes, longitudes, "BOB", "g C/m2/yr")

    return cube # read_data

#************************************************************************
def read_csv(filenameroot, name, make_annual = False):
    """
    Read user supplied CSV for LST into Timeseries object
    """

    indata = np.genfromtxt("{}_{}.dat".format(filenameroot, name), delimiter = ' ', dtype=(str), skip_header = 1)
   
    years = np.array([i[:4] for i in indata[:,0]]).astype(int)
    months = np.array([i[-2:] for i in indata[:,0]]).astype(int)

    times = years + (months - 1)/12.

    if name == "GFAS1p3":
        name = "GFASv1.3"
    elif name == "GFAS1p0":
        name = "GFASv1.0"

    if make_annual:

        monthly = indata[:,1].astype(float)

        annual = np.sum(monthly.reshape(-1,12), axis = 1)/1000.

        years = (years.reshape(-1,12)[:,0])

        if years[-1] == int(settings.YEAR) + 1:
            print "skipping final year {} ".format(years[-1])
            timeseries = utils.Timeseries(name, years[:-1], annual[:-1])
        else:
            timeseries = utils.Timeseries(name, years, annual)
    else:

        timeseries = utils.Timeseries(name, times, indata[:,1].astype(float))

    return timeseries # read_csv

#************************************************************************
def run_all_plots():


    #************************************************************************
    # Timeseries

    fig, (ax1, ax2) = plt.subplots(2, figsize = (10,8), sharex=True)

    gfed = read_csv(data_loc + "timeseries_NAme", "GFED3.1")
    gfas3 = read_csv(data_loc + "timeseries_NAme", "GFAS1p3")
    gfas0 = read_csv(data_loc + "timeseries_NAme", "GFAS1p0")

    utils.plot_ts_panel(ax1, [gfed, gfas3, gfas0], "-", "land_surface", loc = LEGEND_LOC)

    ax1.set_ylabel("Tg(C) per month", fontsize = settings.FONTSIZE)
    ax1.set_ylim([0,100])

    gfed = read_csv(data_loc + "timeseries_TAsi", "GFED3.1")
    gfas3 = read_csv(data_loc + "timeseries_TAsi", "GFAS1p3")
    gfas0 = read_csv(data_loc + "timeseries_TAsi", "GFAS1p0")

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

    plt.savefig(image_loc+"BOB_ts{}".format(settings.OUTFMT))

    plt.close()

    #************************************************************************
    # Global Anomalies

#    cube = read_data(data_loc + "anomaly{}from2003-2015_smooth.txt".format(settings.YEAR), 2.1021584473146504e-05)
    cube_list = iris.load(data_loc + "anomaly{}from2003-2015_smooth.grb".format(settings.YEAR))

    # 2016 - these grib files have 6 identical fields (email from J Kaiser 10/3/2017)
    cube = cube_list[0]
    # convert from kg(C) m-2 s-1 to g(C) m-2 a-1
    cube.data = cube.data * 1000 * 60 * 60 * 24 * 365.

    bounds = [-1000, -160, -80, -40, -10, -5, 5, 10, 40, 80, 160, 1000]

    utils.plot_smooth_map_iris(image_loc + "BOB_anomalies{}".format(settings.YEAR), cube, settings.COLOURMAP_DICT["land_surface"], bounds, "Anomalies from 2003-15 (g C m"+r'$^{-2}$'+" yr"+r'$^{-1}$'+")")
    utils.plot_smooth_map_iris(image_loc + "p2.1_BOB_anomalies{}".format(settings.YEAR), cube, settings.COLOURMAP_DICT["land_surface"], bounds, "Anomalies from 2003-15 (g C m"+r'$^{-2}$'+" yr"+r'$^{-1}$'+")", figtext = "(ad) Carbon Emissions from Biomass Burning")

    #************************************************************************
    # Global actuals

#    cube = read_data(data_loc + "year2016_smooth.txt", 0)
    cube_list = iris.load(data_loc + "year{}_smooth.grb".format(settings.YEAR))
    cube = cube_list[0]
    # convert from kg(C) m-2 s-1 to g(C) m-2 a-1
    cube.data = cube.data * 1000 * 60 * 60 * 24 * 365.

    bounds = [0,5,10,40,80,120,160,200,240,500]

    print "using hard coded sequential colourmap"
    utils.plot_smooth_map_iris(image_loc + "BOB_actuals{}".format(settings.YEAR), cube, plt.cm.YlOrBr, bounds, "g C m"+r'$^{-2}$'+" yr"+r'$^{-1}$'+"")

    #************************************************************************
    # Global climatology

    cube_list = iris.load(data_loc + "clim2003-2015_smooth.grb")
    cube = cube_list[0]
    # convert from kg(C) m-2 s-1 to g(C) m-2 a-1
    cube.data = cube.data * 1000 * 60 * 60 * 24 * 365.

    bounds = [0,5,10,40,80,120,160,200,240,500]

    print "using hard coded sequential colourmap"
    utils.plot_smooth_map_iris(image_loc + "BOB_climatology{}".format(settings.YEAR), cube, plt.cm.YlOrBr, bounds, "g C m"+r'$^{-2}$'+" yr"+r'$^{-1}$'+"")



    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
