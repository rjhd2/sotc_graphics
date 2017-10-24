#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for Surface Albedo (ABD) section.
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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import matplotlib.cm as mpl_cm
import matplotlib as mpl

import iris
import datetime as dt
import struct

import utils # RJHD utilities
import settings


data_loc = "/data/local/rdunn/SotC/{}/data/ABD/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

LEGEND_LOC = "lower left"
minor_tick_interval = 1
minorLocator = MultipleLocator(minor_tick_interval)

print "faking time axis - 16 day period"
# http://modis-atmos.gsfc.nasa.gov/ALBEDO/

start = dt.datetime(2003,1,1)
# get difference in days, 16 day period, and add a couple.
DURATION = int(round(((dt.datetime(int(settings.YEAR),12,31) - start).days)/16.))+2

dates=[start + dt.timedelta(days=i*16) for i in range(DURATION)]

times = np.ma.array([(d - start).days/365.  for d in dates]) + 2003
times.mask = np.zeros(times.shape)

#************************************************************************
def read_binary(filename):
    '''
    Read from binary file - using "struct"

    :param str filename: file to read
    :returns: np.array
    '''

    with open(filename, "rb") as f:
        fileContent = f.read()

    data = struct.unpack("d" * (len(fileContent) // 8), fileContent)

    return np.array(data) # read_binary

#************************************************************************
def read_binary_ts(filename):
    '''
    Read from binary file - using "struct"

    :param str filename: file to read
    :returns: np.array
    '''

    with open(filename, "rb") as f:
        fileContent = f.read()

    globe = np.array(struct.unpack("d" * (DURATION), fileContent[:DURATION*8]))
    north = np.array(struct.unpack("d" * (DURATION), fileContent[DURATION*8:DURATION*8*2]))
    south = np.array(struct.unpack("d" * (DURATION), fileContent[DURATION*8*2:DURATION*8*3]))
    aver = struct.unpack("f" * (DURATION*3), fileContent[DURATION*8*3:])

    aver = np.array(aver).reshape(3, DURATION)

    globe = utils.Timeseries("Globe", times, globe)
    north = utils.Timeseries("N. Hemisphere", times, north)
    south = utils.Timeseries("S. Hemisphere", times, south)

    globe_sm = utils.Timeseries("Globe Smoothed", times, aver[0,:])
    north_sm = utils.Timeseries("N. Hemisphere Smoothed", times, aver[1,:])
    south_sm = utils.Timeseries("S. Hemisphere Smoothed", times, aver[2,:])

    return  [globe, north, south, globe_sm, north_sm, south_sm] # read_binary_ts

#************************************************************************
def run_all_plots():

    #************************************************************************
    # Timeseries

    IRdata = read_binary_ts(data_loc + "TimeseriesBHRNIRpost_{}.bin".format(int(settings.YEAR)+1))
    Vdata = read_binary_ts(data_loc + "TimeseriesBHRVpost_{}.bin".format(int(settings.YEAR)+1))

    fig, (ax1, ax2) = plt.subplots(2, figsize = (10, 8), sharex=True)
    COLOURS = settings.COLOURS["land_surface"]

    for dataset in Vdata:
        print dataset.name
        locs, = np.where(dataset.data != 0) # remove zero bits (mainly for smoothed)
        ls = "--"
        lw = 1
        if dataset.name.split(" ")[-1] == "Smoothed":
            ls = "-"
            lw = 2
        ax1.plot(dataset.times[locs], dataset.data[locs], c = COLOURS[dataset.name], ls = ls, label = dataset.name, lw = lw)

    ax1.axhline(0, c = '0.5', ls = '--')
    utils.thicken_panel_border(ax1)

    for dataset in IRdata:
        print dataset.name
        locs, = np.where(dataset.data != 0)
        ls = "--"
        lw = 1
        if dataset.name.split(" ")[-1] == "Smoothed":
            ls = "-"
            lw = 2            
        ax2.plot(dataset.times[locs], dataset.data[locs], c = COLOURS[dataset.name], ls = ls, label = dataset.name, lw = lw)

    ax2.legend(loc = LEGEND_LOC, ncol = 2, frameon = False, prop={'size':settings.LEGEND_FONTSIZE * 0.8}, labelspacing = 0.1, columnspacing = 0.5)
    ax2.axhline(0, c = '0.5', ls = '--')
    utils.thicken_panel_border(ax2)

    #*******************
    # prettify

    fig.text(0.01, 0.5, "Normalised Anomalies (%)", va='center', rotation='vertical', fontsize = settings.FONTSIZE)

    plt.xlim([2003,2017])

    for ax in [ax1, ax2]:
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        ax.xaxis.set_minor_locator(minorLocator)
    for tick in ax2.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)

    ax1.text(0.02, 0.9, "(a) Visible", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE)
    ax2.text(0.02, 0.9, "(b) Near Infrared", transform = ax2.transAxes, fontsize = settings.LABEL_FONTSIZE)

    fig.subplots_adjust(right = 0.95, top = 0.95, bottom = 0.05, hspace = 0.001)

    plt.savefig(image_loc + "ABD_ts{}".format(settings.OUTFMT))
    plt.close()

    #************************************************************************
    # Hovmullers
    IRdata = read_binary(data_loc + "HovMullerBHRNIRpost_v{}.bin".format(int(settings.YEAR)+1))
    Vdata = read_binary(data_loc + "HovMullerBHRVpost_v{}.bin".format(int(settings.YEAR)+1))

    # reshape - from Readme
    IRdata = IRdata.reshape(360, DURATION)
    Vdata = Vdata.reshape(360, DURATION)

    IRdata = np.ma.masked_where(IRdata <= -100., IRdata) * 100 # from readme
    Vdata = np.ma.masked_where(Vdata <= -100., Vdata) * 100

    lats = np.arange(-90,90, 0.5)

    # curtail times
    locs, = np.where(times > int(settings.YEAR) +1 )
    IRdata.mask[:, locs] = True
    Vdata.mask[:, locs] = True
    times.mask[locs] = True

    bounds = [-100, -20, -15, -10, -5, -3, 3, 5, 10, 15, 20, 100]
    utils.plot_hovmuller(image_loc + "ABD_NIR_hovmuller", times, lats, IRdata, settings.COLOURMAP_DICT["land_surface_r"], bounds, "Normalised Anomalies (%)", figtext = "(b)")
    utils.plot_hovmuller(image_loc + "ABD_V_hovmuller", times, lats, Vdata, settings.COLOURMAP_DICT["land_surface_r"], bounds, "Normalised Anomalies (%)", figtext = "(a)")

    print "Still need to combine into single file output"

    #************************************************************************
    # Anomalies
    IRdata = read_binary(data_loc + "Annual_BHRNIRpost_v{}.bin".format(int(settings.YEAR)+1))
    Vdata = read_binary(data_loc + "Annual_BHRVpost_v{}.bin".format(int(settings.YEAR)+1))

    # reshape - from Readme
    IRdata = IRdata.reshape(360, 720)
    Vdata = Vdata.reshape(360, 720)

    IRdata = np.ma.masked_where(IRdata <= -9999., IRdata)
    Vdata = np.ma.masked_where(Vdata <= -9999., Vdata)

    delta = 0.5
    lons = np.arange(-180 + (delta/2.),180,delta)
    lats = np.arange(-90 + (delta/2.),90, delta)

    ir_cube = utils.make_iris_cube_2d(IRdata, lats, lons, "IR Albedo", "%")
    v_cube = utils.make_iris_cube_2d(Vdata, lats, lons, "V Albedo", "%")

    bounds = [-100, -20, -15, -10, -5, -3, 3, 5, 10, 15, 20, 100]

#   Section to show masks - for future reference if required.
#    import copy
#    this_cmap = copy.copy(settings.COLOURMAP_DICT["land_surface_r"])
#    this_cmap.set_bad("0.5", 1.)

    utils.plot_smooth_map_iris(image_loc + "p2.1_ABD_V_{}".format(settings.YEAR), v_cube, settings.COLOURMAP_DICT["land_surface_r"], bounds, "Anomalies from 2003-{} (%)".format(settings.YEAR), figtext = "(ab) Land Surface Albedo in the Visible")
    utils.plot_smooth_map_iris(image_loc + "ABD_V_{}".format(settings.YEAR), v_cube, settings.COLOURMAP_DICT["land_surface_r"], bounds, "Anomalies from 2003-{} (%)".format(settings.YEAR))

    utils.plot_smooth_map_iris(image_loc + "p2.1_ABD_NIR_{}".format(settings.YEAR), ir_cube, settings.COLOURMAP_DICT["land_surface_r"], bounds, "Anomalies from 2003-{} (%)".format(settings.YEAR), figtext = "(ac) Land Surface Albedo in the Near Infrared")
    utils.plot_smooth_map_iris(image_loc + "ABD_NIR_{}".format(settings.YEAR), ir_cube, settings.COLOURMAP_DICT["land_surface_r"], bounds, "Anomalies from 2003-{} (%)".format(settings.YEAR))

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************