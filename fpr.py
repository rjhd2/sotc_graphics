#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for Fraction of Absorbed Photosynthetic Radiation (FPR) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 30                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2021-06-15 10:41:02 +0100 (Tue, 15 Jun #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import struct
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import utils # RJHD utilities
import settings

DATALOC = "{}/{}/data/FPR/".format(settings.ROOTLOC, settings.YEAR)

LEGEND_LOC = 'lower left'
LW = 2
DURATION = (int(settings.YEAR)-1998+1)*12
minor_tick_interval = 1
minorLocator = MultipleLocator(minor_tick_interval)

print("faking time axis - monthly")
times = np.arange(0, DURATION/12., 1/12.) + 1998

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
    aver = struct.unpack("f" * (DURATION*4), fileContent[DURATION*8*3:])

    aver = np.array(aver).reshape(4, DURATION)

    globe = utils.Timeseries("Globe", times, globe)
    north = utils.Timeseries("N. Hemisphere", times, north)
    south = utils.Timeseries("S. Hemisphere", times, south)

    aver = np.ma.masked_where(aver == 0, aver) # remove first/last bits and line should span any inadvertent gaps

    globe_sm = utils.Timeseries("Globe Smoothed", times, aver[0, :])
    north_sm = utils.Timeseries("N. Hemisphere Smoothed", times, aver[1, :])
    south_sm = utils.Timeseries("S. Hemisphere Smoothed", times, aver[2, :])

    return  [globe, north, south, globe_sm, north_sm, south_sm] # read_binary_ts

#************************************************************************
def run_all_plots():
    #************************************************************************
    # Timeseries

    data = read_binary_ts(DATALOC + "TimeSeries_faparanomaliesglobal_bams_v2018_C6_{}.bin".format(int(settings.YEAR)+1))

    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes([0.16, 0.07, 0.82, 0.88])

    COLOURS = settings.COLOURS["land_surface"]

    for dataset in data:
        print(dataset.name)
        ls = "-."
        lw = 1
        if dataset.name.split(" ")[-1] == "Smoothed":
            ls = "-"
            lw = LW
        ax.plot(dataset.times, dataset.data, c=COLOURS[dataset.name], ls=ls, label=dataset.name, lw=lw)

    ax.axhline(0, c='0.5', ls='--')
    utils.thicken_panel_border(ax)

    ax.legend(loc=LEGEND_LOC, ncol=2, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, labelspacing=0.1, columnspacing=0.5)

    #*******************
    # prettify

    fig.text(0.01, 0.5, "Anomaly (FAPAR)", va='center', rotation='vertical', fontsize=settings.FONTSIZE)

    plt.xlim([1998-1, int(settings.YEAR)+2])
    plt.ylim([-0.022, None])

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)
    ax.xaxis.set_minor_locator(minorLocator)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)

    plt.savefig(settings.IMAGELOC + "FPR_ts{}".format(settings.OUTFMT))
    plt.close()

    #************************************************************************
    # Hovmullers
    data = read_binary(DATALOC + "Hovmuller_Global_lat_fapar1998_2010_bams_trois_C6.eps_{}.bin".format(int(settings.YEAR)+1))

    # reshape - from Readme
    data = data.reshape(360, DURATION)

    data = np.ma.masked_where(data == 0., data)

    lats = np.arange(-90, 90, 0.5)

    bounds = [-20, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 20]
    utils.plot_hovmuller(settings.IMAGELOC + "FPR_hovmuller", times, lats, data, settings.COLOURMAP_DICT["phenological"], bounds, "Anomaly (FAPAR)")


    #************************************************************************
    # Anomalies
    data = read_binary(DATALOC + "DataXFigure1fapar1998_2010_bams_trois_C6_v2.eps_v{}.bin".format(int(settings.YEAR)+1))
    data = data.reshape(360, 720)

    data = np.ma.masked_where(data == -100, data) # land/ocean mask
    data = np.ma.masked_where(data <= -9999., data) # missing data
    data = np.ma.masked_where(np.logical_and(data > -0.001, data < 0.001), data) # missing data

    lons = np.arange(-180, 180, 0.5)
    lats = np.arange(-90, 90, 0.5)

    cube = utils.make_iris_cube_2d(data, lats, lons, "FAPAR", "%")

    bounds = [-20, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 20]
    utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_FPR_{}".format(settings.YEAR), cube, settings.COLOURMAP_DICT["phenological"], bounds, "Anomalies from 1998-{} (FAPAR)".format(2010), figtext="(af) Fraction of Absorbed Photosynthetically Active Radiation")
    utils.plot_smooth_map_iris(settings.IMAGELOC + "FPR_{}".format(settings.YEAR), cube, settings.COLOURMAP_DICT["phenological"], bounds, "Anomalies from 1998-{} (FAPAR)".format(2010))


    return # run_all_plots


#************************************************************************
if __name__ == "__main__":

    run_all_plots()
#************************************************************************
#                                 END
#************************************************************************
