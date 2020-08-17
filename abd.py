#!/usr/env/bin python
#************************************************************************
#
#  Plot figures and output numbers for Surface Albedo (ABD) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 28                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2020-04-09 11:37:08 +0100 (Thu, 09 Apr #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
# python3
from __future__ import absolute_import
from __future__ import print_function
import datetime as dt
import struct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import utils # RJHD utilities
import settings


DATALOC = "{}/{}/data/ABD/".format(settings.ROOTLOC, settings.YEAR)


LEGEND_LOC = "lower left"
minor_tick_interval = 1
minorLocator = MultipleLocator(minor_tick_interval)

# http://modis-atmos.gsfc.nasa.gov/ALBEDO/

start = dt.datetime(2003, 1, 1)
# get difference in days, 10 day period, and add a few.
PERIOD = 10.
print("faking time axis - {} day period".format(PERIOD))
print("Built-in round() might cause issues in Python3")
DURATION = int(round(((dt.datetime(int(settings.YEAR), 12, 31) - start).days)/PERIOD))-9

dates = [start + dt.timedelta(days=i*PERIOD) for i in range(DURATION)]

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

    globe_sm = utils.Timeseries("Globe Smoothed", times, aver[0, :])
    north_sm = utils.Timeseries("N. Hemisphere Smoothed", times, aver[1, :])
    south_sm = utils.Timeseries("S. Hemisphere Smoothed", times, aver[2, :])

    return  [globe, north, south, globe_sm, north_sm, south_sm] # read_binary_ts

#************************************************************************
def run_all_plots():

    #************************************************************************
    # Timeseries
    if True:

        IRdata = read_binary_ts(DATALOC + "TimeseriesBHRNIR_C7_poids_{}.bin".format(int(settings.YEAR)+1))
        Vdata = read_binary_ts(DATALOC + "TimeseriesBHRV_C7_poids_{}.bin".format(int(settings.YEAR)+1))

        fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 6.5), sharex=True)
        COLOURS = settings.COLOURS["land_surface"]

        for dataset in Vdata:
            print(dataset.name)
            locs, = np.where(dataset.data != 0)
            ls = "--"
            lw = 1
            if dataset.name.split(" ")[-1] == "Smoothed":
                ls = "-"
                lw = 2            
            ax1.plot(dataset.times[locs], dataset.data[locs], c=COLOURS[dataset.name], ls=ls, label=dataset.name, lw=lw)

        ax1.axhline(0, c='0.5', ls='--')
        utils.thicken_panel_border(ax1)

        for dataset in IRdata:
            print(dataset.name)
            locs, = np.where(dataset.data != 0) # remove zero bits (mainly for smoothed)
            ls = "--"
            lw = 1
            if dataset.name.split(" ")[-1] == "Smoothed":
                ls = "-"
                lw = 2
            ax2.plot(dataset.times[locs], dataset.data[locs], c=COLOURS[dataset.name], ls=ls, label=dataset.name, lw=lw)

        ax2.legend(loc=LEGEND_LOC, ncol=2, frameon=False, prop={'size':settings.LEGEND_FONTSIZE * 0.8}, labelspacing=0.1, columnspacing=0.5)
        ax2.axhline(0, c='0.5', ls='--')
        utils.thicken_panel_border(ax2)

        #*******************
        # prettify

        fig.text(0.01, 0.5, "Normalized Anomalies (%)", va='center', rotation='vertical', fontsize=settings.FONTSIZE)

        plt.xlim([2002.5, int(settings.YEAR)+1.5])
        ax2.set_ylim([-5.8, None])
        for ax in [ax1, ax2]:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)
            ax.xaxis.set_minor_locator(minorLocator)
            ax.set_yticks(ax.get_yticks()[1:-1])
            ax.yaxis.set_ticks_position('left')
        for tick in ax2.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)

        ax1.text(0.02, 0.9, "(a) Visible", transform=ax1.transAxes, fontsize=settings.LABEL_FONTSIZE)
        ax2.text(0.02, 0.9, "(b) Near Infrared", transform=ax2.transAxes, fontsize=settings.LABEL_FONTSIZE)

        fig.subplots_adjust(right=0.95, top=0.95, bottom=0.05, hspace=0.001)

        plt.savefig(settings.IMAGELOC + "ABD_ts{}".format(settings.OUTFMT))
        plt.close()

    #************************************************************************
    # Hovmullers
    if True:
        IRdata = read_binary(DATALOC + "HovMullerBHRNIR_C7_{}.bin".format(int(settings.YEAR)+1))
        Vdata = read_binary(DATALOC + "HovMullerBHRV_C7_{}.bin".format(int(settings.YEAR)+1))

        # reshape - from Readme
        IRdata = IRdata.reshape(360, DURATION)
        Vdata = Vdata.reshape(360, DURATION)

        IRdata = np.ma.masked_where(IRdata <= -100., IRdata)# * 100 # from readme
        Vdata = np.ma.masked_where(Vdata <= -100., Vdata)# * 100

        lats = np.arange(-90, 90, 0.5)

        # curtail times
        locs, = np.where(times > int(settings.YEAR) + 1)
        IRdata.mask[:, locs] = True
        Vdata.mask[:, locs] = True
        times.mask[locs] = True

        # cant use "cmap.set_bad" for contour as ignored when contouring

        bounds = [-100, -20, -15, -10, -5, -1, 1, 5, 10, 15, 20, 100]
        utils.plot_hovmuller(settings.IMAGELOC + "ABD_NIR_hovmuller", times, lats, IRdata, settings.COLOURMAP_DICT["land_surface_r"], bounds, "Normalized Anomalies (%)", figtext="(b)", background="0.7")
        utils.plot_hovmuller(settings.IMAGELOC + "ABD_V_hovmuller", times, lats, Vdata, settings.COLOURMAP_DICT["land_surface_r"], bounds, "Normalized Anomalies (%)", figtext="(a)", background="0.7")

    #************************************************************************
    # Anomalies
    if True:
        IRdata = read_binary(DATALOC + "Annual_BHRNIR_C7_{}.bin".format(int(settings.YEAR)+1))
        Vdata = read_binary(DATALOC + "Annual_BHRV_C7_{}.bin".format(int(settings.YEAR)+1))

        # reshape - from Readme
        IRdata = IRdata.reshape(360, 720)
        Vdata = Vdata.reshape(360, 720)

        IRdata = np.ma.masked_where(IRdata <= -9999., IRdata)
        Vdata = np.ma.masked_where(Vdata <= -9999., Vdata)

        IRdata = np.ma.masked_where(IRdata == 0., IRdata)
        Vdata = np.ma.masked_where(Vdata == 0., Vdata)

        delta = 0.5
        lons = np.arange(-180 + (delta/2.), 180, delta)
        lats = np.arange(-90 + (delta/2.), 90, delta)

        ir_cube = utils.make_iris_cube_2d(IRdata, lats, lons, "IR Albedo", "%")
        v_cube = utils.make_iris_cube_2d(Vdata, lats, lons, "V Albedo", "%")

        bounds = [-100, -20, -15, -10, -5, 0, 5, 10, 15, 20, 100]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_ABD_V_{}".format(settings.YEAR), v_cube, settings.COLOURMAP_DICT["land_surface_r"], bounds, "Anomalies from 2003-{} (%)".format(2010), figtext="(ac) Land Surface Albedo in the Visible")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "ABD_V_{}".format(settings.YEAR), v_cube, settings.COLOURMAP_DICT["land_surface_r"], bounds, "Anomalies from 2003-{} (%)".format(2010))

        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_ABD_NIR_{}".format(settings.YEAR), ir_cube, settings.COLOURMAP_DICT["land_surface_r"], bounds, "Anomalies from 2003-{} (%)".format(2010), figtext="(ad) Land Surface Albedo in the Near Infrared")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "ABD_NIR_{}".format(settings.YEAR), ir_cube, settings.COLOURMAP_DICT["land_surface_r"], bounds, "Anomalies from 2003-{} (%)".format(2010))

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
