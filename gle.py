#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for global land evaporation (GLE) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 31                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2021-09-06 09:52:46 +0100 (Mon, 06 Sep #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************

import sys
import numpy as np
import matplotlib.pyplot as plt

import utils # RJHD utilities
import settings

DATALOC = "{}/{}/data/GLE/".format(settings.ROOTLOC, settings.YEAR)

LEGEND_LOC = 'lower center'

#************************************************************************
def read_time(filename):
    """
    Read user supplied CSV for GLE into Timeseries object
    """

    times = np.arange(1980, int(settings.YEAR)+1, 1)

    indata = np.genfromtxt(filename, delimiter=',', dtype=(float))

    globe = utils.Timeseries("Globe", times, indata[:, 0])
    NH = utils.Timeseries("N. Hemisphere", times, indata[:, 1])
    SH = utils.Timeseries("S. Hemisphere", times, indata[:, 2])
    SOI = utils.Timeseries("SOI", times, indata[:, 3])

    return globe, NH, SH, SOI # read_time

#************************************************************************
def read_hovmuller(filename):
    """
    Read user supplied CSV for GLE for Hovmuller
    """ 

    indata = np.genfromtxt(filename, delimiter=",", dtype=(float), \
                               missing_values="NaN", filling_values=-99.9)

    indata = np.ma.masked_where(indata == -99.9, indata)

    # Jan 1980 to December of report year
    times = np.arange(1980, int(settings.YEAR)+1, 1./12)

    # from email - 83.5N to 56 S at 0.25 degree resoln
    delta_lat = 0.25
    latitudes = np.arange(83.5 - (delta_lat/2.), -56, -delta_lat)
    latitudes = np.arange(90 - (delta_lat/2.), -90, -delta_lat)

    return times, latitudes, indata # read_hovmuller

#************************************************************************
def read_map(filename):
    """
    Read user supplied CSV for GLE for anomaly map
    """

    indata = np.genfromtxt(filename, delimiter=",", dtype=(float), \
                               missing_values="NaN", filling_values=-99.9)

    # 21 Feb 2018 - the missing_value/filling_value doesn't seem to work for some reason
    indata = np.ma.masked_where(indata == -99.9, indata)
    indata = np.ma.masked_where(indata == 0, indata)
    indata = np.ma.masked_where(indata != indata, indata) # catch any final NaNs

    # from email - 83.5N to 56 S at 0.25 degree resoln
    delta = 0.25
    lats = np.arange(90 - (delta/2.), -90, -delta)
    lons = np.arange(-180 + (delta/2.), 180, delta)

    cube = utils.make_iris_cube_2d(indata, lats, lons, "GLE", "mm/yr")

    return cube # read_map

#************************************************************************
def run_all_plots():

    #************************************************************************
    # Global Map
    if True:

        bounds = [-1000, -160, -120, -80, -40, 0, 40, 80, 120, 160, 1000]
        bounds = [-1000, -200, -150, -100, -50, 0, 50, 100, 150, 200, 1000]
        cube = read_map(DATALOC + "map")

        utils.plot_smooth_map_iris(settings.IMAGELOC + "GLE_{}_anoms".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (mm year"+r'$^{-1}$'+")")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_GLE_{}_anoms".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (mm year"+r'$^{-1}$'+")", figtext="(t) Land Evaporation")

    # Evaporation Map
    if False:

        bounds = [-1000, -160, -120, -80, -40, 0, 40, 80, 120, 160, 1000]
        cube = read_map(DATALOC + "map_transp")

        utils.plot_smooth_map_iris(settings.IMAGELOC + "GLE_{}_transp".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (mm year"+r'$^{-1}$'+")")
 

    #************************************************************************
    # Timeseries figures
    if True:
        fig, ax1 = plt.subplots(figsize=(8, 5))

        globe, NH, SH, SOI = read_time(DATALOC + "timeseries")
        globe.zorder=10
        NH.zorder=10
        SH.zorder=10
        SOI.zorder=1

        ax2 = ax1.twinx()

        # SOI - plot first so under the other lines

        interpTimes = np.linspace(SOI.times[0], SOI.times[-1], 1000)
        interpData = np.interp(interpTimes, SOI.times, SOI.data)
        interpSOI = utils.Timeseries("SOI", interpTimes, interpData)

        ax1.fill_between(interpSOI.times, interpSOI.data, where=interpSOI.data >= 0, \
                             color='lightskyblue', zorder=1)
        ax1.fill_between(interpSOI.times, interpSOI.data, where=interpSOI.data <= 0, \
                             color='lightcoral', zorder=1)

        # then plot the GLE timeseries
        utils.plot_ts_panel(ax2, [globe, NH, SH], "-", "hydrological", loc=LEGEND_LOC, ncol=3)

        for data in [globe, NH, SH]:
            slope, dummy, dummy = utils.median_pairwise_slopes(data.times, data.data, -99.9, 1.)

            fit_years, fit_values = utils.mpw_plot_points(slope, data.times, data.data)

            ax2.plot(fit_years, fit_values, c=settings.COLOURS["hydrological"][data.name], \
                         lw=2, ls="--", zorder=10)

        ax2.set_ylabel("Anomalies (mm year"+r'$^{-1}$'+")", fontsize=settings.FONTSIZE)    
        ax2.set_ylim([-29, 29])

        for tick in ax1.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 

        ax1.set_xlim([1978, int(settings.YEAR)+2])
        ax1.set_ylim([-3, 3])
        ax1.set_ylabel("SOI", fontsize=settings.FONTSIZE)

        # apply to both axes
        utils.thicken_panel_border(ax2)
        utils.thicken_panel_border(ax1)

        # and swap the labels and ticks around from standard
        ax1.yaxis.tick_right()
        ax2.yaxis.tick_left()
        ax1.yaxis.set_label_position("right")
        ax2.yaxis.set_label_position("left")
        
        for tick in ax2.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
        for tick in ax1.yaxis.get_major_ticks():
            tick.label2.set_fontsize(settings.FONTSIZE) 

        fig.subplots_adjust(bottom=0.1, top=0.95, hspace=0.001)
        plt.savefig(settings.IMAGELOC+"GLE_ts{}".format(settings.OUTFMT))

        plt.close()

    #************************************************************************
    # Hoemuller figure
    if True:
        bounds = [-100, -10, -8, -4, -2, 0, 2, 4, 8, 10, 100]
        times, latitudes, indata = read_hovmuller(DATALOC + "latitudinal")
        
        utils.plot_hovmuller(settings.IMAGELOC + "GLE_hovmuller", times, latitudes, indata, \
                             settings.COLOURMAP_DICT["hydrological"], bounds, \
                             "Anomaly (mm month"+r'$^{-1}$'+")")


    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
