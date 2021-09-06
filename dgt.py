#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for Drought (DGT) section.
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
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator

import iris
import sys
import datetime as dt
import numpy as np

import utils # RJHD utilities
import settings

DATALOC = "{}/{}/data/DGT/".format(settings.ROOTLOC, settings.YEAR)


#************************************************************************
def read_data(filename):

    indata = np.genfromtxt(filename, skip_header=2, delimiter=",")
    
    # process the years and months to give decimals
    years = indata[:, 0]
    months = indata[:, 1]

    times = years + (months - 1)/12.

    # make timeseries objects
    moderate = utils.Timeseries("Moderate", times, indata[:, 2])
    severe = utils.Timeseries("Severe", times, indata[:, 3])
    extreme = utils.Timeseries("Extreme", times, indata[:, 4])

    return moderate, severe, extreme # read_data


#************************************************************************
def run_all_plots():

    if False:
        #************************************************************************
        # Filled Timeseries
        #************************************************************************


        moderate, severe, extreme = read_data(DATALOC + "dry_areas_bams_2018_series.csv")

        plt.clf()

        # main panel
        ax1 = plt.axes([0.1, 0.1, 0.85, 0.8])

        ax1.fill_between(moderate.times, 0, moderate.data, color="yellow", label="Moderate (<-2)")
        ax1.fill_between(severe.times, 0, severe.data, color="orange", label="Severe (<-3)")
        ax1.fill_between(extreme.times, 0, extreme.data, color="red", label="Extreme (<-4)")

        ax1.set_ylabel("% Area", fontsize=settings.FONTSIZE)

        ax1.legend(loc="upper left", ncol=1, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, labelspacing=0.1, columnspacing=0.5)

        utils.thicken_panel_border(ax1)
        for tick in ax1.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
        for tick in ax1.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 

        minor_tick_interval = 1
        minorLocator = MultipleLocator(minor_tick_interval)
        ax1.xaxis.set_minor_locator(minorLocator)
        ax1.set_xlim([1950, int(settings.YEAR)+2])
        ax1.set_ylim([None, 44])

        # inset panel
        ax2 = plt.axes([0.7, 0.7, 0.29, 0.29])

        ax2.fill_between(moderate.times[-12:], 0, moderate.data[-12:], color="yellow", label="Moderate (<-2)")
        ax2.fill_between(severe.times[-12:], 0, severe.data[-12:], color="orange", label="Severe (<-3)")
        ax2.fill_between(extreme.times[-12:], 0, extreme.data[-12:], color="red", label="Extreme (<-4)")

        utils.thicken_panel_border(ax2)
        for tick in ax2.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE*0.7) 
        for tick in ax2.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE*0.7) 

        ax2.set_xlim([int(settings.YEAR), int(settings.YEAR)+1])
        ax2.set_xticks(np.arange(2017, 2018, 1./12.))
        ax2.set_xticklabels(["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"])
        ax2.set_ylim([None, 34])
        ax2.set_xlim([float(settings.YEAR)-0.09, float(settings.YEAR)+0.99])
        ax2.text(int(settings.YEAR)+0.1, 27, settings.YEAR)

        plt.savefig(settings.IMAGELOC+"DGT_ts{}".format(settings.OUTFMT))
        plt.close()


    #************************************************************************
    # Anomaly Map
    #************************************************************************
    # Maps for selected years

    if True:

        year = int(settings.YEAR)
        cube_list = iris.load(DATALOC + "BAMS.scPDSI.cru.4.05.{}.early.GLOBAL.annual.mean.{}.{}.nc".format(year+1, year-1, year))

        cube = cube_list[0]
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()

        bounds = [-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100]

        # select the year - plot 2020 (cube[1])

        utils.plot_smooth_map_iris(settings.IMAGELOC + "DGT_{}_rel_1901-{}".format(settings.YEAR, settings.YEAR), cube[1], settings.COLOURMAP_DICT["hydrological"], bounds, "Categories relative to 1901-{} (self-calibrating PDSI)".format(settings.YEAR), cb_extra=["Dry", "Wet"], contour=True)
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_DGT_{}_rel_1901-{}".format(settings.YEAR, settings.YEAR), cube[1], settings.COLOURMAP_DICT["hydrological"], bounds, "Categories relative to 1901-{} (self-calibrating PDSI)".format(settings.YEAR), figtext="(s) Drought (self-calibrating PDSI)", cb_extra=["Dry", "Wet"], contour=True, save_netcdf_filename="{}DGT_for_NOAA_{}.nc".format(DATALOC, dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y")))

        # plot 2020-2019

        newcube = cube[1] - cube[0]
        utils.plot_smooth_map_iris(settings.IMAGELOC + "DGT_{}-{}_rel_1901-{}".format(settings.YEAR, int(settings.YEAR)-1, settings.YEAR), newcube, settings.COLOURMAP_DICT["hydrological"], bounds, "Change in Categories from {} to {} (self-calibrating PDSI)".format(int(settings.YEAR)-1, settings.YEAR), cb_extra=["Dry", "Wet"])


    if False:
        # from 2015 report 

        # select the last year
        cube_list = []
        YEARS = [1982, 1997]
        for year in YEARS:

            loc, = np.where(cube.coord("time").points == year)
            cube_list += [cube[loc][0]]

        utils.plot_smooth_map_iris_multipanel(settings.IMAGELOC + "DGR_{}_{}_{}".format(settings.YEAR, YEARS[0], YEARS[1]), cube_list, settings.COLOURMAP_DICT["hydrological"], bounds, "Categories relative to 1901-2015 (self-calibrating PSDI)", shape=(2, 1), title=[str(y) for y in YEARS], figtext=["(a)", "(b)"])

        # select the last year
        cube_list = []
        YEARS = [1985, 1987]
        for year in YEARS:

            loc, = np.where(cube.coord("time").points == year)
            cube_list += [cube[loc][0]]

        utils.plot_smooth_map_iris_multipanel(settings.IMAGELOC + "DGR_{}_{}_{}".format(settings.YEAR, YEARS[0], YEARS[1]), cube_list, settings.COLOURMAP_DICT["hydrological"], bounds, "Categories relative to 1901-2015 (self-calibrating PSDI)", shape=(2, 1), title=[str(y) for y in YEARS], figtext=["(a)", "(b)"])

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()
#************************************************************************
#                                 END
#************************************************************************

