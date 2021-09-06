#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for NMAT section.
#       For BAMS SotC 2020
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


DATALOC = "{}/{}/data/NMT/".format(settings.ROOTLOC, settings.YEAR)
LEGEND_LOC = 'lower right'

#************************************************************************
def read_csv_1(filename, dataset, region):

    all_data = np.genfromtxt(filename, delimiter=",", skip_header=1, dtype=(str))

    locs, = np.where(np.logical_and(all_data[:, 1] == dataset, all_data[:, 3] == region))

    times = all_data[locs, 0].astype(float)
    data = all_data[locs, 2]

    locs, = np.where(data == "NA")
    data[locs] = -99.9
    data = np.ma.array(data).astype(float)
    data.mask = np.zeros(data.shape)
    data.mask[locs] = True

    return utils.Timeseries(dataset, times, data) # read_csv_1

#************************************************************************
def read_csv_2(filename, dataset, region):

    all_data = np.genfromtxt(filename, delimiter=",", skip_header=1, dtype=(str))

    locs, = np.where(np.logical_and(all_data[:, 2] == dataset, all_data[:, 1] == region))

    times = all_data[locs, 0].astype(float)
    data = all_data[locs, 3].astype(float)

    return utils.Timeseries(dataset, times, data) # read_csv_2


#************************************************************************
def run_all_plots():

    #************************************************************************
    # Timeseries 1 figure

    if True:
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(8, 10), sharex=True)

        filename =  "fig_1_data_v2.csv"
        filename =  "fig_1_data.csv"
        # Globe
        classnmat = read_csv_1(DATALOC + filename, "CLASSnmat v2", "Global")
        uahnmat = read_csv_1(DATALOC + filename, "UAHNMAT v1", "Global")
        hadsst = read_csv_1(DATALOC + filename, "HadSST4", "Global")
        utils.plot_ts_panel(ax1, [classnmat, uahnmat, hadsst], "-", "temperature", loc=LEGEND_LOC)
        ax1.text(0.02, 0.88, "(a) Global", transform=ax1.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)

        # Northern Extra-Tropics
        classnmat = read_csv_1(DATALOC + filename, "CLASSnmat v2", "Northern Extra-Tropics")
        uahnmat = read_csv_1(DATALOC + filename, "UAHNMAT v1", "Northern Extra-Tropics")
        hadsst = read_csv_1(DATALOC + filename, "HadSST4", "Northern Extra-Tropics")
        utils.plot_ts_panel(ax2, [classnmat, uahnmat, hadsst], "-", "temperature", loc="")
        ax2.text(0.02, 0.88, "(b) Northern Extra-Tropics", transform=ax2.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)

        # Tropics
        classnmat = read_csv_1(DATALOC + filename, "CLASSnmat v2", "Tropics")
        uahnmat = read_csv_1(DATALOC + filename, "UAHNMAT v1", "Tropics")
        hadsst = read_csv_1(DATALOC + filename, "HadSST4", "Tropics")
        utils.plot_ts_panel(ax3, [classnmat, uahnmat, hadsst], "-", "temperature", loc="")
        ax3.text(0.02, 0.88, "(c) Tropics", transform=ax3.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)

        # Southern Extra-Tropics
        classnmat = read_csv_1(DATALOC + filename, "CLASSnmat v2", "Southern Extra-Tropics")
        uahnmat = read_csv_1(DATALOC + filename, "UAHNMAT v1", "Southern Extra-Tropics")
        hadsst = read_csv_1(DATALOC + filename, "HadSST4", "Southern Extra-Tropics")
        utils.plot_ts_panel(ax4, [classnmat, uahnmat, hadsst], "-", "temperature", loc="")
        ax4.text(0.02, 0.88, "(d) Southern Extra-Tropics", transform=ax4.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)


        plt.xlim([1897, int(settings.YEAR)+3])

        for tick in ax4.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        for ax in [ax1, ax2, ax3, ax4]:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)
            ax.xaxis.set_minor_locator(MultipleLocator(5))
        
            ax.set_ylim([-1.3, 1.3])

        fig.text(0.01, 0.55, "Anomaly ("+r'$^\circ$'+"C)", fontsize=settings.FONTSIZE, rotation="vertical")
        fig.subplots_adjust(right=0.98, top=0.98, bottom=0.04, hspace=0.001)

        plt.savefig(settings.IMAGELOC+"NMT_ts{}".format(settings.OUTFMT))

        plt.close()

    #************************************************************************
    # Timeseries 1 figure
    if True:
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(8, 10), sharex=True)

        # Globe
        sst = read_csv_2(DATALOC + "fig_2_data.csv", "SST", "Global")
        at = read_csv_2(DATALOC + "fig_2_data.csv", "Air Temperature", "Global")
        utils.plot_ts_panel(ax1, [sst, at], "-", "temperature", loc=LEGEND_LOC)
        ax1.text(0.02, 0.88, "(a) Global", transform=ax1.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)

        # Northern Extra-Tropics
        sst = read_csv_2(DATALOC + "fig_2_data.csv", "SST", "Northern Extra-Tropics")
        at = read_csv_2(DATALOC + "fig_2_data.csv", "Air Temperature", "Northern Extra-Tropics")
        utils.plot_ts_panel(ax2, [sst, at], "-", "temperature", loc="")
        ax2.text(0.02, 0.88, "(b) Northern Extra-Tropics", transform=ax2.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)

        # Tropics
        sst = read_csv_2(DATALOC + "fig_2_data.csv", "SST", "Tropics")
        at = read_csv_2(DATALOC + "fig_2_data.csv", "Air Temperature", "Tropics")
        utils.plot_ts_panel(ax3, [sst, at], "-", "temperature", loc="")
        ax3.text(0.02, 0.88, "(c) Tropics", transform=ax3.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)

        # Southern Extra-Tropics
        sst = read_csv_2(DATALOC + "fig_2_data.csv", "SST", "Southern Extra-Tropics")
        at = read_csv_2(DATALOC + "fig_2_data.csv", "Air Temperature", "Southern Extra-Tropics")
        utils.plot_ts_panel(ax4, [sst, at], "-", "temperature", loc="")
        ax4.text(0.02, 0.88, "(d) Southern Extra-Tropics", transform=ax4.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)

        plt.xlim([1948, int(settings.YEAR)+2])

        for tick in ax4.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        for ax in [ax1, ax2, ax3, ax4]:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)
            ax.xaxis.set_minor_locator(MultipleLocator(5))
        
            ax.set_ylim([-0.7, 0.9])

        fig.text(0.01, 0.55, "Anomaly ("+r'$^\circ$'+"C)", fontsize=settings.FONTSIZE, rotation="vertical")
        fig.subplots_adjust(right=0.98, top=0.98, bottom=0.04, hspace=0.001)

        plt.savefig(settings.IMAGELOC+"NMT_ts_sst{}".format(settings.OUTFMT))

        plt.close()

    #************************************************************************
    # Anomaly figure

    if True:
        # Read in ERA anomalies

        cube_list = iris.load(os.path.join(DATALOC, "CLASSnmat_2.0.0.0_anomaly_1981_2010_ANNMEAN_{}.nc".format(settings.YEAR)))
        for cube in cube_list:
            if cube.var_name == "t10m_anomaly":
                break

        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()

        bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "NMT_CLASS_anomalies_{}".format(settings.YEAR), cube[0], settings.COLOURMAP_DICT["temperature"], bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", title="CLASSNMAT")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_NMT_CLASS_anomalies_{}".format(settings.YEAR), cube[0], settings.COLOURMAP_DICT["temperature"], bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", figtext="(f) Night Marine Air Temperature")

    #************************************************************************
    # Anomaly figure

    if True:
        # Read in ERA anomalies

        cube_list = iris.load(os.path.join(DATALOC, "fig_3_data.nc".format(settings.YEAR)))
        cube = cube_list[0]

        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()

        bounds = [-100, -0.5, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.5, 100]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "NMT_trend_diff_{}".format(settings.YEAR), cube, settings.COLOURMAP_DICT["temperature"], bounds, "Trend ("+r'$^{\circ}$'+"C decade "+r'$^{-1}$'+")", title="MAT-SST trend")

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()


#************************************************************************
#                                 END
#************************************************************************
