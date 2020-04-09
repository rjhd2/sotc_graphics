#!/usr/bin/env python
# python3
from __future__ import absolute_import
from __future__ import print_function
#************************************************************************
#
#  Plot figures and output numbers for Cloudiness (CLD) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 26                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2019-04-17 15:34:18 +0100 (Wed, 17 Apr #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import os
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator

import scipy.io

import utils # RJHD utilities
import settings


DATALOC = "{}/{}/data/CLD/".format(settings.ROOTLOC, settings.YEAR)

LEGEND_LOC = 'lower left'

CLIMSTART = 2003
CLIMEND = 2015

#************************************************************************
def read_ts(filename, anomaly=False, fullbase=False):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read
    :param bool anomaly: run anomalies
    :param bool fullbase: run anomalies using full baseperiod

    :returns: Timeseries object s
    '''

    raw_data = np.genfromtxt(filename, dtype=(float), skip_header=1)
    
    raw_data = np.ma.masked_where(raw_data <= -999.0, raw_data)    

    times = raw_data[:, 0]

    cstart, = np.where(times == CLIMSTART)
    cend, = np.where(times == CLIMEND)

    if fullbase:
        clims = np.ma.mean(raw_data, axis=0)
    else:
        clims = np.ma.mean(raw_data[cstart[0]:cend[0]], axis=0)

    if anomaly:
        data = raw_data - clims
    else:
        data = raw_data

    patmosx = utils.Timeseries("PATMOS-x/AVHRR", times, data[:, 1])
    hirs = utils.Timeseries("HIRS", times, data[:, 2])
    misr = utils.Timeseries("MISR", times, data[:, 3])
    modis = utils.Timeseries("AQUA MODIS C6", times, data[:, 4])
    calipso = utils.Timeseries("CALIPSO", times, data[:, 5])
    ceres = utils.Timeseries("CERES", times, data[:, 6])
    satcorps = utils.Timeseries("SatCORPS", times, data[:, 7])
    clara_a2 = utils.Timeseries("CLARA-A2", times, data[:, 8])
    patmosdx = utils.Timeseries("PATMOS-x/AQUA MODIS", times, data[:, 9])
    cci = utils.Timeseries("Cloud CCI AVHRR-PMv3", times, data[:, 10])

    return patmosx, hirs, misr, modis, calipso, ceres, satcorps, clara_a2, patmosdx, cci # read_ts


#************************************************************************
def run_all_plots():

    #************************************************************************
    # Cloudiness timeseries
    if True:
        plt.clf()
        fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 6.5), sharex=True)

        infilename = os.path.join(DATALOC, "{}_global_cloudiness_timeseries_v2.txt".format(settings.YEAR))

        # anomalies
        patmosx, hirs, misr, modis, calipso, ceres, satcorps, clara_a2, patmosdx, cci = \
            read_ts(infilename, anomaly=True)

        utils.plot_ts_panel(ax1, [patmosx, hirs, misr, modis, calipso, ceres, satcorps, clara_a2, patmosdx, cci], "-", "hydrological", loc="")

        ax1.text(0.02, 0.9, "(a) Satellite - Anomalies", transform=ax1.transAxes, fontsize=settings.FONTSIZE)

        # actuals
        patmosx, hirs, misr, modis, calipso, ceres, satcorps, clara_a2, patmosdx, cci = \
            read_ts(infilename)

        utils.plot_ts_panel(ax2, [patmosx, hirs, misr, modis, calipso, ceres, satcorps, clara_a2, patmosdx, cci], "-", "hydrological", loc=LEGEND_LOC, ncol=3)

        ax2.text(0.02, 0.9, "(b) Satellite - Actual", transform=ax2.transAxes, fontsize=settings.FONTSIZE)

        #*******************
        # prettify
        ax1.set_ylabel("Anomaly (%)", fontsize=settings.FONTSIZE)
        ax2.set_ylabel("(%)", fontsize=settings.FONTSIZE)

        for tick in ax2.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 

        minorLocator = MultipleLocator(1)
        for ax in [ax1, ax2]:
            utils.thicken_panel_border(ax)
            ax.set_yticks(ax.get_yticks()[1:])
            ax.xaxis.set_minor_locator(minorLocator)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE) 

        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        fig.subplots_adjust(right=0.96, top=0.99, bottom=0.04, left=0.09, hspace=0.001)

        plt.xlim([hirs.times[0]-1, hirs.times[-1]+1])
        ax1.set_ylim([-4.4, 6.9])
        ax2.set_ylim([0, 95])

        plt.savefig(settings.IMAGELOC+"CLD_ts{}".format(settings.OUTFMT))
        plt.close()

    #************************************************************************
    # make a version using the full base period of the data and only the pre2000 datasets
    if False:
        plt.clf()
        fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 6.5), sharex=True)

        # anomalies
        patmosx, hirs, misr, modis, calipso, ceres, satcorps, clara_a2, patmosdx, cci = \
            read_ts(DATALOC + "{}_global_cloudiness_timeseries.txt".format(settings.YEAR), \
                        anomaly=True, fullbase=True)

        utils.plot_ts_panel(ax1, [patmosx, hirs, satcorps, clara_a2, cci], "-", "hydrological", loc="")

        ax1.text(0.02, 0.9, "(a) Satellite - Anomalies", transform=ax1.transAxes, fontsize=settings.FONTSIZE)

        # actuals
        patmosx, hirs, misr, modis, calipso, ceres, satcorps, clara_a2, patmosdx, cci = \
            read_ts(DATALOC + "{}_global_cloudiness_timeseries.txt".format(settings.YEAR))

        utils.plot_ts_panel(ax2, [patmosx, hirs, satcorps, clara_a2, cci], "-", "hydrological", \
                                loc=LEGEND_LOC, ncol=3)

        ax2.text(0.02, 0.9, "(b) Satellite - Actual", transform=ax2.transAxes, fontsize=settings.FONTSIZE)

        #*******************
        # prettify
        ax1.set_ylabel("Anomaly (%)", fontsize=settings.FONTSIZE)
        ax2.set_ylabel("(%)", fontsize=settings.FONTSIZE)

        for tick in ax2.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 

        minorLocator = MultipleLocator(1)
        for ax in [ax1, ax2]:
            utils.thicken_panel_border(ax)
            ax.set_yticks(ax.get_yticks()[1:])
            ax.xaxis.set_minor_locator(minorLocator)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE) 

        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        fig.subplots_adjust(right=0.95, top=0.95, hspace=0.001)

        plt.xlim([hirs.times[0]-1, hirs.times[-1]+1])
        ax1.set_ylim([-4.9, 4.4])
        ax2.set_ylim([0, 95])


        plt.savefig(settings.IMAGELOC+"CLD_ts_fullbaseperiod{}".format(settings.OUTFMT))
        plt.close()

    #************************************************************************
    # Cloudiness map
    if True:
        mapfile_dict = scipy.io.readsav(DATALOC + "patmosx_global_monthly_cloudiness_anomaly_map_{}.sav".format(settings.YEAR))

        annual_anoms = mapfile_dict["annual_anom"]*100.
        djf_anoms = mapfile_dict["djf_anom"]*100.
        jja_anoms = mapfile_dict["jja_anom"]*100.
        mam_anoms = mapfile_dict["mam_anom"]*100.
        son_anoms = mapfile_dict["son_anom"]*100.

        lats = mapfile_dict["lat"]
        lons = mapfile_dict["lon"]

        cube = utils.make_iris_cube_2d(annual_anoms, lats[:, 0], lons[0], "CLD_anom", "%")

        bounds = [-100, -15, -10, -5, -2.5, 0, 2.5, 5, 10, 15, 100]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "CLD_{}_anoms".format(settings.YEAR), cube, \
                                       settings.COLOURMAP_DICT["hydrological"], bounds, \
                                       "Anomalies from 1981-2010 (%)")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_CLD_{}_anoms".format(settings.YEAR), \
                                       cube, settings.COLOURMAP_DICT["hydrological"], bounds, \
                                       "Anomalies from 1981-2010 (%)", figtext="(n) Cloudiness")


    #************************************************************************
    # Cloudiness Seasonal Map
    if True:
        cubelist = []
        for season in [djf_anoms, mam_anoms, jja_anoms, son_anoms]:

            cube = utils.make_iris_cube_2d(season, lats[:, 0], lons[0], "CLD_anom", "%")
            cubelist += [cube]


        utils.plot_smooth_map_iris_multipanel(settings.IMAGELOC + "CLD_{}_anoms_seasons".format(settings.YEAR), \
                                                  cubelist, settings.COLOURMAP_DICT["hydrological"], \
                                                  bounds, "Anomaly (%)", shape=(2, 2), \
                                                  title=["DJF", "MAM", "JJA", "SON"], \
                                                  figtext=["(a)", "(b)", "(c)", "(d)"])


    #************************************************************************
    # Cloudiness Hovmoller
    if True:
        data_dict = scipy.io.readsav(DATALOC + "patmosx_global_monthly_cloudiness_hovmuller_{}.sav".format(settings.YEAR))

        lats = data_dict["latitude"]
        times = data_dict["hov_time"]
        anoms = data_dict["hov_anom"]*100.

        utils.plot_hovmuller(settings.IMAGELOC + "CLD_hovmuller", times, lats, anoms, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomaly (%)")

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
