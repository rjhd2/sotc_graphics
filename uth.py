#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for Upper Tropospheric Humidity (UTH) section.
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
from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

import utils # RJHD utilities
import settings


DATALOC = "{}/{}/data/UTH/".format(settings.ROOTLOC, settings.YEAR)

DECIMAL_MONTHS = np.arange(12)/12.

LEGEND_LOC = 'upper left'

#************************************************************************
def read_ts(filename, start, name, smooth=False):
    '''
    Read the timeseries data, and use the hard-coded start and end years
    to make up the time axis for returning Timeseries objects.

    :param str filename: file to read
    :param int start: year of start of data 
    :param str name: name of Timeseries object
    :param bool smooth: if True, then == value of smoothing

    :returns: Timeseries object
    '''

    data = np.genfromtxt(filename, dtype=(float), skip_header=5)

    # create the times
    times = []
    year = start
    while year <= int(settings.YEAR):

        times += [list(year + DECIMAL_MONTHS)]

        year += 1

    times = np.array([val for sublist in times for val in sublist])

    assert len(times) == len(data)

    # smooth if required
    if smooth:
        data = utils.boxcar(data, smooth)


    ts = utils.Timeseries(name, times, data)

    return ts # read_ts

#************************************************************************
def read_map(data_loc, name):
    '''
    Read the map data

    :param str data_loc: location where to find files
    :param str name: name of measuring platform

    :returns: cube
    '''

    lons = np.genfromtxt(data_loc + "{}_lon_map.aa".format(name), dtype=(float), skip_header=5)
    lats = np.genfromtxt(data_loc + "{}_lat_map.aa".format(name), dtype=(float), skip_header=5)

    data = np.zeros((len(lats), len(lons)))

    this_lat = []
    lon_ctr = 0
    lat_ctr = 0

    with open(data_loc + "{}_{}_anom_map.aa".format(name, settings.YEAR), 'r') as infile:
        for ll, line in enumerate(infile):
            if ll > 4: # skip the first 5 lines
                for ls in line.split():
                    this_lat.append(ls)

                if len(this_lat) == len(lons):
                    # read in all for one longitude
                    data[lat_ctr, :] = this_lat
                    lat_ctr += 1

                    this_lat = []
                    
    cube = utils.make_iris_cube_2d(data, lats, lons, "UTH_anom", "%")

    return cube # read_map

#************************************************************************
def run_all_plots():



    #************************************************************************
    # Upper Tropospheric Humidity timeseries figure
    if True:
        HIRSSTART = 1979
        MWSTART = 1999
        ERASTART = 1979

        # smooth by 3 months
        hirs = read_ts(DATALOC + "hirs_data.aa", HIRSSTART, "HIRS", smooth=3)
        mw = read_ts(DATALOC + "mw_data.aa", MWSTART, "Microwave", smooth=3)
        era5 = read_ts(DATALOC + "era5_data.aa", ERASTART, "ERA5", smooth=3)


        fig = plt.figure(figsize=(8, 5))
        ax = plt.axes([0.11, 0.08, 0.86, 0.90])

        utils.plot_ts_panel(ax, [hirs, mw, era5], "-", "hydrological", loc=LEGEND_LOC, ncol=3)

        #*******************
        # prettify

        fig.text(0.01, 0.5, "Anomalies (% rh)", va='center', rotation='vertical', fontsize=settings.FONTSIZE)

        plt.ylim([-2.0, 2.0])
        plt.xlim([1979, int(settings.YEAR)+1.5])

        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)

        plt.savefig(settings.IMAGELOC+"UTH_ts{}".format(settings.OUTFMT))
        plt.close()

    #************************************************************************
    # HIRS map
    if True:
        cube = read_map(DATALOC, "hirs")

        bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]

    #    utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_UTH_{}_anoms_hirs".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 2001-2010 (% rh)", figtext = "(n) Upper Tropospheric Humidity")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "UTH_{}_anoms_hirs".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 2001-2010 (% rh)")


    #************************************************************************
    # MW map
    if True:
        cube = read_map(DATALOC, "mw")

        bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "UTH_{}_anoms_mw".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 2001-2010 (% rh)")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_UTH_{}_anoms_mw".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 2001-2010 (% rh)", figtext="(j) Upper Tropospheric Humidity")

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
