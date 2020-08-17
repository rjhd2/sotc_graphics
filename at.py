#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for lower stratosphere temperature (LST) section.
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

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import utils # RJHD utilities
import settings

DATALOC = "{}/{}/data/AT/".format(settings.ROOTLOC, settings.YEAR)

LEGEND_LOC = "lower left"
LW = 3

#************************************************************************
def read_csv(filename):
    """
    Read user supplied CSV for LST into Timeseries object
    """

    indata = np.genfromtxt(filename, delimiter=",", dtype=(float), skip_header = 1, encoding="latin-1")

    times = indata[:,0] 

    data = indata[:,1]
    fit = indata[:,2]

    return utils.Timeseries("AT", times, data), utils.Timeseries("AT", times, fit)


#************************************************************************
def make_smoothed_ts(ts, points):
    '''
    Make a smoothed timeseries array.  

    :param Timeseries ts: input timeseries
    :param int points: for boxcar
    '''

    smoothed = np.ma.zeros(len(ts.data))
    smoothed.mask = np.zeros(smoothed.shape)
    
    for d in range(len(ts.data)):

        # if a running centred mean:
        if d < points/2.:
            smoothed.mask[d] = 1
        elif d > (len(ts.data) - points/2.):
            smoothed.mask[d] = 1
        else:
            smoothed[d] = np.mean(ts.data[d - points/2: d +  points/2])

        # if a running previous N-month mean
        # if d < points:
        #     smoothed[d] = np.mean(ts.data[:d])
        # elif d > (len(ts.data) - points):
        #     smoothed.mask[d] = 1
        # else:
        #     smoothed[d] = np.mean(ts.data[d - points: d])
    
    smoothed_ts = utils.Timeseries(ts.name, ts.times, smoothed) 

    return smoothed_ts # make_masked_smoothed_ts


#************************************************************************
def run_all_plots():
    #************************************************************************
    # Timeseries figures


    at, at24 = read_csv(DATALOC + "mlo_trans_csv.txt")

    # get 6 and 24 month smoothed curves
    # at6 = make_smoothed_ts(at, 6)
    # at24 = make_smoothed_ts(at, 24)

    at_mean = np.mean(at.data[at.times < 1962])


    # fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 6.5), sharex=True)
    fig = plt.figure(figsize=(8, 5))
    ax1 = plt.axes([0.15, 0.1, 0.8, 0.85])

    minor_tick_interval = 1
    minorLocator = MultipleLocator(minor_tick_interval)
    COLOURS = settings.COLOURS["radiation"]

    # full y - as scatter, do here
    ax1.plot(at.times, at.data, c = COLOURS[at.name], marker = ".", label = at.name, ls = "")
    ax1.axhline(at_mean, c = '0.5', ls = '--')
    ax1.xaxis.set_minor_locator(minorLocator)
    utils.thicken_panel_border(ax1)
    ax1.set_ylim([0.8,0.95])
    ax1.plot(at24.times, at24.data, c = "0.3", ls = "-")

    # # zoomed y
    # ax2.plot(at.times, at.data, c = COLOURS[at.name], marker = ".", label = at.name, ls = "")

    # ax2.plot(at6.times, at6.data, c = "r", ls = "-")
    # ax2.plot(at24.times, at24.data, c = "b", ls = "-")

    # ax2.axhline(at_mean, c = '0.5', ls = '--')
    # ax2.xaxis.set_minor_locator(minorLocator)
    # utils.thicken_panel_border(ax2)
    # ax2.set_ylim([0.9,0.949])

    # prettify

    fig.text(0.03, 0.5, "Apparent Transmission", va='center', rotation='vertical', fontsize = settings.FONTSIZE)

    ax1.text(1961, 0.9, "Agung", fontsize=settings.FONTSIZE)
    ax1.text(1970, 0.88, "El Chichon", fontsize=settings.FONTSIZE)
    ax1.text(1993, 0.88, "Pinatubo", fontsize=settings.FONTSIZE)

    ax1.text(2010, 0.9, "Nabro", fontsize=settings.FONTSIZE)
    ax1.plot([2012.5, 2012], [0.908, 0.92], "k-")

    ax1.fill_between([2002, 2012], [0.7, 0.7], [1.0, 1.0], color="0.7")
    ax1.fill_between([2017, 2020], [0.7, 0.7], [1.0, 1.0], color="0.7")

    plt.xlim([1954,int(settings.YEAR)+2])
    for tick in ax1.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
    for tick in ax1.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)

    # for tick in ax2.yaxis.get_major_ticks():
    #         tick.label.set_fontsize(settings.FONTSIZE)
    # for tick in ax2.xaxis.get_major_ticks():
    #     tick.label.set_fontsize(settings.FONTSIZE)

    fig.subplots_adjust(right=0.96, top=0.98, bottom=0.04, hspace=0.001)

    plt.savefig(settings.IMAGELOC + "AT_ts{}".format(settings.OUTFMT))
    plt.close()

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()
#************************************************************************
#                                 END
#************************************************************************
