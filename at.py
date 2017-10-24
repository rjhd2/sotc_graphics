#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for lower stratosphere temperature (LST) section.
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

import utils # RJHD utilities
import settings

data_loc = "/data/local/rdunn/SotC/{}/data/AT/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

LEGEND_LOC = "lower left"
LW = 3

#************************************************************************
def read_csv(filename):
    """
    Read user supplied CSV for LST into Timeseries object
    """

    indata = np.genfromtxt(filename, dtype=(float), skip_header = 1)

    years = indata[:,0] 
    months = indata[:,1]

    data = indata[:,2]
    stdev = indata[:,3]

    times = years + (months - 1.)/12.

    return utils.Timeseries("AT", times, data)


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


    at = read_csv(data_loc + "apparentTransmissionMLO_monthly.txt")

    # get 6 and 24 month smoothed curves
    at6 = make_smoothed_ts(at, 6)
    at24 = make_smoothed_ts(at, 24)

    at_mean = np.mean(at.data[at.times < 1973])


    fig, (ax1, ax2) = plt.subplots(2, figsize = (10,8), sharex=True)

    minor_tick_interval = 1
    minorLocator = MultipleLocator(minor_tick_interval)
    COLOURS = settings.COLOURS["radiation"]

    # full y - as scatter, do here
    ax1.plot(at.times, at.data, c = COLOURS[at.name], marker = ".", label = at.name, ls = "")
    ax1.axhline(at_mean, c = '0.5', ls = '--')
    ax1.xaxis.set_minor_locator(minorLocator)
    utils.thicken_panel_border(ax1)
    ax1.set_ylim([0.8,0.95])


    # zoomed y
    ax2.plot(at.times, at.data, c = COLOURS[at.name], marker = ".", label = at.name, ls = "")

    ax2.plot(at6.times, at6.data, c = "r", ls = "-")
    ax2.plot(at24.times, at24.data, c = "b", ls = "-")

    ax2.axhline(at_mean, c = '0.5', ls = '--')
    ax2.xaxis.set_minor_locator(minorLocator)
    utils.thicken_panel_border(ax2)
    ax2.set_ylim([0.9,0.949])

    # prettify

    fig.text(0.03, 0.5, "Apparent Transmission", va='center', rotation='vertical', fontsize = settings.FONTSIZE)

    plt.xlim([1954,int(settings.YEAR)+1])
    for ax in [ax1, ax2]:
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
    for tick in ax2.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)

    fig.subplots_adjust(right = 0.95, top = 0.95, bottom = 0.05, hspace = 0.001)

    plt.savefig(image_loc + "AT_ts{}".format(settings.OUTFMT))
    plt.close()

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()
#************************************************************************
#                                 END
#************************************************************************
