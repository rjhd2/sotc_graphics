#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for alpine glaciers (AGL) section.
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
# python3
from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.cm as mpl_cm
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator

import utils # RJHD utilities
import settings

DATALOC = "{}/{}/data/AGL/".format(settings.ROOTLOC, settings.YEAR)


COLOURS = settings.COLOURS["cryosphere"]
LW = 2
SCALE = 1000.

#************************************************************************
def read_glacier(filename):
    """
    Read user supplied CSV for AGL into Timeseries object
    """

    indata = np.genfromtxt(filename, delimiter=",", dtype=(int), skip_header=1, filling_values=-9999)
   

    balance = utils.Timeseries("Balance", indata[:, 0], indata[:, 3]/SCALE)
    cumul_balance = utils.Timeseries("Cumulative Balance", indata[:, 0], np.ma.masked_where(indata[:,4] == -9999, indata[:, 4]/SCALE))

    return balance, cumul_balance # read_glacier


#************************************************************************
def run_all_plots():

    #************************************************************************
    # Glaciers figure

    balance, cumul_balance = read_glacier(DATALOC + "global_mass_balance_2019.csv")

    fig = plt.figure(figsize=(8, 5.5))
    plt.clf()
    ax1 = plt.axes([0.165, 0.09, 0.7, 0.88])

    ax1.bar(balance.times[:-1], balance.data[:-1], color=COLOURS[balance.name], label=balance.name, align="center")
    ax1.bar(balance.times[-1], balance.data[-1], color=COLOURS[balance.name], label=balance.name, align="center", alpha=0.5)

    ax2 = ax1.twinx()

    ax2.plot(cumul_balance.times, cumul_balance.data, marker="o", ms=5, ls='-', c=COLOURS[cumul_balance.name], \
                 label=cumul_balance.name, lw=LW)

    ax1.axhline(0, c='0.5', ls='--')


    ax1.set_xlim([1980, int(settings.YEAR)+2])
    ax1.set_ylim([-2, 0.2])
    ax1.set_ylabel("Mean Specific Annual\nBalance (1000mm w.e.)", fontsize=settings.FONTSIZE, color="r")
    ax2.set_ylim([-22, 2])
    ax2.set_ylabel("Cumulative Mean Specific\nAnnual Balance (1000mm w.e.)", fontsize=settings.FONTSIZE)
    minorLocator = MultipleLocator(1)
    ax1.xaxis.set_minor_locator(minorLocator)

    utils.thicken_panel_border(ax1)
    utils.thicken_panel_border(ax2)
    ax1.yaxis.set_tick_params(right=False)
    # red tickmarks
    ax1.tick_params(axis='y', colors='red')

    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)
    for tick in ax2.yaxis.get_major_ticks():
        tick.label2.set_fontsize(settings.FONTSIZE)

    plt.savefig(settings.IMAGELOC+"AGL_ts{}".format(settings.OUTFMT))

    plt.close()

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 End
#************************************************************************
