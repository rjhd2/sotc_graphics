#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for snow cover (SNW) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 22                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2018-04-06 15:34:21 +0100 (Fri, 06 Apr #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.cm as mpl_cm
import matplotlib as mpl

import struct

import utils # RJHD utilities
import settings


data_loc = "/data/local/rdunn/SotC/{}/data/SNW/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

LEGEND_LOC = 'upper left'

SCALE = 1.e6


#************************************************************************
def read_snow(filename):
    """
    Read user supplied fixed field file into array
    """


    # fieldwidths = (8,8,8,8,8,8,8,8)
    # fmtstring = ''.join('%ds' % f for f in fieldwidths)

    # parse = struct.Struct(fmtstring).unpack_from

    # indata = []
    # with open(filename, 'r') as infile:
    #     for ll, line in enumerate(infile):

    #         if ll != 0: # skip the first
    #             fields = parse(line)
                
    #             indata += [[int(f) for f in fields]]
            
    indata = np.genfromtxt(filename, skip_header = 1, delimiter = ",")

    indata = np.ma.array(indata)
    indata = np.ma.masked_where(indata == 0, indata)

    # process the years and months to give decimals
    years = indata[:,0]
    months = indata[:,1]

    times = years + (months - 1)/12.

    # make timeseries objects
    NH = utils.Timeseries("N Hemisphere", times, indata[:,3]/SCALE)
    Eurasia = utils.Timeseries("Eurasia", times, indata[:,5]/SCALE)
    NAmer = utils.Timeseries("N America", times, indata[:,7]/SCALE)
   
    return NH, Eurasia, NAmer # read_snow

#************************************************************************
def run_all_plots():

    #************************************************************************
    # Snow cover figure

    NH, Eurasia, NAmer = read_snow(data_loc + "rutgers-sce-anom12-31-{}_rd.csv".format(settings.YEAR))


    fig = plt.figure(figsize = (10,6))
    plt.clf()
    ax = plt.axes([0.10, 0.10, 0.87, 0.87])

    utils.plot_ts_panel(ax, [NH, Eurasia, NAmer], "-", "cryosphere", loc = LEGEND_LOC)

    # sort formatting
    plt.xlim([1966,int(settings.YEAR)+1])
    plt.ylim([-1.9,3.3])
    ax.set_ylabel("Anomaly (Million km"+r'$^2$'+")", fontsize = settings.FONTSIZE)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 

    plt.savefig(image_loc+"SNW_ts{}".format(settings.OUTFMT))

    plt.close()

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()


#************************************************************************
#                                 End
#************************************************************************
