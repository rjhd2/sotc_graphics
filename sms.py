#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for Soil Moisture (SMS) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 23                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2018-06-05 17:55:11 +0100 (Tue, 05 Jun #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.cm as mpl_cm
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator

import iris
import calendar
import datetime as dt

import utils # RJHD utilities
import settings


data_loc = "/data/local/rdunn/SotC/{}/data/SMS/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

LEGEND_LOC = 'lower right'
LW = 3

MONTHS = [calendar.month_name[i][:3] for i in range(1,13)]

#************************************************************************
def convert_times(cube):
    '''
    Convert the netcdf times in "x from y" into decimal years

    :param cube cube: input cube

    :returns: times
    '''

    # extract the time data
    timeUnits = cube.coord("time").units
    dt_time = timeUnits.num2date(cube.coord("time").points)

    times = np.array([(date.year + (date.month - 1)/12.)  for date in dt_time ])
    
    return times # convert_times
    

#************************************************************************
def run_all_plots():


    #************************************************************************
    # Annual Map

    cube_list = iris.load(data_loc + "yearlyAnomalies_{}.nc".format(settings.YEAR))

    cube = cube_list[0]
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()

    bounds = [-100, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 100]

    utils.plot_smooth_map_iris(image_loc + "SMS_{}_esa_cci".format(settings.YEAR), cube[-1], settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1991-2014 (m"+r'$^{3}$'+"m"+r'$^{-3}$'+")")

    utils.plot_smooth_map_iris(image_loc + "p2.1_SMS_{}_esa_cci".format(settings.YEAR), cube[-1], settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1991-2014 (m"+r'$^{3}$'+"m"+r'$^{-3}$'+")", figtext = "(g) Soil Moisture")
    

    #************************************************************************
    # Seasonal Map

   
    cube_list = iris.load(data_loc + "anomalyMaps{}.nc".format(settings.YEAR))

    cube = cube_list[0]
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()

    month_list = [cube[i] for i in range(-12,0,1)]

    utils.plot_smooth_map_iris_multipanel(image_loc + "SMS_{}_anoms_seasons".format(settings.YEAR), month_list, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomaly (m"+r'$^{3}$'+"m"+r'$^{-3}$'+")", shape = (6,2), title = MONTHS, figtext = ["(a)","(b)","(c)","(d)", "(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"])


    #************************************************************************
    # Timeseries

    cube_list = iris.load(data_loc + "monthAnomaliesPerHemisphere.nc")

    south_obs = cube_list[0]
    north_obs = cube_list[1]
    north = cube_list[2]
    south = cube_list[3]
    glob = cube_list[4]
    glob_obs = cube_list[5]

    # as date in days-from - let Iris do the heavy lifting and do this one manually

    fig = plt.figure(figsize = (8,6))
    ax1 = plt.axes([0.14,0.18,0.85,0.8])

    iris.plot.plot(south, 'b', label = "S. Hemisphere", lw = LW)
    iris.plot.plot(north, "0.5", label = "N. Hemisphere", lw = LW)
    iris.plot.plot(glob, "k", label = "Global", lw = LW)
    ax1.text(0.02, 0.9, "ESA CCI SM", transform = ax1.transAxes, fontsize = settings.FONTSIZE*0.8)

    # number of observations
    ax2 = plt.axes([0.14,0.05,0.85,0.13], sharex = ax1)

    south_obs.data = 100 * south_obs.data / 244243.
    north_obs.data = 100 * north_obs.data / 244243.
    glob_obs.data = 100 * glob_obs.data / 244243.

    iris.plot.plot(south_obs, 'b', lw = LW)
    iris.plot.plot(north_obs, "0.5", lw = LW)
    iris.plot.plot(glob_obs, "k", lw = LW)


    #*******************
    # prettify
    ax1.set_ylim([-0.021,0.021])
    ax1.set_ylabel("Anomaly (m"+r'$^{3}$'+"m"+r'$^{-3}$'+")", fontsize = settings.FONTSIZE*0.8)
    ax1.axhline(0, c = '0.5', ls = '--')

    ax1.set_yticks([-0.02,-0.01,0,0.01,0.02])

    ax2.set_ylim([0,100])
    ax2.set_ylabel("% \n obs.", fontsize = settings.FONTSIZE*0.8)
    ax2.set_yticks([0,25,50,75])


    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE*0.8) 
    for ax in [ax1, ax2]:
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE*0.8) 
        ax.yaxis.set_ticks_position('left')

    ax1.legend(loc = LEGEND_LOC, ncol = 1, frameon = False, prop = {'size':settings.LEGEND_FONTSIZE}, labelspacing = 0.1, columnspacing = 0.5)
    
    utils.thicken_panel_border(ax1)
    utils.thicken_panel_border(ax2)

    minorLocator = MultipleLocator(365.242199) # in days since
    ax1.xaxis.set_minor_locator(minorLocator)

    lims = ax1.get_xlim()
    ax1.set_xlim([lims[0]-100, lims[1]+100])
#    ax2.set_xticklabels("")

    plt.savefig(image_loc+"SMS_ts_esa_cci{}".format(settings.OUTFMT))
    plt.close()

    #************************************************************************
    # Hovmuller

    cube_list = iris.load(data_loc + "hovmoeller_diagram.nc")

    for cube in cube_list:
        if cube.name() == "hovmoeller": break

    latitudes = cube.coord("latitude").points
    anoms = cube.data * 0.01 # added for offset in 2017 report

    bounds = [-100, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 100]

    # extract the time data
    timeUnits = cube.coord("time").units
    dt_time = timeUnits.num2date(cube.coord("time").points)

    times = np.array([(date.year + (date.month - 1)/12.)  for date in dt_time ])

    utils.plot_hovmuller(image_loc + "SMS_hovmuller", times, latitudes, anoms, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomaly (m"+r'$^{3}$'+"m"+r'$^{-3}$'+")")


    return # run_all_plots


#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
