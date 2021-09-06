#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for Soil Moisture (SMS) section.
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
import calendar
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator

import iris

import utils # RJHD utilities
import settings


DATALOC = "{}/{}/data/SMS/".format(settings.ROOTLOC, settings.YEAR)

LEGEND_LOC = 'lower right'
LW = 2

MONTHS = [calendar.month_name[i][:3] for i in range(1, 13)]

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

    times = np.array([(date.year + (date.month - 1)/12.)  for date in dt_time])
    
    return times # convert_times
    

#************************************************************************
def run_all_plots():


    #************************************************************************
    # Timeseries
    if True:

        cube_list = np.array(iris.load(DATALOC + "monthAnomaliesPerHemisphere.nc"))

        names = np.array([c.var_name for c in cube_list])

        north_obs = cube_list[names == "nObs_north"][0]
        south_obs = cube_list[names == "nObs_south"][0]
        glob_obs = cube_list[names == "nObs_global"][0]
        north = cube_list[names == "Anomalies_north"][0]
        south = cube_list[names == "Anomalies_south"][0]
        glob = cube_list[names == "Anomalies_global"][0]

        # as date in days-from - let Iris do the heavy lifting and do this one manually

        fig = plt.figure(figsize=(8, 6))
        ax1 = plt.axes([0.14, 0.2, 0.84, 0.79])

        iris.plot.plot(south, 'r', label="S. Hemisphere", lw=LW)
        iris.plot.plot(north, "b", label="N. Hemisphere", lw=LW)
        iris.plot.plot(glob, "k", label="Global", lw=3)
        ax1.text(0.02, 0.9, "ESA CCI SM", transform=ax1.transAxes, fontsize=settings.FONTSIZE)

        # number of observations
        ax2 = plt.axes([0.14, 0.07, 0.84, 0.13], sharex=ax1)

        south_obs.data = 100 * south_obs.data / 244243.
        north_obs.data = 100 * north_obs.data / 244243.
        glob_obs.data = 100 * glob_obs.data / 244243.

        iris.plot.plot(south_obs, 'r', lw=LW)
        iris.plot.plot(north_obs, "b", lw=LW)
        iris.plot.plot(glob_obs, "k", lw=3)


        #*******************
        # prettify
        ax1.set_ylim([-0.021, 0.021])
        ax1.set_ylabel("Anomaly (m"+r'$^{3}$'+"m"+r'$^{-3}$'+")", fontsize=settings.FONTSIZE)
        ax1.axhline(0, c='0.5', ls='--')

        ax1.set_yticks([-0.02, -0.01, 0, 0.01, 0.02])

        ax2.set_ylim([0, 100])
        ax2.set_ylabel("% \n obs.", fontsize=settings.FONTSIZE)
        ax2.set_yticks([0, 25, 50, 75])


        for tick in ax2.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
        for ax in [ax1, ax2]:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE) 
            ax.yaxis.set_ticks_position('left')


        ax1.legend(loc=LEGEND_LOC, ncol=1, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, labelspacing=0.1, columnspacing=0.5)

        utils.thicken_panel_border(ax1)
        utils.thicken_panel_border(ax2)

        minorLocator = MultipleLocator(365.242199) # in days since
        ax1.xaxis.set_minor_locator(minorLocator)

        lims = ax1.get_xlim()
        ax1.set_xlim([lims[0]-100, lims[1]+100])
    #    ax2.set_xticklabels("")

        plt.savefig(settings.IMAGELOC+"SMS_ts_esa_cci{}".format(settings.OUTFMT))
        plt.close()

    #************************************************************************
    # Hovmuller
    if True:

        cube_list = iris.load(DATALOC + "hovmoeller_diagram.nc")

        for cube in cube_list:
            if cube.name() == "sm": break

        data = cube.data[:]
        data.mask = np.zeros(data.shape)
        data.fill_value = -9999

        data = np.ma.masked_where(data == data.fill_value, data)
        cube.data = data

        latitudes = cube.coord("latitude").points
        anoms = cube.data

        bounds = [-100, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 100]

        # extract the time data
        timeUnits = cube.coord("time").units
        dt_time = timeUnits.num2date(cube.coord("time").points)

        times = np.array([(date.year + (date.month - 1)/12.)  for date in dt_time])

        utils.plot_hovmuller(settings.IMAGELOC + "SMS_hovmuller", times, latitudes, anoms, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomaly (m"+r'$^{3}$'+"m"+r'$^{-3}$'+")")

    #************************************************************************
    # Annual Map
    if True:
        cube_list = iris.load(DATALOC + "anomalyMaps_yearly.nc")

        cube = cube_list[0]
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()

        bounds = [-100, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 100]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "SMS_{}_esa_cci".format(settings.YEAR), cube[-1], settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1991-2010 (m"+r'$^{3}$'+"m"+r'$^{-3}$'+")")

        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_SMS_{}_esa_cci".format(settings.YEAR), cube[-1], settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1991-2010 (m"+r'$^{3}$'+"m"+r'$^{-3}$'+")", figtext="(r) Soil Moisture")

        # 2020 only special request
        utils.plot_smooth_map_iris(settings.IMAGELOC + "SMS_{}_esa_cci_for_RvdS".format(settings.YEAR), cube[-1], settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies for 2020 (m"+r'$^{3}$'+"m"+r'$^{-3}$'+")")
    #************************************************************************
    # Seasonal Map
    if True:
  
        cube_list = iris.load(DATALOC + "anomalyMaps_monthly.nc")

        cube = cube_list[0]
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()

        month_list = [cube[i] for i in range(-12, 0, 1)]

        utils.plot_smooth_map_iris_multipanel(settings.IMAGELOC + "SMS_{}_anoms_seasons".format(settings.YEAR), month_list, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomaly (m"+r'$^{3}$'+"m"+r'$^{-3}$'+")", shape=(6, 2), title=MONTHS, figtext=["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"])


    return # run_all_plots


#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
