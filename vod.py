#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for Surface Humidity (HUM) section.
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
import os
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import calendar

import iris

import utils # RJHD utilities
import settings

DATALOC = "{}/{}/data/VOD/".format(settings.ROOTLOC, settings.YEAR)
LEGEND_LOC = 'upper left'



#
def process_cube(cube, name):
    timeUnits = cube.coord("time").units
    dt_time = timeUnits.num2date(cube.coord("time").points)
    times = np.array([d.year for d in dt_time])
    return utils.Timeseries(name, times, cube.data)


#************************************************************************
def run_all_plots():

    #*******************
    # anomaly map
    if True:
        cube = iris.load_cube(os.path.join(DATALOC, "vodca_Ku-band_bams_6.1_era5flagged_anomalyMaps_yearly.nc"), "anomalies")
        cube.coord("latitude").guess_bounds()
        cube.coord("longitude").guess_bounds()

        date_constraint = utils.periodConstraint(cube, dt.datetime(int(settings.YEAR), 1, 1), dt.datetime(int(settings.YEAR)+1, 1, 1)) 
        
        year_cube = cube.extract(date_constraint)

        bounds = [-100, -0.1, -0.05, -0.02, -0.01, 0, 0.01, 0.02, 0.05, 0.1, 100]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "VOD_anomalies", year_cube, settings.COLOURMAP_DICT["phenological"], \
                                   bounds, "Anomalies from 1991-2010", figtext="", title="VOD anomalies")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_VOD_anomalies", year_cube, settings.COLOURMAP_DICT["phenological"], \
                                   bounds, "Anomalies from 1991-2010", figtext="(ah) Vegetation Optical Depth", title="")

    #*******************
    # monthly anomaly map
    if True:
        cube = iris.load_cube(os.path.join(DATALOC, "vodca_Ku-band_bams_6.1_era5flagged_anomalyMaps_monthly.nc"), "anomalies")
        cube.coord("latitude").guess_bounds()
        cube.coord("longitude").guess_bounds()

        month_list = [cube[i] for i in range(-12, 0, 1)]

        MONTHS = [calendar.month_name[i][:3] for i in range(1, 13)]

        utils.plot_smooth_map_iris_multipanel(settings.IMAGELOC + "VOD_{}_anoms_months".format(settings.YEAR), month_list, settings.COLOURMAP_DICT["phenological"], bounds, "Anomaly (m"+r'$^{3}$'+"m"+r'$^{-3}$'+")", shape=(6, 2), title=MONTHS, figtext=["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"])


    #*******************
    # timeseries + SOI
    if True:
        cubelist = iris.load(os.path.join(DATALOC, "vodca_Ku-band_bams_6.1_era5flagged_yearAnomaliesPerHemisphere.nc"))

        for cube in cubelist:
            if cube.name() == "northern_hemisphere_coverage":
                NH_cover = process_cube(cube, "N. Hemisphere")
            elif cube.name() == "southern_hemisphere_coverage":
                SH_cover = process_cube(cube, "S. Hemisphere")
            elif cube.name() == "global_coverage":
                G_cover = process_cube(cube, "Globe")

            elif cube.name() == "northern_hemisphere_anom":
                NH = process_cube(cube, "N. Hemisphere")
            elif cube.name() == "southern_hemisphere_anom":
                SH = process_cube(cube, "S. Hemisphere")
            elif cube.name() == "global_anom":
                globe = process_cube(cube, "Globe")

            elif cube.name() == "linearmodel_north":
                NH_fit = process_cube(cube, "N. Hemisphere")
                NH_fit.ls="--"
            elif cube.name() == "linearmodel_south":
                SH_fit = process_cube(cube, "S. Hemisphere")
                SH_fit.ls = "--"
            elif cube.name() == "linearmodel_glob":
                globe_fit = process_cube(cube, "Globe")
                globe_fit.ls = "--"

            elif cube.name() == "soi":
                SOI = process_cube(cube, "SOI")

        fig = plt.figure(figsize=(8, 6))
        ax1 = plt.axes([0.15, 0.2, 0.74, 0.79])
        ax2 = ax1.twinx()

        globe.zorder=10
        NH.zorder=10
        SH.zorder=10

        # SOI - plot first so under the other lines
  
        interpTimes = np.linspace(SOI.times[0], SOI.times[-1], 1000)
        interpData = np.interp(interpTimes, SOI.times, SOI.data)
        interpSOI = utils.Timeseries("SOI", interpTimes, interpData)

        ax1.fill_between(interpSOI.times, interpSOI.data, where=interpSOI.data >= 0, \
                             color='lightskyblue', zorder=-1)
        ax1.fill_between(interpSOI.times, interpSOI.data, where=interpSOI.data <= 0, \
                             color='lightcoral', zorder=-1)

        # then plot the VOD
        utils.plot_ts_panel(ax2, [globe, NH, SH], "-", "vod", loc=LEGEND_LOC)

        for fit in [globe_fit, NH_fit, SH_fit]:
            ax2.plot(fit.times, fit.data, c=settings.COLOURS["vod"][fit.name], \
                         lw=2, ls="--", zorder=10)

        ax2.set_ylabel("VOD Anomalies", fontsize=settings.FONTSIZE)    
        ax2.set_ylim([-0.016, 0.016])


        ax1.set_xlim([1986, int(settings.YEAR)+2])
        ax1.set_ylim([-1.4, 1.4])
        ax1.set_ylabel("SOI", fontsize=settings.FONTSIZE)

        for tick in ax2.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
        for tick in ax1.yaxis.get_major_ticks():
            tick.label2.set_fontsize(settings.FONTSIZE) 

        # apply to both axes
        utils.thicken_panel_border(ax2)
        utils.thicken_panel_border(ax1)

        # and swap the labels and ticks around from standard
        ax1.yaxis.tick_right()
        ax2.yaxis.tick_left()
        ax1.yaxis.set_label_position("right")
        ax2.yaxis.set_label_position("left")

        # coverage
        ax3 = plt.axes([0.15, 0.07, 0.74, 0.13], sharex=ax1)

        utils.plot_ts_panel(ax3, [G_cover, NH_cover, SH_cover], "-", "vod", loc="")        
        ax3.set_ylim([61,74])
        ax3.set_ylabel("Valid\nObs(%)", fontsize=settings.FONTSIZE)    

        for tick in ax3.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
        for ax in [ax1, ax3]:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE) 
        for tick in ax2.yaxis.get_major_ticks():
            tick.label2.set_fontsize(settings.FONTSIZE) 

        plt.savefig(settings.IMAGELOC+"VOD_ts{}".format(settings.OUTFMT))

        plt.close()
        
    #*******************
    # Other maps
    if True:
        # ESA Biomass
        cube = iris.load(os.path.join(DATALOC, "ESACCI-BIOMASS-L4-AGB-MERGED-0d25-2017-fv1.0.nc"))[0]
        cube.coord("latitude").guess_bounds()
        cube.coord("longitude").guess_bounds()

        # mask zeros
        cube.data = np.ma.masked_where(cube.data == 0, cube.data)

        bounds = [0, 25, 50, 75, 100, 125, 150, 200, 300, 400]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "VOD_ESACCI_biomass", cube, settings.COLOURMAP_DICT["land_surface_sequential"], \
                                   bounds, "AGB (Mg/ha)", figtext="(c)", title="ESA CCI Biomass")

        
        # MODIS LAI
        cube = iris.load(os.path.join(DATALOC, "MODIS_LAI_average.nc"))[0]
        cube.coord("latitude").guess_bounds()
        cube.coord("longitude").guess_bounds()

        bounds = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "VOD_MODIS_LAI", cube, settings.COLOURMAP_DICT["land_surface_sequential"], \
                                   bounds, "", figtext="(d)", title="MODIS LAI")

        # L-band
        cube = iris.load(os.path.join(DATALOC, "vod_L_average.nc"))[0]
        cube.coord("latitude").guess_bounds()
        cube.coord("longitude").guess_bounds()

        bounds = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "VOD_Lband", cube, settings.COLOURMAP_DICT["land_surface_sequential"], \
                                   bounds, "", figtext="(a)", title="VOD-L")

        # X-band
        cube = iris.load(os.path.join(DATALOC, "vod_X_average.nc"))[0]
        cube.coord("latitude").guess_bounds()
        cube.coord("longitude").guess_bounds()

        bounds = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "VOD_Xband", cube, settings.COLOURMAP_DICT["land_surface_sequential"], \
                                   bounds, "", figtext="(b)", title="VOD-X")

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
