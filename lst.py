#!/usr/bin/env python
# python3
from __future__ import absolute_import
from __future__ import print_function
##************************************************************************
#
#  Plot figures and output numbers for lower stratosphere temperature (LST) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 29                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2020-08-05 12:12:39 +0100 (Wed, 05 Aug #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************

import datetime as dt
import calendar
import numpy as np
import matplotlib.pyplot as plt

import iris

import utils # RJHD utilities
import settings

DATALOC = "{}/{}/data/LST/".format(settings.ROOTLOC, settings.YEAR)

CLIM_PERIOD = "8110"
LEGEND_LOC = 'lower left'

MONTHS = [calendar.month_name[i][:3] for i in range(1, 13)]

#************************************************************************
def read_csv(filename):
    """
    Read user supplied CSV for LST into Timeseries object
    """

    indata = np.genfromtxt(filename, delimiter=',', dtype=(float), \
                               skip_header=5, missing_values="", filling_values=-99.9)

    indata = np.ma.masked_where(indata == -99.9, indata)

    raobcore = utils.Timeseries("RAOBCORE v1.7", indata[:, 0], indata[:, 1])
    rich = utils.Timeseries("RICH v1.7", indata[:, 0], indata[:, 2])
    ratpac = utils.Timeseries("RATPAC A2", indata[:, 0], indata[:, 3])
 #   unsw = utils.Timeseries("UNSW v1.0", indata[:, 0], indata[:, 4])
    UAH = utils.Timeseries("UAH v6.0", indata[:, 0], indata[:, 4])
    rss = utils.Timeseries("RSS v4.0", indata[:, 0], indata[:, 5])
    noaa = utils.Timeseries("NOAA v4.1", indata[:, 0], indata[:, 6])

#    era5 = utils.Timeseries("ERA5", indata[:, 0], indata[:, 8])
    jra = utils.Timeseries("JRA-55", indata[:, 0], indata[:, 7])
    merra = utils.Timeseries("MERRA-2", indata[:, 0], indata[:, 8])

#    cmip = utils.Timeseries("CMIP5", indata[:, 0], indata[:, 11])
#    ssu3 = utils.Timeseries("SSU-3", indata[:, 0], indata[:, 12])

    return UAH, rss, ratpac, raobcore, rich, noaa, jra, merra # read_csv

#************************************************************************
def read_ssu_csv(filename):
    """
    Read user supplied CSV for LST into Timeseries object
    """

    indata = np.genfromtxt(filename, delimiter=',', dtype=(float), \
                               skip_header=5, missing_values="", filling_values=-99.9)

    indata = np.ma.masked_where(indata == -99.9, indata)
    ssu2 = utils.Timeseries("NOAA", indata[:, 0], indata[:, 2])
    ncar = utils.Timeseries("NCAR", indata[:, 0], indata[:, 4])

    return ssu2, ncar # read_ssu_csv

#************************************************************************
def read_ssu(filename):
    """
    Read user supplied CSV for LST into Timeseries object
    """

    indata = np.genfromtxt(filename, dtype=(float), \
                               skip_header=2, missing_values="", filling_values=-99.9)

    indata = np.ma.masked_where(indata == -99.9, indata)
    ssu1 = utils.Timeseries("SSU+MLS", indata[:, 0], indata[:, 1])
    ssu1_2 = utils.Timeseries("SSU+AMSU", indata[:, 0], indata[:, 2])
    ssu2 = utils.Timeseries("SSU+MLS", indata[:, 0], indata[:, 3])
    ssu2_2 = utils.Timeseries("SSU+AMSU", indata[:, 0], indata[:, 4])
    ssu3 = utils.Timeseries("SSU+MLS", indata[:, 0], indata[:, 5])
    ssu3_2 = utils.Timeseries("SSU+AMSU", indata[:, 0], indata[:, 6])

    return ssu1, ssu1_2, ssu2, ssu2_2, ssu3, ssu3_2 # read_ssu

#************************************************************************
def read_qbo_csv(filename):
    """
    Read user supplied CSV for QBO etc into Timeseries object
    """

    indata = np.genfromtxt(filename, delimiter=',', dtype=(float), skip_header=2, missing_values="", filling_values=-99.9)

    indata = np.ma.masked_where(indata == -99.9, indata)

    years = indata[:, 0]
    pentad = indata[:, 1]

    decimal_years = years + ((pentad - 1) * 5.) / 365.

    north = utils.Timeseries("North", decimal_years, indata[:, 2])
    south = utils.Timeseries("South", decimal_years, indata[:, 3])
    qbo = utils.Timeseries("QBO", decimal_years, indata[:, 4])


    return north, south, qbo # read_qbo_csv

#************************************************************************
def read_qbo_ncdf(filename, name):
    """
    Read user supplied netcdf for QBO into Timeseries object
    """

    cube_list = iris.load(filename)

    cube = cube_list[0]

    # get times and start time
    times = cube.coord("time")
    start = dt.datetime.strptime(times.units.origin.split()[-1], "%Y-%m-%d")

    # extract decimal years
    decimal_years = start.year + times.points.reshape(-1, 12)/12.
    decimal_years = decimal_years.reshape(-1)

    qbo = utils.Timeseries(name, decimal_years, cube.data)

    return qbo # read_qbo_ncdf

#************************************************************************
def read_merra_monthly(filename):
    """
    Read user supplied CSV for LST into Timeseries object
    """
    indata = np.genfromtxt(filename, dtype=(float), skip_header=1)

    # process the years and months to give decimals
    years = indata[:, 0]
    months = indata[:, 1]

    times = years + (months - 1)/12.

    fourh = utils.Timeseries("400", times, indata[:, 2])
    threeh = utils.Timeseries("300", times, indata[:, 3])
    twofiveh = utils.Timeseries("250", times, indata[:, 4])
    twoh = utils.Timeseries("200", times, indata[:, 5])
    onefiveh = utils.Timeseries("`50", times, indata[:, 6])
    fourtwoav = utils.Timeseries("Average", times, indata[:, 13])
    oneh = utils.Timeseries("100", times, indata[:, 7])
    seventy = utils.Timeseries("70", times, indata[:, 8])
    fifty = utils.Timeseries("50", times, indata[:, 9])
    thirty = utils.Timeseries("30", times, indata[:, 10])
    twenty = utils.Timeseries("20", times, indata[:, 11])
    ten = utils.Timeseries("10", times, indata[:, 12])
    seventytwentyav = utils.Timeseries("Average", times, indata[:, 14])

    return fourh, threeh, twofiveh, twoh, onefiveh, fourtwoav, oneh, \
        seventy, fifty, thirty, twenty, ten, seventytwentyav # read_merra_monthly

#************************************************************************
def read_zonal(filename, platform):
    """
    Read user supplied CSV for LST into Timeseries object
    """
    indata = np.genfromtxt(filename, dtype=(float), skip_header=1)

    print("15/6/17 - error in latitudes, needing to be flipped - Craig Long")
    if platform == "R":
        cfsr = utils.Timeseries("CFSR", indata[:, 1], -indata[:, 0])
        merra = utils.Timeseries("MERRA-2", indata[:, 2], -indata[:, 0])
        era = utils.Timeseries("ERA-Interim", indata[:, 3], -indata[:, 0])
        jra = utils.Timeseries("JRA-55", indata[:, 4], -indata[:, 0])

        return cfsr, merra, era, jra

    if platform == "S":
        star = utils.Timeseries("NOAA v3.0", indata[:, 1], -indata[:, 0])
        uah = utils.Timeseries("UAH v6.0", indata[:, 2], -indata[:, 0])
        rss = utils.Timeseries("RSS v3.3", indata[:, 3], -indata[:, 0])

        return star, uah, rss # read_zonal

#************************************************************************
def read_polar(filename):
    """
    Read user supplied data for LST into Timeseries object
    """
    indata = np.genfromtxt(filename, dtype=(float), skip_header=1)

    year = indata[:, 0]
    
    north = indata[:, 1]
    south = indata[:, 2]

    north = utils.Timeseries("North Pole", year, north)
    south = utils.Timeseries("South Pole -15"+r'$^\circ$'+"C", year, south)

    return north, south # read_polar

#************************************************************************
def read_era5(filename):
    """
    Read user supplied data for LST into Timeseries object
    """
    alldata = np.genfromtxt(filename, dtype=(float), skip_header=6, encoding="latin-1")

    years = alldata[:, 0]
    months = alldata[:, 1]
    land = alldata[:, 2]
    ocean = alldata[:, 3]
    combined = alldata[:, 4]

    times = years + (months - 1)/12.

    eral = utils.Timeseries("ERA5", times, land)
    erao = utils.Timeseries("ERA5", times, ocean)
    eralo = utils.Timeseries("ERA5", times, combined)

    return eral, erao, eralo # read_era5
#************************************************************************
def run_all_plots():

    #************************************************************************
    # Timeseries figures
    if True:

        fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 10), sharex=True)

        UAH, rss, ratpac, raobcore, rich, noaa, jra, merra = read_csv(DATALOC + "SotC_AnnTemps_2020_0220_LSTGL.csv")
#        ssu2, ncar = read_ssu_csv(DATALOC + "2018_LTT_LST_SSU_date0401_SSU.csv") 

        eral, erao, eralo = read_era5(DATALOC + "ERA5_TLS_GLOBAL")
        eralo_ann = utils.Timeseries("ERA5", np.reshape(eralo.times, [-1, 12])[:,0], utils.annual_average(eralo.data))

        # Sondes [no RATPAC for 2019]
        utils.plot_ts_panel(ax1, [raobcore, rich], "-", "temperature", loc=LEGEND_LOC)

        # satellites
        utils.plot_ts_panel(ax2, [UAH, noaa, rss], "-", "temperature", loc=LEGEND_LOC)

        # reanalyses
        jra_actuals, jra_anoms = utils.read_jra55(settings.REANALYSISLOC + "JRA-55_MSUch4_global_ts.txt", "temperature")
        merra_actuals, merra_anoms = utils.read_merra_LT_LS(settings.REANALYSISLOC + "MERRA2_MSU_Tanom_ann_{}.dat".format(settings.YEAR), LS=True)
        utils.plot_ts_panel(ax3, [eralo_ann, jra_anoms, merra_anoms], "-", "temperature", loc=LEGEND_LOC)

    #    ax3.set_ylabel("Anomaly ("+r'$^\circ$'+"C)", fontsize=settings.FONTSIZE)

        # Upper Stratosphere
#        utils.plot_ts_panel(ax4, [ssu2, ncar], "-", "temperature", loc=LEGEND_LOC)

        # sort formatting
        plt.xlim([1957, raobcore.times[-1]+1])

        for tick in ax3.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        for ax in [ax1, ax2, ax3]:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)

            ax.set_ylim([-1.2, 2.2])
 
        # sort labelling
        ax1.text(0.02, 0.88, "(a) Radiosondes", transform=ax1.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)
        ax2.text(0.02, 0.88, "(b) Satellites", transform=ax2.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)
        ax3.text(0.02, 0.88, "(c) Reanalyses", transform=ax3.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)
 
        fig.text(0.01, 0.45, "Anomaly ("+r'$^\circ$'+"C)", fontsize=settings.FONTSIZE, rotation="vertical")
        fig.subplots_adjust(right=0.98, top=0.98, bottom=0.04, hspace=0.001)

        plt.savefig(settings.IMAGELOC+"LST_ts{}".format(settings.OUTFMT))

        plt.close()


    #************************************************************************
    # Timeseries figures
    if True:

        fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 10), sharex=True)

        
        ssu1, ssu1_2, ssu2, ssu2_2, ssu3, ssu3_2 = read_ssu(DATALOC + "SSU.dat") 

        # SSU3
        utils.plot_ts_panel(ax1, [ssu3, ssu3_2], "-", "temperature", loc="")
        # SSU2
        utils.plot_ts_panel(ax2, [ssu2, ssu2_2], "-", "temperature", loc="")
        # SSU1
        utils.plot_ts_panel(ax3, [ssu1, ssu1_2], "-", "temperature", loc=LEGEND_LOC)

        # sort formatting
        plt.xlim([1957, ssu1.times[-1]+2])

        for tick in ax3.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        for ax in [ax1, ax2, ax3]:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)

            ax.set_ylim([-1.3, 2.2])
            
        # sort labelling
        ax1.text(0.02, 0.88, "(a) SSU3", transform=ax1.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)
        ax2.text(0.02, 0.88, "(b) SSU2", transform=ax2.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)
        ax3.text(0.02, 0.88, "(c) SSU1", transform=ax3.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)

        fig.text(0.01, 0.45, "Anomaly ("+r'$^\circ$'+"C)", fontsize=settings.FONTSIZE, rotation="vertical")
        fig.subplots_adjust(right=0.98, top=0.98, bottom=0.04, hspace=0.001)

        plt.savefig(settings.IMAGELOC+"LST_SSU_ts{}".format(settings.OUTFMT))

        plt.close()

    #************************************************************************
    # Combined Timeseries figures
    if True:

        fig = plt.figure(figsize=(12, 8))

        # manually set up the 10 axes
        w = 0.42
        h = 0.31
        c = 0.51
        ax1 = plt.axes([c-w, 0.99-h, w, h])
        ax2 = plt.axes([c, 0.99-h, w, h])
        ax3 = plt.axes([c-w, 0.99-(2*h), w, h], sharex=ax1)
        ax4 = plt.axes([c, 0.99-(2*h), w, h], sharex=ax2)
        ax5 = plt.axes([c-w, 0.99-(3*h), w, h], sharex=ax1)
        ax6 = plt.axes([c, 0.99-(3*h), w, h], sharex=ax2)

        UAH, rss, ratpac, raobcore, rich, noaa, jra, merra = read_csv(DATALOC + "SotC_AnnTemps_2020_0220_LSTGL.csv")
#        ssu2, ncar = read_ssu_csv(DATALOC + "2018_LTT_LST_SSU_date0401_SSU.csv") 

        eral, erao, eralo = read_era5(DATALOC + "ERA5_TLS_GLOBAL")
        eralo_ann = utils.Timeseries("ERA5", np.reshape(eralo.times, [-1, 12])[:,0], utils.annual_average(eralo.data))
        
        # update to ERA5.1 during revisions for SotC2019
        era51 = np.genfromtxt(DATALOC + "ERA5_1_update.dat")
        eralo_ann = utils.Timeseries("ERA5", era51[:, 0], era51[:, 1])

        # Sondes [no RATPAC for 2019]
        utils.plot_ts_panel(ax2, [raobcore, rich, ratpac], "-", "temperature", loc=LEGEND_LOC)

        # satellites
        utils.plot_ts_panel(ax4, [UAH, noaa, rss], "-", "temperature", loc=LEGEND_LOC)

        # reanalyses
        jra_actuals, jra_anoms = utils.read_jra55(settings.REANALYSISLOC + "JRA-55_MSUch4_global_ts.txt", "temperature")
        merra_actuals, merra_anoms = utils.read_merra_LT_LS(settings.REANALYSISLOC + "MERRA2_MSU_Tanom_ann_{}.dat".format(settings.YEAR), LS=True)
        utils.plot_ts_panel(ax6, [eralo_ann, jra_anoms, merra_anoms], "-", "temperature", loc=LEGEND_LOC)

        
        ssu1, ssu1_2, ssu2, ssu2_2, ssu3, ssu3_2 = read_ssu(DATALOC + "SSU.dat") 

        # SSU3
        utils.plot_ts_panel(ax1, [ssu3, ssu3_2], "-", "temperature", loc="")
        # SSU2
        utils.plot_ts_panel(ax3, [ssu2, ssu2_2], "-", "temperature", loc="")
        # SSU1
        utils.plot_ts_panel(ax5, [ssu1, ssu1_2], "-", "temperature", loc=LEGEND_LOC)
 
        # sort labelling
        for ax in [ax2, ax4, ax6]:
            ax.yaxis.tick_right()
            for tick in ax.yaxis.get_major_ticks():
                tick.label2.set_fontsize(settings.FONTSIZE)
        for ax in [ax1, ax3, ax5]:
            ax.yaxis.tick_right()
            ax.yaxis.set_ticks_position('left')

        ax2.text(0.02, 0.88, "(d) Radiosondes", transform=ax2.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)
        ax4.text(0.02, 0.88, "(e) Satellites", transform=ax4.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)
        ax6.text(0.02, 0.88, "(f) Reanalyses", transform=ax6.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)
 
        ax1.text(0.02, 0.88, "(a) SSU3", transform=ax1.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)
        ax3.text(0.02, 0.88, "(b) SSU2", transform=ax3.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)
        ax5.text(0.02, 0.88, "(c) SSU1", transform=ax5.transAxes, \
                     fontsize=settings.LABEL_FONTSIZE)

        plt.setp([a.get_xticklabels() for a in fig.axes[:-2]], visible=False)

        # sort formatting
        ax1.set_xlim([1957, raobcore.times[-1]+3])
        ax2.set_xlim([1957, raobcore.times[-1]+3])

        for ax in [ax5, ax6]:
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)
        for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)

            ax.set_ylim([-1.2, 2.2])         

        fig.text(0.01, 0.55, "Anomaly ("+r'$^\circ$'+"C)", fontsize=settings.FONTSIZE, rotation="vertical")
        fig.subplots_adjust(right=0.98, top=0.98, bottom=0.04, hspace=0.001)

        plt.savefig(settings.IMAGELOC+"LST_combined_ts{}".format(settings.OUTFMT))

        plt.close()

    #************************************************************************
    # Polar Figure
    if False:
        fig, ax1 = plt.subplots(figsize=(8, 5))

        north, south = read_polar(DATALOC + "polar_data_{}.dat".format(settings.YEAR))

        ax1.plot(north.times, north.data, c="b", ls="-", label=north.name)
        ax1.plot(south.times, south.data, c="r", ls="-", label=south.name)

        ax1.legend(loc=LEGEND_LOC, frameon=False, prop={'size':settings.FONTSIZE})
        ax1.set_xlim([1978, 2018+2])
        ax1.set_ylabel("Anomaly ("+r'$^\circ$'+"C)", fontsize=settings.FONTSIZE)

        utils.thicken_panel_border(ax1)

        plt.savefig(settings.IMAGELOC+"LST_polar_ts{}".format(settings.OUTFMT))

        for tick in ax1.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        for tick in ax1.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)

        plt.close()

    #************************************************************************
    # Polar and QBO Timeseries figures

    # fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 10), sharex=True)

    # north, south, qbo = read_qbo_csv(DATALOC + "SOC_Strat_Data_QBO.csv")

    # noaa_qbo = read_qbo_ncdf(DATALOC + "qbo_noaa.nc", "NOAA v4.0")
    # uah_qbo = read_qbo_ncdf(DATALOC + "qbo_uah.nc", "UAH v6.0")
    # rss_qbo = read_qbo_ncdf(DATALOC + "qbo_rss.nc", "RSS v3.3")


    # utils.plot_ts_panel(ax1, [north], "-", "temperature", loc="")
    # utils.plot_ts_panel(ax2, [south], "-", "temperature", loc="")
    # utils.plot_ts_panel(ax3, [noaa_qbo, uah_qbo, rss_qbo], "-", "temperature", loc="", lw=1)

    # lines3, labels3 = ax3.get_legend_handles_labels()

    # plt.figlegend(lines3, labels3, "lower center", frameon=False, ncol=3, fontsize=settings.FONTSIZE)

    # ax1.set_ylabel("Anomaly ("+r'$^\circ$'+"C)", fontsize=settings.FONTSIZE)
    # ax2.set_ylabel("Anomaly ("+r'$^\circ$'+"C)", fontsize=settings.FONTSIZE)
    # ax3.set_ylabel("QBO Index", fontsize=settings.FONTSIZE)

    # ax1.text(0.02, 0.88, "(a) North Polar Pentad Anomalies", transform=ax1.transAxes, \
    #              fontsize=settings.LABEL_FONTSIZE)
    # ax2.text(0.02, 0.88, "(b) South Polar Pentad Anomalies", transform=ax2.transAxes, \
    #              fontsize=settings.LABEL_FONTSIZE)
    # ax3.text(0.02, 0.88, "(c) QBO from Lower Stratospheric Temperature", transform=ax3.transAxes, \
    #              fontsize=settings.LABEL_FONTSIZE)

    # ax3.text(-0.15, 0.88, "West", transform=ax3.transAxes, fontsize=settings.LABEL_FONTSIZE)
    # ax3.text(-0.15, 0.01, "East", transform=ax3.transAxes, fontsize=settings.LABEL_FONTSIZE)

    # # sort formatting
    # plt.xlim([1978, int(settings.YEAR)+1])

    # for tick in ax3.xaxis.get_major_ticks():
    #     tick.label.set_fontsize(settings.FONTSIZE)
    # for ax in [ax1, ax2, ax3]:
    #     for tick in ax.yaxis.get_major_ticks():
    #         tick.label.set_fontsize(settings.FONTSIZE)

    # ax1.set_ylim([-10, 19])
    # ax2.set_ylim([-10, 23])
    # ax3.set_ylim([-1.1, 1.4])

    # fig.subplots_adjust(right=0.95, top=0.95, hspace=0.001)

    # plt.savefig(settings.IMAGELOC+"LST_qbo_ts{}".format(settings.OUTFMT))

    # plt.close()

    #************************************************************************
    # MERRA Trop and Strat Timeseries figures

    # fourh, threeh, twofiveh, twoh, onefiveh, fourtwoav, oneh, seventy, fifty, thirty, twenty, ten, seventytwentyav = read_merra_monthly(DATALOC + "MERRA2_400_20_temp_anom_1980-{}.txt".format(settings.YEAR))

    # fourtwoav.data = np.ma.masked_where(fourtwoav.times < 1994, fourtwoav.data)
    # seventytwentyav.data = np.ma.masked_where(seventytwentyav.times < 1994, seventytwentyav.data)
    # fourtwoav.times = np.ma.masked_where(fourtwoav.times < 1994, fourtwoav.times)
    # seventytwentyav.times = np.ma.masked_where(seventytwentyav.times < 1994, seventytwentyav.times)


    # fig, (ax1, ax2) = plt.subplots(2, figsize = (8, 8), sharex=True)

    # # trop
    # fit = utils.fit_plot_points(0.025, -50.01, fourtwoav.times)
    # ax1.plot(fourh.times, fit, c="0.5", lw = 2, ls = "-", label = "fit", zorder = 10)
    # utils.plot_ts_panel(ax1, [fourh, threeh, twofiveh, twoh, fourtwoav], "-", "lst", loc = LEGEND_LOC)
    # # strat
    # fit = utils.fit_plot_points(-0.0014, 2.4527, seventytwentyav.times)
    # ax2.plot(fourh.times, fit, c="0.5", lw = 2, ls = "-", label = "fit", zorder = 10)
    # utils.plot_ts_panel(ax2, [seventy, fifty, thirty, twenty, seventytwentyav], "-", "lst", loc = LEGEND_LOC)

    # # sort formatting
    # plt.xlim([1978, fourh.times[-1]+1])
    # ax1.set_ylim([-1.7,1.5])
    # ax2.set_ylim([-2.0,2.7])

    # for tick in ax2.xaxis.get_major_ticks():
    #     tick.label.set_fontsize(settings.FONTSIZE)
    # for ax in [ax1, ax2]:
    #     for tick in ax.yaxis.get_major_ticks():
    #         tick.label.set_fontsize(settings.FONTSIZE)

    # # sort labelling
    # ax1.text(0.02, 0.9, "(a) MERRA-2 Troposphere", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE)
    # ax2.text(0.02, 0.9, "(b) MERRA-2 Stratosphere", transform = ax2.transAxes, fontsize = settings.LABEL_FONTSIZE)

    # ax1.text(0.7, 0.9, "y=0.025x - 50.01", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE*0.8)
    # ax2.text(0.7, 0.9, "y=-0.0014x + 2.4527", transform = ax2.transAxes, fontsize = settings.LABEL_FONTSIZE*0.8)

    # ax1.set_ylabel("Anomaly ("+r'$^\circ$'+"C)", fontsize = settings.FONTSIZE)
    # ax2.set_ylabel("Anomaly ("+r'$^\circ$'+"C)", fontsize = settings.FONTSIZE)

    # fig.subplots_adjust(right = 0.95, top = 0.95, hspace = 0.001)

    # plt.savefig(settings.IMAGELOC+"LST_merra_ts{}".format(settings.OUTFMT))

    # plt.close()

    #************************************************************************
    # Zonal figures

    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (8, 6.5), sharey=True)

    # cfsr, merra, era, jra = read_zonal(DATALOC + "TLS_Reanal_zonal_trends_1994-{}.txt".format(settings.YEAR), "R")
    # star, uah, rss = read_zonal(DATALOC + "TLS_Satellite_zonal_trends_1994-{}.txt".format(settings.YEAR), "S")

    # # Reanalyses
    # utils.plot_ts_panel(ax1, [cfsr, merra, era, jra], "--", "temperature", loc = "lower left", ncol = 1)

    # # Satellite
    # utils.plot_ts_panel(ax2, [star, uah, rss], "--", "temperature", loc = "center right", ncol = 1)

    # ax1.axvline(0,color = "0.5", ls = "--")
    # ax2.axvline(0,color = "0.5", ls = "--")

    # # sort formatting
    # plt.ylim([-90,90])
    # ax1.set_ylabel("Latitude", fontsize = settings.LABEL_FONTSIZE)
    # ax1.set_xlabel("Trend ("+r'$^\circ$'+"C decade"+r'$^{-1}$'+")", fontsize = settings.FONTSIZE)
    # ax2.set_xlabel("Trend ("+r'$^\circ$'+"C decade"+r'$^{-1}$'+")", fontsize = settings.FONTSIZE)


    # for tick in ax1.yaxis.get_major_ticks():
    #     tick.label.set_fontsize(settings.FONTSIZE)
    # for ax in [ax1, ax2]:
    #     ax.set_xlim([-0.3,0.6])
    #     ax.set_xticks(np.arange(-0.2,0.8,0.2))
    #     ax.set_yticks(np.arange(-90,120,30))
    #     for tick in ax.xaxis.get_major_ticks():
    #         tick.label.set_fontsize(settings.FONTSIZE)

    # # sort labelling
    # ax1.text(0.02, 0.9, "(a) Radiosondes", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE)
    # ax2.text(0.02, 0.9, "(b) Satellites", transform = ax2.transAxes, fontsize = settings.LABEL_FONTSIZE)

    # fig.subplots_adjust(right = 0.95, top = 0.95, hspace = 0.001)

    # plt.savefig(settings.IMAGELOC+"LST_profiles{}".format(settings.OUTFMT))

    # plt.close()


 
    #************************************************************************
    # ERA5 Anomaly figure
    if True:
        # Read in ERA anomalies

        cube_list = iris.load(settings.REANALYSISLOC + "era5_tls_{}01-{}12_ann_ano.nc".format(settings.YEAR, settings.YEAR))
        
        cube = cube_list[0]
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
        
        bounds=[-4, -1.2, -0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.8, 1.2, 4]
        
        utils.plot_smooth_map_iris(settings.IMAGELOC + "LST_{}_anoms_era5".format(settings.YEAR), cube[0], \
                                   settings.COLOURMAP_DICT["temperature"], bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", title="ERA5")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_LST_{}_anoms_era5".format(settings.YEAR), cube[0], \
                                   settings.COLOURMAP_DICT["temperature"], bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", \
                                   figtext="(f) Lower Stratosphere Temperature")



    #************************************************************************
    # merra Anomaly figure
    if False:
        import netCDF4 as ncdf

        ncfile = ncdf.Dataset(settings.REANALYSISLOC + "merra2_tls_ANNUAL_anom.nc")

        var=ncfile.variables["TLS ANOM"][:] # this is a masked array
        nlons = ncfile.variables["LONGITUDES"][:]
        nlats = ncfile.variables["LATITUDES"][:]

        cube = utils.make_iris_cube_2d(var, nlats, nlons, "LST", "C")

        utils.plot_smooth_map_iris(settings.IMAGELOC + "LST_{}_anoms_merra".format(settings.YEAR), cube, settings.COLOURMAP_DICT["temperature"], bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)")


    #************************************************************************
    # 2015 MERRA seasonal figure
    
    # import netCDF4 as ncdf

    # month_list = []

    # for month in MONTHS:
    #     print month

    #     # IRIS doesn't like the whitespace in the "TLS ANOM"
    #     ncfile=ncdf.Dataset(settings.REANALYSISLOC + "merra2_tls_{}_anom.nc".format(month.upper()),'r')

    #     var=ncfile.variables["TLS ANOM"][:] # this is a masked array
    #     lons = ncfile.variables["LONGITUDES"][:]
    #     lats = ncfile.variables["LATITUDES"][:]

    #     ncfile.close()

    #     cube = utils.make_iris_cube_2d(var, lats, lons, "TLS_ANOM", "C")

    #     month_list += [cube]

    # # pass to plotting routine
    # utils.plot_smooth_map_iris_multipanel(settings.IMAGELOC + "LST_{}_monthly_merra".format(settings.YEAR), month_list, \
    #                                           settings.COLOURMAP_DICT["temperature"], bounds, \
    #                                           "Anomaly ("+r'$^{\circ}$'+"C)", shape = (6,2), \
    #                                           title = MONTHS, \
    #                                           figtext = ["(a)","(b)","(c)","(d)", "(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"])

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
