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
# python3
from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import settings
import utils

DATALOC = "{}/{}/data/HUM/".format(settings.ROOTLOC, settings.YEAR)

LEGEND_LOC = 'upper left'

LW = 3
BBOX = (0, 0.9)


#*********************************************
def read_ts(filename, var, domain):

    indata = np.genfromtxt(filename, skip_header=1, dtype=float)

    years = indata[:, 0]

    indata = np.ma.masked_where(indata <= -99.9, indata)

    if domain == "L":
        if var == "q":
            off = 0
        elif var == "rh":
            off = 22

        hadisdh = utils.Timeseries("HadISDH", years, indata[:, 1+off])
        hadcruh = utils.Timeseries("HadCRUH", years, indata[:, 2+off])
        hadcruhext = utils.Timeseries("HadCRUHExt", years, indata[:, 3+off])
        dai = utils.Timeseries("Dai", years, indata[:, 4+off])
        era5_msk = utils.Timeseries("ERA5 mask", years, indata[:, 5+off])
        merra_msk = utils.Timeseries("MERRA-2 mask", years, indata[:, 6+off])
        jra_msk = utils.Timeseries("JRA-55 mask", years, indata[:, 7+off])
        era5 = utils.Timeseries("ERA5", years, indata[:, 8+off])
        erai = utils.Timeseries("ERA-Interim", years, indata[:, 9+off])
        merra = utils.Timeseries("MERRA-2", years, indata[:, 10+off])
        jra = utils.Timeseries("JRA-55", years, indata[:, 11+off])
        cr20 = utils.Timeseries("20CRv3", years, indata[:, 12+off])

        hadcruhext.ls = "--"
        era5_msk.ls = "--"
        erai.ls = "--"
        merra_msk.ls = "--"
        jra_msk.ls = "--"

    elif domain == "M":
        if var == "q":
            off = 12
            extra = 0
        elif var == "rh":
            off = 34
            extra = -2

        hadisdh = utils.Timeseries("HadISDH", years, indata[:, 1+off])
        hadcruh = utils.Timeseries("HadCRUH", years, indata[:, 2+off])
        dai = utils.Timeseries("Dai", years, indata[:, 3+off])
        nocs = utils.Timeseries("NOCS v2.0", years, indata[:, 4+off]) # not in RH
        hoaps = utils.Timeseries("HOAPS", years, indata[:, 5+off]) #  not in RH
        era5 = utils.Timeseries("ERA5", years, indata[:, 6+off+extra])
        erai = utils.Timeseries("ERA-Interim", years, indata[:, 7+off+extra])
        merra = utils.Timeseries("MERRA-2", years, indata[:, 8+off+extra])
        jra = utils.Timeseries("JRA-55", years, indata[:, 9+off+extra])
        cr20 = utils.Timeseries("20CRv3", years, indata[:, 10+off+extra])

    # mask out early ERA data
    pre79, = np.where(years < 1979)
    erai.data.mask[pre79] = True
    erai.ls = "--"


    if domain == "L":
        return [hadisdh, hadcruh, hadcruhext, dai, erai, era5, merra, jra, era5_msk, merra_msk, cr20] # read_ts
    elif domain == "M":
        return [hadisdh, hadcruh, dai, nocs, hoaps, erai, era5, merra, jra, cr20] # read_ts

#*********************************************
def read_ts_unc(filename, var, domain):

    indata = np.genfromtxt(filename, skip_header=1, dtype=float)

    years = indata[:, 0]

    indata = np.ma.masked_where(indata <= -99.9, indata)

    if domain == "L":
        if var == "q":
            off = 0
        elif var == "rh":
            off = 18

        hadisdh = utils.Timeseries("HadISDH", years, indata[:, 1+off])
        hadisdh_l = utils.Timeseries("HadISDH", years, indata[:, 2+off])
        hadisdh_u = utils.Timeseries("HadISDH", years, indata[:, 3+off])
        era5_msk = utils.Timeseries("ERA5 mask", years, indata[:, 4+off])
        merra_msk = utils.Timeseries("MERRA-2 mask", years, indata[:, 5+off])
        jra_msk = utils.Timeseries("JRA-55 mask", years, indata[:, 6+off])
        era5 = utils.Timeseries("ERA5", years, indata[:, 7+off])
        merra = utils.Timeseries("MERRA-2", years, indata[:, 8+off])
        jra = utils.Timeseries("JRA-55", years, indata[:, 9+off])
        cr20 = utils.Timeseries("20CRv3", years, indata[:, 10+off])

        era5_msk.ls = "--"
        merra_msk.ls = "--"
        jra_msk.ls = "--"

    elif domain == "M":
        if var == "q":
            off = 10
            extra = 0
        elif var == "rh":
            off = 28
            extra = -1

        hadisdh = utils.Timeseries("HadISDH", years, indata[:, 1+off])
        hadisdh_l = utils.Timeseries("HadISDH", years, indata[:, 2+off])
        hadisdh_u = utils.Timeseries("HadISDH", years, indata[:, 3+off])
        nocs = utils.Timeseries("NOCS v2.0", years, indata[:, 4+off]) # not in RH
        era5 = utils.Timeseries("ERA5", years, indata[:, 5+off+extra])
        merra = utils.Timeseries("MERRA-2", years, indata[:, 6+off+extra])
        jra = utils.Timeseries("JRA-55", years, indata[:, 7+off+extra])
        cr20 = utils.Timeseries("20CRv3", years, indata[:, 8+off+extra])

    hadisdh.zorder = 10

    if domain == "L":
        return [hadisdh, hadisdh_l, hadisdh_u, era5, merra, jra, era5_msk, merra_msk, cr20] # read_ts_unc
    elif domain == "M":
        return [hadisdh, hadisdh_l, hadisdh_u, nocs, era5, merra, jra, cr20] # read_ts_unc


#*********************************************
def read_maps(filename, name, units, footer=False):

    if footer:
        indata = np.genfromtxt(filename, dtype=(float), skip_footer=2)
    else:
        indata = np.genfromtxt(filename, dtype=(float))

    indata = np.ma.masked_where(indata <= -99.999, indata)

    indata = indata[::-1, :]


    nlat, nlon = indata.shape

    delta_lat = 180./nlat
    delta_lon = 360./nlon
    
    # presume -90 --> 90, -180  --> 180

    LATS = np.arange(-90. + (delta_lat/2.), 90. + (delta_lat/2.), delta_lat)
    LONS = np.arange(-180. + (delta_lon/2.), 180. + (delta_lon/2.), delta_lon)

    cube = utils.make_iris_cube_2d(indata, LATS, LONS, name, units)

    return cube # read_maps

#************************************************************************
def run_all_plots():

    #*********************************************
    # Timeseries plot
    if False:
        (hadisdhLQ, hadcruhLQ, hadcruhextLQ, daiLQ, eraiLQ, era5LQ, merraLQ, jraLQ, era5_mskLQ, merra_mskLQ, cr20LQ) = \
            read_ts(DATALOC + "HUM_timeseries_ALL{}.txt".format(settings.YEAR), "q", "L")
        (hadisdhMQ, hadcruhMQ, daiMQ, nocsMQ, hoapsMQ, eraiMQ, era5MQ, merraMQ, jraMQ, cr20MQ) = \
            read_ts(DATALOC + "HUM_timeseries_ALL{}.txt".format(settings.YEAR), "q", "M")
        (hadisdhLR, hadcruhLR, hadcruhextLR, daiLR, eraiLR, era5LR, merraLR, jraLR, era5_mskLR, merra_mskLR, cr20LR) = \
            read_ts(DATALOC + "HUM_timeseries_ALL{}.txt".format(settings.YEAR), "rh", "L")
        (hadisdhMR, hadcruhMR, daiMR, nocsMR, hoapsMR, eraiMR, era5MR, merraMR, jraMR, cr20MR) = \
            read_ts(DATALOC + "HUM_timeseries_ALL{}.txt".format(settings.YEAR), "rh", "M")

        COLOURS = settings.COLOURS["hydrological"]
        fig = plt.figure(figsize=(14, 12))

        # manually set up the 8 axes
        w = 0.45 # width
        h = 0.24 # height
        c = 0.53 # centre line of plots
        ax1 = plt.axes([c-w, 0.99-h, w, h])
        ax2 = plt.axes([c, 0.99-h, w, h])
        ax3 = plt.axes([c-w, 0.99-(2*h), w, h], sharex=ax1)
        ax4 = plt.axes([c, 0.99-(2*h), w, h], sharex=ax2)
        ax5 = plt.axes([c-w, 0.99-(3*h), w, h], sharex=ax1)
        ax6 = plt.axes([c, 0.99-(3*h), w, h], sharex=ax2)
        ax7 = plt.axes([c-w, 0.99-(4*h), w, h], sharex=ax1)
        ax8 = plt.axes([c, 0.99-(4*h), w, h], sharex=ax2)

        if int(settings.YEAR) < 2019:
            # in situ
            utils.plot_ts_panel(ax1, [hadisdhLQ, hadcruhLQ, hadcruhextLQ, daiLQ, era5_mskLQ], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)
            utils.plot_ts_panel(ax3, [hadisdhMQ, hadcruhMQ, daiMQ, nocsMQ, hoapsMQ], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)
            utils.plot_ts_panel(ax5, [hadisdhLR, hadcruhLR, hadcruhextLR, daiLR, era5_mskLR], "-", "hydrological", loc="")
            utils.plot_ts_panel(ax7, [hadisdhMR, hadcruhMR, daiMR], "-", "hydrological", loc="")

            # reanalyses
            utils.plot_ts_panel(ax2, [eraiLQ, era5LQ, merraLQ, jraLQ], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)
            utils.plot_ts_panel(ax4, [eraiMQ, era5MQ, merraMQ, jraMQ], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)
            utils.plot_ts_panel(ax6, [eraiLR, era5LR, jraLR], "-", "hydrological", loc="")
            utils.plot_ts_panel(ax8, [eraiMR, era5MR, jraMR], "-", "hydrological", loc="")
        else:
            # in situ
            utils.plot_ts_panel(ax1, [hadisdhLQ], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)
            utils.plot_ts_panel(ax3, [hadisdhMQ, nocsMQ], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)
            utils.plot_ts_panel(ax5, [hadisdhLR], "-", "hydrological", loc="")
            utils.plot_ts_panel(ax7, [hadisdhMR], "-", "hydrological", loc="")

            # hadisdh uncertainties in due course

            # reanalyses
            utils.plot_ts_panel(ax2, [era5LQ, merraLQ, jraLQ, cr20LQ], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)
            utils.plot_ts_panel(ax4, [era5MQ, merraMQ, jraMQ, cr20MQ], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)
            utils.plot_ts_panel(ax6, [era5LR, jraLR, cr20LR], "-", "hydrological", loc="")
            utils.plot_ts_panel(ax8, [era5MR, jraMR, cr20MR], "-", "hydrological", loc="")

        # prettify
        ax1.set_xlim([1957, int(settings.YEAR)+2])
        ax2.set_xlim([1957, int(settings.YEAR)+2])
        ax1.set_xticklabels(["", "1960", "1970", "1980", "1990", "2000", "2010", ""])

        for ax in [ax1, ax2, ax3, ax4]:
            ax.set_ylim([-0.39, 0.8])
        for ax in [ax5, ax6, ax7, ax8]:
            ax.set_ylim([-1.8, 1.5])

        for ax in [ax1, ax3, ax5, ax7]:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)
        for tick in ax7.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        for tick in ax8.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)

        plt.setp([a.get_xticklabels() for a in fig.axes[:-2]], visible=False)
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.setp(ax4.get_yticklabels(), visible=False)
        plt.setp(ax6.get_yticklabels(), visible=False)
        plt.setp(ax8.get_yticklabels(), visible=False)


        ax1.text(0.02, 0.85, "(a) In Situ Land q", transform=ax1.transAxes, fontsize=settings.FONTSIZE)
        ax2.text(0.02, 0.85, "(b) Reanalyses Land q", transform=ax2.transAxes, fontsize=settings.FONTSIZE)
        ax3.text(0.02, 0.85, "(c) In Situ Ocean q", transform=ax3.transAxes, fontsize=settings.FONTSIZE)
        ax4.text(0.02, 0.85, "(d) Reanalyses Ocean q", transform=ax4.transAxes, fontsize=settings.FONTSIZE)
        ax5.text(0.02, 0.85, "(e) In Situ Land RH", transform=ax5.transAxes, fontsize=settings.FONTSIZE)
        ax6.text(0.02, 0.85, "(f) Reanalyses Land RH", transform=ax6.transAxes, fontsize=settings.FONTSIZE)
        ax7.text(0.02, 0.85, "(g) In Situ Ocean RH", transform=ax7.transAxes, fontsize=settings.FONTSIZE)
        ax8.text(0.02, 0.85, "(h) Reanalyses Ocean RH", transform=ax8.transAxes, fontsize=settings.FONTSIZE)


        plt.figtext(0.01, 0.75, "Specific Humidity (g kg"+r'$^{-1}$'+")", va='center', rotation='vertical', fontsize=settings.FONTSIZE)
        plt.figtext(0.01, 0.25, "Relative Humidity (%rh)", va='center', rotation='vertical', fontsize=settings.FONTSIZE)

        fig.subplots_adjust(right=0.98, top=0.95, bottom=0.05, hspace=0.001)

        plt.savefig(settings.IMAGELOC + "HUM_ts{}".format(settings.OUTFMT))
        plt.close()

    #*********************************************
    # Timeseries uncertainty plot
    if True:
        (hadisdhLQ, hadisdhLQ_l, hadisdhLQ_u, era5LQ, merraLQ, jraLQ, era5_mskLQ, merra_mskLQ, cr20LQ) = \
            read_ts_unc(DATALOC + "HUM_timeseries_ALL{}_unc.txt".format(settings.YEAR), "q", "L")
        (hadisdhMQ, hadisdhMQ_l, hadisdhMQ_u, nocsMQ, era5MQ, merraMQ, jraMQ, cr20MQ) = \
            read_ts_unc(DATALOC + "HUM_timeseries_ALL{}_unc.txt".format(settings.YEAR), "q", "M")
        (hadisdhLR, hadisdhLR_l, hadisdhLR_u,  era5LR, merraLR, jraLR, era5_mskLR, merra_mskLR, cr20LR) = \
            read_ts_unc(DATALOC + "HUM_timeseries_ALL{}_unc.txt".format(settings.YEAR), "rh", "L")
        (hadisdhMR, hadisdhMR_l, hadisdhMR_u, nocsMR, era5MR, merraMR, jraMR, cr20MR) = \
            read_ts_unc(DATALOC + "HUM_timeseries_ALL{}_unc.txt".format(settings.YEAR), "rh", "M")

        COLOURS = settings.COLOURS["hydrological"]
        fig = plt.figure(figsize=(14, 12))

        # manually set up the 8 axes
        w = 0.45 # width
        h = 0.24 # height
        c = 0.53 # centre line of plots
        ax1 = plt.axes([c-w, 0.99-h, w, h])
        ax2 = plt.axes([c, 0.99-h, w, h])
        ax3 = plt.axes([c-w, 0.99-(2*h), w, h], sharex=ax1)
        ax4 = plt.axes([c, 0.99-(2*h), w, h], sharex=ax2)
        ax5 = plt.axes([c-w, 0.99-(3*h), w, h], sharex=ax1)
        ax6 = plt.axes([c, 0.99-(3*h), w, h], sharex=ax2)
        ax7 = plt.axes([c-w, 0.99-(4*h), w, h], sharex=ax1)
        ax8 = plt.axes([c, 0.99-(4*h), w, h], sharex=ax2)

        # in situ
        utils.plot_ts_panel(ax1, [hadisdhLQ], "-", "hydrological", loc="",)
        ax1.fill_between(hadisdhLQ.times, hadisdhLQ_u.data, hadisdhLQ_l.data, color='0.8', label="")
        utils.plot_ts_panel(ax3, [hadisdhMQ, nocsMQ], "-", "hydrological", loc="")
        ax3.fill_between(hadisdhMQ.times, hadisdhMQ_u.data, hadisdhMQ_l.data, color='0.8', label="")
        utils.plot_ts_panel(ax5, [hadisdhLR], "-", "hydrological", loc="")
        ax5.fill_between(hadisdhLR.times, hadisdhLR_u.data, hadisdhLR_l.data, color='0.8', label="")
        utils.plot_ts_panel(ax7, [hadisdhMR], "-", "hydrological", loc="")
        ax7.fill_between(hadisdhMR.times, hadisdhMR_u.data, hadisdhMR_l.data, color='0.8', label="")

    
        # reanalyses
        utils.plot_ts_panel(ax2, [era5LQ, merraLQ, jraLQ, cr20LQ], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)
        utils.plot_ts_panel(ax4, [era5MQ, merraMQ, jraMQ, cr20MQ], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)
        utils.plot_ts_panel(ax6, [era5LR, jraLR, cr20LR], "-", "hydrological", loc="")
        utils.plot_ts_panel(ax8, [era5MR, jraMR, cr20MR], "-", "hydrological", loc="")

        # fix legend
        unc_patch = ax1.fill(np.NaN, np.NaN, '0.8', zorder = 1, alpha=0.7)

        for ax in [ax1, ax3]:
            lines = []
            labels = []
            for line in ax.get_lines():
                if line.get_label() == "HadISDH":
                    lines += [(line, unc_patch[0])]
                    labels += [line.get_label()]                   
                else:
                    lines += [line]
                    labels += [line.get_label()]


            ax.legend(lines, labels, \
                       loc=LEGEND_LOC, ncol=2, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, \
                       labelspacing=0.1, columnspacing=0.5, bbox_to_anchor=BBOX)

        # prettify
        ax1.set_xlim([1957, int(settings.YEAR)+2])
        ax2.set_xlim([1957, int(settings.YEAR)+2])
        ax1.set_xticklabels(["", "1960", "1970", "1980", "1990", "2000", "2010", ""])

        for ax in [ax1, ax2, ax3, ax4]:
            ax.set_ylim([-0.39, 0.8])
        for ax in [ax5, ax6, ax7, ax8]:
            ax.set_ylim([-1.8, 1.5])

        for ax in [ax1, ax3, ax5, ax7]:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)
        for tick in ax7.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        for tick in ax8.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)

        plt.setp([a.get_xticklabels() for a in fig.axes[:-2]], visible=False)
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.setp(ax4.get_yticklabels(), visible=False)
        plt.setp(ax6.get_yticklabels(), visible=False)
        plt.setp(ax8.get_yticklabels(), visible=False)


        ax1.text(0.02, 0.85, "(a) In Situ Land q", transform=ax1.transAxes, fontsize=settings.FONTSIZE)
        ax2.text(0.02, 0.85, "(b) Reanalyses Land q", transform=ax2.transAxes, fontsize=settings.FONTSIZE)
        ax3.text(0.02, 0.85, "(c) In Situ Ocean q", transform=ax3.transAxes, fontsize=settings.FONTSIZE)
        ax4.text(0.02, 0.85, "(d) Reanalyses Ocean q", transform=ax4.transAxes, fontsize=settings.FONTSIZE)
        ax5.text(0.02, 0.85, "(e) In Situ Land RH", transform=ax5.transAxes, fontsize=settings.FONTSIZE)
        ax6.text(0.02, 0.85, "(f) Reanalyses Land RH", transform=ax6.transAxes, fontsize=settings.FONTSIZE)
        ax7.text(0.02, 0.85, "(g) In Situ Ocean RH", transform=ax7.transAxes, fontsize=settings.FONTSIZE)
        ax8.text(0.02, 0.85, "(h) Reanalyses Ocean RH", transform=ax8.transAxes, fontsize=settings.FONTSIZE)


        plt.figtext(0.01, 0.75, "Specific Humidity (g kg"+r'$^{-1}$'+")", va='center', rotation='vertical', fontsize=settings.FONTSIZE)
        plt.figtext(0.01, 0.25, "Relative Humidity (%rh)", va='center', rotation='vertical', fontsize=settings.FONTSIZE)

        fig.subplots_adjust(right=0.98, top=0.95, bottom=0.05, hspace=0.001)

        plt.savefig(settings.IMAGELOC + "HUM_ts_unc{}".format(settings.OUTFMT))
        plt.close()

    input("stOP")
    #*********************************************
    # Map plots
    ## RH
    bounds = [-50, -12, -9, -6, -3, 0, 3, 6, 9, 12, 50]
    # RH HadISDH
    if False:
        cube = read_maps(DATALOC + "HUMrh_anomalymap_HADISDHland{}.txt".format(settings.YEAR), "HadISDH RH", None, footer=True)

        utils.plot_smooth_map_iris(settings.IMAGELOC + "HUM_RH_hadisdh_land", cube, settings.COLOURMAP_DICT["hydrological"], \
                                       bounds, "Anomalies from 1981-2010 (%rh)", figtext="", title="")

    # RH HadISDH Land and Marine
    if True:
   
        cube = read_maps(DATALOC + "HUMrh_anomalymap_HADISDH{}.txt".format(settings.YEAR), "HadISDH RH", None, footer=True)

        utils.plot_smooth_map_iris(settings.IMAGELOC + "HUM_RH_hadisdh_combined", cube, settings.COLOURMAP_DICT["hydrological"], \
                                   bounds, "Anomalies from 1981-2010 (%rh)", figtext="", title="")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_HUM_RH_hadisdh_combined", cube, settings.COLOURMAP_DICT["hydrological"], \
                                   bounds, "Anomalies from 1981-2010 (%rh)", figtext="(h) Surface Relative Humidity", title="")

    # RH ERA 
    if True:
        cube = read_maps(DATALOC + "HUMrh_anomalymap_ERA5{}.txt".format(settings.YEAR), "ERA-I RH", None, footer=True)
        
        utils.plot_smooth_map_iris(settings.IMAGELOC + "HUM_RH_era5", cube, settings.COLOURMAP_DICT["hydrological"], \
                                   bounds, "Anomalies from 1981-2010 (%rh)", figtext="", title="")
        #    utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_HUM_RH_era5", cube, settings.COLOURMAP_DICT["hydrological"], \
            #                                   bounds, "Anomalies from 1981-2010 (%rh)", figtext="(o) Surface Relative Humidity", title="")

    # RH MERRA 
    if True:
        cube = read_maps(DATALOC + "HUMrh_anomalymap_MERRA2{}.txt".format(settings.YEAR), "MERRA2 RH", None, footer=True)
   
        utils.plot_smooth_map_iris(settings.IMAGELOC + "HUM_RH_merra", cube, settings.COLOURMAP_DICT["hydrological"], \
                                   bounds, "Anomalies from 1981-2010 (%rh)", figtext="", title="")

    ## Q
    bounds = [-20., -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 20]
    # q HadISDH
    if False:
        cube = read_maps(DATALOC + "HUMq_anomalymap_HADISDHland{}.txt".format(settings.YEAR), "HadISDH q", "g/kg", footer=True)

        utils.plot_smooth_map_iris(settings.IMAGELOC + "HUM_q_hadisdh_land", cube, settings.COLOURMAP_DICT["hydrological"], \
                                   bounds, "Anomalies from 1981-2010 (g kg"+r'$^{-1}$'+")", figtext="", title="")
        #    utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_HUM_q_hadisdh_land", cube, settings.COLOURMAP_DICT["hydrological"], \
            #                                   bounds, "Anomalies from 1981-2010 (g kg"+r'$^{-1}$'+")", figtext="(n) Surface Specific Humidity", title="")

    # q HadISDH Land and Marine
    if True:
        cube = read_maps(DATALOC + "HUMq_anomalymap_HADISDH{}.txt".format(settings.YEAR), "HadISDH q", "g/kg", footer=True)

        utils.plot_smooth_map_iris(settings.IMAGELOC + "HUM_q_hadisdh_combined", cube, settings.COLOURMAP_DICT["hydrological"], \
                                   bounds, "Anomalies from 1981-2010 (g kg"+r'$^{-1}$'+")", figtext="", title="")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_HUM_q_hadisdh_combined", cube, settings.COLOURMAP_DICT["hydrological"], \
                                   bounds, "Anomalies from 1981-2010 (g kg"+r'$^{-1}$'+")", figtext="(g) Surface Specific Humidity", title="")

    # q ERA
    if True:
        cube = read_maps(DATALOC + "HUMq_anomalymap_ERA5{}.txt".format(settings.YEAR), "ERA-I q", "g/kg", footer=True)

        utils.plot_smooth_map_iris(settings.IMAGELOC + "HUM_q_era5", cube, settings.COLOURMAP_DICT["hydrological"], \
                                   bounds, "Anomalies from 1981-2010 (g kg"+r'$^{-1}$'+")", figtext="", title="")

    # MERRA?
    if True:
        cube = read_maps(DATALOC + "HUMq_anomalymap_MERRA2{}.txt".format(settings.YEAR), "MERRA2 q", "g/kg", footer=True)

        utils.plot_smooth_map_iris(settings.IMAGELOC + "HUM_q_merra", cube, settings.COLOURMAP_DICT["hydrological"], \
                                   bounds, "Anomalies from 1981-2010 (g kg"+r'$^{-1}$'+")", figtext="", title="")



    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
