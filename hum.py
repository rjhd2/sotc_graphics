#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for Surface Humidity (HUM) section.
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
from mpl_toolkits.basemap import Basemap

import matplotlib.cm as mpl_cm
import matplotlib as mpl

import settings
import utils

data_loc = "/data/local/rdunn/SotC/{}/data/HUM/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

LEGEND_LOC = 'upper left'

LW = 3
BBOX = (0,0.9)


#*********************************************
def read_ts(filename, var, domain):

    indata = np.genfromtxt(filename, skip_header = 1, dtype = float)

    years = indata[:,0]

    indata = np.ma.masked_where(indata <= -99.9, indata)

    if domain == "L":
        if var == "q":
            off = 0
        elif var == "rh":
            off = 19

        hadisdh = utils.Timeseries("HadISDH", years, indata[:,1+off])
        hadcruh = utils.Timeseries("HadCRUH", years, indata[:,2+off])
        hadcruhext = utils.Timeseries("HadCRUHExt", years, indata[:,3+off])
        dai = utils.Timeseries("Dai", years, indata[:,4+off])
        era = utils.Timeseries("ERA-Interim", years, indata[:,5+off])
        merra = utils.Timeseries("MERRA-2", years, indata[:,6+off])
#        ncep = utils.Timeseries("NCEP", years, indata[:,7+off])
        jra = utils.Timeseries("JRA-55", years, indata[:,8+off])
#        cr20 = utils.Timeseries("20CR", years, indata[:,9+off])

        hadcruhext.ls = "--"

    elif domain == "M":
        if var == "q":
            off = 9
            extra = 2
        elif var == "rh":
            off = 28
            extra = 0

        hadisdh = utils.Timeseries("HadISDH", years, indata[:,1+off])
        hadcruh = utils.Timeseries("HadCRUH", years, indata[:,2+off])
        dai = utils.Timeseries("Dai", years, indata[:,3+off])
        nocs = utils.Timeseries("NOCS v2.0", years, indata[:,4+off])
        hoaps = utils.Timeseries("HOAPS", years, indata[:,5+off])
        era = utils.Timeseries("ERA-Interim", years, indata[:,4+off+extra])
        merra = utils.Timeseries("MERRA-2", years, indata[:,5+off+extra])
#        ncep = utils.Timeseries("NCEP", years, indata[:,6+off+extra])
        jra = utils.Timeseries("JRA-55", years, indata[:,7+off+extra])
#        cr20 = utils.Timeseries("20CR", years, indata[:,8+off+extra])

    # mask out early ERA data
    pre79, = np.where(years < 1979)
    era.data.mask[pre79] = True


    if domain == "L":
        return [hadisdh, hadcruh, hadcruhext, dai, era, merra, jra] # read_ts
    elif domain == "M":
        return [hadisdh, hadcruh, dai, nocs, hoaps, era, merra, jra] # read_ts


#*********************************************
def read_maps(filename, name, units, footer = False):

    if footer:
        indata = np.genfromtxt(filename, dtype=(float), skip_footer = 2)
    else:
        indata = np.genfromtxt(filename, dtype=(float))

    indata = np.ma.masked_where(indata <= -99.999, indata)

    indata = indata[::-1,:]


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

    land_q = read_ts(data_loc + "HUM_timeseries_ALL{}.txt".format(settings.YEAR), "q", "L")
    marine_q = read_ts(data_loc + "HUM_timeseries_ALL{}.txt".format(settings.YEAR), "q", "M")
    land_rh = read_ts(data_loc + "HUM_timeseries_ALL{}.txt".format(settings.YEAR), "rh", "L")
    marine_rh = read_ts(data_loc + "HUM_timeseries_ALL{}.txt".format(settings.YEAR), "rh", "M")


    COLOURS = settings.COLOURS["hydrological"]
    fig = plt.figure(figsize = (14,10))

    # manually set up the 8 axes
    w=0.42
    h=0.24
    ax1 = plt.axes([0.50-w,0.99-h,w,h])
    ax2 = plt.axes([0.50,  0.99-h,w,h])
    ax3 = plt.axes([0.50-w,0.99-(2*h),w,h],sharex=ax1)
    ax4 = plt.axes([0.50,  0.99-(2*h),w,h],sharex=ax2)
    ax5 = plt.axes([0.50-w,0.99-(3*h),w,h],sharex=ax1)
    ax6 = plt.axes([0.50,  0.99-(3*h),w,h],sharex=ax2)
    ax7 = plt.axes([0.50-w,0.99-(4*h),w,h],sharex=ax1)
    ax8 = plt.axes([0.50,  0.99-(4*h),w,h],sharex=ax2)

    utils.plot_ts_panel(ax1, land_q[:4], "-", "hydrological", loc = LEGEND_LOC, bbox = BBOX)
    utils.plot_ts_panel(ax3, [marine_q[1],marine_q[2],marine_q[3],marine_q[4]], "-", "hydrological", loc = LEGEND_LOC, bbox = BBOX)
    utils.plot_ts_panel(ax5, land_rh[:4], "-", "hydrological", loc = "")
    utils.plot_ts_panel(ax7, [marine_rh[1],marine_rh[2]], "-", "hydrological", loc = "")

    utils.plot_ts_panel(ax2, land_q[-3:], "-", "hydrological", loc = LEGEND_LOC, bbox = BBOX)
    utils.plot_ts_panel(ax4, marine_q[-3:], "-", "hydrological", loc = LEGEND_LOC, bbox = BBOX)
    utils.plot_ts_panel(ax6, land_rh[-3:], "-", "hydrological", loc = "")
    utils.plot_ts_panel(ax8, marine_rh[-3:], "-", "hydrological", loc = "")

    # prettify
    ax1.set_xlim([1958,2019])
    ax2.set_xlim([1958,2019])

    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_ylim([-0.5,0.5])
    for ax in [ax5, ax6, ax7, ax8]:
        ax.set_ylim([-1.5, 1.5])

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


    ax1.text(0.02, 0.85, "(a) In Situ Land q", transform = ax1.transAxes, fontsize = settings.FONTSIZE)
    ax2.text(0.02, 0.85, "(b) Reanalyses Land q", transform = ax2.transAxes, fontsize = settings.FONTSIZE)
    ax3.text(0.02, 0.85, "(c) In Situ & Satellite Ocean q", transform = ax3.transAxes, fontsize = settings.FONTSIZE)
    ax4.text(0.02, 0.85, "(d) Reanalysis Ocean q", transform = ax4.transAxes, fontsize = settings.FONTSIZE)
    ax5.text(0.02, 0.85, "(e) In Situ Land RH", transform = ax5.transAxes, fontsize = settings.FONTSIZE)
    ax6.text(0.02, 0.85, "(f) Reanalysis Land RH", transform = ax6.transAxes, fontsize = settings.FONTSIZE)
    ax7.text(0.02, 0.85, "(g) In Situ Ocean RH", transform = ax7.transAxes, fontsize = settings.FONTSIZE)
    ax8.text(0.02, 0.85, "(h) Reanalysis Ocean RH", transform = ax8.transAxes, fontsize = settings.FONTSIZE)


    plt.figtext(0.01, 0.75, "Specific Humidity (g kg"+r'$^{-1}$'+")", va='center', rotation='vertical', fontsize = settings.FONTSIZE)
    plt.figtext(0.01, 0.25, "Relative Humidity (%rh)", va='center', rotation='vertical', fontsize = settings.FONTSIZE)

    fig.subplots_adjust(right = 0.95, top = 0.95, bottom = 0.05, hspace = 0.001)

    plt.savefig(image_loc + "HUM_ts{}".format(settings.OUTFMT))
    plt.close()

    #*********************************************
    # Map plots

    # RH HadISDH
    cube = read_maps(data_loc + "HUMrh_anomalymap_HADISDHland{}.txt".format(settings.YEAR), "HadISDH RH", None, footer = False)

    bounds = [-20., -4, -3, -2, -1, 0, 1, 2, 3, 4, 20]

    utils.plot_smooth_map_iris(image_loc + "HUM_RH_hadisdh_land", cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (%rh)", figtext = "", title = "")

    # RH ERA 
    cube = read_maps(data_loc + "HUMrh_anomalymap_ERAI{}.txt".format(settings.YEAR), "ERA-I RH", None, footer = True)

    utils.plot_smooth_map_iris(image_loc + "HUM_RH_era", cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (%rh)", figtext = "", title = "")
    utils.plot_smooth_map_iris(image_loc + "p2.1_HUM_RH_era", cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (%rh)", figtext = "(m) Surface Relative Humidity", title = "")

    # q HadISDH
    cube = read_maps(data_loc + "HUMq_anomalymap_HADISDHland{}.txt".format(settings.YEAR), "HadISDH q", "g/kg", footer = False)

    bounds = [-20., -1.2, -0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9, 1.2, 20]

    utils.plot_smooth_map_iris(image_loc + "HUM_q_hadisdh_land", cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (g kg"+r'$^{-1}$'+")", figtext = "", title = "")
    utils.plot_smooth_map_iris(image_loc + "p2.1_HUM_q_hadisdh_land", cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (g kg"+r'$^{-1}$'+")", figtext = "(l) Surface Specific Humidity", title = "")

    # q ERA
    cube = read_maps(data_loc + "HUMq_anomalymap_ERAI{}.txt".format(settings.YEAR), "ERA-I q", "g/kg", footer = True)

    utils.plot_smooth_map_iris(image_loc + "HUM_q_era", cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (g kg"+r'$^{-1}$'+")", figtext = "", title = "")

    # MERRA?



    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
