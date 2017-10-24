#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for Cloudiness (CLD) section.
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

import iris
import iris.quickplot as qplt
import cartopy.crs as ccrs

import scipy.io

import utils # RJHD utilities
import settings


data_loc = "/data/local/rdunn/SotC/{}/data/CLD/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

LEGEND_LOC = 'upper left'
BBOX = (0.2,1.0)

CLIMSTART = 2003
CLIMEND = 2015

#************************************************************************
def read_ts(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read

    :returns: Timeseries object s
    '''

    raw_data = np.genfromtxt(filename, dtype = (float), skip_header = 1)
    
    raw_data = np.ma.masked_where(raw_data <= -999.0, raw_data)    

    times = raw_data[:,0]

    cstart, = np.where(times == CLIMSTART)
    cend, = np.where(times == CLIMEND)

    clims = np.ma.mean(raw_data[cstart:cend], axis = 0)

    data = raw_data - clims


    patmosx = utils.Timeseries("PATMOS-x/AVHRR", times, data[:,1])
    hirs = utils.Timeseries("HIRS", times, data[:,2])
    misr = utils.Timeseries("MISR", times, data[:,3])
    modis = utils.Timeseries("PATMOS-x/AQUA MODIS C6", times, data[:,4])
    calipso = utils.Timeseries("CALIPSO", times, data[:,5])
    ceres = utils.Timeseries("CERES", times, data[:,6])
    satcorps = utils.Timeseries("SatCORPS", times, data[:,7])
    clara_a2 = utils.Timeseries("CLARA-A2", times, data[:,8])

    return patmosx, hirs, misr, modis, calipso, ceres, satcorps, clara_a2 # read_ts


#************************************************************************
def run_all_plots():

    #************************************************************************
    # Cloudiness timeseries

    patmosx, hirs, misr, modis, calipso, ceres, satcorps, clara_a2 = read_ts(data_loc + "{}_global_cloudiness_timeseries.txt".format(settings.YEAR))

    fig = plt.figure(figsize = (12,6))

    ax1 = plt.axes([0.1,0.1,0.88,0.88])

    utils.plot_ts_panel(ax1, [patmosx, hirs, misr, modis, calipso, ceres, satcorps, clara_a2], "-", "hydrological", loc = LEGEND_LOC, bbox = BBOX, ncol = 3)

    ax1.text(0.02, 0.9, "Satellite", transform = ax1.transAxes, fontsize = settings.FONTSIZE)

    #*******************
    # prettify
    ax1.set_ylim([-4.4,6.9])
    ax1.set_xlim([hirs.times[0]-1,hirs.times[-1]+1])
    ax1.set_ylabel("Anomaly (%)", fontsize = settings.FONTSIZE)

    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 

    plt.savefig(image_loc+"CLD_ts{}".format(settings.OUTFMT))
    plt.close()

    #************************************************************************
    # Cloudiness map

    mapfile_dict = scipy.io.readsav(data_loc + "patmosx_global_monthly_cloudiness_anomaly_map_{}.sav".format(settings.YEAR))

    annual_anoms = mapfile_dict["annual_anom"]*100.
    djf_anoms = mapfile_dict["djf_anom"]*100.
    jja_anoms = mapfile_dict["jja_anom"]*100.
    mam_anoms = mapfile_dict["mam_anom"]*100.
    son_anoms = mapfile_dict["son_anom"]*100.

    lats = mapfile_dict["lat"]
    lons = mapfile_dict["lon"]

    cube = utils.make_iris_cube_2d(annual_anoms, lats[:,0], lons[0]  "CLD_anom", "%")

    bounds=[-100, -15, -10, -5, -2.5, 0, 2.5, 5, 10, 15, 100]

    utils.plot_smooth_map_iris(image_loc + "CLD_{}_anoms".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (%)")
    utils.plot_smooth_map_iris(image_loc + "p2.1_CLD_{}_anoms".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (%)", figtext = "(n) Cloudiness")


    #************************************************************************
    # Cloudiness Seasonal Map

    cubelist = []
    for season in [djf_anoms, mam_anoms, jja_anoms, son_anoms]:

        cube = utils.make_iris_cube_2d(season, lats[:,0], lons[0], "CLD_anom", "%")
        cubelist += [cube]


    utils.plot_smooth_map_iris_multipanel(image_loc + "CLD_{}_anoms_seasons".format(settings.YEAR), cubelist, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomaly (%)", shape = (2,2), title = ["DJF", "MAM", "JJA", "SON"], figtext = ["(a)","(b)", "(c)", "(d)"])


    #************************************************************************
    # Cloudiness Hovmoller

    data_dict = scipy.io.readsav(data_loc + "patmosx_global_monthly_cloudiness_hovmuller_{}.sav".format(settings.YEAR))

    lats = data_dict["latitude"]
    times = data_dict["hov_time"]
    anoms = data_dict["hov_anom"]*100.

    utils.plot_hovmuller(image_loc + "CLD_hovmuller", times, lats, anoms, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomaly (%)")

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
