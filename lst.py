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

import matplotlib.cm as mpl_cm
import matplotlib as mpl

import iris
import iris.quickplot as qplt
import cartopy.crs as ccrs

import struct
import calendar

import utils # RJHD utilities
import settings

data_loc = "/data/local/rdunn/SotC/{}/data/LST/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

CLIM_PERIOD="8110"
LEGEND_LOC = 'lower left'

MONTHS = [calendar.month_name[i][:3] for i in range(1,13)]

#************************************************************************
def read_csv(filename):
    """
    Read user supplied CSV for LST into Timeseries object
    """

    indata = np.genfromtxt(filename, delimiter = ',', dtype=(float), skip_header = 2, missing_values = "", filling_values = -99.9)
   
    indata = np.ma.masked_where(indata == -99.9, indata)

    raobcore = utils.Timeseries("RAOBCORE v1.5", indata[:,0], indata[:,1])
    rich = utils.Timeseries("RICH v1.5", indata[:,0], indata[:,2])
    ratpac = utils.Timeseries("RATPAC A2", indata[:,0], indata[:,3])
    unsw = utils.Timeseries("UNSW v1.0", indata[:,0], indata[:,4])
    UAH = utils.Timeseries("UAH v6.0", indata[:,0], indata[:,5])
    rss = utils.Timeseries("RSS v3.3", indata[:,0], indata[:,6])
    noaa = utils.Timeseries("NOAA v3.0", indata[:,0], indata[:,7])
    cfsr = utils.Timeseries("CFSR", indata[:,0], indata[:,8])

    era = utils.Timeseries("ERA-Interim", indata[:,0], indata[:,10])

    merra = utils.Timeseries("MERRA-2", indata[:,0], indata[:,9])
    jra = utils.Timeseries("JRA-55", indata[:,0], indata[:,11])

    return UAH, rss, ratpac, raobcore, rich, noaa, cfsr, era, unsw # read_csv

#************************************************************************
def read_merra_monthly(filename):
    """
    Read user supplied CSV for LST into Timeseries object
    """
    indata = np.genfromtxt(filename, dtype=(float), skip_header = 1)

    # process the years and months to give decimals
    years = indata[:,0]
    months = indata[:,1]

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

    return fourh, threeh, twofiveh, twoh, onefiveh, fourtwoav, oneh, seventy, fifty, thirty, twenty, ten, seventytwentyav # read_merra_monthly

#************************************************************************
def read_zonal(filename, platform):
    """
    Read user supplied CSV for LST into Timeseries object
    """
    indata = np.genfromtxt(filename, dtype=(float), skip_header = 1)

    print "15/6/17 - error in latitudes, needing to be flipped - Craig Long"
    if platform == "R":
        cfsr = utils.Timeseries("CFSR", indata[:,1], -indata[:,0])
        merra = utils.Timeseries("MERRA-2", indata[:,2], -indata[:,0])
        era = utils.Timeseries("ERA-Interim", indata[:,3], -indata[:,0])
        jra= utils.Timeseries("JRA-55", indata[:,4], -indata[:,0])
        
        return cfsr, merra, era, jra

    if platform == "S":
        star = utils.Timeseries("NOAA v3.0", indata[:,1], -indata[:,0])
        uah = utils.Timeseries("UAH v6.0", indata[:,2], -indata[:,0])
        rss = utils.Timeseries("RSS v3.3", indata[:,3], -indata[:,0])

        return star, uah, rss # read_zonal
        


#************************************************************************
def run_all_plots():

    #************************************************************************
    # Timeseries figures

    fig, (ax1, ax2, ax3) = plt.subplots(3, figsize = (10,12), sharex=True)

    UAH, rss, ratpac, raobcore, rich, noaa, cfsr, era, unsw = read_csv(data_loc + "TLS_annual_anomalies_{}.csv".format(settings.YEAR))

    # Sondes
    utils.plot_ts_panel(ax1, [raobcore, rich, ratpac, unsw], "-", "temperature", loc = LEGEND_LOC)

    # satellites
    utils.plot_ts_panel(ax2, [UAH, noaa, rss], "-", "temperature", loc = LEGEND_LOC)

    ax2.set_ylabel("Anomaly ("+r'$^\circ$'+"C)", fontsize = settings.FONTSIZE)

    # reanalyses
    jra_actuals, jra_anoms = utils.read_jra55(reanalysis_loc + "JRA-55_MSUch4_global_ts.txt", "temperature")
    merra_actuals, merra_anoms = utils.read_merra_LT_LS(reanalysis_loc + "MERRA2_MSU_Tanom_ann.dat", LS = True)
    utils.plot_ts_panel(ax3, [era, jra_anoms, merra_anoms, cfsr], "-", "temperature", loc = LEGEND_LOC)

    # sort formatting
    plt.xlim([1957,era.times[-1]+1])

    for tick in ax3.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 
    for ax in [ax1, ax2, ax3]:
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 

        ax.set_ylim([-1.0,2.4])

    # sort labelling
    ax1.text(0.02, 0.9, "(a) Radiosondes", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE)
    ax2.text(0.02, 0.9, "(b) Satellites", transform = ax2.transAxes, fontsize = settings.LABEL_FONTSIZE)
    ax3.text(0.02, 0.9, "(c) Reanalyses", transform = ax3.transAxes, fontsize = settings.LABEL_FONTSIZE)

    fig.subplots_adjust(right = 0.95, top = 0.95, hspace = 0.001)

    plt.savefig(image_loc+"LST_ts{}".format(settings.OUTFMT))

    plt.close()

    #************************************************************************
    # MERRA Trop and Strat Timeseries figures

    fourh, threeh, twofiveh, twoh, onefiveh, fourtwoav, oneh, seventy, fifty, thirty, twenty, ten, seventytwentyav = read_merra_monthly(data_loc + "MERRA2_400_20_temp_anom_1980-{}.txt".format(settings.YEAR))

    fourtwoav.data = np.ma.masked_where(fourtwoav.times < 1994, fourtwoav.data)
    seventytwentyav.data = np.ma.masked_where(seventytwentyav.times < 1994, seventytwentyav.data)
    fourtwoav.times = np.ma.masked_where(fourtwoav.times < 1994, fourtwoav.times)
    seventytwentyav.times = np.ma.masked_where(seventytwentyav.times < 1994, seventytwentyav.times)

    
    fig, (ax1, ax2) = plt.subplots(2, figsize = (10,10), sharex=True)

    # trop
    fit = utils.fit_plot_points(0.025, -50.01, fourtwoav.times)
    ax1.plot(fourh.times, fit, c="0.5", lw = 2, ls = "-", label = "fit", zorder = 10)
    utils.plot_ts_panel(ax1, [fourh, threeh, twofiveh, twoh, fourtwoav], "-", "lst", loc = LEGEND_LOC)
    # strat
    fit = utils.fit_plot_points(-0.0014, 2.4527, seventytwentyav.times)
    ax2.plot(fourh.times, fit, c="0.5", lw = 2, ls = "-", label = "fit", zorder = 10)
    utils.plot_ts_panel(ax2, [seventy, fifty, thirty, twenty, seventytwentyav], "-", "lst", loc = LEGEND_LOC)

    # sort formatting
    plt.xlim([1978, fourh.times[-1]+1])
    ax1.set_ylim([-1.7,1.5])
    ax2.set_ylim([-2.0,2.7])

    for tick in ax2.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 
    for ax in [ax1, ax2]:
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 

    # sort labelling
    ax1.text(0.02, 0.9, "(a) MERRA-2 Troposphere", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE)
    ax2.text(0.02, 0.9, "(b) MERRA-2 Stratosphere", transform = ax2.transAxes, fontsize = settings.LABEL_FONTSIZE)

    ax1.text(0.7, 0.9, "y=0.025x - 50.01", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE*0.8)
    ax2.text(0.7, 0.9, "y=-0.0014x + 2.4527", transform = ax2.transAxes, fontsize = settings.LABEL_FONTSIZE*0.8)

    ax1.set_ylabel("Anomaly ("+r'$^\circ$'+"C)", fontsize = settings.FONTSIZE)
    ax2.set_ylabel("Anomaly ("+r'$^\circ$'+"C)", fontsize = settings.FONTSIZE)

    fig.subplots_adjust(right = 0.95, top = 0.95, hspace = 0.001)

    plt.savefig(image_loc+"LST_merra_ts{}".format(settings.OUTFMT))

    plt.close()

    #************************************************************************
    # Zonal figures

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (10,8), sharey=True)

    cfsr, merra, era, jra = read_zonal(data_loc + "TLS_Reanal_zonal_trends_1994-{}.txt".format(settings.YEAR), "R")
    star, uah, rss = read_zonal(data_loc + "TLS_Satellite_zonal_trends_1994-{}.txt".format(settings.YEAR), "S")

    # Reanalyses
    utils.plot_ts_panel(ax1, [cfsr, merra, era, jra], "--", "temperature", loc = "lower left", ncol = 1)

    # Satellite
    utils.plot_ts_panel(ax2, [star, uah, rss], "--", "temperature", loc = "center right", ncol = 1)

    ax1.axvline(0,color = "0.5", ls = "--")
    ax2.axvline(0,color = "0.5", ls = "--")

    # sort formatting
    plt.ylim([-90,90])
    ax1.set_ylabel("Latitude", fontsize = settings.LABEL_FONTSIZE)
    ax1.set_xlabel("Trend ("+r'$^\circ$'+"C decade"+r'$^{-1}$'+")", fontsize = settings.FONTSIZE)
    ax2.set_xlabel("Trend ("+r'$^\circ$'+"C decade"+r'$^{-1}$'+")", fontsize = settings.FONTSIZE)


    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 
    for ax in [ax1, ax2]:
        ax.set_xlim([-0.3,0.6])
        ax.set_xticks(np.arange(-0.2,0.8,0.2))
        ax.set_yticks(np.arange(-90,120,30))
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 

    # sort labelling
    ax1.text(0.02, 0.9, "(a) Radiosondes", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE)
    ax2.text(0.02, 0.9, "(b) Satellites", transform = ax2.transAxes, fontsize = settings.LABEL_FONTSIZE)

    fig.subplots_adjust(right = 0.95, top = 0.95, hspace = 0.001)

    plt.savefig(image_loc+"LST_profiles{}".format(settings.OUTFMT))

    plt.close()


    #************************************************************************
    # ERA Anomaly figure

    # Read in ERA anomalies

    cube_list = iris.load(reanalysis_loc + "TLS_anvj_moda_ann{}{}-ann19812010.nc".format(settings.YEAR, settings.YEAR))

    cube = cube_list[0]
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()

    bounds=[-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]

    utils.plot_smooth_map_iris(image_loc + "LST_{}_anoms_era".format(settings.YEAR), cube[0][0], settings.COLOURMAP_DICT["temperature"], bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)")
    utils.plot_smooth_map_iris(image_loc + "p2.1_LST_{}_anoms_era".format(settings.YEAR), cube[0][0], settings.COLOURMAP_DICT["temperature"], bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", figtext = "(f) Lower Stratospheric Temperature")


    #************************************************************************
    # MERRA Anomaly figure
    import netCDF4 as ncdf

    ncfile = ncdf.Dataset(reanalysis_loc + "merra2_tls_ANNUAL_anom.nc")

    var=ncfile.variables["TLS ANOM"][:] # this is a masked array
    nlons = ncfile.variables["LONGITUDES"][:]
    nlats = ncfile.variables["LATITUDES"][:]

    cube = utils.make_iris_cube_2d(var, nlats, nlons, "LST", "C")

    utils.plot_smooth_map_iris(image_loc + "LST_{}_anoms_merra".format(settings.YEAR), cube, settings.COLOURMAP_DICT["temperature"], bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)")


    #************************************************************************
    # 2015 MERRA seasonal figure

    import netCDF4 as ncdf

    month_list = []

    for month in MONTHS:
        print month

        # IRIS doesn't like the whitespace in the "TLS ANOM" 
        ncfile=ncdf.Dataset(reanalysis_loc + "merra2_tls_{}_anom.nc".format(month.upper()),'r')

        var=ncfile.variables["TLS ANOM"][:] # this is a masked array
        lons = ncfile.variables["LONGITUDES"][:]
        lats = ncfile.variables["LATITUDES"][:]

        ncfile.close()

        cube = utils.make_iris_cube_2d(var, lats, lons, "TLS_ANOM", "C")

        month_list += [cube]

    # pass to plotting routine
    utils.plot_smooth_map_iris_multipanel(image_loc + "LST_{}_monthly_merra".format(settings.YEAR), month_list, settings.COLOURMAP_DICT["temperature"], bounds, "Anomaly ("+r'$^{\circ}$'+"C)", shape = (6,2), title = MONTHS, figtext = ["(a)","(b)","(c)","(d)", "(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"])

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
