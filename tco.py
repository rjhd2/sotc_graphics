#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for Tropospheric Ozone (TCO) section.
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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import matplotlib.cm as mpl_cm
import matplotlib as mpl

import iris
import iris.quickplot as qplt
import cartopy.crs as ccrs

import calendar
import datetime as dt

import utils # RJHD utilities
import settings

data_loc = "/data/local/rdunn/SotC/{}/data/TCO/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

LW = 3
LEGEND_LOC = "center right"


#************************************************************************
# def read_data(filename):
# from 2016 report
#     indata = np.genfromtxt(filename, dtype = (float), skip_header = 12, missing_values = "NaN")

#     year = indata[:,0]
#     month = indata[:,1]

#     times = year + (month - 1.)/12.

#     monthly_global = utils.Timeseries("mg", times, indata[:,2])
#     annual_global = utils.Timeseries("ag", times, indata[:,3])
    
#     monthly_NH = utils.Timeseries("mnh", times, indata[:,4])
#     annual_NH = utils.Timeseries("anh", times, indata[:,5])
    
#     monthly_SH = utils.Timeseries("msh", times, indata[:,6])
#     annual_SH = utils.Timeseries("ash", times, indata[:,7])

#     return monthly_global, annual_global, monthly_NH, annual_NH, monthly_SH, annual_SH #  read_data

#************************************************************************
def read_data(filename, name):

    indata = np.genfromtxt(filename, dtype = (str), missing_values = "-999.000")

    times = indata[:,0]
    data = indata[:,1].astype(float)
    data = np.ma.masked_where(data == -999.000, data)

    dt_times = [dt.datetime.strptime(d, "%b%y") for d in times]

    times = [d.year + (d.month - 1.)/12. for d in dt_times]

    timeseries = utils.Timeseries(name, np.array(times), data)

    return timeseries #  read_data
#************************************************************************
def read_map(filename, name, units):
    '''
    Read data for maps and convert to cube.  Given as single list files

    :param str filename: file to read

    :returns: cube
    '''
    
    indata = np.genfromtxt(filename, dtype = (float), skip_header = 1)

    lats = np.arange(-177.5, 177.5+5, 5) # hard coded from header
    lons = np.arange(-57.5, 57.5+5, 5) # hard coded from header
    anoms = indata[:,1:]

    cube = utils.make_iris_cube_2d(anoms.T, lons, lats, name, units)

    return cube


#************************************************************************
def run_all_plots():

    #************************************************************************
    # Global Anomaly map

    cube = read_map(data_loc + "tco_omimls_anomaly_{}.txt".format(settings.YEAR), "TCO_anom", "DU")

    bounds = np.array([-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100])
    bounds = np.array([-100, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0, 100])

    utils.plot_smooth_map_iris(image_loc + "TCO_anomaly_{}".format(settings.YEAR), cube, settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2005-{} (DU)".format(int(settings.YEAR[2:])-1), contour = True)
    utils.plot_smooth_map_iris(image_loc + "p2.1_TCO_anomaly_{}".format(settings.YEAR), cube, settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2005-{} (DU)".format(int(settings.YEAR[2:])-1), figtext = "(x) OMI/MLS Tropospheric Column Ozone", contour = True)


    #************************************************************************
    # Global Trend map
#    cube = read_map(data_loc + "".format(settings.YEAR())

#    bounds = np.array([-10, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 10])
#    utils.plot_smooth_map_iris(image_loc + "TCO_trend_{}".format(settings.YEAR), cube,  settings.COLOURMAP_DICT["composition"], bounds, "(DU per decade)")

    #************************************************************************
    # Timeseries
    # monthly_global, annual_global, monthly_NH, annual_NH, monthly_SH, annual_SH = read_data(data_loc + "OMI_MLS_trop_ozone_burden_2004_2015.txt")

    monthly_global = read_data(data_loc + "BAMS_SOTC_TROPOSPHERIC_OZONE_TG_60Sto60N.txt", "mg")
    monthly_SH = read_data(data_loc + "BAMS_SOTC_TROPOSPHERIC_OZONE_TG_0to60S.txt", "msh")
    monthly_NH = read_data(data_loc + "BAMS_SOTC_TROPOSPHERIC_OZONE_TG_0to60N.txt", "mnh")

    annual_global = read_data(data_loc + "BAMS_SOTC_TROPOSPHERIC_OZONE_TG_RUNNING_MEAN_60Sto60N.txt", "ag")
    annual_SH = read_data(data_loc + "BAMS_SOTC_TROPOSPHERIC_OZONE_TG_RUNNING_MEAN_0to60S.txt", "ash")
    annual_NH = read_data(data_loc + "BAMS_SOTC_TROPOSPHERIC_OZONE_TG_RUNNING_MEAN_0to60N.txt", "anh")

    minor_tick_interval = 1
    minorLocator = MultipleLocator(minor_tick_interval)
    COLOURS = settings.COLOURS["composition"]

    fig = plt.figure(figsize = (10,7))
    ax = plt.axes([0.13, 0.07, 0.75, 0.86])

    plt.plot(monthly_global.times, monthly_global.data,  'k', ls = '-', label = r"1.76$\pm$0.59 Tg yr$^{-1}$", lw = LW)
    plt.plot(annual_global.times, annual_global.data,  'k', ls = '--',  lw = LW)
    plt.text(2004, 250, "60"+r'$^{\circ}$'+"S - 60"+r'$^{\circ}$'+"N", va='center', color = 'k', fontsize = settings.FONTSIZE)

    plt.plot(monthly_NH.times, monthly_NH.data,  'r', ls = '-', label = r"0.93$\pm$0.49 Tg yr$^{-1}$", lw = LW)
    plt.plot(annual_NH.times, annual_NH.data,  'r', ls = '--',  lw = LW)
    plt.text(2004, 180, "0"+r'$^{\circ}$'+" - 60"+r'$^{\circ}$'+"N", va='center', color = 'r', fontsize = settings.FONTSIZE)

    plt.plot(monthly_SH.times, monthly_SH.data,  'c', ls = '-', label = r"0.83$\pm$0.58 Tg yr$^{-1}$", lw = LW)
    plt.plot(annual_SH.times, annual_SH.data,  'c', ls = '--',  lw = LW)
    plt.text(2004, 110, "60"+r'$^{\circ}$'+"S - 0"+r'$^{\circ}$'+"", va='center', color = 'c', fontsize = settings.FONTSIZE)

    ax.legend(loc = LEGEND_LOC, ncol = 1, frameon = False, prop = {'size':settings.LEGEND_FONTSIZE}, labelspacing = 0.1, columnspacing = 0.5)

    # prettify
    ax.xaxis.set_minor_locator(minorLocator)
    utils.thicken_panel_border(ax)

    fig.text(0.035, 0.5, "Tropospheric Ozone Mass\n(Tg)", va='center', rotation='vertical', ha="center", fontsize = settings.FONTSIZE)

    plt.xlim([2003.5, int(settings.YEAR)+1.5])
    plt.ylim([90, 340])
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)

    plt.savefig(image_loc + "TCO_ts{}".format(settings.OUTFMT))
    plt.close()

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
