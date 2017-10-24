#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for Tropospheric Ozone (TCO) section.
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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import matplotlib.cm as mpl_cm
import matplotlib as mpl

import iris
import iris.quickplot as qplt
import cartopy.crs as ccrs

import calendar

import utils # RJHD utilities
import settings

data_loc = "/data/local/rdunn/SotC/{}/data/TCO/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

LW = 3
LEGEND_LOC = "center right"


#************************************************************************
def read_data(filename):

    indata = np.genfromtxt(filename, dtype = (float), skip_header = 12, missing_values = "NaN")

    year = indata[:,0]
    month = indata[:,1]

    times = year + (month - 1.)/12.

    monthly_global = utils.Timeseries("mg", times, indata[:,2])
    annual_global = utils.Timeseries("ag", times, indata[:,3])
    
    monthly_NH = utils.Timeseries("mnh", times, indata[:,4])
    annual_NH = utils.Timeseries("anh", times, indata[:,5])
    
    monthly_SH = utils.Timeseries("msh", times, indata[:,6])
    annual_SH = utils.Timeseries("ash", times, indata[:,7])

    return monthly_global, annual_global, monthly_NH, annual_NH, monthly_SH, annual_SH #  read_data

#************************************************************************
def read_map(filename):
    '''
    Read data for maps and convert to cube.  Given as single list files

    :param str filename: file to read

    :returns: cube
    '''
    
    indata = np.genfromtxt(filename, delimiter = ",", dtype = (float), skip_header = 10)

    lats = indata[:,0]
    lons = indata[:,1]
    anoms = indata[:,2]

    longitudes = np.unique(lons)
    latitudes = np.unique(lats)

    data = np.ma.zeros((len(latitudes), len(longitudes)))
    data.mask = np.ones(data.shape)
    
    for v, val in enumerate(anoms):

        xloc, = np.where(longitudes == lons[v])
        yloc, = np.where(latitudes == lats[v])

        data[yloc[0], xloc[0]] = val    

    data = np.ma.masked_where(data <= -999, data)

    cube = utils.make_iris_cube_2d(data, latitudes, longitudes, "TCO_anom", "DU")

    return cube


#************************************************************************
def run_all_plots():

    #************************************************************************
    # Global Anomaly map

    cube = read_map(data_loc + "TropOzone_anomalies_{}v2_rd.txt".format(settings.YEAR))

    bounds = np.array([-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100])
    bounds = np.array([-100, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0, 100])

    utils.plot_smooth_map_iris(image_loc + "TCO_anomaly_{}".format(settings.YEAR), cube, settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2005-14 (DU)", contour = True)
    utils.plot_smooth_map_iris(image_loc + "p2.1_TCO_anomaly_{}".format(settings.YEAR), cube, settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2005-14 (DU)", figtext = "(w) OMI/MLS Tropospheric Column Ozone", contour = True)


    #************************************************************************
    # Global Actuals map
    cube = read_map(data_loc + "OMI_MLS_TCO_annual_mean_2015.txt")

    bounds = np.array([0, 12, 16, 20, 24, 28, 32, 36, 40, 50])
    print "using hard coded sequential colourmap"
    utils.plot_smooth_map_iris(image_loc + "TCO_mean_{}".format(settings.YEAR), cube, plt.cm.YlOrBr, bounds, "(DU)")

    #************************************************************************
    # Timeseries
    monthly_global, annual_global, monthly_NH, annual_NH, monthly_SH, annual_SH = read_data(data_loc + "OMI_MLS_trop_ozone_burden_2004_2015.txt")

    minor_tick_interval = 1
    minorLocator = MultipleLocator(minor_tick_interval)
    COLOURS = settings.COLOURS["composition"]

    fig = plt.figure(figsize = (10,8))
    ax = plt.axes([0.13, 0.07, 0.75, 0.86])

    plt.plot(monthly_global.times, monthly_global.data,  'k', ls = '-', label = r"2.31$\pm$0.71 Tg yr$^{-1}$", lw = LW)
    plt.plot(annual_global.times, annual_global.data,  'k', ls = '-',  lw = LW)
    plt.text(2004, 300, "60"+r'$^{\circ}$'+"S - 60"+r'$^{\circ}$'+"N", va='center', color = 'k', fontsize = settings.FONTSIZE)

    plt.plot(monthly_NH.times, monthly_NH.data,  'r', ls = '-', label = r"1.18$\pm$0.62 Tg yr$^{-1}$", lw = LW)
    plt.plot(annual_NH.times, annual_NH.data,  'r', ls = '-',  lw = LW)
    plt.text(2004, 170, "0"+r'$^{\circ}$'+" - 60"+r'$^{\circ}$'+"N", va='center', color = 'c', fontsize = settings.FONTSIZE)

    plt.plot(monthly_SH.times, monthly_SH.data,  'c', ls = '--', label = r"1.13$\pm$0.74 Tg yr$^{-1}$", lw = LW)
    plt.plot(annual_SH.times, annual_SH.data,  'c', ls = '--',  lw = LW)
    plt.text(2004, 100, "60"+r'$^{\circ}$'+"S - 0"+r'$^{\circ}$'+"", va='center', color = 'c', fontsize = settings.FONTSIZE)

    ax.legend(loc = LEGEND_LOC, ncol = 1, frameon = False, prop = {'size':settings.LEGEND_FONTSIZE}, labelspacing = 0.1, columnspacing = 0.5)

    # prettify
    ax.xaxis.set_minor_locator(minorLocator)
    utils.thicken_panel_border(ax)

    fig.text(0.03, 0.5, "Tg", va='center', rotation='vertical', fontsize = settings.FONTSIZE)

    plt.xlim([2003.5, 2017.5])
    plt.ylim([80, 330])
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
