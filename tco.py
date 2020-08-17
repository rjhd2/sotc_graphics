#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for Tropospheric Ozone (TCO) section.
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

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import datetime as dt
import numpy as np

import utils # RJHD utilities
import settings

DATALOC = "{}/{}/data/TCO/".format(settings.ROOTLOC, settings.YEAR)

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

    indata = np.genfromtxt(filename, dtype=(str), missing_values="-999.000")

    times = indata[:, 0]
    data = indata[:, 1].astype(float)
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
    
    indata = np.genfromtxt(filename, dtype=(float), skip_header=1)

    lons = np.arange(-177.5, 177.5+5, 5) # hard coded from header
    lats = np.arange(-57.5, 57.5+5, 5) # hard coded from header
    anoms = indata[:, 1:]

    cube = utils.make_iris_cube_2d(anoms.T, lats, lons, name, units)

    return cube

#************************************************************************
def read_significance(filename, signame):
    '''
    Read data for significance and convert to lat/lon lists.  Given as single list files

    :param str filename: datafile to read
    :param str signame: significance file to read

    :returns: lats lons data
    '''
    
    # read in both data (significance is 1/0)
    indata = np.genfromtxt(filename, dtype=(float), skip_header=1)
    sigdata = np.genfromtxt(signame, dtype=(float), skip_header=1)

    lons = np.arange(-177.5, 177.5+5, 5) # hard coded from header
    lats = np.arange(-57.5, 57.5+5, 5) # hard coded from header
    anoms = indata[:, 1:]
    sigs = sigdata[:, 1:]

    # blank lists to store
    sig_lats = []
    sig_lons = []
    sig_data = []
    
    # for each longitude line
    for s, sig in enumerate(sigs):
        
        locs, = np.where(sig == 1.)

        # store matching values
        sig_lons += [lons[s] for l in locs]
        sig_lats += [lats[l] for l in locs]
        sig_data += [anoms[s][l] for l in locs]

    return np.array(sig_lats), np.array(sig_lons), np.array(sig_data) # read_significance

#************************************************************************
def plot_smooth_map_iris(outname, cube, cmap, bounds, cb_label, scatter=[], \
                             figtext="", title="", contour=False, cb_extra="", save_netcdf_filename=""):
    '''
    Standard map - 

    :param str outname: output filename root
    :param array cube: cube to plot
    :param obj cmap: colourmap to use
    :param array bounds: bounds for discrete colormap
    :param str cb_label: colorbar label
    :param str save_netcdf_filename: filename to save the output plot to a cube.
    '''

    import matplotlib as mpl
    import cartopy
    import copy
    import iris

    norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)

    fig = plt.figure(figsize=(8, 5.5))

    plt.clf()
    ax = plt.axes([0.01, 0.12, 0.98, 0.88], projection=cartopy.crs.Robinson())
    ax.gridlines() #draw_labels=True)
    ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
    ax.coastlines()

    ext = ax.get_extent() # save the original extent

    if settings.OUTFMT in [".eps", ".pdf"]:
        if cube.coord("latitude").points.shape[0] > 180 or cube.coord("longitude").points.shape[0] > 360:
            regrid_size = 1.0
            print("Regridding cube for {} output to {} degree resolution".format(settings.OUTFMT, regrid_size))
            print("Old Shape {}".format(cube.data.shape))
            plot_cube = regrid_cube(cube, regrid_size, regrid_size)
            print("New Shape {}".format(plot_cube.data.shape))
        else:
            plot_cube = copy.deepcopy(cube)
    else:
        plot_cube = copy.deepcopy(cube)

    mesh = iris.plot.pcolormesh(plot_cube, cmap=cmap, norm=norm)

    if len(scatter) > 0:
        lons, lats, data = scatter
        plt.scatter(lons, lats, c=data, cmap=cmap, norm=norm, s=10, \
                        transform=cartopy.crs.Geodetic(), edgecolor='0.2', linewidth=1.0)


    cb = plt.colorbar(mesh, orientation='horizontal', pad=0.05, fraction=0.05, \
                        aspect=30, ticks=bounds[1:-1], drawedges=True)
    # thicken border of colorbar and the dividers
    # http://stackoverflow.com/questions/14477696/customizing-colorbar-border-color-on-matplotlib
    cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
    cb.ax.tick_params(axis='x', labelsize=settings.FONTSIZE, direction='in', size=0)
    cb.set_label(label=cb_label, fontsize=settings.FONTSIZE)

#    cb.outline.set_color('k')
    cb.outline.set_linewidth(2)
    cb.dividers.set_color('k')
    cb.dividers.set_linewidth(2)

    if cb_extra != "":
        fig.text(0.04, 0.05, cb_extra[0], fontsize=settings.FONTSIZE * 0.8, ha="left")
        fig.text(0.96, 0.05, cb_extra[1], fontsize=settings.FONTSIZE * 0.8, ha="right")

    ax.set_extent(ext, ax.projection) # fix the extent change from colormesh

    plt.title(title)
    fig.text(0.03, 0.95, figtext, fontsize=settings.FONTSIZE * 0.8)

    plt.savefig(outname + settings.OUTFMT)
    plt.close()

    return # plot_smooth_map_iris


#************************************************************************
def run_all_plots():

    #************************************************************************
    # Global Anomaly map

    cube = read_map(DATALOC + "tco_omimls_anomaly_{}.txt".format(settings.YEAR), "TCO_anom", "DU")

    bounds = np.array([-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100])
    bounds = np.array([-100, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0, 100])

    utils.plot_smooth_map_iris(settings.IMAGELOC + "TCO_anomaly_{}".format(settings.YEAR), cube, settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2004-08 (DU)", contour=True)
    utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_TCO_anomaly_{}".format(settings.YEAR), cube, settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2004-08 (DU)", figtext="(aa) OMI/MLS Tropospheric Column Ozone", contour=True)


    #************************************************************************
    # Global Trend map
    cube = read_map(DATALOC + "tco_omimls_trends_{}.txt".format(settings.YEAR), "TCO_trend", "DU")

    bounds = np.array([-100, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 100])
#    utils.plot_smooth_map_iris(settings.IMAGELOC + "TCO_trend_{}".format(settings.YEAR), cube, settings.COLOURMAP_DICT["composition"], bounds, "(DU per decade)")

    sig_lats, sig_lons, sig_data = read_significance(DATALOC + "tco_omimls_trends_{}.txt".format(settings.YEAR), DATALOC + "tco_omimls_trends_95pct_signficance_{}.txt".format(settings.YEAR))
    # use local adapted routine.
    plot_smooth_map_iris(settings.IMAGELOC + "TCO_trend_significance_{}".format(settings.YEAR), cube, settings.COLOURMAP_DICT["composition"], bounds, "(DU per decade)", scatter = (sig_lons, sig_lats, sig_data))


    #************************************************************************
    # Timeseries
    # monthly_global, annual_global, monthly_NH, annual_NH, monthly_SH, annual_SH = read_data(DATALOC + "OMI_MLS_trop_ozone_burden_2004_2015.txt")

    monthly_global = read_data(DATALOC + "BAMS_SOTC_TROPOSPHERIC_OZONE_TG_60Sto60N_{}.txt".format(settings.YEAR), "mg")
    monthly_SH = read_data(DATALOC + "BAMS_SOTC_TROPOSPHERIC_OZONE_TG_0to60S_{}.txt".format(settings.YEAR), "msh")
    monthly_NH = read_data(DATALOC + "BAMS_SOTC_TROPOSPHERIC_OZONE_TG_0to60N_{}.txt".format(settings.YEAR), "mnh")

    annual_global = read_data(DATALOC + "BAMS_SOTC_TROPOSPHERIC_OZONE_TG_RUNNING_MEAN_60Sto60N_{}.txt".format(settings.YEAR), "ag")
    annual_SH = read_data(DATALOC + "BAMS_SOTC_TROPOSPHERIC_OZONE_TG_RUNNING_MEAN_0to60S_{}.txt".format(settings.YEAR), "ash")
    annual_NH = read_data(DATALOC + "BAMS_SOTC_TROPOSPHERIC_OZONE_TG_RUNNING_MEAN_0to60N_{}.txt".format(settings.YEAR), "anh")

    minor_tick_interval = 1
    minorLocator = MultipleLocator(minor_tick_interval)
    COLOURS = settings.COLOURS["composition"]

    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes([0.14, 0.07, 0.84, 0.90])

    plt.plot(monthly_global.times, monthly_global.data, 'k', ls='-', label=r"1.79$\pm$0.47 Tg yr$^{-1}$", lw=LW)
    plt.plot(annual_global.times, annual_global.data, 'k', ls='--', lw=LW)
    plt.text(2004, 250, "60"+r'$^{\circ}$'+"S - 60"+r'$^{\circ}$'+"N", va='center', color='k', fontsize=settings.FONTSIZE)

    plt.plot(monthly_NH.times, monthly_NH.data, 'r', ls='-', label=r"0.89$\pm$0.39 Tg yr$^{-1}$", lw=LW)
    plt.plot(annual_NH.times, annual_NH.data, 'r', ls='--', lw=LW)
    plt.text(2004, 180, "0"+r'$^{\circ}$'+" - 60"+r'$^{\circ}$'+"N", va='center', color='r', fontsize=settings.FONTSIZE)

    plt.plot(monthly_SH.times, monthly_SH.data, 'c', ls='-', label=r"0.90$\pm$0.46 Tg yr$^{-1}$", lw=LW)
    plt.plot(annual_SH.times, annual_SH.data, 'c', ls='--', lw=LW)
    plt.text(2004, 110, "60"+r'$^{\circ}$'+"S - 0"+r'$^{\circ}$'+"", va='center', color='c', fontsize=settings.FONTSIZE)

    ax.legend(loc=LEGEND_LOC, ncol=1, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, labelspacing=0.1, columnspacing=0.5)

    # prettify
    ax.xaxis.set_minor_locator(minorLocator)
    utils.thicken_panel_border(ax)

    fig.text(0.04, 0.5, "Tropospheric Ozone\n(Tg)", va='center', rotation='vertical', ha="center", fontsize=settings.FONTSIZE)

    plt.xlim([2003.5, int(settings.YEAR)+1.5])
    plt.ylim([90, 340])
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)

    plt.savefig(settings.IMAGELOC + "TCO_ts{}".format(settings.OUTFMT))
    plt.close()

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
