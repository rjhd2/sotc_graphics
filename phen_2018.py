#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for Phenology (PHEN) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 21                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2017-12-22 11:57:17 +0000 (Fri, 22 Dec #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import os
import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
import matplotlib.path as mpath
from matplotlib.ticker import MultipleLocator
import matplotlib.image as mpimg

import iris
import cartopy

import utils # RJHD utilities
import settings


data_loc = "{}/{}/data/PHEN/".format(settings.ROOTLOC, settings.YEAR)
reanalysis_loc = "{}/{}/data/RNL/".format(settings.ROOTLOC, settings.YEAR)
image_loc = "{}/{}/images/".format(settings.ROOTLOC, settings.YEAR)

LEGEND_LOC = 'lower left'



#************************************************************************
def read_modis_ts(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read
    :returns: Timeseries object s
    '''

    raw_data = np.genfromtxt(filename, dtype=(float), skip_header=1, delimiter=",")

    years = raw_data[:, 0]  

   
    # 2018 entries
    sos_nh = utils.Timeseries("SOS", years, raw_data[:, 1])
    sos_na = utils.Timeseries("SOS", years, raw_data[:, 2])
    sos_ea = utils.Timeseries("SOS", years, raw_data[:, 3])

    eos_nh = utils.Timeseries("EOS", years, raw_data[:, 5])
    eos_na = utils.Timeseries("EOS", years, raw_data[:, 6])
    eos_ea = utils.Timeseries("EOS", years, raw_data[:, 7])

    sprt_nh = utils.Timeseries("Spring T", years, raw_data[:, 9])
    sprt_na = utils.Timeseries("Spring T", years, raw_data[:, 10])
    sprt_ea = utils.Timeseries("Spring T", years, raw_data[:, 11])

    falt_nh = utils.Timeseries("Autumn T", years, raw_data[:, 13])
    falt_na = utils.Timeseries("Autumn T", years, raw_data[:, 14])
    falt_ea = utils.Timeseries("Autumn T", years, raw_data[:, 15])

    return sos_na, sos_ea, sprt_na, sprt_ea # read_modis_ts

#************************************************************************
def read_us_phenocam(filename):
    
    raw_data = np.genfromtxt(filename, dtype=(str), skip_header=1)

    lat = raw_data[:, 1].astype(float)
    lon = raw_data[:, 2].astype(float)

    return lat, lon # read_us_phenocam


#************************************************************************
def plot_modis_ts(axl, sos, sprt, dummy, label, anomalies, legend_loc):

    utils.plot_ts_panel(axl, [sos, dummy], "-", "phenological", loc=legend_loc)   


    # make twin
    axr = axl.twinx()
    utils.plot_ts_panel(axr, [sprt], "-", "phenological", loc="")   

    # prettify
    axl.set_ylim([-10, 10])
    axr.set_ylim([3, -3])

    # labels
    axl.text(0.02, 0.83, label, transform=axl.transAxes, fontsize=settings.FONTSIZE*0.8)
    axl.text(0.47, 0.88, anomalies[0], transform=axl.transAxes)
    axl.text(0.47, 0.78, anomalies[1], transform=axl.transAxes)

    # ticks etc
    minorLocator = MultipleLocator(1)
    for ax in [axl]:
        utils.thicken_panel_border(ax)
        ax.set_yticks(ax.get_yticks()[1:-1])
        ax.xaxis.set_minor_locator(minorLocator)

    for ax in [axr]:
        ax.yaxis.tick_right()
        utils.thicken_panel_border(ax)
        ax.set_yticks(ax.get_yticks()[1:-1])
        ax.xaxis.set_minor_locator(minorLocator)

    axl.set_xlim([1999, 2020])

    return # plot_modis_ts

#************************************************************************
def read_uk_oak_csv(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read
    :returns: Timeseries object s
    '''

    raw_data = np.genfromtxt(filename, dtype=(str), skip_header=2, delimiter=",")
    
    indata = raw_data[:, 1].astype(float)
    indata = np.ma.masked_where(indata == -99, indata)
    
    oak = utils.Timeseries("Quercus robur", raw_data[:, 0].astype(int), indata)
    
    return oak  # read_uk_oak_csv

#************************************************************************
def read_windermere_csv(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read
    :returns: Timeseries objects
    '''

    raw_data = np.genfromtxt(filename, dtype=(int), skip_header=1, delimiter=",")
    
    times = raw_data[:, 0]
    north = utils.Timeseries("North Basin", times, raw_data[:, 1])
    south = utils.Timeseries("South Basin", times, raw_data[:, 2])

    return north, south  # read_windermere_csv

#************************************************************************
def read_us_phenocam_csv(filename):


    raw_data = np.genfromtxt(filename, dtype=(str), skip_header=1, delimiter=",")

    # remove "" or NA
    locs = np.where(raw_data == "")
    raw_data[locs] = "-99.9"
    locs = np.where(raw_data == "NA")
    raw_data[locs] = "-99.9"    

    barrow_17 = utils.Timeseries("2017", raw_data[:, 2].astype(int), raw_data[:, 3].astype(float))
    barrow_17_sm = utils.Timeseries("2017", raw_data[:, 2].astype(int), raw_data[:, 4].astype(float))
    barrow_18 = utils.Timeseries("2018", raw_data[:, 2].astype(int), raw_data[:, 6].astype(float))
    barrow_18_sm = utils.Timeseries("2018", raw_data[:, 2].astype(int), raw_data[:, 7].astype(float))

    ozarks_17 = utils.Timeseries("2017", raw_data[:, 12].astype(int), raw_data[:, 13].astype(float))
    ozarks_17_sm = utils.Timeseries("2017", raw_data[:, 12].astype(int), raw_data[:, 14].astype(float))
    ozarks_18 = utils.Timeseries("2018", raw_data[:, 12].astype(int), raw_data[:, 16].astype(float))
    ozarks_18_sm = utils.Timeseries("2018", raw_data[:, 12].astype(int), raw_data[:, 17].astype(float))

    turkey_17 = utils.Timeseries("2017", raw_data[:, 21].astype(int), raw_data[:, 22].astype(float))
    turkey_17_sm = utils.Timeseries("2017", raw_data[:, 21].astype(int), raw_data[:, 23].astype(float))
    turkey_18 = utils.Timeseries("2018", raw_data[:, 21].astype(int), raw_data[:, 25].astype(float))
    turkey_18_sm = utils.Timeseries("2018", raw_data[:, 21].astype(int), raw_data[:, 26].astype(float))

    barrow_17.data = np.ma.masked_where(barrow_17.data == -99.9, barrow_17.data)
    barrow_17_sm.data = np.ma.masked_where(barrow_17_sm.data == -99.9, barrow_17_sm.data)
    barrow_18.data = np.ma.masked_where(barrow_18.data == -99.9, barrow_18.data)
    barrow_18_sm.data = np.ma.masked_where(barrow_18_sm.data == -99.9, barrow_18_sm.data)

    ozarks_17.data = np.ma.masked_where(ozarks_17.data == -99.9, ozarks_17.data)
    ozarks_17_sm.data = np.ma.masked_where(ozarks_17_sm.data == -99.9, ozarks_17_sm.data)
    ozarks_18.data = np.ma.masked_where(ozarks_18.data == -99.9, ozarks_18.data)
    ozarks_18_sm.data = np.ma.masked_where(ozarks_18_sm.data == -99.9, ozarks_18_sm.data)
 
    turkey_17.data = np.ma.masked_where(turkey_17.data == -99.9, turkey_17.data)
    turkey_17_sm.data = np.ma.masked_where(turkey_17_sm.data == -99.9, turkey_17_sm.data)
    turkey_18.data = np.ma.masked_where(turkey_18.data == -99.9, turkey_18.data)
    turkey_18_sm.data = np.ma.masked_where(turkey_18_sm.data == -99.9, turkey_18_sm.data)


    barrow = [barrow_17, barrow_18, barrow_17_sm, barrow_18_sm]
    ozarks = [ozarks_17, ozarks_18, ozarks_17_sm, ozarks_18_sm]
    turkey = [turkey_17, turkey_18, turkey_17_sm, turkey_18_sm]

    return barrow, ozarks, turkey # read_us_phenocam

#************************************************************************
def plot_us_phenocam(ax, data, xlim, label, ylabel = "", legend = False):

    from matplotlib.ticker import MultipleLocator

    minor_tick_interval=10
    minorLocator = MultipleLocator(minor_tick_interval)
    
    ax.plot(data[0].times, data[0].data, c="peru", marker="o", ls=None)
    ax.plot(data[1].times, data[1].data, c="teal", marker="o", ls=None)

    ax.plot(data[2].times, data[2].data, c="peru", lw=2, ls="-", label=data[2].name)
    ax.plot(data[3].times, data[3].data, c="teal", lw=2, ls="-", label=data[3].name)

    if legend:
        ax.legend(loc="lower right", frameon=False, ncol=1, fontsize=settings.FONTSIZE * 0.6)

    ax.xaxis.set_minor_locator(minorLocator)
    # turn of RHS y ticks
    ax.yaxis.set_ticks_position('left')
    ax.set_xlim(xlim)
    ax.set_ylim([0.28, 0.44])
    ax.set_ylabel(ylabel)
    ax.text(0.1, 0.8, label)

    utils.thicken_panel_border(ax)

    return # plot_us_phenocam

#************************************************************************
def plot_images(ax, filename):

    img = mpimg.imread(os.path.join(data_loc, filename))
    ax.imshow(img)
    ax.set_xticks([])
    ax.set_yticks([])

    return # plot_images

#************************************************************************
def run_all_plots():


    #***********************
    # MODIS - centre

    cubelist = iris.load(os.path.join(data_loc, "MODIS.CMG.{}.SOS.EOS.Anomaly.nc".format(settings.YEAR)))

    for c, cube in enumerate(cubelist):
        if cube.name() == "SOS":
            sos_cube = cubelist[c]
       
    # deal with NANS
    sos_cube.data = np.ma.masked_where(sos_cube.data != sos_cube.data, sos_cube.data)

    # read in sites

    us_locations = read_us_phenocam(os.path.join(data_loc, "phenocam_locs.txt"))   
            
    fig = plt.figure(figsize=(8, 12))
    plt.clf()

    # set up plot settings
    BOUNDS = [-100, -20, -10, -5, -2, 0, 2, 5, 10, 20, 100]

    LABELS = "(c) Start of Season (SOS)" #, "(b) End of Season (EOS)"]

    # boundary circle
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    # axes for polar plot
    ax = plt.axes([0.05, 0.01, 0.9, 0.65], projection=cartopy.crs.NorthPolarStereo(central_longitude=300.0))

    plot_cube = sos_cube

    if settings.OUTFMT in [".eps", ".pdf"]:
        if plot_cube.coord("latitude").points.shape[0] > 90 or plot_cube.coord("longitude").points.shape[0] > 360:
            regrid_size = 1.0
            print("Regridding cube for {} output to {} degree resolution".format(settings.OUTFMT, regrid_size))
            print("Old Shape {}".format(plot_cube.data.shape))
            plot_cube = utils.regrid_cube(plot_cube, regrid_size, regrid_size)
            print("New Shape {}".format(plot_cube.data.shape))

    ax.gridlines() #draw_labels=True)
    ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
    ax.coastlines()
    ax.set_boundary(circle, transform=ax.transAxes)


    cmap = settings.COLOURMAP_DICT["phenological_r"]
    norm = mpl.cm.colors.BoundaryNorm(BOUNDS, cmap.N)
    mesh = iris.plot.pcolormesh(plot_cube, cmap=cmap, norm=norm, axes=ax)

    # plot scatter
    COL = "yellow"
    ax.scatter(us_locations[1], us_locations[0], c=COL, s=100, edgecolor="k", transform=cartopy.crs.Geodetic(), zorder=10)
    ax.scatter(-2.9376, 54.3739, c=COL, s=100, edgecolor="k", transform=cartopy.crs.Geodetic(), zorder=10)
    
    # uk box
    region = [-10.0, 49.0, 3.0, 60.0]
    ax.plot([region[0], region[0]], [region[1], region[3]], c=COL, ls='-', lw=4, zorder=10, transform=cartopy.crs.PlateCarree())
    ax.plot([region[2], region[2]], [region[1], region[3]], c=COL, ls='-', lw=4, zorder=10, transform=cartopy.crs.PlateCarree())
    ax.plot([region[0], region[2]], [region[1], region[1]], c=COL, ls='-', lw=4, zorder=10, transform=cartopy.crs.Geodetic())
    ax.plot([region[0], region[2]], [region[3], region[3]], c=COL, ls='-', lw=4, zorder=10, transform=cartopy.crs.Geodetic())
    
    COL = "k"
    ax.plot([region[0], region[0]], [region[1], region[3]], c=COL, ls='-', lw=5, zorder=9, transform=cartopy.crs.PlateCarree())
    ax.plot([region[2], region[2]], [region[1], region[3]], c=COL, ls='-', lw=5, zorder=9, transform=cartopy.crs.PlateCarree())
    ax.plot([region[0], region[2]], [region[1], region[1]], c=COL, ls='-', lw=5, zorder=9, transform=cartopy.crs.Geodetic())
    ax.plot([region[0], region[2]], [region[3], region[3]], c=COL, ls='-', lw=5, zorder=9, transform=cartopy.crs.Geodetic())

    # label axes
    ax.text(-0.1, 1.0, LABELS, fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)

    cb = plt.colorbar(mesh, orientation='horizontal', ticks=BOUNDS[1:-1], label="Anomaly (days)", drawedges=True, fraction=0.1, pad=0.05, aspect=15, shrink=0.8)
    # prettify
    cb.set_ticklabels(["{:g}".format(b) for b in BOUNDS[1:-1]])
    cb.outline.set_linewidth(2)
    cb.dividers.set_color('k')
    cb.dividers.set_linewidth(2)

    ax.set_extent([-180, 180, 30, 90], cartopy.crs.PlateCarree())

    for lat in range(30, 100, 10):
        ax.text(180, lat, '{}$^\circ$N'.format(lat), transform=cartopy.crs.Geodetic())


    fig.subplots_adjust(bottom=0.05, top=0.95, left=0.04, right=0.95, wspace=0.02)

    del sos_cube
    del cubelist

    #***********************
    # MODIS timeseries - 2018
    
    sos_na, sos_ea, sprt_na_orig, sprt_ea_orig = read_modis_ts(os.path.join(data_loc, "MODIS.CMG.{}.SOS.EOS.SPRT.FALT.TS.csv".format(settings.YEAR)))

    dummy, sos_na = utils.calculate_climatology_and_anomalies_1d(sos_na, 2000, 2010)
    dummy, sos_ea = utils.calculate_climatology_and_anomalies_1d(sos_ea, 2000, 2010)

    dummy, sprt_na = utils.calculate_climatology_and_anomalies_1d(sprt_na_orig, 2000, 2010)
    dummy, sprt_ea = utils.calculate_climatology_and_anomalies_1d(sprt_ea_orig, 2000, 2010)

    
    ax = plt.axes([0.1, 0.7, 0.8, 0.15])
    label = "(b) North America"
    anomalies = ["2018 SOS Anomaly = 1.86 days", "2018 Spring T anomaly = -0.75 "+r'$^{\circ}$'+"C"]

    plot_modis_ts(ax, sos_na, sprt_na, sprt_na_orig, label, anomalies, LEGEND_LOC)

    na_trend = -0.64/10
    ax.plot([sos_na.times[0], sos_na.times[-1]], utils.trendline(na_trend, sos_na.times), ls="--", zorder = 10, c="g", lw=2)
    print(utils.trendline(na_trend, sos_na.times))

    ax1 = plt.axes([0.1, 0.85, 0.8, 0.15], sharex=ax)
    label = "(a) Eurasia"
    anomalies = ["2018 SOS Anomaly = 2.01 days", "2018 Spring T anomaly = 0.13 "+r'$^{\circ}$'+"C"]

    plot_modis_ts(ax1, sos_ea, sprt_ea, sprt_ea_orig, label, anomalies, "")

    ea_trend = -1.59/10
    ax1.plot([sos_ea.times[0], sos_ea.times[-1]], utils.trendline(ea_trend, sos_ea.times), ls="--", zorder = 10, c="g", lw=2)
    print(utils.trendline(ea_trend, sos_ea.times))
    
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax.get_xticklabels(), visible=True)

    fig.text(0.05, 0.92, "SOS Anomaly (days)", rotation="vertical")
    fig.text(0.95, 0.92, "Temperature Anomaly ("+r'$^{\circ}$'+"C)", rotation="vertical")

    plt.savefig(image_loc + "PHEN_modis_{}{}".format(settings.YEAR, settings.OUTFMT))

    #***********************
    # US timeseries - 2018

    fig = plt.figure(figsize=(8, 10))
    plt.clf()

    barrow, ozarks, turkey = read_us_phenocam_csv(os.path.join(data_loc, "Richardson PhenoCam plots for Sidebar.csv"))

    ax = plt.axes([0.1, 0.03, 0.35, 0.17])
    plot_us_phenocam(ax, barrow, [160, 229], "Barrow")
    ax.text(0.05, 0.85, "(e) Barrow", transform=ax.transAxes)
    ax = plt.axes([0.1, 0.23, 0.35, 0.17])
    plot_us_phenocam(ax, ozarks, [90, 159], "Ozarks", ylabel="Vegetation Greenness Index", legend = True)
    ax.text(0.05, 0.85, "(d) Ozarks", transform=ax.transAxes)
    ax = plt.axes([0.1, 0.43, 0.35, 0.17])
    plot_us_phenocam(ax, turkey, [60, 129], "Turkey Point")
    ax.text(0.05, 0.85, "(c) Turkey Point", transform=ax.transAxes)

    ax = plt.axes([0.48, 0.02, 0.25, 0.2])
    plot_images(ax, "barrow_2017_06_27_100002.jpg")
    ax = plt.axes([0.48, 0.22, 0.25, 0.2])
    plot_images(ax, "missouriozarks_2017_04_30_163113.jpg")
    ax = plt.axes([0.48, 0.42, 0.25, 0.2])
    plot_images(ax, "turkeypointenf39_2017_04_18_113106.jpg")
    
    ax = plt.axes([0.73, 0.02, 0.25, 0.2])
    plot_images(ax, "barrow_2018_06_27_103003.jpg")
    ax = plt.axes([0.73, 0.22, 0.25, 0.2])
    plot_images(ax, "missouriozarks_2018_04_30_163113.jpg")
    ax = plt.axes([0.73, 0.42, 0.25, 0.2])
    plot_images(ax, "turkeypointenf39_2018_04_18_113302.jpg")
    
    plt.figtext(0.6, 0.02, "2017")
    plt.figtext(0.85, 0.02, "2018")
    

    #***********************
    # UK timeseries - 2018

    ax = plt.axes([0.1, 0.65, 0.85, 0.15])

    oak = read_uk_oak_csv(os.path.join(data_loc, "UK Oak budburst 2000-2018.csv"))
    utils.plot_ts_panel(ax, [oak], "-", "phenological", loc="lower left")   
    ax.set_ylim([94, 124])
    ax.text(0.03, 0.8, "(b) Oak Budburst", transform=ax.transAxes)
    # fix the legend label to be italic
    handles, labels = ax.get_legend_handles_labels()
    labels[0] = "$\mathregular{\mathit{Quercus}}$"+" "+"$\mathregular{\mathit{robur}}$"
    ax.legend(handles, labels, frameon=False, loc="lower left")

    # need legend once have other panel's data
    north, south = read_windermere_csv(os.path.join(data_loc, "Windermere_2000-2018.csv"))
    ax2 = plt.axes([0.1, 0.8, 0.85, 0.15], sharex=ax)
    utils.plot_ts_panel(ax2, [north, south], "-", "phenological", loc="lower right")   
    ax2.set_ylim([79, 159])
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.text(0.03, 0.8, "(a) Windermere", transform=ax2.transAxes)

    ax2.set_xlim([1998, 2020])
    ax.set_xlim([1998, 2020])

    fig.text(0.03, 0.835, "Day of year", rotation = "vertical")

    plt.savefig(image_loc + "PHEN_timeseries_{}{}".format(settings.YEAR, settings.OUTFMT))
    plt.close()


    return # run_all_plots


#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
