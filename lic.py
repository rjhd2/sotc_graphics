#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for Lake Ice Cover (LIC) section.
#       For BAMS SotC 2019
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

DATALOC = "{}/{}/data/LIC/".format(settings.ROOTLOC, settings.YEAR)

LEGEND_LOC = 'lower left'

#************************************************************************
def read_column_csv(filename, era=True):

    raw_data = np.genfromtxt(filename, dtype=(float), skip_header=1, delimiter=",", encoding="latin-1")

    if era:
        name = "ERA5"
    else:
        name = "In Situ"

    times = raw_data[:, 0]
    on = utils.Timeseries(name, times, raw_data[:, 1])
    off = utils.Timeseries(name, times, raw_data[:, 2])
    duration = utils.Timeseries(name, times, raw_data[:, 3])

    return on, off, duration # read_column_csv

#************************************************************************
def read_separate_csv(filename, era=True):

    raw_data = np.genfromtxt(filename, dtype=(float), skip_header=1, delimiter=",", encoding="latin-1")

    times = raw_data[:, 0]
    insitu = utils.Timeseries("In Situ", times, raw_data[:, 1])
    era = utils.Timeseries("ERA5", times, raw_data[:, 2])

    return insitu, era # read_separate_csv

#************************************************************************
def read_continuous_csv(filename):

    all_data = []

    with open(filename, "r") as infile:

        time = []
        data = []
        name = ""
        for line in infile:
            sline = line.split(",")
            if sline[0].lower() == "year":
                continue
            else:
                if name != "" and name != sline[1]:
                    all_data += [utils.Timeseries("Lake", time, data)]
                    time = []
                    data = []
                # else just read in
                name = sline[1]
                time += [int(sline[0])]
                data += [float(sline[2])]
        # and last one
        all_data += [utils.Timeseries("All", time, data)]        
        
    return all_data # read_continuous_csv

#************************************************************************
def read_greatlakes_csv(filename, era=True):

    raw_data = np.genfromtxt(filename, dtype=(str), delimiter=",", encoding="latin-1")

    times = raw_data[1:, 0].astype(int)
    timeseries = []
    for col in range(raw_data.shape[1]):
        if col == 0:
            continue
        # if wanting to autogenerate lake names 
        ts = utils.Timeseries(raw_data[0, col], times, raw_data[1:, col].astype(float))
        if raw_data[0, col] == "Average":
            ts.lw=3
        else:
            ts.lw=2
        timeseries += [ts]
        # if raw_data[0, col] != "Average":
        #     timeseries += [utils.Timeseries("Lake", times, raw_data[1:, col].astype(float))]
        # else:
        #     timeseries += [utils.Timeseries("All", times, raw_data[1:, col].astype(float))]
            
    return timeseries # read_greatlakes_csv

#************************************************************************
def run_all_plots():

    cubelist = iris.load(os.path.join(DATALOC, "Fig1a-c.nc"))
    names = np.array([c.name() for c in cubelist])

    LABELS = {"Ice Start" : "(a) Ice On", "Ice End" : "(b) Ice Off", "Ice Duration" : "(c) Ice Duration"}
    COLORS = {"Ice Start" : plt.cm.RdBu_r, "Ice End" : plt.cm.RdBu, "Ice Duration" : plt.cm.RdBu}

    BOUNDS =  [-100, -20, -10, -5, -2, 0, 2, 5, 10, 20, 100]

    for name in names:
        if name == "Ice Depth":
            continue

        c, = np.where(names == name)[0]

        cube = cubelist[c]
       
        fig = plt.figure(figsize=(8, 9.5))
        plt.clf()

         # boundary circle
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)

        # axes for polar plot
        ax = plt.axes([0.01, 0.02, 0.98, 0.98], projection=cartopy.crs.NorthPolarStereo(central_longitude=300.0))

        plot_cube = cube

        # regrid depending on output format
        if settings.OUTFMT in [".eps", ".pdf"]:
            if plot_cube.coord("latitude").points.shape[0] > 90 or plot_cube.coord("longitude").points.shape[0] > 360:
                regrid_size = 1.0
                print("Regridding cube for {} output to {} degree resolution".format(settings.OUTFMT, regrid_size))
                print("Old Shape {}".format(plot_cube.data.shape))
                plot_cube = utils.regrid_cube(plot_cube, regrid_size, regrid_size)
                print("New Shape {}".format(plot_cube.data.shape))

        # prettify
        ax.gridlines() #draw_labels=True)
        ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
        ax.coastlines()
        ax.set_boundary(circle, transform=ax.transAxes)

        cmap = COLORS[name]
        norm = mpl.cm.colors.BoundaryNorm(BOUNDS, cmap.N)
        mesh = iris.plot.pcolormesh(plot_cube, cmap=cmap, norm=norm, axes=ax)

        # label axes
        ax.text(0.01, 1.0, "{}".format(LABELS[name]), fontsize=settings.FONTSIZE, transform=ax.transAxes)

        cb = plt.colorbar(mesh, orientation='horizontal', ticks=BOUNDS[1:-1], drawedges=True, fraction=0.1, pad=0.01, aspect=20, shrink=0.8)
        # prettify
        cb.set_label(label="Anomaly (days)", fontsize=settings.FONTSIZE)
        cb.ax.tick_params(axis='x', labelsize=settings.FONTSIZE, direction='in', size=0)
        cb.set_ticklabels(["{:g}".format(b) for b in BOUNDS[1:-1]])
        cb.outline.set_linewidth(2)
        cb.dividers.set_color('k')
        cb.dividers.set_linewidth(2)

        ax.set_extent([-180, 180, 30, 90], cartopy.crs.PlateCarree())

        for lat in range(30, 100, 10):
            ax.text(180, lat, '{}$^\circ$N'.format(lat), transform=cartopy.crs.Geodetic())

        plt.savefig(settings.IMAGELOC + "LIC_map_{}_{}{}".format(settings.YEAR, "".join(name.split()), settings.OUTFMT))


    #************************************************************************
    # and temperatures
    BOUNDS =  [-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100]
    cubelist = iris.load(os.path.join(DATALOC, "Fig1d.nc"))

    cube = cubelist[0]

    fig = plt.figure(figsize=(8, 9.5))
    plt.clf()

     # boundary circle
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    # axes for polar plot
    ax = plt.axes([0.01, 0.02, 0.98, 0.98], projection=cartopy.crs.NorthPolarStereo(central_longitude=300.0))

    plot_cube = cube

    # regrid depending on output format
    if settings.OUTFMT in [".eps", ".pdf"]:
        if plot_cube.coord("latitude").points.shape[0] > 90 or plot_cube.coord("longitude").points.shape[0] > 360:
            regrid_size = 1.0
            print("Regridding cube for {} output to {} degree resolution".format(settings.OUTFMT, regrid_size))
            print("Old Shape {}".format(plot_cube.data.shape))
            plot_cube = utils.regrid_cube(plot_cube, regrid_size, regrid_size)
            print("New Shape {}".format(plot_cube.data.shape))

    # prettify
    ax.gridlines() #draw_labels=True)
    ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
    ax.coastlines()
    ax.set_boundary(circle, transform=ax.transAxes)

    cmap = plt.cm.RdBu_r
    norm = mpl.cm.colors.BoundaryNorm(BOUNDS, cmap.N)
    mesh = iris.plot.pcolormesh(plot_cube, cmap=cmap, norm=norm, axes=ax)

    # label axes
    ax.text(0.01, 1.0, "(d) Nov-Apr Air Temperature", fontsize=settings.FONTSIZE, transform=ax.transAxes)

    cb = plt.colorbar(mesh, orientation='horizontal', ticks=BOUNDS[1:-1], drawedges=True, fraction=0.1, pad=0.01, aspect=20, shrink=0.8)
    # prettify
    cb.set_label(label="Anomaly ($^\circ$C)", fontsize=settings.FONTSIZE)
    cb.ax.tick_params(axis='x', labelsize=settings.FONTSIZE, direction='in', size=0)
    cb.set_ticklabels(["{:g}".format(b) for b in BOUNDS[1:-1]])
    cb.outline.set_linewidth(2)
    cb.dividers.set_color('k')
    cb.dividers.set_linewidth(2)

    ax.set_extent([-180, 180, 30, 90], cartopy.crs.PlateCarree())

    for lat in range(30, 100, 10):
        ax.text(180, lat, '{}$^\circ$N'.format(lat), transform=cartopy.crs.Geodetic())

    plt.savefig(settings.IMAGELOC + "LIC_map_{}_{}{}".format(settings.YEAR, "AirT", settings.OUTFMT))

    #************************************************************************
    # Timeseries
#    era_start, era_end, era_length = read_column_csv(os.path.join(DATALOC, "era5_lic_anom_1979_{}_NH.csv".format(settings.YEAR)))
#    situ_start, situ_end, situ_length = read_column_csv(os.path.join(DATALOC, "situ_lic_anom_1979_{}_NH.csv".format(settings.YEAR)), era=False)
    
    situ_start, era_start = read_separate_csv(os.path.join(DATALOC, "Fig2a.csv"))
    situ_end, era_end = read_separate_csv(os.path.join(DATALOC, "Fig2b.csv"))
    situ_length, era_length = read_separate_csv(os.path.join(DATALOC, "Fig2c.csv"))
    

    plt.clf()
    fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 10), sharex=True)
    
    # start
    utils.plot_ts_panel(ax1, [era_start, situ_start], "-", "cryosphere", loc="")

    # end
    utils.plot_ts_panel(ax2, [era_end, situ_end], "-", "cryosphere", loc="")
    ax2.set_ylabel("Anomaly (day)", fontsize=settings.FONTSIZE)

    # duration
    utils.plot_ts_panel(ax3, [era_length, situ_length], "-", "cryosphere", loc="lower left")

    # sort formatting
    for tick in ax3.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)

    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)
    for tick in ax2.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)
    for tick in ax3.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)

    # sort labelling
    ax1.text(0.02, 0.9, "(a) Start of Ice Cover", transform=ax1.transAxes, fontsize=settings.LABEL_FONTSIZE)
    ax2.text(0.02, 0.9, "(b) End of Ice Cover", transform=ax2.transAxes, fontsize=settings.LABEL_FONTSIZE)
    ax3.text(0.02, 0.9, "(c) Duration of Ice Cover", transform=ax3.transAxes, fontsize=settings.LABEL_FONTSIZE)

    fig.subplots_adjust(bottom=0.05, right=0.99, top=0.99, hspace=0.001)
  

    plt.savefig(settings.IMAGELOC + "LIC_ts_{}{}".format(settings.YEAR, settings.OUTFMT))

    #************************************************************************
    # Timeseries
#    indata = read_continuous_csv(os.path.join(DATALOC, "Fig3.csv"))
    indata = read_greatlakes_csv(os.path.join(DATALOC, "Fig3.csv"))

    fig = plt.figure(figsize=(8, 6))
    plt.clf()
    ax = plt.axes([0.1, 0.05, 0.88, 0.9])

    utils.plot_ts_panel(ax, indata, "-", "cryosphere", loc="lower left")

    fig.text(0.02, 0.5, "Maximum ice cover (anomaly, %)", va='center', rotation='vertical', ha="center", fontsize=settings.FONTSIZE)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)

    plt.savefig(settings.IMAGELOC + "LIC_GL_ts_{}{}".format(settings.YEAR, settings.OUTFMT))

    return # run_all_plots

 
#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
