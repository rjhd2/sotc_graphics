#!/usr/bin/env python
# python3
from __future__ import absolute_import
from __future__ import print_function
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


DATALOC = "{}/{}/data/PHEN/".format(settings.ROOTLOC, settings.YEAR)

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

    falt_nh = utils.Timeseries("Fall T", years, raw_data[:, 13])
    falt_na = utils.Timeseries("Fall T", years, raw_data[:, 14])
    falt_ea = utils.Timeseries("Fall T", years, raw_data[:, 15])

    return sos_nh, eos_nh, sprt_nh, falt_nh # read_modis_ts

#************************************************************************
def read_modis_uk_ts(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read
    :returns: Timeseries object s
    '''

    raw_data = np.genfromtxt(filename, dtype=(float), skip_header=1, delimiter=",")

    years = raw_data[:, 0]  

   
    # 2018 entries
    sos_uk = utils.Timeseries("MODIS", years, raw_data[:, 1])
    eos_uk = utils.Timeseries("MODIS", years, raw_data[:, 2])

    return sos_uk, eos_uk # read_modis_uk_ts

#************************************************************************
def read_us_phenocam(filename):
    
    raw_data = np.genfromtxt(filename, dtype=(str), skip_header=1)

    lat = raw_data[:, 1].astype(float)
    lon = raw_data[:, 0].astype(float)

    return lat, lon # read_us_phenocam


#************************************************************************
def plot_modis_ts(axl, sos, sprt, dummy, label, anomalies, legend_loc):

    utils.plot_ts_panel(axl, [sos, dummy], "-", "phenological", loc=legend_loc)   


    # make twin
    axr = axl.twinx()
    utils.plot_ts_panel(axr, [sprt], "-", "phenological", loc="")   

    # prettify
    axl.set_ylim([-10, 10])

    # labels
    axl.text(-0.17, 1.08, label, transform=axl.transAxes, fontsize=settings.FONTSIZE)
    axl.text(0.4, 0.88, anomalies[0], transform=axl.transAxes, fontsize=settings.FONTSIZE)
    axl.text(0.4, 0.78, anomalies[1], transform=axl.transAxes, fontsize=settings.FONTSIZE)

    # ticks etc
    minorLocator = MultipleLocator(1)
    majorLocator = MultipleLocator(5)
    for ax in [axl]:
        if "EOS" in anomalies[0]:
            axr.set_ylim([-3, 3])
            ax.set_ylabel("EOS Anomaly (days)", fontsize=settings.FONTSIZE, color='g')
        elif "SOS" in anomalies[0]:
            axr.set_ylim([3, -3])
            ax.set_ylabel("SOS Anomaly (days)", fontsize=settings.FONTSIZE, color='g')
        ax.tick_params(axis='y', color='g')

    for ax in [axr]:
        ax.yaxis.tick_right()
        ax.set_ylabel("Temperature Anomaly ("+r'$^{\circ}$'+"C)", fontsize=settings.FONTSIZE, color='m')
        ax.tick_params(axis='y', color='m')
        ax.yaxis.set_tick_params(right=True, which="both", width=2, direction="in")
        
    for ax in [axl, axr]:
        utils.thicken_panel_border(ax)
        ax.set_yticks(ax.get_yticks()[1:-1])
        ax.xaxis.set_minor_locator(minorLocator)
        ax.xaxis.set_major_locator(majorLocator)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
            tick.label2.set_fontsize(settings.FONTSIZE)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)

        ax.set_xlim([1999, int(settings.YEAR)+1])

    return # plot_modis_ts

#************************************************************************
def read_uk_oak_csv(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read
    :returns: Timeseries object s
    '''

    raw_data = np.genfromtxt(filename, dtype=(str), skip_header=1, skip_footer=2, delimiter=",")
    
    indata = raw_data[:, 1:].astype(float)
    indata = np.ma.masked_where(indata == -99, indata)
    
    oak_sos = utils.Timeseries("Quercus robur", raw_data[:, 0].astype(int), indata[:, 3])
    oak_eos = utils.Timeseries("Quercus robur", raw_data[:, 0].astype(int), indata[:, 0])
    
    return oak_sos, oak_eos  # read_uk_oak_csv

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

    raw_data = np.genfromtxt(filename, dtype=(str), delimiter=",", encoding="latin-1")

    times = raw_data[0, 2:].astype(int)
    modis_sos = utils.Timeseries("MODIS", times, raw_data[1, 2:].astype(int)) 
    modis_eos = utils.Timeseries("MODIS", times, raw_data[2, 2:].astype(int))

    # find phenocam offset
    start = np.where(raw_data[4, :] == "")[0][-1] + 1

    pheno_sos_10 = utils.Timeseries("PhenoCam 10%", raw_data[0, start:].astype(int), raw_data[4, start:].astype(int)) 
    pheno_sos_25 = utils.Timeseries("PhenoCam 25%", raw_data[0, start:].astype(int), raw_data[5, start:].astype(int)) 
    pheno_sos_50 = utils.Timeseries("PhenoCam 50%", raw_data[0, start:].astype(int), raw_data[6, start:].astype(int)) 
    pheno_sos_m = utils.Timeseries("", raw_data[0, start:].astype(int), raw_data[8, start:].astype(int)) 
    pheno_sos_p = utils.Timeseries("", raw_data[0, start:].astype(int), raw_data[9, start:].astype(int)) 

    pheno_sos = (pheno_sos_10, pheno_sos_25, pheno_sos_50, pheno_sos_m, pheno_sos_p)
    
    pheno_eos_10 = utils.Timeseries("PhenoCam 10%", raw_data[0, start:].astype(int), raw_data[11, start:].astype(int)) 
    pheno_eos_25 = utils.Timeseries("PhenoCam 25%", raw_data[0, start:].astype(int), raw_data[12, start:].astype(int)) 
    pheno_eos_50 = utils.Timeseries("PhenoCam 50%", raw_data[0, start:].astype(int), raw_data[13, start:].astype(int)) 
    pheno_eos_m = utils.Timeseries("", raw_data[0, start:].astype(int), raw_data[15, start:].astype(int)) 
    pheno_eos_p = utils.Timeseries("", raw_data[0, start:].astype(int), raw_data[16, start:].astype(int)) 

    pheno_eos = (pheno_eos_10, pheno_eos_25, pheno_eos_50, pheno_eos_m, pheno_eos_p)

    return modis_sos, modis_eos, pheno_sos, pheno_eos # read_us_phenocam

#************************************************************************
def plot_us_phenocam(ax, modis, pheno, sos=True):

    from matplotlib.ticker import MultipleLocator

    minor_tick_interval=1
    minorLocator = MultipleLocator(minor_tick_interval)
    major_tick_interval=5
    majorLocator = MultipleLocator(major_tick_interval)

    if sos:
        colors = ["#31a354", "#a1d99b", "#e5f5e0"]
    else:
        colors = ["#d95f0e", "#fec44f", "#fff7bc"]

    pheno_10, pheno_25, pheno_50, pheno_m, pheno_p = pheno

    ax.plot(modis.times, modis.data, c="k", ls="-", label=modis.name, lw=3)

    ax.plot(pheno_10.times, pheno_10.data, c=colors[0], ls="-", label=pheno_10.name, lw=3)
    ax.plot(pheno_25.times, pheno_25.data, c=colors[1], ls="-", label=pheno_25.name, lw=3)
    ax.errorbar(pheno_25.times, pheno_25.data, yerr=[pheno_m.data, pheno_p.data], c=colors[1], ls="-", label=None)
    ax.plot(pheno_50.times, pheno_50.data, c=colors[2], ls="-", label=pheno_50.name, lw=3)

    if sos:
        ax.legend(loc="upper left", frameon=False, ncol=1, fontsize=settings.FONTSIZE*0.8)
    else:
        ax.legend(loc="lower left", frameon=False, ncol=1, fontsize=settings.FONTSIZE*0.8)

    ax.xaxis.set_minor_locator(minorLocator)
    ax.xaxis.set_major_locator(majorLocator)
    # turn of RHS y ticks
    ax.yaxis.set_ticks_position('left')
    ax.set_xlim([2000, 2020])
#    ax.set_ylabel("Day of Year", fontsize=settings.FONTSIZE)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE)

    utils.thicken_panel_border(ax)

    return # plot_us_phenocam

#************************************************************************
def plot_images(ax, filename):

    img = mpimg.imread(os.path.join(DATALOC, filename))
    ax.imshow(img)
    ax.set_xticks([])
    ax.set_yticks([])

    return # plot_images

#************************************************************************
def run_all_plots():


    #***********************
    # MODIS - centre
    if True:
        cubelist = iris.load(os.path.join(DATALOC, "MODIS.CMG.{}.SOS.EOS.Anomaly.nc".format(settings.YEAR)))
        names = np.array([c.name() for c in cubelist])
        # set up plot settings
        BOUNDS = [-100, -20, -10, -5, -2, 0, 2, 5, 10, 20, 100]
        LABELS = {"SOS": "(c) Start of Season (SOS)", "EOS": "(d) End of Season (EOS)"}

        for season in ["SOS", "EOS"]:

            c, = np.where(names == season)[0]

            cube = cubelist[c]

            # deal with NANS
            cube.data = np.ma.masked_where(cube.data != cube.data, cube.data)


            fig = plt.figure(figsize=(8, 11))
            plt.clf()

             # boundary circle
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)

            # axes for polar plot
            ax = plt.axes([0.01, 0.02, 0.98, 0.65], projection=cartopy.crs.NorthPolarStereo(central_longitude=300.0))

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

            if season == "SOS":
                cmap = settings.COLOURMAP_DICT["phenological_r"]
            elif season == "EOS":
                cmap = settings.COLOURMAP_DICT["phenological"]

            norm = mpl.cm.colors.BoundaryNorm(BOUNDS, cmap.N)
            mesh = iris.plot.pcolormesh(plot_cube, cmap=cmap, norm=norm, axes=ax)

            # # read in sites
            if season == "EOS":
                pass
            elif season == "SOS":
                lake_locations = read_us_phenocam(os.path.join(DATALOC, "lake_coords.csv"))   
                # scatter
                COL = "chartreuse"
                ax.scatter(lake_locations[1], lake_locations[0], c=COL, s=100, edgecolor="k", transform=cartopy.crs.Geodetic(), zorder=10)
                COL = "m"
                # Harvard Forest - 2019
                ax.scatter(-72.17, 42.54, c=COL, s=100, edgecolor="k", transform=cartopy.crs.Geodetic(), zorder=10)

                # uk box
                COL = "y"
                region = [-10.0, 49.0, 3.0, 60.0]
                ax.plot([region[0], region[0]], [region[1], region[3]], c=COL, ls='-', lw=4, zorder=10, transform=cartopy.crs.PlateCarree())
                ax.plot([region[2], region[2]], [region[1], region[3]], c=COL, ls='-', lw=4, zorder=10, transform=cartopy.crs.PlateCarree())
                ax.plot([region[0], region[2]], [region[1], region[1]], c=COL, ls='-', lw=4, zorder=10, transform=cartopy.crs.Geodetic())
                ax.plot([region[0], region[2]], [region[3], region[3]], c=COL, ls='-', lw=4, zorder=10, transform=cartopy.crs.Geodetic())

            # COL = "k"
            # ax.plot([region[0], region[0]], [region[1], region[3]], c=COL, ls='-', lw=5, zorder=9, transform=cartopy.crs.PlateCarree())
            # ax.plot([region[2], region[2]], [region[1], region[3]], c=COL, ls='-', lw=5, zorder=9, transform=cartopy.crs.PlateCarree())
            # ax.plot([region[0], region[2]], [region[1], region[1]], c=COL, ls='-', lw=5, zorder=9, transform=cartopy.crs.Geodetic())
            # ax.plot([region[0], region[2]], [region[3], region[3]], c=COL, ls='-', lw=5, zorder=9, transform=cartopy.crs.Geodetic())

            # label axes
            ax.text(-0.1, 1.0, LABELS[season], fontsize=settings.FONTSIZE, transform=ax.transAxes)

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

            fig.subplots_adjust(bottom=0.05, top=0.95, left=0.04, right=0.95, wspace=0.02)

            del cube

            #***********************
            # MODIS timeserise
            sos_nh, eos_nh, sprt_nh_orig, falt_nh_orig = read_modis_ts(os.path.join(DATALOC, "MODIS.CMG.{}.SOS.EOS.SPRT.FALT.TS.csv".format(settings.YEAR)))

            dummy, sos_nh = utils.calculate_climatology_and_anomalies_1d(sos_nh, 2000, 2010)
            dummy, eos_nh = utils.calculate_climatology_and_anomalies_1d(eos_nh, 2000, 2010)

            dummy, sprt_nh = utils.calculate_climatology_and_anomalies_1d(sprt_nh_orig, 2000, 2010)
            dummy, falt_nh = utils.calculate_climatology_and_anomalies_1d(falt_nh_orig, 2000, 2010)

            ax = plt.axes([0.15, 0.73, 0.75, 0.23])
            if season == "SOS":
                label = "(a) Start of Season"
                anomalies = ["{} SOS Anomaly = -4.3 days".format(settings.YEAR), "{} Spr. T anomaly = 0.19 ".format(settings.YEAR)+r'$^{\circ}$'+"C"]

                plot_modis_ts(ax, sos_nh, sprt_nh, sprt_nh_orig, label, anomalies, LEGEND_LOC)

            elif season == "EOS":
                label = "(b) End of Season"
                anomalies = ["{} EOS Anomaly = 2.4 days".format(settings.YEAR), "{} Fall T anomaly = -0.53 ".format(settings.YEAR)+r'$^{\circ}$'+"C"]

                plot_modis_ts(ax, eos_nh, falt_nh, falt_nh_orig, label, anomalies, LEGEND_LOC)


            plt.savefig(settings.IMAGELOC + "PHEN_modis_{}_{}{}".format(settings.YEAR, season, settings.OUTFMT))

        del cubelist


   #***********************
    # US timeseries - 2018
    if True:
        fig = plt.figure(figsize=(8, 9.5))
        plt.clf()

        modis_sos, modis_eos, pheno_sos, pheno_eos = read_us_phenocam_csv(os.path.join(DATALOC, "Richardson Data for SOC 2019 Figures.csv"))

        # images
        ax = plt.axes([0.01, 0.66, 0.49, 0.3])
        plot_images(ax, "HarvardForest_20190511.jpg")
        ax = plt.axes([0.5, 0.66, 0.49, 0.3])
        plot_images(ax, "HarvardForest_20191024.jpg")

        # timeseries
        ax = plt.axes([0.11, 0.35, 0.84, 0.3])
        plot_us_phenocam(ax, modis_eos, pheno_eos, sos=False)
#        ax.text(0.05, 0.85, "(a)", transform=ax.transAxes, fontsize=settings.FONTSIZE)
        ax.set_ylim([275, 369])
        plt.setp(ax.get_xticklabels(), visible=False)

        ax = plt.axes([0.11, 0.05, 0.84, 0.3])
        plot_us_phenocam(ax, modis_sos, pheno_sos)
#        ax.text(0.05, 0.85, "(b)", transform=ax.transAxes, fontsize=settings.FONTSIZE)
        ax.set_ylim([101, 144])

        fig.text(0.02, 0.97, "(a)", transform=ax.transAxes, fontsize=settings.FONTSIZE)
        fig.text(0.02, 0.3, "Day of year", rotation = "vertical", fontsize=settings.FONTSIZE)
        plt.savefig(settings.IMAGELOC + "PHEN_UStimeseries_{}{}".format(settings.YEAR, settings.OUTFMT))
        plt.close()



    #***********************
    # US timeseries - 2018
    if False:
        fig = plt.figure(figsize=(8, 7))
        plt.clf()

        modis_sos, modis_eos, pheno_sos, pheno_eos = read_us_phenocam_csv(os.path.join(DATALOC, "Richardson Data for SOC 2019 Figures.csv"))

        # timeseries
        ax = plt.axes([0.11, 0.5, 0.64, 0.45])
        plot_us_phenocam(ax, modis_eos, pheno_eos, sos=False)
        ax.text(0.05, 0.85, "(c)", transform=ax.transAxes, fontsize=settings.FONTSIZE)
        ax.set_ylim([275, 369])
        plt.setp(ax.get_xticklabels(), visible=False)
        ax = plt.axes([0.75, 0.5, 0.25, 0.4])
        plot_images(ax, "HarvardForest_20191024.jpg")

        ax = plt.axes([0.11, 0.05, 0.64, 0.45])
        plot_us_phenocam(ax, modis_sos, pheno_sos)
        ax.text(0.05, 0.85, "(d)", transform=ax.transAxes, fontsize=settings.FONTSIZE)
        ax.set_ylim([101, 144])
        ax = plt.axes([0.75, 0.05, 0.25, 0.4])
        plot_images(ax, "HarvardForest_20190511.jpg")

        fig.text(0.02, 0.4, "Day of year", rotation = "vertical", fontsize=settings.FONTSIZE)
        plt.savefig(settings.IMAGELOC + "PHEN_UStimeseries_{}{}".format(settings.YEAR, settings.OUTFMT))
        plt.close()
   

    #***********************
    # UK timeseries - 2018
    if True:
        from matplotlib.ticker import MultipleLocator
        majorLocator = MultipleLocator(5)

        fig = plt.figure(figsize=(8, 9.5))
        plt.clf()

         # images
        ax = plt.axes([0.01, 0.66, 0.48, 0.3])
        plot_images(ax, "Sarah Burgess first leaf.jpg")
        ax = plt.axes([0.5, 0.66, 0.48, 0.3])
        plot_images(ax, "Judith Garforth oak bare tree 2019.jpg")

        sos_uk, eos_uk = read_modis_uk_ts(os.path.join(DATALOC, "MODIS.CMG.{}.SOS.EOS.SPRT.FALT.TS.UK_DH.csv".format(settings.YEAR)))
        oak_sos, oak_eos = read_uk_oak_csv(os.path.join(DATALOC, "UK_Oakleaf_data.csv"))

        # timeseries
        ax = plt.axes([0.11, 0.35, 0.84, 0.3])
        utils.plot_ts_panel(ax, [oak_eos, eos_uk], "-", "phenological", loc="center left")
#        ax.text(0.05, 0.85, "(c)", transform=ax.transAxes, fontsize=settings.FONTSIZE)
        ax.set_ylim([210, 354])
        plt.setp(ax.get_xticklabels(), visible=False)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
        ax.xaxis.set_major_locator(majorLocator)

        # naughtily, manually tweak the upper oak plot
        for line in ax.get_lines():
            if line.get_color() == "g":
                line.set_color("#d95f0e")
        leg = ax.get_legend()
        for line in leg.get_lines():
            if line.get_color() == "g":
                line.set_color("#d95f0e")
    
        ax = plt.axes([0.11, 0.05, 0.84, 0.3])
        utils.plot_ts_panel(ax, [oak_sos, sos_uk], "-", "phenological", loc="center left")
#        ax.text(0.05, 0.85, "(d)", transform=ax.transAxes, fontsize=settings.FONTSIZE)
        ax.set_ylim([66, 132])
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
        ax.xaxis.set_major_locator(majorLocator)

 
        fig.text(0.02, 0.97, "(b)", transform=ax.transAxes, fontsize=settings.FONTSIZE)
        fig.text(0.02, 0.3, "Day of year", rotation = "vertical", fontsize=settings.FONTSIZE)

        plt.savefig(settings.IMAGELOC + "PHEN_UKtimeseries_{}{}".format(settings.YEAR, settings.OUTFMT))
        plt.close()

    #***********************
    # Lake Boxplot
    if True:

        import pandas as pd

        df = pd.read_csv(DATALOC + "LakeData_forRobert.csv")

        # rename columns
        cols = []
        for col in df.columns:
            if len(col.split()) >= 2:
                df.rename(columns={col: col.split()[0]}, inplace=True)
                cols += [col.split()[0]]

        fig = plt.figure(figsize=(8, 7))
        plt.clf()
        ax = plt.axes([0.1, 0.25, 0.89, 0.74])
        df.boxplot(column=cols, ax=ax, grid=False, )

        # messily pull out 2019
        this_year = df.iloc[-1]
        this_year = this_year.to_frame()
        plt.plot(np.arange(11)+1, this_year[19][1:], "ro")

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
            tick.label.set_rotation("vertical")

        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 

        utils.thicken_panel_border(ax)
 
        plt.ylabel("Day of year", fontsize=settings.FONTSIZE)
        plt.savefig(settings.IMAGELOC + "PHEN_lakes_boxplot_{}{}".format(settings.YEAR, settings.OUTFMT))
        plt.close()

        
    return # run_all_plots


#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
