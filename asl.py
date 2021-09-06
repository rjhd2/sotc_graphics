#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for Aerosols (ASL) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 31                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2021-09-06 09:52:46 +0100 (Mon, 06 Sep #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl

import copy
import numpy as np

import iris
import cartopy

import utils # RJHD utilities
import settings

DATALOC = "{}/{}/data/ASL/".format(settings.ROOTLOC, settings.YEAR)

LW = 3
LEGEND_LOC = "upper left"

#************************************************************************
def read_data(filename, name):
    """
    Read data file and returns Timeseries
    
    :param str filename: file to process
    :param str name: name to give Timeseries

    :returns: Timeseries object
    """

    indata = np.genfromtxt(filename, dtype=(str))

    date = indata[:, 0]
    data = indata[:, 1].astype(float)

    month = [d[5:7] for d in date]
    year = [d[:4] for d in date]

    times = np.array(year).astype(int) + (np.array(month).astype(int) - 1)/12.

    return utils.Timeseries(name, times, data) # read_data


#************************************************************************
#************************************************************************
def run_all_plots():

    #************************************************************************
    # Timeseries plot
    if True:

        monthly = read_data(DATALOC + "monthlist.txt", "AOD monthly")
        annual = read_data(DATALOC + "yearlist.txt", "AOD annual")

        minor_tick_interval = 1
        minorLocator = MultipleLocator(minor_tick_interval)
        major_tick_interval = 5
        majorLocator = MultipleLocator(major_tick_interval)
        COLOURS = settings.COLOURS["composition"]

        fig = plt.figure(figsize=(8, 5))
        ax = plt.axes([0.13, 0.07, 0.85, 0.86])

        plt.plot(monthly.times, monthly.data, COLOURS[monthly.name], ls='-', label=monthly.name, lw=LW)
        plt.plot(annual.times, annual.data, COLOURS[annual.name], ls='-', label=annual.name, lw=LW)

        ax.legend(loc=LEGEND_LOC, ncol=1, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, labelspacing=0.1, columnspacing=0.5)

        # prettify
        ax.xaxis.set_minor_locator(minorLocator)
        ax.xaxis.set_major_locator(majorLocator)
        utils.thicken_panel_border(ax)

        fig.text(0.03, 0.5, "AOD", va='center', rotation='vertical', fontsize=settings.FONTSIZE)

        plt.xlim([2002, int(settings.YEAR)+2])
        plt.ylim([0.09, 0.21])
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)

        plt.savefig(settings.IMAGELOC + "ASL_ts{}".format(settings.OUTFMT))
        plt.close()

    #************************************************************************
    # Global maps

    # actuals
    if True:
        total = iris.load(DATALOC + "AOD_{}.nc".format(settings.YEAR))
        bounds = np.arange(0, 1.1, 0.1)

        utils.plot_smooth_map_iris(settings.IMAGELOC + "ASL_total_actuals_{}".format(settings.YEAR), total[0][0], plt.cm.YlOrBr, bounds, "AOD")

    # total
    if True:
        total = iris.load(DATALOC + "diffAOD{}.nc".format(settings.YEAR[-2:]))
        bounds = np.array([-10, -0.15, -0.10, -0.05, -0.025, -0.01, 0.01, 0.025, 0.05, 0.10, 0.15, 10])

        utils.plot_smooth_map_iris(settings.IMAGELOC + "ASL_total_anomalies_{}".format(settings.YEAR), total[0][0], settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-{} (AOD)".format(settings.YEAR[2:]))
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_ASL_total_anomalies_{}".format(settings.YEAR), total[0][0], settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-20{} (AOD)".format(int(settings.YEAR[2:])-1), figtext="(x) Total Aerosol")

    # ratio AOD
    if True:
        extreme = iris.load(DATALOC + "ratioaod_{}.nc".format(settings.YEAR))
        bounds = np.arange(0, 2.2, 0.2)

        utils.plot_smooth_map_iris(settings.IMAGELOC + "ASL_ratio_{}".format(settings.YEAR), extreme[0][0], settings.COLOURMAP_DICT["composition"], bounds, "Ratio of AOD to 2003-19 average")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_ASL_ratio_{}".format(settings.YEAR), extreme[0][0], settings.COLOURMAP_DICT["composition"], bounds, "Ratio of AOD to 2003-19 average", figtext="(y) AOD Ratio")

    # extreme AOD
    if True:
        extreme = iris.load(DATALOC + "NB999.nc")#.format(settings.YEAR))
        extreme = extreme[0][0]
        bounds = np.array([0, 2.5, 5, 10, 20, 30, 40, 90])
        bounds = np.arange(0, 22, 2)
        extreme.data = np.ma.masked_where(extreme.data < 0, extreme.data)

        utils.plot_smooth_map_iris(settings.IMAGELOC + "ASL_extreme_days_{}".format(settings.YEAR), extreme, plt.cm.YlOrBr, bounds, "Number of days with AOD above the 99.9th percentile (days)")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_ASL_extreme_days_{}".format(settings.YEAR), extreme, plt.cm.YlOrBr, bounds, "Number of days with AOD above the 99.9th percentile (days)", figtext="(z) Extreme Aerosol Days")


    # Biomass Burning
    if False:
        BB = iris.load(DATALOC + "diffBB.nc")#.format(settings.YEAR))

        utils.plot_smooth_map_iris(settings.IMAGELOC + "ASL_BB_anomalies_{}".format(settings.YEAR), BB[0], settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-20{} (AOD)".format(int(settings.YEAR[2:])-1))
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_ASL_BB_anomalies_{}".format(settings.YEAR), BB[0], settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-20{} (AOD)".format(int(settings.YEAR[2:])-1), figtext="(y) Black Carbon & Organic Matter Aerosol")

    # Dust
    if False:
        dust = iris.load(DATALOC + "diffDU.nc")#.format(settings.YEAR))

        utils.plot_smooth_map_iris(settings.IMAGELOC + "ASL_dust_anomalies_{}".format(settings.YEAR), dust[0], settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-20{} (AOD)".format(int(settings.YEAR[2:])-1))
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_ASL_dust_anomalies_{}".format(settings.YEAR), dust[0], settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 2003-20{} (AOD)".format(int(settings.YEAR[2:])-1), figtext="(z) Dust Aerosol")

    #************************************************************************
    # Total, Trend and Event map
    if False:

        # as one colorbar per map, don't use the utils Iris routine

        fig = plt.figure(figsize=(8, 15))

        #total = iris.load(DATALOC + "Aerosol_Fig2a_Total_Averages.nc")
        #trend = iris.load(DATALOC + "Aerosol_Fig2b_Total_Trends.nc")
        #event = iris.load(DATALOC + "Aerosol_Fig2c_ExtremeDayCounts_{}.nc".format(settings.YEAR))

        total = iris.load(DATALOC + "meanAOD.nc")
        trend = iris.load(DATALOC + "trend.nc")
        number = iris.load(DATALOC + "numbermonth.nc")

        bounds_a = np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0])
        bounds_b = np.array([-10, -0.010, -0.005, -0.003, -0.002, -0.001, 0.001, 0.002, 0.003, 0.005, 0.010, 10])
        bounds_c = np.array([0, 1, 2, 3, 4, 5, 6, 7, 10])

        all_cubes = [total, trend, number]
        all_bounds = [bounds_a, bounds_b, bounds_c]

        # do colourmaps by hand
        cmap = [plt.cm.YlOrBr, settings.COLOURMAP_DICT["composition"], plt.cm.YlOrBr]
        PLOTLABELS = ["(a)", "(b)", "(c)"]
        cb_label = ["AOD", "AOD yr"+r'$^{-1}$', "Months"]

        # spin through axes
        for a in range(3):  

            norm = mpl.cm.colors.BoundaryNorm(all_bounds[a], cmap[a].N)

            ax = plt.subplot(3, 1, a+1, projection=cartopy.crs.Robinson())

            ax.gridlines() #draw_labels=True)
            ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
            ax.coastlines()

            ext = ax.get_extent() # save the original extent

            mesh = iris.plot.pcolormesh(all_cubes[a][0], cmap=cmap[a], norm=norm, axes=ax)

            ax.set_extent(ext, ax.projection) # fix the extent change from colormesh
            ax.text(0.0, 1.0, PLOTLABELS[a], fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)

            # sort the colourbar
            cb = fig.colorbar(mesh, ax=ax, orientation='horizontal', \
                                  ticks=all_bounds[a][1:-1], label=cb_label[a], drawedges=True, pad=0.05, fraction=0.07, aspect=30)
            if a == 1:
                cb.set_ticklabels(["{:6.3f}".format(b) for b in all_bounds[a][1:-1]])
            else:
                cb.set_ticklabels(["{:g}".format(b) for b in all_bounds[a][1:-1]])

            cb.outline.set_linewidth(2)
            cb.dividers.set_color('k')
            cb.dividers.set_linewidth(2)


        fig.subplots_adjust(bottom=0.05, top=0.95, left=0.04, right=0.95, wspace=0.02)

        plt.savefig(settings.IMAGELOC + "ASL_Trends{}".format(settings.OUTFMT))

        plt.close()

    #************************************************************************
    # Total, and two trends maps
    if True:

        # as one colorbar per map, don't use the utils Iris routine

        fig = plt.figure(figsize=(8, 15))

#        total = iris.load(DATALOC + "AOD2003-2009.nc")[0]
        total = iris.load(DATALOC + "AOD_03{}.nc".format(settings.YEAR[-2:]))[0]
#        trend1 = iris.load(DATALOC + "trend2003-2019.nc")[0]
        trend1 = iris.load(DATALOC + "medianaod18.nc")[0]
        trend1.data = np.ma.masked_where(trend1.data < -50, trend1.data)
        trend2 = iris.load(DATALOC + "medianaod12-20sig.nc")[0]
        trend2.data = np.ma.masked_where(trend2.data == -99, trend2.data)

#        bounds_a = np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0])
        bounds_a = np.arange(0, 1.1, 0.1)
        bounds_b = np.array([-10, -0.020, -0.010, -0.006, -0.004, -0.002, 0.002, 0.004, 0.006, 0.010, 0.020, 10])
        bounds_c = np.array([-10, -0.020, -0.010, -0.006, -0.004, -0.002, 0.002, 0.004, 0.006, 0.010, 0.020, 10])

        all_cubes = [total, trend1, trend2]
        all_bounds = [bounds_a, bounds_b, bounds_c]

        # do colourmaps by hand
        cmap = [plt.cm.YlOrBr, settings.COLOURMAP_DICT["composition"], settings.COLOURMAP_DICT["composition"]]
        PLOTLABELS = ["(a) Mean 2003-{}".format(settings.YEAR[-2:]), \
                      "(b) Trend 2003-{}".format(settings.YEAR[-2:]), \
                      "(c) Trend 2012-{}".format(settings.YEAR[-2:])]
        cb_label = ["AOD", "AOD yr"+r'$^{-1}$', "AOD yr"+r'$^{-1}$']

        # spin through axes
        for a in range(3):  

            norm = mpl.cm.colors.BoundaryNorm(all_bounds[a], cmap[a].N)

            ax = plt.subplot(3, 1, a+1, projection=cartopy.crs.Robinson())

            ax.gridlines() #draw_labels=True)
            ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
            ax.coastlines()

            ext = ax.get_extent() # save the original extent

            cube = all_cubes[a][0]
            if settings.OUTFMT in [".eps", ".pdf"]:
                if cube.coord("latitude").points.shape[0] > 180 or cube.coord("longitude").points.shape[0] > 360:
                    regrid_size = 1.0
                    print("Regridding cube for {} output to {} degree resolution".format(settings.OUTFMT, regrid_size))
                    print("Old Shape {}".format(cube.data.shape))
                    plot_cube = utils.regrid_cube(cube, regrid_size, regrid_size)
                    print("New Shape {}".format(plot_cube.data.shape))
                else:
                    plot_cube = copy.deepcopy(cube)
            else:
                plot_cube = copy.deepcopy(cube)


            mesh = iris.plot.pcolormesh(plot_cube, cmap=cmap[a], norm=norm, axes=ax)

            ax.set_extent(ext, ax.projection) # fix the extent change from colormesh
            ax.text(0.0, 1.0, PLOTLABELS[a], fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)

            # sort the colourbar
            cb = fig.colorbar(mesh, ax=ax, orientation='horizontal', \
                                  ticks=all_bounds[a][1:-1], label=cb_label[a], drawedges=True, pad=0.05, fraction=0.07, aspect=30)
            cb.set_ticklabels(["{:g}".format(b) for b in all_bounds[a][1:-1]])

            cb.outline.set_linewidth(2)
            cb.dividers.set_color('k')
            cb.dividers.set_linewidth(2)


        fig.subplots_adjust(bottom=0.05, top=0.95, left=0.04, right=0.95, wspace=0.02)

        plt.savefig(settings.IMAGELOC + "ASL_Trends{}".format(settings.OUTFMT))

        plt.close()

    #************************************************************************
    # Forcing map
    if True:

        bounds = np.array([-10, -1.4, -1.0, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 1.0, 1.4, 10])


        RFari = iris.load(DATALOC + "aerorf_camsra_rfari_2003-{}.nc".format(settings.YEAR))

        utils.plot_smooth_map_iris(settings.IMAGELOC + "ASL_RFari_mean_{}".format(settings.YEAR), RFari[0][0], settings.COLOURMAP_DICT["composition"], bounds, "2003-{} average of SW forcing from aerosol-cloud interactions".format(settings.YEAR), figtext="(c) CAMSRA: RFari (SW). Mean -0.45"+r'$\pm$'+"0.13 Wm"+r'$^{-2}$')

        RFaci = iris.load(DATALOC + "aerorf_camsra_rfaci_2003-{}.nc".format(settings.YEAR))

        utils.plot_smooth_map_iris(settings.IMAGELOC + "ASL_RFaci_mean_{}".format(settings.YEAR), RFaci[0][0], settings.COLOURMAP_DICT["composition"], bounds, "2003-{} average of SW forcing from aerosol-radiation interactions".format(settings.YEAR), figtext="(e) CAMSRA: RFaci (SW). Mean -0.63"+r'$\pm$'+"0.42 Wm"+r'$^{-2}$')

        bounds = np.arange(-0.05, 0.6, 0.05)
        AnthAOD = iris.load(DATALOC + "aerorf_camsra_od550anth_2003-{}.nc".format(settings.YEAR))

        utils.plot_smooth_map_iris(settings.IMAGELOC + "ASL_anthAOD_yearly_mean_{}".format(settings.YEAR), AnthAOD[0][0], plt.cm.YlOrBr, bounds, "Average of {} anthropogenic AOD".format(settings.YEAR), figtext="(a) CAMSRA: Anthropogenic AOD. Mean 0.065"+r'$\pm$'+"0.012")

        # and timeseries
        cubelist = iris.load(DATALOC + "aerorf_camsra_timeseries_2003-{}_v2.nc".format(settings.YEAR))

        LABELS = {"Anthropogenic AOD" : "(b)", "RFari" : "(d)", "RFaci" : "(f)"}
        ERRORS = {"Anthropogenic AOD" : 0.178, "RFari" : 0.286, "RFaci" : 0.667}

        for cube in cubelist:

            name = str(cube.var_name)
            if name[0] == "r":
                name = name[:2].upper()+name[2:]
                units = " (W m" + r'$^{-2}$'+")"
            else:
                name = "Anthropogenic AOD"
                units = ""

            fig = plt.figure(figsize=(8, 6))
            ax = plt.axes([0.13, 0.10, 0.8, 0.86])
            iris.plot.plot(cube, 'k', label=cube.var_name, lw=LW)

            timeUnits = cube.coord("time").units
            dt_time = timeUnits.num2date(cube.coord("time").points)
            upper = cube.data*(1 + ERRORS[name])
            lower = cube.data*(1 - ERRORS[name])
            ax.fill_between(dt_time, upper, lower, color="0.7", label=None)

            # prettify
            ax.set_ylabel("{}{}".format(name, units), fontsize=settings.FONTSIZE)

            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE) 
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE) 
            ax.yaxis.set_ticks_position('left')
            utils.thicken_panel_border(ax)

            minorLocator = MultipleLocator(365.242199) # in days since
            ax.xaxis.set_minor_locator(minorLocator)
            utils.thicken_panel_border(ax)

            # set labels
            if cube.var_name == "rfari":
                ax.set_ylim([-0.8, -0.3])
            elif cube.var_name == "rfaci":
                ax.set_ylim([-1.6, 0.])
            elif cube.var_name == "od550anth":
                ax.set_ylim([0.05, 0.085])

            ax.text(0.02, 0.9, LABELS[name], transform=ax.transAxes, fontsize=settings.FONTSIZE)

            plt.savefig(settings.IMAGELOC+"ASL_ts_{}{}".format(cube.var_name, settings.OUTFMT))
            plt.close()
 
    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
