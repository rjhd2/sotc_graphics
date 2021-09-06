#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for lower troposphere temperature (LTT) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 30                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2021-06-15 10:41:02 +0100 (Tue, 15 Jun #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import os
import datetime as dt
import copy
import numpy as np
import matplotlib.pyplot as plt

import iris
import cf_units

import utils # RJHD utilities
import settings



DATALOC = "{}/{}/data/LTT/".format(settings.ROOTLOC, settings.YEAR)

CLIMSTART = 1981
CLIMEND = 2010

DECIMAL_MONTHS = np.arange(12)/12.

LEGEND_LOC = 'lower right'


#************************************************************************
def get_times(time_coord):

    # read in data and extract start time
    points = time_coord.points
    unit = str(time_coord.units)


    start = dt.datetime.strptime(unit.split(" ")[-1], "%Y-%M-%d")
    if "months" in unit:
        # calculate the years
        years = start.year + np.floor(time_coord.points/12)
        months = ((time_coord.points%12)+1)/12
        return years + months # get_times

    else:
        # this one is easier
        years = (start.year + time_coord.points) + 0.5 # plot in centre of year

        return years # get_times


#************************************************************************
def read_csv(filename):
    """
    Read user supplied CSV for LTT into Timeseries object
    """

    indata = np.genfromtxt(filename, delimiter=',', dtype=(float), skip_header=7)

    indata = np.ma.masked_where(indata == "", indata)

    raobcore = utils.Timeseries("RAOBCORE v1.7", indata[:, 0], indata[:, 1])
    rich = utils.Timeseries("RICH v1.7", indata[:, 0], indata[:, 2])
    ratpac = utils.Timeseries("RATPAC A2", indata[:, 0], indata[:, 3])
#    unsw = utils.Timeseries("UNSW v1.0", indata[:, 0], indata[:, 4])

    UAH = utils.Timeseries("UAH v6.0", indata[:, 0], indata[:, 4])
    rss = utils.Timeseries("RSS v4.0", indata[:, 0], indata[:, 5])

    era5 = utils.Timeseries("ERA5", indata[:, 0], indata[:, 6])
    jra = utils.Timeseries("JRA-55", indata[:, 0], indata[:, 7])
    merra = utils.Timeseries("MERRA-2", indata[:, 0], indata[:, 8])

    return raobcore, rich, ratpac, UAH, rss, era5, jra, merra # read_csv


#************************************************************************
def moving_average(data, N=12):
    cumsum = np.cumsum(np.insert(data, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N) # moving_average

#************************************************************************
def read_netcdf_ts(filename, smooth=True):

    cubelist = iris.load(os.path.join(DATALOC, filename))

    for cube in cubelist:

        if "ratpac" in cube.name():
            # fix cube for 2020 run
            cube.coords()[0].standard_name="time"

        times = get_times(cube.coord("time"))
        data = np.ma.array(cube.data)
        data.mask = np.zeros(data.shape[0])

        # deal with NaNs
        data = np.ma.masked_where(np.isnan(data), data)

        goods, = np.where(data.mask == False)
        times = times[goods]
        data = data[goods]

        if "ratpac" not in cube.name() and smooth:

            N = 12
            newdata = moving_average(data, N)
            newtimes = times[N-2:-1]
        else:
            newdata = data
            newtimes = times
            
        if cube.name() ==  "tlt_global_ratpacva2":
            ratpac = utils.Timeseries("RATPAC vA2", newtimes, newdata)
        elif cube.name() ==  "tlt_global_rssv4":
            rss = utils.Timeseries("RSS v4.0", newtimes, newdata)
        elif cube.name() ==  "tlt_global_richv17":
            rich = utils.Timeseries("RICH v1.7", newtimes, newdata)
        elif cube.name() ==  "tlt_global_raobcorev17":
            raobcore = utils.Timeseries("RAOBCORE v1.7", newtimes, newdata)
        elif cube.name() ==  "tlt_global_uahv6":
            uah = utils.Timeseries("UAH v6.0", newtimes, newdata)
        elif cube.name() ==  "tlt_global_era51":
            era = utils.Timeseries("ERA5", newtimes, newdata)
        elif cube.name() ==  "tlt_global_merra2":
            merra = utils.Timeseries("MERRA-2", newtimes, newdata)
        elif cube.name() ==  "tlt_global_jra55":
            jra = utils.Timeseries("JRA-55", newtimes, newdata)

    return raobcore, rich, ratpac, uah, rss, era, merra, jra # read_netcdf_ts

#************************************************************************
def read_mei(filename):
    """
    Read MEI file into Timeseries object

    http://www.esrl.noaa.gov/psd/enso/mei/table.html

    :param str filename: input file to read

    :returns: Timeseries

    """

    with open(filename, "r") as infile:

        all_data = []

        for lc, line in enumerate(infile):
            # issues with np.genfromtxt for html, so do long way around
            if lc <= 12:
                continue

            elif line.strip() == "":
                break

            else:
                all_data += [line.split()]

    all_data = np.array(all_data).astype(float)

    years = all_data[:, 0]

    # create the times
    times = []
    year = years[0]
    while year <= years[-1]:

        times += list(year + DECIMAL_MONTHS)
        year += 1

    # unravel the data
    data = all_data[:, 1:]
    data = data.reshape(-1)

    assert len(times) == len(data)

    mei = utils.Timeseries("MEI", np.array(times), data)

    return mei # read_mei

#************************************************************************
def read_hilo(filename):

    indata = np.genfromtxt(filename, dtype=(float))

    times = indata[:, 0]

    high = utils.Timeseries("Highest", times, indata[:, 1])
    low = utils.Timeseries("Lowest", times, indata[:, 2])

    return high, low # read_hilo

#************************************************************************
def run_all_plots():

    #************************************************************************
    # Timeseries figure (3 panels)
    if True:
        for region in ["global"]:

            if region == "global":
#                raobcore, rich, ratpac, UAH, rss, era5, merra, jra = \
#                    read_csv(DATALOC + "SotC_AnnTemps_2020_0520_LTTGL.csv")

                raobcore, rich, ratpac, UAH, rss, era5, merra, jra = \
                    read_netcdf_ts(DATALOC + "dunn_timeseries.nc", smooth=True)
                

            plt.clf()
            fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 10), sharex=True)

            # sondes
            utils.plot_ts_panel(ax1, [raobcore, ratpac, rich], "-", "temperature", loc=LEGEND_LOC)

            # satellites
            utils.plot_ts_panel(ax2, [UAH, rss], "-", "temperature", loc=LEGEND_LOC)

            ax2.set_ylabel("Anomaly ("+r'$^\circ$'+"C)", fontsize=settings.FONTSIZE)


            # reanalyses
            if region == "global":
    #            jra_actuals, jra_anoms = utils.read_jra55(settings.REANALYSISLOC + "JRA-55_MSUch2LT_global_ts.txt", "temperature")
    #            merra_actuals, merra_anoms = utils.read_merra_LT_LS(settings.REANALYSISLOC + "MERRA2_MSU_Tanom_ann_{}.dat".format(settings.YEAR), LT=True)
    #            utils.plot_ts_panel(ax3, [erai, era5, jra_anoms, merra_anoms], "-", "temperature", loc=LEGEND_LOC)
    #            twenty_cr_actuals = utils.read_20cr(settings.REANALYSISLOC + "tlt.global.txt", "temperature")
    #            dummy, twenty_cr_anoms = utils.calculate_climatology_and_anomalies_1d(twenty_cr_actuals, 1981, 2010)


                utils.plot_ts_panel(ax3, [era5, jra, merra], "-", "temperature", loc=LEGEND_LOC)
            else:
                utils.plot_ts_panel(ax3, [era5], "-", "temperature", loc=LEGEND_LOC)


            # sort formatting
            plt.xlim([jra.times[0]-1, jra.times[-1]+1])
            ax1.set_ylim([-0.89, 0.99])
            ax2.set_ylim([-0.89, 0.99])
            ax3.set_ylim([-0.89, 0.99])

            for tick in ax3.xaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)

            for tick in ax1.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)
            for tick in ax2.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)
            for tick in ax3.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)

            # sort labelling
            ax1.text(0.02, 0.9, "(a) Radiosondes", transform=ax1.transAxes, fontsize=settings.LABEL_FONTSIZE)
            ax2.text(0.02, 0.9, "(b) Satellites", transform=ax2.transAxes, fontsize=settings.LABEL_FONTSIZE)
            ax3.text(0.02, 0.9, "(c) Reanalyses", transform=ax3.transAxes, fontsize=settings.LABEL_FONTSIZE)

            fig.subplots_adjust(right=0.99, top=0.98, bottom=0.05, hspace=0.001, left=0.15)

            plt.savefig(settings.IMAGELOC+"LTT_ts_{}{}".format(region, settings.OUTFMT))

            plt.close()

    # #************************************************************************
    # # Land Fraction - 2018

    # high, low = read_hilo(DATALOC + "high_low.dat")
    # fig, ax1 = plt.subplots(1, figsize=(8, 5))

    # ax1.plot(high.times, high.data, c="r", ls="-", lw=2, label=high.name)
    # ax1.plot(low.times, low.data, c="b", ls="-", lw=2, label=low.name)
    
    # ax1.legend(loc="upper right", ncol=2, frameon=False, prop={'size':settings.FONTSIZE})
    # ax1.set_xlim([1978, int(settings.YEAR)+2])
    # ax1.set_ylabel("Percentage of Global Area (%)", fontsize=settings.FONTSIZE)

    # for tick in ax1.yaxis.get_major_ticks():
    #     tick.label.set_fontsize(settings.FONTSIZE)
    # for tick in ax1.xaxis.get_major_ticks():
    #     tick.label.set_fontsize(settings.FONTSIZE)
    # utils.thicken_panel_border(ax1)

    # plt.savefig(settings.IMAGELOC+"LTT_land_area_ts{}".format(settings.OUTFMT))

    # plt.close()

    #************************************************************************
    # Read in ERA5 anomalies
    if False:
        cube_list = iris.load(DATALOC + "2019TLTAnom.nc")
        names = np.array([c.var_name for c in cube_list])

        loc, = np.where(names == "tltAnnualAnom")[0]

        cube = cube_list[loc]

        bounds = np.array([-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100])
        bounds = np.array([-100, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 100])
        bounds = np.array([-100, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 100])

        cmap = settings.COLOURMAP_DICT["temperature"]
        utils.plot_smooth_map_iris(settings.IMAGELOC + "LTT_{}_anoms_era5".format(settings.YEAR), cube, cmap, \
                                       bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", title="ERA5")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_LTT_{}_anoms_era5".format(settings.YEAR), cube, \
                                       cmap, bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", \
                                       figtext="(e) Lower Tropospheric Temperature")


    #************************************************************************
    # Read in Annual anomalies with record pixels
    if True:
        cube_list = iris.load(DATALOC + "dunn_records.nc")
        names = np.array([c.var_name for c in cube_list])

        loc, = np.where(names == "annual_anomaly_map")[0]

        cube = cube_list[loc]

        # get the RSS data to scatter
        loc, = np.where(names == "rss_record_mask")[0]
        record = cube_list[loc]

        lats, lons, data = [], [], []
        for t, lat in enumerate(record.coord("latitude").points):
            for n, lon in enumerate(record.coord("longitude").points):
                if record.data[t, n] == 1:
                    data += [cube.data[t, n]]
                    lats += [np.mean(record.coord("latitude").bounds[t])]
                    lons += [np.mean(record.coord("longitude").bounds[n])]

        # get the UAH data to hatch
        hatching = copy.deepcopy(cube) # store data from the map
        loc, = np.where(names == "uah_record_mask")[0]
        record = cube_list[loc]
        # apply the mask from UAH
        hatching.data = np.ma.masked_where(record.data == 0, hatching.data)
        
        bounds = np.array([-100, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 100])

        cmap = settings.COLOURMAP_DICT["temperature"]
        utils.plot_smooth_map_iris(settings.IMAGELOC + "LTT_{}_anoms_satellite".format(settings.YEAR), cube, cmap, \
                                       bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", title="Satellite", \
                                   scatter=(lons, lats, data), smarker="dots", hatch=hatching)
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_LTT_{}_anoms_satellite".format(settings.YEAR), cube, cmap, \
                                       bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", \
                                   figtext="(e) Lower Tropospheric Temperature", \
                                   scatter=(lons, lats, data), smarker="dots", hatch=hatching)


    #************************************************************************
    # Read in ERA anomalies with record pixels
    if False:
        cube_list = iris.load(DATALOC + "2019TLTAnom_warmest.nc")
        names = np.array([c.var_name for c in cube_list])

        loc, = np.where(names == "tltAnnualAnom")[0]

        cube = cube_list[loc]

        loc, = np.where(names == "recordMask")[0]
        record = cube_list[loc]
        record.coord("latitude").guess_bounds()
        record.coord("longitude").guess_bounds()
 
        lats, lons, data = [], [], []
        for t, lat in enumerate(record.coord("latitude").points):
            for n, lon in enumerate(record.coord("longitude").points):
                if record.data[t, n] == 1:
                    data += [cube.data[t, n]]
                    lats += [np.mean(record.coord("latitude").bounds[t])]
                    lons += [np.mean(record.coord("longitude").bounds[n])]
        
        bounds = np.array([-100, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 100])

        cmap = settings.COLOURMAP_DICT["temperature"]
        utils.plot_smooth_map_iris(settings.IMAGELOC + "LTT_{}_anoms_era5".format(settings.YEAR), cube, cmap, \
                                       bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", title="ERA5", \
                                   scatter=(lons, lats, data), smarker="dots")


        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_LTT_{}_anoms_era5".format(settings.YEAR), cube, \
                                       cmap, bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", \
                                       figtext="(e) Lower Tropospheric Temperature", \
                                   scatter=(lons, lats, data), smarker="dots")

    #************************************************************************
    # ERA-I Hovmuller

    # times, latitudes, data = utils.erai_2dts_read(settings.REANALYSISLOC, "ltt")

    # bounds = np.array([-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100])
    # bounds = np.array([-100, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 100])

    # # sort the time axis
    # start = (times[0] - 101)/10000.
    # end = (times[-1] - 1201)/10000.

    # year = start
    # new_times = []

    # months = (np.arange(12))/12.

    # while year <= end:
    #     new_times += [np.array([year for i in range(12)]) + months]
    #     year += 1

    # new_times = np.array(new_times).reshape(-1)
    # # sort the climatology to 1981-2010 (monthly!)

    # # reshape to deal with monthly
    # data = data.reshape(-1, 12, data.shape[-1])
    # new_times = new_times.reshape(-1, 12)

    # # extract climatology period
    # start_loc, = np.where(new_times[:, 0] == CLIMSTART)
    # end_loc, = np.where(new_times[:, 0] == CLIMEND)
    # clim_data = data[start_loc[0]:end_loc[0] + 1, :, :]
    # climatology = np.mean(clim_data, axis=0)

    # # make anomalies
    # data = np.array([data[i, :, :] - climatology for i in range(data.shape[0])])

    # # return to original shapes
    # data = data.reshape(-1, data.shape[-1])
    # new_times = new_times.reshape(-1)

    # utils.plot_hovmuller(settings.IMAGELOC + "LTT_hovmuller_erai", new_times, latitudes, data.T, \
    #                          settings.COLOURMAP_DICT["temperature"], bounds, \
    #                          "Anomaly ("+r'$^{\circ}$'+"C)", cosine=True)


    # version with MEI on top
    # mei = read_mei(os.path.join(DATALOC, "MEI.dat"))


    # utils.plot_hovmuller(settings.IMAGELOC + "LTT_hovmuller_era_MEI", new_times, latitudes, data.T, \
    #                          settings.COLOURMAP_DICT["temperature"], bounds, \
    #                          "Anomaly ("+r'$^{\circ}$'+"C)", cosine=True, extra_ts=mei)



    #************************************************************************
    # Timeseries of record high/low.
    if False:
        cubelist = iris.load(os.path.join(DATALOC, "dunn_records.nc"))

        names = np.array([c.var_name for c in cube_list])

        loc, = np.where(names == "rss_cold")[0]
        times = get_times(cube_list[loc].coord("time"))
        low1 = utils.Timeseries("Coldest", times, cube_list[loc].data)
        low1.lw = 1

        loc, = np.where(names == "rss_hot")[0]
        times = get_times(cube_list[loc].coord("time"))
        high1 = utils.Timeseries("Warmest", times, cube_list[loc].data)
        high1.lw = 1

        loc, = np.where(names == "uah_cold")[0]
        times = get_times(cube_list[loc].coord("time"))
        low2 = utils.Timeseries("Coldest", times, cube_list[loc].data)

        loc, = np.where(names == "uah_hot")[0]
        times = get_times(cube_list[loc].coord("time"))
        high2 = utils.Timeseries("Warmest", times, cube_list[loc].data)

        plt.clf()
        fig, (ax1) = plt.subplots(1, figsize=(8, 6))
                
        utils.plot_ts_panel(ax1, [high1, low1], "-", "temperature", loc="upper right")
        ax1.fill_between(high1.times, high1.data, high2.data, color="r")
        ax1.fill_between(low1.times, low1.data, low2.data, color="b")

        # manually plot to avoid getting legend entries
        ax1.plot(high2.times, high2.data, c=settings.COLOURS["temperature"][high2.name], lw=1)
        ax1.plot(low2.times, low2.data, c=settings.COLOURS["temperature"][low2.name], lw=1)


        for tick in ax1.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        for tick in ax1.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        ax1.set_ylabel("Global Area (%)", fontsize=settings.LABEL_FONTSIZE)
        fig.subplots_adjust(right=0.98, top=0.98, bottom=0.05, hspace=0.001)

        plt.savefig(settings.IMAGELOC+"LTT_records{}".format(settings.OUTFMT))

        plt.close()

    #************************************************************************
    # MEI/N3.4
    if False:

        cubelist = iris.load(os.path.join(DATALOC, "dunn_records.nc"))

        names = np.array([c.var_name for c in cube_list])

        loc, = np.where(names == "n34")[0]
        times = get_times(cube_list[loc].coord("time"))
        n34 = utils.Timeseries(u"Niño 3.4", times, cube_list[loc].data)

        fig, ax1 = plt.subplots(1, figsize=(8, 6))

#        ax1.plot(n34.times, n34.data, c="k", ls="-", lw=2, label="Ni\~no 3.4")
        ax1.fill_between(n34.times, n34.data, color="b", where=n34.data<0)
        ax1.fill_between(n34.times, n34.data, color="r", where=n34.data>0)

#        ax1.legend(loc="upper right", ncol=2, frameon=False, prop={'size':settings.FONTSIZE})
        ax1.set_xlim([1978, int(settings.YEAR)+2])
        ax1.set_ylabel(u"Niño 3.4", fontsize=settings.FONTSIZE)

        for tick in ax1.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        for tick in ax1.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        utils.thicken_panel_border(ax1)
        fig.subplots_adjust(right=0.98, top=0.98, bottom=0.05, hspace=0.001)

        plt.savefig(settings.IMAGELOC+"LTT_n34_ts{}".format(settings.OUTFMT))

        plt.close()

    #************************************************************************
    # MEI & records combined
    if True:

        plt.clf()
        fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 8), sharex=True)

        cubelist = iris.load(os.path.join(DATALOC, "dunn_records.nc"))

        names = np.array([c.var_name for c in cube_list])

        # N3.4
        loc, = np.where(names == "n34")[0]
        times = get_times(cube_list[loc].coord("time"))
        n34 = utils.Timeseries(u"Niño 3.4", times, cube_list[loc].data)

        ax1.fill_between(n34.times, n34.data, color="b", where=n34.data<0)
        ax1.fill_between(n34.times, n34.data, color="r", where=n34.data>0)

        ax1.set_xlim([1978, int(settings.YEAR)+2])
        ax1.set_ylabel(u"Niño 3.4 ("+r'$^{\circ}$'+"C)", fontsize=settings.FONTSIZE)

        for tick in ax1.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        for tick in ax1.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        utils.thicken_panel_border(ax1)
        ax1.text(0.02, 0.9, "(a)", transform=ax1.transAxes, fontsize=settings.FONTSIZE)
        # Records

        loc, = np.where(names == "rss_cold")[0]
        times = get_times(cube_list[loc].coord("time"))
        low1 = utils.Timeseries("Coldest", times, cube_list[loc].data)
        low1.lw = 1

        loc, = np.where(names == "rss_hot")[0]
        times = get_times(cube_list[loc].coord("time"))
        high1 = utils.Timeseries("Warmest", times, cube_list[loc].data)
        high1.lw = 1

        loc, = np.where(names == "uah_cold")[0]
        times = get_times(cube_list[loc].coord("time"))
        low2 = utils.Timeseries("Coldest", times, cube_list[loc].data)

        loc, = np.where(names == "uah_hot")[0]
        times = get_times(cube_list[loc].coord("time"))
        high2 = utils.Timeseries("Warmest", times, cube_list[loc].data)
              
        utils.plot_ts_panel(ax2, [high1, low1], "-", "temperature", loc="upper right")
        ax2.fill_between(high1.times, high1.data, high2.data, color="r")
        ax2.fill_between(low1.times, low1.data, low2.data, color="b")

        # manually plot to avoid getting legend entries
        ax2.plot(high2.times, high2.data, c=settings.COLOURS["temperature"][high2.name], lw=1)
        ax2.plot(low2.times, low2.data, c=settings.COLOURS["temperature"][low2.name], lw=1)


        for tick in ax2.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        for tick in ax2.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        ax2.set_ylabel("Global Area (%)", fontsize=settings.FONTSIZE)
        ax2.text(0.02, 0.9, "(b)", transform=ax2.transAxes, fontsize=settings.FONTSIZE)


        fig.subplots_adjust(right=0.98, top=0.98, bottom=0.05, hspace=0.001)

        plt.savefig(settings.IMAGELOC+"LTT_n34_records_ts{}".format(settings.OUTFMT))

        plt.close()
   
    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
