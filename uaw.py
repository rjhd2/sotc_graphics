#!/usr/local/sci/python
# python3
from __future__ import absolute_import
from __future__ import print_function
#************************************************************************
#
#  Plot figures and output numbers for Upper Air Winds (UAW) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 26                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2019-04-17 15:34:18 +0100 (Wed, 17 Apr #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************

import datetime as dt
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker
import netCDF4 as ncdf
import iris

import utils # RJHD utilities
import settings

data_loc = "{}/{}/data/UAW/".format(settings.ROOTLOC, settings.YEAR)
reanalysis_loc = "{}/{}/data/RNL/".format(settings.ROOTLOC, settings.YEAR)
image_loc = "{}/{}/images/".format(settings.ROOTLOC, settings.YEAR)

CLIM_PERIOD = "8110"
LEGEND_LOC = "upper right"



#************************************************************************
def make_masked_smoothed_ts(name, time, data, points):
    '''
    Make a smoothed timeseries of a masked array.  Only use the unmasked points

    :param str name: name for timesereis
    :param array time: time data
    :param array data: data
    :param int points: for boxcar
    '''

    smoothed = np.ma.zeros(len(data))
    smoothed.mask = data.mask

    goods, = np.where(smoothed.mask == False)
    smoothed[goods] = utils.boxcar(data.compressed(), points)
    
    ts = utils.Timeseries(name, np.ma.array(time, mask=data.mask), smoothed) 

    return ts # make_masked_smoothed_ts

#************************************************************************
def make_masked_annual_ts(name, time, data):
    '''
    Make annual timeseries of a masked monthly array.  Only use the unmasked points

    :param str name: name for timesereis
    :param array time: time data
    :param array data: data
    '''

    data = data.reshape(-1, 12)

    counts = np.ma.count(data, axis=1)

    annual = np.ma.mean(data, axis=1)
    
    annual = np.ma.masked_where(counts < 12, annual)

    years = time.reshape(-1, 12)[:, 0]

    ts = utils.Timeseries(name, np.ma.array(years, mask=annual.mask), annual) 

    return ts # make_masked_smoothed_ts

#************************************************************************
def read_uaw_ts(filename, smooth=False, annual=False):

    # IRIS doesn't like the Conventions attribute
    ncfile = ncdf.Dataset(filename, 'r')
    
    time = ncfile.variables["time"][:]
    # time = YYYYMM
    times = np.array([int(str(t)[:4])+((float(str(t)[4:])-1)/12.)  for t in time])
    
    # smooth using 12 point boxcar
    if smooth:
        points = 12
        merra = make_masked_smoothed_ts("MERRA-2", times, ncfile.variables["MERRA2"][:], points)
#        grasp = make_masked_smoothed_ts("GRASP", times, ncfile.variables["GRASP"][:], points)
        cera20c = make_masked_smoothed_ts("CERA20C", times, ncfile.variables["CERA20C"][:], points)
        jra55 = make_masked_smoothed_ts("JRA-55", times, ncfile.variables["JRA55"][:], points)
        erai = make_masked_smoothed_ts("ERA-Interim", times, ncfile.variables["ERAI"][:], points)
        era5 = make_masked_smoothed_ts("ERA5", times, ncfile.variables["ERA5"][:], points)
   
    if annual:
        merra = make_masked_annual_ts("MERRA-2", times, ncfile.variables["MERRA2"][:])
#        grasp = make_masked_annual_ts("GRASP", times, ncfile.variables["GRASP"][:])
        cera20c = make_masked_annual_ts("CERA20C", times, ncfile.variables["CERA20C"][:])
        jra55 = make_masked_annual_ts("JRA-55", times, ncfile.variables["JRA55"][:])
        erai = make_masked_annual_ts("ERA-Interim", times, ncfile.variables["ERAI"][:])
        era5 = make_masked_annual_ts("ERA5", times, ncfile.variables["ERA5"][:])

    if not smooth and not annual:
        merra = utils.Timeseries("MERRA-2", times, ncfile.variables["MERRA2"][:])
#        grasp = utils.Timeseries("GRASP", times, ncfile.variables["GRASP"][:])
        cera20c = utils.Timeseries("CERA20C", times, ncfile.variables["CERA20C"][:])
        jra55 = utils.Timeseries("JRA-55", times, ncfile.variables["JRA55"][:])
        erai = utils.Timeseries("ERA-Interim", times, ncfile.variables["ERAI"][:])
        era5 = utils.Timeseries("ERA5", times, ncfile.variables["ERA5"][:])
    

    ncfile.close()

    # set ERA-Interim linestyle
    erai.ls = "--"

    return era5, erai, cera20c, merra, jra55 # read_uaw_ts

#************************************************************************
def read_QBO(filename):

    indata = np.genfromtxt(filename, dtype=(float), skip_header=1, missing_values="NA", filling_values=-999.)

    # process the years and months to give decimals 
    years = indata[:, 0]
    months = indata[:, 1]
    times = years + (months - 1)/12.

    pre1960, = np.where(times <= 1960)

    indata = np.ma.masked_where(indata <= -999., indata)

    # levels of 3,4,5,6,8,10,12,15,20,25,30,35,40,45,50,60,70,80,90
    qbo = utils.Timeseries("GIUB", times[pre1960], indata[pre1960, 16])

    return qbo # read_QBO



#************************************************************************
def run_all_plots():

    #************************************************************************
    # Timeseries - 2016
    # grasp, erai, era_presat, merra, jra55 = read_uaw_ts(data_loc + "20N-40N300.nc", smooth = True)
    # qbo = read_QBO(data_loc + "qbo_1908_2015_REC_ERA40_ERAINT.txt")

    # fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, figsize = (10,15), sharex=True)

    # # Observations
    # utils.plot_ts_panel(ax1, [grasp], "-", "circulation", loc = LEGEND_LOC, ncol=2, extra_labels = [" (0.02)"])

    # # Reanalyses
    # utils.plot_ts_panel(ax2, [erai, era_presat, jra55, merra], "-", "circulation", loc = LEGEND_LOC, ncol=2, extra_labels = [" (-0.20)", " (0.33)", " (-0.13)", " (-0.07)"])

    # grasp, erai, era_presat, merra, jra55 = read_uaw_ts(data_loc + "10S-10N50.nc", smooth = False)

    # # Observations
    # utils.plot_ts_panel(ax3, [qbo, grasp], "-", "circulation", loc = LEGEND_LOC, ncol=2, extra_labels = ["", " (-0.31)"])
    # ax3.set_ylabel("Zonal Anomaly (m s"+r'$^{-1}$'+")", fontsize = settings.FONTSIZE)

    # # Reanalyses
    # utils.plot_ts_panel(ax4, [erai, era_presat, jra55, merra], "-", "circulation", loc = LEGEND_LOC, ncol=2, extra_labels = [" (-0.33)"," (-0.16)", " (-0.37)", " (0.30)"])

    # fig.subplots_adjust(left = 0.11, right = 0.99, top = 0.99, hspace = 0.001)

    # # turn on 4th axis ticks

    # for tick in ax4.get_xticklabels():
    #     tick.set_visible(True)

    # # delete the 5th axis and recreate - to break the sharex link
    # fig.delaxes(ax5)
    # ax5 = fig.add_subplot(515)
    # pos = ax5.get_position()
    # new_pos = [pos.x0, pos.y0 - 0.05, pos.width, pos.height]
    # ax5.set_position(new_pos)

    # # Obs & Reanalyses
    # utils.plot_ts_panel(ax5, [grasp, erai, jra55, merra], "-", "circulation", loc = LEGEND_LOC, ncol=2)

    # # sort formatting
    # for ax in [ax4, ax5]:
    #     for tick in ax.xaxis.get_major_ticks():
    #         tick.label.set_fontsize(settings.FONTSIZE) 

    # for ax in [ax1, ax2, ax3, ax4, ax5]:
    #     for tick in ax.yaxis.get_major_ticks():
    #         tick.label.set_fontsize(settings.FONTSIZE) 

    # # x + y limit
    # ax1.set_xlim([1930,2017.9])
    # ax1.set_ylim([-13,6])
    # ax1.yaxis.set_ticks([-10, -5, 0])
    # ax2.set_ylim([-3.8,3.8])
    # ax2.yaxis.set_ticks([-2, 0, 2, 4])
    # ax3.set_ylim([-34,24])
    # ax4.set_ylim([-34,24])
    # ax5.set_ylim([-34,24])

    # ax5.set_xlim([2000,2017.9])

    # # sort labelling
    # ax1.text(0.02, 0.87, "(a) Observations 20"+r'$^\circ$'+" - 40"+r'$^\circ$'+"N 300hPa", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE)
    # ax2.text(0.02, 0.87, "(b) Reanalyses 20"+r'$^\circ$'+" - 40"+r'$^\circ$'+"N 300hPa", transform = ax2.transAxes, fontsize = settings.LABEL_FONTSIZE)
    # ax3.text(0.02, 0.87, "(c) Observations & Reconstructions 10"+r'$^\circ$'+"S - 10"+r'$^\circ$'+"N 50hPa", transform = ax3.transAxes, fontsize = settings.LABEL_FONTSIZE)
    # ax4.text(0.02, 0.87, "(d) Reanalyses 10"+r'$^\circ$'+"S - 10"+r'$^\circ$'+"N 50hPa", transform = ax4.transAxes, fontsize = settings.LABEL_FONTSIZE)
    # ax5.text(0.02, 0.87, "(e) Observations & Reanalyses 10"+r'$^\circ$'+"S - 10"+r'$^\circ$'+"N 50hPa", transform = ax5.transAxes, fontsize = settings.LABEL_FONTSIZE)


    # plt.savefig(image_loc+"UAW_ts{}".format(settings.OUTFMT))

    #************************************************************************
    # Timeseries - 2018

    plt.figure(figsize=(10, 6))
    plt.clf()
    ax = plt.axes([0.12, 0.10, 0.87, 0.87])
 
    # Globe
#    grasp, erai, cera, merra, jra55 = read_uaw_ts(data_loc + "Globe850.nc", annual=True)
    era5, erai, cera, merra, jra55 = read_uaw_ts(data_loc + "Globe850.nc", annual=True)
    utils.plot_ts_panel(ax, [merra, erai, era5, jra55, cera], "-", "circulation", \
                        loc=LEGEND_LOC, ncol=2, extra_labels=[" (0.03)", " (0.07)", \
                                                                  " (0.03)", " (0.06)", ""])

    # sort formatting
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 

    # x + y limit
    ax.set_xlim([1958, int(settings.YEAR)+0.9])
    ax.set_ylim([-0.39, 1.0])
    ax.yaxis.set_ticks_position('left')
    ax.set_ylabel("Wind Anomaly (m s"+r'$^{-1}$'+")", fontsize=settings.LABEL_FONTSIZE)
 
    # # sort labelling
    ax.text(0.02, 0.87, "Globe 850hPa", transform=ax.transAxes, fontsize=settings.LABEL_FONTSIZE)

    plt.savefig(image_loc+"UAW_globe_ts{}".format(settings.OUTFMT))

    # #*******
    # # Tropics timeseries

    # fig = plt.figure(figsize=(10, 6))
    # plt.clf()
    # ax = plt.axes([0.10, 0.10, 0.87, 0.87])
 
    # # 10N to 10S
    # grasp, erai, cera, merra, jra55 = read_uaw_ts(data_loc + "10S-10N50.nc")
    # utils.plot_ts_panel(ax, [merra, erai, jra55, grasp], "-", "circulation",\
    #                         loc=LEGEND_LOC, ncol=2, extra_labels=[" (0.17)", " (-0.40)", \
    #                                                                   " (-0.51)", " (-0.30)"])

    # # sort formatting
    # for tick in ax.xaxis.get_major_ticks():
    #     tick.label.set_fontsize(settings.FONTSIZE) 

    # for tick in ax.yaxis.get_major_ticks():
    #     tick.label.set_fontsize(settings.FONTSIZE) 

    # # x + y limit
    # ax.set_xlim([2000, int(settings.YEAR)+2])
    # ax.set_ylim([-28, 18])
    # ax.yaxis.set_ticks_position('left')
    # ax.set_ylabel("Wind Anomaly (m s"+r'$^{-1}$'+")", fontsize=settings.LABEL_FONTSIZE)
 
    # # # sort labelling
    # ax.text(0.02, 0.87, "10"+r'$^\circ$'+"S - 10"+r'$^\circ$'+"N 50hPa", \
    #             transform=ax.transAxes, fontsize=settings.LABEL_FONTSIZE)

    # plt.savefig(image_loc+"UAW_tropics_ts{}".format(settings.OUTFMT))

    #************************************************************************
    # Global Map - ERA-I Anomaly figure

    # Read in ERA-I anomalies

    # IRIS doesn't like the Conventions attribute
    ncfile = ncdf.Dataset(data_loc + "ERAI_850.nc", 'r')

    var = ncfile.variables["ws"][:] # this is a masked array
    lons = ncfile.variables["longitude"][:]
    lats = ncfile.variables["latitude"][:]

    ncfile.close()

    # monthly data, so take mean
    mean = np.mean(var, axis=0)

    cube = utils.make_iris_cube_2d(mean, lats, lons, "UAW_ANOM", "m/s")

    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]

    utils.plot_smooth_map_iris(image_loc + "p2.1_UAW_{}_anoms_erai".format(settings.YEAR), \
                                   cube, settings.COLOURMAP_DICT["circulation"], bounds, \
                                   "Anomalies from 1981-2010 (m s"+r'$^{-1}$'+")", \
                                   figtext="(w) Upper Air (850-hPa) Winds")
    utils.plot_smooth_map_iris(image_loc + "UAW_{}_anoms_erai".format(settings.YEAR), \
                                   cube, settings.COLOURMAP_DICT["circulation"], bounds, \
                                   "Anomalies from 1981-2010 (m s"+r'$^{-1}$'+")")

    #************************************************************************
    # Global Map - ERA5 Anomaly figure

    # Read in ERA5 anomalies

    # IRIS doesn't like the Conventions attribute
    ncfile = ncdf.Dataset(data_loc + "ERA5_850.nc", 'r')

    var = ncfile.variables["ws"][:] # this is a masked array
    lons = ncfile.variables["longitude"][:]
    lats = ncfile.variables["latitude"][:]

    ncfile.close()

    # monthly data, so take mean
    mean = np.mean(var, axis=0)

    cube = utils.make_iris_cube_2d(mean, lats, lons, "UAW_ANOM", "m/s")

    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]

    utils.plot_smooth_map_iris(image_loc + "p2.1_UAW_{}_anoms_era5".format(settings.YEAR), \
                                   cube, settings.COLOURMAP_DICT["circulation"], bounds, \
                                   "Anomalies from 1981-2010 (m s"+r'$^{-1}$'+")", \
                                   figtext="(w) Upper Air (850-hPa) Winds")
    utils.plot_smooth_map_iris(image_loc + "UAW_{}_anoms_era5".format(settings.YEAR), \
                                   cube, settings.COLOURMAP_DICT["circulation"], bounds, \
                                   "Anomalies from 1981-2010 (m s"+r'$^{-1}$'+")")

    #************************************************************************
    # QBO plot - https://www.geo.fu-berlin.de/en/met/ag/strat/produkte/qbo/index.html
    
    levels = np.array([70., 50., 40., 30., 20., 15., 10.])
    times = []
    dttimes = []
    data = np.zeros((levels.shape[0], 13))
    factor = 0.1
    j = 0
    
    with open(data_loc + "qbo.dat", "r") as infile:

        for line in infile:
            line = line.split()

            if len(line) > 0:
                    
                # get current year
                try:
                    if int(line[1][:2]) >= int(settings.YEAR[-2:]) and int(line[1][:2]) <= int(settings.YEAR[-2:])+1:
                        month = int(line[1][-2:])
                        times += [month]
                        dttimes += [dt.datetime(int(settings.YEAR), month, 1)]
                        data[:, j] = [float(i)*factor for i in line[2:]]
                        j += 1
                except ValueError:
                    pass

    data = np.array(data)
    times = np.array(times)
    times[-1] += 12

    # And now plot
    cmap = settings.COLOURMAP_DICT["circulation"]
    bounds = [-100., -45., -30., -15., -10., -5., 0., 5., 10., 15., 30., 45., 100]
    norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)

    fig = plt.figure(figsize=(8, 8))
    plt.clf()
    ax = plt.axes([0.12, 0.07, 0.8, 0.9])

    times, levels = np.meshgrid(times, levels)

    con = plt.contourf(times, levels, data, bounds, cmap=cmap, norm=norm, vmax=bounds[-1], vmin=bounds[1])

    plt.ylabel("Pressure (hPa)", fontsize=settings.FONTSIZE)
    plt.xlabel(settings.YEAR, fontsize=settings.FONTSIZE)
    plt.xticks(times[0], [dt.datetime.strftime(d, "%b") for d in dttimes], fontsize=settings.FONTSIZE*0.8)

    plt.xlim([1, 13])
    plt.ylim([70, 10])

    ax.set_yscale("log", subsy=[])
    plt.gca().yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(10))
    plt.gca().yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    plt.gca().yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())

    plt.yticks(np.arange(70, 0, -10), ["{}".format(l) for l in np.arange(70, 0, -10)], fontsize=settings.FONTSIZE)

    # colourbar and prettify
    cb = plt.colorbar(con, orientation='horizontal', pad=0.1, fraction=0.05, aspect=30, \
                          ticks=bounds[1:-1], label="zonal wind (m/s)", drawedges=True)

    cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
    cb.ax.tick_params(axis='x', labelsize=settings.FONTSIZE*0.6, direction='in')

    cb.set_label(label="zonal wind (m/s)", fontsize=settings.FONTSIZE*0.6)
#    cb.outline.set_color('k')
    cb.outline.set_linewidth(2)
    cb.dividers.set_color('k')
    cb.dividers.set_linewidth(2)
                
    utils.thicken_panel_border(ax)

    plt.savefig(image_loc+"UAW_QBO_levels{}".format(settings.OUTFMT))    


    #************************************************************************
    # 200hPa winds in 1980 and 2018    

    cube_list = iris.load(data_loc + "ws200_spread_197901_201801.nc")
    names = np.array([str(cube.var_name) for cube in cube_list])

    # hard coded labels
    mu = {"1980":"1.8", "2018": "1.0"}
    rms = {"1980":"2.0", "2018": "1.1"}
    label = {"1980":"(a)", "2018": "(b)"}
    bounds = [0, 0.5, 1, 1.5, 2, 2.5, 3.0, 100]
    for name in names:
        print(name)
        cube_index, = np.where(names == name)
        cube = cube_list[cube_index[0]]
        
        year = name.split("_")[-1][:4]
        
        utils.plot_smooth_map_iris(image_loc + "UAW_200hPa_Jan{}".format(year), cube, plt.cm.BuPu, bounds, "m/s", figtext="{} January {}, mean={}, RMS={}".format(label[year], year, mu[year], rms[year]))
       


    #************************************************************************
    # Plots
    label = {"1980":"(c)", "2018": "(d)"}
    for year in ["1980", "2018"]:

        cube_list = iris.load(data_loc + "v200_zonal_{}01.nc".format(year))

        for cube in cube_list:
            if cube.var_name == "products":
                names = cube
            elif cube.var_name == "v200_array":
                data_array = cube

        latitudes = data_array.coord("latitude").points
        plt.figure()
        ax=plt.axes([0.13, 0.13, 0.85, 0.85])

        COLOURS = {"ERA5 ens. mean": "red", "ERA5 ensemble": "orange", "ERA5 HRES": "orange", "JRA55": "c", "MERRA-2": "m", "ERA-Interim": "orange"}

        for name, data in zip(names.data, data_array.data):

            str_name = "".join(str(name.compressed(), "latin-1").rstrip())
            # manually fix names
            if str_name == "ERAI":
                str_name = "ERA-Interim"
            elif str_name == "ERA5 ens mean":
                str_name = "ERA5 ens. mean"

            if str_name[:3] == "mem":
                plt.plot(latitudes[1:], data[1:], c="orange")
            elif str_name == "ERA-Interim":
                plt.plot(latitudes[1:], data[1:], c=COLOURS[str_name], label=str_name, lw=2, ls="--")
            else:
                plt.plot(latitudes[1:], data[1:], c=COLOURS[str_name], label=str_name, lw=2)

        plt.legend(loc="upper right", ncol=1, frameon=False)
        plt.xlabel("Latitude", fontsize=settings.FONTSIZE*0.8)
        plt.ylabel("m/s", fontsize=settings.FONTSIZE*0.8)
        plt.text(0.03, 0.92, "{} January {}".format(label[year], year), transform=ax.transAxes, fontsize=settings.FONTSIZE*0.8)

        plt.xlim([-90, 90])
        plt.xticks(np.arange(-90, 120, 30))
        plt.ylim([-1, 4])
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE*0.8) 
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE*0.8) 
        utils.thicken_panel_border(ax)

        plt.savefig(image_loc+"UAW_200hPa_Jan{}_ts{}".format(year, settings.OUTFMT))       

    return

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
