#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for total column water vapour (TCW) section.
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
# python3
from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.cm as mpl_cm
import matplotlib as mpl

import iris
import iris.quickplot as qplt
import cartopy.crs as ccrs

import struct

import utils # RJHD utilities
import settings

data_loc = "{}/{}/data/TCW/".format(settings.ROOTLOC, settings.YEAR)
reanalysis_loc = "{}/{}/data/RNL/".format(settings.ROOTLOC, settings.YEAR)
image_loc = "{}/{}/images/".format(settings.ROOTLOC, settings.YEAR)

LEGEND_LOC = 'upper left'
BBOX = (0, 0.9)
YLIM = [-1, 1.6]

#************************************************************************
def read_csv(filename, domain="L"):
    """
    Read user supplied CSV for LST into Timeseries object
    """

    # fieldwidths = (5, 8, 12, 10, 12, 14, 10, 8)
    # fmtstring = ''.join('%ds' % f for f in fieldwidths)

    # parse = struct.Struct(fmtstring).unpack_from

    # indata = []
    # with open(filename, 'r') as infile:
    #     for ll, line in enumerate(infile):

    #         if ll > 3: # skip the first
    #             fields = parse(line)
                
    #             indata += [[float(f.strip()) for f in fields]]

    indata = np.genfromtxt(filename, skip_header=2, dtype=str, encoding="latin-1")
            
    indata = np.ma.array(indata)
    indata[indata =="NaN"] = "-99.9"
    indata = indata.astype(float)
    indata = np.ma.masked_where(indata == -99.9, indata)

    # process the years and months to give decimals
    times = indata[:, 0]

    if domain == "L":
        jra = utils.Timeseries("JRA-55", times, indata[:, 1])
        merra2 = utils.Timeseries("MERRA-2", times, indata[:, 2])
        erai = utils.Timeseries("ERA-Interim", times, indata[:, 3])
        era5 = utils.Timeseries("ERA5", times, indata[:, 4])
        cosmic = utils.Timeseries("COSMIC RO", times, indata[:, 5])
        obs = utils.Timeseries("GNSS", times, indata[:, 6])

    elif domain == "O":
        obs = utils.Timeseries("RSS", times, indata[:, 1])
        jra = utils.Timeseries("JRA-55", times, indata[:, 2])
        merra2 = utils.Timeseries("MERRA-2", times, indata[:, 3])
        erai = utils.Timeseries("ERA-Interim", times, indata[:, 4])
        era5 = utils.Timeseries("ERA5", times, indata[:, 5])
        cosmic = utils.Timeseries("COSMIC RO", times, indata[:, 6])

    return merra2, erai, era5, jra, cosmic, obs  # read_csv



#************************************************************************
def read_hovmuller(filename):

    all_data = np.genfromtxt(filename, dtype=(float))

    # convert the years and months to give decimals
    years = all_data[:, 0]
    months = all_data[:, 1]
    all_times = years + (months - 1)/12.

    all_latitudes = all_data[:, 2]
    all_data = np.ma.masked_where(all_data[:, 3] == -99.99, all_data[:, 3])
    
    # have columnar data = convert into an array

    times = np.unique(all_times)
    latitudes = np.unique(all_latitudes)

    data = np.ma.zeros((len(times), len(latitudes)))
    data.mask = np.ones(data.shape)

    for v, val in enumerate(all_data):

        xloc, = np.where(times == all_times[v])
        yloc, = np.where(latitudes == all_latitudes[v])

        data[xloc[0], yloc[0]] = val     
        

    return times, latitudes, data # read_hovmuller


#************************************************************************
def read_map_data(filename):
    '''
    Read data for maps and convert to cube.  Given as single list files

    :param str filename: file to read

    :returns: cube
    '''
    

    anoms, lats, lons = read_scatter_data(filename)

    longitudes = np.unique(lons)
    latitudes = np.unique(lats)

    data = np.ma.zeros((len(latitudes), len(longitudes)))
    data.mask = np.ones(data.shape)
    
    for v, val in enumerate(anoms):

        xloc, = np.where(longitudes == lons[v])
        yloc, = np.where(latitudes == lats[v])

        data[yloc[0], xloc[0]] = val     

    cube = utils.make_iris_cube_2d(data, latitudes, longitudes, "TCW_anom", "mm")

    return cube

#************************************************************************
def read_scatter_data(filename):
    '''
    Read data for maps and return data, latitudes and longitudes

    :param str filename: file to read

    :returns: anomalies, latitudes, longitudes
    '''

    all_jra = np.genfromtxt(filename, dtype=(float), skip_header=2)

    lats = all_jra[:, 0]
    lons = all_jra[:, 1]
    anoms = all_jra[:, 2]

    return anoms, lats, lons

#************************************************************************
def run_all_plots():

    #************************************************************************
    # Total Column Water figure

    # land
    merra2_land, erai_land, era5_land, jra_land, cosmic_land, gnss_land=read_csv(data_loc + "time_series_tpw_land.txt", domain="L")
    gnss_land.name = "GNSS (Ground Based)"
    erai_land.ls = "--"

    # ocean
    merra2_ocean, erai_ocean, era5_ocean, jra_ocean, cosmic_ocean, radiometer_ocean = read_csv(data_loc + "time_series_tpw_ocean.txt", domain="O")
    radiometer_ocean.name = "RSS Satellite"
    erai_ocean.ls = "--"


    #************************************************************************
    # Timeseries figures
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(10, 12), sharex=True)

    # Obs - ocean
    utils.plot_ts_panel(ax1, [radiometer_ocean, cosmic_ocean], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)

    # Reanalyses - ocean
    utils.plot_ts_panel(ax2, [erai_ocean, era5_ocean, jra_ocean, merra2_ocean], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)

    # Obs - land
    utils.plot_ts_panel(ax3, [gnss_land, cosmic_land], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)

    # Reanalyses - land
    utils.plot_ts_panel(ax4, [erai_land, era5_land, jra_land, merra2_land], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)


    # prettify
    fig.subplots_adjust(right=0.95, top=0.95, hspace=0.001)

    for tick in ax4.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 
    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_ylim(YLIM)
        ax.yaxis.set_ticks([-1, -0.5, 0, 0.5, 1])
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
    plt.xlim([1979,int(settings.YEAR)+1])

    # sort labelling
    ax1.text(0.02, 0.88, "(a) Observations Ocean", transform=ax1.transAxes, fontsize=settings.LABEL_FONTSIZE)
    ax2.text(0.02, 0.88, "(b) Reanalyses Ocean", transform=ax2.transAxes, fontsize=settings.LABEL_FONTSIZE)
    ax3.text(0.02, 0.88, "(c) Observations Land", transform=ax3.transAxes, fontsize=settings.LABEL_FONTSIZE)
    ax4.text(0.02, 0.88, "(d) Reanalyses Land", transform=ax4.transAxes, fontsize=settings.LABEL_FONTSIZE)

    fig.text(0.03, 0.5, "Anomalies (mm)", va='center', rotation='vertical', fontsize=settings.FONTSIZE)

    plt.savefig(image_loc+"TCW_ts{}".format(settings.OUTFMT))

    plt.close()


    #************************************************************************
    # JRA Hovmuller figure

    times, latitudes, data = read_hovmuller(data_loc + "data_for_tpw_hofmueller_{}_JRA-55.txt".format(settings.YEAR))

    bounds = np.array([-100, -6, -3, -1.5, -0.5, 0, 0.5, 1.5, 3, 6, 100])

    utils.plot_hovmuller(image_loc + "TCW_hovmuller_jra", times, latitudes, data.T, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomaly (mm)")

    #************************************************************************
    # ERA5 Hovmuller figure

    times, latitudes, data = read_hovmuller(data_loc + "data_for_tpw_hofmueller_{}_ERA5.txt".format(settings.YEAR))

    bounds = np.array([-100, -6, -3, -1.5, -0.5, 0, 0.5, 1.5, 3, 6, 100])

    utils.plot_hovmuller(image_loc + "TCW_hovmuller_era5", times, latitudes, data.T, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomaly (mm)")

    #************************************************************************
    # Microwave, Cosmic and GNSS figure - plate 2.1

    bounds = np.array([-100, -6, -3, -1.5, -0.5, 0, 0.5, 1.5, 3, 6, 100])
    cosmic = read_map_data(data_loc + "rss_cosmic_data_for_tpw_map_{}.1981_2010.txt".format(settings.YEAR))

    gnss_anoms, gnss_lats, gnss_lons = read_scatter_data(data_loc + "GNSS_ground_stations_tpw_{}.1981_2010.txt".format(settings.YEAR))

    utils.plot_smooth_map_iris(image_loc + "p2.1_TCW_{}_anoms_cosmic".format(settings.YEAR), cosmic, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (mm)", figtext="(r) Total Column Water Vapour", scatter=(gnss_lons, gnss_lats, gnss_anoms))
    utils.plot_smooth_map_iris(image_loc + "TCW_{}_anoms_cosmic".format(settings.YEAR), cosmic, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (mm)", scatter=(gnss_lons, gnss_lats, gnss_anoms))


    #************************************************************************
    # ERA and GNSS figure

    era5 = read_map_data(data_loc + "ERA5_data_for_tpw_map_{}.1981_2010.txt".format(settings.YEAR))

    utils.plot_smooth_map_iris(image_loc + "TCW_{}_anoms_era5".format(settings.YEAR), era5, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (mm)", scatter=(gnss_lons, gnss_lats, gnss_anoms))
    utils.plot_smooth_map_iris(image_loc + "p2.1_TCW_{}_anoms_era5".format(settings.YEAR), era5, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (mm)", figtext="(i) Total Column Water Vapour", scatter=(gnss_lons, gnss_lats, gnss_anoms))


    #************************************************************************
    # JRA-55 1997 and 2015 figure

    # plotyears = ["1997", "2015"]

    # cubelist = []
    # for year in plotyears:

    #     cube = read_map_data(data_loc + "data_for_tpw_map_jra55_late_{}.txt".format(year))

    #     cubelist += [cube]

    # # pass to plotting routine
    # bounds = np.array([-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100])

    # utils.plot_smooth_map_iris_multipanel(image_loc + "TCW_{}_year_jra".format(settings.YEAR), cubelist, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomaly (mm)", shape=(2,1), title=plotyears, figtext=["(a)","(b)"])


    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()
