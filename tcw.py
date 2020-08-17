#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for total column water vapour (TCW) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 29                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2020-08-05 12:12:39 +0100 (Wed, 05 Aug #$:  Date of last commit
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

DATALOC = "{}/{}/data/TCW/".format(settings.ROOTLOC, settings.YEAR)

LEGEND_LOC = 'upper left'
BBOX = (0.05, 0.9)
YLIM = [-1, 1.6]

#************************************************************************
def read_csv(filename, domain="L"):
    """
    Read user supplied CSV for LST into Timeseries object
    """

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
def read_ncdf_ts(filename, domain="L"):
    """
    Read user supplied CSV for LST into Timeseries object
    """

    cube_list = iris.load(filename)

    names = np.array([c.var_name for c in cube_list])

    loc, = np.where(names == "Time")[0]
    times = cube_list[loc].data

    for name in names:
        # skip wrong domains
        if domain == "L" and "Land" not in name:
            continue
        elif domain == "O" and "Ocean" not in name:
            continue

        loc, = np.where(names == name)[0]
        cube = cube_list[loc]

        # remove NANs
        cube.data = np.ma.masked_where(cube.data != cube.data, cube.data)

        if "JRA" in name:
            jra = utils.Timeseries("JRA-55", times, cube.data)
        elif "MERRA" in name:
            merra2 = utils.Timeseries("MERRA-2", times, cube.data)
        elif "ERA" in name:
            era5 = utils.Timeseries("ERA5", times, cube.data)
        elif "COSMIC" in name:
             cosmic = utils.Timeseries("COSMIC RO", times, cube.data)
        elif "GNSS" in name:
             obs = utils.Timeseries("GNSS", times, cube.data)
        elif "RAD" in name:
             obs = utils.Timeseries("RSS", times, cube.data)
        
    return merra2, era5, jra, cosmic, obs  # read_ncdf_ts


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
def read_ncdf_hovmuller(filename):

    START = 1900

    cube_list = iris.load(filename)

    for cube in cube_list:

        if cube.var_name == "Latitude":
            latitudes = cube.data
        elif cube.var_name == "Time":
            times = cube.data
            # months since 1900
            times = times/12.
            times += START
        else:
            data = cube.data.T

    return times, latitudes, data # read_ncdf_hovmuller

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
def read_ncdf_map_data(filename, variable):
    '''
    Read data for maps and convert to cube.  Given as single list files

    :param str filename: file to read

    :returns: cube
    '''
    cube = iris.load(filename, "{}_TPW_anomaly_map_{}".format(variable, settings.YEAR))[0]

    lons = iris.load(filename, "Longitude")[0].data
    lats = iris.load(filename, "Latitude")[0].data

    latcoord = iris.coords.DimCoord(lats, standard_name='latitude', units='degrees')
    loncoord = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
 
    cube.add_dim_coord(latcoord, 0)
    cube.add_dim_coord(loncoord, 1)

    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()

    return cube # read_ncdf_map_data

#************************************************************************
def read_ncdf_scatter_data(filename):
    '''
    Read data for maps and return data, latitudes and longitudes

    :param str filename: file to read

    :returns: anomalies, latitudes, longitudes
    '''

    lons = iris.load(filename, "GNSS_TPW_longitude")[0].data
    lats = iris.load(filename, "GNSS_TPW_latitude")[0].data
    anoms = iris.load(filename, "GNSS_TPW_anomaly_{}".format(settings.YEAR))[0].data

    return anoms, lats, lons

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

# #************************************************************************
def read_COSMIC_late(filename, era, annual=True, do_era=False):

    all_data = np.genfromtxt(filename, dtype=(float))

    all_data = np.ma.masked_where(all_data == -99.9, all_data)

    # convert the years and months to give decimals
    all_years = all_data[:, 0]
    months = all_data[:, 1]
    all_times = all_years + (months - 1)/12.

    if annual:
        years = all_years[::12]
        yearly = utils.annual_average(all_data[:, 2])

        if do_era:
            locs = np.in1d(era.times, years)
            era_mean = np.ma.mean(era.data[locs])
            input(era_mean)
            return utils.Timeseries("COSMIC RO", years, yearly+era_mean) 
        return utils.Timeseries("COSMIC RO", years, yearly) 
    else:
        return utils.Timeseries("COSMIC RO", all_times, all_data[:, 2]) # read_COSMIC_late

#************************************************************************
def run_all_plots():

    #************************************************************************
    # Total Column Water figure
    if True:
        # land
#        merra2_land, erai_land, era5_land, jra_land, cosmic_land, gnss_land=read_csv(DATALOC + "time_series_tpw_land.txt", domain="L")
        merra2_land, era5_land, jra_land, cosmic_land, gnss_land = read_ncdf_ts(DATALOC + "TPW_2019_anom_TS.v2.nc", domain="L")
        dummy, dummy, dummy, cosmic_land, dummy = read_ncdf_ts(DATALOC + "TPW_2019_anom_TS.COSMIC_WRONG.nc", domain="L")
        gnss_land.name = "GNSS (Ground Based)"
        

        # ocean
#        merra2_ocean, erai_ocean, era5_ocean, jra_ocean, cosmic_ocean, radiometer_ocean = read_csv(DATALOC + "time_series_tpw_ocean.txt", domain="O")
        merra2_ocean, era5_ocean, jra_ocean, cosmic_ocean, radiometer_ocean = read_ncdf_ts(DATALOC + "TPW_2019_anom_TS.v2.nc", domain="O")
        dummy, dummy, dummy, cosmic_ocean, dummy = read_ncdf_ts(DATALOC + "TPW_2019_anom_TS.COSMIC_WRONG.nc", domain="O")
        radiometer_ocean.name = "RSS Satellite"


        # updated file June 2020
        print("remove this from June 2020")
        alldata = np.genfromtxt(DATALOC + "TCWV_2020_ts_updated.dat")
        alldata = np.ma.masked_where(alldata == 0.0, alldata)
        cosmic_land = utils.Timeseries("COSMIC RO", alldata[:, 0], alldata[:, 7])
        cosmic_ocean = utils.Timeseries("COSMIC RO", alldata[:, 0], alldata[:, 2])


        #************************************************************************
        # Timeseries figures
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(8, 10), sharex=True)

        # Obs - ocean
        utils.plot_ts_panel(ax1, [radiometer_ocean, cosmic_ocean], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)

        # Reanalyses - ocean
        utils.plot_ts_panel(ax2, [era5_ocean, jra_ocean, merra2_ocean], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)

        # Obs - land
        utils.plot_ts_panel(ax3, [gnss_land, cosmic_land], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)

        # Reanalyses - land
        utils.plot_ts_panel(ax4, [era5_land, jra_land, merra2_land], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)


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

        plt.savefig(settings.IMAGELOC+"TCW_ts_v3{}".format(settings.OUTFMT))
        plt.close()
        
        if False:
            with open("TCW_2020_ts.dat", "w") as outfile:
                outfile.write("year \t RadiometerO \t CosmicO \t ERA5O   \t JRA555O \t MERRA2O \t GNSSL   \t CosmicL \t ERA5L   \t JRA55L \t MERRA2L\n")
                for t, time in enumerate(era5_ocean.times):
                    outfile.write("{} \t {:6.4f} \t {:6.4f} \t {:6.4f} \t {:6.4f} \t {:6.4f} \t {:6.4f} \t {:6.4f} \t {:6.4f} \t {:6.4f} \t {:6.4f}\n".format(time, radiometer_ocean.data[t], cosmic_ocean.data[t], era5_ocean.data[t], jra_ocean.data[t], merra2_ocean.data[t], gnss_land.data[t], cosmic_land.data[t], era5_land.data[t], jra_land.data[t], merra2_land.data[t]))

 
    #************************************************************************
    # JRA Hovmuller figure
    if False:
        times, latitudes, data = read_hovmuller(DATALOC + "data_for_tpw_hofmueller_{}_JRA-55.txt".format(settings.YEAR))
        
        bounds = np.array([-100, -6, -3, -1.5, -0.5, 0, 0.5, 1.5, 3, 6, 100])
        
        utils.plot_hovmuller(settings.IMAGELOC + "TCW_hovmuller_jra", times, latitudes, data.T, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomaly (mm)")

    #************************************************************************
    # ERA5 Hovmuller figure
    if True:
#        times, latitudes, data = read_hovmuller(DATALOC + "data_for_tpw_hofmueller_{}_ERA5.txt".format(settings.YEAR))
        times, latitudes, data = read_ncdf_hovmuller(DATALOC + "TPW_2019_anom_time_lat.nc")

        bounds = np.array([-100, -6, -3, -1.5, -0.5, 0, 0.5, 1.5, 3, 6, 100])

        utils.plot_hovmuller(settings.IMAGELOC + "TCW_hovmuller_era5", times, latitudes, data.T, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomaly (mm)")

    #************************************************************************
    # Microwave, Cosmic and GNSS figure - plate 2.1
    bounds = np.array([-100, -6, -3, -1.5, -0.5, 0, 0.5, 1.5, 3, 6, 100])
    if True:
#        cosmic = read_map_data(DATALOC + "rss_cosmic_data_for_tpw_map_{}.1981_2010.txt".format(settings.YEAR))
        cosmic = read_ncdf_map_data(DATALOC + "TPW_2019_anom_maps.nc", "RSS_COSMIC")

#        gnss_anoms, gnss_lats, gnss_lons = read_scatter_data(DATALOC + "GNSS_ground_stations_tpw_{}.1981_2010.txt".format(settings.YEAR))
        gnss_anoms, gnss_lats, gnss_lons = read_ncdf_scatter_data(DATALOC + "TPW_2019_anom_maps.nc")

        # March 2020 - updated file
        all_gnss = np.genfromtxt(DATALOC + "GNSS_TCWV_anomaly_2019_ref_2006_2014.txt")
        gnss_lats = all_gnss[:, 1].astype(float)
        gnss_lons = all_gnss[:, 2].astype(float)
        gnss_anoms = all_gnss[:, 3].astype(float)
        

        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_TCW_{}_anoms_cosmic".format(settings.YEAR), cosmic, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (mm)", figtext="(r) Total Column Water Vapour", scatter=(gnss_lons, gnss_lats, gnss_anoms))
        utils.plot_smooth_map_iris(settings.IMAGELOC + "TCW_{}_anoms_cosmic".format(settings.YEAR), cosmic, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (mm)", scatter=(gnss_lons, gnss_lats, gnss_anoms))


    #************************************************************************
    # ERA and GNSS figure
    if True:
#        era5 = read_map_data(DATALOC + "ERA5_data_for_tpw_map_{}.1981_2010.txt".format(settings.YEAR))
        era5 = read_ncdf_map_data(DATALOC + "TPW_2019_anom_maps.nc", "ERA5")

        utils.plot_smooth_map_iris(settings.IMAGELOC + "TCW_{}_anoms_era5".format(settings.YEAR), era5, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (mm)", scatter=(gnss_lons, gnss_lats, gnss_anoms))
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_TCW_{}_anoms_era5".format(settings.YEAR), era5, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2010 (mm)", figtext="(i) Total Column Water Vapour", scatter=(gnss_lons, gnss_lats, gnss_anoms))


    #************************************************************************
    # JRA-55 1997 and 2015 figure

    # plotyears = ["1997", "2015"]

    # cubelist = []
    # for year in plotyears:

    #     cube = read_map_data(DATALOC + "data_for_tpw_map_jra55_late_{}.txt".format(year))

    #     cubelist += [cube]

    # # pass to plotting routine
    # bounds = np.array([-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100])

    # utils.plot_smooth_map_iris_multipanel(settings.IMAGELOC + "TCW_{}_year_jra".format(settings.YEAR), cubelist, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomaly (mm)", shape=(2,1), title=plotyears, figtext=["(a)","(b)"])


    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()
