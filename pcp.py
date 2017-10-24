#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for Precipitation (PCP) section.
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
from mpl_toolkits.basemap import Basemap

import matplotlib.cm as mpl_cm
import matplotlib as mpl

import iris
import iris.quickplot as qplt
import cartopy.crs as ccrs

import utils # RJHD utilities
import settings


data_loc = "/data/local/rdunn/SotC/{}/data/PCP/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

LEGEND_LOC = 'upper left'
BBOX = (0.3,1.0)

#************************************************************************
def read_land(filename):

    indata = np.genfromtxt(filename, dtype=(float), skip_header = 1)

    indata = np.ma.masked_where(indata == -999.9999, indata)

    ghcn = utils.Timeseries("GHCN", indata[:,0], indata[:,1])
    gpcc = utils.Timeseries("GPCC", indata[:,0], indata[:,2])
    gpcp = utils.Timeseries("GPCPv23", indata[:,0], indata[:,3])
    ghcn2 = utils.Timeseries("GHCNv2", indata[:,0], indata[:,4])

    return ghcn, gpcc, gpcp, ghcn2 # read_land

def read_map(filename):

    resolution = 2.5

    lats = np.arange(88.75, -88.75 - resolution, -resolution)
    lons = np.arange(-178.75, 178.75 + resolution, resolution)

    data = np.zeros((len(lats), len(lons)))

    lat_ctr = 0
    this_lat = []
    with open(filename, 'r') as infile:
     
        for ll, line in enumerate(infile):
            for ls in line.split():
                this_lat.append(ls)

            if len(this_lat) == len(lons):
                # read in all for one longitude
                data[lat_ctr, :] = this_lat
                lat_ctr += 1
                
                this_lat = []

    data = np.ma.masked_where(data <= -99999.99, data)

    cube = utils.make_iris_cube_2d(data, lats, lons, "PCP_anom", "mm")

    return cube


#************************************************************************
def run_all_plots():


    #************************************************************************
    # Precipitation Timeseries

    ghcn, gpcc, gpcp, ghcn2 = read_land(data_loc + "Fig2.1h_land_in_situ.dat")

    fig = plt.figure(figsize = (10,6))

    ax1 = plt.axes([0.1,0.1,0.85,0.85])

    # Land
    utils.plot_ts_panel(ax1, [ghcn, gpcc, gpcp, ghcn2], "-", "hydrological", loc = LEGEND_LOC, bbox = BBOX)

    ax1.set_ylim([-60,100])
    ax1.yaxis.set_ticks([-50,0,50,100])
    ax1.set_ylabel("Anomaly (mm)", fontsize = settings.FONTSIZE)

    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 

    # ax1.text(0.02, 0.9, "(a) Land in Situ", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE)

    plt.savefig(image_loc+"PCP_ts{}".format(settings.OUTFMT))
    plt.close()

    # old 3-panel plot for 2015 report

    # fig = plt.figure(figsize = (10,10))

    # ax1 = plt.axes([0.13,0.68,0.84,0.3])
    # ax2 = plt.axes([0.13,0.38,0.84,0.3], sharex = ax1)
    # ax3 = plt.axes([0.13,0.04,0.84,0.3])

    # # Land
    # utils.plot_ts_panel(ax1, [ghcn, gpcc, gpcp, ghcn2], "-", "hydrological", loc = LEGEND_LOC, bbox = BBOX)

    # # Ocean
    # cube = iris.load(data_loc + "ocean_precipitation_sotc2015.nc", "anomaly_time_series")[0]
    # time = iris.load(data_loc + "ocean_precipitation_sotc2015.nc", "time")[0]

    # ocean_month = utils.Timeseries("GPCP", time.data, cube.data)

    # # and make the annual averages
    # annuals = cube.data.reshape(-1,12)
    # annuals = np.mean(annuals, axis = 1)

    # ocean_year = utils.Timeseries("GPCP", np.arange(int(time.data[0]), int(time.data[-1])+1), annuals)

    # utils.plot_ts_panel(ax2, [ocean_year], "-", "hydrological", loc = LEGEND_LOC, bbox = BBOX)


    # utils.plot_ts_panel(ax3, [ocean_month], "-", "hydrological", loc = LEGEND_LOC, bbox = (0.5,1.0))

    # #*******************
    # # prettify

    # ax1.set_ylim([-60,100])
    # ax1.yaxis.set_ticks([-50,0,50,100])
    # ax2.set_ylim([-60,100])
    # ax2.set_ylabel("Anomaly (mm yr"+r'$^{-1}$'+")", fontsize = settings.FONTSIZE)
    # ax2.yaxis.set_ticks([-50,0,50,100])
    # ax3.set_ylim([-140,150])
    # ax3.yaxis.set_ticks([-100,-50,0,50,100,150])
    # ax3.set_xlim([1988,int(settings.YEAR)])

    # for ax in [ax1, ax2, ax3]:
    #     for tick in ax.yaxis.get_major_ticks():
    #         tick.label.set_fontsize(settings.FONTSIZE) 
    #     for tick in ax.xaxis.get_major_ticks():
    #         tick.label.set_fontsize(settings.FONTSIZE) 

    # plt.setp(ax1.get_xticklabels(), visible=False)

    # # sort labelling
    # ax1.text(0.02, 0.9, "(a) Land in Situ", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE)
    # ax2.text(0.02, 0.9, "(b) Ocean Satellite", transform = ax2.transAxes, fontsize = settings.LABEL_FONTSIZE)
    # ax3.text(0.02, 0.9, "(c) Ocean Satellite Expanded", transform = ax3.transAxes, fontsize = settings.LABEL_FONTSIZE)


    # plt.savefig(image_loc+"PCP_ts{}".format(settings.OUTFMT))
    # plt.close()

    #************************************************************************
    # GPCP map

    cube = read_map(data_loc + "GPCP_{}anomaly_base_1981-2000.txt".format(settings.YEAR))
    bounds=[-2000, -400, -300, -200, -100, 0, 100, 200, 300, 400, 2000]

    utils.plot_smooth_map_iris(image_loc + "PCP_{}_anoms_gpcp".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2000 (mm)")
    utils.plot_smooth_map_iris(image_loc + "p2.1_PCP_{}_anoms_gpcp".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2000 (mm)", figtext = "(i) Precipitation")

    #************************************************************************
    # GPCP map - Ocean only


    # # Read Oceans
    # print "FIX YEAR"
    # cube = iris.load(data_loc + "ocean_precipitation_sotc2015.nc", "anomaly_map")[0]

    # # add in the coordinates which haven't been stored sensibly
    # for c, coord in enumerate(["latitude", "longitude"]):
    #     coord_cube = iris.load(data_loc + "ocean_precipitation_sotc2015.nc", coord)[0]
    #     if c == 0:
    #         iris_coord = iris.coords.DimCoord(coord_cube.data[:,0], standard_name=coord, units='degrees')
    #     elif c == 1:
    #         iris_coord = iris.coords.DimCoord(coord_cube.data[0], standard_name=coord, units='degrees')

    #     cube.add_dim_coord(iris_coord,c)
    #     cube.coord(coord).guess_bounds()

    # bounds=[-2000, -400, -300, -200, -100, 0, 100, 200, 300, 400, 2000]

    # utils.plot_smooth_map_iris(image_loc + "PCP_{}_anoms_gpcp_ocean".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2000 (mm)")


    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 End
#************************************************************************
