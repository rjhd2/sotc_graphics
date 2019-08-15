#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for Precipitation (PCP) section.
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
from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import utils # RJHD utilities
import settings


data_loc = "{}/{}/data/PCP/".format(settings.ROOTLOC, settings.YEAR)
reanalysis_loc = "{}/{}/data/RNL/".format(settings.ROOTLOC, settings.YEAR)
image_loc = "{}/{}/images/".format(settings.ROOTLOC, settings.YEAR)

LEGEND_LOC = 'upper left'
BBOX = (0.3, 1.0)

#************************************************************************
def read_land(filename):

    indata = np.genfromtxt(filename, dtype=(float), skip_header=1)

    indata = np.ma.masked_where(indata == -99.99, indata)

    ghcn = utils.Timeseries("GHCN", indata[:, 0], indata[:, 1])
    ghcn2 = utils.Timeseries("GHCNv2", indata[:, 0], indata[:, 2])
    gpcc = utils.Timeseries("GPCC", indata[:, 0], indata[:, 3])
    gpcp = utils.Timeseries("GPCPv23", indata[:, 0], indata[:, 4])

#    cfsr = utils.Timeseries("CFSR", indata[:, 0], indata[:, 5])
    erai = utils.Timeseries("ERA-Interim", indata[:, 0], indata[:, 5])
    merra = utils.Timeseries("MERRA-2", indata[:, 0], indata[:, 6])

    return ghcn, gpcc, gpcp, ghcn2, erai, merra # read_land

#************************************************************************
def read_domain(filename, name):

    years = np.arange(1979, int(settings.YEAR)+1)

    if name not in ["Nino 3.4"]:
        data = np.genfromtxt(filename, dtype=(float)) * 365. # convert /day to /year
    else:
        data = np.genfromtxt(filename, dtype=(float))

    return utils.Timeseries(name, years, data) # read_domain

#************************************************************************
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

    ghcn, gpcc, gpcp, ghcn2, erai, merra2 = read_land(data_loc + "Land_insitu_timeseries-1979.dat")

    fig = plt.figure(figsize=(10, 6))

    ax1 = plt.axes([0.1, 0.1, 0.85, 0.85])

    # Land
    utils.plot_ts_panel(ax1, [ghcn, gpcc, gpcp, ghcn2, erai, merra2], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)

    ax1.set_ylim([-60, 100])
    ax1.yaxis.set_ticks([-50, -25, 0, 25, 50, 75, 100])
    ax1.set_ylabel("Anomaly (mm)", fontsize=settings.FONTSIZE)

    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 

    # ax1.text(0.02, 0.9, "(a) Land in Situ", transform=ax1.transAxes, fontsize=settings.LABEL_FONTSIZE)

    plt.savefig(image_loc+"PCP_ts{}".format(settings.OUTFMT))
    plt.close()

    # old 3-panel plot for 2015 report

    fig = plt.figure(figsize=(10, 10))

    ax1 = plt.axes([0.11, 0.68, 0.78, 0.3])
    ax2 = plt.axes([0.11, 0.38, 0.78, 0.3],  sharex=ax1)
    ax3 = plt.axes([0.11, 0.08, 0.78, 0.3],  sharex=ax1)

    # Land
    utils.plot_ts_panel(ax1, [ghcn, gpcc, gpcp, ghcn2, erai, merra2], "-", "hydrological", loc=LEGEND_LOC, bbox=BBOX)

    # Ocean
    ocean_year = read_domain(data_loc + "gpcp_v23_globalocean_ann_1979_2018", "GPCP")

    o_clim, o_anoms = utils.calculate_climatology_and_anomalies_1d(ocean_year, 1981, 2010)

    ax2.plot(o_anoms.times, o_anoms.data, c="b", ls="-", label=o_anoms.name, lw=2)

    # Nino and others
    land_year = read_domain(data_loc + "gpcp_v23_globalland_ann_1979_2018", "GPCP")
    l_clim, l_anoms = utils.calculate_climatology_and_anomalies_1d(land_year, 1981, 2010)
    combined_year = read_domain(data_loc + "gpcp_v23_globallandocean_ann_1979_2018", "GPCP")
    c_clim, c_anoms = utils.calculate_climatology_and_anomalies_1d(combined_year, 1981, 2010)

    nino = read_domain(data_loc + "nino34_1979_2018_ann_anomaly", "Nino 3.4")
    
    ax3.plot(o_anoms.times, o_anoms.data, c="b", ls="-", label="{} Ocean".format(o_anoms.name), lw=2, zorder=10)
    ax3.plot(l_anoms.times, l_anoms.data, c="lime", ls="-", label="{} Land".format(l_anoms.name), lw=2, zorder=10)
    ax3.plot(c_anoms.times, c_anoms.data, c="r", ls="-", label="{} Land + Ocean".format(c_anoms.name), lw=2, zorder=10)
    ax3.plot([1960, 1961], [0, 0], c="k", label="Nino 3.4")

    ax4 = ax3.twinx()
    ax4.plot(nino.times, nino.data, c="k", zorder=1)
    ax4.fill_between(nino.times, nino.data, 0, color="0.5", label=nino.name, zorder=1)
    ax4.set_ylim([-1.9, 1.9])

    ax3.patch.set_visible(False)
    ax3.set_zorder(ax4.get_zorder() + 1)

    #*******************
    # prettify
    minorLocator = MultipleLocator(1)
    ax1.set_ylim([-60, 100])
    ax1.yaxis.set_ticks([-50, 0, 50, 100])
    for ax in [ax2, ax3]:
        ax.set_ylim([-50, 50])
        ax.axhline(0, c='0.5', ls='--')
        ax.legend(loc=LEGEND_LOC, ncol=2, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, \
                          labelspacing=0.1, columnspacing=0.5, bbox_to_anchor=BBOX)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_ticks_position('left')
        utils.thicken_panel_border(ax)
    
    ax2.set_ylabel("Anomaly (mm yr"+r'$^{-1}$'+")", fontsize=settings.FONTSIZE)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)


    ax4.set_xlim([1979, int(settings.YEAR)+1])
    ax4.yaxis.set_label_position("right")
    ax3.yaxis.set_tick_params(right=False)
    utils.thicken_panel_border(ax4)

    for ax in [ax1, ax2, ax3]:
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
    for tick in ax3.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 
    for tick in ax4.yaxis.get_major_ticks():
        tick.label2.set_fontsize(settings.FONTSIZE)


    # sort labelling
    ax1.text(0.02, 0.9, "(a) Land in Situ", transform=ax1.transAxes, fontsize=settings.LABEL_FONTSIZE)
    ax2.text(0.02, 0.9, "(b) Ocean", transform=ax2.transAxes, fontsize=settings.LABEL_FONTSIZE)
    ax3.text(0.02, 0.9, "(c) ", transform=ax3.transAxes, fontsize=settings.LABEL_FONTSIZE)


    plt.savefig(image_loc+"PCP_ts_3panel{}".format(settings.OUTFMT))
    plt.close()

    #************************************************************************
    # GPCP map

    cube = read_map(data_loc + "GPCP_{}anomaly_base_1981-2000.txt".format(settings.YEAR))
    bounds = [-2000, -400, -300, -200, -100, 0, 100, 200, 300, 400, 2000]

    utils.plot_smooth_map_iris(image_loc + "PCP_{}_anoms_gpcp".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2000 (mm)")
    utils.plot_smooth_map_iris(image_loc + "p2.1_PCP_{}_anoms_gpcp".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomalies from 1981-2000 (mm)", figtext="(k) Precipitation")

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
