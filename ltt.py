#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for lower troposphere temperature (LTT) section.
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



data_loc = "/data/local/rdunn/SotC/{}/data/LTT/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

CLIMSTART=1981
CLIMEND=2010

LEGEND_LOC = 'lower right'

#************************************************************************
def read_csv(filename):
    """
    Read user supplied CSV for LTT into Timeseries object
    """

    indata = np.genfromtxt(filename, delimiter = ',', dtype=(float), skip_header = 1)
   
    indata = np.ma.masked_where(indata == -99.9, indata)

    raobcore = utils.Timeseries("RAOBCORE v1.5", indata[:,0], indata[:,1])
    rich = utils.Timeseries("RICH v1.5", indata[:,0], indata[:,2])
    ratpac = utils.Timeseries("RATPAC A2", indata[:,0], indata[:,3])
    unsw = utils.Timeseries("UNSW v1.0", indata[:,0], indata[:,4])

    UAH = utils.Timeseries("UAH v6.0", indata[:,0], indata[:,5])
    rss = utils.Timeseries("RSS v3.3", indata[:,0], indata[:,6])

    era = utils.Timeseries("ERA-Interim", indata[:,0], indata[:,7])
    jra = utils.Timeseries("JRA-55", indata[:,0], indata[:,8])
    merra = utils.Timeseries("MERRA-2", indata[:,0], indata[:,9])
    
    return raobcore, rich, ratpac, unsw, UAH, rss, era, jra, merra # read_csv



#************************************************************************
def run_all_plots():

    #************************************************************************
    # Timeseries figure (2 panels)

    for region in ["global"]:

        if region == "global":
            raobcore, rich, ratpac, unsw, UAH, rss, era, merra, jra = read_csv(data_loc + "TLT_GL_{}_RD_v2.csv".format(settings.YEAR))

        plt.clf()
        fig, (ax1, ax2, ax3) = plt.subplots(3, figsize = (10,12), sharex=True)

        # sondes
        utils.plot_ts_panel(ax1, [raobcore, rich, ratpac, unsw], "-", "temperature", loc = LEGEND_LOC)

        # satellites   
        utils.plot_ts_panel(ax2, [UAH, rss], "-", "temperature", loc = LEGEND_LOC)

        ax2.set_ylabel("Anomaly ("+r'$^\circ$'+"C)", fontsize = settings.FONTSIZE)


        # reanalyses
        if region == "global":
            jra_actuals, jra_anoms = utils.read_jra55(reanalysis_loc + "JRA-55_MSUch2LT_global_ts.txt", "temperature")
            merra_actuals, merra_anoms = utils.read_merra_LT_LS(reanalysis_loc + "MERRA2_MSU_Tanom_ann.dat", LT = True)

            utils.plot_ts_panel(ax3, [era, jra_anoms, merra_anoms], "-", "temperature", loc = LEGEND_LOC)
        else:
            utils.plot_ts_panel(ax3, [era], "-", "temperature", loc = LEGEND_LOC)


        # sort formatting
        plt.xlim([era.times[0]-1,era.times[-1]+1])
        ax1.set_ylim([-0.89,0.89])
        ax2.set_ylim([-0.89,0.89])
        ax3.set_ylim([-0.89,0.89])

        for tick in ax3.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 

        for tick in ax1.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
        for tick in ax2.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
        for tick in ax3.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 

        # sort labelling
        ax1.text(0.02, 0.9, "(a) Radiosondes", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE)
        ax2.text(0.02, 0.9, "(b) Satellites", transform = ax2.transAxes, fontsize = settings.LABEL_FONTSIZE)
        ax3.text(0.02, 0.9, "(c) Reanalyses", transform = ax3.transAxes, fontsize = settings.LABEL_FONTSIZE)

        fig.subplots_adjust(right = 0.95, top = 0.95, hspace = 0.001)

        plt.savefig(image_loc+"LTT_ts_{}{}".format(region, settings.OUTFMT))

        plt.close()

    #************************************************************************
    # Read in ERA anomalies

    cube_list = iris.load(reanalysis_loc + "TLT_anvj_moda_ann{}{}-ann19812010.nc".format(settings.YEAR, settings.YEAR))

    cube = cube_list[0]
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()

    bounds=[-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
    cmap=settings.COLOURMAP_DICT["temperature"]
    utils.plot_smooth_map_iris(image_loc + "LTT_{}_anoms_era".format(settings.YEAR), cube[0][0], cmap, bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)")
    utils.plot_smooth_map_iris(image_loc + "p2.1_LTT_{}_anoms_era".format(settings.YEAR), cube[0][0], cmap, bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", figtext = "(e) Lower Tropospheric Temperature")

    #************************************************************************
    # ERA Hovmuller

    times, latitudes, data = utils.era_2dts_read(reanalysis_loc, "ltt")

    bounds = np.array([-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100])

    # sort the time axis
    start = (times[0] - 101)/10000.
    end = (times[-1] - 1201)/10000.

    year = start
    new_times = []

    months = (np.arange(12))/12.

    while year <= end:
        new_times += [np.array([year for i in range (12)]) + months]
        year += 1

    new_times = np.array(new_times).reshape(-1)
    # sort the climatology to 1981-2010 (monthly!)

    # reshape to deal with monthly
    data = data.reshape(-1, 12, data.shape[-1]) 
    new_times = new_times.reshape(-1,12)

    # extract climatology period
    start_loc, = np.where(new_times[:,0] == CLIMSTART)
    end_loc, = np.where(new_times[:,0] == CLIMEND)
    clim_data = data[start_loc[0]:end_loc[0] + 1, :, :]
    climatology = np.mean(clim_data, axis = 0)

    # make anomalies
    data = np.array([data[i,:,:] - climatology for i in range(data.shape[0])])

    # return to original shapes
    data = data.reshape(-1, data.shape[-1])
    new_times = new_times.reshape(-1)

    utils.plot_hovmuller(image_loc + "LTT_hovmuller_era", new_times, latitudes, data.T, settings.COLOURMAP_DICT["temperature"], bounds, "Anomaly ("+r'$^{\circ}$'+"C)")

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
