#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for Upper Air Winds (UAW) section.
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
import netCDF4 as ncdf

import matplotlib.cm as mpl_cm
import matplotlib as mpl

import iris
import iris.quickplot as qplt
import cartopy.crs as ccrs

import calendar

import utils # RJHD utilities
import settings

data_loc = "/data/local/rdunn/SotC/{}/data/UAW/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

CLIM_PERIOD="8110"
LEGEND_LOC = "lower left"



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
    
    ts = utils.Timeseries(name, np.ma.array(time, mask = data.mask), smoothed) 

    return ts # make_masked_smoothed_ts


#************************************************************************
def read_uaw_ts(filename, smooth = True):

    # IRIS doesn't like the Conventions attribute
    ncfile=ncdf.Dataset(filename,'r')
    
    time = ncfile.variables["time"][:]
    # time = YYYYMM
    times = np.array([int(str(t)[:4])+(float(str(t)[4:])/12.)  for t in time])
    
    era20c = ncfile.variables["ERA20C"][:] 
    cera20c = ncfile.variables["CERA20C"][:]


    # smooth using 12 point boxcar
    if smooth:
        points = 12
        merra = make_masked_smoothed_ts("MERRA-2", times, ncfile.variables["MERRA2"][:], points)
        grasp = make_masked_smoothed_ts("GRASP", times, ncfile.variables["GRASP"][:], points)
        era_presat = make_masked_smoothed_ts("ERApreSAT", times, ncfile.variables["ERApreSAT"][:], points)
        jra55 = make_masked_smoothed_ts("JRA-55", times, ncfile.variables["JRA55"][:], points)
        erai = make_masked_smoothed_ts("ERA-Interim", times, ncfile.variables["ERAI"][:], points)
    else:
        merra = utils.Timeseries("MERRA-2", times, ncfile.variables["MERRA2"][:])
        grasp = utils.Timeseries("GRASP", times, ncfile.variables["GRASP"][:])
        era_presat = utils.Timeseries("ERApreSAT", times, ncfile.variables["ERApreSAT"][:])
        jra55 = utils.Timeseries("JRA-55", times, ncfile.variables["JRA55"][:])
        erai = utils.Timeseries("ERA-Interim", times, ncfile.variables["ERAI"][:])

    

    ncfile.close()

    return grasp, erai, era_presat, merra, jra55 # read_uaw_ts

#************************************************************************
def read_QBO(filename):

    indata = np.genfromtxt(filename, dtype = (float), skip_header = 1, missing_values = "NA", filling_values = -999.)

    # process the years and months to give decimals 
    years = indata[:,0]
    months = indata[:,1]
    times = years + (months - 1)/12.

    pre1960, = np.where(times <= 1960)

    indata = np.ma.masked_where(indata <= -999., indata)

    # levels of 3,4,5,6,8,10,12,15,20,25,30,35,40,45,50,60,70,80,90
    qbo = utils.Timeseries("GIUB", times[pre1960], indata[pre1960,16])

    return qbo # read_QBO



#************************************************************************
def run_all_plots():

    #************************************************************************
    # Timeseries
    grasp, erai, era_presat, merra, jra55 = read_uaw_ts(data_loc + "20N-40N300.nc", smooth = True)
    qbo = read_QBO(data_loc + "qbo_1908_2015_REC_ERA40_ERAINT.txt")

    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, figsize = (10,15), sharex=True)

    # Observations
    utils.plot_ts_panel(ax1, [grasp], "-", "circulation", loc = LEGEND_LOC, ncol=2, extra_labels = [" (0.02)"])

    # Reanalyses
    utils.plot_ts_panel(ax2, [erai, era_presat, jra55, merra], "-", "circulation", loc = LEGEND_LOC, ncol=2, extra_labels = [" (-0.20)", " (0.33)", " (-0.13)", " (-0.07)"])

    grasp, erai, era_presat, merra, jra55 = read_uaw_ts(data_loc + "10S-10N50.nc", smooth = False)

    # Observations
    utils.plot_ts_panel(ax3, [qbo, grasp], "-", "circulation", loc = LEGEND_LOC, ncol=2, extra_labels = ["", " (-0.31)"])
    ax3.set_ylabel("Zonal Anomaly (m s"+r'$^{-1}$'+")", fontsize = settings.FONTSIZE)

    # Reanalyses
    utils.plot_ts_panel(ax4, [erai, era_presat, jra55, merra], "-", "circulation", loc = LEGEND_LOC, ncol=2, extra_labels = [" (-0.33)"," (-0.16)", " (-0.37)", " (0.30)"])

    fig.subplots_adjust(left = 0.11, right = 0.99, top = 0.99, hspace = 0.001)

    # turn on 4th axis ticks

    for tick in ax4.get_xticklabels():
        tick.set_visible(True)

    # delete the 5th axis and recreate - to break the sharex link
    fig.delaxes(ax5)
    ax5 = fig.add_subplot(515)
    pos = ax5.get_position()
    new_pos = [pos.x0, pos.y0 - 0.05, pos.width, pos.height]
    ax5.set_position(new_pos)

    # Obs & Reanalyses
    utils.plot_ts_panel(ax5, [grasp, erai, jra55, merra], "-", "circulation", loc = LEGEND_LOC, ncol=2)

    # sort formatting
    for ax in [ax4, ax5]:
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 

    for ax in [ax1, ax2, ax3, ax4, ax5]:
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 

    # x + y limit
    ax1.set_xlim([1930,2017.9])
    ax1.set_ylim([-13,6])
    ax1.yaxis.set_ticks([-10, -5, 0])
    ax2.set_ylim([-3.8,3.8])
    ax2.yaxis.set_ticks([-2, 0, 2, 4])
    ax3.set_ylim([-34,24])
    ax4.set_ylim([-34,24])
    ax5.set_ylim([-34,24])

    ax5.set_xlim([2000,2017.9])

    # sort labelling
    ax1.text(0.02, 0.87, "(a) Observations 20"+r'$^\circ$'+" - 40"+r'$^\circ$'+"N 300hPa", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE)
    ax2.text(0.02, 0.87, "(b) Reanalyses 20"+r'$^\circ$'+" - 40"+r'$^\circ$'+"N 300hPa", transform = ax2.transAxes, fontsize = settings.LABEL_FONTSIZE)
    ax3.text(0.02, 0.87, "(c) Observations & Reconstructions 10"+r'$^\circ$'+"S - 10"+r'$^\circ$'+"N 50hPa", transform = ax3.transAxes, fontsize = settings.LABEL_FONTSIZE)
    ax4.text(0.02, 0.87, "(d) Reanalyses 10"+r'$^\circ$'+"S - 10"+r'$^\circ$'+"N 50hPa", transform = ax4.transAxes, fontsize = settings.LABEL_FONTSIZE)
    ax5.text(0.02, 0.87, "(e) Observations & Reanalyses 10"+r'$^\circ$'+"S - 10"+r'$^\circ$'+"N 50hPa", transform = ax5.transAxes, fontsize = settings.LABEL_FONTSIZE)


    plt.savefig(image_loc+"UAW_ts{}".format(settings.OUTFMT))

    #************************************************************************
    # Global Map - ERA Anomaly figure

    # Read in ERA anomalies

    # IRIS doesn't like the Conventions attribute
    ncfile=ncdf.Dataset(data_loc + "ERAI_850.nc",'r')

    var=ncfile.variables["u"][:] # this is a masked array
    lons = ncfile.variables["longitude"][:]
    lats = ncfile.variables["latitude"][:]

    ncfile.close()

    # monthly data, so take mean
    mean = np.ma.mean(var, axis = 0)

    cube = utils.make_iris_cube_2d(mean, lats, lons, "UAW_ANOM", "m/s")

    bounds=[-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]

    utils.plot_smooth_map_iris(image_loc + "p2.1_UAW_{}_anoms_era".format(settings.YEAR), cube, settings.COLOURMAP_DICT["circulation"], bounds, "Anomalies from 1981-2010 (m s"+r'$^{-1}$'+")", figtext = "(t) Upper Air (850-hPa) Winds")
    utils.plot_smooth_map_iris(image_loc + "UAW_{}_anoms_era".format(settings.YEAR), cube, settings.COLOURMAP_DICT["circulation"], bounds, "Anomalies from 1981-2010 (m s"+r'$^{-1}$'+")")


#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
