#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for Terrestrial Water Storage (TWS) section.
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

import matplotlib.cm as mpl_cm
import matplotlib as mpl

import iris
import datetime as dt
import glob

import utils # RJHD utilities
import settings


data_loc = "/data/local/rdunn/SotC/{}/data/TWS/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

LEGEND_LOC = 'upper left'

START = dt.datetime(2002,8,1)
END = dt.datetime(2015,9,1)
MDI = -999.9

NOGRACE = [dt.datetime(2002,6,15), dt.datetime(2002,7,15), dt.datetime(2003,6,15), dt.datetime(2011,1,15),\
               dt.datetime(2011,6,15), dt.datetime(2012,5,15), dt.datetime(2012,10,15), dt.datetime(2013,3,15),\
               dt.datetime(2013,8,15), dt.datetime(2013,9,15), dt.datetime(2014,2,15), dt.datetime(2014,7,15),\
               dt.datetime(2014,12,15), dt.datetime(2015,6,15), dt.datetime(2015,10,15), dt.datetime(2015,11,15),\
               dt.datetime(2016,4,15), dt.datetime(2016,9,15), dt.datetime(2016,10,15)]

#************************************************************************
def decimal_to_dt(date):
    '''
    http://stackoverflow.com/questions/20911015/decimal-years-to-datetime-in-python

    '''
    year = int(date)
    rem = date - year

    base = dt.datetime(year, 1, 1)
    result = base + dt.timedelta(seconds=(base.replace(year=base.year + 1) - base).total_seconds() * rem)

    return result # decimal_to_dt

#************************************************************************
def test_if_missing(time):
    '''
    Test if a missing month is one of the listed ones
    '''

    result = False

    dt_time = decimal_to_dt(time)

    for ng in NOGRACE:
        # test each of the listed missing months
        if ng.year == dt_time.year and ng.month == dt_time.month:
            result = True
            print "testing {} {}".format(dt.datetime.strftime(dt_time, "%Y-%m-%d"), dt.datetime.strftime(ng, "%Y-%m-%d"))
           
            break

    return result # test_if_missing
 
#************************************************************************
def read_ts(filename):
    '''
    Read the GRACE timeseries and return a Timeseries object

    :param str filename: file to read

    :returns: Timeseries object
    '''

    indata = np.genfromtxt(filename, dtype = (float))

    data = np.ma.masked_where(indata[:,1] == 998.0, indata[:,1])

    return utils.Timeseries("GRACE", indata[:,0], data) # read_ts

#************************************************************************
def read_hovmuller_2015(data_loc):
    '''
    Retired for 2016 SotC but retained in case of value in future

    Scrappy subroutine to read in the files - one per month (some missing)
    Has latitudes and data as two columns. Appends to lists and then sorts out

    :param str data_loc: where to find the files

    :returns: times, latitudes and data
    '''

    # initialise the counters
    year = START.year
    month = START.month
    
    # initialise the holding lists
    data = []
    lats = []
    time = []

    # start at the beginning
    today = dt.datetime(year, month, 1)
    while today <= END:

        # if file exists, read in, else read in dummy data of same size
        try:
            indata = np.genfromtxt(data_loc + "avg_lat_csr05_ds_{:04g}{:02g}.txt".format(today.year, today.month), dtype = (float))

            data += [indata[:,1]]
            lats += [indata[:,0]]     
        
        except IOError:
            data += [[MDI for i in data[-1]]]
            lats += [lats[-1]]

        # convert to decimal years
        time += [year + (month-1)/12.]

        # sort out the counters if at the end of a year
        month += 1
        if month > 12:
            month = 1
            year += 1

        # and increment the while loop
        today = dt.datetime(year, month, 1)

    # conversion to array should ensure failure if change in coverage over time
    latitudes = np.array(lats)[0]
    # mask out the months where there weren't any data
    data = np.ma.masked_where(np.array(data) == MDI, np.array(data))

    return np.array(time), latitudes, data # read_hovmuller_2015

#************************************************************************
def read_hovmuller(data_loc):
    '''
    Scrappy subroutine to read in the files - one per month (some missing)
    Has latitudes and data as two columns. Appends to lists and then sorts out

    :param str data_loc: where to find the files

    :returns: times, latitudes and data
    '''

    # initialise the counters
    year = START.year
    month = START.month
    
    # initialise the holding lists
    data = []
    lats = []
    time = []

    # start at the beginning
    infiles = glob.glob(data_loc + "avg_lat_JPLM05_*.txt")
    infiles.sort()

    print "SORT MISSING MONTHS IN HOVMULLER - ACCOUNT FOR SATELLITE DRIFT"

    for filename in infiles:

        # convert to decimal years
        year = filename.split("/")[-1].split(".")[0].split("_")[-1]
        decimal = filename.split("/")[-1].split(".")[1]
        file_time = int(year) + float(decimal)/100.

        if file_time == 2012.0:
            print "file time adjusted as this is really December 2011"
            file_time = 2011.95

        dt_time = decimal_to_dt(file_time)        

        # test for missing months
        # as long as there's more than one time already added
        if len(time) > 1:
            # get the difference between current and previous in days
            ndays = (dt_time - decimal_to_dt(time[-1])).days

            print file_time, dt.datetime.strftime(dt_time, "%Y-%m-%d"), ndays

            while ndays > 45 : # if more than 45 days (expecting roughly every 30 or so)

                inserted_time = time[-1] + 0.1 # add increment of 36 days to test

                # only add blank data if the missing month is one of the listed ones
                if test_if_missing(inserted_time):

                    time += [inserted_time]
                    print "missing month added at {}".format(time[-1])

                    data += [[MDI for i in data[-1]]]
                    lats += [lats[-1]]
                    
                    ndays = (dt_time - decimal_to_dt(time[-1])).days

        time += [int(year) + float(decimal)/100.]

        # if file exists, read in, else read in dummy data of same size
        try:
            indata = np.genfromtxt(filename, dtype = (float))

            data += [indata[:,1]]
            lats += [indata[:,0]]     
        
        except IOError:
            data += [[MDI for i in data[-1]]]
            lats += [lats[-1]]
            
    # conversion to array should ensure failure if change in coverage over time
    latitudes = np.array(lats)[0]
    # mask out the months where there weren't any data
    data = np.ma.masked_where(np.array(data) == MDI, np.array(data))

    return np.array(time), latitudes, data # read_hovmuller

#************************************************************************
def read_map_data(filename):
    '''
    Read data for maps and convert to cube.  Given as single list files

    :param str filename: file to read

    :returns: cube
    '''
    
    grace = np.genfromtxt(filename, dtype = (float))

    lats = grace[:,1]
    lons = grace[:,0]
    anoms = grace[:,2]

    longitudes = np.unique(lons)
    latitudes = np.unique(lats)

    data = np.ma.zeros((len(latitudes), len(longitudes)))
    data.mask = np.ones(data.shape)
    
    for v, val in enumerate(anoms):

        xloc, = np.where(longitudes == lons[v])
        yloc, = np.where(latitudes == lats[v])

        data[yloc[0], xloc[0]] = val     

    cube = utils.make_iris_cube_2d(data, latitudes, longitudes, "TWS_anom", "mm")

    return cube


#************************************************************************
def run_all_plots():

    #************************************************************************
    # GRACE timeseries

    grace = read_ts(data_loc + "avg_JPLM05_land.txt")

    fig = plt.figure(figsize = (10,6))

    ax1 = plt.axes([0.1,0.1,0.88,0.88])

    utils.plot_ts_panel(ax1, [grace], "-", "hydrological", loc = LEGEND_LOC)

    #*******************
    # prettify
    ax1.set_ylim([-3.2,2.3])
    ax1.set_xlim([2002,int(settings.YEAR)+1.3])
    ax1.set_ylabel("Anomaly (cm)", fontsize = settings.FONTSIZE)

    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 

    plt.savefig(image_loc+"TWS_ts{}".format(settings.OUTFMT))
    plt.close()


    #************************************************************************
    # GRACE Hovmuller


    times, latitudes, data = read_hovmuller(data_loc)

    bounds = np.array([-200, -12, -9, -6, -3, 0, 3, 6, 9, 12, 200])

    utils.plot_hovmuller(image_loc + "TWS_hovmuller_grace", times, latitudes, data.T, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomaly (cm)")


    #************************************************************************
    # Difference Map

    cube = read_map_data(data_loc + "JPLM05_2016-2015.txt")

    utils.plot_smooth_map_iris(image_loc + "p2.1_TWS_{}_diffs".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Difference between 2016 and 2015 Equivalent Depth of Water (cm)", figtext = "(h) Terrestrial Water Storage")
    utils.plot_smooth_map_iris(image_loc + "TWS_{}_diffs".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Difference between 2016 and 2015 Equivalent Depth of Water (cm)")

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
