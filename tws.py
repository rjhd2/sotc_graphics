#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for Terrestrial Water Storage (TWS) section.
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

import datetime as dt
import glob
import numpy as np
import matplotlib.pyplot as plt

import utils # RJHD utilities
import settings


data_loc = "{}/{}/data/TWS/".format(settings.ROOTLOC, settings.YEAR)
reanalysis_loc = "{}/{}/data/RNL/".format(settings.ROOTLOC, settings.YEAR)
image_loc = "{}/{}/images/".format(settings.ROOTLOC, settings.YEAR)

LEGEND_LOC = 'lower left'

START = dt.datetime(2002, 8, 1)
END = dt.datetime(2017, 12, 1)
MDI = -999.9

NOGRACE = [dt.datetime(2002, 6, 15), dt.datetime(2002, 7, 15), dt.datetime(2003, 6, 15), dt.datetime(2011, 1, 15), \
               dt.datetime(2011, 6, 15), dt.datetime(2012, 5, 15), dt.datetime(2012, 10, 15), dt.datetime(2013, 3, 15), \
               dt.datetime(2013, 8, 15), dt.datetime(2013, 9, 15), dt.datetime(2014, 2, 15), dt.datetime(2014, 7, 15), \
               dt.datetime(2014, 12, 15), dt.datetime(2015, 6, 15), dt.datetime(2015, 10, 15), dt.datetime(2015, 11, 15), \
               dt.datetime(2016, 4, 15), dt.datetime(2016, 9, 15), dt.datetime(2016, 10, 15), dt.datetime(2016, 11, 3), \
               dt.datetime(2017, 2, 28)]

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
            print("testing {} {}".format(dt.datetime.strftime(dt_time, "%Y-%m-%d"), dt.datetime.strftime(ng, "%Y-%m-%d")))
           
            break

    return result # test_if_missing
 
#************************************************************************
def read_ts(filename, name):
    '''
    Read the GRACE timeseries and return a Timeseries object

    :param str filename: file to read
    :param str name: name of timeseries

    :returns: Timeseries object
    '''

    indata = np.genfromtxt(filename, dtype=(float))

    if name == "GRACE":
        data = np.ma.masked_where(indata[:, 1] == 998.0, indata[:, 1])
    times = indata[:, 0].astype(str)

    if name == "Model":
        years = np.array([str(y)[:4] for y in times]).astype(float)
        months = np.array([str(y)[4:6] for y in times]).astype(float)

        time = years + (months-1)/12. + (1/24.) # place in centre of month

        # mask first 3 months (Matt Rodell email, 14-2-2019)
        data = np.ma.array(indata[:, 2])
        data.mask = np.zeros(data.shape[0])
        data.mask[:3] = True
    else:
        time = times.astype(float)        

    return utils.Timeseries(name, time, data) # read_ts

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
            indata = np.genfromtxt(data_loc + "avg_lat_csr05_ds_{:04g}{:02g}.txt".format(today.year, today.month), dtype=(float))

            data += [indata[:, 1]]
            lats += [indata[:, 0]]     
        
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
def read_hovmuller_2017(data_loc):
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

    print("SORT MISSING MONTHS IN HOVMULLER - ACCOUNT FOR SATELLITE DRIFT")

    for filename in infiles:
        print(filename)

        # convert to decimal years
        year = filename.split("/")[-1].split(".")[0].split("_")[-1]
        decimal = filename.split("/")[-1].split(".")[1]
        file_time = int(year) + float(decimal)/100.

        if file_time == 2012.0:
            print("file time adjusted as this is really December 2011")
            file_time = 2011.95

        dt_time = decimal_to_dt(file_time)        

        # test for missing months
        # as long as there's more than one time already added
        if len(time) > 1:
            # get the difference between current and previous in days
            ndays = (dt_time - decimal_to_dt(time[-1])).days

            print(file_time, dt.datetime.strftime(dt_time, "%Y-%m-%d"), ndays)

            while ndays > 45: # if more than 45 days (expecting roughly every 30 or so)

                inserted_time = time[-1] + 0.1 # add increment of 36 days to test

                # only add blank data if the missing month is one of the listed ones
                if test_if_missing(inserted_time):

                    time += [inserted_time]
                    print("missing month added at {}".format(time[-1]))

                    data += [[MDI for i in data[-1]]]
                    lats += [lats[-1]]
                    
                    ndays = (dt_time - decimal_to_dt(time[-1])).days

        time += [int(year) + float(decimal)/100.]

        # if file exists, read in, else read in dummy data of same size
        try:
            indata = np.genfromtxt(filename, dtype=(float))

            data += [indata[:, 1]]
            lats += [indata[:, 0]]     
        
        except IOError:
            data += [[MDI for i in data[-1]]]
            lats += [lats[-1]]
            
    # conversion to array should ensure failure if change in coverage over time
    latitudes = np.array(lats)[0]
    # mask out the months where there weren't any data
    data = np.ma.masked_where(np.array(data) == MDI, np.array(data))

    return np.array(time), latitudes, data # read_hovmuller_2017

#************************************************************************
def read_hovmuller(filename):
    '''
    Scrappy subroutine to read in file
    Has latitudes, dates and data as three columns. 

    :param str filename: file to read

    :returns: times, latitudes and data
    '''

    indata = np.genfromtxt(filename, dtype=(float))

    latitudes = np.unique(indata[:, 0])
    latitudes.sort()
    times = np.unique(indata[:, 1]).astype(int)
    times.sort()

    data = np.zeros([len(latitudes), len(times)])

    for d in indata:

        lat, = np.where(d[0] == latitudes)
        tim, = np.where(d[1] == times)

        data[lat, tim] = d[2]

    # mask out the months where there weren't any data
    data = np.ma.masked_where(np.array(data) == MDI, np.array(data))

    years = np.array([str(y)[:4] for y in times]).astype(float)
    months = np.array([str(y)[4:6] for y in times]).astype(float)

    time = years + (months-1)/12.
    
    return np.array(time), latitudes, data # read_hovmuller

#************************************************************************
def read_map_data(filename):
    '''
    Read data for maps and convert to cube.  Given as single list files

    :param str filename: file to read

    :returns: cube
    '''
    
    grace = np.genfromtxt(filename, dtype=(float))

    lats = grace[:, 0]
    lons = grace[:, 1]
    anoms = grace[:, 2]

    longitudes = np.unique(lons)
    latitudes = np.unique(lats)

    data = np.ma.zeros((len(latitudes), len(longitudes)))
    data.mask = np.ones(data.shape)
    
    for v, val in enumerate(anoms):

        xloc, = np.where(longitudes == lons[v])
        yloc, = np.where(latitudes == lats[v])

        data[yloc[0], xloc[0]] = val     

    data = np.ma.masked_where(data == -9999.0, data)

    cube = utils.make_iris_cube_2d(data, latitudes, longitudes, "TWS_anom", "mm")

    return cube


#************************************************************************
def run_all_plots():

    #************************************************************************
    # GRACE timeseries

    model = read_ts(data_loc + "glb_avg_tws_2003-2018_model.txt", "Model")
    model.ls = "--"
    grace = read_ts(data_loc + "glb_avg_tws_2002-2017_GRACE.txt", "GRACE")

    fig = plt.figure(figsize=(10, 6))

    ax1 = plt.axes([0.1, 0.1, 0.88, 0.88])

    utils.plot_ts_panel(ax1, [grace, model], "-", "hydrological", loc=LEGEND_LOC)

    #*******************
    # prettify
    ax1.set_ylim([-5.2, 2.3])
    ax1.set_xlim([2003, int(settings.YEAR)+1.3])
    ax1.set_ylabel("Anomaly (cm)", fontsize=settings.FONTSIZE)

    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 

    plt.savefig(image_loc+"TWS_ts{}".format(settings.OUTFMT))
    plt.close()


    #************************************************************************
    # GRACE Hovmuller


    times, latitudes, data = read_hovmuller(data_loc + "zonal_mean_tws.txt")

    bounds = np.array([-200, -12, -9, -6, -3, 0, 3, 6, 9, 12, 200])

    utils.plot_hovmuller(image_loc + "TWS_hovmuller_grace", times, latitudes, data, settings.COLOURMAP_DICT["hydrological"], bounds, "Anomaly (cm)")


    #************************************************************************
    # Difference Map

    cube = read_map_data(data_loc + "tws_changes_{}-{}.txt".format(settings.YEAR, int(settings.YEAR)-1))

    utils.plot_smooth_map_iris(image_loc + "p2.1_TWS_{}_diffs".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Difference between {} and {} Equivalent Depth of Water (cm)".format(settings.YEAR, int(settings.YEAR)-1), figtext="(q) Terrestrial Water Storage")

    utils.plot_smooth_map_iris(image_loc + "TWS_{}_diffs".format(settings.YEAR), cube, settings.COLOURMAP_DICT["hydrological"], bounds, "Difference between {} and {} Equivalent Depth of Water (cm)".format(settings.YEAR, int(settings.YEAR)-1))

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
