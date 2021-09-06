#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for lower stratosphere temperature (LST) section.
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
import copy
import subprocess
#import urllib2
import urllib.request, urllib.error, urllib.parse
import gc
import datetime as dt

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import matplotlib as mpl
import matplotlib.path as mpath
from matplotlib.ticker import MultipleLocator

import iris
import cartopy

import utils # RJHD utilities
import settings

DATALOC = "{}/{}/data/SLP/".format(settings.ROOTLOC, settings.YEAR)

LEGEND_LOC = 'lower left'
LW = 2

# data sources:
# AAO Obtained at PSD http://www.esrl.noaa.gov/psd/data/correlation/aao.data
# or http://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/aao.shtml
# AO http://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/ao.shtml
# NAO https://climatedataguide.ucar.edu/climate-data/hurrell-north-atlantic-oscillation-nao-index-station-based
# SOI BOM ftp://ftp.bom.gov.au/anon/home/ncc/www/sco/soi/soiplaintext.html <- manually fix header

# Winter NAO data calculated using IDL SLP_specialNAOtimeseries.pro script
# uses TIDL, give three-month range to cover last 3 winters.
# HadSLP from www.metoffice.gov.uk/hadobs/hadslp2

#*********************************************
def http(host, remote_loc, filename, local_loc, diagnostics=False):
    """
    Use HTTP handler via urllib2

    :param str host: HTTP hostename
    :param str remote_loc: remote directory
    :param str filename: remote filename
    :param str local_loc: local directory
    :param bool diagnostics: extra output
    """
    from socket import error as SocketError
    import errno

    success = 0
    try:
        if diagnostics:
            print("{}/{}/{}".format(host, remote_loc, filename))

        remote_file = urllib.request.urlopen("{}/{}/{}".format(host, remote_loc, filename))
        data = remote_file.read()

        local_file = open(os.path.join(local_loc, filename), "w")
        local_file.write(str(data, encoding = "latin-1"))
        local_file.close()
        if diagnostics:
            print("     Downloaded")
        success = 1

    # handle the error
    except urllib.error.URLError as e:
        print("Some sort of urllib2 error for {}".format("{}/{}/{}".format(host, remote_loc, filename)))
        if e.reason == "Not Found":
            # server doesn't have the file, so waiting won't help
            print("File {} Not Found - no further attempts".format(filename))
            success = 1

    except urllib.error.HTTPError as e:
        print(e.code)
        print(e.read())

    except SocketError as e:
        if e.errno != errno.ECONNRESET:
            raise # Not the error we are looking for
        
        # likely that too many calls at once.  
        # Sleep for a while and try again.  Handled in caller

    return success # http

#*********************************************
def doftp(host, remote_loc, filename, local_loc, user=None, passw=None, diagnostics=False):
    """
    Interface to Met Office doftp command

    :param str host: FTP hostename
    :param str remote_loc: remote directory
    :param str filename: remote filename
    :param str local_loc: local directory
    :param bool diagnostics: extra output
    """

    try:
        if user != None:
            if passw != None:
                # user and passw
                if diagnostics:
                    print("doftp -host {} -user {} -pass{} -cwd {} -get {}={}/{}".format(host, user, passw, remote_loc, filename, local_loc, filename))
                    subprocess.check_call(["doftp", '-host', host, '-user', user, '-pass', passw, '-cwd', remote_loc, '-get', "{}={}/{}".format(filename, local_loc, filename)])
            else:
                # user, no passw
                if diagnostics:
                    print("doftp -host {} -user {} -cwd {} -get {}={}/{}".format(host, user, remote_loc, filename, local_loc, filename))
                    subprocess.check_call(["doftp", '-host', host, '-user', user, '-cwd', remote_loc, '-get', "{}={}/{}".format(filename, local_loc, filename)])
        else:
            # no user or passw
            if diagnostics:
                print("doftp -host {} -cwd {} -get {}={}/{}".format(host, remote_loc, filename, local_loc, filename))
                subprocess.check_call(["doftp", '-host', host, '-cwd', remote_loc, '-get', "{}={}/{}".format(filename, local_loc, filename)])

        if diagnostics:
            print("     Downloaded")
                   
    # handle the error
    except subprocess.CalledProcessError:
        print("waiting 10 sec and trying again")
        import time
        time.sleep(10) # wait 10 seconds and onto next while loop

    except OSError:
        # executable not found
        print("Issue with doftp")
        raise OSError("doftp not found")

    return # doftp

#************************************************************************
def read_hadslp(filename):
    '''
    Read monthly HadSLP2 fields, returns cube

    '''


    years = []
    months = []

    # process each line, and separate each month using year/month note
    all_months = []
    with open(filename, 'r') as infile:
        for ll, line in enumerate(infile):

            if len(line.split()) == 2:
                years += [int(line.split()[0])]
                months += [int(line.split()[1])]

                if ll != 0:
                    all_months += [this_month]
                this_month = []

            else:
                this_month += [line.split()]
    # and get the final one!
    all_months += [this_month]
                
    all_months = np.array(all_months).astype(int) / 100. # convert to hPa
    all_months = np.swapaxes(all_months, 1, 2) # swap lats and lons around

    all_months = all_months[:, :, ::-1] # invert latitudes

    delta = 5.
    latitudes = np.arange(-90, 90 + delta, delta)
    longitudes = np.arange(-180 + delta/2., 180 + delta/2., delta)

    times = np.array([(dt.datetime(years[i], months[i], 1, 0, 0) - dt.datetime(years[0], months[0], 1, 0, 0)).days for i in range(len(years))])

    cube = utils.make_iris_cube_3d(all_months, times, "days since {}-{}-01 00:00:00".format(years[0], months[0]), longitudes, latitudes, "HadSLP", "bar")

    return cube # read_hadslp


#************************************************************************
def read_era5_slp():
    '''
    Read daily SLP from ERA5 downloaded locally, make monthly fields, and single cube.
    '''

    import os, fnmatch
    from iris.experimental.equalise_cubes import equalise_attributes
    from iris.util import unify_time_units

    era5_loc = "/data/users/hadjj/ERA5/"
    
    
    yearly_cubes = iris.cube.CubeList([])
    for year in range(1979, int(settings.YEAR)+1):
        print("processing ERA5 {}".format(year))
        monthly_cubes = iris.cube.CubeList([])
        for month in range(1, 13):

            this_month = iris.load(os.path.join(era5_loc, str(year), "{:02d}".format(month), "ERA5_MSLP_{}-{:02d}-*.nc".format(year, month)))

            equalise_attributes(this_month)
            unify_time_units(this_month)
            
            this_month = this_month.merge_cube() # these don't have a time axis, so merge
            
            monthly = this_month.collapsed(["time"], iris.analysis.MEAN)

            monthly = utils.regrid_cube(monthly, 1.0, 1.0) # regrid to 1x1 degree
    
            monthly_cubes += [monthly]

        # merge the 12 months into single year - try to split up this process
        print("making year cube {}".format(year))
        yearly_cubes += [monthly_cubes.merge_cube()] # these don't have a time axis so merge

    print("making single cube")
    unify_time_units(yearly_cubes)
    cube = yearly_cubes.concatenate_cube() #  these have a time axis, so concatenate

    # convert to hPa
    cube.data = cube.data/100.

    iris.save(cube, os.path.join(DATALOC, "era5_slp.nc"))
    print("finished ERA5 processing")

    return cube # read_era5_slp


#************************************************************************
def read_daily_aao(filename, name):
    '''
    Read the AAO data, returns Timeseries

    '''

    all_data = np.genfromtxt(filename, dtype=(float))

    years = all_data[:, 0].astype(int)
    months = all_data[:, 1].astype(int)
    days = all_data[:, 2].astype(int)

    times = np.array([dt.datetime(years[i], months[i], days[i]) for i in range(len(days))])

    data = all_data[:, 3]

    smoothed = utils.boxcar(data, 3)

    return utils.Timeseries(name, times, data), utils.Timeseries(name, times, smoothed)  # read_ao


#************************************************************************
def read_a_ao(filename, name, skip):
    '''
    Read the AO and AAO data, returns Timeseries

    '''

    all_data = np.genfromtxt(filename, dtype=(float), skip_header=skip)

    years = all_data[:, 0]
    months = all_data[:, 1]
    data = all_data[:, 2]

    times = np.array(years) + ((np.array(months)-1.)/12.)

    return utils.Timeseries(name, times, data) # read_a_ao

#************************************************************************
def read_soi(filename):
    '''
    Read the SOI data, returns Timeseries

    '''
    try:
        all_data = np.genfromtxt(filename, dtype=(float), skip_header=12, skip_footer=3)
    except ValueError:
        # presume last year has incomplete months
        all_data = np.genfromtxt(filename, dtype=(float), skip_header=12, skip_footer=4)
    
    years = all_data[:, 0]
    data = all_data[:, 1:]

    times = np.arange(years[0], years[-1]+1, 1/12.)

    return utils.Timeseries("SOI", times, data.reshape(-1)) # read_a_ao

#************************************************************************
def read_snao(filename):
    '''
    Read the SNAO data, returns Timeseries

    '''

    all_data = np.genfromtxt(filename, dtype=(float), skip_header=0)

    years = all_data[:, 0]
    data = all_data[:, 1]

    return utils.Timeseries("SNAO", years, data) # read_snao

#************************************************************************
def read_nao(filename):
    '''
    Read the NAO data, returns Timeseries

    '''

    all_data = np.genfromtxt(filename, dtype=(float), skip_header=9)

    years = all_data[:, 0]
    data = all_data[:, 1] # just get DJF column

    return utils.Timeseries("NAO", years, data) # read_nao

#************************************************************************
def plt_bars(ax, ts, text, w=1./12., invert=False, label=True):
    '''
    Plot the bars on an axes

    :param Axes ax: axes to plot on
    :param Timeseries ts: timeseries object
    :param str text: extra text for plotlabel
    :param float w: width of the bars
    :param bool invert: change order of colours

    '''

    if invert:
        plus = ts.data < 0.0
        minus = ts.data >= 0.0
    else:
        plus = ts.data > 0.0
        minus = ts.data <= 0.0

    ax.bar(ts.times[plus], ts.data[plus], color='r', ec="r", width=w, align="center")
    ax.bar(ts.times[minus], ts.data[minus], color='b', ec="b", width=w, align="center")
    ax.axhline(0, c='0.5', ls='--')
    if label:
        ax.text(0.03, 0.8, "{} {}".format(text, ts.name), transform=ax.transAxes, fontsize=settings.FONTSIZE)

    return # plt_months

#************************************************************************
def read_winter_nao(DATALOC, years):
    '''
    Read the NAO data, returns Timeseries and smoothed version

    '''

    ts_list = []
    smoothed = []
    
    for y in years:

        all_data = np.genfromtxt(DATALOC + "SLP_WinterNAOtimeseries_{}.txt".format(y), dtype=(float), skip_header=1)

        days = all_data[:, 1].astype(int)
        months = all_data[:, 0].astype(int)
        years = np.array([y for i in days]).astype(int)
        years[months < 12] = years[months < 12] + 1

        times = np.array([dt.datetime(years[i], months[i], days[i]) for i in range(len(days))])

        data = np.ma.masked_where(all_data[:, 2] <= -99.99, all_data[:, 2])
                                  
        ts_list += [utils.Timeseries("{}/{}".format(y, int(y[2:])+1), times, data)]

        smoothed += [np.ma.masked_where(all_data[:, 3] <= -99.99, all_data[:, 3])]

    return ts_list, smoothed # read_winter_nao

#************************************************************************
def run_all_plots(download=False):

    if download:
        # download if required

        # AAO
        http("https://www.esrl.noaa.gov/", "psd/data/correlation", "aao.data", DATALOC, diagnostics=True)
        http("https://www.cpc.ncep.noaa.gov", "products/precip/CWlink/daily_ao_index/aao", "monthly.aao.index.b79.current.ascii", DATALOC, diagnostics=True)
        doftp("ftp://ftp.cpc.ncep.noaa.gov", "/cwlinks/", "norm.daily.aao.index.b790101.current.ascii", DATALOC, diagnostics=True)
        # AO
        http("https://www.cpc.ncep.noaa.gov", "products/precip/CWlink/daily_ao_index", "monthly.ao.index.b50.current.ascii", DATALOC, diagnostics=True)
        # NAO
        http("https://climatedataguide.ucar.edu", "sites/default/files", "nao_station_seasonal.txt", DATALOC, diagnostics=True)
        print("HadSLP2r & SOI likely to fail")
        # HadSLP2r
        http("http://www.metoffice.gov.uk", "hadobs/hadslp2/data", "hadslp2r.asc.gz", DATALOC, diagnostics=True)
        # SOI 
        doftp("ftp://ftp.bom.gov.au", "/anon/home/ncc/www/sco/soi/", "soiplaintext.html", DATALOC, diagnostics=True)
        # https://metoffice.sharepoint.com/sites/NetworksCommsSite/SitePages/FTP-Proxy.aspx

        print("Now check downloads")
        sys.exit()


    #************************************************************************
    # Timeseries figures - different indices
    if True:
        SOI = read_soi(DATALOC + "soiplaintext.html")
        AO = read_a_ao(DATALOC + "monthly.ao.index.b50.current.ascii", "AO", 3)
        AAO = read_a_ao(DATALOC + "monthly.aao.index.b79.current.ascii", "AAO", 5)
    #    SNAO = read_snao(DATALOC + "JA STANDARDISED SNAO.txt") # from Chris Folland - unavailable in 2019
        NAO = read_nao(DATALOC + "nao_station_seasonal.txt")

        fig = plt.figure(figsize=(8, 6))

        # manually set up the 10 axes
        w = 0.42
        h = 0.23
        c = 0.51
        ax1 = plt.axes([c-w, 0.99-h, w, h])
        ax2 = plt.axes([c, 0.99-h, w, h])
        ax3 = plt.axes([c-w, 0.99-(2*h), w, h], sharex=ax1)
        ax4 = plt.axes([c, 0.99-(2*h), w, h], sharex=ax2)
        ax5 = plt.axes([c-w, 0.99-(3*h), w, h], sharex=ax1)
        ax6 = plt.axes([c, 0.99-(3*h), w, h], sharex=ax2)
        ax7 = plt.axes([c-w, 0.99-(4*h), w, h], sharex=ax1)
        ax8 = plt.axes([c, 0.99-(4*h), w, h], sharex=ax2)
    #    ax9 = plt.axes([c-w, 0.99-(5*h), w, h], sharex=ax1)
    #    ax10= plt.axes([c, 0.99-(5*h), w, h], sharex=ax2)

        plt_bars(ax1, SOI, "(a)", w=1./12, invert=True)
        plt_bars(ax3, AO, "(c)", w=1./12)
        plt_bars(ax5, AAO, "(e)", w=1./12)
        plt_bars(ax7, NAO, "(g)", w=1.)
    #    plt_bars(ax9, SNAO, "(i)", w=1.)

        ax1.set_xlim([1860, int(settings.YEAR)+1])

        plt_bars(ax2, SOI, "(b)", w=1./12, invert=True)
        plt_bars(ax4, AO, "(d)", w=1./12)
        plt_bars(ax6, AAO, "(f)", w=1./12)
        plt_bars(ax8, NAO, "(h)", w=0.9)
    #    plt_bars(ax10, SNAO, "(j)", w=0.9)

        ax2.set_xlim([2006, int(settings.YEAR)+2])
        for ax in [ax2, ax4, ax6, ax8]:#, ax10]:
            ax.yaxis.tick_right()

        # prettify
        for ax in [ax1, ax3, ax5, ax7]:#, ax9]:
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)
            utils.thicken_panel_border(ax)
            ax.yaxis.set_ticks_position('left')

        for ax in [ax2, ax4, ax6, ax8]:#, ax10]:
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)
            ax.yaxis.tick_right()
            for tick in ax.yaxis.get_major_ticks():
                tick.label2.set_fontsize(settings.FONTSIZE)
            utils.thicken_panel_border(ax)

        plt.setp([a.get_xticklabels() for a in fig.axes[:-2]], visible=False)

        ax1.xaxis.set_ticks([1880, 1920, 1960, 2000])
        ax2.xaxis.set_ticks([2007, 2010, 2013, 2016, 2019])
        ax1.yaxis.set_ticks([-40, 0, 40])
        ax2.yaxis.set_ticks([-40, 0, 40])
        ax3.yaxis.set_ticks([-4, 0, 4])
        ax4.yaxis.set_ticks([-4, 0, 4])
        ax5.yaxis.set_ticks([-2, 0, 2])
        ax6.yaxis.set_ticks([-2, 0, 2])
        ax7.yaxis.set_ticks([-4, 0, 4])
        ax8.yaxis.set_ticks([-4, 0, 4])
    #    ax9.yaxis.set_ticks([-2, 0, 2])
    #    ax10.yaxis.set_ticks([-2, 0, 2])

        minorLocator = MultipleLocator(1)
        ax2.xaxis.set_minor_locator(minorLocator)

        ax1.set_ylim([-40, 45])
        ax2.set_ylim([-40, 45])
        ax3.set_ylim([-4.5, 4.9])
        ax4.set_ylim([-4.5, 4.9])

        ax5.set_ylabel("Standard Units", fontsize=settings.FONTSIZE)

        plt.savefig(settings.IMAGELOC+"SLP_ts{}".format(settings.OUTFMT))

        plt.close()

    #************************************************************************
    # Timeseries figures - winter NAO
    #   initial run usually only Dec + Jan for recent year
    if True:
        # tried minor locator

        YEARS = ["2018", "2019", "2020"]
        plot_data, smoothed_data = read_winter_nao(DATALOC, YEARS)

        fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 6.5))

        TEXTS = ["(d)", "(e)", "(f)"]
        for a, ax in enumerate([ax1, ax2, ax3]):
            year = int(YEARS[a])

            plt_bars(ax, plot_data[a], TEXTS[a], w=1, invert=False)
            utils.thicken_panel_border(ax)
            ax.xaxis.set_ticks([dt.datetime(year, 12, 1), dt.datetime(year+1, 1, 1), dt.datetime(year+1, 2, 1), dt.datetime(year+1, 3, 1)])

            ax.plot(plot_data[a].times+dt.timedelta(hours=12), smoothed_data[a], "k-", lw=LW)
            ax.set_xlim([dt.datetime(year, 11, 27), dt.datetime(year+1, 2, 4) + dt.timedelta(days=30)])

        ax3.xaxis.set_ticklabels(["Dec", "Jan", "Feb", "Mar"], fontsize=settings.FONTSIZE)
        ax2.set_ylabel("Winter NAO Index (hPa)", fontsize=settings.FONTSIZE)

        fig.subplots_adjust(right=0.98, top=0.98, bottom=0.05, hspace=0.001)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

        for ax in [ax1, ax2, ax3]:
            ax.set_ylim([-69, 99])
#            ax.yaxis.set_ticks_position('left')

            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)

        plt.savefig(settings.IMAGELOC+"SLP_ts_winter_nao{}".format(settings.OUTFMT))
        plt.close()

    #************************************************************************
    # Timeseries figures - AAO
    #   initial run usually only Dec + Jan for recent year
    if True:
        # download latest file

        YEARS = ["2018", "2019", "2020"]
        plot_data, smoothed_data = read_daily_aao(DATALOC+"norm.daily.aao.index.b790101.current.ascii", "AAO")

        fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 6.5))

        TEXTS = ["(j)", "(k)", "(l)"]
        for a, ax in enumerate([ax1, ax2, ax3]):
            year = int(YEARS[a])

            locs, = np.where(np.logical_and(plot_data.times >= dt.datetime(year,12,1), plot_data.times < dt.datetime(year+1, 3, 1)))

            plt_bars(ax, utils.Timeseries("", plot_data.times[locs], plot_data.data[locs]), TEXTS[a], w=1, invert=False)
            ax.plot(plot_data.times[locs]+dt.timedelta(hours=12), smoothed_data.data[locs], "k-", lw=LW)

            utils.thicken_panel_border(ax)
            ax.xaxis.set_ticks([dt.datetime(year, 12, 1), dt.datetime(year+1, 1, 1), dt.datetime(year+1, 2, 1), dt.datetime(year+1, 3, 1)])

            ax.set_xlim([dt.datetime(year, 11, 27), dt.datetime(year+1, 2, 4) + dt.timedelta(days=30)])

        ax3.xaxis.set_ticklabels(["Dec", "Jan", "Feb", "Mar"], fontsize=settings.FONTSIZE)
        ax2.set_ylabel("AAO Index", fontsize=settings.FONTSIZE)

        fig.subplots_adjust(right=0.98, top=0.98, bottom=0.05, hspace=0.001)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

        for ax in [ax1, ax2, ax3]:
            ax.set_ylim([-5, 5])
#            ax.yaxis.set_ticks_position('left')

            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)

        plt.savefig(settings.IMAGELOC+"SLP_ts_winter_aao{}".format(settings.OUTFMT))
        plt.close()


    #************************************************************************
    # Global map
    if True:
        if not os.path.exists(os.path.join(DATALOC, "era5_slp.nc")):
            input("Is this running on SPICE (salloc --time=120 --mem=15G) as it will fail on VLD?")
            cube = read_era5_slp()
        else:
            cube = iris.load_cube(os.path.join(DATALOC, "era5_slp.nc"))
        # cube = read_hadslp(DATALOC + "hadslp2r.asc")

        # restrict to 1900 to last full year
        date_constraint = utils.periodConstraint(cube, dt.datetime(1900, 1, 1), dt.datetime(int(settings.YEAR)+1, 1, 1)) 
        cube = cube.extract(date_constraint)

        # convert to 1981-2010 climatology.
        clim_constraint = utils.periodConstraint(cube, dt.datetime(1981, 1, 1), dt.datetime(2011, 1, 1)) 
        clim_cube = cube.extract(clim_constraint)

        clim_data = clim_cube.data.reshape(-1, 12, clim_cube.data.shape[-2], clim_cube.data.shape[-1])

        # more than 15 years present
        climatology = np.ma.mean(clim_data, axis=0)
        nyears = np.ma.count(clim_data, axis=0)
        climatology = np.ma.masked_where(nyears <= 15, climatology) # Kate keeps GT 15.

        # extract final year
        final_year_constraint = utils.periodConstraint(cube, dt.datetime(int(settings.YEAR), 1, 1), dt.datetime(int(settings.YEAR)+1, 1, 1)) 
        final_year_cube = cube.extract(final_year_constraint)

        final_year_cube.data = final_year_cube.data - climatology

        # more than 6 months present
        annual_cube = final_year_cube.collapsed(['time'], iris.analysis.MEAN)
        nmonths = np.ma.count(final_year_cube.data, axis=0)
        annual_cube.data = np.ma.masked_where(nmonths <= 6, annual_cube.data)

        bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "SLP_{}_anoms_era5".format(settings.YEAR), annual_cube, settings.COLOURMAP_DICT["circulation"], bounds, "Anomalies from 1981-2010 (hPa)")
        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_SLP_{}_anoms_era5".format(settings.YEAR), annual_cube, settings.COLOURMAP_DICT["circulation"], bounds, "Anomalies from 1981-2010 (hPa)", figtext="(u) Sea Level Pressure")

        plt.close()

        del annual_cube
        del final_year_cube
        gc.collect()

    #************************************************************************
    # Polar Figures (1x3)
    if True:
        # apply climatology - incomplete end year, so repeat climatology and then truncate
        climatology = np.tile(climatology, ((cube.data.shape[0]//12)+1, 1, 1))
        anoms = iris.cube.Cube.copy(cube)

        anoms.data = anoms.data - climatology[0:anoms.data.shape[0], :, :]

        bounds = [-100, -8, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 8, 100]

        # north and south poles
        for hem in ("N", "S"):

            # set up a 1 x 3 set of axes
            fig = plt.figure(figsize=(4, 9.5))
            plt.clf()

            # set up plot settings
            cmap = settings.COLOURMAP_DICT["circulation"]
            norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)
            PLOTYEARS = [2018, 2019, 2020]
            if hem == "N":
                PLOTLABELS = ["(a) 2018/19", "(b) 2019/20", "(c) 2020/21"]
            elif hem == "S":
                PLOTLABELS = ["(g) 2018/19", "(h) 2019/20", "(i) 2020/21"]


            # boundary circle
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)


            # spin through axes
            for a in range(3):  

                if hem == "N":
                    ax = plt.subplot(3, 1, a+1, projection=cartopy.crs.NorthPolarStereo())
                    lat_constraint = utils.latConstraint([3, 80]) 
                elif hem == "S":
                    ax = plt.subplot(3, 1, a+1, projection=cartopy.crs.SouthPolarStereo())
                    lat_constraint = utils.latConstraint([-80, -3]) 

                plot_cube = iris.cube.Cube.copy(anoms)

                # extract 3 winter months
                date_constraint = utils.periodConstraint(anoms, dt.datetime(PLOTYEARS[a], 12, 1), dt.datetime(PLOTYEARS[a]+1, 3, 1)) 
                plot_cube = plot_cube.extract(date_constraint)

                # plot down to (almost) equator
                plot_cube = plot_cube.extract(lat_constraint)

                # take the mean
                try:
                    plot_cube = plot_cube.collapsed(['time'], iris.analysis.MEAN)
                except iris.exceptions.CoordinateCollapseError:
                    pass

                ax.gridlines() #draw_labels=True)
                ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
                ax.coastlines()
                ax.set_boundary(circle, transform=ax.transAxes)
                if hem == "N":
                    ax.set_extent([-180, 180, 0, 90], cartopy.crs.PlateCarree())
                elif hem == "S":
                    ax.set_extent([-180, 180, 0, -90], cartopy.crs.PlateCarree())

                ext = ax.get_extent() # save the original extent

                mesh = iris.plot.pcolormesh(plot_cube, cmap=cmap, norm=norm, axes=ax)

                ax.set_extent(ext, ax.projection) # fix the extent change from colormesh
                ax.text(0.0, 1.02, PLOTLABELS[a], fontsize=settings.FONTSIZE, transform=ax.transAxes)

            # add a colourbar for the figure
            cbar_ax = fig.add_axes([0.77, 0.07, 0.04, 0.9])
            cb = plt.colorbar(mesh, cax=cbar_ax, orientation='vertical', ticks=bounds[1:-1], drawedges=True)
            cb.ax.tick_params(axis='y', labelsize=settings.FONTSIZE, direction='in', size=0)
            cb.set_label(label="Anomaly (hPa)", fontsize=settings.FONTSIZE)

            # prettify
            cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
            cb.outline.set_linewidth(2)
            cb.dividers.set_color('k')
            cb.dividers.set_linewidth(2)


            fig.subplots_adjust(bottom=0.01, top=0.97, left=0.01, right=0.8, wspace=0.02)

            plt.title("")
            fig.text(0.03, 0.95, "", fontsize=settings.FONTSIZE * 0.8)

            plt.savefig(settings.IMAGELOC + "SLP_polar_era5_{}Hem{}".format(hem, settings.OUTFMT))
            plt.close()

            del plot_cube
        # ----
        del climatology
        gc.collect()

        #************************************************************************
        # SNAO figures
        if False:
            date_constraint = utils.periodConstraint(anoms, dt.datetime(int(settings.YEAR), 7, 1), dt.datetime(int(settings.YEAR), 9, 1)) 
            ja_cube = anoms.extract(date_constraint)
            ja_cube = ja_cube.collapsed(['time'], iris.analysis.MEAN)

            bounds = [-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100]
            cmap = settings.COLOURMAP_DICT["circulation"]
            norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)
            PLOTLABEL = "(a) July - August"

            fig = plt.figure(figsize=(8, 10))
            plt.clf()
            # make axes by hand
            axes = ([0.05, 0.45, 0.9, 0.55], [0.1, 0.05, 0.88, 0.35])

            ax = plt.axes(axes[0], projection=cartopy.crs.Robinson())

            ax.gridlines() #draw_labels=True)
            ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
            ax.coastlines()

            mesh = iris.plot.pcolormesh(ja_cube, cmap=cmap, norm=norm, axes=ax)

            ax.text(-0.05, 1.05, PLOTLABEL, fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)

            # colorbar
            cb = plt.colorbar(mesh, orientation='horizontal', ticks=bounds[1:-1], label="Anomaly (hPa)", drawedges=True, pad=0.05)

            # prettify
            cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
            cb.outline.set_linewidth(2)
            cb.dividers.set_color('k')
            cb.dividers.set_linewidth(2)

            # and the timeseries
            ax = plt.axes(axes[1])

            snao = np.genfromtxt(DATALOC+"{} DAILY SNAO.txt".format(settings.YEAR), dtype=(float))

            data = snao[:, 1]
            times = np.array([dt.datetime(int(settings.YEAR), 7, 1) + dt.timedelta(days=i) for i in range(len(data))])

            SNAO = utils.Timeseries("Summer NAO Index", times, data)

            plt_bars(ax, SNAO, "", w=0.7, label=False)
            ax.text(-0.1, 1.05, "(b) Summer NAO Index", fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)
            ax.yaxis.set_label("Summer NAO Index")
            # label every 14 days
            ticks = [dt.datetime(int(settings.YEAR), 7, 1) + dt.timedelta(days=i) for i in range(0, len(data), 14)]
            ax.xaxis.set_ticks(ticks)
            ax.xaxis.set_ticklabels([dt.datetime.strftime(t, "%d %b %Y") for t in ticks])
            minorLocator = MultipleLocator(1)
            ax.xaxis.set_minor_locator(minorLocator)

            utils.thicken_panel_border(ax3)

            plt.savefig(settings.IMAGELOC + "SLP_SNAO{}".format(settings.OUTFMT))
            plt.close()

            del ja_cube
            gc.collect()

        #************************************************************************
        # EU figure
        if False:
            contour = True

            date_constraint = utils.periodConstraint(anoms, dt.datetime(int(settings.YEAR), 7, 1), dt.datetime(int(settings.YEAR), 8, 1)) 
            jul_cube = anoms.extract(date_constraint)
            date_constraint = utils.periodConstraint(anoms, dt.datetime(int(settings.YEAR), 8, 1), dt.datetime(int(settings.YEAR), 9, 1)) 
            aug_cube = anoms.extract(date_constraint)

            bounds = [-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100]
            cmap = settings.COLOURMAP_DICT["circulation"]
            norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)
            PLOTLABELS = ["(a) July", "(b) August"]

            fig = plt.figure(figsize=(8, 14))
            plt.clf()
            # make axes by hand
            axes = ([0.1, 0.65, 0.8, 0.3], [0.1, 0.35, 0.8, 0.3], [0.1, 0.05, 0.8, 0.25])

            cube_list = [jul_cube, aug_cube]
            for a in range(2):
                ax = plt.axes(axes[a], projection=cartopy.crs.PlateCarree(central_longitude=-10))

                ax.gridlines() #draw_labels=True)
                ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
                ax.coastlines()
                ax.set_extent([-70, 40, 30, 80], cartopy.crs.PlateCarree())

                if contour:
                    from scipy.ndimage.filters import gaussian_filter
                    sigma = 0.5
                    cube = cube_list[a]

                    cube.data = gaussian_filter(cube.data, sigma)
                    mesh = iris.plot.contourf(cube, bounds, cmap=cmap, norm=norm, axes=ax)

                else:
                    mesh = iris.plot.pcolormesh(cube_list[a], cmap=cmap, norm=norm, axes=ax)

                ax.text(-0.05, 1.05, PLOTLABELS[a], fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)

            # colorbar
            cb = plt.colorbar(mesh, orientation='horizontal', ticks=bounds[1:-1], label="Anomaly (hPa)", drawedges=True)

            # prettify
            cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
            cb.outline.set_linewidth(2)
            cb.dividers.set_color('k')
            cb.dividers.set_linewidth(2)

            # and the timeseries
            ax3 = plt.axes(axes[2])

            snao = np.genfromtxt(DATALOC+"{} DAILY SNAO.txt".format(settings.YEAR), dtype=(float))

            data = snao[:, 1]
            times = np.array([dt.datetime(int(settings.YEAR), 7, 1) + dt.timedelta(days=i) for i in range(len(data))])

            SNAO = utils.Timeseries("Summer NAO Index", times, data)

            plt_bars(ax3, SNAO, "", w=0.7, label=False)
            ax3.text(-0.05, 1.05, "(c) Summer NAO Index", fontsize=settings.FONTSIZE * 0.8, transform=ax3.transAxes)
            ax3.yaxis.set_label("Summer NAO Index")
            # label every 14 days
            ticks = [dt.datetime(int(settings.YEAR), 7, 1) + dt.timedelta(days=i) for i in range(0, len(data), 14)]
            ax3.xaxis.set_ticks(ticks)
            ax3.xaxis.set_ticklabels([dt.datetime.strftime(t, "%d %b %Y") for t in ticks])
            minorLocator = MultipleLocator(1)
            ax3.xaxis.set_minor_locator(minorLocator)
            ax3.yaxis.set_ticks_position('left')

            utils.thicken_panel_border(ax3)

            plt.savefig(settings.IMAGELOC + "SLP_NAtlantic{}".format(settings.OUTFMT))
            plt.close()

            del cube
            del anoms
            gc.collect()


    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--download', dest='download', action='store_true', default=False,
                        help='Force download of data')

    args = parser.parse_args()         
    run_all_plots(download=args.download)

#************************************************************************
#                                 END
#************************************************************************
