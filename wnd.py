#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for winds (WND) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 23                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2018-06-05 17:55:11 +0100 (Tue, 05 Jun #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.cm as mpl_cm
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator

import iris
import iris.quickplot as qplt
import cartopy.crs as ccrs

import netCDF4 as ncdf

import utils # RJHD utilities
import settings


data_loc = "/data/local/rdunn/SotC/{}/data/WND/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

CLIM_PERIOD="1981-2010"

LEGEND_LOC = 'upper left'
LW = 3
BBOX = (0,0.9)
YLIM = [-0.28, 0.45]


#************************************************************************
def combine_arrays(tup):
    """
    Manually combine N masked arrays and masks

    :param tuple tup: input tuple containing arrays
    :returns: array - final combined array
    """


    return np.ma.masked_array(np.hstack(tup), \
                                  mask=np.hstack([arr.mask for arr in tup]))
    

#************************************************************************
def read_australia(filename, mean_then_clim = True):
    """
    Read in the Australian data from Tim McVicar

    Resaved .xls as .csv and removed some rows and columns

    :param str filename: full path and name of input file
    :param bool mean_then_clim: how to calculate the anomalies

    :returns: lat, lon, anomalies and trends

    """

    print "save .xls as .csv from LibreOffice using fixed space format.\n Ignore station name column\n change trend calculation dates if necessary"

    indata = np.genfromtxt(filename, dtype=(str), skip_header = 3, skip_footer = 3)

    aus_lat = np.array([float(x) for x in indata[:,3]])
    aus_lon = np.array([float(x) for x in indata[:,2]])
    aus_anom_8110 = np.array([float(x) for x in indata[:,-3]])
    aus_clim_8110 = np.array([float(x) for x in indata[:,-4]])
    aus_anom_8810 = np.array([float(x) for x in indata[:,-1]])
    aus_clim_8810 = np.array([float(x) for x in indata[:,-2]])
    aus_trend_79_pres = np.array([float(x) for x in indata[:,-6]])

    all_timeseries = np.array(indata[:, 5:-6], dtype=np.float)


    years = np.arange(1974, int(settings.YEAR) + 1, 1)
    clim_start, = np.where(years == 1981)
    clim_end, = np.where(years == 2010+1)

    # two methods of calculating regional anomalies
    if mean_then_clim:
        # take mean of all station series to form regional actuals
        # then subtract climatology taken from regional mean series
        mean_timeseries = np.mean(all_timeseries, axis=0)
        anomalies = mean_timeseries - (np.mean(mean_timeseries[clim_start[0] : clim_end[0]]))

    else:
        # calculate climatology for each station, get anomalies for each station
        # then calculate the mean
        clims = np.mean(all_timeseries[:, clim_start : clim_end], axis = 1)
        all_anomalies = all_timeseries - np.tile(clims, (all_timeseries.shape[1], 1)).transpose()
        anomalies = np.mean(all_anomalies, axis = 0)

    return  aus_lat, aus_lon, np.ma.array(aus_anom_8110, mask = np.zeros(len(aus_anom_8110))), \
        np.ma.array(aus_trend_79_pres, mask = np.zeros(len(aus_trend_79_pres))), \
        years, anomalies, aus_clim_8110 # read_australia

#************************************************************************
def read_hadisd_annual_anomalies(region):
    '''
    Read the annual anomalies (also at 3 and 10m) from HadISD output
    Also read the slope, number of stations and mean speed

    :param obj region: region object
    :returns: arrays of years, anomalies, 3m and 10m anomalies
    '''


    indata = np.genfromtxt("{}/{}_Wind_annual_anomalies_{}.dat".format(data_loc, region.fname, CLIM_PERIOD), dtype=(str), skip_header = 13)
    years = np.array([float(x) for x in indata[:,0]])
    anomaly = np.array([float(x) for x in indata[:,1]])
    meter_3 = np.array([float(x) for x in indata[:,2]])
    meter_10 = np.array([float(x) for x in indata[:,3]])

    # also extract data for table

    with file("{}/{}_Wind_annual_anomalies_{}.dat".format(data_loc, region.fname, CLIM_PERIOD), "r") as infile:

        for line in infile:
            if len(line.strip()) == 0:
                continue

            # if data there, convert to float and store
            if line.split()[0] == "slope":
                region.slope = float(line.split()[2])
            if line.split()[0] == "Nstations":
                region.nstat = float(line.split()[1])
            if line.split()[0] == "Mean":
                region.mean = float(line.split()[3])

    return years, anomaly, meter_3, meter_10 # read_hadisd_annual_anomalies



#************************************************************************
def read_hadisd_global_summary(australia = False):
    '''
    Read the station by station data.  Trends increased to per decade

    :param bool australia: include australia or not
    :returns: ids, lats, lons, means, anomalies and trends
    '''


    if australia:
        indata = np.genfromtxt("{}/Global_sotc_summary_clim_{}_trend_1979.csv".format(data_loc, CLIM_PERIOD), dtype=(str), skip_header = 10, skip_footer = 6)
    else:
        indata = np.genfromtxt("{}/GlobalNoOz_sotc_summary_clim_{}_trend_1979.csv".format(data_loc, CLIM_PERIOD), dtype=(str), skip_header = 10, skip_footer = 6)

    stn_id = indata[:,0]
    lon = np.array([float(x) for x in indata[:,1]])
    lat = np.array([float(x) for x in indata[:,2]])
    low = np.array([float(x) for x in indata[:,3]])

    mean8110 = np.array([float(x) for x in indata[:,4]])
    anomaly8110 = np.array([float(x) for x in indata[:,5]])
    trend79_pres = np.array([float(x) for x in indata[:,6]]) 

    mean8110 = np.ma.masked_where(mean8110 == -999.0, mean8110)
    anomaly8110 = np.ma.masked_where(anomaly8110 == -999.0, anomaly8110)
    trend79_pres = np.ma.masked_where(trend79_pres == -999.0, trend79_pres)
    
    if len(mean8110.mask.shape) == 0: mean8110.mask = np.zeros(mean8110.shape[0])
    if len(anomaly8110.mask.shape) == 0: anomaly8110.mask = np.zeros(anomaly8110.shape[0])
    if len(trend79_pres.mask.shape) == 0: trend79_pres.mask = np.zeros(trend79_pres.shape[0])


    return stn_id, lon, lat,  mean8110, anomaly8110, trend79_pres # read_hadisd_global_summary


#************************************************************************
def read_map_data(filename):
    '''
    Read data for maps and convert to cube.  Given as single list files

    :param str filename: file to read

    :returns: cube
    '''
    
    era = np.genfromtxt(filename, dtype = (float), skip_header = 3)

    lats = era[:,0]
    lons = era[:,1]
    anoms = era[:,2]

    longitudes = np.unique(lons)
    latitudes = np.unique(lats)

    data = np.ma.zeros((len(latitudes), len(longitudes)))
    data.mask = np.ones(data.shape)
    
    for v, val in enumerate(anoms):

        xloc, = np.where(longitudes == lons[v])
        yloc, = np.where(latitudes == lats[v])

        data[yloc[0], xloc[0]] = val     

    cube = utils.make_iris_cube_2d(data, latitudes, longitudes, "WND_anom", "m/s")

    return cube # read_map_data

#************************************************************************
def read_radiometer(filename):
    '''
    Read radiometer data - may include ERA - check columns

    '''
    
    mdi = -99.9999
    indata = np.genfromtxt(filename, dtype=(float), skip_header = 2, missing_values = "-NaN", filling_values = mdi)

    ts = utils.Timeseries("SSM/I+SSMIS",indata[:,0], np.ma.masked_where(indata[:,3] <= mdi, indata[:,3]))
   
    return ts # read_radiometer

#************************************************************************
def read_ts_cube(filename, variable, name):

    cube = iris.load(filename, variable)[0]

    monthly = cube.data
    monthly = monthly.reshape(-1, 12)
    annual = np.ma.mean(monthly, axis = 1)

    start = int(cube[0].coord("time").units.name.split(" ")[-1])
    times = np.array([start + i for i in range(len(annual))])

    return utils.Timeseries(name, times, annual) # read_ts_cube

#************************************************************************
class Region(object):
    '''
    Class for region
    '''
    
    def __init__(self, name, fname, color):
        self.name = name
        self.color = color
        self.fname = fname

    def __str__(self):     
        return "Region: {}".format(self.name)

    __repr__ = __str__
   

#************************************************************************
#************************************************************************


Australia = Region("Australia", "Australia", "green")
Europe = Region("Europe", "Europe", "blue")
EastAsia = Region("East Asia", "EastAsia", "brown")
CentAsia = Region("Central Asia", "CentAsia", "red")
NorthUSA = Region("North America", "NorthAmer", "grey")
Globe = Region("Globe (excl Austr)", "GlobalNoOz", "black")

all_regions = [Globe, NorthUSA, Europe, CentAsia, EastAsia, Australia]

#************************************************************************
def run_all_plots():

    #************************************************************************
    # Read in Australian Data
    aus_lat, aus_lon, aus_anom_8110, aus_trend_79_pres, aus_years, aus_anomalies, aus_clim_8110 = read_australia("{}/u_Extracted_ANN_1974_{}_RD.csv".format(data_loc, settings.YEAR))

    # Timeseries figure

    # ERA data
    era_globe, era_oceans, era_land, eratropics = utils.era_ts_read(reanalysis_loc, "wnd", annual = True)

    land_era_clim, land_era_anoms = utils.calculate_climatology_and_anomalies_1d(era_land, 1981, 2010)

    land_merra_anoms = utils.read_merra(reanalysis_loc + "MERRA-2_SfcAnom_ts_{}.dat".format(settings.YEAR), "wind", "L", anomalies = True)

    jra_actuals, jra_anoms = utils.read_jra55(reanalysis_loc + "JRA-55_ws10m_globalland_ts.txt", "windspeed")

    # Plot timeseries figure
    fig, (ax1, ax2, ax3) = plt.subplots(3, figsize = (10,12), sharex=True)

    print "{} {} {} {} {}".format("name","mean","anomaly","trend/dec","N station")
    for region in all_regions:

        if region.name == "Australia":
            # use the Australian data
 
            years = aus_years
            anomalies = aus_anomalies

            print anomalies

            region.nstat = len(aus_lat)
            region.slope = np.mean(aus_trend_79_pres) # just do the mean of all station slopes
            region.mean = np.mean(aus_clim_8110)

        else:
            # Read in the HadISD annual anomalies
            years, anomalies, m3, m10 = read_hadisd_annual_anomalies(region)

        # Print data for table
        print "{} {} {} {} {}".format(region.name, region.mean, anomalies[-1], region.slope * 10., region.nstat)

        order = anomalies.argsort()
        ranks = order.argsort()

        print "{} highest {} lowest".format(len(ranks) - ranks[-1], ranks[-1])

        # plot data
        if region.fname == "GlobalNoOz":
            ax1.plot(years, anomalies, c = region.color, label = region.name, lw = 3)
        else:
            ax1.plot(years, anomalies, c = region.color, label = region.name, lw = 2)

        if region.fname == "GlobalNoOz":
            ax2.plot(years, m3, c = region.color, lw = 3)
            ax3.plot(years, m10, c = region.color, lw = 3)
        elif region.name != "Australia":
            ax2.plot(years, m3, c = region.color, lw = 2)
            ax3.plot(years, m10, c = region.color, lw = 2)

        if region.fname == "GlobalNoOz":
            ax1.plot(land_era_anoms.times, land_era_anoms.data, c = settings.COLOURS["circulation"]["ERA-Interim"], label = "ERA-Interim (land only)", lw = 2)
            ax1.plot(land_merra_anoms.times, land_merra_anoms.data, c = settings.COLOURS["circulation"]["MERRA-2"], label = "MERRA-2 (land only)", lw = 2)
            ax1.plot(jra_anoms.times, jra_anoms.data, c = settings.COLOURS["circulation"]["JRA-55"], label = "JRA-55 (land only)", lw = 2)

    # finish off plot
    ax1.axhline(0, c = '0.5', ls = '--')

    ax1.text(0.02, 0.9, "(a) In Situ all Speeds", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE)
    ax2.text(0.02, 0.9, "(b) In Situ >3 m s"+r'$^{-1}$'+" Winds", transform = ax2.transAxes, fontsize = settings.LABEL_FONTSIZE)
    ax3.text(0.02, 0.9, "(c) In Situ >10 m s"+r'$^{-1}$'+" Winds", transform = ax3.transAxes, fontsize = settings.LABEL_FONTSIZE)

    fig.text(0.03, 0.82, "Wind Anomaly (m s"+r'$^{-1}$'+")", va='center', rotation='vertical', fontsize = settings.LABEL_FONTSIZE)
    fig.text(0.03, 0.4, "Wind Frequency (% yr"+r'$^{-1}$'+")", va='center', rotation='vertical', fontsize = settings.LABEL_FONTSIZE)

    ax1.legend(loc="upper right", ncol=2, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, labelspacing=0.1, columnspacing=0.5, bbox_to_anchor=(1.0,0.95))

    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

    fig.subplots_adjust(right = 0.95, top = 0.95, hspace = 0.001)

    for tick in ax3.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 

    minorLocator = MultipleLocator(1)
    for ax in [ax1, ax2, ax3]:
        utils.thicken_panel_border(ax)
        ax.set_yticks(ax.get_yticks()[1:])
        ax.xaxis.set_minor_locator(minorLocator)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
        ax.yaxis.set_ticks_position('left')

    plt.xlim([1970,int(settings.YEAR)+1])

    plt.savefig(image_loc+"WND_land_ts{}".format(settings.OUTFMT))
    plt.close()

    #************************************************************************
    # 2015 HadISD Anomaly figure

    # Read in HadISD station anomalies
    stn_id, hadisd_lons, hadisd_lats, mean8110, hadisd_anomaly_8110, hadisd_trend_79_pres = read_hadisd_global_summary()

    # combine together 
    #lats = hadisd_lats
    #lons = hadisd_lons
    #anom = hadisd_anomaly_8110
    #trend = hadisd_trend_79_pres

    lats = np.append(hadisd_lats, aus_lat)
    lons = np.append(hadisd_lons, aus_lon)
    anom = combine_arrays((hadisd_anomaly_8110, aus_anom_8110))
    trend = combine_arrays((hadisd_trend_79_pres, aus_trend_79_pres)) * 10.

    bounds=[-100, -0.8, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4, 0.8, 100]
    bounds=[-100, -0.4, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.4, 100]
    utils.scatter_plot_map(image_loc + "WND_{}_obs_trend".format(settings.YEAR), trend, lons, lats, settings.COLOURMAP_DICT["circulation"], bounds, "(m s"+r'$^{-1}$'+" decade"+r'$^{-1}$)')

    bounds=[-100, -1.2, -0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.8, 1.2, 100]
    utils.scatter_plot_map(image_loc + "WND_{}_obs_anomaly".format(settings.YEAR), anom, lons, lats, settings.COLOURMAP_DICT["circulation"], bounds, "(m s"+r'$^{-1}$)')

    total = float(len(anom.compressed()))
    pos, = np.ma.where(anom > 0)
    neg, = np.ma.where(anom < 0)
    print "positive {} negative {}".format(len(pos)/total, len(neg)/total)
    pos, = np.ma.where(anom > 0.5)
    neg, = np.ma.where(anom < -0.5)
    print "positive {} negative {} (than 0.5)".format(len(pos)/total, len(neg)/total)
    pos, = np.ma.where(anom > 1.0)
    neg, = np.ma.where(anom < -1.0)
    print "positive {} negative {} (than 1.0)".format(len(pos)/total, len(neg)/total)


    #************************************************************************
    # ERA + HadISD Anomaly figure

    # Read in ERA anomalies

    #cube_list = iris.load(reanalysis_loc + "SI10_ansf_moda_ann{}{}-ann19812010.nc".format(settings.YEAR, settings.YEAR))

    #cube = cube_list[0]
    #cube.coord('latitude').guess_bounds()
    #cube.coord('longitude').guess_bounds()

    #bounds=[-4, -1.2, -0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.8, 1.2, 4]

    #utils.plot_smooth_map_iris(image_loc + "WND_{}_era_obs_anomaly".format(settings.YEAR), cube[0][0], settings.COLOURMAP_DICT["circulation_r"], bounds, "Anomalies from 1981-2010 (m s"+r'$^{-1}$)', scatter = (lons, lats, anom))
    #utils.plot_smooth_map_iris(image_loc + "WND_{}_era_anomaly".format(settings.YEAR), cube[0][0], settings.COLOURMAP_DICT["circulation_r"], bounds, "Anomalies from 1981-2010 (m s"+r'$^{-1}$)')

    #************************************************************************
    # MERRA + HadISD Anomaly figure

    # Read in MERRA anomalies

    cube_list = iris.load(reanalysis_loc + "MERRA-2_SfcAnom_{}.nc".format(settings.YEAR), "10m Wind Speed Anomaly (1981-2010)")

    cube = cube_list[0]
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()

    bounds=[-4, -1.2, -0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.8, 1.2, 4]

    utils.plot_smooth_map_iris(image_loc + "WND_{}_merra_obs_anomaly".format(settings.YEAR), cube[0], settings.COLOURMAP_DICT["circulation_r"], bounds, "Anomalies from 1981-2010 (m s"+r'$^{-1}$)', scatter = (lons, lats, anom))
    utils.plot_smooth_map_iris(image_loc + "p2.1_WND_{}_merra_obs_anomaly".format(settings.YEAR), cube[0], settings.COLOURMAP_DICT["circulation_r"], bounds, "Anomalies from 1981-2010 (m s"+r'$^{-1}$)', figtext = "(v) Surface Winds", scatter = (lons, lats, anom))
    utils.plot_smooth_map_iris(image_loc + "WND_{}_merra_anomaly".format(settings.YEAR), cube[0], settings.COLOURMAP_DICT["circulation_r"], bounds, "Anomalies from 1981-2010 (m s"+r'$^{-1}$)')

    #************************************************************************
    # ERA + HadISD Trend figure

    # Read in ERA anomalies

    cube_list = iris.load(reanalysis_loc + "ERA-Int_trendsb_10SI_ann_1979_{}.nc".format(settings.YEAR))

    cube = cube_list[0]
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()

    bounds=[-4, -0.8, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4, 0.8, 4]
    bounds=[-100, -0.4, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.4, 100]
    utils.plot_smooth_map_iris(image_loc + "WND_{}_era_trend".format(settings.YEAR), cube[0][0], settings.COLOURMAP_DICT["circulation_r"], bounds, "Trend from 1979-{} (m s".format(settings.YEAR)+r'$^{-1}$'+" decade"+r'$^{-1}$)')
    utils.plot_smooth_map_iris(image_loc + "WND_{}_era_obs_trend".format(settings.YEAR), cube[0][0], settings.COLOURMAP_DICT["circulation_r"], bounds, "Trend from 1979-{} (m s".format(settings.YEAR)+r'$^{-1}$'+" decade"+r'$^{-1}$)', scatter = (lons, lats, trend))


    #************************************************************************
    # Ocean timeseries

    satellite = read_radiometer(data_loc + "wind_data_for_global_ocean_time_series.annual.txt")

#    print "NO IN SITU OCEAN DATA FOR 2016, using 2015 data"
#    nocs = read_ts_cube(data_loc + "NOCSv2.0_oceanW_5by5_8110anoms_areaTS_FEB2016.nc", "Globally Average 70S-70N", "NOCSv2.0")

#    WASwind = read_ts_cube(data_loc + "waswind_v1_0_1.monthly_areaTS_19502011.nc","Globally Averaged Anomalies 70S-70N", "WASwind")

#    print "FIXING WASWIND TIMES - DATAFILE HAS WRONG DESCRIPTOR"
#    WASwind.times = WASwind.times - (1973-1950)

    jra_actuals, jra_anoms = utils.read_jra55(reanalysis_loc + "JRA-55_ws10m_globalocean_ts.txt", "wind")

    era_globe, era_ocean, era_land, eratropics = utils.era_ts_read(reanalysis_loc, "wnd", annual = True)
    ocean_era_clim, ocean_era_anoms = utils.calculate_climatology_and_anomalies_1d(era_ocean, 1981, 2010)

    merra_anoms = utils.read_merra(reanalysis_loc + "MERRA-2_SfcAnom_ts_{}.dat".format(settings.YEAR), "wind", "O", anomalies = True)

#    fig, (ax1, ax2, ax3) = plt.subplots(3, figsize = (10,12), sharex=True)
    fig, (ax1, ax3) = plt.subplots(2, figsize = (10,8), sharex=True)

    # Satellite
    utils.plot_ts_panel(ax1, [satellite], "-", "circulation", loc = LEGEND_LOC, bbox = BBOX)

    # In Situ
#    utils.plot_ts_panel(ax2, [nocs, WASwind], "-", "circulation", loc = LEGEND_LOC, bbox = BBOX)

#    ax2.set_ylabel("Anomaly (m s"+r'$^{-1}$'+")", fontsize = settings.FONTSIZE)

    # Reanalyses
    utils.plot_ts_panel(ax3, [ocean_era_anoms, jra_anoms, merra_anoms], "-", "temperature", loc = LEGEND_LOC, bbox = BBOX)

    #*******************
    # prettify

    # sort formatting
    plt.xlim([1979, int(settings.YEAR) + 1])

    for tick in ax3.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 
#    for ax in [ax1, ax2, ax3]:
    for ax in [ax1, ax3]:
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 
        ax.set_ylim(YLIM)
        ax.yaxis.set_ticks([-0.2,0,0.2,0.4])
        ax.yaxis.set_ticks_position('left')

    # sort labelling
    ax1.text(0.02, 0.9, "(a) Satellites", transform = ax1.transAxes, fontsize = settings.LABEL_FONTSIZE)
#    ax2.text(0.02, 0.9, "(b) In Situ", transform = ax2.transAxes, fontsize = settings.LABEL_FONTSIZE)
    ax3.text(0.02, 0.9, "(b) Reanalyses", transform = ax3.transAxes, fontsize = settings.LABEL_FONTSIZE)

    fig.subplots_adjust(right = 0.95, top = 0.95, hspace = 0.001)

    plt.savefig(image_loc+"WND_ocean_ts{}".format(settings.OUTFMT))

    plt.close()

    #************************************************************************
    # Ocean maps


    anoms = read_map_data(data_loc + "data_for_{}_wind_anom_map.merra2.1981.2010.txt".format(settings.YEAR))
    bounds=[-40, -1.2, -0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.8, 1.2, 40]

    utils.plot_smooth_map_iris(image_loc + "WND_ocean_merra2", anoms, settings.COLOURMAP_DICT["circulation_r"], bounds, "Anomaly (m s"+r'$^{-1}$'+")")


    #************************************************************************
    # MERRA/RSS ocean + HadISD Trend figure

    # Read in MERRA/RSS anomalies

    anoms = read_map_data(data_loc + "data_for_1988_{}_wind_trend_map.rss_ocean_merra2_land.txt".format(settings.YEAR))

    bounds=[-4, -0.8, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4, 0.8, 4]
    bounds=[-100, -0.4, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.4, 100]
    utils.plot_smooth_map_iris(image_loc + "WND_{}_merra-rss_trend".format(settings.YEAR), anoms, settings.COLOURMAP_DICT["circulation_r"], bounds, "Trend from 1988-{} (m s".format(settings.YEAR)+r'$^{-1}$'+" decade"+r'$^{-1}$)')
    utils.plot_smooth_map_iris(image_loc + "WND_{}_merra-rss_obs_trend".format(settings.YEAR), anoms, settings.COLOURMAP_DICT["circulation_r"], bounds, "Trend from 1988-{} (land/ice and ocean) and 1979-{} (land, points) (m s".format(settings.YEAR, settings.YEAR)+r'$^{-1}$'+" decade"+r'$^{-1}$)', scatter = (lons, lats, trend))

    #*******************************************************
    # twin plot

    late_1997 = read_map_data(data_loc + "data_for_late_1997_wind_anom_map.era_int.1988.2010.txt")
    late_2015 = read_map_data(data_loc + "data_for_late_2015_wind_anom_map.era_int.1988.2010.txt")

    bounds=[-40, -1.2, -0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.8, 1.2, 40]

    utils.plot_smooth_map_iris_multipanel(image_loc + "WND_ocean_1997_2015", [late_1997, late_2015], settings.COLOURMAP_DICT["circulation_r"], bounds, "Anomaly (m s"+r'$^{-1}$'+")", shape = (2,1))

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 End
#************************************************************************
