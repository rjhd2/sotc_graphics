#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for P extremes (PEX) section.
#       For BAMS SotC 2017 
#
#************************************************************************
#                    SVN Info
# $Rev:: 21                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2017-12-22 11:57:17 +0000 (Fri, 22 Dec #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import os
import copy
import struct
import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
from matplotlib.ticker import MultipleLocator

import iris
import cartopy

import utils # RJHD utilities
import settings

DATALOC = "{}/{}/data/PEX/".format(settings.ROOTLOC, settings.YEAR)

DATASTART=1951
CLIMSTART = 1961
CLIMEND = 1990
SELECTED_YEAR = int(settings.YEAR)
THRESHOLD = 0.9

LEGEND_LOC = 'upper left'
BBOX = (0, 0.9)
LW = 3


SEASON_DICT = {"MAM" : ["Mar", "Apr", "May"], "JJA" : ["Jun", "Jul", "Aug"], \
                   "SON" : ["Sep", "Oct", "Nov"], "DJF" : ["Dec", "Jan", "Feb"]}

SEASON_DICT_ERA = {"MAM" : ["March", "April", "May"], "JJA" : ["June", "July", "August"], \
                   "SON" : ["September", "October", "November"], "DJF" : ["December", "January", "February"]}

SEASONS = ["DJF", "MAM", "JJA", "SON"]

ETCCDI_INDICES = ["PRCPTOT", "Rx1day", "Rx5day", "R20mm", "R95p"] # "R10mm"

ETCCDI_LABELS = {"PRCPTOT" : "Total Precipitation", "Rx1day" : "Maximum 1 day precipitation total", "Rx5day" : "Maximum 5 day precipitation total", "R10mm" : "Number of heavy precipitation days", "R20mm" : "Number of very heavy precipitation days", "R95p" : "Precipitation from very wet days"}

ETCCDI_UNITS = {"PRCPTOT" : "mm", "Rx1day" : "mm", "Rx5day" : "mm", "R10mm" : "days", "R20mm" : "days", "R95p" : "mm"}


DWD_INDICES = ["PD10", "PD20", "RX1", "RX5", "R95P", "DI"]
DWD_NAMES = {"PD10" : "R10mm", "PD20" : "R20mm", "RX1" : "Rx1day", "RX5" : "Rx5day", "R95P" : "R95p", "DI" : "Drought"}

DWD_LABELS = {"CDD" : "Consecutive dry days", "CWD" : "Consecutive wet days", "DI" : "Drought Index", "DD" : "Number of dry days", "PD" : "Number of wet days", "RX1" : "Maximum 1 day precipitation total", "RX5" : "Maximum 5 day precipitation total", "PD10" : "Number of heavy precipitation days (>10mm)", "PD20" : "Number of very heavy precipitation days (>20mm)", "R95P" : "Average precipitation from very wet days", "SDII" : "Specific Daily Intensity Index"}

DWD_UNITS = {"CDD" : "days", "CWD" : "days", "DI" : "DI", "DD" : "days", "PD" : "days", "RX1" : "mm", "RX5" : "mm", "PD10" : "days", "PD20" : "days", "R95P" : "mm/day", "SDII" : "mm/day"}

#***************************************
def binomial(n, k): 
    bc = [1 for i in range(0, k+1)] 
    for j in range(1, n-k+1): 
        for i in range(1, k+1): 
            bc[i] = bc[i-1] + bc[i] 
    return bc[k] #binomial


#***************************************
def binomialfilter(tosmooth, mdi, n, pad=False):
    '''

    data - data array
    mdi - missing data indicator
    n - filter size
    '''
    

    smoothed = np.zeros(len(tosmooth))
    smoothed.fill(mdi)

    weights = []
    for k in range(n):
        weights += [binomial(n-1, k)]

    weights = np.array(weights)

    for o, obs in enumerate(tosmooth):
        if (o >= int(n/2)) and (o <= len(tosmooth)-int(n/2)):
            
            chunk = tosmooth[o-int(n/2): o+int(n/2)+1]
            good = np.where(chunk != mdi)
            
            if len(good[0]) >= int(n/2)+1:
                
                norm = sum(weights[good])
                weighted = sum(chunk[good]*weights[good]/norm)

                smoothed[o] = weighted

        elif pad == True:

            if (o < int(n/2)):
                smoothed[o] = np.mean(tosmooth[:int(n/2)])
            elif (o > len(tosmooth)-n/2):
                smoothed[o] = np.mean(tosmooth[int(-n/2):])
            

    return smoothed #binomialfilter

#************************************************************************
# def GetYears(cube):

#     times = cube.coord('time').points
#     years = np.round(np.array([(t - 101)/10000 for t in times]))

#     return years

#************************************************************************
def GetYears(cube, is_era5 = False):

    if is_era5:

        timeUnits = cube.coord("time").units
        dt_time = timeUnits.num2date(cube.coord("time").points)

        years = np.array([d.year for d in dt_time])

    else:
        times = cube.coord('time').points
        years = times + DATASTART
    
    return years

#************************************************************************
def ApplyClimatology(cube):

    years = GetYears(cube)

    clim_years=np.where((years >= CLIMSTART) & (years <= CLIMEND))
    
    climatology = np.ma.mean(cube.data[clim_years], axis=0)
    
    cube.data = cube.data - climatology                                
        
    return cube

#************************************************************************
def obtain_timeseries(filename, cube_name, ts_name):

    # re-read the file to overwrite memory of what was done to this cube.
    cube_list = iris.load(filename)
    names = np.array([cube.name() for cube in cube_list])

    selected_cube, = np.where(names == cube_name)

    cube = cube_list[selected_cube]
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()  

    nyears = np.ma.count(cube.data, axis=0)

    # repeat so selection will work across all years
    nyears = np.array([nyears for i in range(cube.data.shape[0])])

    cube.data = np.ma.masked_where(nyears < THRESHOLD * cube.data.shape[0], cube.data)

    weight_areas = iris.analysis.cartography.area_weights(cube)
    ts = cube.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=weight_areas)

    ts = utils.Timeseries(ts_name, GetYears(ts), ts.data * 3.65)

    return ts

#*********************************************
def parse(string, fields):
    """
    Replacement for struct and parser in Py 3
    """

    cumulative = 0
    outlist = []

    for field in fields:
        outlist += [string[cumulative: cumulative + field]]
        cumulative += field

    return np.array(outlist) # parse

#************************************************************************
def read_ghcnd(filename):

    fieldwidths = (11, 5, 3, 6, 9, 8, 5, 6, 3, 8, 8, 10, 10, 7, 25)

    indata = []
    with open(filename, 'r', encoding="latin-1") as infile:
        for ll, line in enumerate(infile):
            fields = parse(line, fieldwidths)
               
            indata += [[f.strip() for f in fields]]

    indata = np.array(indata)

    ids = indata[:, 0]
    record_y = indata[:, 1].astype(int)
    record_m = indata[:, 2].astype(int)
    ratio = indata[:, 3].astype(float)
    value_mm = indata[:, 4].astype(float)
    value_in = indata[:, 5].astype(float)
    data_length = indata[:, 6].astype(int)
    prev_record_yr = indata[:, 7].astype(int)
    prev_record_mth = indata[:, 8].astype(int)
    prev_value_mm = indata[:, 9].astype(float)
    prev_value_in = indata[:, 10].astype(float)
    lats = indata[:, 11].astype(float)
    lons = indata[:, 12].astype(float)
    elev = indata[:, 13].astype(float)
    name = indata[:, 14]


    return lats[::-1], lons[::-1], ratio[::-1], value_mm[::-1], prev_value_mm[::-1] # read_ghcnd

#************************************************************************
def calc_rank(data):
    
    order = np.argsort(data)

#    return np.argsort(order) # calc_rank

    # derive plotting array

    # get lowest/highest 4 and mark with -4 --> -1, and 1 --> 4 for plotting
    N = 5
    plot_order = np.argsort(order)
    if len(data) > 1:

        for rank in range(len(plot_order)):

            loc, = np.where(plot_order == rank)

            if rank in np.arange(0, N+1):
                plot_order[loc] = rank-N+1
                
            elif rank in np.arange(len(data)-N, len(data)):
                plot_order[loc] = rank - (len(data)-N)

            else:
                plot_order[loc] = 0

    return plot_order # calc_rank



#************************************************************************
def get_ranks(incube):


    lats = incube.coord('latitude').points
    lons = incube.coord('longitude').points
    time = incube.coord('time').points

    data = incube.data
    data = np.swapaxes(data, 1, 2) # swap lats and lons around

    rank = np.ma.zeros(data.shape)
    rank.mask = data.mask
    
    for lt in range(len(lats)):
        for ln in range(len(lons)):
            clean_loc, = np.where(data.mask[:, ln, lt] == False)
            rank[clean_loc, ln, lt] = calc_rank(data[clean_loc, ln, lt])

    cube = utils.make_iris_cube_3d(rank, time, "unknown", lons, lats, "TEX rank", "1")

    return cube # get_ranks

#************************************************************************
def plot_rank_map(outname, cube, cmap, bounds, cb_label, scatter=[], figtext="", title=""):
    '''
    Standard scatter map

    :param str outname: output filename root
    :param array cube: cube to plot
    :param obj cmap: colourmap to use
    :param array bounds: bounds for discrete colormap
    :param str cb_label: colorbar label
    '''

    norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)

    fig = plt.figure(figsize=(8, 5.5))

    plt.clf()
    ax = plt.axes([0.01, 0.12, 0.98, 0.88], projection=cartopy.crs.Robinson())
    ax.gridlines() #draw_labels=True)
    ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
    ax.coastlines()

    ext = ax.get_extent() # save the original extent

    mesh = iris.plot.pcolormesh(cube, cmap=cmap, norm=norm)

    if len(scatter) > 0:

        lons, lats, data = scatter

        plt.scatter(lons, lats, c=data, cmap=cmap, norm=norm, s=25, \
                        transform=cartopy.crs.Geodetic(), edgecolor='0.5', linewidth=0.5)


    cb = plt.colorbar(mesh, orientation='horizontal', pad=0.05, fraction=0.05, \
                        aspect=30, ticks=bounds[1:-1], drawedges=True)

    # thicken border of colorbar and the dividers
    # http://stackoverflow.com/questions/14477696/customizing-colorbar-border-color-on-matplotlib
    cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
    cb.ax.tick_params(axis='x', labelsize=settings.FONTSIZE, direction='in', size=0)

    cb.set_label(label=cb_label, fontsize=settings.FONTSIZE)

#    cb.outline.set_color('k')
    cb.outline.set_linewidth(2)
    cb.dividers.set_color('k')
    cb.dividers.set_linewidth(2)

    # label colorbar with sensibly placed and named labels
    cb.ax.get_xaxis().set_ticks([])
    cb.ax.get_xaxis().set_ticklabels([""])
    for j, lab in enumerate(['lowest', '2nd lowest', '3rd lowest', '', '3rd highest', '2nd highest', 'highest']):
        cb.ax.text((2 * j + 1) / 14.0, -0.5, lab, ha='center', va='center', fontsize=settings.FONTSIZE*0.8)

    ax.set_extent(ext, ax.projection) # fix the extent change from colormesh

    plt.title(title)
    fig.text(0.03, 0.95, figtext, fontsize=settings.FONTSIZE)

    plt.savefig(outname + settings.OUTFMT)
    plt.close()

    return # plot_rank_map_iris


#************************************************************************
def read_dwd_percentile_old(filename):
    """
    Read data from .txt file into Iris cube

    :param str filename: file to process

    :returns: cube
    """
    # use header to hard code the final array shapes
 
    longitudes = np.arange(-179.5, 180.5, 1.)
    latitudes = np.arange(89.5, -90.5, -1.)
    
    data = np.ma.zeros((latitudes.shape[0], longitudes.shape[0]))

    # read in the dat
    indata = np.genfromtxt(filename, dtype=(float))

    this_lat = []
    tl = 0
    # process each row, append until have complete latitude band
    for row in indata:
        this_lat += [row]

        if len(this_lat) == longitudes.shape[0]:
            # copy into final array and reset
            data[tl, :] = this_lat
            tl += 1
            this_lat = []

    # mask the missing values
    data = np.ma.masked_where(data <= -999.000, data)

    cube = utils.make_iris_cube_2d(data, latitudes, longitudes, "R90p", "%")

    return cube # read_dwd_percentile_old

#************************************************************************
def read_dwd_percentile(filename):
    """
    Read data from .txt file into Iris cube

    :param str filename: file to process

    :returns: cube
    """
    # use header to hard code the final array shapes
 
    # read in the dat
    indata = np.genfromtxt(filename, dtype=(float))

    longitudes = np.unique(indata[:, 0])
    latitudes = np.unique(indata[:, 1])

    data = np.ma.zeros((latitudes.shape[0], longitudes.shape[0]))
    data.mask = np.ones(data.shape)

    for val in indata:

        yloc, = np.where(longitudes == val[0])
        xloc, = np.where(latitudes == val[1])

        data[xloc[0], yloc[0]] = val[2]

    # mask the missing values
    data = np.ma.masked_where(data <= -99999.99, data)

    cube = utils.make_iris_cube_2d(data, latitudes, longitudes, "R90p", "%")

    return cube # read_dwd_percentile

#************************************************************************
def read_era5(filename, widths=7):

    if widths == 7:
        fieldwidths = (7, 7, 7, 7, 7, 7, 7, 7, 7, 7)
        longitudes = np.arange(0, 360, 0.28125)
    elif widths == 6:
        fieldwidths = (6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6)
        longitudes = np.arange(0, 360, 0.28302) #  wrong in the header

    try:
        with open(DATALOC + filename, "r", encoding="latin-1") as infile:

            latitudes = []
            data = []
            for lc, line in enumerate(infile):
                # skip header
                if lc < 2:
                    continue
                
                sp_line = line.split()
                if sp_line[0] == "latitude":
                    # read in last data
                    if lc > 3:
                        data += [values]

                    latitudes += [float(sp_line[1])]
                    values = []
                else:
                    values += [float(t) for t in parse(line, fieldwidths)]

            # and final one
            data += [values]    
            data = np.array(data)
            latitudes = np.array(latitudes)

    except IOError:
        print("{} doesn't exist".format(DATALOC + filename))


    cube = utils.make_iris_cube_2d(data, latitudes, longitudes, "max Precipitation", "mm")

    return cube

#************************************************************************
def read_cei(filename):

    indata = np.genfromtxt(DATALOC + filename, delimiter=',', encoding='latin-1')

    return utils.Timeseries("CEI", indata[:, 0], indata[:, 1]) # read_cei


#
#************************************************************************
def run_all_plots():

    # # CEI timeseries
    if False:
        # from https://www.ncdc.noaa.gov/extremes/cei/graph/us/03-05/4 (Spring, Step4 indicator)
        cei = read_cei("CEI_step4_figure_data.csv")
        smoothed = binomialfilter(cei.data, -99.9, 9, pad = False)
        smoothed = np.ma.masked_where(smoothed == -99.9, smoothed)

        fig = plt.figure(figsize=(8, 5))
        plt.clf()
        ax1 = plt.axes([0.11, 0.08, 0.86, 0.90])

        ax1.bar(cei.times, cei.data, color="g", label=cei.name, align="center", width=1, edgecolor="darkgreen")
        ax1.plot(cei.times, smoothed, "r", lw = LW)

        ax1.plot([cei.times[0], cei.times[-1]], [np.mean(cei.data), np.mean(cei.data)], "k", lw=1)

        ax1.set_xlim([1910, int(settings.YEAR)+2])
        ax1.set_ylabel("%", fontsize=settings.FONTSIZE)
        minorLocator = MultipleLocator(1)
        ax1.xaxis.set_minor_locator(minorLocator)

        utils.thicken_panel_border(ax1)
        ax1.yaxis.set_tick_params(right=False)
        for tick in ax1.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        for tick in ax1.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)

#        ax1.text(0.02, 0.9, "(e)", transform=ax1.transAxes, fontsize=settings.FONTSIZE)
        plt.savefig(settings.IMAGELOC+"PEX_CEI_ts{}".format(settings.OUTFMT))

        plt.close()

    #*************************
    # GHCNDEX indices
    rank_bounds = [-4.5, -3.5, -2.5, -1.5, 1.5, 2.5, 3.5, 4.5]

    if True:

        for index in ETCCDI_INDICES:


            cube_list = iris.load(DATALOC + "GHCND_{}_1951-{}_RegularGrid_global_2.5x2.5deg_LSmask.nc".format(index, int(settings.YEAR) + 1))
            names = np.array([cube.name() for cube in cube_list])

            #*************
            # plot annual map

            selected_cube, = np.where(names == "Ann")

            total_cube = cube_list[selected_cube[0]]
            total_cube.coord('latitude').guess_bounds()
            total_cube.coord('longitude').guess_bounds()  

            anoms = copy.deepcopy(total_cube)
            anoms = ApplyClimatology(anoms)

            # select the year to plot
            years = GetYears(total_cube)
            loc, = np.where(years == SELECTED_YEAR)

            # sort the bounds and colourbars
            if index in ["Rx1day"]:
                bounds = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
                cmap = settings.COLOURMAP_DICT["precip_sequential"]
            elif index in ["Rx5day"]:
                bounds = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
                cmap = settings.COLOURMAP_DICT["precip_sequential"]
            elif index in ["R10mm"]:
                bounds = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
                cmap = settings.COLOURMAP_DICT["precip_sequential"]
            elif index in ["R20mm"]:
                bounds = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
                cmap = settings.COLOURMAP_DICT["precip_sequential"]
            elif index in ["R95p"]:
                bounds = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225]
                cmap = settings.COLOURMAP_DICT["precip_sequential"]
            elif index in ["PRCPTOT"]:
                bounds = [0, 25, 50, 75, 100, 125, 150, 200, 300, 400]
                cmap = settings.COLOURMAP_DICT["precip_sequential"]

            utils.plot_smooth_map_iris(settings.IMAGELOC + "PEX_{}_{}_ghcndex".format(index, settings.YEAR), total_cube[loc[0]], cmap, bounds, "{} ({})".format(index, ETCCDI_UNITS[index]), title="GHCNDEX {} - {}".format(index, ETCCDI_LABELS[index]))

            # sort the bounds and colourbars
            if index in ["Rx1day"]:
                bounds = [-100, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 100]
                bounds = [-100, -40, -30, -20, -10, -5, 0, 5, 10, 20, 30, 40, 100]
                cmap = settings.COLOURMAP_DICT["hydrological"]
            elif index in ["Rx5day"]:
                bounds = [-100, -40, -30, -20, -10, -5, 0, 5, 10, 20, 30, 40, 100]
                cmap = settings.COLOURMAP_DICT["hydrological"]
            elif index in ["R10mm"]:
                bounds = [-20, -8, -6, -4, -2, 0, 2, 4, 6, 8, 20]
                cmap = settings.COLOURMAP_DICT["hydrological"]
            elif index in ["R20mm"]:
                bounds = [-20, -4, -3, -2, -1, 0, 1, 2, 3, 4, 20]
                cmap = settings.COLOURMAP_DICT["hydrological"]
            elif index in ["R95p"]:
                bounds = [-1000, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 1000]
                cmap = settings.COLOURMAP_DICT["hydrological"]
            elif index in ["PRCPTOT"]:
                bounds = [-1000, -150, -80, -60, -40, -20, 0, 20, 40, 60, 80, 150, 1000]
                cmap = settings.COLOURMAP_DICT["hydrological"]

            utils.plot_smooth_map_iris(settings.IMAGELOC + "PEX_{}_{}_anoms_ghcndex".format(index, settings.YEAR), anoms[loc[0]], cmap, bounds, "Anomalies from 1961-90 ({})".format(ETCCDI_UNITS[index]), title="GHCNDEX {} - {}".format(index, ETCCDI_LABELS[index]))
            if index == "Rx1day":
                utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_PEX_{}_{}_anoms_ghcndex".format(index, settings.YEAR), anoms[loc[0]], cmap, bounds, "Anomalies from 1961-90 ({})".format(ETCCDI_UNITS[index]), figtext="(k) Maximum 1 Day Precipitation Amount")

            rank_cube = get_ranks(anoms)

            plot_rank_map(settings.IMAGELOC + "PEX_{}_{}_rank_ghcndex".format(index, settings.YEAR), rank_cube[loc[0]], cmap, rank_bounds, "Rank", title="GHCNDEX {} - {}".format(index, ETCCDI_UNITS[index]))
            #*************
            # plot season maps (2x2)

            if index in ["Rx1day", "Rx5day"]:

                for sc, cube in enumerate([total_cube, anoms]):

                    season_list = []
                    for season in SEASONS:

                        # extract each month
                        month_data = []
                        months = SEASON_DICT[season]
                        for month in months:

                            selected_cube, = np.where(names == month)
                            cube = cube_list[selected_cube[0]]

                            if month == "Dec":
                                # need to extract from previous year - cheat by rolling data around
                                cube.data = np.roll(cube.data, 1, axis=0)
                                cube.data.mask[0, :, :] = True # and mask out the previous years'

                            month_data += [cube.data]

                        # finished getting all months, make a dummy cube to populate
                        month_data = np.ma.array(month_data)
                        season_cube = copy.deepcopy(cube)

                        # take appropriate seasonal value
                        season_cube.data = np.ma.max(month_data, axis=0)

                        # mask if fewer that 2 months present
                        nmonths_locs = np.ma.count(month_data, axis=0)
                        season_cube.data = np.ma.masked_where(nmonths_locs < 2, season_cube.data)

                        # make anomalies
                        if sc == 1:
                            season_cube = ApplyClimatology(season_cube)

                        # fix for plotting
                        season_cube.coord('latitude').guess_bounds()
                        season_cube.coord('longitude').guess_bounds()

                        # select the year to plot
                        years = GetYears(cube)
                        loc, = np.where(years == SELECTED_YEAR)

                        # add to list
                        season_list += [season_cube[loc[0]]]

                    # sort the bounds and colourbars
                    if sc == 0:
                        if index in ["Rx1day", "Rx5day"]:
                            bounds = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
                            cmap = settings.COLOURMAP_DICT["precip_sequential"]

                        # pass to plotting routine
                        utils.plot_smooth_map_iris_multipanel(settings.IMAGELOC + "PEX_{}_{}_seasons_ghcndex".format(index, settings.YEAR), season_list, cmap, bounds, "GHCNDEX {} ({})".format(index, ETCCDI_UNITS[index]), shape=(2, 2), title=SEASONS, figtext=["(a)", "(b)", "(c)", "(d)"], figtitle="{} - {}".format(index, ETCCDI_LABELS[index]))
                    elif sc == 1:
                        cmap = settings.COLOURMAP_DICT["hydrological"]
                        if index in ["Rx1day"]:
                            bounds = [-100, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 100]
                        elif index in ["Rx5day"]:
                            bounds = [-100, -20, -16, -12, -8, -4, 0, 4, 8, 12, 16, 20, 100]

                        # pass to plotting routine
                        utils.plot_smooth_map_iris_multipanel(settings.IMAGELOC + "PEX_{}_{}_anoms_seasons_ghcndex".format(index, settings.YEAR), season_list, cmap, bounds, "Anomalies from 1961-90 ({})".format(ETCCDI_UNITS[index]), shape=(2, 2), title=SEASONS, figtext=["(a)", "(b)", "(c)", "(d)"], figtitle="GHCNDEX {} - {}".format(index, ETCCDI_LABELS[index]))

    
    #*************************
    # DWD indices
    if True:
        for index in DWD_INDICES:
            print(index)
            if not os.path.exists(DATALOC + "First_Guess_Daily_{}_{}.nc".format(settings.YEAR, index)):
                print("File {} missing".format("First_Guess_Daily_{}_{}.nc".format(settings.YEAR, index)))
                continue
            cube_list = iris.load(DATALOC + "First_Guess_Daily_{}_{}.nc".format(settings.YEAR, index))

            if len(cube_list) == 1:
                cube = cube_list[0]
            else:
                # these have two fields
                for cube in cube_list:
                    if cube.units == "mm per 5 day" and index == "RX5":
                        break
                    elif cube.var_name == "consecutive_wet_days_index_per_time_period" and index == "CWD":
                        break
                    elif cube.var_name == "consecutive_dry_days_index_per_time_period" and index == "CDD":
                        break
                    else:
                        print("Check cube for {} for extra fields".format(index))

            cube = cube[0] # take only single slice
            cube.coord('latitude').guess_bounds()
            cube.coord('longitude').guess_bounds()  

            if index in ["CDD"]:
                bounds = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
                cmap = settings.COLOURMAP_DICT["precip_sequential_r"]
            elif index in ["CWD"]:
                bounds = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
                bounds = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
                cmap = settings.COLOURMAP_DICT["precip_sequential"]
            elif index in ["DD"]:
                bounds = [0, 40, 80, 120, 160, 200, 240, 280, 320, 360]
                cmap = settings.COLOURMAP_DICT["precip_sequential_r"]
            elif index in ["PD"]:
                bounds = [0, 40, 80, 120, 160, 200, 240, 280, 320, 360]
                cmap = settings.COLOURMAP_DICT["precip_sequential"]
            elif index in ["RX1", "PD10"]:
                bounds = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
                cmap = settings.COLOURMAP_DICT["precip_sequential"]
            elif index in ["RX5"]:
                bounds = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450]
                bounds = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
                cmap = settings.COLOURMAP_DICT["precip_sequential"]
            elif index in ["R10", "R20", "R95P", "SDII", "PD20"]:
                bounds = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
                cmap = settings.COLOURMAP_DICT["precip_sequential"]

            utils.plot_smooth_map_iris(settings.IMAGELOC + "PEX_{}_{}_dwd".format(index, settings.YEAR), cube, cmap, bounds, "({})".format(DWD_UNITS[index]), title="GPCC {} - {}".format(DWD_NAMES[index], DWD_LABELS[index]))


    #*************************
    # DWD differences indices
    if True:
        for index in DWD_INDICES:
            print(index)
            if not os.path.exists(DATALOC + "Diff_{}-Mean_{}.nc".format(index, settings.YEAR)):
                print("File {} missing".format("Diff_{}-Mean_{}.nc".format(index, settings.YEAR)))
                continue
            cube_list = iris.load(DATALOC + "Diff_{}-Mean_{}.nc".format(index, settings.YEAR))

            if len(cube_list) == 1:
                cube = cube_list[0]
            else:
                # these have two fields
                for cube in cube_list:
                    if cube.units == "mm per 5 day" and index == "RX5":
                        break

            cube = cube[0] # take only single slice
            cube.coord('latitude').guess_bounds()
            cube.coord('longitude').guess_bounds()  

            if index in ["RX1"]:
                bounds = [-100, -50, -25, -10, -5, 0, 5, 10, 25, 50, 100]
                cmap = settings.COLOURMAP_DICT["hydrological"]
            elif index in ["RX5"]:
                bounds = [-50, -30, -20, -10, -5, 0, 5, 10, 20, 30, 50]
                bounds = [-100, -50, -25, -10, -5, 0, 5, 10, 25, 50, 100]
                cmap = settings.COLOURMAP_DICT["hydrological"]
            elif index in ["PD10", "PD20", "R95P"]:
                bounds = [-50, -30, -20, -10, -5, 0, 5, 10, 20, 30, 50]
                cmap = settings.COLOURMAP_DICT["hydrological"]

            utils.plot_smooth_map_iris(settings.IMAGELOC + "PEX_{}_{}_diff_dwd".format(index, settings.YEAR), cube, cmap, bounds, "Anomalies from 1982-2016 ({})".format(DWD_UNITS[index]), title="GPCC {} - {}".format(DWD_NAMES[index], DWD_LABELS[index]))
            if index == "PD10":
                utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_PEX_{}_{}_diff_dwd".format("R10mm", settings.YEAR), cube, cmap, bounds, "Anomalies from 1982-2016 ({})".format(DWD_UNITS[index]), figtext="(l) {} anomalies".format(DWD_NAMES[index]))
            else:
                utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_PEX_{}_{}_diff_dwd".format(index, settings.YEAR), cube, cmap, bounds, "Anomalies from 1982-2016 ({})".format(DWD_UNITS[index]), figtext="(l) {} anomalies".format(DWD_NAMES[index]))


    #*************************
    # MERRA map Rx1day
    if True:
        index = "Rx1day"
        cube = iris.load(DATALOC + "MERRA2_ann_{}_rx1day_gl_anom.nc".format(settings.YEAR))[0]
        
        bounds = [-100, -60, -30, -10, -5, 0, 5, 10, 30, 60, 100]
        cmap = settings.COLOURMAP_DICT["hydrological"]
        
        utils.plot_smooth_map_iris(settings.IMAGELOC + "PEX_Rx1day_{}_anoms_merra2".format(settings.YEAR), cube[0], cmap, bounds, "Anomalies from 1981-2010 ({})".format(ETCCDI_UNITS[index]), title="MERRA-2 {} - {}".format(index, ETCCDI_LABELS[index]))
        # utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_PEX_Rx1day_{}_anoms_merra2".format(settings.YEAR), cube[0], cmap, bounds, "Anomalies from 1981-2010 ({})".format(ETCCDI_UNITS[index]), figtext="(l) Rx1day anomalies")

    #*************************
    # MERRA map R10mm
    if True:
        index = "R10mm"
        cube = iris.load(DATALOC + "MERRA2_ann_{}_r10mm_gl_anom.nc".format(settings.YEAR))[0]
        
        bounds = [-400, -30, -20, -10, -5, 0, 5, 10, 20, 30, 400]
        cmap = settings.COLOURMAP_DICT["hydrological"]
        
        utils.plot_smooth_map_iris(settings.IMAGELOC + "PEX_R10mm_{}_anoms_merra2".format(settings.YEAR), cube[0], cmap, bounds, "Anomalies from 1981-2010({})".format(ETCCDI_UNITS[index]), title="MERRA-2 {} - {}".format(index, ETCCDI_LABELS[index]))

     #*************************
    # MERRA map R20mm
    if True:
        index = "R20mm"
        cube = iris.load(DATALOC + "MERRA2_ann_{}_r20mm_gl_anom.nc".format(settings.YEAR))[0]
        
        bounds = [-400, -30, -20, -10, -5, 0, 5, 10, 20, 30, 400]
        cmap = settings.COLOURMAP_DICT["hydrological"]
        
        utils.plot_smooth_map_iris(settings.IMAGELOC + "PEX_R20mm_{}_anoms_merra2".format(settings.YEAR), cube[0], cmap, bounds, "Anomalies from 1981-2010 ({})".format(ETCCDI_UNITS[index]), title="MERRA-2 {} - {}".format(index, ETCCDI_LABELS[index]))

   #*************************
    # ERA5 map - R10mm anomaly map
    if True:
        index = "R10mm"
        cube = read_era5("precip_ndays_ge10_anomaly_{}.txt".format(settings.YEAR), widths=6)

        bounds = [-400, -30, -20, -10, -5, 0, 5, 10, 20, 30, 400]
        cmap = settings.COLOURMAP_DICT["hydrological"]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "PEX_R10mm_{}_anoms_era5".format(settings.YEAR), cube, cmap, bounds, "Anomalies from 1981-2010 ({})".format(ETCCDI_UNITS[index]), title="ERA5 {} - {}".format(index, ETCCDI_LABELS[index]))

    #*************************
    # ERA5 map - R20mm anomaly map
    if True:
        index = "R20mm"
        cube = read_era5("precip_ndays_ge20_anomaly_{}.txt".format(settings.YEAR), widths=6)

        bounds = [-400, -30, -20, -10, -5, 0, 5, 10, 20, 30, 400]
        cmap = settings.COLOURMAP_DICT["hydrological"]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "PEX_R20mm_{}_anoms_era5".format(settings.YEAR), cube, cmap, bounds, "Anomalies from 1981-2010 ({})".format(ETCCDI_UNITS[index]), title="ERA5 {} - {}".format(index, ETCCDI_LABELS[index]))

    #*************************
    # ERA5 map - Rx1day map
    if False:
        index = "Rx1day"
        cube = read_era5("prmax_{}.txt".format(settings.YEAR))
        
        bounds = [0, 2, 5, 10, 20, 40, 80, 160, 300, 450]
        cmap = settings.COLOURMAP_DICT["precip_sequential"]
        
        utils.plot_smooth_map_iris(settings.IMAGELOC + "PEX_Rx1day_{}_era5".format(settings.YEAR), cube, cmap, bounds, "{} ({})".format(index, ETCCDI_UNITS[index]), title="ERA5 {} - {}".format(index, ETCCDI_LABELS[index]))
 
    #*************************
    # ERA5 map - Rx1day anomaly map
    if False:
        index = "Rx1day"
        cube = read_era5("prmax1d_anomaly_for_{}_wrt_1981-2010.txt".format(settings.YEAR))

        bounds = [-400, -100, -75, -50, -25, 0, 25, 50, 75, 100, 400]
        bounds = [-400, -100, -50, -20, -10, 0, 10, 20, 50, 100, 400]
        bounds = [-100, -40, -30, -20, -10, -5, 0, 5, 10, 20, 30, 40, 100]
        cmap = settings.COLOURMAP_DICT["hydrological"]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "PEX_Rx1day_{}_anoms_era5".format(settings.YEAR), cube, cmap, bounds, "{} ({})".format(index, ETCCDI_UNITS[index]), title="ERA5 {} - {}".format(index, ETCCDI_LABELS[index]))
    #    utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_PEX_Rx1day_{}_anoms_era5".format(settings.YEAR), cube, cmap, bounds, "{} ({})".format(index, ETCCDI_UNITS[index]), figtext="(i) Rx1day anomalies")

    #*************************
    # ERA5 map - Rx1day anomalies as percent
    if True:
        index = "Rx1day"
        cube = read_era5("prmax1d_anomaly_percent_{}.txt".format(settings.YEAR))

        bounds = [0, 10, 25, 50, 75, 100, 150, 200, 250, 300, 1000]
        cmap = settings.COLOURMAP_DICT["hydrological"]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "PEX_Rx1day_{}_anoms_percent_era5".format(settings.YEAR), cube, cmap, bounds, "Anomalies from 1981-2010 (%)", title="ERA5 - Rx1day %")
        # utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_PEX_Rx1day_{}_anoms_percent_era5".format(settings.YEAR), cube, cmap, bounds, "Anomalies from 1981-2010 (%)", figtext="(m) Rx1day anomalies")

    #*************************
    # DWD percentile
    if False:
     
#        dwd_cube = read_dwd_percentile(DATALOC + "GPCC_perzentile_{}.xyzras".format(settings.YEAR))

        cube_list = iris.load(DATALOC + "Quantile_12month_{}01-{}12.nc".format(settings.YEAR, settings.YEAR))
        dwd_cube = cube_list[0][0]

        bounds = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        cmap = settings.COLOURMAP_DICT["hydrological"]
        utils.plot_smooth_map_iris(settings.IMAGELOC + "PEX_{}_{}_dwd".format("R90", settings.YEAR), dwd_cube, cmap, bounds, "GPCC {} ({})".format("percentile", "%"), title="Percentile of the annual precipitation total")
        # utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_PEX_{}_{}_dwd".format("R90", settings.YEAR), dwd_cube, cmap, bounds, "GPCC {} ({})".format("percentile", "%"), figtext="(l) Percentile of the Annual Precipitation Total")

    #*************************
    # DWD drought
    if False:
     
#        dwd_cube = read_dwd_percentile(DATALOC + "GPCC_perzentile_{}.xyzras".format(settings.YEAR))

        cube_list = iris.load(DATALOC + "GPCC_DI_201912_12.nc")
        dwd_cube = cube_list[0][0]

        bounds = [-10, -5, -3, -2, -1, 0, 1, 2, 3, 5, 10]
        cmap = settings.COLOURMAP_DICT["hydrological"]
        utils.plot_smooth_map_iris(settings.IMAGELOC + "PEX_{}_{}_dwd".format("DI", settings.YEAR), dwd_cube, cmap, bounds, "{} ({})".format("DI", "-"), title="12 month Drought Index")

    #*************************
    # GHCND totals and ratios
    if False:
        import cartopy.feature as cfeature

        NAMES = {"japan" : "Japan", "neaus" : "Australia", "hawaii" : "Hawaii"}
        EXTENTS = {"japan" : (132, 33, [126, 139, 30, 36]), "neaus" : (145, -17, [139, 151, -26, -8]), "hawaii" : (-158, 20, [-161, -154, 18.5, 22.5])}

        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor='face',
                                                facecolor=cfeature.COLORS['land'])
        land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                                edgecolor='face',
                                                facecolor=cfeature.COLORS['land'])

        for index in ["Rx5day", "Rx1day"]:
            for state in ["03-neaus", "07-japan", "08-hawaii"]:

                print("{} - {}".format(state, index))
                lats, lons, value_mm, prev_value_mm = read_ghcnd("{}/{}-{}-non0-sorted-{}.txt".format(DATALOC, index, settings.YEAR, state))

                ratio = value_mm/prev_value_mm

                # set up figure
                if "neaus" in state:
                    fig = plt.figure(figsize=(8, 5))
                else:
                    fig = plt.figure(figsize=(8, 3.5))


                if index == "Rx5day":
                    bounds = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900]
                elif index == "Rx1day":
                    bounds = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450]
                cmap = settings.COLOURMAP_DICT["precip_sequential"]
                norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)

    #            ratio_cmap = plt.cm.YlOrRd
    #            ratio_bounds = [1, 1.05, 1.1, 1.15, 1.2, 1.3, 1.5, 1.75, 2.0 ]
    #            ratio_cmap = settings.COLOURMAP_DICT["hydrological"]
    #            ratio_bounds = [0.01, 0.5, 0.75, 0.9, 0.95, 1.0, 1.05, 1.1, 1.25, 2.0, 10.0]
                ratio_cmap = plt.cm.Blues
                ratio_bounds = [0.01, 0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.5, 2.0, 10.0]
                ratio_norm = mpl.cm.colors.BoundaryNorm(ratio_bounds, ratio_cmap.N)

                plt.clf()

                # set up axes
                ax0 = plt.axes([0.01, 0.12, 0.45, 0.85], projection=cartopy.crs.Stereographic(central_longitude=EXTENTS[state[3:]][0], central_latitude=EXTENTS[state[3:]][1]))
                ax1 = plt.axes([0.51, 0.12, 0.45, 0.85], projection=cartopy.crs.Stereographic(central_longitude=EXTENTS[state[3:]][0], central_latitude=EXTENTS[state[3:]][1]))

                for ax in [ax0, ax1]:
                    ax.set_extent(EXTENTS[state[3:]][2], cartopy.crs.PlateCarree())

                    states_provinces = cfeature.NaturalEarthFeature(
                        category='cultural',
                        name='admin_1_states_provinces_lines',
                        scale='50m',
                        facecolor='none')
                    ax.add_feature(states_provinces, edgecolor='gray')
                    ax.coastlines(resolution="10m", linewidth=0.5)
                    ax.add_feature(land_50m, zorder=0, facecolor="0.9", edgecolor="k")

                    # add other features
                    ax.gridlines() #draw_labels=True)
                    ax.add_feature(cartopy.feature.BORDERS, zorder=0, facecolor="0.9", edgecolor="k")


                #  plot 
                for ax, data, cm, bnds, nrm, label in zip([ax0, ax1], [value_mm, ratio], [cmap, ratio_cmap], [bounds, ratio_bounds], [norm, ratio_norm], ["{} (mm)".format(index), "Ratio to previous record"]):

                    scatter = ax.scatter(lons, lats, c=data, cmap=cm, norm=nrm, s=50, \
                                             transform=cartopy.crs.Geodetic(), edgecolor='0.5', linewidth=0.5)


                    # thicken border of colorbar and the dividers
                    # http://stackoverflow.com/questions/14477696/customizing-colorbar-border-color-on-matplotlib
                    if "Ratio" in label:
                        cb = fig.colorbar(scatter, ax=ax, orientation='horizontal', pad=0.05, fraction=0.05, \
                                        aspect=30, ticks=bnds[1:-1], label=label, drawedges=True)
                        cb.set_ticklabels(["{:g}".format(b) for b in bnds[1:-1]])                
                    else:
                        cb = fig.colorbar(scatter, ax=ax, orientation='horizontal', pad=0.05, fraction=0.05, \
                                        aspect=30, ticks=bnds[1:-1], label=label, drawedges=True)
                        cb.set_ticklabels(["{:g}".format(b) for b in bnds])

                    cb.outline.set_linewidth(2)
                    cb.dividers.set_color('k')
                    cb.dividers.set_linewidth(2)
                    cb.ax.tick_params(axis='x', labelsize=settings.FONTSIZE*0.6, direction='in')

                if index == "Rx5day" and state[3:] == "hawaii":                
                    ax0.text(0.05, 1.05, "(a)", transform=ax0.transAxes, fontsize=settings.FONTSIZE*0.8)
                    ax1.text(0.05, 1.05, "(b)", transform=ax1.transAxes, fontsize=settings.FONTSIZE*0.8)
                elif index == "Rx1day" and state[3:] == "hawaii":
                    ax0.text(0.05, 1.05, "(c)", transform=ax0.transAxes, fontsize=settings.FONTSIZE*0.8)
                    ax1.text(0.05, 1.05, "(d)", transform=ax1.transAxes, fontsize=settings.FONTSIZE*0.8)


                plt.savefig(settings.IMAGELOC + "PEX_{}_{}-{}".format(index, NAMES[state[3:]], state[:2]) + settings.OUTFMT)
                plt.close()

    #*************************
    # GHCND totals and ratios
    if True:
        import cartopy.feature as cfeature

        NAMES = {"NEU" : "NEurope"}
        EXTENT = {"NEU" : (10, 57, [-10, 30, 45, 70])}

        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor='face',
                                                facecolor=cfeature.COLORS['land'])
        land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                                edgecolor='face',
                                                facecolor=cfeature.COLORS['land'])

        for index in ["Rx5day", "Rx1day"]:
            for state in ["NEU"]:

                print("{} - {}".format(state, index))

                lats, lons, ratio, value_mm, prev_value_mm = read_ghcnd("{}/{}-europe-{}-06.txt".format(DATALOC, index, settings.YEAR))

                # set up figure
                fig = plt.figure(figsize=(8, 5))

                if index == "Rx5day":
                    bounds = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900]
                    bounds = [0, 10, 20, 40, 60, 80, 100, 150, 200, 450]
                elif index == "Rx1day":
                    bounds = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450]
                    bounds = [0, 10, 20, 30, 40, 50, 60, 80, 100, 300]
                cmap = settings.COLOURMAP_DICT["precip_sequential"]
                norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)

                ratio_cmap = plt.cm.Blues
                ratio_bounds = [0.01, 0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.5, 2.0, 10.0]
                ratio_norm = mpl.cm.colors.BoundaryNorm(ratio_bounds, ratio_cmap.N)

                plt.clf()

                # set up axes
                ax0 = plt.axes([0.01, 0.12, 0.45, 0.85], projection=cartopy.crs.Stereographic(central_longitude=EXTENT[state][0], central_latitude=EXTENT[state][1]))
                ax1 = plt.axes([0.51, 0.12, 0.45, 0.85], projection=cartopy.crs.Stereographic(central_longitude=EXTENT[state][0], central_latitude=EXTENT[state][1]))

                for ax in [ax0, ax1]:
                    ax.set_extent(EXTENT[state][2], cartopy.crs.PlateCarree())

                    states_provinces = cfeature.NaturalEarthFeature(
                        category='cultural',
                        name='admin_1_states_provinces_lines',
                        scale='50m',
                        facecolor='none')
                    ax.add_feature(states_provinces, edgecolor='gray')
                    ax.coastlines(resolution="10m", linewidth=0.5)
                    ax.add_feature(land_50m, zorder=0, facecolor="0.9", edgecolor="k")

                    # add other features
                    ax.gridlines() #draw_labels=True)
                    ax.add_feature(cartopy.feature.BORDERS, zorder=0, facecolor="0.9", edgecolor="k")


                #  plot 
                for ax, data, cm, bnds, nrm, label in zip([ax0, ax1], [value_mm, ratio], [cmap, ratio_cmap], [bounds, ratio_bounds], [norm, ratio_norm], ["{} (mm)".format(index), "Ratio to previous record"]):

                    scatter = ax.scatter(lons, lats, c=data, cmap=cm, norm=nrm, s=50, \
                                             transform=cartopy.crs.Geodetic(), edgecolor='0.5', linewidth=0.5)


                    # thicken border of colorbar and the dividers
                    # http://stackoverflow.com/questions/14477696/customizing-colorbar-border-color-on-matplotlib
                    if "Ratio" in label:
                        cb = fig.colorbar(scatter, ax=ax, orientation='horizontal', pad=0.05, fraction=0.05, \
                                        aspect=30, ticks=bnds[1:-1], label=label, drawedges=True)
                        cb.set_ticklabels(["{:g}".format(b) for b in bnds[1:-1]])                
                    else:
                        cb = fig.colorbar(scatter, ax=ax, orientation='horizontal', pad=0.05, fraction=0.05, \
                                        aspect=30, ticks=bnds[1:-1], label=label, drawedges=True)
                        cb.set_ticklabels(["{:g}".format(b) for b in bnds])

                    cb.outline.set_linewidth(2)
                    cb.dividers.set_color('k')
                    cb.dividers.set_linewidth(2)
                    cb.ax.tick_params(axis='x', labelsize=settings.FONTSIZE*0.6, direction='in')

                if index == "Rx5day":                
                    ax0.text(0.05, 1.05, "(c)", transform=ax0.transAxes, fontsize=settings.FONTSIZE*0.8)
                    ax1.text(0.05, 1.05, "(d)", transform=ax1.transAxes, fontsize=settings.FONTSIZE*0.8)
                elif index == "Rx1day":
                    ax0.text(0.05, 1.05, "(a)", transform=ax0.transAxes, fontsize=settings.FONTSIZE*0.8)
                    ax1.text(0.05, 1.05, "(b)", transform=ax1.transAxes, fontsize=settings.FONTSIZE*0.8)


                plt.savefig(settings.IMAGELOC + "PEX_{}-{}".format(index, NAMES[state]) + settings.OUTFMT)
                plt.close()

    # correlation between Rx1day and R20mm - quick look 2020 - likely delete this snipped
    if False:
        rx_cube = iris.load(DATALOC + "MERRA2_ann_{}_rx1day_gl_anom.nc".format(settings.YEAR))[0]
        r2_cube = iris.load(DATALOC + "MERRA2_ann_{}_r20mm_gl_anom.nc".format(settings.YEAR))[0]
        
        correlations = np.corrcoef(rx_cube.data.reshape(-1), r2_cube.data.reshape(-1))

        plt.clf()
        plt.plot(rx_cube.data.reshape(-1), r2_cube.data.reshape(-1), marker=".", ls="", alpha=0.3)
        plt.text(0.8, 0.9, r'$r={:5.3f}$'.format(correlations[0,1]), transform=plt.gca().transAxes)
        plt.xlabel("Rx1day (mm)")
        plt.ylabel("R20mm (days)")
        plt.savefig(settings.IMAGELOC + "PEX_MERRA2_Rx1day_R20mm_correl.png")
        plt.close()

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
