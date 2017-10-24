#!/usr/local/sci/python
#************************************************************************
#
#  Plot figures and output numbers for T extremes (TEX) section.
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

import copy
import string

import iris
import iris.quickplot as qplt
import cartopy

import datetime as dt

import utils # RJHD utilities
import settings

data_loc = "/data/local/rdunn/SotC/{}/data/TEX/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

CLIMSTART=1961
CLIMEND=1990
SELECTED_YEAR = 2016
THRESHOLD = 0.9

LEGEND_LOC = 'upper left'
BBOX = (0,0.9)
LW = 3


SEASON_DICT = {"MAM":["Mar","Apr","May"], "JJA":["Jun","Jul","Aug"], \
                   "SON":["Sep","Oct","Nov"], "DJF":["Dec","Jan","Feb"]}

SEASON_DICT_ERA = {"MAM":["March","April","May"], "JJA":["June","July","August"], \
                   "SON":["September","October","November"], "DJF":["December","January","February"]}

SEASONS = ["DJF", "MAM", "JJA","SON"]

INDICES = ["TX90p", "TX10p", "TN90p", "TN10p"] #, "TXx", "TXn", "TNx", "TNn"]

INDEX_LABELS = {"TX90p" : "Warm Days", "TX10p" : "Cool Days", "TN90p" : "Warm Nights", "TN10p" : "Cool Nights"}

#***************************************
def binomial(n,k): 
    bc = [1 for i in range(0,k+1)] 
    for j in range(1,n-k+1): 
        for i in range(1,k+1): 
            bc[i] = bc[i-1]+bc[i] 
    return bc[k] #binomial


#***************************************
def binomialfilter(tosmooth,mdi,n,pad=False):
    '''

    data - data array
    mdi - missing data indicator
    n - filter size
    '''
    

    smoothed=np.zeros(len(tosmooth))
    smoothed.fill(mdi)

    weights=[]
    for k in range(n):
        weights+=[binomial(n-1,k)]

    weights=np.array(weights)

    for o,obs in enumerate(tosmooth):
        if (o >= n/2) and (o <= len(tosmooth)-n/2):
            
            chunk=tosmooth[o-(n/2):o+(n/2)+1]
            good=np.where(chunk != mdi)
            
            if len(good[0]) >= (n/2)+1:
                
                norm=sum(weights[good])
                weighted=sum(chunk[good]*weights[good]/norm)

                smoothed[o]=weighted

        elif pad==True:

            if (o < n/2):
                smoothed[o]=np.mean(tosmooth[:n/2])
            elif (o > len(tosmooth)-n/2):
                smoothed[o]=np.mean(tosmooth[-n/2:])
            

    return smoothed #binomialfilter

#************************************************************************
def GetYears(cube):

    times = cube.coord('time').points
    years = np.round(np.array([(t - 101)/10000 for t in times]))

    return years

#************************************************************************
def ApplyClimatology(cube):

    years = GetYears(cube)

    clim_years=np.where((years >= CLIMSTART) & (years <= CLIMEND))
    
    climatology = np.ma.mean(cube.data[clim_years], axis = 0)
    
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

    nyears = np.ma.count(cube.data, axis = 0)

    # repeat so selection will work across all years
    nyears = np.array([nyears for i in range(cube.data.shape[0])])

    cube.data = np.ma.masked_where(nyears < THRESHOLD * cube.data.shape[0], cube.data)

    weight_areas = iris.analysis.cartography.area_weights(cube)
    ts = cube.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=weight_areas)

    ts = utils.Timeseries(ts_name, GetYears(ts), ts.data * 3.65)

    return ts

#************************************************************************
def calc_rank(data):
    
    order = np.argsort(data)

#    return np.argsort(order) # calc_rank

    # derive plotting array

    # get lowest/highest 4 and mark with -4 --> -1, and 1 --> 4 for plotting
    N=5
    plot_order = np.argsort(order)
    if len(data) > 1:

        for rank in range(len(plot_order)):

            loc, = np.where(plot_order == rank)

            if rank in range(0, N+1):
                plot_order[loc] = rank-N+1
                
            elif rank in range(len(data)-N, len(data)):
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
    data = np.swapaxes(data, 1,2) # swap lats and lons around

    rank = np.ma.zeros(data.shape)
    rank.mask = data.mask
    
    for lt in range(len(lats)):
        for ln in range(len(lons)):
            clean_loc, = np.where(data.mask[:,ln,lt] == False)
            rank[clean_loc,ln,lt] = calc_rank(data[clean_loc,ln,lt])

    cube = utils.make_iris_cube_3d(rank, time, "unknown",lons, lats, "TEX rank", "1")

    return cube # get_ranks

#************************************************************************
def plot_rank_map(outname, cube, cmap, bounds, cb_label, scatter = [], figtext = "", title = ""):
    '''
    Standard scatter map

    :param str outname: output filename root
    :param array cube: cube to plot
    :param obj cmap: colourmap to use
    :param array bounds: bounds for discrete colormap
    :param str cb_label: colorbar label
    '''

    norm=mpl.cm.colors.BoundaryNorm(bounds,cmap.N)

    fig = plt.figure(figsize =(10,6.5))

    plt.clf()
    ax = plt.axes([0.05, 0.10, 0.90, 0.90], projection=cartopy.crs.Robinson())
    ax.gridlines() #draw_labels=True)
    ax.add_feature(cartopy.feature.LAND, zorder = 0, facecolor = "0.9", edgecolor = "k")
    ax.coastlines()

    ext = ax.get_extent() # save the original extent

    mesh = iris.plot.pcolormesh(cube, cmap=cmap, norm = norm)

    if len(scatter) > 0:

        lons, lats, data = scatter

        plt.scatter(lons, lats, c = data, cmap = cmap, norm = norm, s=25, \
                        transform = cartopy.crs.Geodetic(), edgecolor = '0.5', linewidth='0.5')


    cb=plt.colorbar(mesh, orientation = 'horizontal', pad = 0.05, fraction = 0.05, \
                        aspect = 30, ticks = bounds[1:-1], drawedges=True)

    # thicken border of colorbar and the dividers
    # http://stackoverflow.com/questions/14477696/customizing-colorbar-border-color-on-matplotlib
    cb.outline.set_color('k')
    cb.outline.set_linewidth(2)
    cb.dividers.set_color('k')
    cb.dividers.set_linewidth(2)
    cb.set_label(cb_label, labelpad=20)


    # label colorbar with sensibly placed and named labels
    cb.ax.get_xaxis().set_ticks([])
    cb.ax.get_xaxis().set_ticklabels([""])
    for j, lab in enumerate(['lowest','2nd lowest','3rd lowest','','3rd highest','2nd highest','highest']):
        cb.ax.text((2 * j + 1) / 14.0, -0.5, lab, ha='center', va='center')



    ax.set_extent(ext, ax.projection) # fix the extent change from colormesh

    plt.title(title)
    fig.text(0.03, 0.95, figtext, fontsize = settings.FONTSIZE * 0.8)

    plt.savefig(outname + settings.OUTFMT)
    plt.close()

    return # plot_rank_map_iris

#************************************************************************
def run_all_plots():

    for index in INDICES:

        # dummy so far
#        NYEARS = 66
#        rank_bounds = [-1,1,2,3,NYEARS-2,NYEARS-1,NYEARS,100]
        rank_bounds = [-4.5,-3.5,-2.5,-1.5,1.5,2.5,3.5,4.5]

        # sort the bounds and colourbars
        if index in ["TX90p", "TN90p"]:
#            bounds = [-100, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 100]
            bounds = [-100, -30, -20, -10, -5, 0, 5, 10, 20, 30, 100]
            cmap=settings.COLOURMAP_DICT["temperature"]
        elif index in ["TX10p", "TN10p"]:
#            bounds = [-100, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 100]
            bounds = [-100, -30, -20, -10, -5, 0, 5, 10, 20, 30, 100]
            cmap=settings.COLOURMAP_DICT["temperature_r"]
        elif index in ["TXx", "TNx", "TXn", "TNn"]:
            bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
            cmap=settings.COLOURMAP_DICT["temperature"]

        cube_list = iris.load(data_loc + "GHCND_{}_1951-{}_RegularGrid_global_2.5x2.5deg_LSmask.nc".format(index, int(settings.YEAR) + 1))
        names = np.array([cube.name() for cube in cube_list])

        #*************
        # plot annual map

        selected_cube, = np.where(names == "Ann")

        cube = cube_list[selected_cube]
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()  

        if index in ["TX90p", "TN90p", "TX10p", "TN10p"]:
            # change from % to days
            cube.data = cube.data * 3.65


        cube = ApplyClimatology(cube)

        # select the year to plot
        years = GetYears(cube)
        loc, = np.where(years == SELECTED_YEAR)

        utils.plot_smooth_map_iris(image_loc + "TEX_{}_{}_anoms_ghcndex".format(index, settings.YEAR), cube[loc[0]], cmap, bounds, "Anomalies from 1961-90 (Days)", title = "{} - {}".format(index, INDEX_LABELS[index]))

        if index == "TX90p":
            utils.plot_smooth_map_iris(image_loc + "p2.1_TEX_{}_{}_anoms_ghcndex".format(index, settings.YEAR), cube[loc[0]], cmap, bounds, "Anomalies from 1961-90 (Days)", title = "", figtext = "(c) Warm Days", save_netcdf_filename = "{}{}_for_NOAA_{}.nc".format(data_loc, index, dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y")))
        if index == "TX10p":
            utils.plot_smooth_map_iris(image_loc + "p2.1_TEX_{}_{}_anoms_ghcndex".format(index, settings.YEAR), cube[loc[0]], cmap, bounds, "Anomalies from 1961-90 (Days)", title = "", figtext = "(d) Cool Days", save_netcdf_filename = "{}{}_for_NOAA_{}.nc".format(data_loc, index, dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y")))


        rank_cube = get_ranks(cube)

        plot_rank_map(image_loc + "TEX_{}_{}_rank_ghcndex".format(index, settings.YEAR), rank_cube[loc[0]], cmap, rank_bounds, "Rank", title = "{} - {}".format(index, INDEX_LABELS[index]))


        #*************
        # plot season maps (2x2)

        season_list = []
        for season in SEASONS:

            # extract each month
            month_data = []
            months = SEASON_DICT[season]
            for month in months:

                selected_cube, = np.where(names == month)
                cube = cube_list[selected_cube]

                if month  == "Dec":
                    # need to extract from previous year - cheat by rolling data around
                    cube.data = np.roll(cube.data, 1, axis = 0)
                    cube.data.mask[0,:,:] = True # and mask out the previous years'

                if index in ["TX90p", "TN90p", "TX10p", "TN10p"]:
                    # change from % to days
                    cube.data = cube.data * (3.65/4.) # assume a season is 1/4 of a year

                month_data += [cube.data]

            # finished getting all months, make a dummy cube to populate
            month_data = np.ma.array(month_data)
            season_cube = copy.deepcopy(cube)

            # take appropriate seasonal value
            if index in ["TX90p", "TN90p", "TX10p", "TN10p"]:
                season_cube.data = np.ma.mean(month_data, axis = 0)
            elif index in ["TXx", "TNx"]:
                season_cube.data = np.ma.max(month_data, axis = 0)
            elif index in ["TXn", "TNn"]:
                season_cube.data = np.ma.min(month_data, axis = 0)

            # mask if fewer that 2 months present
            nmonths_locs = np.ma.count(month_data, axis = 0)
            season_cube.data = np.ma.masked_where(nmonths_locs < 2, season_cube.data)

            # make anomalies
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
        if index in ["TX90p", "TN90p"]:
            bounds = [-100, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 100]
        elif index in ["TX10p", "TN10p"]:
            bounds = [-100, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 100]
        elif index in ["TXx", "TNx", "TXn", "TNn"]:
            bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]

        # pass to plotting routine
        utils.plot_smooth_map_iris_multipanel(image_loc + "TEX_{}_{}_seasons_ghcndex".format(index, settings.YEAR), season_list, cmap, bounds, "Anomalies from 1961-90 (Days)", shape = (2,2), title = SEASONS, figtext = ["(a)","(b)","(c)","(d)"], figtitle = "{} - {}".format(index, INDEX_LABELS[index]))

    #*************
    # timeseries obs

    for limit in ["X", "N"]:

        fig, axes = plt.subplots(2, figsize = (10, 8), sharex=True)

        for ix, index in enumerate(["T{}90p".format(limit), "T{}10p".format(limit)]):

            index_ts = obtain_timeseries(data_loc + "GHCND_{}_1951-{}_RegularGrid_global_2.5x2.5deg_LSmask.nc".format(index, int(settings.YEAR) + 1), "Ann", "GHCNDEX")

            utils.plot_ts_panel(axes[ix], [index_ts], "-", "temperature", loc = LEGEND_LOC, bbox = BBOX)
            axes[ix].text(0.02, 0.9, "({}) {}".format(string.lowercase[ix], index), transform = axes[ix].transAxes, fontsize = settings.FONTSIZE)

            # and smoothed
            index_ts.data.fill_value = -99.9
            smoothed = binomialfilter(index_ts.data.filled(), -99.9, 5, pad = False)
            smoothed = np.ma.masked_where(smoothed == -99.9, smoothed)

            axes[ix].plot(index_ts.times, smoothed, "r--", lw = LW)
            
            # print the index, the current anomaly and the rank information
            print index, index_ts.data.compressed()[-1]/36.5, np.argsort(np.argsort(index_ts.data.compressed())) + 1

        # prettify
        plt.xlim([1950,2019])
        axes[0].set_ylim([15,None])
        if limit == "X":
            axes[1].set_ylim([19,59])
        elif limit == "N":
            axes[1].set_ylim([11,59])


        fig.text(0.03, 0.5, "Number of Days", va='center', rotation='vertical', fontsize = settings.FONTSIZE)

        for ax in axes:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)
        for tick in axes[1].xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)

        fig.subplots_adjust(right = 0.95, top = 0.95, bottom = 0.05, hspace = 0.001)

        plt.savefig(image_loc+"TEX_T{}_ts_ghcndex{}".format(limit,settings.OUTFMT))
        plt.close()



    #*************
    # timeseries ERA


    for limit in ["X", "N"]:

        fig, axes = plt.subplots(2, figsize = (10, 8), sharex=True)

        for ix, index in enumerate(["T{}90p".format(limit), "T{}10p".format(limit)]):

            index_ts = obtain_timeseries(data_loc + "ERA-Int_1979-{}_{}_LSmask.nc".format(settings.YEAR, index), "Annual", "ERA-Interim")

            utils.plot_ts_panel(axes[ix], [index_ts], "-", "temperature", loc = LEGEND_LOC, bbox = BBOX)

            axes[ix].text(0.02, 0.9, "({}) {}".format(string.lowercase[ix], index), transform = axes[ix].transAxes, fontsize = settings.FONTSIZE)

            # and smoothed
            index_ts.data.fill_value = -99.9
            smoothed = binomialfilter(index_ts.data.filled(), -99.9, 5, pad = False)
            smoothed = np.ma.masked_where(smoothed == -99.9, smoothed)

            axes[ix].plot(index_ts.times, smoothed, "r--", lw = LW)

        # prettify
        plt.xlim([1950,2019])
        axes[0].set_ylim([15,None])
        axes[1].set_ylim([16,59])

        fig.text(0.03, 0.5, "Number of Days", va='center', rotation='vertical', fontsize = settings.FONTSIZE)

        for ax in axes:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)
        for tick in axes[1].xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)

        fig.subplots_adjust(right = 0.95, top = 0.95, bottom = 0.05, hspace = 0.001)

        plt.savefig(image_loc+"TEX_T{}_ts_era{}".format(limit,settings.OUTFMT))
        plt.close()


    #*************
    # ERA Maps (annual only)

    for index in INDICES:

        # sort the bounds and colourbars
        if index in ["TX90p", "TN90p"]:
#           bounds = [-100, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 100]
            bounds = [-100, -30, -20, -10, -5, 0, 5, 10, 20, 30, 100]
            cmap=settings.COLOURMAP_DICT["temperature"]
        elif index in ["TX10p", "TN10p"]:
#            bounds = [-100, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 100]
            bounds = [-100, -30, -20, -10, -5, 0, 5, 10, 20, 30, 100]
            cmap=settings.COLOURMAP_DICT["temperature_r"]
        elif index in ["TXx", "TNx", "TXn", "TNn"]:
            bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
            cmap=settings.COLOURMAP_DICT["temperature"]

        cube_list = iris.load(data_loc + "ERA-Int_1979-{}_{}_LSmask.nc".format(settings.YEAR, index))
        names = np.array([cube.name() for cube in cube_list])

        #*************
        # plot annual map

        selected_cube, = np.where(names == "Annual")

        cube = cube_list[selected_cube]
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()  

        if index in ["TX90p", "TN90p", "TX10p", "TN10p"]:
            # change from % to days
            cube.data = cube.data * 3.65

        cube = ApplyClimatology(cube)

        # select the year to plot
        years = GetYears(cube)
        loc, = np.where(years == SELECTED_YEAR)

        utils.plot_smooth_map_iris(image_loc + "TEX_{}_{}_anoms_era".format(index, settings.YEAR), cube[loc[0]], cmap, bounds, "Anomalies from 1961-90 (Days)", title = "ERA-Interim {} - {}".format(index, INDEX_LABELS[index]))

        #*************
        # plot season maps (2x2)

        season_list = []
        for season in SEASONS:

            # extract each month
            month_data = []
            months = SEASON_DICT_ERA[season]
            for month in months:

                selected_cube, = np.where(names == month)
                cube = cube_list[selected_cube]

                if month  == "December":
                    # need to extract from previous year - cheat by rolling data around
                    cube.data = np.roll(cube.data, 1, axis = 0)
                    cube.data.mask[0,:,:] = True # and mask out the previous years'

                if index in ["TX90p", "TN90p", "TX10p", "TN10p"]:
                    # change from % to days
                    cube.data = cube.data * (3.65/4.) # assume a season is 1/4 of a year

                month_data += [cube.data]

            # finished getting all months, make a dummy cube to populate
            month_data = np.ma.array(month_data)
            season_cube = copy.deepcopy(cube)

            # take appropriate seasonal value
            if index in ["TX90p", "TN90p", "TX10p", "TN10p"]:
                season_cube.data = np.ma.mean(month_data, axis = 0)
            elif index in ["TXx", "TNx"]:
                season_cube.data = np.ma.max(month_data, axis = 0)
            elif index in ["TXn", "TNn"]:
                season_cube.data = np.ma.min(month_data, axis = 0)

            # mask if fewer that 2 months present
            nmonths_locs = np.ma.count(month_data, axis = 0)
            season_cube.data = np.ma.masked_where(nmonths_locs < 2, season_cube.data)

            # make anomalies
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
        if index in ["TX90p", "TN90p"]:
            bounds = [-100, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 100]
        elif index in ["TX10p", "TN10p"]:
            bounds = [-100, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 100]
        elif index in ["TXx", "TNx", "TXn", "TNn"]:
            bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]

        # pass to plotting routine
        utils.plot_smooth_map_iris_multipanel(image_loc + "TEX_{}_{}_seasons_era".format(index, settings.YEAR), season_list, cmap, bounds, "Anomalies from 1961-90 (Days)", shape = (2,2), title = SEASONS, figtext = ["(a)","(b)","(c)","(d)"], figtitle = "{} - {}".format(index, INDEX_LABELS[index]))

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************