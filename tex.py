#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for T extremes (TEX) section.
#       For BAMS SotC 2016 
#
#************************************************************************
#                    SVN Info
# $Rev:: 31                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2021-09-06 09:52:46 +0100 (Mon, 06 Sep #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.cm as mpl_cm
import matplotlib as mpl

import copy
import string

import iris
import iris.coord_categorisation
import iris.quickplot as qplt
import cartopy
import cf_units

import datetime as dt
import calendar

import utils # RJHD utilities
import settings

DATALOC = "{}/{}/data/TEX/".format(settings.ROOTLOC, settings.YEAR)

DATASTART=1951
CLIMSTART=1961
CLIMEND=1990
SELECTED_YEAR = int(settings.YEAR)
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

INDEX_LABELS = {"TX90p" : "Warm Days", "TX10p" : "Cool Days", "TN90p" : "Warm Nights", "TN10p" : "Cool Nights", "TXx" : "max Tmax", "TXn" : "min Tmax", "TNx" : "max Tmin", "TNn": "min Tmin"}

INDEX_UNITS = {"TX90p" : "Days", "TX10p" : "Days", "TN90p" : "Days", "TN10p" : "Days", "TXx" : '$^{\circ}$'+"C", "TXn" : '$^{\circ}$'+"C", "TNx" : '$^{\circ}$'+"C", "TNn": '$^{\circ}$'+"C"}

FIGURE_LABELS = {"TX90p" : "(a)", "TX10p" : "(c)", "TN90p" : "(d)", "TN10p" : "(b)", "TXx" : "(a)", "TXn" : "(c)", "TNx" : "(b)", "TNn": "(d)"}

SEASON_LABELS = {"TX90p" : ["(a)","(b)","(c)","(d)"], "TX10p" : ["(e)","(f)","(g)","(h)"], "TN90p" : ["(i)","(j)","(h)","(l)"], "TN10p" : ["(m)","(n)","(o)","(p)"], "TXx" : ["(a)","(b)","(c)","(d)"], "TXn" : ["(e)","(f)","(g)","(h)"], "TNx" : ["(a)","(b)","(c)","(d)"] , "TNn" : ["(e)","(f)","(g)","(h)"]}

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
        weights+=[binomial(n-1, k)]

    weights=np.array(weights)

    # this needs testing under python3 given apparent integer division
    for o,obs in enumerate(tosmooth):
        if (o >= n/2) and (o <= len(tosmooth)-n//2):
            
            chunk=tosmooth[o-(n//2):o+(n//2)+1]
            good=np.where(chunk != mdi)
            
            if len(good[0]) >= (n//2)+1:
                
                norm=sum(weights[good])
                weighted=sum(chunk[good]*weights[good]/norm)

                smoothed[o]=weighted

        elif pad==True:

            if (o < n//2):
                smoothed[o]=np.mean(tosmooth[:n//2])
            elif (o > len(tosmooth)-n/2):
                smoothed[o]=np.mean(tosmooth[-n//2:])
            

    return smoothed #binomialfilter

# #************************************************************************
# def GetYears(cube, is_era5 = False):

#     if is_era5:

#         timeUnits = cube.coord("time").units
#         dt_time = timeUnits.num2date(cube.coord("time").points)

#         years = np.array([d.year for d in dt_time])

#     else:
#         times = cube.coord('time').points
#         years = np.round(np.array([(t - 101)/10000 for t in times]))

#     return years

#************************************************************************
def GetYears(cube, is_era5 = False):

    if is_era5:

        timeUnits = cube.coord("time").units
        dt_time = timeUnits.num2date(cube.coord("time").points)

        years = np.array([d.year for d in dt_time])

    else:
        times = cube.coord('time').points
        if times[0] < 19500000:
            years = times + DATASTART
        else:
            years = np.round(np.array([(t - 101)/10000 for t in times]))
    return years

#************************************************************************
def ApplyClimatology(cube, is_era5=False):

    years = GetYears(cube, is_era5=is_era5)

    clim_years=np.where((years >= CLIMSTART) & (years <= CLIMEND))
    
    climatology = np.ma.mean(cube.data[clim_years], axis = 0)
    
    cube.data = cube.data - climatology                                
        
    return cube

#************************************************************************
def obtain_timeseries(filename, cube_name, ts_name, index, is_era5 = False):

    # re-read the file to overwrite memory of what was done to this cube.
    cube_list = iris.load(filename)
    names = np.array([cube.name() for cube in cube_list])

    selected_cube, = np.where(names == cube_name)
    if len(selected_cube) == 0:
        names = np.array([cube.var_name for cube in cube_list])
        selected_cube, = np.where(names == cube_name)
        

    cube = cube_list[selected_cube[0]]
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()  

    nyears = np.ma.count(cube.data, axis = 0)

    # repeat so selection will work across all years
    nyears = np.array([nyears for i in range(cube.data.shape[0])])

    cube.data = np.ma.masked_where(nyears < THRESHOLD * cube.data.shape[0], cube.data)

    weight_areas = iris.analysis.cartography.area_weights(cube)
    ts = cube.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=weight_areas)

    times = GetYears(ts, is_era5=is_era5)

    if index in ["TX90p", "TN90p", "TX10p", "TN10p"]:
        ts = utils.Timeseries(ts_name, times, ts.data * 3.65)

    else:
        ts = utils.Timeseries(ts_name, times, ts.data)

    # get coverage
    cover_cube = copy.deepcopy(cube)
    # replace data with 1
    cover_cube.data.data[:] = 1
    # and get timeseries
    weight_areas = iris.analysis.cartography.area_weights(cover_cube)
    cover_ts = cover_cube.collapsed(['latitude', 'longitude'], iris.analysis.SUM, weights=weight_areas)

    # normalise - obtain max area possible with each gridbox = 1
    norm_cube = copy.deepcopy(cover_cube)[0]
    norm_cube.data.mask = False
    weight_areas = iris.analysis.cartography.area_weights(norm_cube)
    norm = norm_cube.collapsed(['latitude', 'longitude'], iris.analysis.SUM, weights=weight_areas)
   
    # apply normalisation and convert to %.  Also convert to %land.
    cover_ts = (cover_ts/norm)*100/0.3
    cover_ts = utils.Timeseries("Coverage", GetYears(cover_ts), cover_ts.data)
        
    return ts, cover_ts # obtain_timeseries

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

    # if have following year already present
    while time[-1] > (int(settings.YEAR)*10000+101):
        time = time[: -1]
        data = data[: -1, :, :]

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

    fig = plt.figure(figsize=(8, 5.5))

    plt.clf()
    ax = plt.axes([0.01, 0.12, 0.98, 0.88], projection=cartopy.crs.Robinson())
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
#    cb.outline.set_color('k')
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

#***************************************
def fix_time_coord(incube):

    try:
        newtimes = np.array([dt.datetime.strptime("{}".format(int(d)), "%Y%m%d") for d in incube.coord("time").points])
    except ValueError:
        # 2019 run
        try:
            newtimes = np.array([dt.datetime.strptime("{}".format(int(d)), "%Y") for d in incube.coord("time").points])
        except ValueError:
            # 2020 run
            newtimes = np.array([dt.datetime.strptime("{}".format(int(d+1951)), "%Y") for d in incube.coord("time").points])
        
        
    newdiff = np.array([(t - newtimes[0]).days for t in newtimes])
    
    # replace cube time coordinate with new one
    time_unit = cf_units.Unit('days since ' + dt.datetime.strftime(newtimes[0], "%Y-%m-%d %H:%M"), calendar=cf_units.CALENDAR_GREGORIAN)   
    timecoord = iris.coords.DimCoord(newdiff, standard_name='time', units=time_unit, var_name="time") # add bounds?
    incube.remove_coord('time')
    incube.add_dim_coord(timecoord, 0)

    return incube # fix_time_coord

#*******************************
def compute_coverage_error(observations, reanalysis):
    '''
    Calculate the coverage error on a monthly basis

    Takes each month in observations, find all the corresponding calendar months in reanalysis
    (i.e. for Jan 1973, finds all Jans in ERA).  Then mask by data coverage, and get
    range of residuals in global average.  Use these residuals to estimate error

    From HadEX3 code

    '''
    
    offset = np.zeros(len(observations.coord("time").points))
    st_dev = np.zeros(len(observations.coord("time").points))

    # add month names into data
    iris.coord_categorisation.add_month(reanalysis, 'time', name='month')
    iris.coord_categorisation.add_month(observations, 'time', name='month')

    for m, month  in enumerate(observations.coord("time").points):

        # get weightings for cubes
        grid_areas = iris.analysis.cartography.cosine_latitude_weights(observations[m]) 
        area_average = observations[m].collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)

        if area_average == observations.data.fill_value:
            # no data for this month
            offset[m] = observations.data.fill_value
            st_dev[m] = observations.data.fill_value
            
        else:
            # make a copy so can update the mask without worrying
            reanal = copy.deepcopy(reanalysis)

            # get a clean-mean
            grid_areas = iris.analysis.cartography.cosine_latitude_weights(reanal)              

            total_mean = reanal.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)

            # apply the observation coverage mask to all months
            # combine the lsm with the data mask
            combined_mask = np.ma.mask_or(np.array([observations[m].data.mask for i in range(reanal.data.shape[0])]), reanal.data.mask)

            reanal.data.mask = combined_mask

            # get a masked mean
            grid_areas = iris.analysis.cartography.cosine_latitude_weights(reanal)              
            masked_mean = reanal.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)
            
            # calculate residuals and find mean and st_dev
            residuals = masked_mean.data - total_mean.data

            offset[m] = np.mean(residuals)
            st_dev[m] = np.std(residuals, ddof=1) # to match IDL

    return offset, st_dev # compute_coverage_error

#************************************************************************
def run_all_plots():

    if True:
        #*************
        # plot annual map

        for index in INDICES:

            # sort the bounds and colourbars
            if index in ["TX90p", "TN90p"]:
                bounds = [-100, -40, -30, -20, -10, 0, 10, 20, 30, 40, 100]
                cmap=settings.COLOURMAP_DICT["temperature"]
            elif index in ["TX10p", "TN10p"]:
                bounds = [-100, -40, -30, -20, -10, 0, 10, 20, 30, 40, 100]
                cmap=settings.COLOURMAP_DICT["temperature_r"]
            elif index in ["TXx", "TNx", "TXn", "TNn"]:
                bounds = [-100, -6, -4, -2, -1, 0, 1, 2, 4, 6, 100]
                cmap=settings.COLOURMAP_DICT["temperature"]

            cube_list = iris.load(DATALOC + "GHCND_{}_1951-{}_RegularGrid_global_2.5x2.5deg_LSmask.nc".format(index, int(settings.YEAR) + 1))
            names = np.array([cube.name() for cube in cube_list])

            selected_cube, = np.where(names == "Ann")[0]

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

            utils.plot_smooth_map_iris(settings.IMAGELOC + "TEX_{}_{}_anoms_ghcndex".format(index, settings.YEAR), cube[loc[0]], cmap, bounds, "Anomalies from 1961-90 ({})".format(INDEX_UNITS[index]), title = "{} - {}".format(index, INDEX_LABELS[index]), figtext = FIGURE_LABELS[index])

            if index == "TX90p":
                utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_TEX_{}_{}_anoms_ghcndex".format(index, settings.YEAR), cube[loc[0]], cmap, bounds, "Anomalies from 1961-90 ({})".format(INDEX_UNITS[index]), title = "", figtext = "(c) Warm Days", save_netcdf_filename = "{}{}_for_NOAA_{}.nc".format(DATALOC, index, dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y")))
            if index == "TN10p":
                utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_TEX_{}_{}_anoms_ghcndex".format(index, settings.YEAR), cube[loc[0]], cmap, bounds, "Anomalies from 1961-90 ({})".format(INDEX_UNITS[index]), title = "", figtext = "(d) Cool Nights", save_netcdf_filename = "{}{}_for_NOAA_{}.nc".format(DATALOC, index, dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y")))

    if False:                
        #*************
        # plot rank annual map
        for index in INDICES:

            cube_list = iris.load(DATALOC + "GHCND_{}_1951-{}_RegularGrid_global_2.5x2.5deg_LSmask.nc".format(index, int(settings.YEAR) + 1))
            names = np.array([cube.name() for cube in cube_list])

            selected_cube, = np.where(names == "Ann")[0]

            cube = cube_list[selected_cube]
            cube.coord('latitude').guess_bounds()
            cube.coord('longitude').guess_bounds()  

            if index in ["TX90p", "TN90p", "TX10p", "TN10p"]:
                # change from % to days
                cube.data = cube.data * 3.65


            cube = ApplyClimatology(cube)

            rank_bounds = [-4.5,-3.5,-2.5,-1.5,1.5,2.5,3.5,4.5]

            rank_cube = get_ranks(cube)

            # select the year to plot
            years = GetYears(cube)
            loc, = np.where(years == SELECTED_YEAR)

            plot_rank_map(settings.IMAGELOC + "TEX_{}_{}_rank_ghcndex".format(index, settings.YEAR), rank_cube[loc[0]], cmap, rank_bounds, "Rank", title = "{} - {}".format(index, INDEX_LABELS[index]))


    if True:                
        #*************
        # plot seasonal map 2x2
        for index in INDICES:

            cube_list = iris.load(DATALOC + "GHCND_{}_1951-{}_RegularGrid_global_2.5x2.5deg_LSmask.nc".format(index, int(settings.YEAR) + 1))
            names = np.array([cube.name() for cube in cube_list])

            selected_cube, = np.where(names == "Ann")[0]

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

            season_list = []
            for season in SEASONS:

                # extract each month
                month_data = []
                months = SEASON_DICT[season]
                for month in months:

                    selected_cube, = np.where(names == month)[0]
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
                bounds = [-100, -6, -4, -2, -1, 0, 1, 2, 4, 6, 100]

            # pass to plotting routine
            utils.plot_smooth_map_iris_multipanel(settings.IMAGELOC + "TEX_{}_{}_seasons_ghcndex".format(index, settings.YEAR), season_list, cmap, bounds, "Anomalies from 1961-90 ({})".format(INDEX_UNITS[index]), shape = (2,2), title = SEASONS, figtext = SEASON_LABELS[index], figtitle = "{} - {}".format(index, INDEX_LABELS[index]))

    #*************
    # timeseries obs as a pair
    if True:
        for index_pair in [["TX90p", "TN10p"], ["TN90p", "TX10p"]]:

            fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 6.5), sharex=True)

            ax3 = ax1.twinx()
            ax4 = ax2.twinx()

            axes = (ax1, ax2, ax3, ax4)

            for ix, index in enumerate(index_pair):

                index_ts, cover_ts = obtain_timeseries(DATALOC + "GHCND_{}_1951-{}_RegularGrid_global_2.5x2.5deg_LSmask.nc".format(index, int(settings.YEAR) + 1), "Ann", "GHCNDEX", index)

                utils.plot_ts_panel(axes[ix], [index_ts], "-", "temperature", loc = "", bbox = BBOX) # no legend as single line
                axes[ix].text(0.02, 0.9, "({}) {}".format(string.ascii_lowercase[ix], index), transform = axes[ix].transAxes, fontsize = settings.FONTSIZE)

                # red tickmarks
                axes[ix].tick_params(axis='y', colors='red', direction="in")
                axes[ix].yaxis.tick_left() # none on right

                # and smoothed
                index_ts.data.fill_value = -99.9
                smoothed = binomialfilter(index_ts.data.filled(), -99.9, 5, pad = False)
                smoothed = np.ma.masked_where(smoothed == -99.9, smoothed)

                axes[ix].plot(index_ts.times, smoothed, "r--", lw = LW)

                # print the index, the current anomaly and the rank information
                # print(index, index_ts.data.compressed()[-1]/36.5, np.argsort(np.argsort(index_ts.data.compressed())) + 1)

                axes[ix+2].plot(cover_ts.times, cover_ts.data, "k:", lw = 2)
                axes[ix+2].yaxis.set_label_position("right")
                axes[ix+2].yaxis.set_ticks_position('right')


            # prettify
            plt.xlim([1950-2, int(settings.YEAR)+2])
            if "X" in index:
                axes[0].set_ylim([21,79])
                axes[1].set_ylim([11,59])
            elif "N" in index:
                axes[0].set_ylim([21,79])
                axes[1].set_ylim([11,59])

            ax3.set_ylim([0,98])
            ax4.set_ylim([0,98])

            ax3.tick_params(axis='y', direction="in", width=2, length=10)
            ax4.tick_params(axis='y', direction="in", width=2, length=10)


            fig.text(0.03, 0.5, "Number of Days", va='center', rotation='vertical', fontsize = settings.FONTSIZE, color="r")
            fig.text(0.97, 0.5, "% land covered", va='center', rotation='vertical', fontsize = settings.FONTSIZE)

            for ax in axes:
                for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(settings.FONTSIZE)
                    tick.label2.set_fontsize(settings.FONTSIZE)
            for tick in axes[1].xaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)

            fig.subplots_adjust(left=0.1, right = 0.9, top = 0.98, bottom = 0.05, hspace = 0.001)

            plt.savefig(settings.IMAGELOC+"TEX_{}+{}_ts_ghcndex{}".format(index_pair[0], index_pair[1], settings.OUTFMT))
            plt.close()



    # #*************
    # # timeseries ERA-Int


    # for index_pair in [["TX90p", "TN10p"], ["TN90p", "TX10p"]]:

    #     fig, axes = plt.subplots(2, figsize=(8, 6.5), sharex=True)

    #     for ix, index in enumerate(index_pair):

    #         index_ts, cover_ts = obtain_timeseries(DATALOC + "ERA-Int_1979-{}_{}_LSmask.nc".format(settings.YEAR, index), "Annual", "ERA-Interim", index)

    #         utils.plot_ts_panel(axes[ix], [index_ts], "-", "temperature", loc = LEGEND_LOC, bbox = BBOX)

    #         axes[ix].text(0.02, 0.9, "({}) {}".format(string.ascii_lowercase[ix], index), transform = axes[ix].transAxes, fontsize = settings.FONTSIZE)

    #         # and smoothed
    #         index_ts.data.fill_value = -99.9
    #         smoothed = binomialfilter(index_ts.data.filled(), -99.9, 5, pad = False)
    #         smoothed = np.ma.masked_where(smoothed == -99.9, smoothed)

    #         axes[ix].plot(index_ts.times, smoothed, "r--", lw = LW)

    #     # prettify
    #     plt.xlim([1950,2019])
    #     axes[0].set_ylim([15,None])
    #     axes[1].set_ylim([16,59])

    #     fig.text(0.03, 0.5, "Number of Days", va='center', rotation='vertical', fontsize = settings.FONTSIZE)

    #     for ax in axes:
    #         for tick in ax.yaxis.get_major_ticks():
    #             tick.label.set_fontsize(settings.FONTSIZE)
    #     for tick in axes[1].xaxis.get_major_ticks():
    #         tick.label.set_fontsize(settings.FONTSIZE)

    #     fig.subplots_adjust(right = 0.95, top = 0.95, bottom = 0.05, hspace = 0.001)

    #     plt.savefig(settings.IMAGELOC+"TEX_{}+{}_ts_erai{}".format(index_pair[0], index_pair[1], settings.OUTFMT))
    #     plt.close()


    # #*************
    # # ERA-Int Maps (annual only)

    # for index in INDICES:

    #     # sort the bounds and colourbars
    #     if index in ["TX90p", "TN90p"]:
    #         bounds = [-100, -30, -20, -10, -5, 0, 5, 10, 20, 30, 100]
    #         cmap=settings.COLOURMAP_DICT["temperature"]
    #     elif index in ["TX10p", "TN10p"]:
    #         bounds = [-100, -30, -20, -10, -5, 0, 5, 10, 20, 30, 100]
    #         cmap=settings.COLOURMAP_DICT["temperature_r"]
    #     elif index in ["TXx", "TNx", "TXn", "TNn"]:
    #         bounds = [-100, -6, -4, -2, -1, 0, 1, 2, 4, 6, 100]
    #         cmap=settings.COLOURMAP_DICT["temperature"]

    #     cube_list = iris.load(DATALOC + "ERA-Int_1979-{}_{}_LSmask.nc".format(settings.YEAR, index))
    #     names = np.array([cube.name() for cube in cube_list])

    #     #*************
    #     # plot annual map

    #     selected_cube, = np.where(names == "Annual")[0]

    #     cube = cube_list[selected_cube]
    #     cube.coord('latitude').guess_bounds()
    #     cube.coord('longitude').guess_bounds()  

    #     if index in ["TX90p", "TN90p", "TX10p", "TN10p"]:
    #         # change from % to days
    #         cube.data = cube.data * 3.65

    #     cube = ApplyClimatology(cube)

    #     # select the year to plot
    #     years = GetYears(cube)
    #     loc, = np.where(years == SELECTED_YEAR)

    #     utils.plot_smooth_map_iris(settings.IMAGELOC + "TEX_{}_{}_anoms_erai".format(index, settings.YEAR), cube[loc[0]], cmap, bounds, "Anomalies from 1981-2010 ({})".format(INDEX_UNITS[index]), title = "ERA-Interim {} - {}".format(index, INDEX_LABELS[index]), figtext = FIGURE_LABELS[index])

    #     #*************
    #     # plot season maps (2x2)

    #     season_list = []
    #     for season in SEASONS:

    #         # extract each month
    #         month_data = []
    #         months = SEASON_DICT_ERA[season]
    #         for month in months:

    #             selected_cube, = np.where(names == month)[0]
    #             cube = cube_list[selected_cube]

    #             if month  == "December":
    #                 # need to extract from previous year - cheat by rolling data around
    #                 cube.data = np.roll(cube.data, 1, axis = 0)
    #                 cube.data.mask[0,:,:] = True # and mask out the previous years'

    #             if index in ["TX90p", "TN90p", "TX10p", "TN10p"]:
    #                 # change from % to days
    #                 cube.data = cube.data * (3.65/4.) # assume a season is 1/4 of a year

    #             month_data += [cube.data]

    #         # finished getting all months, make a dummy cube to populate
    #         month_data = np.ma.array(month_data)
    #         season_cube = copy.deepcopy(cube)

    #         # take appropriate seasonal value
    #         if index in ["TX90p", "TN90p", "TX10p", "TN10p"]:
    #             season_cube.data = np.ma.mean(month_data, axis = 0)
    #         elif index in ["TXx", "TNx"]:
    #             season_cube.data = np.ma.max(month_data, axis = 0)
    #         elif index in ["TXn", "TNn"]:
    #             season_cube.data = np.ma.min(month_data, axis = 0)

    #         # mask if fewer that 2 months present
    #         nmonths_locs = np.ma.count(month_data, axis = 0)
    #         season_cube.data = np.ma.masked_where(nmonths_locs < 2, season_cube.data)

    #         # make anomalies
    #         season_cube = ApplyClimatology(season_cube)

    #         # fix for plotting
    #         season_cube.coord('latitude').guess_bounds()
    #         season_cube.coord('longitude').guess_bounds()

    #         # select the year to plot
    #         years = GetYears(cube)
    #         loc, = np.where(years == SELECTED_YEAR)

    #         # add to list
    #         season_list += [season_cube[loc[0]]]

    #     # sort the bounds and colourbars
    #     if index in ["TX90p", "TN90p"]:
    #         bounds = [-100, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 100]
    #     elif index in ["TX10p", "TN10p"]:
    #         bounds = [-100, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 100]
    #     elif index in ["TXx", "TNx", "TXn", "TNn"]:
    #         bounds = [-100, -6, -4, -2, -1, 0, 1, 2, 4, 6, 100]

    #     # pass to plotting routine
    #     utils.plot_smooth_map_iris_multipanel(settings.IMAGELOC + "TEX_{}_{}_seasons_erai".format(index, settings.YEAR), season_list, cmap, bounds, "Anomalies from 1981-2010 ({})".format(INDEX_UNITS[index]), shape = (2,2), title = SEASONS, figtext = SEASON_LABELS[index], figtitle = "{} - {}".format(index, INDEX_LABELS[index]))

    #*************
    # timeseries ERA5 as a pair
    if True:
 
        for index_pair in [["TX90p", "TN10p"], ["TN90p", "TX10p"]]:

            fig, axes = plt.subplots(2, figsize=(8, 6.5), sharex=True)

            for ix, index in enumerate(index_pair):

                index_ts, cover_ts = obtain_timeseries(settings.ERA5LOC_TEX + "ERA5_{}_1979-{}_land.nc".format(index, settings.YEAR), "Ann", "ERA5", index, is_era5 = True)

                utils.plot_ts_panel(axes[ix], [index_ts], "-", "temperature", loc = LEGEND_LOC, bbox = BBOX)

                axes[ix].text(0.02, 0.9, "({}) {}".format(string.ascii_lowercase[ix+2], index), transform = axes[ix].transAxes, fontsize = settings.FONTSIZE)

                # and smoothed
                index_ts.data.fill_value = -99.9
                smoothed = binomialfilter(index_ts.data.filled(), -99.9, 5, pad = False)
                smoothed = np.ma.masked_where(smoothed == -99.9, smoothed)

                axes[ix].plot(index_ts.times, smoothed, "m--", lw = LW)

            # prettify
            plt.xlim([1950-2, int(settings.YEAR)+2])
            axes[0].set_ylim([15,None])
            axes[1].set_ylim([16,59])

            fig.text(0.03, 0.5, "Number of Days", va='center', rotation='vertical', fontsize = settings.FONTSIZE)

            for ax in axes:
                for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(settings.FONTSIZE)
            for tick in axes[1].xaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)

            fig.subplots_adjust(right = 0.95, top = 0.98, bottom = 0.05, hspace = 0.001)

            plt.savefig(settings.IMAGELOC+"TEX_{}+{}_ts_era5{}".format(index_pair[0], index_pair[1], settings.OUTFMT))
            plt.close()


    #*************
    # ERA Maps
    if True:
        for index in INDICES:

            # sort the bounds and colourbars
            if index in ["TX90p", "TN90p"]:
#                bounds = [-100, -30, -20, -10, -5, 0, 5, 10, 20, 30, 100]
                bounds = [-100, -40, -30, -20, -10, 0, 10, 20, 30, 40, 100]
                cmap=settings.COLOURMAP_DICT["temperature"]
            elif index in ["TX10p", "TN10p"]:
#                bounds = [-100, -30, -20, -10, -5, 0, 5, 10, 20, 30, 100]
                bounds = [-100, -40, -30, -20, -10, 0, 10, 20, 30, 40, 100]
                cmap=settings.COLOURMAP_DICT["temperature_r"]
            elif index in ["TXx", "TNx", "TXn", "TNn"]:
                bounds = [-100, -6, -4, -2, -1, 0, 1, 2, 4, 6, 100]
                cmap=settings.COLOURMAP_DICT["temperature"]

            cube_list = iris.load(settings.ERA5LOC_TEX + "ERA5_{}_1979-{}_land.nc".format(index, settings.YEAR))
            names = np.array([cube.var_name for cube in cube_list])

            #*************
            # plot annual map

            selected_cube, = np.where(names == "Ann")[0]

            cube = cube_list[selected_cube]
            cube.coord('latitude').guess_bounds()
            cube.coord('longitude').guess_bounds()  

            if index in ["TX90p", "TN90p", "TX10p", "TN10p"]:
                # change from % to days
                cube.data = cube.data * 3.65

            cube = ApplyClimatology(cube, is_era5=True)

            # select the year to plot
            years = GetYears(cube, is_era5 = True)
            loc, = np.where(years == SELECTED_YEAR)

            utils.plot_smooth_map_iris(settings.IMAGELOC + "TEX_{}_{}_anoms_era5".format(index, settings.YEAR), cube[loc[0]], cmap, bounds, "Anomalies from 1981-2010 ({})".format(INDEX_UNITS[index]), title = "ERA5 {} - {}".format(index, INDEX_LABELS[index]), figtext = FIGURE_LABELS[index])


            #*************
            # plot season maps (2x2)

            season_list = []
            for season in SEASONS:

                # extract each month
                month_data = []
                months = SEASON_DICT[season]
                for month in months:

                    selected_cube, = np.where(names == month)[0]
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
                season_cube = ApplyClimatology(season_cube, is_era5=True)

                # fix for plotting
                season_cube.coord('latitude').guess_bounds()
                season_cube.coord('longitude').guess_bounds()

                # select the year to plot
                years = GetYears(cube, is_era5 = True)
                loc, = np.where(years == SELECTED_YEAR)

                # add to list
                season_list += [season_cube[loc[0]]]

            # sort the bounds and colourbars
            if index in ["TX90p", "TN90p"]:
                bounds = [-100, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 100]
            elif index in ["TX10p", "TN10p"]:
                bounds = [-100, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 100]
            elif index in ["TXx", "TNx", "TXn", "TNn"]:
                bounds = [-100, -6, -4, -2, -1, 0, 1, 2, 4, 6, 100]

            # pass to plotting routine
            utils.plot_smooth_map_iris_multipanel(settings.IMAGELOC + "TEX_{}_{}_seasons_era5".format(index, settings.YEAR), season_list, cmap, bounds, "Anomalies from 1981-2010 ({})".format(INDEX_UNITS[index]), shape = (2,2), title = SEASONS, figtext = SEASON_LABELS[index], figtitle = "ERA5 {} - {}".format(index, INDEX_LABELS[index]))

    #*************
    # timeseries obs with uncertainties
    if True:
        for index_pair in [["TX90p", "TN10p"], ["TN90p", "TX10p"]]:

            fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 6.5), sharex=True)

            ax3 = ax1.twinx()
            ax4 = ax2.twinx()

            axes = (ax1, ax2, ax3, ax4)

            for ix, index in enumerate(index_pair):

                index_ts, cover_ts = obtain_timeseries(DATALOC + "GHCND_{}_1951-{}_RegularGrid_global_2.5x2.5deg_LSmask.nc".format(index, int(settings.YEAR) + 1), "Ann", "GHCNDEX", index)

                utils.plot_ts_panel(axes[ix], [index_ts], "-", "temperature", loc = "", bbox = BBOX) # no legend as single line
                axes[ix].text(0.02, 0.9, "({}) {}".format(string.ascii_lowercase[ix], index), transform = axes[ix].transAxes, fontsize = settings.FONTSIZE)

                # obs cube
                cube_list = iris.load(DATALOC + "GHCND_{}_1951-{}_RegularGrid_global_2.5x2.5deg_LSmask.nc".format(index, int(settings.YEAR) + 1))
                names = np.array([cube.var_name for cube in cube_list])
                selected_cube, = np.where(names == "Ann")[0]

                ghcndex_cube = cube_list[selected_cube]
                ghcndex_cube.coord('latitude').guess_bounds()
                ghcndex_cube.coord('longitude').guess_bounds() 
                ghcndex_cube = fix_time_coord(ghcndex_cube)

                if index in ["TX90p", "TN90p", "TX10p", "TN10p"]:
                    # change from % to days
                    ghcndex_cube.data = ghcndex_cube.data * 3.65

                # era5 cube
                cube_list = iris.load(settings.ERA5LOC_TEX + "ERA5_{}_1979-{}_land.nc".format(index, settings.YEAR))
                names = np.array([cube.var_name for cube in cube_list])
                selected_cube, = np.where(names == "Ann")[0]

                era5_cube = cube_list[selected_cube]
                era5_cube.coord('latitude').guess_bounds()
                era5_cube.coord('longitude').guess_bounds()  

                if index in ["TX90p", "TN90p", "TX10p", "TN10p"]:
                    # change from % to days
                    era5_cube.data = era5_cube.data * 3.65

                # need to regrid
                era5_cube = era5_cube.regrid(ghcndex_cube, iris.analysis.Linear(extrapolation_mode="mask"))

                coverage_offset, coverage_stdev = compute_coverage_error(ghcndex_cube, era5_cube)
                coverage_stdev *= 2. # 90%, 2s.d.
                axes[ix].fill_between(index_ts.times, index_ts.data-coverage_stdev, index_ts.data+coverage_stdev, color='mistyrose', label="ERA5 coverage uncertainty")


                # red tickmarks
                axes[ix].tick_params(axis='y', colors='red', direction="in")

                # and smoothed
                index_ts.data.fill_value = -99.9
                smoothed = binomialfilter(index_ts.data.filled(), -99.9, 5, pad = False)
                smoothed = np.ma.masked_where(smoothed == -99.9, smoothed)

                axes[ix].plot(index_ts.times, smoothed, "r--", lw = LW)

                # print the index, the current anomaly and the rank information
                # print(index, index_ts.data.compressed()[-1]/36.5, np.argsort(np.argsort(index_ts.data.compressed())) + 1)

                # plot the coverage
                axes[ix+2].plot(cover_ts.times, cover_ts.data, "k:", lw = 2)
                axes[ix+2].yaxis.set_label_position("right")
                axes[ix+2].yaxis.set_ticks_position('right')


            # prettify
            plt.xlim([1950-2, int(settings.YEAR)+2])
            axes[0].set_ylim([15,None])
            if "X" in index:
                axes[0].set_ylim([21,79])
                axes[1].set_ylim([19,59])
            elif "N" in index:
                axes[0].set_ylim([21,79])
                axes[1].set_ylim([11,59])

            for ax in [ax1, ax2]:
                utils.thicken_panel_border(ax)
                ax.yaxis.set_ticks_position("left")

            for ax in [ax3, ax4]:         
                ax.set_ylim([0,98])
                utils.thicken_panel_border(ax)



            fig.text(0.03, 0.5, "Number of Days", va='center', rotation='vertical', fontsize = settings.FONTSIZE, color="r")
            fig.text(0.97, 0.5, "% land covered", va='center', rotation='vertical', fontsize = settings.FONTSIZE)

            for ax in axes:
                for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(settings.FONTSIZE)
                    tick.label2.set_fontsize(settings.FONTSIZE)
            for tick in axes[1].xaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)

            fig.subplots_adjust(left=0.1, right = 0.9, top = 0.98, bottom = 0.05, hspace = 0.001)

            plt.savefig(settings.IMAGELOC+"TEX_{}+{}_ts_ghcndex_uncertainties{}".format(index_pair[0], index_pair[1], settings.OUTFMT))
            plt.close()


    #*************
    # MERRA Siberian heatwave timeseries 
    if True:
        
        indata = np.genfromtxt(DATALOC + "merra_siberia.dat", delimiter=",", skip_header=2)
        
        times = indata[:, 0]
        hwf = utils.Timeseries("HWF", times, indata[:, 1])
        hwm = utils.Timeseries("HWM", times, indata[:, 2])
        plt.clf()
        fig, ax1 = plt.subplots(1, figsize=(8, 6))
        ax2 = ax1.twinx()
        utils.plot_ts_panel(ax1, [hwf, hwm], "-", "temperature", loc=LEGEND_LOC)
        ax1.set_ylim([0, 6])
        utils.plot_ts_panel(ax2, [hwm], "-", "temperature", loc="")
        ax2.set_ylim([6.4, 8.8])

        plt.xlim([1978, int(settings.YEAR)+2])
        
        for tick in ax1.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)

        for tick in ax1.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)
        for tick in ax2.yaxis.get_major_ticks():
            tick.label2.set_fontsize(settings.FONTSIZE)

        ax1.yaxis.set_label_position("left")
        ax1.yaxis.set_ticks_position('left')
        ax1.tick_params(axis='y', direction="in", width=2, length=10, color="b")
        ax1.set_ylabel("HWF", color="b", fontsize = settings.FONTSIZE)

        ax2.yaxis.set_label_position("right")
        ax2.yaxis.set_ticks_position('right')
        ax2.tick_params(axis='y', direction="in", width=2, length=10, color="r")
        ax2.set_ylabel('HWM ($^{\circ}$'+"C)", color="r", fontsize = settings.FONTSIZE)

        plt.title("MERRA-2 over 60-160E, 50-80N, AMJ, land only", fontsize = settings.FONTSIZE)
        plt.savefig(settings.IMAGELOC+"TEX_merra2_siberia_ts{}".format(settings.OUTFMT))

        plt.close()


    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
