#!/usr/local/sci/python
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

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.cm as mpl_cm
import matplotlib as mpl

import copy
import string

import iris
import iris.quickplot as qplt
import cartopy

import datetime as dt
import struct

import utils # RJHD utilities
import settings

data_loc = "/data/local/rdunn/SotC/{}/data/PEX/".format(settings.YEAR)
reanalysis_loc = "/data/local/rdunn/SotC/{}/data/RNL/".format(settings.YEAR)
image_loc = "/data/local/rdunn/SotC/{}/images/".format(settings.YEAR)

CLIMSTART=1961
CLIMEND=1990
SELECTED_YEAR = 2017
THRESHOLD = 0.9

LEGEND_LOC = 'upper left'
BBOX = (0,0.9)
LW = 3


SEASON_DICT = {"MAM":["Mar","Apr","May"], "JJA":["Jun","Jul","Aug"], \
                   "SON":["Sep","Oct","Nov"], "DJF":["Dec","Jan","Feb"]}

SEASON_DICT_ERA = {"MAM":["March","April","May"], "JJA":["June","July","August"], \
                   "SON":["September","October","November"], "DJF":["December","January","February"]}

SEASONS = ["DJF", "MAM", "JJA","SON"]

ETCCDI_INDICES = ["PRCPTOT", "Rx1day", "Rx5day", "R10mm", "R20mm", "R95p"]

ETCCDI_LABELS = {"PRCPTOT" : "Total Precipitation", "Rx1day" : "Maximum 1 day precipitation total", "Rx5day" : "Maximum 5 day precipitation total", "R10mm" : "Number of heavy precipitation days", "R20mm" : "Number of very heavy precipitation days", "R95p" : "Precipitation from very wet days"}

ETCCDI_UNITS = {"PRCPTOT" : "mm", "Rx1day" : "mm", "Rx5day" : "mm", "R10mm" : "days", "R20mm" : "days", "R95p" : "mm"}


DWD_INDICES = ["RX1", "RX5", "R10", "R20", "R95P"]

DWD_LABELS = {"RX1" : "Maximum 1 day precipitation total", "RX5" : "Maximum 5 day precipitation total", "R10" : "Number of heavy precipitation days", "R20" : "Number of very heavy precipitation days", "R95P" : "Average precipitation from very wet days"}

DWD_UNITS = {"RX1" : "mm", "RX5" : "mm", "R10" : "days", "R20" : "days", "R95P" : "mm/day"}

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
def read_ghcnd(filename):


    fieldwidths = (11,5,3,8,8,5,6,3,8,8,10,10,7,3,25)
    fmtstring = ''.join('%ds' % f for f in fieldwidths)

    parse = struct.Struct(fmtstring).unpack_from

    indata = []
    with open(filename, 'r') as infile:
        for ll, line in enumerate(infile):
            fields = parse(line)
               
            indata += [[f.strip() for f in fields]]

    indata = np.array(indata)

    ids = indata[:,0]
    record_y = indata[:,1].astype(int)
    record_m = indata[:,2].astype(int)
    value_mm = indata[:,3].astype(float)
    value_in = indata[:,4].astype(float)
    data_length = indata[:,5].astype(int)
    prev_record_yr = indata[:,6].astype(int)
    prev_record_mth = indata[:,7].astype(int)
    prev_value_mm = indata[:,8].astype(float)
    prev_value_in = indata[:,9].astype(float)
    lats = indata[:,10].astype(float)
    lons = indata[:,11].astype(float)
    elev = indata[:,12].astype(float)
    state = indata[:,13]
    name = indata[:,14]

    return lats, lons, value_mm, prev_value_mm # read_ghcnd

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


#************************************************************************
def read_dwd_percentile(filename):
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
    indata = np.genfromtxt(filename, dtype = (float))

    this_lat = []
    tl = 0
    # process each row, append until have complete latitude band
    for row in indata:
        this_lat += [row]

        if len(this_lat) == longitudes.shape[0]:
            # copy into final array and reset
            data[tl, :] = this_lat
            tl += 1
            this_lat=[]

    # mask the missing values
    data = np.ma.masked_where(data <= -999.000, data)

    cube = utils.make_iris_cube_2d(data, latitudes, longitudes, "R90p", "%")

    return cube # read_dwd_percentile

#
#************************************************************************
def run_all_plots():

    rank_bounds = [-4.5,-3.5,-2.5,-1.5,1.5,2.5,3.5,4.5]

    #*************************
    # GHCNDEX indices

    for index in ETCCDI_INDICES:


        cube_list = iris.load(data_loc + "GHCND_{}_1951-{}_RegularGrid_global_2.5x2.5deg_LSmask.nc".format(index, int(settings.YEAR) + 1))
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
            cmap=settings.COLOURMAP_DICT["precip_sequential"]
        elif index in ["Rx5day"]:
            bounds = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
            cmap=settings.COLOURMAP_DICT["precip_sequential"]
        elif index in ["R10mm"]:
            bounds = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
            cmap=settings.COLOURMAP_DICT["precip_sequential"]
        elif index in ["R20mm"]:
            bounds = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
            cmap=settings.COLOURMAP_DICT["precip_sequential"]
        elif index in ["R95p"]:
            bounds = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225]
            cmap=settings.COLOURMAP_DICT["precip_sequential"]
        elif index in ["PRCPTOT"]:
            bounds = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225]
            cmap=settings.COLOURMAP_DICT["precip_sequential"]

        utils.plot_smooth_map_iris(image_loc + "PEX_{}_{}_ghcndex".format(index, settings.YEAR), total_cube[loc[0]], cmap, bounds, "{} ({})".format(index, ETCCDI_UNITS[index]), title = "{} - {}".format(index, ETCCDI_LABELS[index]))

        # sort the bounds and colourbars
        if index in ["Rx1day"]:
            bounds = [-100, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 100]
            cmap=settings.COLOURMAP_DICT["hydrological"]
        elif index in ["Rx5day"]:
            bounds = [-100, -40, -30, -20, -10, -5, 0, 5, 10, 20, 30, 40, 100]
            cmap=settings.COLOURMAP_DICT["hydrological"]
        elif index in ["R10mm"]:
            bounds = [-20, -8, -6, -4, -2, 0, 2, 4, 6, 8, 20]
            cmap=settings.COLOURMAP_DICT["hydrological"]
        elif index in ["R20mm"]:
            bounds = [-20, -4, -3, -2, -1, 0, 1, 2, 3, 4, 20]
            cmap=settings.COLOURMAP_DICT["hydrological"]
        elif index in ["R95p"]:
            bounds = [-1000, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 1000]
            cmap=settings.COLOURMAP_DICT["hydrological"]
        elif index in ["PRCPTOT"]:
            bounds = [-1000, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 1000]
            cmap=settings.COLOURMAP_DICT["hydrological"]

        utils.plot_smooth_map_iris(image_loc + "PEX_{}_{}_anoms_ghcndex".format(index, settings.YEAR), anoms[loc[0]], cmap, bounds, "{} ({})".format(index, ETCCDI_UNITS[index]), title = "{} - {}".format(index, ETCCDI_LABELS[index]))
        if index == "Rx1day":
            utils.plot_smooth_map_iris(image_loc + "p2.1_PEX_{}_{}_anoms_ghcndex".format(index, settings.YEAR), anoms[loc[0]], cmap, bounds, "{} ({})".format(index, ETCCDI_UNITS[index]), figtext = "(k) Maximum 1 Day Precipitation Amount")

        rank_cube = get_ranks(anoms)

        plot_rank_map(image_loc + "PEX_{}_{}_rank_ghcndex".format(index, settings.YEAR), rank_cube[loc[0]], cmap, rank_bounds, "Rank", title = "{} - {}".format(index, ETCCDI_UNITS[index]))
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

                        if month  == "Dec":
                            # need to extract from previous year - cheat by rolling data around
                            cube.data = np.roll(cube.data, 1, axis = 0)
                            cube.data.mask[0,:,:] = True # and mask out the previous years'

                        month_data += [cube.data]

                    # finished getting all months, make a dummy cube to populate
                    month_data = np.ma.array(month_data)
                    season_cube = copy.deepcopy(cube)

                    # take appropriate seasonal value
                    season_cube.data = np.ma.max(month_data, axis = 0)

                    # mask if fewer that 2 months present
                    nmonths_locs = np.ma.count(month_data, axis = 0)
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
                        cmap=settings.COLOURMAP_DICT["precip_sequential"]

                    # pass to plotting routine
                    utils.plot_smooth_map_iris_multipanel(image_loc + "PEX_{}_{}_seasons_ghcndex".format(index, settings.YEAR), season_list, cmap, bounds, "{} ({})".format(index, ETCCDI_UNITS[index]), shape = (2,2), title = SEASONS, figtext = ["(a)","(b)","(c)","(d)"], figtitle = "{} - {}".format(index, ETCCDI_LABELS[index]))
                elif sc == 1:
                    cmap=settings.COLOURMAP_DICT["hydrological"]
                    if index in ["Rx1day"]:
                        bounds = [-100, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 100]
                    elif index in [ "Rx5day"]:
                        bounds = [-100, -20, -16, -12, -8, -4, 0, 4, 8, 12, 16, 20, 100]

                    # pass to plotting routine
                    utils.plot_smooth_map_iris_multipanel(image_loc + "PEX_{}_{}_anoms_seasons_ghcndex".format(index, settings.YEAR), season_list, cmap, bounds, "{} ({})".format(index, ETCCDI_UNITS[index]), shape = (2,2), title = SEASONS, figtext = ["(a)","(b)","(c)","(d)"], figtitle = "{} - {}".format(index, ETCCDI_LABELS[index]))


    #*************************
    # DWD indices
            
    for index in DWD_INDICES:
        cube_list = iris.load(data_loc + "first_guess_daily_{}_{}.nc".format(settings.YEAR, index))

        if len(cube_list) == 1:
            cube = cube_list[0]
        else:
            for cube in cube_list:
                if cube.units == "mm per 5 day" and index == "RX5":
                    break

        cube = cube[0] # take only single slice
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()  
        
        if index in ["RX1"]:
            bounds = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
            cmap=settings.COLOURMAP_DICT["precip_sequential"]
        if index in ["RX5"]:
            bounds = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450]
            cmap=settings.COLOURMAP_DICT["precip_sequential"]
        elif index in ["R10", "R20", "R95P"]:
            bounds = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
            cmap=settings.COLOURMAP_DICT["precip_sequential"]

        utils.plot_smooth_map_iris(image_loc + "PEX_{}_{}_dwd".format(index, settings.YEAR), cube, cmap, bounds, "{} ({})".format(index, DWD_UNITS[index]), title = "{} - {}".format(index, DWD_LABELS[index]))
    

    #*************************
    # DWD percentile

    dwd_cube = read_dwd_percentile(data_loc + "perzentile_201701-201712.ras")
    
    bounds = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    cmap=settings.COLOURMAP_DICT["hydrological"]
    utils.plot_smooth_map_iris(image_loc + "PEX_{}_{}_dwd".format("R90", settings.YEAR), dwd_cube, cmap, bounds, "{} ({})".format("percentile", "%"), title = "Percentile of the annual precipitation total")
    utils.plot_smooth_map_iris(image_loc + "p2.1_PEX_{}_{}_dwd".format("R90", settings.YEAR), dwd_cube, cmap, bounds, "{} ({})".format("percentile", "%"), figtext = "(j) Percentile of the Annual Precipitation Total")


    #*************************
    # GHCND totals and ratios
    import cartopy.feature as cfeature

    NAMES = {"tx" : "Texas", "as" : "Australia", "rq" : "Puerto Rico"}
    EXTENTS = {"tx" : (-95, 30, [-100, -90, 27, 32]), "as" : (150, -25, [140, 160, -35, -15]), "rq" : (-65, 30, [-70, -60, 15, 20])}

    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                            edgecolor='face',
                                            facecolor=cfeature.COLORS['land'])
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                            edgecolor='face',
                                            facecolor=cfeature.COLORS['land'])
    
    for index in ["Rx5day", "Rx1day"]:
        for state in ["tx-201708", "as-201703", "rq-201709"]:

            lats, lons, value_mm, prev_value_mm = read_ghcnd("{}/{}-{}.txt".format(data_loc, state, index))

            ratio = value_mm/prev_value_mm

            # set up figure
            if "as" in state:
                fig = plt.figure(figsize =(10,6))
            else:
                fig = plt.figure(figsize =(10,4))


            if index == "Rx5day":
                bounds = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900]
            elif index == "Rx1day":
                bounds = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450]
            cmap = settings.COLOURMAP_DICT["precip_sequential"]
            norm=mpl.cm.colors.BoundaryNorm(bounds,cmap.N)

#            ratio_cmap = plt.cm.YlOrRd
#            ratio_bounds = [1, 1.05, 1.1, 1.15, 1.2, 1.3, 1.5, 1.75, 2.0 ]
#            ratio_cmap = settings.COLOURMAP_DICT["hydrological"]
#            ratio_bounds = [0.01, 0.5, 0.75, 0.9, 0.95, 1.0, 1.05, 1.1, 1.25, 2.0, 10.0]
            ratio_cmap = plt.cm.Blues
            ratio_bounds = [0.01, 0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.5, 2.0, 10.0]
            ratio_norm=mpl.cm.colors.BoundaryNorm(ratio_bounds,ratio_cmap.N)

            plt.clf()

            # set up axes
            ax0 = plt.axes([0.01, 0.12, 0.45, 0.85], projection=cartopy.crs.Stereographic(central_longitude=EXTENTS[state[:2]][0], central_latitude = EXTENTS[state[:2]][1]))
            ax1 = plt.axes([0.51, 0.12, 0.45, 0.85], projection=cartopy.crs.Stereographic(central_longitude=EXTENTS[state[:2]][0], central_latitude = EXTENTS[state[:2]][1]))

            for ax in [ax0, ax1]:
                ax.set_extent(EXTENTS[state[:2]][2], cartopy.crs.PlateCarree())
                
                states_provinces = cfeature.NaturalEarthFeature(
                    category='cultural',
                    name='admin_1_states_provinces_lines',
                    scale='50m',
                    facecolor='none')
                ax.add_feature(states_provinces, edgecolor='gray')
                ax.coastlines(resolution = "50m")
                ax.add_feature(land_50m, zorder = 0, facecolor = "0.9", edgecolor = "k")

                # add other features
                ax.gridlines() #draw_labels=True)
                ax.add_feature(cartopy.feature.BORDERS, zorder = 0, facecolor = "0.9", edgecolor = "k")


            #  plot 
            for ax, data, cm, bnds, nrm, label in zip([ax0, ax1], [value_mm, ratio], [cmap, ratio_cmap], [bounds, ratio_bounds], [norm, ratio_norm], [index, "Ratio to previous record"]):

                scatter = ax.scatter(lons, lats, c = data, cmap = cm, norm = nrm, s=50, \
                                         transform = cartopy.crs.Geodetic(), edgecolor = '0.5', linewidth='0.5')


                # thicken border of colorbar and the dividers
                # http://stackoverflow.com/questions/14477696/customizing-colorbar-border-color-on-matplotlib
                if "Ratio" in label:
                    cb=fig.colorbar(scatter, ax = ax, orientation = 'horizontal', pad = 0.05, fraction = 0.05, \
                                    aspect = 30, ticks = bnds[1:-1], label = label, drawedges=True)
                    cb.set_ticklabels(["{:g}".format(b) for b in bnds[1:-1]])                
                else:
                    cb=fig.colorbar(scatter, ax = ax, orientation = 'horizontal', pad = 0.05, fraction = 0.05, \
                                    aspect = 30, ticks = bnds[1:-1], label = label, drawedges=True)
                    cb.set_ticklabels(["{:g}".format(b) for b in bnds])

                cb.outline.set_linewidth(2)
                cb.dividers.set_color('k')
                cb.dividers.set_linewidth(2)


            plt.savefig(image_loc + "PEX_{}_{}-{}".format(index, NAMES[state[:2]], state[-6:]) + settings.OUTFMT)
            plt.close()


    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************