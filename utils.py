#!/usr/bin/env python
#************************************************************************
#
#  Utility scripts for BAMS SotC 2015
#
#************************************************************************
#                    SVN Info
# $Rev:: 31                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2021-09-06 09:52:46 +0100 (Mon, 06 Sep #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import math
import os
import sys
import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
import scipy.stats

import iris
import iris.plot
import cartopy
import cf_units

import settings # RJHD settings


#************************************************************************
class Timeseries(object):
    '''
    Class for timeseries
    '''

    def __init__(self, name, times, data):
        self.name = name
        self.times = times
        self.data = data


    def __str__(self):
        return "timeseries of {}".format(self.name)

    __repr__ = __str__




#************************************************************************
def scatter_plot_map(outname, data, lons, lats, cmap, bounds, cb_label, title="", figtext=""):
    '''
    Standard scatter map

    :param str outname: output filename root
    :param array data: data to plot
    :param array lons: longitudes
    :param array lats: latitudes
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

    scatter = plt.scatter(lons, lats, c=data, cmap=cmap, norm=norm, s=25, \
                        transform=cartopy.crs.Geodetic(), edgecolor='0.5', linewidth=0.5)


    cb=plt.colorbar(scatter, orientation='horizontal', pad=0.05, fraction=0.05, \
                        aspect=30, ticks=bounds[1:-1], label=cb_label, drawedges=True)

    # thicken border of colorbar and the dividers
    # http://stackoverflow.com/questions/14477696/customizing-colorbar-border-color-on-matplotlib
    cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
    cb.ax.tick_params(axis='x', labelsize=settings.FONTSIZE, direction='in', size=0)

    cb.set_label(label=cb_label, fontsize=settings.FONTSIZE)

#    cb.outline.set_color('k')
    cb.outline.set_linewidth(2)
    cb.dividers.set_color('k')
    cb.dividers.set_linewidth(2)

    ax.set_extent(ext, ax.projection) # fix the extent change from colormesh

    plt.title(title)
    fig.text(0.03, 0.95, figtext, fontsize=settings.FONTSIZE * 0.8)

    plt.savefig(outname + settings.OUTFMT)
    plt.close()

    return # scatter_plot_map

#************************************************************************
def annual_average(indata):
    '''
    Calculates the annual average of monthly data

    Checks for array shape if necessary

    :param array indata: input data

    :returns: annual mean array
    '''

    # if a 1-D array, reshape
    if len(indata.shape) == 1:
        newdata = np.ma.reshape(indata, [-1, 12])

    return np.ma.mean(newdata, axis=1) # annual_average

#************************************************************************
def erai_ts_read(data_loc, variable, annual=False):
    '''
    read ERAI data and return

    :param str data_loc: data_loc from where to read file (determined in script)
    :param str variable: variable name to be processed
    :param bool annual: return annual means

    :returns" years, globe, oceans, land, tropics - arrays of timeseries.

    '''
    # set up filenames
    if variable == "wnd":
        filename = "tsmm1d_10SI_197901-{}12.txt".format(settings.YEAR)
    elif variable == "sat":
        filename = "tsmm1d_2T_197901-{}12.txt".format(settings.YEAR)

    all_era = np.genfromtxt(data_loc + filename, skip_header=2, skip_footer=29, dtype=(float))

    # extract the data
    date = all_era[:, 1]
    globe = all_era[:, 2]
    oceans = all_era[:, 3]
    land = all_era[:, 4]
    tropics = all_era[:, 5]

    # make annual averages
    if annual:
        globe = annual_average(globe)
        land = annual_average(land)
        oceans = annual_average(oceans)
        tropics = annual_average(tropics)

        date = (np.reshape(date, [-1, 12])[:, 0] / 100.) - 0.01

    return Timeseries("ERA-Interim", date, globe), \
        Timeseries("ERA-Interim", date, oceans), \
        Timeseries("ERA-Interim", date, land),\
        Timeseries("ERA-Interim", date, tropics) # erai_ts_read


#************************************************************************
def erai_2dts_read(data_loc, variable):
    '''
    read ERA data and returns hovmuller data

    :param str variable: variable name to be processed

    :returns:

    '''
    # set up filenames
    if variable == "wnd":
        filename = "timtw_ERA-Int_10SI_197901-{}12.txt".format(settings.YEAR)

    elif variable == "sat":
        filename = "timtw_ERA-Int_2T_197901-{}12.txt".format(settings.YEAR)

    elif variable == "ltt":
        filename = "timtw_ERA-Int_TLT_197901-{}12.txt".format(settings.YEAR)

    all_era = np.genfromtxt(data_loc + filename, dtype=(float))

    times = all_era[:, 0]
    data = all_era[:, 1:]

    n_grids = data.shape[1]
    grid_cell = 180. / (n_grids - 1)
    print("presuming 90 - -90 - swapped Feb2018")

    latitudes = np.arange(90.0, -90.0 - grid_cell, -grid_cell)

    return times, latitudes, data # era_2dts_read

#************************************************************************
def era5_ts_read(data_loc, variable, annual=False):
    '''
    read ERA data and return

    :param str data_loc: data_loc from where to read file (determined in script)
    :param str variable: variable name to be processed
    :param bool annual: return annual means

    :returns" years, globe, oceans, land, tropics - arrays of timeseries.

    '''
    # set up filenames
    if variable == "wnd":
        var = "WS10"
        var = "10SI"
        loc = "ERA5/SFCWIND"
        start = 1979
    elif variable == "sat":
        var = "T2M"
        var = "2T"
        loc = "ERA5/SFCTEMP"
        start = 1967

    cube_list = iris.load(os.path.join(data_loc, loc, "era5_{}_{}01-{}12_mon_area_avg_series.nc".format(var.lower(), start, settings.YEAR)))

    names = np.array([cube.name() for cube in cube_list])

    for cube in cube_list:

        if cube.var_name == "{}_global".format(var.lower(),):
            # global
            globe = cube.data[:].squeeze()
            timeUnits = cube.coord("time").units
            dt_time = timeUnits.num2date(cube.coord("time").points)
        elif cube.var_name == "{}_land".format(var.lower(),):
            # land
            land = cube.data[:].squeeze()
        elif cube.var_name == "{}_ocean".format(var.lower(),):
            # ocean
            oceans = cube.data[:].squeeze()
        elif cube.var_name == "{}_20S20N".format(var.lower(),):
            # tropics
            tropics = cube.data[:].squeeze()

    # make annual averages
    if annual:
        globe = annual_average(globe)
        land = annual_average(land)
        oceans = annual_average(oceans)
        tropics = annual_average(tropics)
        
        # and extract the date
        date = np.array([y.year for y in np.reshape(dt_time, [-1, 12])[:, 0]])
        
    else:
        # extract the date
        date = np.array([(date.year + (date.month - 1)/12.)  for date in dt_time ])

    return Timeseries("ERA5", date, globe), \
        Timeseries("ERA5", date, oceans), \
        Timeseries("ERA5", date, land),\
        Timeseries("ERA5", date, tropics) # era5_ts_read

#************************************************************************
def calculate_climatology_and_anomalies_1d(indata, start, end):
    '''
    Calculate the climatology and anomalies for a 1-d timeseries

    :param array indata: input data - Timeseries object
    :param float start: start year
    :param float end: end year

    :returns: climatology and anomaly Timeseries object
    '''
    # find range to use
    start_loc, = np.where(indata.times == start)
    end_loc, = np.where(indata.times == end + 1)

    if len(start_loc) == 0:
        # series starts too late
        start_loc = [0]
    if len(end_loc) == 0:
        # series ends too early
        end_loc = [-1]

    # get the mean and subtract to get anomalies
    try:
        climatology = np.mean(indata.data[start_loc : end_loc]) # include last year
    except TypeError:
        climatology = np.mean(indata.data[start_loc[0] : end_loc[0]]) # include last year


    outdata = Timeseries(indata.name, indata.times, indata.data - climatology)

    return climatology, outdata # apply_climatology

#***************************************
def median_pairwise_slopes(xdata, ydata, mdi, sigma=1.0):
    '''
    Calculate the median of the pairwise slopes - assumes no missing values

    :param array xdata: x array
    :param array ydata: y array
    :param float mdi: missing data indicator
    :param float sigma: std range for upper/lower
    :returns: float of slope
    '''
    # run through all pairs, and obtain slope
    slopes = []
    for i in range(len(xdata)):
        for j in range(i+1, len(xdata)):
            if ydata[i] > mdi and ydata[j] > mdi:
                slopes += [(ydata[j] - ydata[i]) / (xdata[j] - xdata[i])]

    mpw = np.median(np.array(slopes))

    # copied from median_pairwise.pro methodology (Mark McCarthy)
    slopes.sort()

    good_data = np.where(ydata != mdi)[0]

    n = len(ydata[good_data])

    dof = n * (n - 1) // 2
    w = math.sqrt(n * (n - 1) * ((2. * n) + 5.) / 18.)

    percentile_point = scipy.stats.norm.cdf(sigma)*2.

    rank_upper = ((dof + percentile_point * w) / 2.) + 1
    rank_lower = ((dof - percentile_point * w) / 2.) + 1

    if rank_upper >= len(slopes): rank_upper = len(slopes)-1
    if rank_upper < 0: rank_upper = 0
    if rank_lower < 0: rank_lower = 0

    upper = slopes[int(rank_upper)]
    lower = slopes[int(rank_lower)]

    return mpw, lower, upper      # MedianPairwiseSlopes


#************************************************************************
def fit_plot_points(slope, intercept, years):
    """
    Calculate start and end points for line describing best fit

    :param float slope: line slope
    :param float intercept: line intercept

    :returns: yvalues
    """

    yvalues = []

    for y in years:

        yvalues += [slope * y + intercept]

    return yvalues # fit_plot_points

#************************************************************************
def mpw_plot_points(slope, years, values):
    """
    Calculate start and end points for line describing MPW best fit

    :param float slope: line slope
    :param float array years: x-coordinates
    :param float array values: y-coordinates

    :returns: [[x1,x2], [y1,y1]]
    """

    mu_x = np.ma.mean(years)
    mu_y = np.ma.mean(values)

    y1 = slope * (years[0] - mu_x) + mu_y
    y2 = slope * (years[-1] - mu_x) + mu_y

    return [years[0], years[-1]], [y1, y2] # mpw_plot_points


#************************************************************************
def plot_smooth_map_iris(outname, cube, cmap, bounds, cb_label, scatter=[], smarker="o", ssize=25,\
                             figtext="", title="", contour=False, cb_extra="", save_netcdf_filename="", \
                         hatch=None, tall=False):
    '''
    Standard scatter map

    :param str outname: output filename root
    :param array cube: cube to plot
    :param obj cmap: colourmap to use
    :param array bounds: bounds for discrete colormap
    :param str cb_label: label for colorbar
    :param str smarker: scatter marker type
    :param int ssize: scatter marker sizes
    :param str figtext: add text to figure label (top left)
    :param str title: add title
    :param bool contour: plot smooth rather than gridded data
    :param str cb_label: colorbar label
    :param str save_netcdf_filename: filename to save the output plot to a cube.
    :param cube hatch: cube to plot as hatching
    :param bool tall: make a taller map for longer colorbar label (with line-break)
    '''

    norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)

    if tall:
        fig = plt.figure(figsize=(8, 5.7))
        plt.clf()
        ax = plt.axes([0.01, 0.18, 0.98, 0.78], projection=cartopy.crs.Robinson())
    else:
        fig = plt.figure(figsize=(8, 5.5))
        plt.clf()
        ax = plt.axes([0.01, 0.12, 0.98, 0.88], projection=cartopy.crs.Robinson())
        
    ax.gridlines() #draw_labels=True)
    ax.add_feature(cartopy.feature.LAND, zorder=-1, facecolor="0.9", edgecolor="k")
    ax.coastlines()

    ext = ax.get_extent() # save the original extent

    if contour:
        # http://stackoverflow.com/questions/12274529/how-to-smooth-matplotlib-contour-plot
        from scipy.ndimage.filters import gaussian_filter
        sigma = 0.5
        con_cube = copy.deepcopy(cube)

        """
        gaussian_filter doesn't work with masked arrays
        As this messes up things next to masked values,
        go through each value in the array, if it is masked
        replace with mean of surrounding 8 values, as long
        as this itself is not masked

        3-Aug-2017
        """

        try:
            new_data = copy.deepcopy(con_cube.data)
            orig_mask = copy.deepcopy(con_cube.data.mask)

            # spin through each grid-cell
            for lt in range(con_cube.data.shape[0]):
                for ln in range(con_cube.data.shape[1]):
                    # if missing
                    if con_cube.data.data[lt, ln] == con_cube.data.fill_value:
                        # take square of values (hope wrapping isn't an issue)
                        square = con_cube.data[lt-1:lt+2, ln-1:ln+2]
                        # if there is data to mean, then do and overwrite
                        if len(square.compressed()) > 0:
                            new_data[lt, ln] = np.ma.mean(square)

            # replace the data
            con_cube.data = new_data
            con_cube.data = gaussian_filter(con_cube.data, sigma)
            # re-create the mask
            con_cube.data.mask = orig_mask
        except AttributeError:
            # no masked data
            con_cube.data = gaussian_filter(con_cube.data, sigma)


        mesh = iris.plot.contourf(con_cube, bounds, cmap=cmap, norm=norm)

        if save_netcdf_filename != "":
            save_cube_as_netcdf(con_cube, save_netcdf_filename)

    else:
        if settings.OUTFMT in [".eps", ".pdf"]:
            if cube.coord("latitude").points.shape[0] > 180 or cube.coord("longitude").points.shape[0] > 360:
                regrid_size = 1.0
                print("Regridding cube for {} output to {} degree resolution".format(settings.OUTFMT, regrid_size))
                print("Old Shape {}".format(cube.data.shape))
                plot_cube = regrid_cube(cube, regrid_size, regrid_size)
                print("New Shape {}".format(plot_cube.data.shape))
            else:
                plot_cube = copy.deepcopy(cube)
        else:
            plot_cube = copy.deepcopy(cube)
        
        mesh = iris.plot.pcolormesh(plot_cube, cmap=cmap, norm=norm)

        if save_netcdf_filename != "":
            save_cube_as_netcdf(plot_cube, save_netcdf_filename)

    if hatch != None:
        iris.plot.pcolor(hatch, hatch='///', cmap=cmap, norm=norm)

    if len(scatter) > 0:
        lons, lats, data = scatter
        if smarker == "o":
            plt.scatter(lons, lats, c=data, cmap=cmap, norm=norm, s=ssize, marker="o",\
                        transform=cartopy.crs.Geodetic(), edgecolor='0.2', linewidth=0.5, zorder=10)
        elif smarker == "dots":
            plt.scatter(lons, lats, c="0.1", cmap=cmap, norm=norm, s=2, marker="o",\
                        transform=cartopy.crs.Geodetic())


    cb = plt.colorbar(mesh, orientation='horizontal', pad=0.05, fraction=0.05, \
                        aspect=30, ticks=bounds[1:-1], label=cb_label, drawedges=True)

    # thicken border of colorbar and the dividers
    # http://stackoverflow.com/questions/14477696/customizing-colorbar-border-color-on-matplotlib
    cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
    cb.ax.tick_params(axis='x', labelsize=settings.FONTSIZE, direction='in', size=0)

    cb.set_label(label=cb_label, fontsize=settings.FONTSIZE)

#    cb.outline.set_color('k')
    cb.outline.set_linewidth(2)
    cb.dividers.set_color('k')
    cb.dividers.set_linewidth(2)

    if cb_extra != "":
        fig.text(0.04, 0.07, cb_extra[0], fontsize=settings.FONTSIZE, ha="left")
        fig.text(0.96, 0.07, cb_extra[1], fontsize=settings.FONTSIZE, ha="right")

    ax.set_extent(ext, ax.projection) # fix the extent change from colormesh

    plt.title(title, fontsize=settings.FONTSIZE)
    fig.text(0.01, 0.95, figtext, fontsize=settings.FONTSIZE)

    plt.savefig(outname + settings.OUTFMT)
    plt.close()

    return # plot_smooth_map_iris


#************************************************************************
def plot_smooth_map_iris_multipanel(outname, cube_list, cmap, bounds, cb_label, \
                                        shape=(1, 1), scatter=[], figtext=[], title=[], figtitle=""):
    '''
    Standard scatter map

    :param str outname: output filename root
    :param array cube: data to plot
    :param obj cmap: colourmap to use
    :param array bounds: bounds for discrete colormap
    :param str cb_label: colorbar label
    '''
    number_of_panels = len(cube_list)

    assert shape[0] * shape[1] == number_of_panels
        
    norm=mpl.cm.colors.BoundaryNorm(bounds, cmap.N)

    height = 8
    if shape[0] == 2:
        height = 6
    if shape[0] == 4:
        height = 12
    if shape[0] == 6:
        height = 16
    width = 8
#    if shape[1] == 2:
#        width = 10
#    if shape[1] == 1:
#        width = 5

    fig = plt.figure(figsize=(width, height))


    # spin through panels - rows then columns
    for panel in range(number_of_panels):

        ax = plt.subplot(shape[0], shape[1], panel + 1, projection=cartopy.crs.Robinson())

        ax.gridlines() #draw_labels=True)
        ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
        ax.coastlines()

        ext = ax.get_extent() # save the original extent

        cube = cube_list[panel]

        if settings.OUTFMT in [".eps", ".pdf"]:
            if cube.coord("latitude").points.shape[0] > 180 or cube.coord("longitude").points.shape[0] > 360:
                print("Regridding cube for {} output to 1.0 degree resolution".format(settings.OUTFMT))
                print("Old Shape {}".format(cube.data.shape))
                plot_cube = regrid_cube(cube, 1.0, 1.0)
                print("New Shape {}".format(plot_cube.data.shape))
            else:
                plot_cube = copy.deepcopy(cube)
        else:
            plot_cube = copy.deepcopy(cube)

        mesh = iris.plot.pcolormesh(plot_cube, cmap=cmap, norm=norm, axes=ax)

        if len(scatter) > 0:
            if len(scatter[panel]) > 0:

                lons, lats, data = scatter[panel]

                ax.scatter(lons, lats, c=data, cmap=cmap, norm=norm, s=25, \
                            transform=cartopy.crs.Geodetic(), edgecolor='0.5', linewidth=0.5)

        if len(title) > 0:
            ax.title.set_text(title[panel])
            ax.title.set_fontsize(settings.FONTSIZE)
        if len(figtext) > 0:
            ax.text(0.03, 0.95, figtext[panel], fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)

        ax.set_extent(ext, ax.projection) # fix the extent change from colormesh

    if height == 16:
        cbar_ax = fig.add_axes([0.05, 0.05, 0.9, 0.02])
    elif height == 12:
        cbar_ax = fig.add_axes([0.05, 0.07, 0.9, 0.03])
    elif height == 6:
        cbar_ax = fig.add_axes([0.05, 0.09, 0.9, 0.04])

    cb = fig.colorbar(mesh, cax=cbar_ax, orientation='horizontal', \
                         ticks=bounds[1:-1], label=cb_label, drawedges=True)

    # thicken border of colorbar and the dividers
    # http://stackoverflow.com/questions/14477696/customizing-colorbar-border-color-on-matplotlib
    cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
    cb.ax.tick_params(axis='x', labelsize=settings.FONTSIZE, direction='in', size=0)

    cb.set_label(label=cb_label, fontsize=settings.FONTSIZE)

#    cb.outline.set_color('k')
    cb.outline.set_linewidth(2)
    cb.dividers.set_color('k')
    cb.dividers.set_linewidth(2)

    plt.figtext(0.5, 0.95, figtitle, fontsize=settings.FONTSIZE, ha="center")

    if height == 16:
        fig.subplots_adjust(right=0.99, top=0.98, bottom=0.10, left=0.01, hspace=0.3, wspace=0.01)
    else:
        fig.subplots_adjust(right=0.99, top=0.92, bottom=0.15, left=0.01, hspace=0.01, wspace=0.01)

    plt.savefig(outname + settings.OUTFMT)
    plt.close()

    return # plot_smooth_map_iris_multipanel

#************************************************************************
def read_merra(filename, variable, domain, anomalies=True):
    """
    Read Merra from Mike B's large file

    :param str filename: data file name
    :param str variable: which variable to read
    :param str domain: which domain (LO = Land+Ocean, L = Land only, O = Ocean only)
    :param bool anomaly: read in anomalies
    """
    if domain not in ['LO', 'L', 'O']: sys.exit("Specify domain as L, O or LO")

    # set up which column to read
    if variable == "temperature":
        col = "T2m"
    elif variable == "wind":
        col = "WS10m"
    elif variable == "humidity":
        col = "Qv2m"
    else:
        sys.exit("badly specified variable")

    if anomalies:
        col += "a"
    else:
        col += "m"

    if domain == "L":
        col += "L"
    elif domain == "O":
        col += "O"

    # read the data
    indata = np.genfromtxt(filename, delimiter=',', dtype=(str), skip_header=1)

    headings = np.array([x.strip() for x in indata[0, :]])
    locs, = np.where(headings == col)

    if len(locs) != 1:
        # then only take first
        locs = locs[0]

    months = indata[1:, 0]

    # extract what's necessary - reread the file to get as floats
    indata = np.genfromtxt(filename, delimiter=',', dtype=(float), skip_header=2)
    indata = np.ma.masked_where(indata == -9.99000000e+08, indata)

    values = indata[:, locs]

    # convert to annuals
    months = np.array([int(x[-4:]) for x in months])

    months = months.reshape([-1, 12])
    years = months[:, 1]

    values = values.reshape([-1, 12])

    annuals = np.mean(values, axis=1)

    return Timeseries("MERRA-2", years, annuals) # read_merra

#************************************************************************
def read_jra55(filename, variable):
    """
    JRA-55 data to be read (from Shinya)

    :param str filename: file to read

    :returns years, actuals and anomalies
    """

    name = "JRA-55"

    indata = np.genfromtxt(filename, dtype=(float))

    if variable == "temperature":
        actuals = Timeseries(name, indata[:, 0], indata[:, 1] - 273.1)
    else:
        actuals = Timeseries(name, indata[:, 0], indata[:, 1])

    anomalies = Timeseries(name, indata[:, 0], indata[:, 2])

    return actuals, anomalies

#************************************************************************
def read_20cr(filename, variable):
    """
    20CR data to be read (from Gil & Cathy)

    :param str filename: file to read

    :returns years, actuals and anomalies
    """

    name = "20CRv3"

    indata = np.genfromtxt(filename, dtype=(float))
    
    times = indata[:, 0].astype(int).astype(str)
    years = np.array([t[:4] for t in times]).astype(int)
    months = np.array([t[4:] for t in times]).astype(int)
    
    times = years + (months - 1)/12.    

    monthly = np.ma.masked_where(indata[:, 1] == -9999, indata[:, 1])
    annuals = annual_average(monthly)

    if variable == "temperature":
        actuals = Timeseries(name, times[::12], annuals - 273.1)
    else:
        actuals = Timeseries(name, times[::12], annuals)

    return actuals # read_20cr


#************************************************************************
def read_merra_LT_LS(filename, LT=False, LS=False):
    """
    MERRA Lower Trop/Strat in separate files and different formats

    :param str filename: file to read
    :param bool LT: return lower troposphere
    :param bool LS: reutrn lower stratosphere

    :returns years, LT and LS actuals and anomalies
    """

    name = "MERRA-2"

    indata = np.genfromtxt(filename, delimiter=',', dtype=(float), skip_header=2, skip_footer=12)

    LT_actuals = Timeseries(name, indata[:, 0], indata[:, 1] - 273.1)
    LT_anoms = Timeseries(name, indata[:, 0], indata[:, 2])
    LS_actuals = Timeseries(name, indata[:, 0], indata[:, 3] - 273.1)
    LS_anoms = Timeseries(name, indata[:, 0], indata[:, 4])

    if LT and LS or (not LT and not LS):
        input("Choose one of LT or LS")
    elif LT:
        return LT_actuals, LT_anoms
    elif LS:
        return LS_actuals, LS_anoms # read_merra_LT_LS


#************************************************************************
def thicken_panel_border(ax):
    """
    Thicken the border around the timeseries panels.
    """


    for loc in ["top", "bottom", "left", "right"]:
        ax.spines[loc].set_linewidth(2)

    ax.xaxis.set_tick_params(top=True, which="both", width=2, direction="in")
#    ax.yaxis.set_tick_params(right=False, which="both", width=2, direction="in")
    ax.yaxis.set_tick_params(right=True, which="both", width=2, direction="in")
   
    ax.xaxis.set_tick_params(which="minor", length=5)
    ax.yaxis.set_tick_params(which="minor", length=5)

    ax.xaxis.set_tick_params(which="major", length=10)
    ax.yaxis.set_tick_params(which="major", length=10)

    return

#************************************************************************
def plot_ts_panel(ax, datasets, ls, section, loc="lower right", bbox=(), \
                      minor_tick_interval=1, lw=3, ncol=2, extra_labels=[]):
    """
    Plot panel of a timeseries plot - can be a single one

    :param obj ax: axes object
    :param list datasets: list of Timeseries objects to plot
    :param int ls: linestyle
    :param str loc: legend location
    :param str section: which section of BAMS to look up colour for.
    :param tuple bbox: bounding box to anchor argument
    """
    minorLocator = MultipleLocator(minor_tick_interval)

    COLOURS = settings.COLOURS[section]

    if extra_labels == []:
        extra_labels = ["" for d in datasets]

    assert len(extra_labels) == len(datasets)

    for d, dataset in enumerate(datasets):
        # overwrite
        try:
            LS = dataset.ls
#            print("setting line-style for {} = {}".format(dataset.name, LS))
        except AttributeError:
            LS = ls
            
        try:
            LW = dataset.lw
#            print("setting line-width for {} = {}".format(dataset.name, LW))
        except AttributeError:
            LW = lw

        try:
            ZO = dataset.zorder
#            print("setting zorder for {}".format(dataset.name))
        except AttributeError:
            ZO = 5
        
        ax.plot(dataset.times, dataset.data, c=COLOURS[dataset.name], ls=LS, \
                    label=dataset.name + extra_labels[d], lw=LW, zorder=ZO)

    ax.axhline(0, c='0.5', ls='--')

    if loc != "":
        if len(bbox) in [2, 4]:
            ax.legend(loc=loc, ncol=ncol, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, \
                          labelspacing=0.1, columnspacing=0.5, bbox_to_anchor=bbox)
        else:
            ax.legend(loc=loc, ncol=ncol, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, \
                          labelspacing=0.1, columnspacing=0.5)

    ax.xaxis.set_minor_locator(minorLocator)

    # turn off RHS y ticks
    ax.yaxis.set_ticks_position('left')

    thicken_panel_border(ax)

    return # plot_ts_panel

#************************************************************************
def make_iris_cube_2d(data, lats, lons, name, units):
    """
    Make an Iris cube of a single year of data from arrays of lat, lon and data

    :param array data: data (lat x lon)
    :param array lons: longitudes
    :param array lats: latitudes
    :param str name: name for the cube
    :param str units: units of the data

    :returns: cube - iris cube
    """

    # create the iris cube
    cube = iris.cube.Cube(data)
    cube.rename(name)
    cube.units = units

    # single field, so no time dimension needed

    latcoord = iris.coords.DimCoord(lats, standard_name='latitude', units='degrees')
    loncoord = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')

    cube.add_dim_coord(latcoord, 0)
    cube.add_dim_coord(loncoord, 1)

    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()

    # set coordinate system - 07-08-2014
    cs = iris.coord_systems.GeogCS(6371229)

    cube.coord('latitude').coord_system = cs
    cube.coord('longitude').coord_system = cs

    print("check output as lat/lon order changed 3-Aug-2017")

    return cube # make_iris_cube_2d

#************************************************************************
def make_iris_cube_3d(data, times, time_units, lons, lats, name, units):
    """
    Make an Iris cube of a single year of data from arrays of lat, lon and data

    :param array data: data (lon x lat)
    :param array lons: longitudes
    :param array lats: latitudes
    :param array times: times
    :param str name: name for the cube
    :param str units: units of the data

    :returns: cube - iris cube
    """

    # create the iris cube
    cube = iris.cube.Cube(data)
    cube.rename(name)
    cube.units = units

    time_unit = cf_units.Unit(time_units, calendar=cf_units.CALENDAR_GREGORIAN)
    timecoord = iris.coords.DimCoord(times, standard_name='time', units=time_unit, var_name="time") # add bounds?

    loncoord = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    latcoord = iris.coords.DimCoord(lats, standard_name='latitude', units='degrees')

    cube.add_dim_coord(timecoord, 0)
    cube.add_dim_coord(loncoord, 1)
    cube.add_dim_coord(latcoord, 2)

    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()
    cube.coord('time').guess_bounds()

    # set coordinate system - 07-08-2014
    cs = iris.coord_systems.GeogCS(6371229)

    cube.coord('latitude').coord_system = cs
    cube.coord('longitude').coord_system = cs

    return cube # make_iris_cube_3d

#************************************************************************
def plot_hovmuller(outname, times, latitudes, data, cmap, bounds, cb_label, figtext="", title="", \
                       cosine=False, extra_ts=Timeseries("BLANK", [0], [0]), background=""):
    '''
    Plot a Hovmuller plot (latitude versus time).

    :param str outname: output filename root
    :param array times: time data
    :param array latitudes: latitude data
    :param array data: data to contour
    :param obj cmap: colourmap to use
    :param array bounds: bounds for discrete colormap
    :param str cb_label: colorbar label
    :param bool cosine: cosine scale the latitudes
    :param Timeseries extra_ts: plot extra TS on top
    '''

    minorLocator = MultipleLocator(1)

    # set up the figure

    if extra_ts.name == "BLANK":
        fig = plt.figure(figsize=(8, 6))
        plt.clf()
        ax1 = plt.axes([0.10, 0.10, 0.87, 0.87])
    else:
        fig = plt.figure(figsize=(8, 8))
        plt.clf()
        ax1 = plt.axes([0.10, 0.10, 0.87, 0.67])


    if cosine:
        latitudes = np.sin(np.deg2rad(latitudes))


    # add background to show data coverage
    if background != "":
        ax1.patch.set_facecolor(background)
    # make the basic contour plot
    norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)

    contour = plt.contourf(times, latitudes, data, bounds, cmap=cmap, norm=norm, vmax=bounds[-2], vmin=bounds[1])

    # colourbar and prettify
    cb = plt.colorbar(orientation='horizontal', pad=0.07, fraction=0.05, aspect=30, \
                          ticks=bounds[1:-1], label=cb_label, drawedges=True)

    cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
    cb.ax.tick_params(axis='x', labelsize=settings.FONTSIZE, direction='in', size=0)

    cb.set_label(label=cb_label, fontsize=settings.FONTSIZE)

#    cb.outline.set_color('k')
    cb.outline.set_linewidth(2)
    cb.dividers.set_color('k')
    cb.dividers.set_linewidth(2)


    # plot the timeseries above if needed
    if extra_ts.name != "BLANK":
        # add axis afterwards to make sure original plot works
        ax2 = fig.add_axes([0.10, 0.80, 0.87, 0.17])

        ax2.plot(extra_ts.times, extra_ts.data, c="k", ls="-", label=extra_ts.name, lw=2)
        ax2.axhline(0, ls=":", lw=1, c="k")
#        ax2.set_ylim([-2.5, 2.5])
        ax2.set_xticklabels([""])
        ax2.set_ylabel(extra_ts.name)
        ax2.yaxis.set_ticks_position('left')
        # if two axes, need to add labels
        ax2.text(-0.11, 1.0, "(a)", fontsize=settings.FONTSIZE * 0.8, transform=ax2.transAxes)

        ax1.text(-0.11, 1.0, "(b)", fontsize=settings.FONTSIZE * 0.8, transform=ax1.transAxes)

    else:
        # to make the prettifying easier
        ax2 = fig.add_axes([0.99, 0.99, 0.00, 0.00])
        ax2.set_axis_off()



    for ax in [ax1]:
        if cosine:
            ax.set_ylim(np.sin(np.deg2rad(np.array([-90, 90]))))
            ax.set_yticks(np.sin(np.deg2rad(np.array([-90, -60, -30, 0, 30, 60, 90]))))
            ax.set_yticklabels(["-90"+r'$^{\circ}$'+"S", "-60"+r'$^{\circ}$'+"S", \
                                    "-30"+r'$^{\circ}$'+"S", "0"+r'$^{\circ}$'+"", \
                                    "30"+r'$^{\circ}$'+"N", "60"+r'$^{\circ}$'+"N", "90"+r'$^{\circ}$'+"N"], fontsize=settings.FONTSIZE)
        else:
            ax.set_ylim([-90, 90])
            ax.set_yticks([-60, -30, 0, 30, 60])
            ax.set_yticklabels(["-60"+r'$^{\circ}$'+"S", "-30"+r'$^{\circ}$'+"S", \
                                    "0"+r'$^{\circ}$'+"", "30"+r'$^{\circ}$'+"N", "60"+r'$^{\circ}$'+"N"], fontsize=settings.FONTSIZE)

    # prettify the plot
    for ax in [ax1, ax2]:
        try:
            ax.set_xlim([int(times.compressed()[0]-1), int(times.compressed()[-1]+2)])
        except AttributeError:
            ax.set_xlim([int(times[0]-1), int(times[-1]+2)])

        ax.xaxis.set_minor_locator(minorLocator)
        thicken_panel_border(ax)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)

        if "TWS" in outname:
            majorLocator = MultipleLocator(5)
            ax.xaxis.set_major_locator(majorLocator)
        if "ABD" in outname:
            majorLocator = MultipleLocator(5)
            ax.xaxis.set_major_locator(majorLocator)


    # finish off
    plt.title(title)
    fig.text(0.03, 0.95, figtext, fontsize=settings.FONTSIZE)

    plt.savefig(outname + settings.OUTFMT)
    plt.close()

    return # plot_hovmuller

#*********************************************************
def periodConstraint(cube, t1, t2):
    import iris

    # Time constraint now works on datetime values
#    timeUnits = cube.coord('time').units
#    t1 = timeUnits.date2num(t1)
#    t2 = timeUnits.date2num(t2)
    return iris.Constraint(time=lambda cell: t1 <= cell.point < t2) # periodConstraint

#*********************************************************
def latConstraint(lats):
    return iris.Constraint(latitude=lambda cell: lats[0] <= cell <= lats[1]) # latConstraint

#*********************************************************
def boxcar(data, kernel):
    return np.convolve(data, np.ones((kernel,))/kernel, mode="same") # boxcar

#*********************************************************
def regrid_cube(cube, delta_lat, delta_lon):

    lat, lon = cube.coord('latitude'), cube.coord('longitude')
    lat_min, lon_min = lat.points.min(), lon.points.min()
    lat_max, lon_max = lat.points.max(), lon.points.max()
    newlat = np.arange(lat_min, lat_max, delta_lat)
    newlon = np.arange(lon_min, lon_max, delta_lon)

    return cube.interpolate([('latitude', newlat), ('longitude', newlon)], \
                           iris.analysis.Linear()) # regrid_cube

#*********************************************************
def save_cube_as_netcdf(cube, filename):

    iris.save(cube, filename)

    return # save_cube_as_netcdf

#*********************************************************
def trendline(trend, years):

    x_mean = (years[0]+years[-1])/2.

    c = -trend * x_mean

    y0 = trend*years[0] + c
    y1 = trend*years[-1] + c

    return y0, y1 # trendline
