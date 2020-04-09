#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for surface temperature (SAT) section.
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 27                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2019-08-15 16:09:25 +0100 (Thu, 15 Aug #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************

# python3
from __future__ import absolute_import
from __future__ import print_function
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import iris

import utils # RJHD utilities
import settings

DATALOC = "{}/{}/data/SAT/".format(settings.ROOTLOC, settings.YEAR)
IS_timeseries_root = "globalDatasets2019"
LEGEND_LOC = 'upper left'

LW = 3
BBOX = (0.05, 0.9)
YLIM = [-1.4, 1.1]

# Colin to provide 1981-2010 HadCRUT4 timeseries with uncertainty bounds

# Map files from HadOBS (use N+S averaged),
#            Ahira for MLOST
#            /project/earthobs/GLOBAL_SURFACE_TEMPERATURE/GISTEMP/GISS_1200_blend_1x1.pp fof GISS



#************************************************************************
def read_global_t(filename):
    # Ahira's global temperature files (observed)

    indata = np.genfromtxt(filename, delimiter=',', dtype=(float), skip_header=1)

    indata = np.ma.masked_where(indata == -99.9, indata)

    hadley = utils.Timeseries("Hadley", indata[:, 0], indata[:, 1]) # for completeness in 2016
    noaa = utils.Timeseries("NOAA/NCEI", indata[:, 0], indata[:, 2])
    nasa = utils.Timeseries("NASA/GISS", indata[:, 0], indata[:, 3])
#    jma = utils.Timeseries("JMA", indata[:, 0], indata[:, 4])
    jma = utils.Timeseries("JMA", [0], [0])

    try:
        berkeley = utils.Timeseries("Berkeley", indata[:, 0], indata[:, 5])
        return  noaa, nasa, jma, berkeley
    except IndexError:
        return  noaa, nasa, jma # read_global_t

#************************************************************************
def read_nasa_giss(filename):
    """
    Read the NASA GISS data and returns a cube of the year.

    :param str filename: filename to read

    :returns: cube of 1 year of temperature anomalies
    """

    all_giss = np.genfromtxt(filename, dtype=(float), skip_header=2)

    # 2 degree, but just in case
    # read from i, j columns to get the size of the array
    latitudes = np.zeros(np.max(all_giss[:, 1]).astype(int))
    longitudes = np.zeros(np.max(all_giss[:, 0]).astype(int))

    # set up a masked data array
    data = np.ma.zeros((np.max(all_giss[:, 1]).astype(int), np.max(all_giss[:, 0]).astype(int)))
    data.mask = np.ones(data.shape)

    # spin through each line
    for line in all_giss:
        # use the indexing provided
        i = line[1]
        j = line[0]

        data[i-1, j-1] = line[4]

        # and read in the coordinates too
        if j == 1:
            longitudes[i-1] = line[2]
        if i == 1:
            latitudes[j-1] = line[3]

    # mask the missing data
    data = np.ma.masked_where(data > 1000, data)

    cube = utils.make_iris_cube_2d(data, latitudes, longitudes, "temperature", "C")

    return cube # read_nasa_giss

#************************************************************************
def read_noaa_mlost(filename, year):
    """
    Read the NOAA MLOST data and returns a cube of the year.

    :param str filename: filename to read
    :param int year: year to extract

    :returns: cube of 1 year of temperature anomalies
    """

    all_mlost = np.genfromtxt(filename, dtype=(float))


    DELTA = 5
    # read from i, j columns to get the size of the array
    longitudes = np.arange(-180 + (DELTA/2.), 180 + (DELTA/2.), DELTA)
    latitudes = np.arange(-90 + (DELTA/2.), 90 + (DELTA/2.), DELTA)

    # set up a masked data array
    data = np.ma.zeros((len(latitudes), len(longitudes)))
    data.mask = np.ones(data.shape)

    # spin through each line
    for line in all_mlost:
        if line[0] == year:

            lat_loc, = np.where(latitudes == line[1])
            lon_loc, = np.where(longitudes == line[2])

            data[lat_loc, lon_loc] = line[3]
            data.mask[lat_loc, lon_loc] = False

    cube = utils.make_iris_cube_2d(data, latitudes, longitudes, "temperature", "C")

    return cube # read_noaa_mlost

#************************************************************************
def read_hadcrut_crutem(filename, adjust_clim=False):
    """
    Read data from HadCRUT, HadSST and CRUTEM

    :param str filename: infile to read
    :param bool adjust_clim: adjust climatology if required

    :returns: Timeseries object with upper and lower bounds.
    """


    #******************************************
    def apply_clim(data, upper, lower, start):
        '''Calculate and apply climatology'''

        offset = np.mean(data[start:start+30])

        return data-offset, upper-offset, lower-offset # apply_clim
    #******************************************

    # Colin's global temperature files (observed)

    indata = np.genfromtxt(filename, dtype=(float), delimiter=",", skip_header=3)

    indata = np.ma.masked_where(indata == -99.9, indata)

    years = indata[:, 0]
    mean = indata[:, 1]
    # order can be different for different datasets.
    if "CRUTEM" in filename or "crutem" in filename:
        lower = indata[:, -2]
        upper = indata[:, -1]
        name = "CRUTEM"
    elif "HadSST" in filename or "hadsst" in filename:
        lower = indata[:, -2]
        upper = indata[:, -1]
        name = "HadSST"
    elif "hadcrut4" in filename:
        lower = indata[:, -2]
        upper = indata[:, -1]
        name = "HadCRUT4"

    if adjust_clim:
        # these curves are 1961-1990
        # need to adjust to 1981-2010

        locs, = np.where(years == 1981)

        mean, upper, lower = apply_clim(mean, upper, lower, locs) #  does 30 years from start point

    if years[-1] >= dt.datetime.now().year:
        while years[-1] != dt.datetime.now().year - 1:
            years = years[:-1]
            mean = mean[:-1]
            upper = upper[:-1]
            lower = lower[:-1]

    hadcrut = utils.Timeseries(name, years, mean)
    hadcrut.upper = upper
    hadcrut.lower = lower

    return hadcrut # read_hadcrut_crutem



#************************************************************************
def run_all_plots():

    if True:
        # multipanel timeseries

        COLOURS = settings.COLOURS["temperature"]
        fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, figsize=(8, 19), sharex=True)

        # ERA5
        era5_globe, era5_ocean, era5_land, era5tropics = utils.era5_ts_read(settings.REANALYSISLOC, "sat", annual=True)
        land_era5_clim, land_era5_anoms = utils.calculate_climatology_and_anomalies_1d(era5_land, 1981, 2010)
        ocean_era5_clim, ocean_era5_anoms = utils.calculate_climatology_and_anomalies_1d(era5_ocean, 1981, 2010)
        global_era5_clim, global_era5_anoms = utils.calculate_climatology_and_anomalies_1d(era5_globe, 1981, 2010)

        #*******************
        # in situ L+O
        noaa, nasa, jma = read_global_t(DATALOC + "{}_LO.csv".format(IS_timeseries_root))

        hadcrut = read_hadcrut_crutem(DATALOC+"hadcrut4.1981-2010.csv")

        p0 = ax1.plot(noaa.times, noaa.data, c=COLOURS[noaa.name], ls='-', label=noaa.name, lw=LW)
        p1 = ax1.plot(nasa.times, nasa.data, c=COLOURS[nasa.name], ls='-', label=nasa.name, lw=LW)
   #     p2 = ax1.plot(jma.times, jma.data, c=COLOURS[jma.name], ls='-', label=jma.name, lw=LW)
        p3 = ax1.plot(hadcrut.times, hadcrut.data, c=COLOURS[hadcrut.name], ls='-', label=hadcrut.name, lw=LW)
        ax1.fill_between(hadcrut.times, hadcrut.lower, hadcrut.upper, \
                             where=hadcrut.upper > hadcrut.lower, color='0.5', alpha=0.7)
        p4 = ax1.fill(np.NaN, np.NaN, '0.5', alpha=0.7)

        ax1.axhline(0, c='0.5', ls='--')

        ax1.legend([p0[0], p1[0], (p3[0], p4[0])], [noaa.name, nasa.name, hadcrut.name], \
                       loc=LEGEND_LOC, ncol=2, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, \
                       labelspacing=0.1, columnspacing=0.5, bbox_to_anchor=BBOX)

        ax1.text(0.02, 0.9, "(a) In Situ Land and Ocean", transform=ax1.transAxes, fontsize=settings.FONTSIZE)

        utils.thicken_panel_border(ax1)
        # ax1.yaxis.set_ticks_position('left')

        #*******************
        # reanalysis L+O

        merra = utils.read_merra(settings.REANALYSISLOC + "MERRA-2_SfcAnom{}.dat".format(settings.YEAR), "temperature", "LO")
        jra_actuals, jra_anoms = utils.read_jra55(settings.REANALYSISLOC + "JRA-55_tmp2m_global_ts.txt", "temperature")
        twenty_cr_actuals = utils.read_20cr(settings.REANALYSISLOC + "global.2mt.skt.txt", "temperature")
        dummy, twenty_cr_anoms = utils.calculate_climatology_and_anomalies_1d(twenty_cr_actuals, 1981, 2010)
        twenty_cr_anoms.zorder=-1
        
        # 2018 no MERRA
        utils.plot_ts_panel(ax2, [jra_anoms, global_era5_anoms], "-", "temperature", loc=LEGEND_LOC, bbox=BBOX)

        ax2.text(0.02, 0.9, "(b) Reanalysis Land and Ocean", transform=ax2.transAxes, fontsize=settings.FONTSIZE)

        #*******************
        # in situ L

        noaa, nasa, jma = read_global_t(DATALOC +"{}_L.csv".format(IS_timeseries_root))
        crutem = read_hadcrut_crutem(DATALOC + "crutem4_new_logo.1981-2010.csv")

        p0 = ax3.plot(noaa.times, noaa.data, ls='-', c=COLOURS[noaa.name], label=noaa.name, lw=LW)
        p1 = ax3.plot(nasa.times, nasa.data, ls='-', c=COLOURS[nasa.name], label=nasa.name, lw=LW)
    #    p2 = ax3.plot(jma.times, jma.data, ls='-', c=COLOURS[jma.name], label=jma.name, lw=LW)
    #    p3 = ax3.plot(berkeley.times, berkeley.data, ls = '-', c = COLOURS[berkeley.name], label = berkeley.name, lw = LW)
        p4 = ax3.plot(crutem.times, crutem.data, ls='-', c=COLOURS[crutem.name], label=crutem.name, lw=LW)
        ax3.fill_between(crutem.times, crutem.lower, crutem.upper, \
                             where=crutem.upper > crutem.lower, color='0.5', alpha=0.7)
        p5 = ax1.fill(np.NaN, np.NaN, '0.5', alpha=0.7)

        ax3.axhline(0, c='0.5', ls='--')

        ax3.legend([p0[0], p1[0], (p4[0], p5[0])], [noaa.name, nasa.name, crutem.name], \
                       loc=LEGEND_LOC, ncol=2, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, \
                       labelspacing=0.1, columnspacing=0.5, bbox_to_anchor=BBOX)

        ax3.text(0.02, 0.9, "(c) In Situ Land only", transform=ax3.transAxes, fontsize=settings.FONTSIZE)

        utils.thicken_panel_border(ax3)
    #    ax3.yaxis.set_ticks_position('left')

        #*******************
        # reanalysis L

        merra = utils.read_merra(settings.REANALYSISLOC + "MERRA-2_SfcAnom{}.dat".format(settings.YEAR), "temperature", "L")
        jra_actual, jra_anoms = utils.read_jra55(settings.REANALYSISLOC + "JRA-55_tmp2m_globalland_ts.txt", "temperature")
        twenty_cr_actuals = utils.read_20cr(settings.REANALYSISLOC + "air2mland.txt", "temperature")
        dummy, twenty_cr_anoms = utils.calculate_climatology_and_anomalies_1d(twenty_cr_actuals, 1981, 2010)
        twenty_cr_anoms.zorder=-1
 
        # 2018 - No MERRA
        utils.plot_ts_panel(ax4, [jra_anoms, land_era5_anoms], "-", "temperature", loc=LEGEND_LOC, bbox=BBOX)

        ax4.text(0.02, 0.9, "(d) Reanalysis Land only", transform=ax4.transAxes, fontsize=settings.FONTSIZE)

        #*******************
        # in situ O

        noaa, nasa, jma = read_global_t(DATALOC + "{}_O.csv".format(IS_timeseries_root))
        hadsst = read_hadcrut_crutem(DATALOC+"hadsst3_new_logo.1981-2010.csv")

        p0 = ax5.plot(noaa.times, noaa.data, ls='-', c=COLOURS[noaa.name], label=noaa.name, lw=LW)
        p1 = ax5.plot(nasa.times, nasa.data, ls='-', c=COLOURS[nasa.name], label=nasa.name, lw=LW)
   #     p2 = ax5.plot(jma.times, jma.data, ls='-', c=COLOURS[jma.name], label=jma.name, lw=LW)
        p3 = ax5.plot(hadsst.times, hadsst.data, ls='-', c=COLOURS[hadsst.name], label=hadsst.name, lw=LW)
        ax5.fill_between(hadsst.times, hadsst.lower, hadsst.upper, \
                             where=hadsst.upper > hadsst.lower, color='0.5', alpha=0.7)
        p4 = ax1.fill(np.NaN, np.NaN, '0.5', alpha=0.7)

        ax5.axhline(0, c='0.5', ls='--')

        ax5.legend([p0[0], p1[0], (p3[0], p4[0])], [noaa.name, nasa.name, hadsst.name], \
                       loc=LEGEND_LOC, ncol=2, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, \
                       labelspacing=0.1, columnspacing=0.5, bbox_to_anchor=BBOX)

        ax5.text(0.02, 0.9, "(e) In Situ Ocean only", transform=ax5.transAxes, fontsize=settings.FONTSIZE)

        utils.thicken_panel_border(ax5)
    #    ax5.yaxis.set_ticks_position('left')

        #*******************
        # reanalysis O

        merra = utils.read_merra(settings.REANALYSISLOC + "MERRA-2_SfcAnom{}.dat".format(settings.YEAR), "temperature", "O")
        jra_actual, jra_anoms = utils.read_jra55(settings.REANALYSISLOC + "JRA-55_tmp2m_globalocean_ts.txt", "temperature")
        twenty_cr_actuals = utils.read_20cr(settings.REANALYSISLOC + "airsktocean.txt", "temperature")
        dummy, twenty_cr_anoms = utils.calculate_climatology_and_anomalies_1d(twenty_cr_actuals, 1981, 2010)
        twenty_cr_anoms.zorder=-1
 
        # 2018 no MERRA
        utils.plot_ts_panel(ax6, [jra_anoms, ocean_era5_anoms], "-", "temperature", loc=LEGEND_LOC, bbox=BBOX)

        ax6.text(0.02, 0.9, "(f) Reanalysis Ocean only", transform=ax6.transAxes, fontsize=settings.FONTSIZE)

        #*******************
        # prettify

        fig.text(0.03, 0.5, "Anomalies ("+r'$^{\circ}$'+"C)", va='center', rotation='vertical', fontsize=settings.FONTSIZE)


        plt.xlim([1900, int(settings.YEAR)+2])

        minorLocator = MultipleLocator(5)
        for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
            ax.set_ylim(YLIM)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(settings.FONTSIZE)
            ax.xaxis.set_minor_locator(minorLocator)

        for tick in ax6.xaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE)

        fig.subplots_adjust(right=0.96, top=0.995, bottom=0.02, hspace=0.001)

        plt.savefig(settings.IMAGELOC+"SAT_ts{}".format(settings.OUTFMT))

        plt.close()

    #************************************************************************
    # ERA5 Anomaly figure

    if True:
        # Read in ERA anomalies

        cube_list = iris.load(settings.REANALYSISLOC + "era5_t2m_{}01-{}12_ann_ano.nc".format(settings.YEAR, settings.YEAR))
        for cube in cube_list:
            if cube.var_name == "T2M":
                break

        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()

        bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "SAT_{}_anoms_era5".format(settings.YEAR), cube[0], settings.COLOURMAP_DICT["temperature"], bounds, "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", title="ERA5")

    #************************************************************************
    # MERRA2 Anomaly figure
    if True:
        cube_list = iris.load(settings.REANALYSISLOC + "MERRA-2_SfcAnom_{}.nc".format(settings.YEAR))
        for cube in cube_list:
            if cube.var_name == "t2ma":
                break

        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()

        bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "SAT_{}_anoms_merra".format(settings.YEAR), cube[0], \
                                       settings.COLOURMAP_DICT["temperature"], bounds, \
                                       "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", title="MERRA-2")


    #************************************************************************
    # HadCRUT4 Anomaly figure
    if True:
        cube_list = iris.load(DATALOC + "HadCRUT.4.6.0.0.median.nc")
        for cube in cube_list:
            if cube.var_name == "temperature_anomaly":
                break

        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()

        # restrict to 1851 to last full year
        date_constraint = utils.periodConstraint(cube, dt.datetime(1850, 1, 1), dt.datetime(int(settings.YEAR)+1, 1, 1))
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

        bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "SAT_{}_anoms_hadcrut4".format(settings.YEAR), annual_cube, \
                                       settings.COLOURMAP_DICT["temperature"], bounds, \
                                       "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", title="HadCRUT 4.6")

    #************************************************************************
    # NOAA data Anomaly figure - incl plate 2.1
    if True:
        cube = read_noaa_mlost(DATALOC + "mlost-box.ytd.12.1981-2010bp.txt", int(settings.YEAR))

        bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_SAT_{}_anoms_noaa".format(settings.YEAR), cube, \
                                       settings.COLOURMAP_DICT["temperature"], bounds, \
                                       "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", \
                                       figtext="(a) Surface Temperature", \
                                       save_netcdf_filename="{}MLOST_for_NOAA_{}.nc".format(DATALOC, dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y")))

        utils.plot_smooth_map_iris(settings.IMAGELOC + "SAT_{}_anoms_noaa".format(settings.YEAR), cube, \
                                       settings.COLOURMAP_DICT["temperature"], bounds, \
                                       "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", title="NOAAGlobalTemp")

    #************************************************************************
    # NASA GISS Anomaly figure
    if True:
        #cube = read_nasa_giss(DATALOC + "nasa-2015-anomalies-wrt1981-2010bp")

        cube = iris.load(DATALOC + "gistemp1200_GHCNv4_ERSSTv5.nc")[0]

        # convert to 1981-2010 climatology.
        clim_constraint = utils.periodConstraint(cube, dt.datetime(1981, 1, 1), dt.datetime(2011, 1, 1))
        clim_cube = cube.extract(clim_constraint)

        clim_data = clim_cube.data.reshape(-1, 12, clim_cube.data.shape[-2], clim_cube.data.shape[-1])

        # more than 15 years present
        climatology = np.ma.mean(clim_data, axis=0)
        nyears = np.ma.count(clim_data, axis=0)
        climatology = np.ma.masked_where(nyears <= 15, climatology) # Kate keeps GT 15.

        # extract final year
        final_year_constraint = utils.periodConstraint(cube, dt.datetime(int(settings.YEAR), 1, 1), \
                                                           dt.datetime(int(settings.YEAR)+1, 1, 1))
        final_year_cube = cube.extract(final_year_constraint)

        final_year_cube.data = final_year_cube.data - climatology

        # more than 6 months present
        annual_cube = final_year_cube.collapsed(['time'], iris.analysis.MEAN)
        nmonths = np.ma.count(final_year_cube.data, axis=0)
        annual_cube.data = np.ma.masked_where(nmonths <= 6, annual_cube.data)

        bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]

        utils.plot_smooth_map_iris(settings.IMAGELOC + "SAT_{}_anoms_nasa".format(settings.YEAR), annual_cube, \
                                       settings.COLOURMAP_DICT["temperature"], bounds, \
                                       "Anomalies from 1981-2010 ("+r'$^{\circ}$'+"C)", title="NASA GISS")

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()


#************************************************************************
#                                 END
#************************************************************************
