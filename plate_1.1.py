#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for Plate 1.1 .
#       For BAMS SotC 2016
#
#************************************************************************
#                    SVN Info
# $Rev:: 29                                       $:  Revision of last commit
# $Author:: rdunn                                 $:  Author of last commit
# $Date:: 2020-08-05 12:12:39 +0100 (Wed, 05 Aug #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
from __future__ import absolute_import
from __future__ import print_function

import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
import iris

import utils # RJHD utilities
import settings


DATALOC = "{}/{}/data/".format(settings.ROOTLOC, settings.YEAR)
image_loc = "{}/{}/images/".format(settings.ROOTLOC, settings.YEAR)

LW = 2

print("Standardise individual plotting scripts \nThen use main() and import statements so \ncan call the read routines only once")

insitu = "k"
satellite = "r"
reanalyses = "b"


#************************************************************************
def make_plot(ax, ts_list, color, plot_label="", ylabel="", ylim=[], ls="-", scatter=False):
    """
    Plot the lines and other details on the axis.

    """ 
    

    assert isinstance(ts_list, list)

    zorder = 1
    if color == "k": zorder = 10
    if color == "r": zorder = 5

    for ts in ts_list:
    
        if scatter:
            ax.plot(ts.times, ts.data, color=color, ls="", marker=".", zorder=zorder)
        else:
            ax.plot(ts.times, ts.data, color=color, ls=ls, lw=LW, zorder=zorder)

    if ylabel != "":
        ax.set_ylabel(ylabel)
    if ylim != []:
        ax.set_ylim(ylim)
    if plot_label != "":
        ax.text(0.02, 0.75, plot_label, transform=ax.transAxes)

    # autogenerate text for caption
    if color == insitu:
        print("In Situ: {}".format(len(ts_list)))
    elif color == satellite:
        print("Satellite: {}".format(len(ts_list)))
    elif color == reanalyses:
        print("Reanalyses: {}".format(len(ts_list)))

    return # make_plot

  
#************************************************************************
# have separate read script for each section.
'''
black - in situ
blue - reanalyses
red - satellite
'''

#************************************************************************
def read_sie(filename):

    indata = np.genfromtxt(filename, dtype=(float), skip_header=3)
     
    arctic_max = utils.Timeseries("SIE", indata[:, 0], indata[:, 2])
    arctic_min = utils.Timeseries("SIE", indata[:, 0], indata[:, 4])
    antarctic_min = utils.Timeseries("SIE", indata[:, 0], indata[:, 6])
    antarctic_max = utils.Timeseries("SIE", indata[:, 0], indata[:, 8])

    return arctic_max, arctic_min, antarctic_min, antarctic_max # read_sie

#************************************************************************
def read_ohc(filename):

    indata = np.genfromtxt(filename, dtype=(float), skip_header=1)
     
    hadley = utils.Timeseries("OHC", indata[:, 0], indata[:, 1])
    csiro = utils.Timeseries("OHC", indata[:, 0], indata[:, 3])
    pmel = utils.Timeseries("OHC", indata[:, 0], indata[:, 5])
    ncei = utils.Timeseries("OHC", indata[:, 0], indata[:, 7])
    mri = utils.Timeseries("OHC", indata[:, 0], indata[:, 9])
    iap = utils.Timeseries("OHC", indata[:, 0], indata[:, 9])

    # reapply updated climatology
    dummy, hadley = utils.calculate_climatology_and_anomalies_1d(hadley, 1993, 2016)
    dummy, csiro = utils.calculate_climatology_and_anomalies_1d(csiro, 1993, 2016)
    dummy, pmel = utils.calculate_climatology_and_anomalies_1d(pmel, 1993, 2016)
    dummy, ncei = utils.calculate_climatology_and_anomalies_1d(ncei, 1993, 2016)
    dummy, mri = utils.calculate_climatology_and_anomalies_1d(mri, 1993, 2016)
    dummy, iap = utils.calculate_climatology_and_anomalies_1d(iap, 1993, 2016)

    return hadley, csiro, pmel, ncei, mri, iap # read_ohc

#************************************************************************
def read_slr(filename):

    indata = np.genfromtxt(filename, dtype=(float))
     
    slr = utils.Timeseries("SLR", indata[:, 0], indata[:, 1])

    return slr # read_SLR

#************************************************************************
def read_swv(filename):

    indata = np.genfromtxt(filename, dtype=(float), delimiter=",")
     
    swv = utils.Timeseries("SWV", indata[:, 0], indata[:, 1])

    return swv # read_SWV

#************************************************************************
def read_arct(filename):
    """ Read the Arctic Temperatures"""

    indata = np.genfromtxt(filename, delimiter=",", skip_header=1)

    arct = utils.Timeseries("ARCT", indata[:, 0], indata[:, 1])

    return arct # read_arct


#************************************************************************
def toYearFraction(date):
    import time
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

#************************************************************************
def annual_from_monthly(indata):

    try:
        times = indata.times.reshape(-1, 12)
        monthly = indata.data.reshape(-1, 12)
        annuals = np.ma.mean(monthly, axis=1)

        return utils.Timeseries(indata.name, times[:, 0], annuals) # annual_from_monthly
    except:

        annuals = []
        times = []
        year = int(indata.times[0])
        
        while year < int(settings.YEAR)+1:
            times += [year]

            locs, = np.where(np.logical_and(indata.times >= year, indata.times < year + 1))

            annuals += [np.ma.mean(indata.data[locs])]

            year += 1

        return utils.Timeseries(indata.name, times, annuals) # annual_from_monthly






#************************************************************************
#************************************************************************
# 3 col 11 row plot

WIDTH = 0.80/3. # leave 0.2 for axes labels
HEIGHT = 0.95/11.

W_OFFSET = 0.2/3 * 0.9 # start 90% of way through spare.
H_OFFSET = 0.04

# 
plt.figure(figsize = (16, 20))

# col 1
ax_a = plt.axes([W_OFFSET, H_OFFSET + (10 * HEIGHT), WIDTH, HEIGHT])
ax_b = plt.axes([W_OFFSET, H_OFFSET + (9 * HEIGHT), WIDTH, HEIGHT], sharex=ax_a)
ax_c = plt.axes([W_OFFSET, H_OFFSET + (8 * HEIGHT), WIDTH, HEIGHT], sharex=ax_a)
ax_d = plt.axes([W_OFFSET, H_OFFSET + (7 * HEIGHT), WIDTH, HEIGHT], sharex=ax_a)
ax_e = plt.axes([W_OFFSET, H_OFFSET + (6 * HEIGHT), WIDTH, HEIGHT], sharex=ax_a)
ax_f = plt.axes([W_OFFSET, H_OFFSET + (5 * HEIGHT), WIDTH, HEIGHT], sharex=ax_a)
ax_g = plt.axes([W_OFFSET, H_OFFSET + (4 * HEIGHT), WIDTH, HEIGHT], sharex=ax_a)
ax_h = plt.axes([W_OFFSET, H_OFFSET + (3 * HEIGHT), WIDTH, HEIGHT], sharex=ax_a)
ax_i = plt.axes([W_OFFSET, H_OFFSET + (2 * HEIGHT), WIDTH, HEIGHT], sharex=ax_a)
ax_j = plt.axes([W_OFFSET, H_OFFSET + (1 * HEIGHT), WIDTH, HEIGHT], sharex=ax_a)
ax_k = plt.axes([W_OFFSET, H_OFFSET + (0 * HEIGHT), WIDTH, HEIGHT], sharex=ax_a)

# col 2
ax_l = plt.axes([(2 * W_OFFSET) + (1 * WIDTH), H_OFFSET + (10 * HEIGHT), WIDTH, HEIGHT])
ax_m = plt.axes([(2 * W_OFFSET) + (1 * WIDTH), H_OFFSET + (9 * HEIGHT), WIDTH, HEIGHT], sharex=ax_l)
ax_n = plt.axes([(2 * W_OFFSET) + (1 * WIDTH), H_OFFSET + (8 * HEIGHT), WIDTH, HEIGHT], sharex=ax_l)
ax_o = plt.axes([(2 * W_OFFSET) + (1 * WIDTH), H_OFFSET + (7 * HEIGHT), WIDTH, HEIGHT], sharex=ax_l)
ax_p = plt.axes([(2 * W_OFFSET) + (1 * WIDTH), H_OFFSET + (6 * HEIGHT), WIDTH, HEIGHT], sharex=ax_l)
ax_q = plt.axes([(2 * W_OFFSET) + (1 * WIDTH), H_OFFSET + (5 * HEIGHT), WIDTH, HEIGHT], sharex=ax_l)
ax_r = plt.axes([(2 * W_OFFSET) + (1 * WIDTH), H_OFFSET + (4 * HEIGHT), WIDTH, HEIGHT], sharex=ax_l)
ax_s = plt.axes([(2 * W_OFFSET) + (1 * WIDTH), H_OFFSET + (3 * HEIGHT), WIDTH, HEIGHT], sharex=ax_l)
ax_t = plt.axes([(2 * W_OFFSET) + (1 * WIDTH), H_OFFSET + (2 * HEIGHT), WIDTH, HEIGHT], sharex=ax_l)
ax_u = plt.axes([(2 * W_OFFSET) + (1 * WIDTH), H_OFFSET + (1 * HEIGHT), WIDTH, HEIGHT], sharex=ax_l)
ax_v = plt.axes([(2 * W_OFFSET) + (1 * WIDTH), H_OFFSET + (0 * HEIGHT), WIDTH, HEIGHT], sharex=ax_l)

# col 3
ax_w = plt.axes([(3 * W_OFFSET) + (2 * WIDTH), H_OFFSET + (10 * HEIGHT), WIDTH, HEIGHT])
ax_x = plt.axes([(3 * W_OFFSET) + (2 * WIDTH), H_OFFSET + (9 * HEIGHT), WIDTH, HEIGHT], sharex=ax_w)
ax_y = plt.axes([(3 * W_OFFSET) + (2 * WIDTH), H_OFFSET + (8 * HEIGHT), WIDTH, HEIGHT], sharex=ax_w)
ax_z = plt.axes([(3 * W_OFFSET) + (2 * WIDTH), H_OFFSET + (7 * HEIGHT), WIDTH, HEIGHT], sharex=ax_w)
ax_aa = plt.axes([(3 * W_OFFSET) + (2 * WIDTH), H_OFFSET + (6 * HEIGHT), WIDTH, HEIGHT], sharex = ax_w)
ax_ab = plt.axes([(3 * W_OFFSET) + (2 * WIDTH), H_OFFSET + (5 * HEIGHT), WIDTH, HEIGHT], sharex = ax_w)
ax_ac = plt.axes([(3 * W_OFFSET) + (2 * WIDTH), H_OFFSET + (4 * HEIGHT), WIDTH, HEIGHT], sharex = ax_w)
ax_ad = plt.axes([(3 * W_OFFSET) + (2 * WIDTH), H_OFFSET + (3 * HEIGHT), WIDTH, HEIGHT], sharex = ax_w)
ax_ae = plt.axes([(3 * W_OFFSET) + (2 * WIDTH), H_OFFSET + (2 * HEIGHT), WIDTH, HEIGHT], sharex = ax_w)
ax_af = plt.axes([(3 * W_OFFSET) + (2 * WIDTH), H_OFFSET + (1 * HEIGHT), WIDTH, HEIGHT], sharex = ax_w)
ax_ag = plt.axes([(3 * W_OFFSET) + (2 * WIDTH), H_OFFSET + (0 * HEIGHT), WIDTH, HEIGHT], sharex = ax_w)

#***************************
# A - NH POLAR STRATOSPHERIC
print("NH Polar Stratospheric Ozone")

import soz
nh = soz.read_ts(DATALOC + "{}/NH_polar_ozone.txt".format("SOZ"))

make_plot(ax_a, [nh], insitu, plot_label="(a) N. Hemisphere Polar Stratospheric Ozone (Mar)\n(actual values)", ylabel="DU", ylim=[340, 530])

#***************************
# B - SH POLAR STRATOSPHERIC
print("SH Polar Stratospheric Ozone")

sh = soz.read_ts(DATALOC + "{}/SH_polar_ozone.txt".format("SOZ"))

make_plot(ax_b, [sh], insitu, plot_label="(b) S. Hemisphere Polar Stratospheric Ozone (Oct)\n(actual values)", ylabel="DU", ylim=[180, 449])

#***************************
# C - APPARENT TRANSMISSION
print("Apparent Transmission")

import at
at, at24 = at.read_csv(DATALOC + "{}/mlo_trans_csv.txt".format("AT"))
annuals = annual_from_monthly(at)

make_plot(ax_c, [annuals], insitu, plot_label="(c) Apparent Transmission (Mauna Loa)\n(actual values)", ylabel="AT", ylim=[0.81,0.99])

print("Using Arctic Temperature - Jessica/Deke")

#arctic = read_arct(DATALOC + "{}/arctic-60N-temps.csv".format("PLT1_1"))

#make_plot(ax_c, [arctic], insitu, plot_label="(c) Arctic Temperature (60-90N)\n(1981-2010)", ylabel='$^{\circ}$'+"C", ylim=[-1.9, 2.4])
#***************************
# D - SURFACE TEMPERATURE
print("Surface Air Temperature")

import sat
# in situ
noaa, nasa, jma = sat.read_global_t(DATALOC + "{}/{}_LO.csv".format("SAT", sat.IS_timeseries_root))
hadcrut = sat.read_hadcrut_crutem(DATALOC+"{}/hadcrut4.1981-2010.csv".format("SAT"))

make_plot(ax_d, [noaa, nasa, jma, hadcrut], insitu)

# reanalyses
merra = utils.read_merra(settings.REANALYSISLOC + "MERRA-2_SfcAnom{}.dat".format(settings.YEAR), "temperature", "LO")
jra_actuals, jra_anoms = utils.read_jra55(settings.REANALYSISLOC + "JRA-55_tmp2m_global_ts.txt", "temperature")

era5_globe, era5_ocean, era5_land, era5tropics = utils.era5_ts_read(settings.REANALYSISLOC, "sat", annual=True)
global_era5_clim, global_era5_anoms = utils.calculate_climatology_and_anomalies_1d(era5_globe, 1981, 2010)

make_plot(ax_d, [merra, jra_anoms, era5_globe], reanalyses, plot_label="(d) Surface Temperature\n(1981-2010)", ylabel='$^{\circ}$'+"C", ylim=[-0.9, 1.2])

#***************************
# E - LOWER TROP TEMPERATURE
print("Lower Tropospheric Temperature")

import ltt

raobcore, rich, ratpac, UAH, rss, era5, jra,  merra= ltt.read_csv(DATALOC + "{}/SotC_AnnTemps_2020_0220_LTTGL.csv".format("LTT"))

# in situ
make_plot(ax_e, [raobcore, rich, ratpac], insitu)
# satellite
make_plot(ax_e, [UAH, rss], satellite)

jra_actuals, jra_anoms = utils.read_jra55(settings.REANALYSISLOC + "JRA-55_MSUch2LT_global_ts.txt", "temperature")
merra_actuals, merra_anoms = utils.read_merra_LT_LS(settings.REANALYSISLOC + "MERRA2_MSU_Tanom_ann_{}.dat".format(settings.YEAR), LT=True)

# reanalyses
make_plot(ax_e, [era5, jra_anoms, merra_anoms], reanalyses, plot_label="(e) Lower Tropospheric Temperature\n(1981-2010)", ylabel='$^{\circ}$'+"C", ylim=[-0.9, 1.2])

#***************************
# F - LOWER STRAT TEMPERATURE
print("Lower Stratospheric Temperature")

import lst

UAH, rss, ratpac, raobcore, rich, noaa, jra, merra = lst.read_csv(DATALOC + "{}/SotC_AnnTemps_2020_0220_LSTGL.csv".format("LST"))

# upper stratosphere for completeness
ssu1, ssu1_2, ssu2, ssu2_2, ssu3, ssu3_2 = lst.read_ssu(DATALOC + "{}/SSU.dat".format("LST")) 

eral, erao, eralo = lst.read_era5(DATALOC + "{}/ERA5_TLS_GLOBAL".format("LST"))
era5 = utils.Timeseries("ERA5", np.reshape(eralo.times, [-1, 12])[:,0], utils.annual_average(eralo.data))

# in situ
make_plot(ax_f, [raobcore, rich, ratpac], insitu)
# satellite
make_plot(ax_f, [UAH, noaa, rss], satellite)

jra_actuals, jra_anoms = utils.read_jra55(settings.REANALYSISLOC + "JRA-55_MSUch4_global_ts.txt", "temperature")
merra_actuals, merra_anoms = utils.read_merra_LT_LS(settings.REANALYSISLOC + "MERRA2_MSU_Tanom_ann_{}.dat".format(settings.YEAR), LS=True)

# reanalyses
make_plot(ax_f, [era5, jra_anoms, merra_anoms], reanalyses, plot_label="(f) Lower Stratospheric Temperature\n(1981-2010)", ylabel='$^{\circ}$'+"C", ylim=[-0.9, 2.7])

#***************************
# G - Extreme Warm/Cool Days
print("Temperature Extremes")

import tex

tx90p, cover = tex.obtain_timeseries(DATALOC + "{}/GHCND_{}_1951-{}_RegularGrid_global_2.5x2.5deg_LSmask.nc".format("TEX", "TX90p", int(settings.YEAR) + 1), "Ann", "GHCNDEX", "TX90p")
tx10p, cover = tex.obtain_timeseries(DATALOC + "{}/GHCND_{}_1951-{}_RegularGrid_global_2.5x2.5deg_LSmask.nc".format("TEX", "TX10p", int(settings.YEAR) + 1), "Ann", "GHCNDEX", "TX10p")

print("Axis in % days, so resetting back to original values (x3.65 applied in tex.py)")

tx90p.data = tx90p.data / 3.65
tx10p.data = tx10p.data / 3.65

# in situ
make_plot(ax_g, [tx90p], insitu)
make_plot(ax_g, [tx10p], insitu, plot_label="(g) Extremes [Warm Days and Cool days (dotted)]\n(1961-1990)", ylabel="% days", ylim=[5, 25], ls=":")

#***************************
# H - Arctic Sea Ice Extent - Jessica/Deke give data
print("Arctic Sea Ice - Jessica/Deke")

arctic_max, arctic_min, antarctic_min, antarctic_max = read_sie(DATALOC +"{}/SIE.dat".format("PLT1_1"))

make_plot(ax_h, [arctic_max], insitu)
make_plot(ax_h, [arctic_min], insitu, plot_label="(h) Arctic Sea Ice Extent [max and min (dotted)]\n1981-2010", ylabel="x 10"+r'$^6$'+" km"+r'$^2$', ylim=[-2.9, 2.9], ls=":")

#***************************
# I - Antarctic Sea Ice Extent - Jessica/Deke give data
print("Antarctic Sea Ice - Jessica/Deke")

make_plot(ax_i, [antarctic_max], insitu)
make_plot(ax_i, [antarctic_min], insitu, plot_label="(i) Antarctic Sea Ice Extent [max and min (dotted)]\n1981-2010", ylabel="x 10"+r'$^6$'+" km"+r'$^2$', ylim=[-1.2, 3.0], ls=":")

#***************************
# J - Glacier Mass Balance
print("Glacier Mass Balance")

import agl
balance, cumul_balance = agl.read_glacier(DATALOC + "{}/global_mass_balance_{}.csv".format("AGL", settings.YEAR))

# in situ
make_plot(ax_j, [cumul_balance], insitu, plot_label="(j) Glacier Cumulative Mean Specific Balance\n(actual values)", ylabel="equivalent depth\n in water (m)", ylim=[-22, 9])

#***************************
# K - Snow Cover
print("NH Snow Cover")

import snw

NH, Eurasia, NAmer = snw.read_snow(DATALOC + "{}/Robinson-snow-cover-{}.csv".format("SNW", settings.YEAR))

# satellite
make_plot(ax_k, [annual_from_monthly(NH)], satellite, plot_label="(k) Northern Hemisphere Snow Cover Extent\n(1966-{})".format(settings.YEAR), ylabel="x 10"+r'$^6$'+" km"+r'$^2$', ylim=[-1.9, 3.6])

#***************************
# L - Lower Stratospheric Water Vapour

swv = read_swv(DATALOC + "{}/Plate1_StratWV_83hPa_BLD_{}.txt".format("SWV", settings.YEAR))

# in situ
make_plot(ax_l, [swv], insitu, plot_label="(l) Lower Stratospheric Water Vapor\n(actual values)", ylabel="ppmv", ylim=[2.01, 6.99], scatter=True)


#***************************
# M - CLOUDINESS
print("Cloudiness")

import cld

patmosx, hirs, misr, modis, calipso, ceres, satcorps, clara_a2, patmosdx, cci = cld.read_ts(DATALOC + "{}/{}_global_cloudiness_timeseries.txt".format("CLD", settings.YEAR), anomaly=True)

# satellite
make_plot(ax_m, [patmosx, hirs, misr, modis, calipso, ceres, satcorps, clara_a2, patmosdx, cci], satellite, plot_label="(m) Cloudiness\n(2003-2015)", ylabel="%", ylim=[-6, 9])


#***************************
# N - Total Column Water - Land
print("Total Column Water Vapour - Land")

import tcw

#merra2_land, erai_land, era5_land, jra_land, cosmic_land, gnss_land = tcw.read_csv(DATALOC + "{}/time_series_tpw_land.txt".format("TCW"), domain="L")
merra2_land, era5_land, jra_land, cosmic_land, gnss_land=tcw.read_ncdf_ts(DATALOC + "{}/TPW_{}_anom_TS.v2.nc".format("TCW", settings.YEAR), domain="L")
gnss_land.name = "GNSS (Ground Based)"


# updated file June 2020
print("remove TCW from June 2020")
alldata = np.genfromtxt(DATALOC + "TCW/TCWV_2020_ts_updated.dat")
alldata = np.ma.masked_where(alldata == 0.0, alldata)
cosmic_land = utils.Timeseries("COSMIC RO", alldata[:, 0], alldata[:, 7])

# satellite
make_plot(ax_n, [cosmic_land], satellite)
# in situ
make_plot(ax_n, [gnss_land], insitu)
# reanalyses
make_plot(ax_n, [era5_land, jra_land, merra2_land], reanalyses, plot_label="(n) Total Column Water Vapour - Land\n(1981-2010)", ylabel="mm", ylim=[-1.2, 1.9])

#***************************
# O - Total Column Water - Ocean
print("Total Column Water Vapour - Marine")

#merra2_ocean, era5_ocean, jra_ocean, cosmic_ocean, radiometer_ocean = tcw.read_csv(DATALOC + "{}/time_series_tpw_ocean.txt".format("TCW"), domain="O")
merra2_ocean, era5_ocean, jra_ocean, cosmic_ocean, radiometer_ocean = tcw.read_ncdf_ts(DATALOC + "{}/TPW_{}_anom_TS.v2.nc".format("TCW", settings.YEAR), domain="O")
radiometer_ocean.name = "RSS Satellite"

# updated file June 2020
print("remove TCW from June 2020")
alldata = np.genfromtxt(DATALOC + "TCW/TCWV_2020_ts_updated.dat")
alldata = np.ma.masked_where(alldata == 0.0, alldata)
cosmic_ocean = utils.Timeseries("COSMIC RO", alldata[:, 0], alldata[:, 2])

# satellite
make_plot(ax_o, [radiometer_ocean, cosmic_ocean], satellite)
# reanalyses
make_plot(ax_o, [era5_ocean, jra_ocean, merra2_ocean], reanalyses, plot_label="(o) Total Column Water Vapour - Ocean\n(1981-2010)", ylabel="mm", ylim=[-1.2, 2.6])

#***************************
# P - Upper Tropospheric Humidity
print("Upper Tropospheric Humidity")

import uth

HIRSSTART = 1979
MWSTART = 1999
ERASTART = 1979
hirs = uth.read_ts(DATALOC + "{}/hirs_data.aa".format("UTH"), HIRSSTART, "HIRS", smooth=12)
mw = uth.read_ts(DATALOC + "{}/mw_data.aa".format("UTH"), MWSTART, "Microwave", smooth=12)
era5 = uth.read_ts(DATALOC + "{}/era5_data.aa".format("UTH"), ERASTART, "ERA5", smooth=3)

# satellite
make_plot(ax_p, [annual_from_monthly(hirs), annual_from_monthly(mw)], satellite,)
# reanalyses
make_plot(ax_p, [annual_from_monthly(era5)], reanalyses, plot_label="(p) Upper Tropospheric Humidity\n(2001-2010)", ylabel="% rh", ylim=[-0.59, 0.89])

#***************************
# Q - Specific Humidity - Land
print("Specific Humidity - Land")

import hum

hadisdhLQ, hadcruhLQ, hadcruhextLQ, daiLQ, eraiLQ, era5LQ, merraLQ, jraLQ, era5_mskLQ, merra_mskLQ, cr20LQ = hum.read_ts(DATALOC + "{}/HUM_timeseries_ALL{}.txt".format("HUM", settings.YEAR), "q", "L")

# in situ
make_plot(ax_q, [hadisdhLQ] , insitu)
# reanalyses
make_plot(ax_q, [ era5LQ, merraLQ, jraLQ, cr20LQ], reanalyses, plot_label="(q) Specific Humidity - Land\n(1979-2003)", ylabel="g kg"+r'$^{-1}$', ylim=[-0.5, 0.8])

#***************************
# R - Specific Humidity - Ocean
print("Specific Humidity - Marine")

hadisdhMQ, hadcruhMQ, daiMQ, nocsMQ, hoapsMQ, eraiMQ, era5MQ, merraMQ, jraMQ, cr20MQ = hum.read_ts(DATALOC + "{}/HUM_timeseries_ALL{}.txt".format("HUM", settings.YEAR), "q", "M")

# in situ
make_plot(ax_r, [hadisdhMQ, nocsMQ], insitu)
# satellite
# make_plot(ax_r, [hoapsMQ], satellite) # not for 2019
# reanalyses
make_plot(ax_r, [era5MQ, merraMQ, jraMQ, cr20MQ], reanalyses, plot_label="(r) Specific Humidity - Ocean\n(1979-2003)", ylabel="g kg"+r'$^{-1}$', ylim=[-0.5,0.8])


#***************************
# S - Relative Humidity - Land
print("Relative Humidity - Land")

hadisdhLR, hadcruhLR, hadcruhextLR, daiLR, eraiLR, era5LR, merraLR, jraLR, era5_mskLR, merra_mskLR, cr20LR = hum.read_ts(DATALOC + "{}/HUM_timeseries_ALL{}.txt".format("HUM", settings.YEAR), "rh", "L")

# in situ
make_plot(ax_s, [hadisdhLR, hadcruhLR, hadcruhextLR, daiLR,], insitu)
# reanalyses
make_plot(ax_s, [era5LR, jraLR, cr20LR], reanalyses, plot_label="(s) Relative Humidity - Land\n(1979-2003)", ylabel="% rh", ylim=[-1.5, 2.5])


#***************************
# T - Relative Humidity - Ocean
print("Relative Humidity - Marine")

hadisdhMR, hadcruhMR, daiMR, nocsMR, hoapsMR, eraiMR, era5MR, merraMR, jraMR, cr20MR = hum.read_ts(DATALOC + "{}/HUM_timeseries_ALL{}.txt".format("HUM", settings.YEAR), "rh", "M")

# in situ
make_plot(ax_t, [hadisdhMR], insitu)
# reanalyses
make_plot(ax_t, [era5MR, jraMR, cr20MR], reanalyses, plot_label="(t) Relative Humidity - Ocean\n(1979-2003)", ylabel="% rh", ylim=[-0.7, 1.5])


#***************************
# U - Precipitation - Land
print("Precipitation - Land")

import pcp
ghcn, gpcc, gpcp = pcp.read_land(DATALOC + "{}/Land_insitu_timeseries-1979.dat".format("PCP"))

# in situ
make_plot(ax_u, [ghcn, gpcc, gpcp], insitu, plot_label="(u) Precipitation - Land\n(1981-2010)", ylabel="u")
# reanalyses - not for 2019
# make_plot(ax_u, [erai, merra], reanalyses, plot_label="(u) Precipitation - Land\n(1981-2010)", ylabel="u)

#***************************
# V - Precipitation - Ocean
gpcp = pcp.read_ocean(DATALOC + "{}/Ocean_insitu_timeseries-1979.dat".format("PCP"))

# in situ
make_plot(ax_v, [gpcp], insitu, plot_label="(v) Precipitation - Ocean\n(1981-2010)", ylabel="v")

#***************************
# V - Southern Oscillation Index
print("Southern Oscillation Index")

import slp

soi = slp.read_soi(DATALOC + "{}/soiplaintext.html".format("SLP"))

#make_plot(ax_v, [soi], insitu, plot_label="(v) Southern Oscillation Index\n ", ylabel="Standard Units", ylim=[-40, 49])

#***************************
# W - OHC - Jessica & Deke provide
print("Ocean Heat- Jessica/Deke")

hadley, csiro, pmel, ncei, mri, iap = read_ohc(DATALOC + "{}/OHC.dat".format("PLT1_1"))

make_plot(ax_w, [hadley, csiro, pmel, ncei, mri, iap], insitu, plot_label="(w) Ocean Heat Content (0-700m)\n(1983-{})".format(settings.YEAR), ylabel='$10^{21}$'+"J", ylim=[-140, 140])

#***************************
# X - Sea Level Rise - Jessica & Deke provide
print("Sea Level Rise - Jessica/Deke")

slr = read_slr(DATALOC + "{}/SLR.dat".format("PLT1_1"))

make_plot(ax_x, [slr], insitu, plot_label="(x) Sea Level Rise\n(actual values)", ylabel="mm", ylim=[-21, 120])

#***************************
# Y - Tropospheric Ozone
print("Tropospheric Ozone - need actuals")

import tco

tco = tco.read_data(DATALOC + "{}/BAMS_SOTC_TROPOSPHERIC_OZONE_TG_60Sto60N_{}.txt".format("TCO", settings.YEAR), "TCO")

# satellite
make_plot(ax_y, [annual_from_monthly(tco)], satellite, plot_label="(y) Tropospheric Ozone\n(actual values)", ylabel="Ozone Burden (Tg)", ylim=[280, 320])


#***************************
# Z - Tropospheric Wind Speed
print("Tropospheric Wind Speed")

import uaw

era5, erai, merra, jra55 = uaw.read_uaw_ts(DATALOC + "{}/Globe850.nc".format("UAW"), smooth=True)

# sonde
# make_plot(ax_z, [annual_from_monthly(grasp)], insitu)
# reanalyses
make_plot(ax_z, [annual_from_monthly(erai), annual_from_monthly(era5), annual_from_monthly(merra), annual_from_monthly(jra55)], reanalyses, plot_label="(z) Tropospheric Wind Speed at 850hPa\n(1981-2010)", ylabel="m s"+r'$^{-1}$', ylim=[-0.5, 0.7])

#***************************
# AA - LAND WIND SPEED
print("Near Surface Wind Speed - Land")

import wnd

# in situ
Globe = wnd.Region("Globe (excl Austr)", "GlobalNoOz", "black")
years, anomalies, m3, m10 = wnd.read_hadisd_annual_anomalies(Globe)
lwnd = utils.Timeseries("HadISD", years, anomalies)

make_plot(ax_aa, [lwnd], insitu, plot_label="(aa) Land Wind Speed\n(1981-2010)", ylabel="m s"+r'$^{-1}$', ylim=[-0.29, 0.39])

#***************************
# AB - OCEAN WIND SPEED
print("Near Surface Wind Speed - Ocean")

ownd = wnd.read_ts_cube(DATALOC + "{}/rss_wind_trend_anomaly_SOTC_{}.nc".format("WND", settings.YEAR), "RSS_wind_global_annual_anom_ts", "Satellite MW Radiometers")
ocean_obs_clim, ocean_obs_anoms = utils.calculate_climatology_and_anomalies_1d(ownd, 1981, 2010)

make_plot(ax_ab, [ocean_obs_anoms], satellite)

# reanalyses
jra_actuals, jra_anoms = utils.read_jra55(settings.REANALYSISLOC + "JRA-55_ws10m_globalocean_ts.txt", "wind")
merra_anoms = utils.read_merra(settings.REANALYSISLOC + "MERRA-2_SfcAnom{}.dat".format(settings.YEAR), "wind", "O", anomalies=True)
era5_globe, era5_ocean, era5_land, era5tropics = utils.era5_ts_read(settings.REANALYSISLOC, "wnd", annual=True)
ocean_era5_clim, ocean_era5_anoms = utils.calculate_climatology_and_anomalies_1d(era5_ocean, 1981, 2010)
twenty_cr_actuals = utils.read_20cr(settings.REANALYSISLOC + "wspd10m.ocean.txt", "wind speed")
ocean_20cr_clim, ocean_20cr_anoms = utils.calculate_climatology_and_anomalies_1d(twenty_cr_actuals, 1981, 2010)

make_plot(ax_ab, [ocean_era5_anoms, merra_anoms, ocean_20cr_anoms], reanalyses, plot_label="(ab) Ocean Wind Speed\n(1981-2010)", ylabel="m s"+r'$^{-1}$', ylim=[-0.29,0.49])

#***************************
# AC - Biomass Burning
print("Biomass Burning")

import bob

gfed = bob.read_gfed_csv(DATALOC + "{}/data4Johannes.txt".format("BOB"), "global", make_annual=True)
gfas = bob.read_csv(DATALOC + "{}/timeseries_glob".format("BOB"), "GFAS1p4", make_annual=True)

make_plot(ax_ac, [gfed, gfas], satellite, plot_label="(ac) Biomass Burning\n(actual values)", ylabel="Pg C yr"+r'$^{-1}$', ylim=[1.25, 3.7])


#***************************
# AD - Soil Moisture
print("Soil Moisture")

import sms

cube_list = np.array(iris.load(DATALOC + "{}/ESA_CCI_SM_COMBINED_monthAnomaliesPerHemisphere.nc".format("SMS")))
names = np.array([c.var_name for c in cube_list])
glob = cube_list[names == "Anomalies_global"][0]
years = sms.convert_times(glob)

annuals = annual_from_monthly(utils.Timeseries("SMS", years, glob.data))
print(annuals)
make_plot(ax_ad, [annuals], satellite, plot_label="(ad) Soil Moisture\n(1991-2010)", ylabel="m"+r'$^{3}$', ylim=[-0.009, 0.009])


#***************************
# AE - Terrestrial Water Storage
print("Terrestrial Water Storage")

import tws

grace = tws.read_ts(DATALOC + "{}/avg_JPLM06v2_land.txt".format("TWS"), "GRACE")
# NOTE, model used for 2018 report alongside GRACE
annuals = annual_from_monthly(grace)
make_plot(ax_ae, [annuals], satellite, plot_label="(ae) Terrestrial Water Storage\n(2005-2010)", ylabel="equivalent depth\nin water (cm)", ylim=[-2.3, 1.4])

#***************************
# AF - FAPAR
print("FAPAR")

import fpr

data = fpr.read_binary_ts(DATALOC + "{}/TimeSeries_faparanomaliesglobal_bams_v2018_C6_2020.bin".format("FPR"))

for dataset in data:
    if dataset.name == "Globe":
        
        annuals = annual_from_monthly(dataset)
        make_plot(ax_af, [annuals], satellite, plot_label="(af) FAPAR\n(1998-{})".format(settings.YEAR), ylabel="FAPAR", ylim=[-0.009, 0.019])


#***************************
# AG - Land Surface Albedo
print("Land Surface Albedo")

import abd

IRdata = abd.read_binary_ts(DATALOC + "{}/TimeseriesBHRNIR_C7_poids_{}.bin".format("ABD", int(settings.YEAR)+1))
Vdata = abd.read_binary_ts(DATALOC + "{}/TimeseriesBHRV_C7_poids_{}.bin".format("ABD", int(settings.YEAR)+1))

for dataset in Vdata:
    if dataset.name == "Globe":
        annuals = annual_from_monthly(dataset)
        make_plot(ax_ag, [annuals], satellite)

for dataset in IRdata:
    if dataset.name == "Globe":
        annuals = annual_from_monthly(dataset)
        make_plot(ax_ag, [annuals], satellite, plot_label="(ag) Land Surface Albedo - visible & infrared (dotted)\n(2003-{})".format(settings.YEAR), ylabel="%", ylim=[-4, 5], ls=":")


#************************************************************************
# tidy up
all_axes = [ax_a, ax_b, ax_c, ax_d, ax_e, ax_f, ax_g, ax_h, ax_i, ax_j, ax_k, ax_l, ax_m, ax_n, ax_o, ax_p, ax_q, ax_r, ax_s, ax_t, ax_u, ax_v, ax_w, ax_x, ax_y, ax_z, ax_aa, ax_ab, ax_ac, ax_ad, ax_ae, ax_af, ax_ag]

for ax in all_axes:
    ax.axhline(0, ls="--", color="0.5")

ax_a.set_xlim([1950, int(settings.YEAR)+2])
ax_l.set_xlim([1960, int(settings.YEAR)+2])
ax_w.set_xlim([1980, int(settings.YEAR)+2])

# remove x-tick labels on all but lowest 3.
plt.setp([a.get_xticklabels() for a in [ax_a, ax_b, ax_c, ax_d, ax_e, ax_f, ax_g, ax_h, ax_i, ax_j, ax_l, ax_m, ax_n, ax_o, ax_p, ax_q, ax_r, ax_s, ax_t, ax_u, ax_w, ax_x, ax_y, ax_z, ax_aa, ax_ab, ax_ac, ax_ad, ax_ae, ax_af]], visible=False)



plt.savefig(settings.IMAGELOC + "plate_1_1{}".format(settings.OUTFMT))

plt.close()
#************************************************************************
#                                 END
#************************************************************************
