#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for Phenology (PHEN) section.
#       For BAMS SotC 2016
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
import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
import matplotlib.path as mpath
from matplotlib.ticker import MultipleLocator
import matplotlib.image as mpimg

import iris
import cartopy

import utils # RJHD utilities
import settings


data_loc = "{}/{}/data/PHEN/".format(settings.ROOTLOC, settings.YEAR)
reanalysis_loc = "{}/{}/data/RNL/".format(settings.ROOTLOC, settings.YEAR)
image_loc = "{}/{}/images/".format(settings.ROOTLOC, settings.YEAR)

LEGEND_LOC = 'lower left'


#************************************************************************
def read_uk_csv(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read
    :returns: Timeseries object s
    '''

    raw_data = np.genfromtxt(filename, dtype=(str), skip_header=6, skip_footer=13, delimiter=",")
    
    alder = utils.Timeseries("A. glutinosa", raw_data[:, 0].astype(int), raw_data[:, 1].astype(int))
    chestnut = utils.Timeseries("A. hippocastanum", raw_data[:, 0].astype(int), raw_data[:, 2].astype(int))
    oak = utils.Timeseries("Q. robur (2000-{})".format(settings.YEAR[2:]), raw_data[:, 0].astype(int), raw_data[:, 3].astype(int))
    beech = utils.Timeseries("F. sylvatica", raw_data[:, 0].astype(int), raw_data[:, 4].astype(int))
    
    return alder, chestnut, oak, beech # read_uk_csv

#************************************************************************
def read_uk_oak_csv(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read
    :returns: Timeseries object s
    '''

    raw_data = np.genfromtxt(filename, dtype=(str), skip_header=1, delimiter=",")
    
    for rd in raw_data:
        for loc in (1, 2):
            if rd[loc] == "": rd[loc] = "-99"

    indata = raw_data[:, 2].astype(float)
    indata = np.ma.masked_where(indata == -99, indata)
    
    oak = utils.Timeseries("Q. robur (1951-99)", raw_data[:, 0].astype(int), indata)
    
    return oak  # read_uk_oak_csv
#************************************************************************
def read_uk_map_csv(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read
    :returns: Timeseries object s
    '''

    raw_data = np.genfromtxt(filename, dtype=(str), skip_header=1, delimiter=",")
    
    species = raw_data[:, 1]
    day = raw_data[:, 4].astype(int)
    latitude = raw_data[:, -2].astype(float)
    longitude = raw_data[:, -1].astype(float)

    return species, day, latitude, longitude # read_uk_map_csv

#************************************************************************
def read_lake_csv(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read
    :returns: Timeseries object s
    '''

    raw_data = np.genfromtxt(filename, dtype=(str), skip_header=1, delimiter=",")
    
    chlorophyll = utils.Timeseries("Chlorophyll-a", raw_data[: 0].astype(int), raw_data[:, 2].astype(int))
    
    return chlorophyll # read_lake_csv

   
#************************************************************************
def read_dwd_betula(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read
    :returns: Timeseries object s
    '''

    raw_data = np.genfromtxt(filename, dtype=(int), skip_header=3, delimiter=",", skip_footer=7)

    betula = utils.Timeseries("B. pendula", raw_data[:, 0], raw_data[:, 1])
    
    return betula # read_dwd_betula

#************************************************************************
def read_dwd_quercus(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.
    
    :param str filename: file to read
    :returns: Timeseries object s
    '''
    
    raw_data = np.genfromtxt(filename, dtype=(int), skip_header=3, delimiter=",", skip_footer=7, filling_values=-1)

    raw_data = np.ma.masked_where(raw_data < 0, raw_data)

    quercus_out = utils.Timeseries("Q. robur", raw_data[:, 0], raw_data[:, 8])
    quercus_fall = utils.Timeseries("Q. robur", raw_data[:, 0], raw_data[:, 22])
    
    return quercus_out, quercus_fall # read_dwd_betula
   
#************************************************************************
def read_dwd_fagus(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read
    :returns: Timeseries object s
    '''

    raw_data = np.genfromtxt(filename, dtype=(int), skip_header=1, delimiter=",")
    
    betula_fall = utils.Timeseries("B. pendula", raw_data[:, 0], raw_data[:, 2])
    fagus_out = utils.Timeseries("F. sylvatica", raw_data[:, 0], raw_data[:, 4])
    fagus_fall = utils.Timeseries("F. sylvatica", raw_data[:, 0], raw_data[:, 6])
    
    return betula_fall, fagus_out, fagus_fall # read_dwd_fagus

#************************************************************************
def read_bartlett_green(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read
    :returns: Timeseries object s
    '''

    raw_data = np.genfromtxt(filename, dtype=(int))

    start = utils.Timeseries("Greenup", raw_data[:, 0], raw_data[:, 1])
    end = utils.Timeseries("Greendown", raw_data[:, 0], raw_data[:, 2])
    diff = utils.Timeseries("Difference ", raw_data[:, 0], raw_data[:, 3])

    return start, end, diff # read_bartlett_green

#************************************************************************
def read_bartlett_gpp(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read
    :returns: Timeseries object s
    '''

    raw_data = np.genfromtxt(filename, dtype=(float), skip_header=1, delimiter=",")

    gpp = utils.Timeseries("Bartlett GPP", raw_data[:, 0], raw_data[:, 6])
    gcc = utils.Timeseries("Bartlett GCC", raw_data[:, 7], raw_data[:, 8])

    return gpp, gcc # read_bartlett_gpp

#************************************************************************
def read_bartlett_duke(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read
    :returns: Timeseries object s
    '''

    raw_data = np.genfromtxt(filename, dtype=(float), skip_header=1, delimiter=",")

    days = raw_data[:, 0]  

    bartlett = utils.Timeseries("Bartlett GCC", days, raw_data[:, 10])
    duke = utils.Timeseries("Duke GCC", days, raw_data[:, 12])

    return bartlett, duke # read_bartlett_duke

#************************************************************************
def read_modis_ts(filename):
    '''
    Read the timeseries data, and returning Timeseries objects.

    :param str filename: file to read
    :returns: Timeseries object s
    '''

    raw_data = np.genfromtxt(filename, dtype=(float), skip_header=1, delimiter=",")

    years = raw_data[:, 0]  

    # 2017 entries
#    sos = utils.Timeseries("SOS", years, raw_data[:, 1])
#    eos = utils.Timeseries("EOS", years, raw_data[:, 2])
#    lai = utils.Timeseries("LAI", years, raw_data[:, 3])
    
    # 2018 entries
    sos_nh = utils.Timeseries("SOS", years, raw_data[:, 1])
    sos_na = utils.Timeseries("SOS", years, raw_data[:, 2])
    sos_ea = utils.Timeseries("SOS", years, raw_data[:, 3])

    eos_nh = utils.Timeseries("EOS", years, raw_data[:, 5])
    eos_na = utils.Timeseries("EOS", years, raw_data[:, 6])
    eos_ea = utils.Timeseries("EOS", years, raw_data[:, 7])

    sprt_nh = utils.Timeseries("Spring T", years, raw_data[:, 9])
    sprt_na = utils.Timeseries("Spring T", years, raw_data[:, 10])
    sprt_ea = utils.Timeseries("Spring T", years, raw_data[:, 11])

    falt_nh = utils.Timeseries("Autumn T", years, raw_data[:, 13])
    falt_na = utils.Timeseries("Autumn T", years, raw_data[:, 14])
    falt_ea = utils.Timeseries("Autumn T", years, raw_data[:, 15])

    return sos_na, sos_ea, sprt_na, sprt_ea # read_modis_ts

#************************************************************************
def run_all_plots():


    # #***********************
    # # Three panel timeseries - UK, DWD, Windermere
    


    # chlorophyll = read_lake_csv(os.path.join(data_loc, "timeseries", "Windermere.csv"))
    # betula_out = read_dwd_betula(os.path.join(data_loc, "timeseries", "dwd-data_per180112_Betula_pendula_leaf_out.csv"))
    # quercus_out, quercus_fall = read_dwd_quercus(os.path.join(data_loc, "timeseries", "dwd-data_per-180112_Quercus robur leave unf.csv"))
    # betula_fall, fagus_out, fagus_fall = read_dwd_fagus(os.path.join(data_loc, "timeseries", "dwd-data_per180112_Fagus sylvatica and Betula pendula leaf data.csv"))

    # alder, chestnut, oak, beech = read_uk_csv(os.path.join(data_loc, "timeseries", "details_for_climate_report.csv"))
    # long_oak = read_uk_oak_csv(os.path.join(data_loc, "timeseries", "NatureCalendar_pedunculate_oak_first_leaf_1753_1999.csv"))

    # fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 7), sharex=True)

    # # DWD data - budburst
    # utils.plot_ts_panel(ax1, [betula_out, fagus_out, quercus_out], "-", "phenological", loc="")

    # lines, labels = ax1.get_legend_handles_labels()
    # newlabels = []
    # for ll in labels:
    #     newlabels += [r'${}$'.format(ll)]
    # ax1.legend(lines, newlabels, loc="lower left", ncol=2, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, labelspacing=0.1, columnspacing=0.5)
    # ax4 = plt.axes([0.84, 0.84, 0.15, 0.15])  
    # img = mpimg.imread(os.path.join(data_loc, 'DWD_oak_leaf.jpg'))
    # ax4.imshow(img)
    # ax4.set_xticks([])
    # ax4.set_yticks([])
    

    # # UK data - budburst
    # utils.plot_ts_panel(ax2, [alder, chestnut, beech, long_oak, oak], "-", "phenological", loc="")
    # ax2.plot(long_oak.times, long_oak.data, ls=None, marker="o", c=settings.COLOURS["phenological"][long_oak.name], markeredgecolor=settings.COLOURS["phenological"][long_oak.name])
    # lines, labels = ax2.get_legend_handles_labels()
    # newlabels = []
    # for ll in labels:
    #     newlabels += [r'${}$'.format(ll)]
    # ax2.legend(lines, newlabels, loc="lower left", ncol=2, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, labelspacing=0.1, columnspacing=0.5)
    # ax5 = plt.axes([0.84, 0.54, 0.15, 0.15])  
    # img = mpimg.imread(os.path.join(data_loc, '61d1c433-b65c-47bc-b831-7e44ab74a895.jpg'))
    # ax5.imshow(img)
    # ax5.set_xticks([])
    # ax5.set_yticks([])

    # # Windermere
    # utils.plot_ts_panel(ax3, [chlorophyll], "-", "phenological", loc="lower left")
    # ax6 = plt.axes([0.84, 0.23, 0.15, 0.15])  
    # img = mpimg.imread(os.path.join(data_loc, 'Asterionella and aulacoseira 05 May 2017.jpg'))
    # ax6.imshow(img)
    # ax6.set_xticks([])
    # ax6.set_yticks([])

    # # prettify
    # ax1.set_ylim([70, 160])
    # ax2.set_ylim([70, 160])
    # ax3.set_ylim([70, 160])

    # ax2.set_ylabel("Day of year", fontsize=settings.FONTSIZE)

    # ax1.text(0.02, 0.88, "(a) Germany - tree leaf unfurling", transform=ax1.transAxes, fontsize=settings.FONTSIZE)
    # ax2.text(0.02, 0.88, "(b) UK - tree bud burst & leaf out", transform=ax2.transAxes, fontsize=settings.FONTSIZE)
    # ax3.text(0.02, 0.88, "(c) UK - Windermere - plankton", transform=ax3.transAxes, fontsize=settings.FONTSIZE)

    # for tick in ax3.xaxis.get_major_ticks():
    #     tick.label.set_fontsize(settings.FONTSIZE) 

    # minorLocator = MultipleLocator(100)
    # for ax in [ax1, ax2, ax3]:
    #     utils.thicken_panel_border(ax)
    #     ax.set_yticks(ax.get_yticks()[1:-1])
    #     ax.xaxis.set_minor_locator(minorLocator)
    #     for tick in ax.yaxis.get_major_ticks():
    #         tick.label.set_fontsize(settings.FONTSIZE) 

    # ax1.set_xlim([1950, 2020])

    # plt.setp([a.get_xticklabels() for a in [ax1, ax2]], visible=False)
    # fig.subplots_adjust(right=0.95, top=0.95, bottom=0.05, hspace=0.001)

    # plt.savefig(image_loc+"PHEN_tree_ts{}".format(settings.OUTFMT))
    # plt.close()


    # #***********************
    # # Barlett start/end timeseries - image inset
    # #   Hand-copied data from "Bartlett Gcc Transition Dates 2008-2017.xlsx"

    # # up, down, diff = read_bartlett_green(os.path.join(data_loc, "timeseries", "Bartlett_greenupdown.dat"))

    # # fig = plt.figure(figsize = (8, 9))
    # # plt.clf()

    # # # make axes by hand
    # # ax1 = plt.axes([0.15,0.3,0.8,0.65])
    # # ax2 = plt.axes([0.15,0.1,0.8,0.2], sharex = ax1)
    
    # # print "need inset axes for images"

    # # # plot data
    # # ax1.plot(up.times, up.data, c = "g", ls = "-", lw = 2, label = "Green up")
    # # ax1.plot(down.times, down.data, c = "brown", ls = "-", lw = 2, label = "Green down")
    # # ax1.fill_between(up.times, up.data, down.data, color = "0.75")

    # # ax2.plot(diff.times, diff.data, c = "c", ls = "-", lw = 2)

    # # ax1.legend(loc="upper right", ncol = 1, frameon = False, prop = {'size':settings.LEGEND_FONTSIZE}, labelspacing = 0.1, columnspacing = 0.5)

    # # # sort axes
    # # ax1.set_xlim([2008, 2018])
    # # ax2.set_ylim([140,180])
    # # ax1.text(0.02, 0.93, "(a) Bartlett seasons", transform = ax1.transAxes, fontsize = settings.FONTSIZE)
    # # ax2.text(0.02, 0.85, "(b) Difference", transform = ax2.transAxes, fontsize = settings.FONTSIZE)

    # # # prettify
    # # ax1.set_ylabel("Day of Year", fontsize = settings.FONTSIZE)
    # # ax2.set_ylabel("Days", fontsize = settings.FONTSIZE)
    # # minorLocator = MultipleLocator(100)
    # # for ax in [ax1, ax2]:
    # #     utils.thicken_panel_border(ax)
    # #     ax.set_yticks(ax.get_yticks()[1:-1])
    # #     ax.xaxis.set_minor_locator(minorLocator)
    # #     for tick in ax.yaxis.get_major_ticks():
    # #         tick.label.set_fontsize(settings.FONTSIZE) 
    # #     for tick in ax.xaxis.get_major_ticks():
    # #         tick.label.set_fontsize(settings.FONTSIZE) 

    # # plt.setp(ax1.get_xticklabels(), visible=False)

    # # plt.savefig(image_loc+"PHEN_bartlett_ts{}".format(settings.OUTFMT))
    # # plt.close()

    # #***********************
    # # Barlett start/end timeseries - image inset
    # #   Hand-copied data from "Bartlett Gcc Transition Dates 2008-2017.xlsx"
    # #   https://matplotlib.org/examples/pylab_examples/broken_axis.html

    # up, down, diff = read_bartlett_green(os.path.join(data_loc, "timeseries", "Bartlett_greenupdown.dat"))

    # fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 9), sharex=True)

    # print("need inset axes for images")

    # # plot same on both
    # ax1.plot(up.times, up.data, c="g", ls="-", lw=2, label="Green up", marker="o")
    # ax1.plot(down.times, down.data, c="brown", ls="-", lw=2, label="Green down", marker="o")
    # ax1.fill_between(up.times, up.data, down.data, color="lightgreen")

    # ax2.plot(up.times, up.data, c="g", ls="-", lw=2, marker="o")
    # ax2.plot(down.times, down.data, c="brown", ls="-", lw=2, marker="o")
    # ax2.fill_between(up.times, up.data, down.data, color="lightgreen")

    # ax1.invert_yaxis()
    # ax2.invert_yaxis()

    # # hide the spines between ax and ax2
    # ax1.spines['bottom'].set_visible(False)
    # ax2.spines['top'].set_visible(False)
    # ax1.xaxis.tick_top()
    # ax1.tick_params(labeltop='off')  # don't put tick labels at the top
    # ax2.xaxis.tick_bottom()

    # # and the difference
    # ax3.plot(diff.times, diff.data, c="c", ls="-", lw=2, marker="o")

    # # plot the legend
    # ax1.legend(loc="upper right", ncol=1, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, labelspacing=0.1, columnspacing=0.5)
    
    # # prettify
    # ax1.set_xlim([2007, int(settings.YEAR)+1])
    # ax2.set_ylabel("Day of Year", fontsize=settings.FONTSIZE)
    # ax3.set_ylabel("Days", fontsize=settings.FONTSIZE)
    # ax1.text(0.02, 0.88, "(a) Bartlett seasons", transform=ax1.transAxes, fontsize=settings.FONTSIZE)
    # ax3.text(0.02, 0.85, "(b) Difference", transform=ax3.transAxes, fontsize=settings.FONTSIZE)

    # # zoom in
    # ax1.set_ylim([150, 100])
    # ax2.set_ylim([300, 250])
    # ax3.set_ylim([136, 179])

    # minorLocator = MultipleLocator(100)
    # for ax in [ax1, ax2, ax3]:
    #     utils.thicken_panel_border(ax)
    #     ax.set_yticks(ax.get_yticks()[1:-1])
    #     ax.xaxis.set_minor_locator(minorLocator)
    #     for tick in ax.yaxis.get_major_ticks():
    #         tick.label.set_fontsize(settings.FONTSIZE) 
    #     for tick in ax.xaxis.get_major_ticks():
    #         tick.label.set_fontsize(settings.FONTSIZE) 

    # plt.setp(ax1.get_xticklabels(), visible=False)
    # fig.subplots_adjust(right=0.95, top=0.95, hspace=0.001)


    # plt.savefig(image_loc+"PHEN_bartlett_ts_gap{}".format(settings.OUTFMT))
    # plt.close()

    # #***********************
    # # Bartlett GPP & Duke separate

    # fig = plt.figure(figsize=(8, 6.5))
    # plt.clf()
    # ax1 = plt.axes([0.12, 0.52, 0.78, 0.44])
    # ax2 = ax1.twinx()
    # ax3 = plt.axes([0.12, 0.08, 0.78, 0.44], sharex=ax1)
    
    # gpp, gcc = read_bartlett_gpp(os.path.join(data_loc, "timeseries", "Bartlett GPP Flux and Canopy Greenness.csv"))


    # bartlett, duke = read_bartlett_duke(os.path.join(data_loc, "timeseries", "Bartlett 2008-2017 and Duke 2017 Camera Greenness.csv"))

    # utils.plot_ts_panel(ax1, [gcc], "-", "phenological", loc="")
    # utils.plot_ts_panel(ax2, [gpp], "-", "phenological", loc="")
    # utils.plot_ts_panel(ax3, [bartlett, duke], "-", "phenological", loc="upper right", ncol=1)
    # ax3.axhline(0.36, c='0.5', ls='--')
    # ax3.axvline(185, c='0.5', ls=":")
    
    
    # # fix the legend
    # lines1, labels1 = ax1.get_legend_handles_labels()
    # lines2, labels2 = ax2.get_legend_handles_labels()
    # ax2.legend(lines1 + lines2, labels1 + labels2, loc="upper right", ncol=1, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, labelspacing=0.1, columnspacing=0.5)

    # # prettify
    # ax2.yaxis.set_label_position("right")
    # ax2.yaxis.set_ticks_position('right')

    # ax1.set_ylim([0.34, 0.48])
    # ax3.set_ylim([0.34, 0.48])
    # ax2.set_ylim([-2, 12])

    # ax1.set_ylabel("Canopy Greeness", fontsize=settings.FONTSIZE)
    # ax2.set_ylabel("Gross Primary Productivity", fontsize=settings.FONTSIZE)
    # ax3.set_ylabel("Canopy Greeness", fontsize=settings.FONTSIZE)

    # ax1.text(0.02, 0.88, "(c)", transform=ax1.transAxes, fontsize=settings.FONTSIZE)
    # ax3.text(0.02, 0.85, "(d)", transform=ax3.transAxes, fontsize=settings.FONTSIZE)

    # for tick in ax3.xaxis.get_major_ticks():
    #     tick.label.set_fontsize(settings.FONTSIZE) 

    # minorLocator = MultipleLocator(100)
    # for ax in [ax1, ax3]:
    #     utils.thicken_panel_border(ax)
    #     ax.set_yticks(ax.get_yticks()[1:-1])
    #     ax.xaxis.set_minor_locator(minorLocator)
    #     for tick in ax.yaxis.get_major_ticks():
    #         tick.label.set_fontsize(settings.FONTSIZE) 

    # utils.thicken_panel_border(ax2)
    # ax2.set_yticks(ax2.get_yticks()[1:-1])
    # for tick in ax2.yaxis.get_major_ticks():
    #     tick.label2.set_fontsize(settings.FONTSIZE)
    
    # plt.setp(ax1.get_xticklabels(), visible=False)
    # plt.setp(ax2.get_xticklabels(), visible=False)
    # fig.subplots_adjust(right=0.95, top=0.95, hspace=0.001)

    # ax4 = plt.axes([0.14, 0.65, 0.15, 0.15])  
    # img = mpimg.imread(os.path.join(data_loc, 'bbc7_2017_01_20_140005.jpg'))
    # ax4.imshow(img)
    # ax4.set_xticks([])
    # ax4.set_yticks([])
    # ax4.set_title("January")
    # ax5 = plt.axes([0.72, 0.65, 0.15, 0.15])  
    # img = mpimg.imread(os.path.join(data_loc, 'bbc7_2017_07_17_130005.jpg'))
    # ax5.imshow(img)
    # ax5.set_xticks([])
    # ax5.set_yticks([])
    # ax5.set_title("July")

    # ax6 = plt.axes([0.14, 0.25, 0.15, 0.15])  
    # img = mpimg.imread(os.path.join(data_loc, 'bbc7_2017_07_04_124505.jpg'))
    # ax6.imshow(img)
    # ax6.set_xticks([])
    # ax6.set_yticks([])
    # ax6.set_title("Bartlett")
    # ax7 = plt.axes([0.72, 0.25, 0.15, 0.15])  
    # img = mpimg.imread(os.path.join(data_loc, 'dukehw_2017_07_04_120110.jpg'))
    # ax7.imshow(img)
    # ax7.set_xticks([])
    # ax7.set_yticks([])
    # ax7.set_title("Duke")


    # plt.savefig(image_loc+"PHEN_Bartlett_GPP_Duke_separate{}".format(settings.OUTFMT))
    # plt.close()

    # #***********************
    # # Bartlett GPP & Duke combined

    # # fig = plt.figure(figsize = (8, 5))
    # # plt.clf()
    # # ax1 = plt.axes([0.12, 0.1, 0.78, 0.8])
    # # ax2 = ax1.twinx()
    
    # # gpp, gcc = read_bartlett_gpp(os.path.join(data_loc, "timeseries", "Bartlett GPP Flux and Canopy Greenness.csv"))

    # # bartlett, duke = read_bartlett_duke(os.path.join(data_loc, "timeseries", "Bartlett 2008-2017 and Duke 2017 Camera Greenness.csv"))

    # # utils.plot_ts_panel(ax1, [gcc, duke], "-", "phenological", loc = "")
    # # utils.plot_ts_panel(ax2, [gpp], "-", "phenological", loc = "")
    
    # # # fix the legend
    # # lines1, labels1 = ax1.get_legend_handles_labels()
    # # lines2, labels2 = ax2.get_legend_handles_labels()
    # # # have to remove duplicate label
    # # ax2.legend(lines1 + lines2, labels1 + labels2, loc="upper right", ncol = 1, frameon = False, prop = {'size':settings.LEGEND_FONTSIZE}, labelspacing = 0.1, columnspacing = 0.5)

    # # # prettify
    # # ax2.yaxis.set_label_position("right")

    # # ax1.set_ylim([0.34, 0.48])
    # # ax2.set_ylim([-2, 12])

    # # ax1.set_ylabel("Canopy Greeness", fontsize = settings.FONTSIZE)
    # # ax2.set_ylabel("Gross Primary Productivity", fontsize = settings.FONTSIZE)
 
    # # for tick in ax1.xaxis.get_major_ticks():
    # #     tick.label.set_fontsize(settings.FONTSIZE) 

    # # minorLocator = MultipleLocator(100)
    # # for ax in [ax1, ax3]:
    # #     utils.thicken_panel_border(ax)
    # #     ax.set_yticks(ax.get_yticks()[1:])
    # #     ax.xaxis.set_minor_locator(minorLocator)
    # #     for tick in ax.yaxis.get_major_ticks():
    # #         tick.label.set_fontsize(settings.FONTSIZE) 

    # # utils.thicken_panel_border(ax2)
    # # ax2.set_yticks(ax2.get_yticks()[1:-1])
    # # for tick in ax2.yaxis.get_major_ticks():
    # #     tick.label2.set_fontsize(settings.FONTSIZE)
    
    # # fig.subplots_adjust(right = 0.95, top = 0.95, hspace = 0.001)

    # # plt.savefig(image_loc+"PHEN_Bartlett_GPP_Duke_combined{}".format(settings.OUTFMT))
    # # plt.close()

    # #***********************
    # # UK Map
    # import cartopy.feature as cfeature
    # land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
    #                                         edgecolor='face',
    #                                         facecolor=cfeature.COLORS['land'])
 
    # species, days, lats, lons = read_uk_map_csv(os.path.join(data_loc, "timeseries", "UK_Leafout_4Trees.csv"))
    
    # cmap = plt.cm.YlGn
    # bounds = np.arange(80, 160, 10)
    # norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)

    # fig = plt.figure(figsize=(8, 10.5))
    # plt.clf()
    # ax = plt.axes([0.05, 0.05, 0.9, 0.9], projection=cartopy.crs.LambertConformal(central_longitude=-7.5))

    # ax.gridlines() #draw_labels=True)
    # ax.add_feature(land_50m, zorder=0, facecolor="0.9", edgecolor="k")
    # ax.set_extent([-10, 4, 48, 60], cartopy.crs.PlateCarree())

    # scat = ax.scatter(lons, lats, c=days, transform=cartopy.crs.PlateCarree(), cmap=cmap, norm=norm, s=25, edgecolor='0.5', linewidth='0.5', zorder=10)

    # utils.thicken_panel_border(ax)

    # cb = plt.colorbar(scat, orientation='horizontal', ticks=bounds[1:-1], label="Day of year", drawedges=True, fraction=0.1, pad=0.05, aspect=15, shrink=0.8)

    # # prettify
    # cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
    # cb.outline.set_linewidth(2)
    # cb.dividers.set_color('k')
    # cb.dividers.set_linewidth(2)

    # plt.savefig(image_loc + "PHEN_UK_map{}".format(settings.OUTFMT))
    # plt.close()


    #***********************
    # MODIS

    cubelist = iris.load(os.path.join(data_loc, "MODIS.CMG.{}.SOS.EOS.Anomaly.nc".format(settings.YEAR)))

    for c, cube in enumerate(cubelist):
        if cube.name() == "EOS":
            eos_cube = cubelist[c]
        elif cube.name() == "SOS":
            sos_cube = cubelist[c]
        elif cube.name() == "MAX":
            max_cube = cubelist[c]
       
    # deal with NANS
#    eos_cube.data = np.ma.masked_where(eos_cube.data != eos_cube.data, eos_cube.data)
    sos_cube.data = np.ma.masked_where(sos_cube.data != sos_cube.data, sos_cube.data)
#    max_cube.data = np.ma.masked_where(max_cube.data != max_cube.data, max_cube.data)


            
    # set up a 1 x 2 set of axes (2018 only needs one panel)
    fig = plt.figure(figsize=(8, 8))
    plt.clf()

    # set up plot settings
    BOUNDS = [[-100, -20, -10, -5, -2, 0, 2, 5, 10, 20, 100]]#, [-100, -20, -10, -5, -2, 0, 2, 5, 10, 20, 100]]
    CMAPS = [settings.COLOURMAP_DICT["phenological_r"]]#, settings.COLOURMAP_DICT["phenological"]]

    CUBES = [sos_cube]# eos_cube]
    LABELS = ["(a) Start of Season (SOS)"]#, "(b) End of Season (EOS)"]

    # boundary circle
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    # spin through axes
    for a in range(1):  

        ax = plt.subplot(1, 1, a+1, projection=cartopy.crs.NorthPolarStereo())

        plot_cube = CUBES[a]

        if settings.OUTFMT in [".eps", ".pdf"]:
            if plot_cube.coord("latitude").points.shape[0] > 90 or plot_cube.coord("longitude").points.shape[0] > 360:
                regrid_size = 1.0
                print("Regridding cube for {} output to {} degree resolution".format(settings.OUTFMT, regrid_size))
                print("Old Shape {}".format(plot_cube.data.shape))
                plot_cube = utils.regrid_cube(plot_cube, regrid_size, regrid_size)
                print("New Shape {}".format(plot_cube.data.shape))

        ax.gridlines() #draw_labels=True)
        ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
        ax.coastlines()
        ax.set_boundary(circle, transform=ax.transAxes)
        ax.set_extent([-180, 180, 45, 90], cartopy.crs.PlateCarree())

        ext = ax.get_extent() # save the original extent

        cmap = CMAPS[a]
        norm = mpl.cm.colors.BoundaryNorm(BOUNDS[a], cmap.N)
        mesh = iris.plot.pcolormesh(plot_cube, cmap=cmap, norm=norm, axes=ax)

        ax.set_extent(ext, ax.projection) # fix the extent change from colormesh
        ax.text(-0.1, 1.0, LABELS[a], fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)

        cb = plt.colorbar(mesh, orientation='horizontal', ticks=BOUNDS[a][1:-1], label="Anomaly (days)", drawedges=True, fraction=0.1, pad=0.05, aspect=15, shrink=0.8)
        # prettify
        cb.set_ticklabels(["{:g}".format(b) for b in BOUNDS[a][1:-1]])
        cb.outline.set_linewidth(2)
        cb.dividers.set_color('k')
        cb.dividers.set_linewidth(2)

        ax.set_extent([-180, 180, 45, 90], cartopy.crs.PlateCarree())

    fig.subplots_adjust(bottom=0.05, top=0.95, left=0.04, right=0.95, wspace=0.02)

    plt.title("")

    plt.savefig(image_loc + "PHEN_modis_polar{}".format(settings.OUTFMT))
    plt.close()

    del eos_cube
    del sos_cube
    del cubelist

    #***********************
    # MODIS LAI
    # max_cube.data = np.ma.masked_where(max_cube.data <= -9000, max_cube.data)
    # if settings.OUTFMT in [".eps", ".pdf"]:
    #     if max_cube.coord("latitude").points.shape[0] > 180 or max_cube.coord("longitude").points.shape[0] > 360:
    #         regrid_size = 1.0
    #         print("Regridding cube for {} output to {} degree resolution".format(settings.OUTFMT, regrid_size))
    #         print("Old Shape {}".format(max_cube.data.shape))
    #         max_cube = utils.regrid_cube(max_cube, regrid_size, regrid_size)
    #         print("New Shape {}".format(max_cube.data.shape))

    # fig = plt.figure(figsize=(8, 8))
    # plt.clf()
    # ax = plt.axes([0.05, 0.05, 0.9, 0.9], projection=cartopy.crs.NorthPolarStereo())

    # ax.gridlines() #draw_labels=True)
    # ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
    # ax.coastlines()
    # ax.set_boundary(circle, transform=ax.transAxes)
    # ax.set_extent([-180, 180, 45, 90], cartopy.crs.PlateCarree())

    # ext = ax.get_extent() # save the original extent

    # cmap = settings.COLOURMAP_DICT["phenological_r"]
    # bounds = [-100, -5, -3, -2, -1, 0, 1, 2, 3, 5, 100]
    # norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)
    # mesh = iris.plot.pcolormesh(max_cube, cmap=cmap, norm=norm, axes=ax)

    # ax.set_extent(ext, ax.projection) # fix the extent change from colormesh
    # ax.text(-0.1, 1.0, "Max Leaf Area Index", fontsize=settings.FONTSIZE * 0.8, transform=ax.transAxes)

    # cb = plt.colorbar(mesh, orientation='horizontal', ticks=bounds[1:-1], label="Anomaly (m"+r'$^2'+"/m"+r'$^2'+")", drawedges=True, fraction=0.1, pad=0.05, aspect=15, shrink=0.8)
    # # prettify
    # cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
    # cb.outline.set_linewidth(2)
    # cb.dividers.set_color('k')
    # cb.dividers.set_linewidth(2)

    # ax.set_extent([-180, 180, 45, 90], cartopy.crs.PlateCarree())

    # plt.savefig(image_loc + "PHEN_modis_lai{}".format(settings.OUTFMT))
    # plt.close()

    #***********************
    # MODIS timeseries - 2017

    # sos, eos, lai = read_modis_ts(os.path.join(data_loc, "MODIS.CMG.{}.SOS.EOS.MAX.AnomalyTS.csv".format(settings.YEAR)))


    # plt.clf()

    # fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 4), sharex=True)
    
    # # plot same on both
    # ax1.plot(sos.times, sos.data, c="g", ls="-", lw=2, label="Start of Season", marker="o")
    # ax1.plot(eos.times, eos.data, c="brown", ls="-", lw=2, label="End of Season", marker="o")
    # ax1.fill_between(eos.times, eos.data, sos.data, color="lightgreen")

    # ax2.plot(sos.times, sos.data, c="g", ls="-", lw=2, marker="o")
    # ax2.plot(eos.times, eos.data, c="brown", ls="-", lw=2, marker="o")
    # ax2.fill_between(eos.times, eos.data, sos.data, color="lightgreen")

    # # invert so green down is on the bottom
    # ax1.invert_yaxis()
    # ax2.invert_yaxis()

    # # hide the spines between ax and ax2
    # ax1.spines['bottom'].set_visible(False)
    # ax2.spines['top'].set_visible(False)
    # ax1.xaxis.tick_top()
    # ax1.tick_params(labeltop='off')  # don't put tick labels at the top
    # ax2.xaxis.tick_bottom()

    # # plot the legend
    # ax1.legend(loc="upper right", ncol=1, frameon=False, prop={'size':settings.LEGEND_FONTSIZE}, labelspacing=0.1, columnspacing=0.5)
    
    # # prettify
    # ax1.set_xlim([1999, 2018])
    # ax2.set_ylabel("Day of Year", fontsize=settings.FONTSIZE)

    # # zoom in
    # ax1.set_ylim([160, 110])
    # ax2.set_ylim([300, 250])

    # # set label
    # ax1.text(-0.1, 0.92, "(c)", transform=ax1.transAxes, fontsize=settings.FONTSIZE)

    # minorLocator = MultipleLocator(100)
    # for ax in [ax1, ax2]:
    #     utils.thicken_panel_border(ax)
    #     ax.set_yticks(ax.get_yticks()[1:-1])
    #     ax.xaxis.set_minor_locator(minorLocator)
    #     for tick in ax.yaxis.get_major_ticks():
    #         tick.label.set_fontsize(settings.FONTSIZE) 
    #     for tick in ax.xaxis.get_major_ticks():
    #         tick.label.set_fontsize(settings.FONTSIZE) 

    # plt.setp(ax1.get_xticklabels(), visible=False)
    # fig.subplots_adjust(right=0.95, top=0.95, hspace=0.001)

    # plt.savefig(image_loc+"PHEN_modis_ts{}".format(settings.OUTFMT))
    # plt.close()

    #***********************
    # MODIS timeseries - 2018
    
    sos_na, sos_ea, sprt_na, sprt_ea = read_modis_ts(os.path.join(data_loc, "MODIS.CMG.{}.SOS.EOS.SPRT.FALT.TS.csv".format(settings.YEAR)))

    dummy, sos_na = utils.calculate_climatology_and_anomalies_1d(sos_na, 2000, 2010)
    dummy, sos_ea = utils.calculate_climatology_and_anomalies_1d(sos_ea, 2000, 2010)

    fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 6.5), sharex=True)

    # North America
    # use un-anomalised spring T to get legend without the data
    utils.plot_ts_panel(ax1, [sos_na, sprt_na], "-", "phenological", loc=LEGEND_LOC)   

    # Eurasia
    utils.plot_ts_panel(ax2, [sos_ea, sprt_ea], "-", "phenological", loc=LEGEND_LOC)

    # make twin axes for Spring T
    dummy, sprt_na = utils.calculate_climatology_and_anomalies_1d(sprt_na, 2000, 2010)
    dummy, sprt_ea = utils.calculate_climatology_and_anomalies_1d(sprt_ea, 2000, 2010)

    ax3 = ax1.twinx()
    utils.plot_ts_panel(ax3, [sprt_na], "-", "phenological", loc="")   

    ax4 = ax2.twinx()
    utils.plot_ts_panel(ax4, [sprt_ea], "-", "phenological", loc="")   

    # prettify
    ax1.set_ylim([-8, 8])
    ax2.set_ylim([-8, 8])
    ax3.set_ylim([2, -2])
    ax4.set_ylim([2, -2])

    # labels
    ax1.text(0.02, 0.88, "(a) North America", transform=ax1.transAxes, fontsize=settings.FONTSIZE)
    ax2.text(0.02, 0.88, "(b) Eurasia", transform=ax2.transAxes, fontsize=settings.FONTSIZE)

    ax1.text(0.47, 0.88, "2018 SOS Anomaly = 1.86 days", transform=ax1.transAxes, fontsize=settings.FONTSIZE*0.8)
    ax1.text(0.47, 0.78, "2018 Spring T anomaly = -0.75 "+r'$^{\circ}$'+"C", transform=ax1.transAxes, fontsize=settings.FONTSIZE*0.8)
    ax2.text(0.47, 0.88, "2018 SOS Anomaly = 2.01 days", transform=ax2.transAxes, fontsize=settings.FONTSIZE*0.8)
    ax2.text(0.47, 0.78, "2018 Spring T anomaly = 0.13 "+r'$^{\circ}$'+"C", transform=ax2.transAxes, fontsize=settings.FONTSIZE*0.8)

    fig.text(0.01, 0.5, "SOS Anomaly (days)", va='center', rotation='vertical', fontsize = settings.FONTSIZE)
    fig.text(0.95, 0.5, "Temperature Anomaly ("+r'$^{\circ}$'+"C)", va='center', rotation='vertical', fontsize = settings.FONTSIZE)

    # ticks and labels
    for tick in ax2.xaxis.get_major_ticks():
        tick.label.set_fontsize(settings.FONTSIZE) 

    minorLocator = MultipleLocator(1)
    for ax in [ax1, ax2]:
        utils.thicken_panel_border(ax)
        ax.set_yticks(ax.get_yticks()[1:-1])
        ax.xaxis.set_minor_locator(minorLocator)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(settings.FONTSIZE) 

    for ax in [ax3, ax4]:
        ax.yaxis.tick_right()
        utils.thicken_panel_border(ax)
        ax.set_yticks(ax.get_yticks()[1:-1])
        ax.xaxis.set_minor_locator(minorLocator)
        for tick in ax.yaxis.get_major_ticks():
            tick.label2.set_fontsize(settings.FONTSIZE) 

    # final settings
    ax1.set_xlim([1999, 2020])

    plt.setp([a.get_xticklabels() for a in [ax1]], visible=False)
    fig.subplots_adjust(left=0.1, right=0.85, top=0.95, bottom=0.05, hspace=0.001)

    plt.savefig(image_loc+"PHEN_modis_ts{}".format(settings.OUTFMT))
    plt.close()



    return # run_all_plots


#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
