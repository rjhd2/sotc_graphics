#!/usr/bin/env python
#************************************************************************
#
#  Global Settings for SotC plots
#  RJHD, Exeter 2017
#
#************************************************************************
#                    SVN Info
# $Rev:: 31                                         $:  Revision of last commit
# $Author:: rdunn                                      $:  Author of last commit
# $Date:: 2021-09-06 09:52:46 +0100 (Mon, 06 Sep #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import matplotlib.pyplot as plt
import numpy as np

import os
import configparser

#************************************************************************
def adjust_RdYlBu():
    '''
    Original RdYlBu has a greenish colour on the blue end, which doesn't contrast
    well with the yellow.  As this is a diverging colourmap, want to have a larger
    distinction for the positives and negatives.

    Take the colourmap, move the blues along a bit, and repeat at the end.
    '''

    print("for 2017/18 - change to match the purple adjustment")

    cmap = plt.cm.RdYlBu

    cmaplist = [cmap(i) for i in range(cmap.N)]

    # ignore first 20 of the blues
    N = 20
    retained_colours = np.array(cmaplist[255//2 + N :])

    stretched = np.ones((len(cmaplist[255//2:]), 4))

    # interpolate out
    for i in range(stretched.shape[1]):
        stretched[:, i] = np.interp(np.linspace(0, retained_colours.shape[0], stretched.shape[0]), \
                                        list(range(retained_colours.shape[0])), retained_colours[:, i])

    # convert to tuple
    new_colours = []
    for c in stretched:
        new_colours += [tuple(c)]

    # copy back
    cmaplist[255//2:] = new_colours

    # enforce white for the centre of the colour range
    for i in [126, 127, 128, 129]:
        cmaplist[i] = (1.0, 1.0, 1.0, 1.0)

    # and make the colour map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    cmap_r = cmap.from_list('Custom cmap', cmaplist[::-1], cmap.N)

    return cmap, cmap_r # adjust_RdYlBu

#************************************************************************
def adjust_PuOr():
    '''
    Original PuOr has the lowest purple a very similar shade to the grey
    of the continents.  As this is a diverging colourmap, want to have a larger
    distinction for the positives and negatives.

    Take the colourmap, move the purples along a bit, and repeat at the end.
    '''

    cmap = plt.cm.PuOr

    cmaplist = [cmap(i) for i in range(cmap.N)]

    # ignore first 20 of the purples
    N = 20    
    retained_colours = np.array(cmaplist[255//2 + N :])

    stretched = np.ones((len(cmaplist[255//2:]), 4))

    # interpolate out
    for i in range(stretched.shape[1]):
        stretched[:, i] = np.interp(np.linspace(0, retained_colours.shape[0], stretched.shape[0]), \
                                    list(range(retained_colours.shape[0])), retained_colours[:, i])

    # convert to tuple
    new_colours = []
    for c in stretched:
        new_colours += [tuple(c)]

    # copy back
    cmaplist[255//2:] = new_colours

    # enforce white for the centre of the colour range
    for i in [126, 127, 128, 129]:
        cmaplist[i] = (1.0, 1.0, 1.0, 1.0)

    # and make the colour map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    cmap_r = cmap.from_list('Custom cmap', cmaplist[::-1], cmap.N)

    return cmap, cmap_r # adjust_RdYlBu

#************************************************************************
def make_BrBu():
    '''
    Use the BrBG and the RdBu colourmaps to create a BrBu map in style of Brewer colours
    '''


    cmap = plt.cm.BrBG
    cmaplist = [cmap(i) for i in range(cmap.N)]

    cmap_R = plt.cm.RdBu
    cmaplist_R = [cmap_R(i) for i in range(cmap_R.N)]

    cmaplist[255//2:] = cmaplist_R[255//2:]

    # enforce white for the centre of the colour range
    for i in [126, 127, 128, 129]:
        cmaplist[i] = (1.0, 1.0, 1.0, 1.0)

    # and make the colour map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    cmap_r = cmap.from_list('Custom cmap', cmaplist[::-1], cmap.N)

    return cmap, cmap_r # make_BrBu

#************************************************************************
def make_BrBG():
    '''
    Enforce white at centre of BrBG (not the case using the default)
    '''


    cmap = plt.cm.BrBG
    cmaplist = [cmap(i) for i in range(cmap.N)]

    # enforce white for the centre of the colour range
    for i in [126, 127, 128, 129]:
        cmaplist[i] = (1.0, 1.0, 1.0, 1.0)

    # and make the colour map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    cmap_r = cmap.from_list('Custom cmap', cmaplist[::-1], cmap.N)

    return cmap, cmap_r # make_BrBu


#************************************************************************
#************************************************************************
CONFIG_FILE = os.path.join(os.getcwd(), "configuration.txt")

if not os.path.exists(CONFIG_FILE):
    print("Configuration file missing - {}".format(CONFIG_FILE))
    sys.exit

# read in configuration file
config = configparser.ConfigParser()
config.read(CONFIG_FILE)

# main settings
ROOTLOC = config.get("Paths", "rootloc")
YEAR = config.get("Misc", "year")
OUTFMT = config.get("Format", "outfmt")
FONTSIZE = config.getint("Format", "fontsize")
ERA5LOC_TEX = config.get("Paths", "tex_era5loc")

# derived settings
LEGEND_FONTSIZE = 0.8 * FONTSIZE
LABEL_FONTSIZE = 0.9 * FONTSIZE

IMAGELOC = "{}/{}/images/".format(ROOTLOC, YEAR)
REANALYSISLOC = "{}/{}/data/RNL/".format(ROOTLOC, YEAR)

#************************************************************************
COLOURS = {"temperature" : {"ERA-Interim" : "purple", \
                                "ERA5" : "purple", \
                                "20CRv3" : "orange", \
                                "MERRA-2" : "lime", \
                                "JRA-55" : "c", \
                                "UAH v6.0" : "b", \
                                "RSS v4.0" : "r", \
                                "NOAA v4.1" : "c", \
                                "UW" : "m", \
                                "RAOBCORE v1.7" : "r", \
                                "RICH v1.7" : "y", \
                                "RATPAC vA2" : "m", \
                                "UNSW v1.0" : "c", \
                                "SSU+AMSU" : "k", \
                                "SSU+MLS" : "r", \
                                "CFSR" : "lime", \
                                "NOAA/NCEI" : "r", \
                                "NASA/GISS" : "b", \
                                "JMA" : "c", \
                                "Berkeley" : "y", \
                                "Hadley" : "k", \
                                "HadCRUT5" : "k", \
                                "CRUTEM5" : "k", \
                                "HadSST4" : "k", \
                                "GHCNDEX" : "r", \
                                "CLASSnmat v2" : "r", \
                                "UAHNMAT v1" : "b", \
                                "SST" : "r", \
                                "Air Temperature" : "b", \
                                "CMIP5" : "k", \
                                "NOAA" : "r", \
                                "NCAR" : "k", \
                                "North" : "k", \
                                "South": "k", \
                                "QBO" : "k", \
                            "HWF" : "b", \
                            "HWM" : "r", \
                                "NOAA v4.0" : "k", \
                                "RSS v3.3" : "r", \
                            "Warmest" : "r", \
                            "Coldest" : "b", \
                            "Lake" : "k"}, \
              "lst-ssu" : {"SSU1+MLS" : "yellow", \
                           "SSU2+MLS" : "cyan", \
                           "SSU3+MLS" : "magenta", \
                           "SSU1+AMSU" : "olive", \
                           "SSU2+AMSU" : "teal", \
                           "SSU3+AMSU" : "purple", \
                           "RSS" : "r", \
                           "UAH" : "b", \
                           "NOAA" : "k"},\
               "cryosphere" : {"N Hemisphere" : 'k', \
                                   "Eurasia" : 'r', \
                                   "N America" : 'b', \
                                   "Cumulative Balance" : "k", \
                                   "Balance" : "r",\
                                   "Lake" : "0.7",\
                                   "Erie" : "b",\
                                   "Huron" : "c",\
                                   "Michigan" : "g",\
                                   "Ontario" : "lime",\
                                   "Superior" : "yellow",\
                                   "Average" : "k",\
                               "ERA5" : "purple",\
                               "In Situ" : "0.3"},\
               "hydrological" : {"ERA-Interim" : "purple", \
                                     "ERA5" : "purple", \
                                     "JRA-55" : "c", \
                                     "MERRA" : "lime", \
                                     "MERRA-2" : "lime", \
                                     "20CRv3" : "orange", \
                                     "Satellite RO" : "b", \
                                     "GNSS (Ground Based)" : "r", \
                                     "RSS Satellite" : '0.5', \
                                     "HIRS" : "k", \
                                     "MSU-HIRS" : "k", \
                                     "Microwave" : "b", \
                                     "GHCN" : "0.5", \
                                     "GPCC" : "r", \
                                     "GPCPv23" : "b", \
                                     "GHCNv2" : "k", \
                                     "GPCP" : "c", \
                                     "PATMOS-x/AVHRR" : "0.5", \
                                     "PATMOS-x/AQUA MODIS" : "c", \
                                     "PATMOS-x/AVHRR+HIRS" : "yellow", \
                                     "MISR" : "brown", \
                                     "AQUA MODIS C6" : "b", \
                                     "CALIPSO" : "r", \
                                     "CERES" : "orange", \
                                     "SatCORPS" : "lime", \
                                     "CLARA-A2": "purple", \
                                     "Cloud CCI AVHRR-PMv3": "g", \
                                     "GRACE FO" : "k", \
                                     "GRACE" : "0.5", \
                                     "Model" : "0.5", \
                                     "HadISDH" : "0.3", \
                                     "HadCRUH" : "k", \
                                     "HadCRUHExt" : "k", \
                                     "Dai" : "r", \
                                     "NOCS v2.0" : "b", \
                                     "HOAPS" : "brown", \
                                     "NCEP" : "w", \
                                     "Globe" : "k", \
                                     "N. Hemisphere" : "b", \
                                     "S. Hemisphere" : "r", \
                                     "ERA5 mask" : "purple", \
                                     "MERRA-2 mask" : "lime"},\
               "circulation" : {"Satellite MW Radiometers" : "k", \
                                    "NOCSv2.0" : "b", \
                                    "WASwind" : "r", \
                                    "ERA-Interim" : "purple", \
                                    "ERA5" : "purple", \
                                    "JRA-55" : "c", \
                                    "MERRA-2" : "lime", \
                                    "MERRA" : "lime", \
                                    "ASCAT" : "r", \
                                    "QuikSCAT" : "brown", \
                                    "20CRv3" : "orange", \
                                    "ERApreSAT" : "purple", \
                                    "CERA20C" : "purple", \
                                    "GRASP" : 'k', \
                                    "GIUB" : "lime"},\
               "radiation" : {"AT" : "r"},\
               "composition" : {"AOD monthly" : "r", \
                                "AOD annual" : "b"},\
               "land_surface" : {"Globe" : "0.7", \
                                     "N. Hemisphere" : "cornflowerblue", \
                                     "S. Hemisphere" : "lightcoral", \
                                     "Globe Smoothed" : "k", \
                                     "N. Hemisphere Smoothed" : "b", \
                                     "S. Hemisphere Smoothed" : "r", \
                                     "GFED4s" : "b", \
                                     "GFASv1.4" : "r", \
                                     "GFASv1.0" : "k"},\
               "vod" : {"Globe" : "k", \
                                     "N. Hemisphere" : "b", \
                                     "S. Hemisphere" : "r"},\
               "lst" : {"400" : "c", \
                            "300" : "m", \
                            "250" : "lime", \
                            "200" : "y", \
                            "150" : "k", \
                            "100" : "orange", \
                            "70" : "c", \
                            "50" : "m", \
                            "30" : "lime", \
                            "20" : "y", \
                            "10" : "k", \
                            "Average" : "k"},\
               "phenological" : {"Chlorophyll-a" : "g", \
                                     "B. pendula" : "b", \
                                     "Q. robur (2000-{})".format(YEAR[2:]) : "r", \
                                     "$SOS_{PO}$" : "g", \
                                     "$EOS_{PO}$" : "g", \
                                     "F. sylvatica" : "k", \
                                     "Q. robur (1951-99)" : "orange", \
                                     "A. hippocastanum" : "c", \
                                     "A. glutinosa": "m", \
                                     "Greenup" : "lime", \
                                     "Greendown" : "orange", \
                                     "Difference" : "k", \
                                     "North Basin" : "g", \
                                     "South Basin" : "m", \
                                     "Duke GCC" : "k",\
                                     "SOS" : "g",\
                                     "Spring T": "m",
                                     "EOS" : "g",\
                                     "Fall T": "m",\
                                     "$SOS_{M}$": "k",\
                                     "$EOS_{M}$": "k"}} #  note space in second Q. robur

# all maps have 10 colours - or 11 with white as central.  
# Can probably make up from Kate's code if better match needed.


RdYlBu, RdYlBu_r = adjust_RdYlBu()
PuOr, PuOr_r = adjust_PuOr()
# Composition & Land Surface - create a version of Brown-Blue brewer colours
BrBu, BrBu_r = make_BrBu()
# adjust BrBG to ensure has central colours equal to white
BrBG, BrBG_r = make_BrBG()

COLOURMAP_DICT = {"temperature" : RdYlBu_r, "temperature_r" : RdYlBu, \
                      "hydrological" : BrBG, "hydrological_r" : BrBG_r, \
                      "precip_sequential" : plt.cm.YlGnBu, "precip_sequential_r" : plt.cm.YlGnBu_r,\
                      "circulation" : PuOr, "circulation_r" : PuOr_r, \
                      "composition" : BrBu_r, "composition_r" : BrBu, \
                      "land_surface" : BrBu_r, "land_surface_r" : BrBu, "land_surface_sequential" : plt.cm.YlGn,\
                      "phenological" : BrBG, "phenological_r" : BrBG_r}



