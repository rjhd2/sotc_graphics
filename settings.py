#!/usr/local/sci/python
#************************************************************************
#
#  Global Settings for SotC plots
#  RJHD, Exeter 2017
#
#************************************************************************
#                    SVN Info
# $Rev:: 22                                         $:  Revision of last commit
# $Author:: rdunn                                      $:  Author of last commit
# $Date:: 2018-04-06 15:34:21 +0100 (Fri, 06 Apr #$:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import matplotlib.pyplot as plt
import numpy as np


#************************************************************************
def adjust_RdYlBu():
    '''
    Original RdYlBu has a greenish colour on the blue end, which doesn't contrast
    well with the yellow.  As this is a diverging colourmap, want to have a larger
    distinction for the positives and negatives.

    Take the colourmap, move the blues along a bit, and repeat at the end.
    '''

    print "for 2017/18 - change to match the purple adjustment"

    cmap = plt.cm.RdYlBu

    cmaplist = [cmap(i) for i in range(cmap.N)]

    # ignore first 20 of the blues
    N=20    
    retained_colours = np.array(cmaplist[255/2 + N :])

    stretched = np.ones((len(cmaplist[255/2:]), 4))

    # interpolate out
    for i in range(stretched.shape[1]):
        stretched[:, i] = np.interp(np.linspace(0, retained_colours.shape[0], stretched.shape[0]), range(retained_colours.shape[0]), retained_colours[:, i])

    # convert to tuple
    new_colours = []
    for c in stretched:
        new_colours += [tuple(c)]

    # copy back
    cmaplist[255/2:] = new_colours

    # enforce white for the centre of the colour range
    for i in [126,127,128,129]:
        cmaplist[i] = (1.0,1.0,1.0,1.0)

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
    N=20    
    retained_colours = np.array(cmaplist[255/2 + N :])

    stretched = np.ones((len(cmaplist[255/2:]), 4))

    # interpolate out
    for i in range(stretched.shape[1]):
        stretched[:, i] = np.interp(np.linspace(0, retained_colours.shape[0], stretched.shape[0]), range(retained_colours.shape[0]), retained_colours[:, i])

    # convert to tuple
    new_colours = []
    for c in stretched:
        new_colours += [tuple(c)]

    # copy back
    cmaplist[255/2:] = new_colours

    # enforce white for the centre of the colour range
    for i in [126,127,128,129]:
        cmaplist[i] = (1.0,1.0,1.0,1.0)

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

    cmaplist[255/2:] = cmaplist_R[255/2:]

    # enforce white for the centre of the colour range
    for i in [126,127,128,129]:
        cmaplist[i] = (1.0,1.0,1.0,1.0)

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
    for i in [126,127,128,129]:
        cmaplist[i] = (1.0,1.0,1.0,1.0)

    # and make the colour map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    cmap_r = cmap.from_list('Custom cmap', cmaplist[::-1], cmap.N)

    return cmap, cmap_r # make_BrBu


#************************************************************************
#************************************************************************
YEAR = "2017"
OUTFMT = ".png"
OUTFMT = ".eps"
#OUTFMT = ".pdf"


FONTSIZE = 20
LEGEND_FONTSIZE = 0.8 * FONTSIZE
LABEL_FONTSIZE = 0.9 * FONTSIZE

COLOURS = {"temperature" : {"ERA-Interim" : "orange", "MERRA-2" : "m", "JRA-55" : "c", "UAH v6.0" : "b", "RSS v4.0" : "r", "NOAA v3.0" : "c", "UW" : "m", "RAOBCORE v1.5" : "r", "RICH v1.5" : "y", "RATPAC A2" : "m", "UNSW v1.0" : "c", "CFSR" : "lime", "NOAA/NCEI" : "r", "NASA/GISS" : "b", "JMA" : "c", "Berkeley" : "y", "HadCRUT4" : "k", "CRUTEM" : "k", "HadSST" : "k", "GHCNDEX" : "r", "CMIP5" : "k", "SSU-3" : "r", "North" : "k", "South": "k", "QBO" : "k", "NOAA v4.0" : "k", "RSS v3.3" : "r"}, \
"cryosphere" : {"N Hemisphere" : 'k', "Eurasia" : 'r', "N America" : 'b', "Cumulative Balance" : "k", "Balance" : "r"},\
"hydrological" : {"ERA-Interim" : "purple", "JRA-55" : "c", "MERRA" : "lime", "MERRA-2" : "lime", "COSMIC RO" : "b", "GNSS (Ground Based)" : "r", "RSS Satellite" : '0.5', "HIRS" : "k", "Microwave" : "b", "GHCN" : "0.5", "GPCC" : "r", "GPCPv23" : "b", "GHCNv2" : "k", "GPCP" : "c", "PATMOS-x/AVHRR" : "0.5", "PATMOS-x/AQUA MODIS" : "c", "MISR" : "brown", "AQUA MODIS C6" : "b", "CALIPSO" : "r", "CERES" : "orange", "SatCORPS" : "lime", "CLARA-A2": "purple", "Cloud CCI": "g", "GRACE" : "k", "HadISDH" : "0.5", "HadCRUH" : "k", "HadCRUHExt" : "k", "Dai" : "r", "NOCS v2.0" : "b", "HOAPS" : "brown", "NCEP" : "w", "20CR" : "w", "Globe" : "k", "N. Hemisphere" : "c", "S. Hemisphere" : "m", "ERA-Interim mask" : "purple", "MERRA-2 mask" : "lime"},\
"circulation" : {"SSM/I+SSMIS" : "k", "NOCSv2.0" : "b", "WASwind" : "r", "ERA-Interim" : "orange", "JRA-55" : "c", "MERRA-2" : "m", "MERRA" : "m", "ERApreSAT" : "purple", "CERA20C" : "purple", "GRASP" : 'k', "GIUB" : "lime"},\
"radiation" : {"AT" : "k"},\
"composition" : {"AOD monthly" : "r", "AOD annual" : "b"},\
"land_surface" : {"Globe" : "k", "N. Hemisphere" : "b", "S. Hemisphere" : "r", "Globe Smoothed" : "k", "N. Hemisphere Smoothed" : "b", "S. Hemisphere Smoothed" : "r", "GFED4s" : "b", "GFASv1.4" : "r", "GFASv1.0" : "k"},\
"lst" : {"400" : "c", "300" : "m", "250" : "lime", "200" : "y", "150" : "k", "100" : "orange", "70" : "c", "50" : "m", "30" : "lime", "20" : "y", "10" : "k", "Average" : "k"},\
"phenological" : {"Chlorophyll-a" : "g", "B. pendula" : "b", "Q. robur (2000-{})".format(YEAR[2:]) : "r", "Q. robur" : "r", "F. sylvatica" : "k", "Q. robur (1951-99)" : "orange", "A. hippocastanum" : "c", "A. glutinosa": "m", "Greenup" : "lime", "Greendown" : "orange", "Difference" : "k", "Bartlett GPP" : "g", "Bartlett GCC" : "b", "Duke GCC" : "k"}} #  note space in second Q. robur



# all maps have 10 colours - or 11 with white as central.  Can probably make up from Kate's code if better match needed.


RdYlBu, RdYlBu_r = adjust_RdYlBu()
PuOr, PuOr_r = adjust_PuOr()
# Composition & Land Surface - create a version of Brown-Blue brewer colours
BrBu, BrBu_r = make_BrBu()
# adjust BrBG to ensure has central colours equal to white
BrBG, BrBG_r = make_BrBG()

COLOURMAP_DICT = {"temperature" : RdYlBu_r, "temperature_r" : RdYlBu, "hydrological" : BrBG, "hydrological_r" : BrBG_r, "precip_sequential" : plt.cm.YlGnBu, "circulation" : PuOr, "circulation_r" : PuOr_r, "composition" : BrBu_r, "composition_r" : BrBu, "land_surface" : BrBu_r, "land_surface_r" : BrBu, "phenological" : BrBG, "phenological_r" : BrBG_r}
