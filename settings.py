#!/usr/local/sci/python
#************************************************************************
#
#  Global Settings for SotC plots
#  RJHD, Exeter 2017
#
#************************************************************************
#                    SVN Info
# $Rev:: 20                                         $:  Revision of last commit
# $Author:: rdunn                                      $:  Author of last commit
# $Date::                                         $:  Date of last commit
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

    # by hand move the blues further down the colourmap
    i = 0
    while i < 255/2 - 20: # visual testing that 20 increments covers the green bit

        cmaplist[255/2 + i] = cmaplist[255/2 + i + 20]
        i += 1

    # enforce white for the centre of the colour range
    for i in [126,127,128,129]:
        cmaplist[i] = (1.0,1.0,1.0,1.0)

    # by hand move repeat the last colour for the remainder
    cmaplist[-1] = (0.14,0.16,0.47,1.0)
    i -= 1
    while i < 255/2:

        cmaplist[255/2 + i] = cmaplist[-1]
        i += 1

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

    # and make the colour map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    cmap_r = cmap.from_list('Custom cmap', cmaplist[::-1], cmap.N)

    return cmap, cmap_r # make_BrBu

YEAR = "2016"
OUTFMT = ".png"
#OUTFMT = ".eps"
#OUTFMT = ".pdf"


FONTSIZE = 20
LEGEND_FONTSIZE = 0.8 * FONTSIZE
LABEL_FONTSIZE = 0.9 * FONTSIZE

COLOURS = {"temperature" : {"ERA-Interim" : "orange", "MERRA-2" : "m", "JRA-55" : "c", "UAH v6.0" : "b", "RSS v3.3" : "r", "NOAA v3.0" : "c", "UW" : "m", "RAOBCORE v1.5" : "r", "RICH v1.5" : "y", "RATPAC A2" : "m", "UNSW v1.0" : "c", "CFSR" : "lime", "NOAA/NCEI" : "r", "NASA/GISS" : "b", "JMA" : "c", "Berkeley" : "y", "HadCRUT4" : "k", "CRUTEM" : "k", "HadSST" : "k", "GHCNDEX" : "r"}, \
"cryosphere" : {"N Hemisphere" : 'k', "Eurasia" : 'r', "N America" : 'b', "Cumulative Balance" : "k", "Balance" : "r"},\
"hydrological" : {"ERA-Interim" : "purple", "JRA-55" : "c", "MERRA" : "lime", "MERRA-2" : "lime", "COSMIC RO" : "b", "GNSS (Ground Based)" : "r", "RSS Satellite" : '0.5', "HIRS" : "k", "Microwave" : "b", "GHCN" : "0.5", "GPCC" : "r", "GPCPv23" : "b", "GHCNv2" : "k", "GPCP" : "c", "PATMOS-x/AVHRR" : "0.5", "MISR" : "brown", "PATMOS-x/AQUA MODIS C6" : "b", "CALIPSO" : "r", "CERES" : "orange", "SatCORPS" : "lime", "CLARA-A2": "purple", "GRACE" : "k", "HadISDH" : "0.5", "HadCRUH" : "k", "HadCRUHExt" : "k", "Dai" : "r", "NOCS v2.0" : "b", "HOAPS" : "brown", "NCEP" : "w", "20CR" : "w", "Globe" : "k", "N. Hemisphere" : "b", "S. Hemisphere" : "r"},\
"circulation" : {"SSM/I+SSMIS" : "k", "NOCSv2.0" : "b", "WASwind" : "r", "ERA-Interim" : "orange", "JRA-55" : "c", "MERRA-2" : "m", "MERRA" : "m", "ERApreSAT" : "purple", "GRASP" : 'k', "GIUB" : "lime"},\
"radiation" : {"AT" : "k"},\
"composition" : {"AOD monthly" : "r", "AOD annual" : "b"},\
"land_surface" : {"Globe" : "k", "N. Hemisphere" : "b", "S. Hemisphere" : "r", "Globe Smoothed" : "k", "N. Hemisphere Smoothed" : "b", "S. Hemisphere Smoothed" : "r", "GFED3.1" : "g", "GFASv1.3" : "b", "GFASv1.0" : "k"},\
"lst" : {"400" : "c", "300" : "m", "250" : "lime", "200" : "y", "150" : "k", "100" : "orange", "70" : "c", "50" : "m", "30" : "lime", "20" : "y", "10" : "k", "Average" : "k"}}




# all maps have 10 colours - or 11 with white as central.  Can probably make up from Kate's code if better match needed.


RdYlBu, RdYlBu_r = adjust_RdYlBu()
PuOr, PuOr_r = adjust_PuOr()
# Composition & Land Surface - create a version of Brown-Blue brewer colours
BrBu, BrBu_r = make_BrBu()



COLOURMAP_DICT = {"temperature" : RdYlBu_r, "temperature_r" : RdYlBu, "hydrological" : plt.cm.BrBG, "hydrological_r" : plt.cm.BrBG_r, "circulation" : PuOr, "circulation_r" : PuOr_r, "composition" : BrBu_r, "composition_r" : BrBu, "land_surface" : BrBu_r, "land_surface_r" : BrBu}
