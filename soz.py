#!/usr/bin/env python
#************************************************************************
#
#  Plot figures and output numbers for Stratospheric Ozone (SOZ) section.
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

import utils # RJHD utilities
import settings

DATALOC = "{}/{}/data/SOZ/".format(settings.ROOTLOC, settings.YEAR)

#************************************************************************
def read_data(filename):
    """
    Read data from .txt file into Iris cube

    :param str filename: file to process

    :returns: cube
    """
    # use header to hard code the final array shapes
    longitudes = np.arange(0.625, 360.625, 1.25)
    latitudes = np.arange(-89.5, 90.5, 1.)
    
    data = np.ma.zeros((latitudes.shape[0], longitudes.shape[0]))

    # read in the dat
    indata = np.genfromtxt(filename, dtype=(float), skip_header=4)

    this_lat = []
    tl = 0
    # process each row, append until have complete latitude band
    for row in indata:
        this_lat += [r for r in row]

        if len(this_lat) == longitudes.shape[0]:
            # copy into final array and reset
            data[tl, :] = this_lat
            tl += 1
            this_lat = []

    # mask the missing values
    data = np.ma.masked_where(data <= -999.000, data)

    cube = utils.make_iris_cube_2d(data, latitudes, longitudes, "SOZ_anom", "DU")

    return cube # read_data

#************************************************************************
def read_ts(filename):

    indata = np.genfromtxt(filename, skip_header=1)

    years = indata[:, 0]
    mean_val = indata[:, 1]

    ts = utils.Timeseries("SOZ", years, mean_val)

    return ts

#************************************************************************
def read_multi_ts(filename):

    all_series = []
    years = []
    values = []

    with open(filename, "r") as infile:
        for line in (infile):
            line = line.split()

            # skip comments
            if line[0] == "#":
                # move on if nothing to process
                if len(years) == 0:
                    continue

                try:
                    # mask bad values
                    values = np.ma.masked_where(np.array(values) < 100, np.array(values))
                    years = np.array(years)

                    # expand years to full range
                    all_years = np.arange(years[0], years[-1]+1)
                    all_values = np.ma.zeros(all_years.shape)
                    all_values.mask = np.ones(all_years.shape)
                    locs = np.in1d(all_years, years)

                    all_values[locs] = values
                    all_values.mask[locs] = False

                    all_series += [utils.Timeseries("SOZ", all_years, all_values)]
                except TypeError:
                    pass

                years = []
                values = []
            else:
                years += [int(line[0])]
                values += [float(line[1])]    
                
        try:
            # mask bad values
            values = np.ma.masked_where(np.array(values) < 100, np.array(values))
            years = np.array(years)

            # expand years to full range
            all_years = np.arange(years[0], years[-1]+1)
            all_values = np.ma.zeros(all_years.shape)
            all_values.mask = np.ones(all_years.shape)
            locs = np.in1d(all_years, years)

            all_values[locs] = values
            all_values.mask[locs] = False

            all_series += [utils.Timeseries("SOZ", all_years, all_values)]
        except TypeError:
            pass

    return all_series # read_multi_ts



#************************************************************************
def run_all_plots():

    #************************************************************************
    # Global Map
    if True:
        cube = read_data(DATALOC + "GSG_{}annual mean_G_ano_ref1998-2008.txt".format(settings.YEAR))

        bounds = [-100, -20, -15, -8, -4, 0, 4, 8, 15, 20, 100]
        bounds = [-200, -90, -40, -15, -5, -1, 1, 5, 15, 40, 90, 200]


        utils.plot_smooth_map_iris(settings.IMAGELOC + "SOZ_{}_anoms".format(settings.YEAR), cube, settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 1998-2008 (DU)")

        utils.plot_smooth_map_iris(settings.IMAGELOC + "p2.1_SOZ_{}_anoms".format(settings.YEAR), cube, settings.COLOURMAP_DICT["composition"], bounds, "Anomalies from 1998-2008 (DU)", figtext="(aa) Stratospheric (Total Column) Ozone")


    #************************************************************************
    # Timeseries (plate 1.1)
    if True:
    #    nh = read_ts(DATALOC + "tozave_60N_90N_3_3.dat")
    #    sh = read_ts(DATALOC + "tozave_90S_60S_10_10.dat")

    #    nh = read_ts(DATALOC + "polar_ozone_NH.dat")
    #    sh = read_ts(DATALOC + "polar_ozone_SH.dat")

        nh_series = read_multi_ts(DATALOC + "polar_ozone_NH.dat")
        sh_series = read_multi_ts(DATALOC + "polar_ozone_SH.dat")

        fig = plt.figure(figsize=(8, 6))
        plt.clf()

        for nh in nh_series:
            plt.plot(nh.times, nh.data, "b-", label="March NH")
        for sh in sh_series:
            plt.plot(sh.times, sh.data, "r-", label="October SH")

        plt.ylabel("total Ozone (DU")
        plt.title("Polar Ozone")

        plt.legend()
        plt.savefig(settings.IMAGELOC + "SOZ_ts{}".format(settings.OUTFMT))

    return # run_all_plots

#************************************************************************
if __name__ == "__main__":

    run_all_plots()

#************************************************************************
#                                 END
#************************************************************************
